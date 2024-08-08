#!/usr/bin/env python
#
import glob
import os
import site
import subprocess
import sys
from os.path import join as pjoin

import rapidtide.io as tide_io

pythonbindir = "/cm/shared/anaconda3/envs/mic/bin"

fsldir = os.environ.get("FSLDIR")
if fsldir is not None:
    fslsubcmd = os.path.join(fsldir, "bin", "fsl_sub")
    flirtcmd = os.path.join(fsldir, "bin", "flirt")
    applywarpcmd = os.path.join(fsldir, "bin", "applywarp")
    fslexists = True
else:
    fslexists = False

freesurferpath = os.path.join(os.environ.get("FREESURFER_HOME"), "bin")
if freesurferpath is not None:
    freesurferexists = True
else:
    freesurferexists = False

antspath = os.environ.get("ANTSPATH")
if antspath is not None:
    antsexists = True
else:
    antsexists = False


def getbatchinfo():
    QSUB = "/cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub"
    SBATCH = "/cm/shared/apps/slurm/current/bin/sbatch"
    if os.path.isfile(QSUB):
        SYSTYPE = "sge"
        SUBMITTER = QSUB
        SINGULARITY = "/usr/bin/singularity"
    else:
        SYSTYPE = "slurm"
        SUBMITTER = SBATCH
        SINGULARITY = "/cm/local/apps/apptainer/current/bin/singularity"
    return SYSTYPE, SUBMITTER, SINGULARITY


def make_runscript(thecommand, jobname="rapidtide", ncpus=8, timelimit="0:02:00", mem="1G", debug=False):
    """
    Create a temporary script file we can submit to qsub.
    """

    import tempfile

    pre = ["#!/bin/bash"]

    # pre += ["#SBATCH -p short"]
    pre += [f"#SBATCH --job-name={jobname}"]
    pre += [f"#SBATCH --output={jobname}.%j.out"]
    pre += [f"#SBATCH --error={jobname}.%j.err"]
    pre += [f"#SBATCH --time={timelimit}"]
    pre += [f"#SBATCH --cpus-per-task={ncpus}"]
    pre += ["#SBATCH --qos=normal"]
    pre += [f"#SBATCH --mem={mem}"]
    pre += ["#SBATCH --mail-type=FAIL # notifications for job fail"]
    pre += ["#SBATCH --mail-user=bbfrederick@mclean.harvard.edu"]

    if debug:
        print(f"thecommand: {thecommand}")
    script = "\n".join(pre) + "\n" + " \\\n    ".join(thecommand) + "\n"

    _, filename = tempfile.mkstemp()
    with open(filename, "w") as fp:
        fp.write(script)
    return filename, script


def runcmd(
    thecmd,
    timelimit="0:02:00",
    mem="1G",
    ncpus=1,
    cluster=False,
    readable=False,
    fake=False,
    waitfor=None,
    debug=False,
):
    SYSTYPE, SUBMITTER, SINGULARITY = getbatchinfo()
    if debug:
        print("RUNCMD:", thecmd)
        print(f"\t{cluster=}")
        print(f"\t{readable=}")
        print(f"\t{fake=}")
        print(f"\t{waitfor=}")
    else:
        if cluster:
            jobname = "log_" + thecmd[0].split("/")[-1]
            scriptfile, thescript = make_runscript(
                thecmd, jobname, ncpus=ncpus, timelimit=timelimit, mem=mem
            )
            waitstr = ""
            if SYSTYPE == "sge":
                if waitfor is not None:
                    waitstr = f"-j {waitfor}"
                sub_cmd = f"{SUBMITTER} -cwd -q fmriprep.q -N {jobname} -w e -R y {waitstr} {thecmd}".split()
            elif SYSTYPE == "slurm":
                if waitfor is not None:
                    waitstr = f"--dependency afterok:{waitfor}"
                sub_cmd = f"{SUBMITTER} {waitstr} {scriptfile}".split()
            if debug:
                print("sub_cmd:", sub_cmd)
            if fake:
                print(sub_cmd)
                print(thescript)
                return None
            else:
                thereturn = subprocess.check_output(sub_cmd).split()
                thepid = (thereturn[-1]).strip()
                thepidstr = str(thepid, "UTF8")
                if debug:
                    print("return value:", thereturn)
                    print("pid value:", thepid)
                    print("pidstr value:", thepidstr)
                return str(thepidstr)
        else:
            if fake:
                if readable:
                    print(thecmd[0])
                    for thearg in thecmd[1:]:
                        print("\t", thearg)
                else:
                    print(" ".join(thecmd))
                print()
            else:
                subprocess.call(thecmd)
            return None


def mriconvert(inputfile, outputfile, cluster=False, fake=False, waitfor=None, debug=False):
    convcmd = []
    convcmd += [f"{freesurferpath}/mri_convert"]
    convcmd += [inputfile]
    convcmd += [outputfile]
    pidnum = runcmd(
        convcmd,
        timelimit="0:02:00",
        mem="1G",
        ncpus=1,
        cluster=cluster,
        fake=fake,
        waitfor=waitfor,
        debug=debug,
    )
    return pidnum


def n4correct(inputfile, outputdir, cluster=False, fake=False, waitfor=None, debug=False):
    thename, theext = tide_io.niftisplitext(inputfile)
    n4cmd = []
    n4cmd += [f"{antspath}/N4BiasFieldCorrection"]
    n4cmd += ["-d", "3"]
    n4cmd += ["-i", inputfile]
    n4cmd += ["-o", pjoin(outputdir, thename + "_n4" + theext)]
    pidnum = runcmd(
        n4cmd,
        timelimit="0:02:00",
        mem="1G",
        ncpus=1,
        cluster=cluster,
        fake=fake,
        waitfor=waitfor,
        debug=debug,
    )
    return pidnum


def antsapply(
    inputname,
    targetname,
    outputroot,
    transforms,
    cluster=False,
    fake=False,
    debug=False,
    waitfor=None,
    interp=None,
):
    applyxfmcmd = []
    applyxfmcmd += [f"{antspath}/antsApplyTransforms"]
    applyxfmcmd += ["--default-value", "0"]
    applyxfmcmd += ["-d", "3"]
    applyxfmcmd += ["-i", inputname]
    applyxfmcmd += ["-o", outputroot]
    applyxfmcmd += ["-r", targetname]
    if interp is not None:
        applyxfmcmd += ["--interpolation", interp]
    for thetransform in transforms:
        applyxfmcmd += ["--transform", thetransform]
    pidnum = runcmd(
        applyxfmcmd,
        timelimit="0:02:00",
        mem="1G",
        ncpus=1,
        cluster=cluster,
        fake=fake,
        waitfor=waitfor,
        debug=debug,
    )
    return pidnum


def fingerprintapply(
    inputfile,
    outputfileroot,
    includemask=None,
    excludemask=None,
    extramask=None,
    atlas="JHU1",
    fitorder=1,
    template="lag",
    limittomask=False,
    nointercept=None,
    cluster=False,
    fake=False,
    waitfor=None,
    debug=False,
):
    fingerprintcmd = []
    fingerprintcmd += [f"{pythonbindir}/fingerprint"]
    fingerprintcmd += [inputfile]
    fingerprintcmd += [outputfileroot]
    fingerprintcmd += ["--atlas", atlas]
    fingerprintcmd += ["--fitorder", str(fitorder)]
    fingerprintcmd += ["--template", template]
    if limittomask:
        fingerprintcmd += ["--limittomask"]
    if nointercept:
        fingerprintcmd += ["--nointercept"]
    if includemask is not None:
        fingerprintcmd += ["--includemask", includemask]
    if excludemask is not None:
        fingerprintcmd += ["--excludemask", excludemask]
    if extramask is not None:
        fingerprintcmd += ["--extramask", extramask]
    pidnum = runcmd(
        fingerprintcmd,
        timelimit="0:02:00",
        mem="1G",
        ncpus=1,
        cluster=cluster,
        fake=fake,
        waitfor=waitfor,
        debug=debug,
    )
    return pidnum


def runqualitycheckapply(
    inputfileroot,
    graymaskspec=None,
    whitemaskspec=None,
    cluster=False,
    fake=False,
    waitfor=None,
    debug=False,
):
    runqualcmd = []
    runqualcmd += [f"{pythonbindir}/runqualitycheck"]
    runqualcmd += [inputfileroot]
    if graymaskspec is not None:
        runqualcmd += ["--graymaskspec", graymaskspec]
    if whitemaskspec is not None:
        runqualcmd += ["--whitemaskspec", whitemaskspec]
    pidnum = runcmd(
        runqualcmd,
        timelimit="0:02:00",
        mem="1G",
        ncpus=1,
        cluster=cluster,
        waitfor=waitfor,
        fake=fake,
        debug=debug,
    )
    return pidnum


def atlasaverageapply(
    inputfile,
    templatefile,
    outputfileroot,
    regionlist=None,
    label=None,
    includemask=None,
    excludemask=None,
    extramask=None,
    nozero=False,
    header=False,
    cluster=False,
    fake=False,
    waitfor=None,
    debug=False,
):
    atlasavgcmd = []
    atlasavgcmd += [f"{pythonbindir}/atlasaverage"]
    atlasavgcmd += [inputfile]
    atlasavgcmd += [templatefile]
    atlasavgcmd += [outputfileroot]
    if nozero:
        atlasavgcmd += ["--ignorezeros"]
    if header:
        atlasavgcmd += ["--headerline"]
    if label is not None:
        atlasavgcmd += ["--datalabel", label]
    if regionlist is not None:
        atlasavgcmd += ["--regionlistfile", regionlist]
    if includemask is not None:
        atlasavgcmd += ["--includemask", includemask]
    if excludemask is not None:
        atlasavgcmd += ["--excludemask", excludemask]
    if extramask is not None:
        atlasavgcmd += ["--extramask", extramask]
    pidnum = runcmd(
        atlasavgcmd,
        timelimit="0:02:00",
        mem="1G",
        ncpus=1,
        cluster=cluster,
        fake=fake,
        waitfor=waitfor,
        debug=debug,
    )
    return pidnum


def readlist(inputfilename):
    inputlist = []
    with open(inputfilename, "r") as thefile:
        lines = thefile.readlines()
        for line in lines:
            inputlist.append(line.strip())
    return inputlist


def makeadir(pathname):
    try:
        os.makedirs(pathname)
    except OSError:
        if os.path.exists(pathname):
            # We are nearly safe
            return True
        else:
            # There was an error on creation, so make sure we know about it
            print("ERROR: ", pathname, " does not exist, and could not create it")
            return False
    return True


def parseconnectomename(thefile, volumeproc=True, debug=False):
    absname = os.path.abspath(thefile)
    therundir, thefmrifile = os.path.split(absname)
    theresultsdir, therun = os.path.split(therundir)
    splitname = thefmrifile.split("_")
    therun = "_".join(splitname[0:2])
    pedir = splitname[2].split(".")[0]
    theMNINonLinDir, dummy = os.path.split(theresultsdir)
    if debug:
        print(f"{absname=}")
        print(f"{therundir=}")
        print(f"{thefmrifile=}")
        print(f"{theresultsdir=}")
        print(f"{theMNINonLinDir=}")
    if volumeproc:
        thesubjdir, dummy = os.path.split(theMNINonLinDir)
        dummy, thesubj = os.path.split(thesubjdir)
    else:
        print("NOT IMPLEMENTED!")
        thesubj = None
    if debug:
        print(f"{thesubj=}")
        print(f"{therun=}")
        print(f"{pedir=}")
    return absname, thesubj, therun, pedir, theMNINonLinDir


def parsecolename(thefile, volumeproc=True):
    absname = os.path.abspath(thefile)
    therundir, thefmrifile = os.path.split(absname)
    theparts = thefmrifile.split("_")
    thesubj = theparts[0]
    therun = theparts[3]
    pedir = theparts[-1][0:2]
    return absname, thesubj, therun, pedir


def findboldfiles_HCP(theroot, thetype, volumeproc, usefixforglm, inputlistfile=None, debug=False):
    # special options depending on whether using volume or grayordinate files
    if inputlistfile is None:
        if volumeproc:
            searchstring = os.path.join(
                theroot,
                "*",
                "preproc",
                "*",
                "MNINonLinear",
                "Results",
                thetype + "_[LR][LR]",
                thetype + "_[LR][LR].nii.gz",
            )
        else:
            searchstring = os.path.join(theroot, "*", thetype + "_[LR][LR]_grayordinate.nii.gz")
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            if volumeproc:
                searchstring = os.path.join(
                    theroot,
                    str(subject),
                    "preproc",
                    "*",
                    "MNINonLinear",
                    "Results",
                    thetype + "_[LR][LR]",
                    thetype + "_[LR][LR].nii.gz",
                )
            else:
                searchstring = os.path.join(
                    theroot, str(subject), thetype + "_[LR][LR]_grayordinate.nii.gz"
                )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def parserecigname(thefile):
    absname = os.path.abspath(thefile)
    therundir, thefilename = os.path.split(absname)
    thetopdir, therunkey = os.path.split(therundir)
    theparts = therunkey.split("_")
    thesubj = theparts[0]
    thetask = theparts[1]
    theses = (theparts[2].split("."))[0]

    return absname, thefilename, thesubj, theses, thetask


def parsebidsname(thefile):
    absname = os.path.abspath(thefile)
    therundir, thefilename = os.path.split(absname)
    theparts = thefilename.split("_")
    nameparts = {}
    for thepart in theparts:
        splitparts = thepart.split("-")
        if len(splitparts) == 2:
            nameparts[splitparts[0]] = splitparts[1]
    try:
        thesubj = nameparts["sub"]
    except KeyError:
        print("no subject key!")
        sys.exit()
    try:
        thetask = nameparts["task"]
    except KeyError:
        thetask = None
    try:
        theses = nameparts["ses"]
    except KeyError:
        theses = None
    try:
        therun = nameparts["run"]
    except KeyError:
        therun = None
    try:
        pedir = nameparts["dir"]
    except KeyError:
        pedir = None
    try:
        thespace = nameparts["space"]
    except KeyError:
        thespace = None

    return absname, thefilename, thesubj, theses, therun, pedir, thetask, thespace


def files_psusleep(
    thetype,
    inputlistfile=None,
    debug=False,
    space="MNI152NLin6Asym",
    bidsroot="/data/ckorponay/Sleep",
):
    if inputlistfile is None:
        searchstring = os.path.join(
            bidsroot,
            "fmriprep",
            "sub*",
            "func",
            f"*{thetype}*smoothAROMAnonaggr_bold.nii.gz",
        )
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            searchstring = os.path.join(
                bidsroot,
                "fmriprep",
                "sub-" + str(subject),
                "func",
                f"*{thetype}*smoothAROMAnonaggr_bold.nii.gz",
            )
            #    f"*{thetype}*{space}*bold.nii.gz",
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def findboldfiles_recig(
    thetype,
    inputlistfile=None,
    altpath=False,
    debug=False,
    bidsroot="/data/ajanes/REcig/fmri",
):
    if (thetype != "resting") and altpath:
        extrapath = "_SpecRegressors"
    else:
        extrapath = ""
    if inputlistfile is None:
        searchstring = os.path.join(
            bidsroot,
            "Recig_*",
            "visit*",
            thetype,
            f"*{thetype}*visit[12]{extrapath}.feat",
            "filtered_func_data.nii.gz",
        )
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            searchstring = os.path.join(
                bidsroot,
                "Recig_" + str(subject),
                "visit*",
                thetype,
                f"*{thetype}*visit[12]{extrapath}.feat",
                "filtered_func_data.nii.gz",
            )
            #    f"*{thetype}*{space}*bold.nii.gz",
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def findboldfiles_fmriprep(
    subjects=None,
    sessions=None,
    tasks=None,
    space="MNI152NLin6Asym",
    bidsroot="/data/frederic/OASIS",
    hassessions=True,
    debug=False,
):
    # find all bold files
    if hassessions:
        pathparts = [
            bidsroot,
            "derivatives",
            "fmriprep",
            "sub*",
            "ses*",
            "func",
            "*bold.nii.gz",
        ]
    else:
        pathparts = [
            bidsroot,
            "derivatives",
            "fmriprep",
            "sub*",
            "func",
            "*bold.nii.gz",
        ]
    searchstring = os.path.join(*pathparts)
    completelist = glob.glob(searchstring)
    if debug:
        print(completelist)

    # now filter
    if subjects is None:
        thesubjs = completelist
    else:
        thesubjs = [s for s in completelist if any(f"sub-{xs}" in s for xs in subjects)]
    if debug:
        print(thesubjs)

    if hassessions:
        if sessions is None:
            thesessions = thesubjs
        else:
            thesessions = [s for s in thesubjs if any(f"ses-{xs}" in s for xs in sessions)]
    if debug:
        print(thesessions)

    if tasks is None:
        thetasks = thesessions
    else:
        thetasks = [s for s in thesessions if any(f"task-{xs}" in s for xs in tasks)]
    if debug:
        print(thetasks)

    thespace = f"space-{space}"
    return [s for s in thetasks if thespace in s]




    """if subject is None:
        pathparts.append("sub-*")
    else:
        # loop over subjects
        for subject in subjects:
            pathparts.append(subject)
            if hassessions:
                # loop over sessions
                if sessions is None:
                    pathparts.append("ses-*")
                else:
                    for session in sessions:
                        pathparts.append(f"ses-*")

        if session is None:
            searchstring = os.path.join(
                bidsroot,
                "derivatives",
                "fmriprep",
                "sub*",
                "ses*",
                "func",
                f"*task-{task}_*{space}*bold.nii.gz",
            )
        searchstring = os.path.join(
            bidsroot,
            "derivatives",
            "fmriprep",
            "sub*",
            "ses*",
            "func",
            f"*task-{task}_*{space}*bold.nii.gz",
        )

        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            searchstring = os.path.join(
                bidsroot,
                "derivatives",
                "fmriprep",
                "sub-" + str(subject),
                "ses*",
                "func",
                f"*task-{task}_*{space}*bold.nii.gz",
            )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]"""


def findaparcaseg_fmriprep(
    inputlistfile=None,
    debug=False,
    bidsroot="/data/frederic/OASIS",
):
    if inputlistfile is None:
        searchstring = os.path.join(
            bidsroot,
            "derivatives",
            "fmriprep",
            "sub*",
            "anat",
            "*desc-aparcaseg_dseg.nii.gz",
        )
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            searchstring = os.path.join(
                bidsroot,
                "derivatives",
                "fmriprep",
                "sub-" + str(subject),
                "anat",
                "*desc-aparcaseg_dseg.nii.gz",
            )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def findboldfiles_BIDS_multisession(
    inputlistfile=None,
    debug=False,
    bidsroot="/data/frederic/OASIS",
):
    if inputlistfile is None:
        searchstring = os.path.join(
            bidsroot,
            "sub*",
            "ses*",
            "func",
            "*bold.nii.gz",
        )
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            searchstring = os.path.join(
                bidsroot,
                "sub-" + str(subject),
                "ses*",
                "func",
                "*bold.nii.gz",
            )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def findboldfiles_cole(
    theroot, thetype, volumeproc, usefixforglm, inputlistfile=None, debug=False
):
    # special options depending on whether using volume or grayordinate files
    if inputlistfile is None:
        if volumeproc:
            searchstring = os.path.join(
                theroot,
                "*",
                "*_Clean_" + thetype + "_[LR][LR].nii.gz",
            )
        else:
            searchstring = os.path.join(theroot, "*", thetype + "_[LR][LR]_grayordinate.nii.gz")
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            if volumeproc:
                searchstring = os.path.join(
                    theroot,
                    str(subject),
                    "preproc",
                    "*",
                    "MNINonLinear",
                    "Results",
                    thetype + "_[LR][LR]",
                    thetype + "_[LR][LR].nii.gz",
                )
                searchstring = os.path.join(
                    theroot,
                    str(subject),
                    "*_Clean_" + thetype + "_[LR][LR].nii.gz",
                )
            else:
                searchstring = os.path.join(
                    theroot, str(subject), thetype + "_[LR][LR]_grayordinate.nii.gz"
                )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def parsertname(thefile, debug=False):
    absname = os.path.abspath(thefile)
    thesessdir, thefmrifile = os.path.split(absname)
    thesubjdir, thesess = os.path.split(thesessdir)
    therootdir, thesubj = os.path.split(thesubjdir)
    therun = thefmrifile.split("_")[-3]
    if debug:
        print(thefmrifile, thesess, thesubj, therun)
    return absname, thesubj, thesess, therun


def findrtfiles(theroot, themaptype, inputlistfile=None, debug=False):
    if inputlistfile is None:
        searchstring = os.path.join(
            theroot,
            "sub-*",
            "ses-*",
            f"*_desc-{themaptype}_map.nii.gz",
        )
        if debug:
            print("searchstring:", searchstring)
        return glob.glob(searchstring)
    else:
        print("using subject list")
        inputlist = readlist(inputlistfile)
        print("subject list is ", inputlist)
        retlist = []
        for subject in inputlist:
            searchstring = os.path.join(
                theroot,
                f"sub-OAS{subject}",
                "ses-*",
                f"*_desc-{themaptype}_map.nii.gz",
            )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def finddatadir():
    # Get the list of directories
    site_packages_dirs = site.getsitepackages()

    # Find the "site-packages" directory in the list
    for dir in site_packages_dirs:
        if dir.endswith("site-packages"):
            sitepackages_dir = dir
            break
        else:
            sitepackages_dir = None
    referencedir = os.path.join(
        sitepackages_dir,
        "micscriptlib",
        "data",
    )
    return referencedir
