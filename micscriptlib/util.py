#!/usr/bin/env python
#
import argparse
import glob
import os
import shutil
import subprocess
import sys
from os.path import join as pjoin

import rapidtide.io as tide_io


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


def make_runscript(thecommand, jobname="rapidtide", ncpus=8, timelimit="0:30:00", mem="16G"):
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

    script = "\n".join(pre) + "\n" + " \\\n    ".join(thecommand) + "\n"

    _, filename = tempfile.mkstemp()
    with open(filename, "w") as fp:
        fp.write(script)
    return filename, script


atlasaveragecommand = "/cm/shared/anaconda3/envs/mic/bin/atlasaverage"
fingerprintcommand = "/cm/shared/anaconda3/envs/mic/bin/fingerprint"
fsldir = os.environ.get("FSLDIR")
if fsldir is not None:
    fslsubcmd = os.path.join(fsldir, "bin", "fsl_sub")
    flirtcmd = os.path.join(fsldir, "bin", "flirt")
    applywarpcmd = os.path.join(fsldir, "bin", "applywarp")
    fslexists = True
else:
    fslexists = False


def runcmd(thecmd, cluster=False, readable=False, fake=False, waitfor=None, debug=False):
    SYSTYPE, SUBMITTER, SINGULARITY = getbatchinfo()
    if debug:
        print(thecmd)
    if fake:
        if readable:
            print(thecmd[0])
            for thearg in thecmd[1:]:
                print("\t", thearg)
        else:
            print(" ".join(thecmd))
        print()
    else:
        if cluster:
            jobname = thecmd[0].split("/")[-1]
            scriptfile, thescript = make_runscript(thecmd, jobname)
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
            thereturn = subprocess.check_output(sub_cmd).split()
            thepid = (thereturn[-1]).strip()
            thepidstr = str(thepid, "UTF8")
            if debug:
                print("return value:", thereturn)
                print("pid value:", thepid)
                print("pidstr value:", thepidstr)
            return str(thepidstr)
        else:
            subprocess.call(thecmd)
            return None


def mriconvert(inputfile, outputfile, cluster=False, fake=False, waitfor=None, debug=False):
    convcmd = []
    convcmd += ["mri_convert"]
    convcmd += [inputfile]
    convcmd += [outputfile]
    pidnum = runcmd(convcmd, cluster=cluster, fake=fake, waitfor=waitfor, debug=debug)
    return pidnum


def n4correct(inputfile, outputdir, cluster=False, fake=False, waitfor=None, debug=False):
    thename, theext = tide_io.niftisplitext(inputfile)
    n4cmd = []
    n4cmd += ["N4BiasFieldCorrection"]
    n4cmd += ["-d", "3"]
    n4cmd += ["-i", inputfile]
    n4cmd += ["-o", pjoin(outputdir, thename + "_n4" + theext)]
    pidnum = runcmd(n4cmd, cluster=cluster, fake=fake, waitfor=waitfor, debug=debug)
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
    applyxfmcmd += ["/cm/shared/apps/ants-2.4.3/bin/antsApplyTransforms"]
    applyxfmcmd += ["--default-value", "0"]
    applyxfmcmd += ["-d", "3"]
    applyxfmcmd += ["-i", inputname]
    applyxfmcmd += ["-o", outputroot]
    applyxfmcmd += ["-r", targetname]
    if interp is not None:
        applyxfmcmd += ["--interpolation", interp]
    for thetransform in transforms:
        applyxfmcmd += ["--transform", thetransform]
    pidnum = runcmd(applyxfmcmd, cluster=cluster, fake=fake, waitfor=waitfor, debug=debug)
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
    fingerprintcmd += [atlasaveragecommand]
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
    pidnum = runcmd(fingerprintcmd, cluster=cluster, fake=fake, waitfor=waitfor, debug=debug)
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
    atlasavgcmd += [atlasaveragecommand]
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
    pidnum = runcmd(atlasavgcmd, cluster=cluster, fake=fake, waitfor=waitfor, debug=debug)
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


def parseconnectomename(thefile, volumeproc=True):
    absname = os.path.abspath(thefile)
    therundir, thefmrifile = os.path.split(absname)
    theresultsdir, therun = os.path.split(therundir)
    if volumeproc:
        theMNINonLinDir, dummy = os.path.split(theresultsdir)
        thesubjdir, dummy = os.path.split(theMNINonLinDir)
        dummy, thesubj = os.path.split(thesubjdir)
        pedir = therun[-2:]
    else:
        splitname = thefmrifile.split("_")
        thesubj = therun
        therun = "_".join(splitname[0:2])
        pedir = splitname[2]
    return absname, thesubj, therun, pedir


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


def findboldfiles_psusleep(
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
    debug=False,
    bidsroot="/data/ajanes/REcig/fmri",
):
    if inputlistfile is None:
        searchstring = os.path.join(
            bidsroot,
            "Recig_*",
            "visit*",
            thetype,
            f"*{thetype}*visit[12].feat",
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
                f"*{thetype}*visit[12].feat",
                "filtered_func_data.nii.gz",
            )
            #    f"*{thetype}*{space}*bold.nii.gz",
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


def findboldfiles_fmriprep(
    inputlistfile=None,
    debug=False,
    space="MNI152NLin6Asym",
    bidsroot="/data/frederic/OASIS",
):
    if inputlistfile is None:
        searchstring = os.path.join(
            bidsroot,
            "derivatives",
            "fmriprep",
            "sub*",
            "ses*",
            "func",
            f"*{space}*bold.nii.gz",
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
                f"*{space}*bold.nii.gz",
            )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


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
