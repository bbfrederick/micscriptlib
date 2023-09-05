#!/usr/bin/env python
#
import argparse
import glob

import sys
import shutil

import os
import subprocess
from os.path import join as pjoin


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
fsldir = os.environ.get("FSLDIR")
if fsldir is not None:
    fslsubcmd = os.path.join(fsldir, "bin", "fsl_sub")
    flirtcmd = os.path.join(fsldir, "bin", "flirt")
    applywarpcmd = os.path.join(fsldir, "bin", "applywarp")
    fslexists = True
else:
    fslexists = False


def runcmd(thecmd, cluster=False, readable=False, fake=False, debug=False):
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
            if SYSTYPE == "sge":
                sub_cmd = f"{QSUB} -cwd -q fmriprep.q -N {jobname} -w e -R y {filename}".split()
            elif SYSTYPE == "slurm":
                sub_cmd = f"{SBATCH} {scriptfile}".split()
                subprocess.call(sub_cmd)
        else:
            subprocess.call(thecmd)


def batchruncmd(thecmd, readable=False, fake=False, debug=False):
    jobname = thecmd[0].split("/")[-1]
    scriptfile, thescript = make_runscript(thecmd, jobname)
    if not fake:
        if SYSTYPE == "sge":
            sub_cmd = f"{QSUB} -cwd -q fmriprep.q -N {jobname} -w e -R y {filename}".split()
        elif SYSTYPE == "slurm":
            sub_cmd = f"{SBATCH} {scriptfile}".split()
            subprocess.call(sub_cmd)
    else:
        print(thescript)


def mriconvert(inputfile, outputfile, cluster=False, fake=False, debug=False):
    convcmd = []
    convcmd += ["mri_convert"]
    convcmd += [inputfile]
    convcmd += [outputfile]
    runcmd(convcmd, cluster=cluster, fake=fake, debug=debug)


def n4correct(inputfile, outputdir, cluster=False, fake=False, debug=False):
    thename, theext = tide_io.niftisplitext(inputfile)
    n4cmd = []
    n4cmd += ["N4BiasFieldCorrection"]
    n4cmd += ["-d", "3"]
    n4cmd += ["-i", inputfile]
    n4cmd += ["-o", pjoin(outputdir, thename + "_n4" + theext)]
    runcmd(n4cmd, cluster=cluster, fake=fake, debug=debug)


def antsapply(
    inputname,
    targetname,
    outputroot,
    transforms,
    cluster=False,
    fake=False,
    debug=False,
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
    runcmd(applyxfmcmd, cluster=cluster, fake=fake, debug=debug)


def atlasaverageapply(
    inputfile,
    templatefile,
    outputfileroot,
    regionlist=None,
    label=None,
    nozero=False,
    header=False,
    cluster=False,
    fake=False,
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
    runcmd(atlasavgcmd, cluster=cluster, fake=fake, debug=debug)


def readlist(inputfilename):
    inputlist = []
    with open(inputfilename, "r") as thefile:
        lines = thefile.readlines()
        for line in lines:
            inputlist.append(int(line))
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
    if volumeproc:
        theparts = thefmrifile.split("_")
        thesubj = theparts[0]
        therun = theparts[3]
        pedir = theparts[-1][0:2]
    else:
        splitname = thefmrifile.split("_")
        thesubj = therun
        therun = "_".join(splitname[0:2])
        pedir = splitname[2]
    return absname, thesubj, therun, pedir


def findboldfiles_HCP(theroot, thetype, volumeproc, usefixforglm, inputlistfile=None, debug=False):
    # special options depending on whether using volume or grayordinate files
    if inputlistfile is None:
        if volumeproc:
            searchstring = os.path.join(
                connectomeroot,
                "*",
                "preproc",
                "*",
                "MNINonLinear",
                "Results",
                thetype + "_[LR][LR]",
                thetype + "_[LR][LR].nii.gz",
            )
        else:
            searchstring = os.path.join(outputroot, "*", thetype + "_[LR][LR]_grayordinate.nii.gz")
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
                    connectomeroot,
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
                    outputroot, str(subject), thetype + "_[LR][LR]_grayordinate.nii.gz"
                )
            if debug:
                print("searchstring:", searchstring)
            retlist.append(glob.glob(searchstring))
        return [val for sublist in retlist for val in sublist]


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

    return absname, thefilename, thesubj, theses, therun, pedir, thespace


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
                connectomeroot,
                "*",
                "*_Clean_" + thetype + "_[LR][LR].nii.gz",
            )
        else:
            searchstring = os.path.join(outputroot, "*", thetype + "_[LR][LR]_grayordinate.nii.gz")
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
                    connectomeroot,
                    str(subject),
                    "preproc",
                    "*",
                    "MNINonLinear",
                    "Results",
                    thetype + "_[LR][LR]",
                    thetype + "_[LR][LR].nii.gz",
                )
                searchstring = os.path.join(
                    connectomeroot,
                    str(subject),
                    "*_Clean_" + thetype + "_[LR][LR].nii.gz",
                )
            else:
                searchstring = os.path.join(
                    outputroot, str(subject), thetype + "_[LR][LR]_grayordinate.nii.gz"
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


def _get_parser():
    # get the command line parameters
    parser = argparse.ArgumentParser(
        prog="extractrtvals",
        description="extracts and averages rapidtide values based on the subject's freesurfer parcellation",
        usage="%(prog)s",
    )
    parser.add_argument(
        "--inputlist",
        metavar="INPUTLIST",
        type=str,
        action="store",
        dest="inputlistfile",
        help="read a list of subjects to analyze (default is to run all subjects)",
        default=None,
    )
    parser.add_argument(
        "--fake",
        action="store_false",
        dest="doforreal",
        help="just print out the commands that will be executed rather than running them",
        default=True,
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        dest="debug",
        help="output extra debugging information",
        default=False,
    )
    return parser


##########################################################################################
##########################################################################################
#
# control flow starts here
#

##########################################################################################


def main():
    theoutputdir = "/data/frederic/tofrapidtide/derivatives"
    try:
        args = _get_parser().parse_args()
    except SystemExit:
        _get_parser().print_help()
        raise

    if args.debug:
        print(args)

    # define some globals
    qspec = ""

    allgoodaparcasegs = {}
    for themaptype in ["maxtime", "maxcorr"]:
        theboldfiles = findrtfiles(
            "/data/saslan/R21_OASIS_new",
            themaptype,
            inputlistfile=args.inputlistfile,
            debug=args.debug,
        )
        if args.debug:
            print("boldfiles:")
            print(theboldfiles)

        if not makeadir(theoutputdir):
            print("cannot initialize output root directory, exiting")
            sys.exit(1)

        for thefile in theboldfiles:
            absname, thesubj, thesess, therun = parsertname(thefile, debug=args.debug)
            absname = os.path.abspath(thefile)
            anatdir = f"/home/frederic/OMG/OASIS/derivatives/fmriprep/{thesubj}/anat"
            mridir = f"/home/frederic/OMG/OASIS/derivatives/freesurfer/{thesubj}/mri"
            badaparcasegfile = f"{anatdir}/{thesubj}_desc-aparcaseg_dseg.nii.gz"
            initaparcasegmgz = f"{mridir}/aparc+aseg.mgz"
            initaparcasegnifti = os.path.join(theoutputdir, thesubj, f"aparcaseg_fs.nii.gz")
            goodaparcasegfile = os.path.join(theoutputdir, thesubj, f"aparcaseg.nii.gz")
            goodaparcasegT1file = os.path.join(
                theoutputdir, thesubj, f"{thesubj}_desc-aparcaseg_dseg.nii.gz"
            )
            goodaparcasegMNIfile = os.path.join(
                theoutputdir,
                thesubj,
                f"{thesubj}_space-MNI152NLin6Asym_desc-aparcaseg_dseg.nii.gz",
            )
            if not makeadir(os.path.join(theoutputdir, thesubj)):
                print("cannot initialize subject directory, exiting")
                sys.exit(1)
            try:
                thisaparcaseg = allgoodaparcasegs[thesubj]
            except KeyError:
                allgoodaparcasegs[thesubj] = goodaparcasegMNIfile
                thisaparcaseg = allgoodaparcasegs[thesubj]

                if not os.path.isfile(thisaparcaseg):
                    print(f"{thisaparcaseg} does not exist - creating")

                    # copy the reference target
                    localbadaparcasegfile = os.path.join(
                        theoutputdir, thesubj, "badaparcaseg.nii.gz"
                    )
                    shutil.copy(badaparcasegfile, localbadaparcasegfile)

                    # get the fs to native transform
                    fstonativexfm = f"{anatdir}/{thesubj}_from-fsnative_to-T1w_mode-image_xfm.txt"
                    localfstonativexfm = os.path.join(theoutputdir, thesubj, "fstonativexfm.txt")
                    shutil.copy(fstonativexfm, localfstonativexfm)

                    # first make a correct aparcaseg file
                    mriconvert(initaparcasegmgz, initaparcasegnifti)

                    T1toMNIfile = (
                        f"{anatdir}/{thesubj}_from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5"
                    )
                    localT1toMNIfile = os.path.join(theoutputdir, thesubj, "T1toMNI.h5")
                    if not os.path.isfile(localT1toMNIfile):
                        shutil.copy(T1toMNIfile, localT1toMNIfile)

                    # now map it to the native T1 space
                    antsapply(
                        initaparcasegnifti,
                        localbadaparcasegfile,
                        goodaparcasegT1file,
                        [localfstonativexfm],
                        cluster=True,
                        fake=(not args.doforreal),
                        debug=args.debug,
                        interp="NearestNeighbor",
                    )

                    # now map it to MNI152NLin6Asym space
                    antsapply(
                        initaparcasegnifti,
                        thefile,
                        thisaparcaseg,
                        [localT1toMNIfile, localfstonativexfm],
                        cluster=True,
                        fake=(not args.doforreal),
                        debug=args.debug,
                        interp="NearestNeighbor",
                    )

            """
            # now warp rapiditide maps to T1 space
            MNItoT1file = (
                f"{anatdir}/{thesubj}_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5"
            )
            localMNItoT1file = os.path.join(theoutputdir, thesubj, "MNItoT1.h5")

            if args.debug:
                print(f"{thisaparcaseg=}")
                print(f"{MNItoT1file=}")
            if not makeadir(os.path.join(theoutputdir, thesubj, thesess)):
                print("cannot initialize output subject session directory, exiting")
                sys.exit(1)
    

            warpeddata = os.path.join(
                theoutputdir, thesubj, thesess, f"{thesubj}_{thesess}_{therun}_space-T1_desc-{themaptype}_map.nii.gz"
            )
    
            if not os.path.isfile(warpeddata):
                print(f"{warpeddata} does not exist - creating")
                if not os.path.isfile(localMNItoT1file):
                    shutil.copy(MNItoT1file, localMNItoT1file)
                antsapply(thefile,
                            thisaparcaseg,
                            warpeddata,
                            [localMNItoT1file],
                            cluster=True,
                            fake=(not args.doforreal),
                            debug=args.debug,
                            interp="NearestNeighbor")

            """
            # average over regions
            # summaryfileroot = os.path.join(
            #    theoutputdir, thesubj, thesess, f"{thesubj}_{thesess}_{therun}_space-T1_desc-{themaptype}_freesurferaverage"
            # )
            summaryfileroot = os.path.join(
                theoutputdir,
                thesubj,
                thesess,
                f"{thesubj}_{thesess}_{therun}_space-MNI152NLin6Asym_desc-{themaptype}_freesurferaverage",
            )
            thedatalabel = "_".join([thesubj, thesess, therun, themaptype])
            regionlistfile = "/data/frederic/tofrapidtide/code/regionlist"
            if not os.path.isfile(summaryfileroot + "_regionsummaries.csv"):
                atlasaverageapply(
                    thefile,
                    thisaparcaseg,
                    summaryfileroot,
                    regionlist=regionlistfile,
                    label=thedatalabel,
                    nozero=True,
                    cluster=True,
                    header=True,
                    fake=(not args.doforreal),
                    debug=args.debug,
                )


if __name__ == "__main__":
    main()
