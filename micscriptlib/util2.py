#!/usr/bin/env python
#
#       $Author: frederic $
#       $Date: 2015/03/06 14:12:29 $
#       $Id: preprocess_spo2.py,v 1.1 2015/03/06 14:12:29 frederic Exp $
#
import argparse
import glob
import os
import subprocess
import sys

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


def make_runscript(thecommand, ncpus=8, timelimit="0:30:00"):
    """
    Create a temporary script file we can submit to qsub.
    """

    import tempfile

    pre = ["#!/bin/bash"]

    # pre += ["#SBATCH -p short"]
    pre += ["#SBATCH --job-name=rapidtide"]
    pre += ["#SBATCH --output=rapidtide.%j.out"]
    pre += ["#SBATCH --error=rapidtide.%j.err"]
    pre += [f"#SBATCH --time={timelimit}"]
    pre += [f"#SBATCH --cpus-per-task={ncpus}"]
    pre += ["#SBATCH --qos=normal"]
    pre += ["#SBATCH --mem=16G"]
    pre += ["#SBATCH --mail-type=FAIL # notifications for job fail"]
    pre += ["#SBATCH --mail-user=bbfrederick@mclean.harvard.edu"]

    script = "\n".join(pre) + "\n" + " \\\n    ".join(thecommand) + "\n"

    _, filename = tempfile.mkstemp()
    with open(filename, "w") as fp:
        fp.write(script)
    return filename, script


# file locations
sourcetype = "HCP"

outputroot = "/data/frederic/connectome"
if sourcetype == "cole":
    connectomeroot = "/data/ckorponay/New_HCP_Cleaned_TomMethod"
    theoutputdir = os.path.join(outputroot, "splitanalysis", "outputtom")
else:
    connectomeroot = "/data2/HCP1200"
    theoutputdir = os.path.join(outputroot, "splitanalysis", "output")

rapidtidecmd = "/cm/shared/anaconda3/envs/mic/bin/rapidtide"


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


def _get_parser():
    # get the command line parameters
    parser = argparse.ArgumentParser(
        prog="splitrapidtide",
        description="runs rapidtide on the connectome data",
        usage="%(prog)s",
    )
    parser.add_argument(
        "--numsections",
        metavar="NUMSECTIONS",
        type=int,
        action="store",
        help="number of sections to split timecourse into (default is 2)",
        default=2,
    )
    parser.add_argument(
        "--startpoint",
        metavar="STARTPOINT",
        type=int,
        action="store",
        help="first TR to use (default is 20)",
        default=20,
    )
    parser.add_argument(
        "--ncpus",
        metavar="NCPUS",
        type=int,
        action="store",
        help="number of cpus per job (default is 8)",
        default=8,
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
        "--noexistcheck",
        action="store_false",
        dest="existcheck",
        help="run rapidtide even if output file already exists",
        default=True,
    )
    parser.add_argument(
        "--swap",
        action="store_true",
        dest="isswapped",
        help="do a null correlation run (analyze with the other PE direction regressor)",
        default=False,
    )
    parser.add_argument(
        "--fake",
        action="store_false",
        dest="doforreal",
        help="just print out the commands that will be executed rather than running them",
        default=True,
    )
    parser.add_argument(
        "--surface",
        action="store_false",
        dest="volumeproc",
        help="process grayordinate file rather than volume file",
        default=True,
    )
    parser.add_argument(
        "--donotusefixforglm",
        action="store_false",
        dest="usefixforglm",
        help="regress lagged timecourses out of original rather than fix cleaned data",
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
    global theoutputdir
    try:
        args = _get_parser().parse_args()
    except SystemExit:
        _get_parser().print_help()
        raise

    if args.debug:
        print(args)

    # define some globals
    extrasuffix = "_spf-2p0"
    spatialfiltwidth = 2.0

    # make the appropriate output directory
    theoutputdir = f"{theoutputdir}_{args.numsections}"

    #    "--nolimitoutput",
    rapidtideopts = [
        "--despecklepasses 4",
        "--filterfreqs 0.009 0.15",
        "--searchrange -7.5 15.0",
        "--lagminthresh 0.25",
        "--lagmaxthresh 2.5",
        f"--nprocs {args.ncpus - 1}",
        "--pickleft",
        "--nofitfilt",
        "--similaritymetric hybrid",
        "--peakfittype gauss",
        "--noprogressbar",
        "--ampthresh 0.15",
    ]

    # set options for volume vs surface
    if args.volumeproc:
        print("setting up for volume processing")
        rapidtideopts.append(f"--spatialfilt {spatialfiltwidth}")
        outputnamesuffix = ""
        qspec = ""
    else:
        print("setting up for grayordinate processing")
        rapidtideopts.append("-c")
        outputnamesuffix = "_grayordinate"
        qspec = ""

    rest1runs = ["rfMRI_REST1"]
    rest2runs = ["rfMRI_REST2"]
    restruns = rest1runs + rest2runs
    emotionruns = ["tfMRI_EMOTION"]
    gamblingruns = ["tfMRI_GAMBLING"]
    languageruns = ["tfMRI_LANGUAGE"]
    motorruns = ["tfMRI_MOTOR"]
    relationalruns = ["tfMRI_RELATIONAL"]
    socialruns = ["tfMRI_SOCIAL"]
    wmruns = ["tfMRI_WM"]

    thetasktypes = (
        emotionruns
        + gamblingruns
        + languageruns
        + motorruns
        + relationalruns
        + socialruns
        + wmruns
    )

    # only run resting state 1 for the time being
    thetypes = rest1runs

    # loop over all run types
    for thetype in thetypes:
        # special options depending on whether using volume or grayordinate files
        if sourcetype == "cole":
            theboldfiles = findboldfiles_cole(
                connectomeroot,
                thetype,
                args.volumeproc,
                args.usefixforglm,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        else:
            theboldfiles = findboldfiles_HCP(
                connectomeroot,
                thetype,
                args.volumeproc,
                args.usefixforglm,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )

        if not makeadir(theoutputdir):
            print("cannot initialize output root directory, exiting")
            sys.exit(1)

        for thefile in theboldfiles:
            if sourcetype == "cole":
                absname, thesubj, therun, pedir = parsecolename(
                    thefile, volumeproc=args.volumeproc
                )
                print(f"{absname=}, {thesubj=}, {therun=}, {pedir=}")
            else:
                absname, thesubj, therun, pedir = parseconnectomename(
                    thefile, volumeproc=args.volumeproc
                )
            absname = os.path.abspath(thefile)
            therundir, thefmrifile = os.path.split(absname)
            theresultsdir, therun = os.path.split(therundir)
            if not makeadir(os.path.join(theoutputdir, thesubj)):
                print("cannot initialize output subject directory, exiting")
                sys.exit(1)

            print("thefmrifile is", thefmrifile)
            if args.isswapped:
                print("doing swapped run")
                swapname = "_swap_"
                if pedir == "LR":
                    regressortouse = "RL"
                else:
                    regressortouse = "LR"
            else:
                print("doing unswapped run")
                swapname = "_noswap_"
                regressortouse = pedir
            if args.volumeproc:
                if args.isswapped:
                    regressoropts = [
                        "--regressor "
                        + os.path.join(
                            theoutputdir,
                            thesubj,
                            thesubj
                            + "_"
                            + thetype
                            + "_noswap_"
                            + regressortouse
                            + outputnamesuffix
                            + extrasuffix
                            + "_reference_fmrires_pass3.txt",
                        )
                    ]
                else:
                    regressoropts = ["--passes 3"]
            else:
                regressoropts = [
                    "--regressor "
                    + os.path.join(
                        theoutputdir,
                        thesubj,
                        thesubj
                        + "_"
                        + thetype
                        + "_noswap_"
                        + regressortouse
                        + outputnamesuffix
                        + extrasuffix
                        + "_reference_fmrires_pass3.txt",
                    )
                ]

            outroot = os.path.join(
                thesubj,
                thesubj + "_" + thetype + swapname + pedir + outputnamesuffix + extrasuffix,
            )
            thecommand = []
            fmrifile = absname
            thecommand.append(rapidtidecmd)
            thecommand.append(fmrifile)
            if args.usefixforglm and (thetype.find("REST") >= 0) and (sourcetype != "cole"):
                cleanspec = "_hp2000_clean"
                # glmname = os.path.join(therundir, thefmrifile[:-7] + '_hp2000_clean.nii.gz')
                glmname = os.path.join(
                    therundir, thefmrifile[:-7] + "_hp2000_clean.nii.gz"
                ).replace("preproc", "fixextended")
                thecommand.append("--glmsourcefile=" + glmname)

            # put in the rapidtide options
            # for theopt in rapidtideopts:
            # splitopt = theopt.split()
            # thecommand.extend(splitopt)
            thecommand += rapidtideopts

            # put in the regressor options
            # thecommand.extend(regressoropts.split())
            thecommand += regressoropts

            # oversample to ~1.0s TR
            osfac = 1
            thetr = 0.72
            while thetr / osfac > 0.7:
                osfac *= 2
            if osfac > 2:
                thecommand.extend(["-O", str(osfac)])

            # before submitting the job, check to see if output file exists
            endpoint = 1199
            numpoints = endpoint - args.startpoint + 1
            pointspersection = numpoints // args.numsections
            print(
                f"dividing timecourse into {args.numsections} sections of {pointspersection} points"
            )
            for section in range(args.numsections):
                sectionname = f"{section + 1}-of-{args.numsections}"
                secstart = args.startpoint + section * pointspersection
                secend = secstart + pointspersection - 1
                inputrange = f"{secstart} {secend}"
                outputname = f"{os.path.join(theoutputdir, outroot)}_{sectionname}"
                dothis = False
                if args.existcheck:
                    if os.path.isfile(outputname + "_DONE.txt"):
                        print(outputname + "_DONE.txt", "exists - skipping")
                        dothis = False
                    else:
                        print(outputname + "_DONE.txt", "does not exist - running")
                        dothis = True
                else:
                    dothis = True

                thiscommand = thecommand + [f"--timerange {inputrange}", outputname]
                scriptfile, thescript = make_runscript(thiscommand)
                if dothis:
                    if args.doforreal:
                        if SYSTYPE == "sge":
                            sub_cmd = f"{QSUB} -cwd -q fmriprep.q -N {args.jobname} -pe fmriprep {args.ncpus} -w e -R y {filename}".split()
                        elif SYSTYPE == "slurm":
                            sub_cmd = f"{SBATCH} {scriptfile}".split()
                            subprocess.call(sub_cmd)
                    else:
                        print(thescript)


if __name__ == "__main__":
    main()
