#!/usr/bin/env python
#
#       $Author: frederic $
#       $Date: 2015/03/06 14:12:29 $
#       $Id: preprocess_spo2.py,v 1.1 2015/03/06 14:12:29 frederic Exp $
#
import argparse
import os
import subprocess
import sys

import micscriptlib.util as micutil

DEFAULT_NUMSECTIONS = 1
DEFAULT_TASK = "rest"
DEFAULT_SOURCETYPE = "HCP"
DEFAULT_TIMELIMIT = "00:5:00"
DEFAULT_MEM = "4G"

fslmeantscmd = os.path.join(os.environ["FSLDIR"], "bin", "fslmeants")


def _get_parser():
    # get the command line parameters
    parser = argparse.ArgumentParser(
        prog="extractgms",
        description="extracts global mean timecourses from the connectome data",
        usage="%(prog)s",
    )
    parser.add_argument(
        "--numsections",
        metavar="NUMSECTIONS",
        type=int,
        action="store",
        help="number of sections to split timecourse into (default is 2)",
        default=DEFAULT_NUMSECTIONS,
    )
    parser.add_argument(
        "--timelimit",
        metavar="TIME",
        type=str,
        action="store",
        help=f"time limit of job in HH:MM:SS. Default is {DEFAULT_TIMELIMIT})",
        default=DEFAULT_TIMELIMIT,
    )
    parser.add_argument(
        "--mem",
        metavar="MEM",
        type=str,
        action="store",
        help=f"Memory limit of job in bytes.  Use the suffix M for megabytes, G for gigabytes. Default is {DEFAULT_MEM})",
        default=DEFAULT_MEM,
    )
    parser.add_argument(
        "--sourcetype",
        dest="sourcetype",
        action="store",
        type=str,
        choices=["HCP", "cole", "recig", "psusleep"],
        help=f"Dataset (default is {DEFAULT_SOURCETYPE})",
        default=DEFAULT_SOURCETYPE,
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
        help="extractgms even if output file already exists",
        default=True,
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
        "--dofix",
        action="store_true",
        dest="dofix",
        help="use the fix file instead of minimally preprocessed",
        default=False,
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


def extractgms_workflow():
    global theoutputdir
    try:
        args = _get_parser().parse_args()
    except SystemExit:
        print("Use --help option for detailed information on options.")
        raise

    if args.debug:
        print(args)

    # file locations
    outputroot = f"/data/frederic/{args.sourcetype}/derivatives"
    if args.sourcetype == "cole":
        inputroot = "/data/ckorponay/New_HCP_Cleaned_TomMethod"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["REST1", "REST2"]
    elif args.sourcetype == "psusleep":
        inputroot = "/data/ckorponay/Sleep/fmriprep"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["rest", "sleep"]
    elif args.sourcetype == "recig":
        inputroot = "/data/ajanes/REcig/fmri"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["cue1", "cue2", "cue3", "cue4", "cue5", "cue6", "resting"]
    else:
        inputroot = "/data2/HCP1200"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["REST1", "REST2"]

    SYSTYPE, SUBMITTER, SINGULARITY = micutil.getbatchinfo()

    if args.debug:
        print(f"sourcetype is {args.sourcetype}")
        print(f"\t{inputroot=}")
        print(f"\t{outputroot=}")
        print(f"\t{theoutputdir=}")
        print(f"\t{thetypes=}")
        print(f"\t{SYSTYPE=}")
        print(f"\t{SUBMITTER=}")
        print(f"\t{SINGULARITY=}")

    # loop over all run types
    for thetype in thetypes:
        # find the bold files
        if args.sourcetype == "cole":
            theboldfiles = micutil.findboldfiles_cole(
                inputroot,
                thetype,
                True,
                False,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        elif args.sourcetype == "psusleep":
            theboldfiles = micutil.findboldfiles_psusleep(
                thetype,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        elif args.sourcetype == "ds001927":
            theboldfiles = micutil.findboldfiles_ds001927(
                thetype,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        elif args.sourcetype == "recig":
            theboldfiles = micutil.findboldfiles_recig(
                thetype,
                inputlistfile=args.inputlistfile,
                altpath=True,
                debug=args.debug,
            )
        else:
            theboldfiles = micutil.findboldfiles_HCP(
                inputroot,
                thetype,
                True,
                False,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )

        if not micutil.makeadir(theoutputdir):
            print("cannot initialize output root directory, exiting")
            sys.exit(1)

        for thefile in theboldfiles:
            if args.sourcetype == "cole":
                absname, thesubj, therun, pedir= micutil.parsecolename(
                    thefile, volumeproc=True
                )
                maskname = os.path.join(therundir, "brainmask_fs.2.nii.gz")
            elif args.sourcetype == "HCP":
                absname, thesubj, therun, pedir, MNIDir = micutil.parseconnectomename(
                    thefile, volumeproc=True, debug=args.debug
                )
                outroot = os.path.join(thesubj, thesubj + "_" + thetype + "_" + pedir)
                maskname = os.path.join(therundir, "brainmask_fs.2.nii.gz")
            elif args.sourcetype == "psusleep":
                (
                    absname,
                    thefmrifilename,
                    thesubj,
                    thesess,
                    therun,
                    pedir,
                    thetask,
                    thespace,
                ) = micutil.parsebidsname(thefile)
                thesubj = f"sub-{thesubj}"
                thetask = f"task-{thetask}"
                outroot = os.path.join(thesubj, thesubj + "_" + thetask)
                maskname = thefile.replace("desc-smoothAROMAnonaggr_bold", "res-2_desc-brain_mask")
            elif args.sourcetype == "recig":
                absname, thefilename, thesubj, thesess, thetask = micutil.parserecigname(thefile)
                thesubj = f"sub-{thesubj}"
                thesess = f"ses-{thesess}"
                thetask = f"task-{thetask}"
                if args.debug:
                    print(f"{thesubj=}, {thesess=}, {thetask=}")
                outroot = os.path.join(thesubj, "_".join([thesubj, thesess, thetask]))
                maskname = thefile.replace("filtered_func_data", "mask")
            else:
                print("illegal sourcetype")
                sys.exit()
            absname = os.path.abspath(thefile)
            therundir, thefmrifile = os.path.split(absname)
            theresultsdir, therun = os.path.split(therundir)
            if not micutil.makeadir(os.path.join(theoutputdir, thesubj)):
                print("cannot initialize output subject directory, exiting")
                sys.exit(1)

            print("thefmrifile is", thefmrifile)

            thecommand = []
            fmrifile = absname
            fixfile = os.path.join(therundir, thefmrifile[:-7] + "_hp2000_clean.nii.gz").replace(
                "preproc", "fixextended"
            )
            thecommand.append(fslmeantscmd)
            thecommand.append("-i")
            if args.dofix:
                thecommand.append(fixfile)
            else:
                thecommand.append(fmrifile)

            # add the mask line
            thecommand.append("-m")
            thecommand.append(maskname)

            fulloutputname = f"{os.path.join(theoutputdir, outroot)}_GMS.txt"
            dothis = False
            if args.existcheck:
                if os.path.isfile(fulloutputname):
                    print(fulloutputname, "exists - skipping")
                    dothis = False
                else:
                    print(fulloutputname, "does not exist - running")
                    dothis = True
            else:
                dothis = True

            thiscommand = thecommand + ["-o"]
            thiscommand.append(fulloutputname)
            scriptfile, thescript = micutil.make_runscript(
                thiscommand, timelimit=args.timelimit, mem=args.mem, ncpus=1, jobname="extractgms"
            )
            if dothis:
                if args.doforreal:
                    if SYSTYPE == "sge":
                        sub_cmd = f"{SUBMITTER} -cwd -q short.q -N fslmeants -w e -R y {thiscommand}".split()
                    elif SYSTYPE == "slurm":
                        sub_cmd = f"{SUBMITTER} {scriptfile}".split()
                    subprocess.call(sub_cmd)
                else:
                    print(thescript)


if __name__ == "__main__":
    extractgms_workflow("psusleep")
