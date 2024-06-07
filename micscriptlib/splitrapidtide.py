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
DEFAULT_STARTPOINT = 12
DEFAULT_NCPUS = 8
DEFAULT_SOURCETYPE = "HCP"
DEFAULT_TIMELIMIT = "00:20:00"
DEFAULT_MEM = "25G"


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
        help=f"number of sections to split timecourse into (default is {DEFAULT_NUMSECTIONS})",
        default=DEFAULT_NUMSECTIONS,
    )
    parser.add_argument(
        "--task",
        metavar="TASK",
        type=str,
        action="store",
        help=f"task to process (default is {DEFAULT_TASK})",
        default=DEFAULT_TASK,
    )
    parser.add_argument(
        "--startpoint",
        metavar="STARTPOINT",
        type=int,
        action="store",
        help=f"first TR to use (default is {DEFAULT_STARTPOINT})",
        default=DEFAULT_STARTPOINT,
    )
    parser.add_argument(
        "--extraargs",
        metavar="ARGS",
        type=str,
        action="store",
        dest="extraargs",
        help="string containing extra arguments to use to invoke rapidtide",
        default=None,
    )
    parser.add_argument(
        "--ncpus",
        metavar="NCPUS",
        type=int,
        action="store",
        help=f"number of cpus per job (default is {DEFAULT_NCPUS})",
        default=DEFAULT_NCPUS,
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
        help="run rapidtide even if output file already exists",
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
        "--donotusefixforglm",
        action="store_false",
        dest="usefixforglm",
        help="regress lagged timecourses out of original rather than fix cleaned data",
        default=True,
    )
    parser.add_argument(
        "--domotion",
        action="store_true",
        dest="domotion",
        help="do motion regression",
        default=False,
    )
    parser.add_argument(
        "--dodesignmat",
        action="store_true",
        dest="dodesignmat",
        help="regress out task design",
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


def splitrapidtide_workflow():
    try:
        args = _get_parser().parse_args()
    except SystemExit:
        _get_parser().print_help()
        raise

    if args.debug:
        print(args)

    # define some globals
    extrasuffix = ""
    spatialfiltwidth = 2.0

    # file locations
    outputroot = f"/data/frederic/{args.sourcetype}/derivatives"
    derivativetype = "rapidtide"
    if args.sourcetype == "cole":
        inputroot = "/data/ckorponay/New_HCP_Cleaned_TomMethod"
        theoutputdir = os.path.join(outputroot, "rapidtide")
        thetypes = ["REST1", "REST2"]
    elif args.sourcetype == "psusleep":
        inputroot = "/data/ckorponay/Sleep/fmriprep"
        theoutputdir = os.path.join(outputroot, "rapidtide")
        thetypes = ["rest", "sleep"]
    elif args.sourcetype == "recig":
        inputroot = "/data/ajanes/REcig/fmri"
        theoutputdir = os.path.join(outputroot, "rapidtide")
        thetypes = ["cue1", "cue2", "cue3", "cue4", "cue5", "resting"]
    else:
        inputroot = "/data2/HCP1200"
        theoutputdir = os.path.join(outputroot, "rapidtide")
        rest1runs = ["rfMRI_REST1"]
        rest2runs = ["rfMRI_REST2"]
        emotionruns = ["tfMRI_EMOTION"]
        gamblingruns = ["tfMRI_GAMBLING"]
        languageruns = ["tfMRI_LANGUAGE"]
        motorruns = ["tfMRI_MOTOR"]
        relationalruns = ["tfMRI_RELATIONAL"]
        socialruns = ["tfMRI_SOCIAL"]
        wmruns = ["tfMRI_WM"]

        thetypes = (
            emotionruns
            + gamblingruns
            + languageruns
            + motorruns
            + relationalruns
            + socialruns
            + wmruns
            + rest1runs
            + rest2runs
        )
    rapidtidecmd = "/cm/shared/anaconda3/envs/mic/bin/rapidtide"
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

    # make the appropriate output directory
    theoutputdir = f"{theoutputdir}_{args.numsections}"

    rapidtideopts = [
        "--despecklepasses 4",
        "--filterfreqs 0.009 0.15",
        "--searchrange -7.5 15.0",
        "--lagminthresh 0.25",
        "--lagmaxthresh 2.5",
        f"--nprocs {args.ncpus}",
        "--pickleft",
        "--nofitfilt",
        "--similaritymetric hybrid",
        "--noprogressbar",
        "--ampthresh 0.15",
        "--outputlevel normal",
    ]

    # put in the extra rapidtide options
    if args.extraargs is not None:
        rapidtideopts.append(args.extraargs)

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

    # loop over all run types
    for thetype in thetypes:
        # special options depending on whether using volume or grayordinate files
        print(f"sourcetype is {args.sourcetype}")
        if args.sourcetype == "cole":
            theboldfiles = micutil.findboldfiles_cole(
                inputroot,
                thetype,
                args.volumeproc,
                args.usefixforglm,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        elif args.sourcetype == "psusleep":
            theboldfiles = micutil.findboldfiles_psusleep(
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
                args.volumeproc,
                args.usefixforglm,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )

        if not micutil.makeadir(theoutputdir):
            print("cannot initialize output root directory, exiting")
            sys.exit(1)

        for thefile in theboldfiles:
            if args.sourcetype == "cole":
                absname, thesubj, therun, pedir = micutil.parsecolename(
                    thefile, volumeproc=args.volumeproc
                )
                print(f"{absname=}, {thesubj=}, {therun=}, {pedir=}")
                outroot = os.path.join(
                    thesubj,
                    thesubj + "_" + thetask + "_" + thesess + outputnamesuffix + extrasuffix,
                )
                motionfile = None
                designfile = None
            elif args.sourcetype == "recig":
                absname, thefmrifilename, thesubj, thesess, thetask = micutil.parserecigname(
                    thefile
                )
                thesubj = f"sub-{thesubj}"
                thetask = f"task-{thetask}"
                thesess = f"ses-{thesess}"
                if args.debug:
                    print(f"{absname=}")
                    print(f"{thefmrifilename=}")
                    print(f"{thesubj=}")
                    print(f"{thesess=}")
                    print(f"{thetask=}")
                outroot = os.path.join(thesubj, thesubj + "_" + thesess + "_" + thetask)
                motionfile = thefile.replace("filtered_func_data.nii.gz", "mc/prefiltered_func_data_mcf.par")
                designfile = thefile.replace("filtered_func_data.nii.gz", "design.mat:col_00,col_01,col_02")
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
                thetask = f"sub-{thetask}"
                if args.debug:
                    print(f"{absname=}")
                    print(f"{thefmrifilename=}")
                    print(f"{thesubj=}")
                    print(f"{thesess=}")
                    print(f"{therun=}")
                    print(f"{pedir=}")
                    print(f"{thetask=}")
                    print(f"{thespace=}")
                outroot = os.path.join(
                    thesubj,
                    thesubj + "_" + thetask + "_" + thesess + outputnamesuffix + extrasuffix,
                )
                motionfile = None
                designfile = None
            else:
                absname, thesubj, therun, pedir = micutil.parseconnectomename(
                    thefile, volumeproc=args.volumeproc
                )
                outroot = os.path.join(
                    thesubj,
                    thesubj + "_" + thetask + "_" + thesess + outputnamesuffix + extrasuffix,
                )
                motiondir, thefmrifile = os.path.split(thefile)
                motionfile = os.path.join(motiondir, "Movement_Regressors.txt:0-5")
                designfile = None
            absname = os.path.abspath(thefile)
            therundir, thefmrifile = os.path.split(absname)
            theresultsdir, therun = os.path.split(therundir)
            if not micutil.makeadir(os.path.join(theoutputdir, thesubj)):
                print("cannot initialize output subject directory, exiting")
                sys.exit(1)

            print("thefmrifile is", thefmrifile)

            thecommand = []
            fmrifile = absname
            thecommand.append(rapidtidecmd)
            thecommand += rapidtideopts
            thecommand.append(fmrifile)
            if args.domotion and motionfile is not None:
                if args.sourcetype == "recig":
                    if thetype == "resting":
                        thecommand.append(f"--motionfile {motionfile}")
                else:
                    thecommand.append(f"--motionfile {motionfile}")
            if args.dodesignmat and (designfile is not None) and (thetype != "resting"):
                thecommand.append(f"--confoundfile {designfile}")
            if args.usefixforglm and (thetype.find("REST") >= 0) and (args.sourcetype == "HCP"):
                cleanspec = "_hp2000_clean"
                # glmname = os.path.join(therundir, thefmrifile[:-7] + '_hp2000_clean.nii.gz')
                glmname = os.path.join(
                    therundir, thefmrifile[:-7] + "_hp2000_clean.nii.gz"
                ).replace("preproc", "fixextended")
                thecommand.append(f"--glmsourcefile {glmname}")


            # before submitting the job, check to see if output file exists
            if args.sourcetype == "psusleep":
                if thetype == "sleep":
                    endpoint = 428
                else:
                    endpoint = 285
            elif args.sourcetype == "recig":
                if thetype == "resting":
                    endpoint = 499
                else:
                    endpoint = 427
            elif args.sourcetype == "HCP":
                endpoint = 1199
            else:
                endpoint = 100
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
                if args.numsections == 1:
                    outputname = f"{os.path.join(theoutputdir, outroot)}"
                else:
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
                scriptfile, thescript = micutil.make_runscript(
                    thiscommand, timelimit=args.timelimit, mem=args.mem, ncpus=args.ncpus, debug=args.debug
                )
                if dothis:
                    if args.doforreal:
                        if SYSTYPE == "sge":
                            sub_cmd = f"{SUBMITTER} -cwd -q fmriprep.q -N {args.jobname} -pe fmriprep {args.ncpus} -w e -R y {thiscommand}".split()
                        elif SYSTYPE == "slurm":
                            sub_cmd = f"{SUBMITTER} {scriptfile}".split()
                        subprocess.call(sub_cmd)
                    else:
                        print(thescript)


if __name__ == "__main__":
    splitrapidtide_workflow()
