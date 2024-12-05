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
DEFAULT_SPATIALFILTWIDTH = 2.0
DEFAULT_SOURCETYPE = "HCPYA"
DEFAULT_TIMELIMIT = "00:20:00"
DEFAULT_MEM = "25G"


def _get_parser():
    # get the command line parameters
    parser = argparse.ArgumentParser(
        prog="runretroglms",
        description="runs retroglm on existing rapidtide analyses",
        usage="%(prog)s",
    )
    parser.add_argument(
        "--inputdir",
        metavar="inputdir",
        type=str,
        action="store",
        help=f"Location of the existing rapidtide runs",
        default=DEFAULT_TASK,
    )
    parser.add_argument(
        "--outputlevel",
        dest="outputlevel",
        action="store",
        type=str,
        choices=["min", "less", "normal", "more", "max"],
        help=(
            "The level of file output produced.  'min' produces only absolutely essential files, 'less' adds in "
            "the GLM filtered data (rather than just filter efficacy metrics), 'normal' saves what you "
            "would typically want around for interactive data exploration, "
            "'more' adds files that are sometimes useful, and 'max' outputs anything you might possibly want. "
            "Selecting 'max' will produce ~3x your input datafile size as output.  "
            f'Default is "normal".'
        ),
        default="normal",
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
        help="string containing extra arguments to use to invoke retroglm",
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
        choices=["HCPYA", "cole", "recig", "psusleep", "ds001927"],
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
        help="run retroglm even if output file already exists",
        default=True,
    )
    parser.add_argument(
        "--alternateoutputdir",
        metavar="DIR",
        type=str,
        action="store",
        help="Override the default output directory and put output here.",
        default=None,
    )
    parser.add_argument(
        "--extrasuffix",
        metavar="SUFFIX",
        type=str,
        action="store",
        help="Extre suffix to add to the end of the output file root (default is an empty string)",
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


def runretroglms_workflow():
    try:
        args = _get_parser().parse_args()
    except SystemExit:
        print("Use --help option for detailed information on options.")
        raise

    if args.debug:
        print(args)

    # set the subject list, if we're using that
    if args.inputlistfile is None:
        print("using all subjects")
        inputlist = None
    else:
        print("using subject list")
        inputlist = micutil.readlist(args.inputlistfile)

    # file locations
    outputroot = f"/data/frederic/{args.sourcetype}"
    derivativetype = "retroglm"
    if args.sourcetype == "cole":
        inputroot = "/data/ckorponay/New_HCP_Cleaned_TomMethod"
        theoutputdir = os.path.join(outputroot, "derivatives", "rapidtide")
        thetypes = ["REST1", "REST2"]
    elif args.sourcetype == "ds001927":
        thebidsroot = "/data/frederic/ds001927"
        hassessions = True
        inputroot = f"{outputroot}/fmriprep"
        theoutputdir = os.path.join(outputroot, "derivatives", "rapidtide")
        thetypes = ["restpre", "restpost"]
    elif args.sourcetype == "psusleep":
        inputroot = "/data/ckorponay/Sleep/fmriprep"
        theoutputdir = os.path.join(outputroot, "derivatives", "rapidtide")
        thetypes = ["rest", "sleep"]
    elif args.sourcetype == "recig":
        inputroot = "/data/ajanes/REcig/fmri"
        theoutputdir = os.path.join(outputroot, "derivatives", "retroglm")
        if args.task == "rest":
            thetypes = ["resting"]
        else:
            thetypes = ["cue1", "cue2", "cue3", "cue4", "cue5", "cue6", "resting"]
    elif args.sourcetype == "HCPA":
        print("using HCPA")
        print("NOT IMPLEMENTED - exiting")
        raise
    elif args.sourcetype == "ABCD":
        print("using ABCD")
        print("NOT IMPLEMENTED - exiting")
        raise
    else:
        # HCPYA
        inputroot = "/data2/HCP1200"
        procroot = args.inputdir
        theoutputdir = os.path.join(outputroot, "derivatives", "retroglm")
        rest1runs = ["rfMRI_REST1"]
        rest2runs = ["rfMRI_REST2"]
        emotionruns = ["tfMRI_EMOTION"]
        gamblingruns = ["tfMRI_GAMBLING"]
        languageruns = ["tfMRI_LANGUAGE"]
        motorruns = ["tfMRI_MOTOR"]
        relationalruns = ["tfMRI_RELATIONAL"]
        socialruns = ["tfMRI_SOCIAL"]
        wmruns = ["tfMRI_WM"]

        '''thetypes = (
            emotionruns
            + gamblingruns
            + languageruns
            + motorruns
            + relationalruns
            + socialruns
            + wmruns
            + rest1runs
            + rest2runs
        )'''

        thetypes = rest1runs + rest2runs
    retroglmcmd = "/cm/shared/miniforge3/envs/mic/bin/retroglm"
    retroglmcmd = "retroglm"
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
    retroglmopts = [
        f"--nprocs {args.ncpus}",
        "--noprogressbar",
    ]


    # put in the extra retroglm options
    if args.extraargs is not None:
        retroglmopts.append(args.extraargs)

    # set options for volume vs surface
    if args.volumeproc:
        print("setting up for volume processing")
        outputnamesuffix = None
        qspec = ""
    else:
        print("setting up for grayordinate processing")
        retroglmopts.append("-c")
        outputnamesuffix = "grayordinate"
        qspec = ""

    retroglmopts.append(f"--outputlevel {args.outputlevel}")

    # loop over all run types
    for thetype in thetypes:
        # special options depending on whether using volume or grayordinate files
        print(f"sourcetype is {args.sourcetype}")
        if args.debug:
            print("about to find files")
            print(f"\t{inputroot=}")
            print(f"\t{procroot=}")
            print(f"\t{thetype=}")
            print(f"\t{args.inputlistfile=}")
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
        elif args.sourcetype == "ds001927":
            theboldfiles = micutil.findboldfiles_fmriprep(
                subjects=inputlist,
                sessions=None,
                tasks=[thetype],
                hassessions=True,
                bidsroot=thebidsroot,
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
            therapidtideruns = micutil.findrapidtideruns_HCPYA(
                procroot,
                thetype,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )

        for datafileroot in therapidtideruns:
            if args.sourcetype == "cole":
                absname, thesubj, therun, pedir = micutil.parsecolename(
                    datafileroot, volumeproc=args.volumeproc
                )
                print(f"{absname=}, {thesubj=}, {therun=}, {pedir=}")
                outroot = os.path.join(
                    thesubj,
                    thesubj + "_" + thetask + "_" + thesess + outputnamesuffix + args.extrasuffix,
                )
            elif args.sourcetype == "recig":
                absname, thefmrifilename, thesubj, thesess, thetask = micutil.parserecigname(
                    datafileroot
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
            elif args.sourcetype == "psusleep" or args.sourcetype == "ds001927":
                (
                    absname,
                    thefmrifilename,
                    thesubj,
                    thesess,
                    therun,
                    pedir,
                    thetask,
                    thespace,
                ) = micutil.parsebidsname(datafileroot)
                thesubj = f"sub-{thesubj}"
                thetask = f"task-{thetask}"
                if thesess is not None:
                    thesess = f"ses-{thesess}"
                if args.debug:
                    print(f"{absname=}")
                    print(f"{thefmrifilename=}")
                    print(f"{thesubj=}")
                    print(f"{thesess=}")
                    print(f"{therun=}")
                    print(f"{pedir=}")
                    print(f"{thetask=}")
                    print(f"{thespace=}")
                if thesess is not None:
                    pathparts = [thesubj, thesess, thesubj + "_" + thesess + "_" +  thetask]
                else:
                    pathparts = [thesubj, thesubj + "_" +  thetask]
                outroot = os.path.join(*pathparts)
            else:
                # HCPYA
                therootname, thesubj, therun, pedir = micutil.parseconnectomerapidtidename(datafileroot, debug=args.debug)
                nameparts = [thesubj, therun, pedir]
                print(f"{therootname=}, {thesubj=}, {therun=}, {pedir=}")
                if outputnamesuffix is not None:
                    nameparts.append(outputnamesuffix)
                if args.extrasuffix is not None:
                    nameparts.append(args.extrasuffix)
                outroot = os.path.join(
                    thesubj,
                    "_".join(nameparts)
                )
                thesess = None

            fmrifile = micutil.findHCPYAsourcefile(inputroot, thesubj, therun, pedir, usefixforglm=args.usefixforglm, debug=args.debug)
            print(f"{fmrifile=}")

            thecommand = []
            thecommand.append(retroglmcmd)
            thecommand.append(datafileroot)
            thecommand.append(fmrifile)
            thecommand += retroglmopts

            # check to see if we should override the output directory
            if args.alternateoutputdir is not None:
                theresultsdir = os.path.join(args.alternateoutputdir, thesubj)
                if args.debug:
                    print(f"alternate output is {theresultsdir=}")
                thecommand.append(f"--alternateoutput {os.path.join(theresultsdir, therootname)}")
                if not micutil.makeadir(theresultsdir):
                    print("cannot initialize output root directory, exiting")
                    sys.exit(1)


            # before submitting the job, check to see if output file exists
            if args.sourcetype == "psusleep":
                if thetype == "sleep":
                    endpoint = 428
                else:
                    endpoint = 285
                thetr = 2.0
            elif args.sourcetype == "recig":
                if thetype == "resting":
                    endpoint = 499
                else:
                    endpoint = 427
                thetr = 0.8
            elif args.sourcetype == "ds001927":
                endpoint = 599
                thetr = 1.25
            elif args.sourcetype == "HCPYA":
                endpoint = 1199
                thetr = 0.72
            else:
                endpoint = 100
                thetr = 3.0
            simcalcstart = int(round(100 * (0.72 / thetr), 0))
            numpoints = endpoint - args.startpoint + 1
            pointspersection = numpoints
            for section in range(1):
                sectionname = f"{section + 1}-of-{1}"
                secstart = args.startpoint + section * pointspersection
                secend = secstart + pointspersection - 1
                inputrange = f"{secstart} {secend}"
                outputname = f"{os.path.join(theoutputdir, outroot)}"

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

                thiscommand = thecommand
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
    runretroglms_workflow()
