#!/usr/bin/env python
#
#       $Author: frederic $
#       $Date: 2015/03/06 14:12:29 $
#       $Id: preprocess_spo2.py,v 1.1 2015/03/06 14:12:29 frederic Exp $
#
from __future__ import print_function

import argparse
import os
import subprocess
import sys

import micscriptlib.util as micutil


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
        help="number of sections to split timecourse into (default is 4)",
        default=4,
    )
    parser.add_argument(
        "--task",
        metavar="TASK",
        type=str,
        action="store",
        help="task to process (default is rest)",
        default="rest",
    )
    parser.add_argument(
        "--startpoint",
        metavar="STARTPOINT",
        type=int,
        action="store",
        help="first TR to use (default is 20)",
        default=12,
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


def splitrapidtide_workflow(sourcetype):
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
    outputroot = f"/data/frederic/{sourcetype}/derivatives"
    derivativetype = "rapidtide"
    if sourcetype == "cole":
        inputroot = "/data/ckorponay/New_HCP_Cleaned_TomMethod"
        theoutputdir = os.path.join(outputroot, "rapidtide")
        thetypes = ["REST1", "REST2"]
    elif sourcetype == "psusleep":
        inputroot = "/data/ckorponay/Sleep/fmriprep"
        theoutputdir = os.path.join(outputroot, "rapidtide")
        thetypes = ["rest", "sleep"]
    elif sourcetype == "recig":
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
        print(f"sourcetype is {sourcetype}")
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
        f"--nprocs {args.ncpus - 1}",
        "--pickleft",
        "--nofitfilt",
        "--similaritymetric hybrid",
        "--peakfittype gauss",
        "--noprogressbar",
        "--ampthresh 0.15",
        "--nolimitoutput",
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

    # loop over all run types
    for thetype in thetypes:
        # special options depending on whether using volume or grayordinate files
        print(f"sourcetype is {sourcetype}")
        if sourcetype == "cole":
            theboldfiles = micutil.findboldfiles_cole(
                inputroot,
                thetype,
                args.volumeproc,
                args.usefixforglm,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        elif sourcetype == "psusleep":
            theboldfiles = micutil.findboldfiles_psusleep(
                thetype,
                inputlistfile=args.inputlistfile,
                debug=args.debug,
            )
        elif sourcetype == "recig":
            theboldfiles = micutil.findboldfiles_recig(
                thetype,
                inputlistfile=args.inputlistfile,
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
            if sourcetype == "cole":
                absname, thesubj, therun, pedir = micutil.parsecolename(
                    thefile, volumeproc=args.volumeproc
                )
                print(f"{absname=}, {thesubj=}, {therun=}, {pedir=}")
                outroot = os.path.join(
                    thesubj,
                    thesubj + "_" + thetask + "_" + thesess + outputnamesuffix + extrasuffix,
                )
            elif sourcetype == "recig":
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
                outroot = os.path.join(thesubj, thesubj + "_" + thetask + "_" + thesess)
            elif sourcetype == "psusleep":
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
            else:
                absname, thesubj, therun, pedir = micutil.parseconnectomename(
                    thefile, volumeproc=args.volumeproc
                )
                outroot = os.path.join(
                    thesubj,
                    thesubj + "_" + thetask + "_" + thesess + outputnamesuffix + extrasuffix,
                )
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
            thecommand.append(fmrifile)
            if args.usefixforglm and (thetype.find("REST") >= 0) and (sourcetype == "HCP"):
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

            # before submitting the job, check to see if output file exists
            if sourcetype == "psusleep":
                if thetype == "sleep":
                    endpoint = 428
                else:
                    endpoint = 285
            elif sourcetype == "recig":
                if thetype == "resting":
                    endpoint = 499
                else:
                    endpoint = 427
            elif sourcetype == "HCP":
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
                    thiscommand, timelimit="0:60:00", ncpus=args.ncpus
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
    splitrapidtide_workflow("recig")
