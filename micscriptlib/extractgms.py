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
        default=2,
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


def extractgms_workflow(sourcetype):
    global theoutputdir
    try:
        args = _get_parser().parse_args()
    except SystemExit:
        _get_parser().print_help()
        raise

    if args.debug:
        print(args)

    # file locations
    outputroot = f"/data/frederic/{sourcetype}/derivatives"
    if sourcetype == "cole":
        inputroot = "/data/ckorponay/New_HCP_Cleaned_TomMethod"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["REST1", "REST2"]
    elif sourcetype == "psusleep":
        inputroot = "/data/ckorponay/Sleep/fmriprep"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["rest", "sleep"]
    elif sourcetype == "recig":
        inputroot = "/data/ajanes/REcig/fmri"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["cue1", "cue2", "cue3", "cue4", "cue5", "resting"]
    else:
        inputroot = "/data2/HCP1200"
        theoutputdir = os.path.join(outputroot, "gms")
        thetypes = ["REST1", "REST2"]

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

    # loop over all run types
    for thetype in thetypes:
        # find the bold files
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
                absname, thesubj, therun, pedir = micutil.parseconnectomename(
                    thefile, volumeproc=args.volumeproc
                )
                maskname = os.path.join(therundir, "brainmask_fs.2.nii.gz")
            elif sourcetype == "HCP":
                absname, thesubj, therun, pedir = micutil.parseconnectomename(
                    thefile, volumeproc=args.volumeproc
                )
                outroot = os.path.join(thesubj, thesubj + "_" + thetype + "_" + pedir)
                maskname = os.path.join(therundir, "brainmask_fs.2.nii.gz")
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
                thetask = f"task-{thetask}"
                outroot = os.path.join(thesubj, thesubj + "_" + thetask)
                maskname = thefile.replace("desc-smoothAROMAnonaggr_bold", "res-2_desc-brain_mask")
            elif sourcetype == "recig":
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
                thiscommand, ncpus=1, jobname="extractgms"
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

            """# before submitting the job, check to see if output file exists
            startpoint = 20
            endpoint = 1199
            numpoints = endpoint - startpoint + 1
            pointspersection = numpoints // args.numsections
            for section in range(args.numsections):
                sectionname = f"{section + 1}-of-{args.numsections}"
                secstart = startpoint + section * pointspersection
                secend = secstart + pointspersection - 1
                inputrange = f"{secstart} {secend}"
                outputname = f"{os.path.join(theoutputdir, outroot)}_{sectionname}.txt"
                dothis = False
                if args.existcheck:
                    if os.path.isfile(outputname):
                        print(outputname, "exists - skipping")
                        dothis = False
                    else:
                        print(outputname, "does not exist - running")
                        dothis = True
                else:
                    dothis = True
    
                thiscommand = thecommand + [f"--timerange {inputrange}", outputname]
                thecommand.append(outputname)
                scriptfile, thescript = micutil.make_runscript(thiscommand)
                if dothis:
                    if args.doforreal:
                        if SYSTYPE == "sge":
                            sub_cmd = f"{QSUB} -cwd -q fmriprep.q -N {args.jobname} -pe fmriprep 1 -w e -R y {filename}".split()
                        elif SYSTYPE == "slurm":
                            sub_cmd = f"{SBATCH} {scriptfile}".split()
                            subprocess.call(sub_cmd)
                    else:
                        print(thescript)
              """


if __name__ == "__main__":
    extractgms_workflow("psusleep")
