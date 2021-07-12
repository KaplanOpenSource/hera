#! /usr/bin/env python
import argparse
import json
import sys
import glob
import shutil
import os


def prepareflowForDispersion(args):
    """
        Prepares the simulation for the dispersion simulation:

        two modes:

        args is a strcuture with the following:
            - args: list of dest, orig.
            - useTimestep: use the timestep as the base
            - toTimestep: the final timestep that is used.
            - copyMesh  : If true, copies the mesh for eac processor and does not create
                          a link



        A. Inplace:
                    - Delete all timesteps except the the requested timestep (useTimestep if none, use the last).
                    - Change the time to the requested time.
                    - Add the hmix/ustar files.
        B. Copy:
                    - Copy the required time step to the requested directory with the right names.
                    - Add the hmix/ustar files.
    """

    baseOFFile = """
    /*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  dev                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      {fieldName};
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      {dimension};

boundaryField
{
  "proc.*" { type processor;} 
}
"""


    if len(args.args) ==0:
        raise ValueError("Must supply original directory")

    orig = args.args[0]

    TS = [os.path.basename(ts) for ts in glob.glob(os.path.join(orig, "processor0", "*")) if
          os.path.basename(ts).replace(".", "").isdigit()]
    TS.sort()

    if args.useTimestep is None:
        # find maximal TS, assume it is parallel:
        uts = TS[-1]
    else:
        # find the closes TS.
        request = float(args.useTimestep)
        uts = TS[min(range(len(TS)), key=lambda i: abs(float(TS[i]) - request))]

    print(f"Using time {uts}")

    dest_time = args.toTimestep

    if len(args.args) == 2:
        # Copy
        dest = args.args[1]
        os.makedirs(args.args[1],exist_ok=True)

        # copy constant, 0 and system.
        print(f"Process general {orig} directory --> {dest}")
        for general in ["constant","system","0"]:
            orig_general = os.path.join(orig,general)
            dest_general = os.path.join(dest,general)
            shutil.copytree(orig_general, dest_general)

        for proc in glob.glob(os.path.join(orig,"processor*")):
            print(f"Process {proc}")
            orig_proc = os.path.join(proc,str(uts))
            dest_proc = os.path.join(dest,os.path.basename(proc),dest_time)
            shutil.copytree(orig_proc,dest_proc)

            if args.copyMesh:
                orig_constant = os.path.join(proc,"constant")
                dest_constant = os.path.join(dest,os.path.basename(proc),"constant")
                shutil.copytree(orig_constant, dest_constant)
            else:
                fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
                destination = os.path.join(dest, os.path.basename(proc), "constant", "polyMesh")
                os.makedirs(os.path.dirname(destination), exist_ok=True)
                os.system(f"ln -s {fullpath} {destination}")

            # Now writing the Hmix and ustar.
            with open(os.path.join(dest_proc,"Hmix"),'w') as HmixFile:
                Hmix = baseOFFile.replace("{fieldName}","Hmix").replace("{dimensions}","[0 1 0 0 0 0 0]")
                HmixFile.write(Hmix)

            with open(os.path.join(dest_proc,"ustar"),'w') as ustarFile:
                ustar = baseOFFile.replace("{fieldName}","ustar").replace("{dimensions}","[0 1 -1 0 0 0 0]")
                ustarFile.write(ustar)


    else:
        print("Not implemented yet")



if __name__ =="__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_flowForDispersion = subparsers.add_parser('flowForDispersion', help='load help')
    #parser_executePipeline = subparsers.add_parser('prepareDisperions', help='executePipeline help')

    ##### Load parameters
    parser_flowForDispersion.add_argument('args',
                             nargs='*',
                             type=str,
                             help='[Original directory] [newdirectory]')

    parser_flowForDispersion.add_argument('--useTimestep',type=str,dest="useTimestep",default=None)
    parser_flowForDispersion.add_argument('--toTimestep', type=str,dest="toTimestep",required=True, default=None)
    parser_flowForDispersion.add_argument('--copyMesh', dest="copyMesh", default=False,action="store_true")

    parser_flowForDispersion.set_defaults(func=prepareflowForDispersion)

    ############################## arg parse
    args = parser.parse_args()
    args.func(args)
