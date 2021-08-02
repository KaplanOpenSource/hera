#! /usr/bin/env python
import argparse
import json
import sys
import pandas
import glob
import shutil
import numpy
import os
from distutils.dir_util import copy_tree
try:
    from hera import toolkitHome
except ImportError:
    print("Hera not found, cannot create sources")

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

dimensions      {dimensions};

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

    fromTime = args.fromTimestep

    print(f"Using Time step {uts} for Steady state")

    if len(args.args) == 2:
        # Copy
        dest = os.path.abspath(args.args[1])
        try:
            os.makedirs(dest,exist_ok=args.overwrite)
        except FileExistsError:
            raise FileExistsError("The case already exists, use --overwrite ")

        # copy constant, 0 and system.
        print(f"Process general {orig} directory --> {dest}")
        for general in ["constant", "system", "0"]:
            orig_general = os.path.join(orig, general)
            dest_general = os.path.join(dest, general)
            if os.path.exists(dest_general):
                copy_tree(orig_general, dest_general)
            else:
                shutil.copytree(orig_general, dest_general)

        for proc in glob.glob(os.path.join(orig,"processor*")):
            for dest_time in [fromTime, args.toTimestep]:
                print(f"time {dest_time}/Process {proc}")
                orig_proc = os.path.join(proc,str(uts))
                dest_proc = os.path.join(dest,os.path.basename(proc),dest_time)
                shutil.copytree(orig_proc,dest_proc)

                # Now writing the Hmix and ustar.
                with open(os.path.join(dest_proc,"Hmix"),'w') as HmixFile:
                    Hmix = baseOFFile.replace("{fieldName}","Hmix").replace("{dimensions}","[0 1 0 0 0 0 0]")
                    HmixFile.write(Hmix)

                with open(os.path.join(dest_proc,"ustar"),'w') as ustarFile:
                    ustar = baseOFFile.replace("{fieldName}","ustar").replace("{dimensions}","[0 1 -1 0 0 0 0]")
                    ustarFile.write(ustar)


            if args.copyMesh:
                orig_constant = os.path.join(proc,"constant")
                dest_constant = os.path.join(dest,os.path.basename(proc),"constant")
                shutil.copytree(orig_constant, dest_constant)
            else:
                fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
                destination = os.path.join(dest, os.path.basename(proc), "constant", "polyMesh")
                os.makedirs(os.path.dirname(destination), exist_ok=True)
                os.system(f"ln -s {fullpath} {destination}")

            # We should look into it more closly, why parallel case doesn't recognize the time steps of the
            # processors. For now, just create these directories in the main root as well.
            for dest_time in [fromTime, args.toTimestep]:
                os.makedirs(os.path.join(dest,dest_time),exist_ok=True)



    else:
        print("Not implemented yet")

def prepareDispersionCase(args):
    """
        Prepares the dispersion case:

         * Copies a template to the new directory,
         * Sets up references to the root directory.


         Assumes the root directory is parallel.


    :param args:
             - args : a list :
                Location:
                   0 : The name of the dispersion template
                   1 : The name of the new dispersion case.
                   2 : The name of the root case.
    :return:
    """

    if len(args.args) < 3:
        raise ValueError("Must supply template dir, case dir and rootCase dir")

    template = args.args[0]
    case     = os.path.abspath(args.args[1])
    flowCase = os.path.abspath(args.args[2])

    shutil.copytree(template,case)
    for proc in glob.glob(os.path.join(flowCase,"processor*")):
        fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
        destination = os.path.join(case, os.path.basename(proc), "constant", "polyMesh")
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        os.system(f"ln -s {fullpath} {destination}")
        os.system(f"ln -s {os.path.abspath(proc)} {os.path.join(case, os.path.basename(proc))}/rootCase")

        # create the 0 directory in all processors.
        os.makedirs(os.path.join(case, os.path.basename(proc), '0'), exist_ok=True)

    # linking the decpomposePar dict from the root.
    root_decomposePar = os.path.abspath(os.path.join(flowCase,"system","decomposeParDict"))
    decompose_dest    = os.path.abspath(os.path.join(case,"system"))
    os.system(f"ln -s {root_decomposePar} {decompose_dest}")

    # linking the rootCase in the root directory of the dispersion case.
    os.system(f"ln -s {flowCase} {os.path.join(case,'rootCase')}")

    # create the 0 directory in the root.
    os.makedirs(os.path.join(case,'0'),exist_ok=True)


def makeSourceCylinder(args):

    center = args.center
    params = dict(x=center[0],
                  y=center[1],
                  z=center[2],
                  radius = args.radius,
                  height = args.height,
                  nParticles=args.particles
                  )
    case   = os.path.abspath(args.case[0])

    LSMtoolkit = toolkitHome.getToolkit(toolkitHome.OF_LSM,"tmpProject",casePath=case)

    LSMtoolkit.makeSource(type="Cylinder", **params)

def makeEscapedMassFile(args):

    case   = os.path.abspath(args.casePath)
    massFileName = f"{args.patch}Mass" if args.massFileName is None else args.massFileName
    dt = args.dt
    LSMtoolkit = toolkitHome.getToolkit(toolkitHome.OF_LSM,"tmpProject",casePath=case)
    data = LSMtoolkit.analysis.getMassFromLog(logFile=args.logFile,solver=args.solver)
    data = data.loc[data.name==args.patch].loc[data.action==args.action]
    data["diffMass"] = data.mass.diff()
    data = data.fillna(0)
    if dt is None:
        timesteps = data.time
    else:
        timesteps = numpy.arange(data.time.min(),data.time.max(),float(dt))
        times = pandas.DataFrame({"time":timesteps})
        data = data.set_index("time").join(times.set_index("time"),how="outer").reset_index()
        data = data.interpolate()
    newstr = "/*--------------------------------*- C++ -*----------------------------------*\n" \
             "| =========                 |                                                 |\n" \
             "| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n" \
             "|  \    /   O peration     | Version:  dev                                   |\n" \
             "|   \  /    A nd           | Web:      www.OpenFOAM.org                      |\n" \
             "|    \/     M anipulation  |                                                 |\n" \
             "\*---------------------------------------------------------------------------*/\n" \
             "FoamFile\n{    version     2.0;\n    format      ascii;\n    class       scalarField;\n    object      kinematicCloudPositions;\n}\n" \
             f"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n{len(data)}\n(\n"
    for time in timesteps:
        newstr += f"{float(data.loc[data.time==time].diffMass)}\n"
    newstr += ")"
    print("saving in ",os.path.join(case,"constant",massFileName))
    f = open(os.path.join(case,"constant",massFileName), "w")
    f.write(newstr)
    f.close()

if __name__ =="__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_flowForDispersion = subparsers.add_parser('flowFieldDispersionSteadyState', help='load help')
    parser_prepareDisperionCase = subparsers.add_parser('prepareDispersionCase', help='executePipeline help')
    parser_makeSourceCylinder = subparsers.add_parser('makeSourceCylinder', help='executePipeline help')
    parser_makeEscapedMassFile = subparsers.add_parser('makeEscapedMassFile', help='makeEscapedMassFile help')

    ##### Prepare flow
    parser_flowForDispersion.add_argument('args',
                             nargs='*',
                             type=str,
                             help='[Original directory] [newdirectory]')

    parser_flowForDispersion.add_argument('--useTimestep',type=str,dest="useTimestep",default=None)
    parser_flowForDispersion.add_argument('--toTimestep', type=str,dest="toTimestep",required=True, default=None)
    parser_flowForDispersion.add_argument('--fromTimestep', type=str, dest="fromTimestep", default="0")
    parser_flowForDispersion.add_argument('--copyMesh', dest="copyMesh", default=False,action="store_true")
    parser_flowForDispersion.add_argument('--overwrite', dest="overwrite", default=False,action="store_true",help="Overwrite directory if exists")

    parser_flowForDispersion.set_defaults(func=prepareflowForDispersion)


    #### Prepare dispersion
    parser_prepareDisperionCase.add_argument('args',
                             nargs='*',
                             type=str,
                             help='[Template] [newdirectory] [rootCase]')

    parser_prepareDisperionCase.set_defaults(func=prepareDispersionCase)

    ############################ Source
    ########## Cylinder

    parser_makeSourceCylinder.add_argument('case',
                             nargs=1,
                             type=str,
                             help='[dispersion case dir]')
    parser_makeSourceCylinder.add_argument('--radius', dest="radius", required=True,type=float)
    parser_makeSourceCylinder.add_argument('--center',nargs=3, dest="center",type=float, required=True)
    parser_makeSourceCylinder.add_argument('--height', dest="height", required=True,type=float)
    parser_makeSourceCylinder.add_argument('--particles', dest="particles", required=True,type=int)
    parser_makeSourceCylinder.set_defaults(func=makeSourceCylinder)

    # makeEscapedMassFile
    parser_makeEscapedMassFile.add_argument('logFile',type=str)
    parser_makeEscapedMassFile.add_argument('casePath', type=str)
    parser_makeEscapedMassFile.add_argument('patch', type=str)
    parser_makeEscapedMassFile.add_argument('-dt', type=str)
    parser_makeEscapedMassFile.add_argument('-massFileName', type=str)
    parser_makeEscapedMassFile.add_argument('-solver', type=str,default="StochasticLagrangianSolver")
    parser_makeEscapedMassFile.add_argument('-action', type=str, default="escape")
    parser_makeEscapedMassFile.set_defaults(func=makeEscapedMassFile)
    ############################## arg parse
    args = parser.parse_args()
    print(args)
    args.func(args)
