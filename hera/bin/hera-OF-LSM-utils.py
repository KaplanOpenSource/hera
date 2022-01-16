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
from hera import toolkitHome
from hera.simulations.openFoam import ofObjectHome

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
    case     = os.path.abspath(args.dispersionCase)
    flowCase = os.path.abspath(args.baseFlow)

    constantDir = os.path.join(case,"constant")
    systemDir   = os.path.join(case, "system")
    os.makedirs(constantDir,exist_ok=True)
    os.makedirs(systemDir, exist_ok=True)

    for fls in glob.glob(os.path.join(flowCase,"constant","*")):
        os.system(f"ln -s {fls} {constantDir}")

    for fls in glob.glob(os.path.join(flowCase,"system","*")):
        os.system(f"ln -s {fls} {systemDir}")

    for proc in glob.glob(os.path.join(flowCase,"processor*")):
        fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
        destination = os.path.join(case, os.path.basename(proc), "constant", "polyMesh")
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        os.system(f"ln -s {fullpath} {destination}")
        os.system(f"ln -s {os.path.abspath(proc)} {os.path.join(case, os.path.basename(proc))}/rootCase")

        # create the 0 directory in all processors.
        os.makedirs(os.path.join(case, os.path.basename(proc), '0'), exist_ok=True)



    # linking the decpomposePar dict from the root.
    #root_decomposePar = os.path.abspath(os.path.join(flowCase,"system","decomposeParDict"))
    #decompose_dest    = os.path.abspath(os.path.join(case,"system"))
    #os.system(f"ln -s {root_decomposePar} {decompose_dest}")

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
    data = data.loc[data.filterType == args.patch].loc[data.action == args.action]
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
    parser_makeSourceCylinder = subparsers.add_parser('makeSourceCylinder', help='executePipeline help')
    parser_makeEscapedMassFile = subparsers.add_parser('makeEscapedMassFile', help='makeEscapedMassFile help')

    #### Prepare dispersion
    parser_prepareDisperionCase = subparsers.add_parser('prepareDispersionCase', help='executePipeline help')
    parser_prepareDisperionCase.add_argument("dispersionCase",type=str,help="The target flow case")
    parser_prepareDisperionCase.add_argument("baseFlow",type=str,help="the flow field to use")
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
    args.func(args)
