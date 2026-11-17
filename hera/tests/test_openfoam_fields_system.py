"""OpenFOAM 10 runs from the fields Hera writes, and Hera reads what it writes back (#1169).

The cavity tutorial (icoFoam) runs three times:

* reference -- the tutorial as shipped, with its own 0/U and 0/p;
* hera      -- the same mesh, with 0/U and 0/p written by Hera
               (OFObjectHome.getEmptyFieldFromCase, addBoundaryField,
               writeToCase -- foamlib underneath);
* parallel  -- the same, decomposed into two processors, with the
               processor fields written by Hera, run with mpirun and
               reconstructed.

Hera's fields describe the same initial and boundary conditions as the
tutorial's, so every run must end with the reference velocity and pressure.
The results are read back through Hera's readFieldAsDataFrame.

OpenFOAM comes from PATH (blockMesh) or, when HERA_OPENFOAM_DOCKER_IMAGE
names an image, from Docker (CI uses openfoam/openfoam10-paraview56).
Without either the module is skipped.
"""
import os
import shutil
import subprocess

import numpy
import pytest

from hera.simulations.openFoam import FLOWTYPE_INCOMPRESSIBLE
from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome
from hera.tests.conftest import compare_outputs

DOCKER_IMAGE = os.environ.get("HERA_OPENFOAM_DOCKER_IMAGE")
NATIVE = shutil.which("blockMesh") is not None

pytestmark = [
    pytest.mark.openfoam,
    pytest.mark.skipif(not (NATIVE or DOCKER_IMAGE),
                       reason="needs OpenFOAM 10: blockMesh on PATH or HERA_OPENFOAM_DOCKER_IMAGE"),
]

TUTORIAL = "$FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity"
END_TIME = "0.5"
DECOMPOSE_PAR_DICT = """FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); }
"""


def _foam(case, command):
    """Run a shell command with the OpenFOAM 10 environment, inside ``case``."""
    script = f"source /opt/openfoam10/etc/bashrc >/dev/null 2>&1; set -e; {command}"
    if NATIVE:
        argv = ["bash", "-c", script]
    else:
        argv = ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
                "-v", f"{case}:/case", "-w", "/case", "--entrypoint", "bash", DOCKER_IMAGE, "-c", script]
    result = subprocess.run(argv, cwd=case, capture_output=True, text=True)
    assert result.returncode == 0, f"{command} failed:\n{result.stdout[-4000:]}\n{result.stderr[-4000:]}"


def _writeHeraFields(case, readParallel):
    home = OFObjectHome()
    velocity = home.getEmptyFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, str(case), readParallel=readParallel)
    velocity.addBoundaryField("movingWall", type="fixedValue", value="uniform (1 0 0)")
    velocity.addBoundaryField("fixedWalls", type="noSlip")
    velocity.writeToCase(str(case), 0)
    if readParallel:
        # The regex route: '".*"' zeroGradient and '"proc.*"' processor.  OpenFOAM
        # matches a patch name against them only if they are written quoted.
        pressure = home.getEmptyField("p", FLOWTYPE_INCOMPRESSIBLE, noOfProc=2, addParallelProc=True)
        pressure.addBoundaryField("frontAndBack", type="empty")
    else:
        pressure = home.getEmptyFieldFromCase("p", FLOWTYPE_INCOMPRESSIBLE, str(case), readParallel=False)
    pressure.writeToCase(str(case), 0)


def _result(case, fieldName, readParallel=False):
    df = OFObjectHome().readFieldAsDataFrame(fieldName, str(case), times=END_TIME, readParallel=readParallel)
    return df.drop(columns="time").reset_index(drop=True)


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    case = tmp_path_factory.mktemp("reference")
    _foam(case, f"cp -r {TUTORIAL}/. . && blockMesh && icoFoam")
    return case


@pytest.fixture(scope="module")
def heraCase(tmp_path_factory):
    case = tmp_path_factory.mktemp("hera")
    _foam(case, f"cp -r {TUTORIAL}/. . && rm 0/U 0/p && blockMesh")
    _writeHeraFields(case, readParallel=False)
    _foam(case, "icoFoam")
    return case


@pytest.fixture(scope="module")
def parallelCase(tmp_path_factory):
    case = tmp_path_factory.mktemp("parallel")
    _foam(case, f"cp -r {TUTORIAL}/. . && rm 0/U 0/p && blockMesh")
    (case / "system" / "decomposeParDict").write_text(DECOMPOSE_PAR_DICT)
    _foam(case, "decomposePar")
    _writeHeraFields(case, readParallel=True)
    _foam(case, "mpirun --oversubscribe -np 2 icoFoam -parallel && reconstructPar -latestTime")
    return case


@pytest.mark.parametrize("fieldName", ["U", "p"])
def test_a_run_from_the_fields_hera_wrote_ends_with_the_reference_result(reference, heraCase, fieldName):
    assert compare_outputs(_result(heraCase, fieldName), _result(reference, fieldName), "dataframe")


@pytest.mark.parametrize("fieldName, columns", [("U", ["Ux", "Uy", "Uz"]), ("p", ["p"])])
def test_a_decomposed_run_from_the_fields_hera_wrote_ends_with_the_reference_result(
        reference, parallelCase, fieldName, columns):
    """Reconstructed, cell for cell.  The parallel linear solvers converge to the
    same tolerance by a different path, so the match is to that tolerance."""
    result, expected = _result(parallelCase, fieldName), _result(reference, fieldName)
    assert list(result["region"]) == list(expected["region"])
    assert numpy.allclose(result[columns], expected[columns], atol=1e-5, equal_nan=True)


def test_hera_reads_the_decomposed_result_cell_for_cell(parallelCase):
    perProcessor = _result(parallelCase, "U", readParallel=True)
    reconstructed = _result(parallelCase, "U")
    internal = perProcessor[perProcessor["region"] == "internalField"]
    assert sorted(internal["processor"].unique()) == [0, 1]
    assert len(internal) == (reconstructed["region"] == "internalField").sum() == 400
    for column in ("Ux", "Uy", "Uz"):
        assert numpy.allclose(numpy.sort(internal[column]),
                              numpy.sort(reconstructed.loc[reconstructed["region"] == "internalField", column]))
