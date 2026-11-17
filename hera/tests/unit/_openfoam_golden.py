"""OpenFOAM field I/O scenarios, frozen against the PyFoam implementation.

Issue #1169 replaced PyFoam with foamlib in ``preprocessOFObjects``.  Every
scenario here goes through the public hera API only (``OFObjectHome``,
``OFField``, ``extractFieldFile``), so the same code runs on both
implementations.  ``hera/tests/expected/openfoam_fields/pyfoam/`` holds what
the last PyFoam implementation (commit 9eb888e6, PyFoam 2023.7) produced for
each scenario: the files it wrote, the DataFrames it returned as CSV, or the
exception it raised (``<name>.error``).  ``test_openfoam_field_golden.py``
runs the scenarios on the current code and compares the values.

The input case under ``openfoam_fields/input_case`` is a 4-cell, 5-patch
mesh plus the same case decomposed into two processors.  It deliberately
contains a scalar patch with exactly 3 faces (foamlib writes any 3-long
array as a uniform vector), patches with no ``value``, a patch with an empty
value list, and processor patches.

The goldens can only be regenerated on a checkout that still has the PyFoam
implementation, with a real PyFoam on the path::

    git checkout 9eb888e6 -- hera/simulations/openFoam/preprocessOFObjects
    PYTHONPATH=.:Hermes:<pyargos>:<dir with PyFoam 2023.7> \\
        python -m hera.tests.unit._openfoam_golden
"""
import io
import os
import shutil

import pandas

from hera.simulations.openFoam import FIELDCOMPUTATION_EULERIAN, FIELDTYPE_TENSOR, FLOWTYPE_INCOMPRESSIBLE
from hera.simulations.openFoam.preprocessOFObjects.OFField import OFField
from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome
from hera.simulations.openFoam.preprocessOFObjects.utils import extractFieldFile

GOLDEN_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "expected", "openfoam_fields"))
INPUT_CASE = os.path.join(GOLDEN_ROOT, "input_case")
PYFOAM_DIR = os.path.join(GOLDEN_ROOT, "pyfoam")

WRITE_SCENARIOS = {}  # name -> fn(home, outputCase): writes files under outputCase
FRAME_SCENARIOS = {}  # name -> fn(home) -> pandas.DataFrame

_U = ["Ux", "Uy", "Uz"]


def _write(fn):
    WRITE_SCENARIOS[fn.__name__] = fn
    return fn


def _frame(fn):
    FRAME_SCENARIOS[fn.__name__] = fn
    return fn


def _input(*parts):
    return os.path.join(INPUT_CASE, *parts)


def _scalarFrame():
    """p on 4 cells, 2 inlet faces and 3 outlet faces (the ambiguous length)."""
    return pandas.DataFrame({
        "region": ["internalField"] * 4 + ["boundaryField"] * 5,
        "boundary": [None] * 4 + ["inlet"] * 2 + ["outlet"] * 3,
        "processorIndex": [0, 1, 2, 3, 0, 1, 0, 1, 2],
        "p": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
    })


def _parallelVectorFrame():
    """U on two processors: 2 cells each, one inlet face on 0, 3 outlet faces on 1."""
    return pandas.DataFrame({
        "processor": [0, 0, 0, 1, 1, 1, 1, 1],
        "region": ["internalField", "internalField", "boundaryField",
                   "internalField", "internalField", "boundaryField", "boundaryField", "boundaryField"],
        "boundary": [None, None, "inlet", None, None, "outlet", "outlet", "outlet"],
        "processorIndex": [0, 1, 0, 0, 1, 0, 1, 2],
        "Ux": [0.1, 0.2, 1.0, 0.3, 0.4, 0.5, 0.6, 0.7],
        "Uy": [0.0, 0.01, 0.0, 0.0, 0.02, 0.0, 0.0, 0.1],
        "Uz": [0.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0, 0.0],
    })


# --------------------------------------------------------------------------- #
#  Scenarios that write files                                                  #
# --------------------------------------------------------------------------- #
@_write
def empty_scalar(home, out):
    home.getEmptyField("p", FLOWTYPE_INCOMPRESSIBLE).writeToCase(out, 0)


@_write
def empty_vector_parallel(home, out):
    home.getEmptyField("U", FLOWTYPE_INCOMPRESSIBLE, noOfProc=2, addParallelProc=True).writeToCase(out, 0)


@_write
def empty_tensor(home, out):
    field = OFField(name="R", fileName="R", fieldType=FIELDTYPE_TENSOR,
                    fieldComputation=FIELDCOMPUTATION_EULERIAN, dimensions={"m": 2, "s": -2})
    field.writeToCase(out, 0)


@_write
def documented_api(home, out):
    """The sequence shown in docs/toolkits/simulations/openfoam.md."""
    field = home.getEmptyField("U", FLOWTYPE_INCOMPRESSIBLE)
    field.setInternalUniformFieldValue([5, 0, 0])
    field.addBoundaryField("inlet", type="fixedValue", value="uniform (5 0 0)")
    field.addBoundaryField("outlet", type="zeroGradient")
    field.addBoundaryField("ground", type="noSlip")
    field.addProcBoundary()
    field.writeToCase(out, 0)


@_write
def scalar_from_dataframe(home, out):
    field = home.getEmptyField("p", FLOWTYPE_INCOMPRESSIBLE)
    field.setFieldFromDataFrame(_scalarFrame())
    field.writeToCase(out, 0)


@_write
def vector_from_dataframe_parallel(home, out):
    field = home.getEmptyField("U", FLOWTYPE_INCOMPRESSIBLE)
    field.setFieldFromDataFrame(_parallelVectorFrame())
    field.writeToCase(out, 0)


@_write
def boundaries_from_case_single(home, out):
    home.getEmptyFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=False).writeToCase(out, 0)


@_write
def boundaries_from_case_single_internal_value(home, out):
    field = home.getEmptyFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE,
                                       internalValue=[5, 0, 0], readParallel=False)
    field.writeToCase(out, 0)


@_write
def boundaries_from_case_parallel(home, out):
    home.getEmptyFieldFromCase("p", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=True).writeToCase(out, 0)


@_write
def roundtrip_single_U(home, out):
    home.readFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=False).writeToCase(out, 0)


@_write
def roundtrip_single_p(home, out):
    home.readFieldFromCase("p", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=False).writeToCase(out, 0)


@_write
def roundtrip_parallel_U(home, out):
    home.readFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=True).writeToCase(out, 0)


@_write
def roundtrip_parallel_p(home, out):
    home.readFieldFromCase("p", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=True).writeToCase(out, 0)


# --------------------------------------------------------------------------- #
#  Scenarios that return a DataFrame                                           #
# --------------------------------------------------------------------------- #
@_frame
def extract_U(home):
    return extractFieldFile(_input("0", "U"), columnNames=_U)


@_frame
def extract_U_patch_list(home):
    return extractFieldFile(_input("0", "U"), columnNames=_U, patchNameList=["inlet", "walls"])


@_frame
def extract_processor0_U_unfiltered(home):
    return extractFieldFile(_input("processor0", "0", "U"), columnNames=_U, filterInternalPatches=False, processor=0)


@_frame
def extract_C(home):
    return extractFieldFile(_input("0", "C"), columnNames=["Cx", "Cy", "Cz"])


@_frame
def extract_T(home):
    return extractFieldFile(_input("0", "T"), columnNames=["T"])


@_frame
def extract_p(home):
    return extractFieldFile(_input("0", "p"), columnNames=["p"])


@_frame
def getDataFrame_single_U(home):
    return home.readFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=False).getDataFrame()


@_frame
def getDataFrame_parallel_U(home):
    return home.readFieldFromCase("U", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=True).getDataFrame()


@_frame
def getDataFrame_single_T(home):
    return home.readFieldFromCase("T", FLOWTYPE_INCOMPRESSIBLE, INPUT_CASE, readParallel=False).getDataFrame()


@_frame
def getDataFrame_scalar_from_dataframe(home):
    field = home.getEmptyField("p", FLOWTYPE_INCOMPRESSIBLE)
    field.setFieldFromDataFrame(_scalarFrame())
    return field.getDataFrame()


@_frame
def readFieldAsDataFrame_single_U(home):
    return home.readFieldAsDataFrame("U", INPUT_CASE, times=0, readParallel=False)


@_frame
def readFieldAsDataFrame_parallel_U(home):
    return home.readFieldAsDataFrame("U", INPUT_CASE, times=0, readParallel=True)


@_frame
def readFieldAsDataFrame_single_cellCenters(home):
    return home.readFieldAsDataFrame("cellCenters", INPUT_CASE, times=0, readParallel=False)


@_frame
def readFieldAsDataFrame_parallel_p(home):
    return home.readFieldAsDataFrame("p", INPUT_CASE, times=0, readParallel=True)


# --------------------------------------------------------------------------- #
#  Normalization shared by the generator and the test                          #
# --------------------------------------------------------------------------- #
_SORT_KEYS = ("time", "processor", "region", "boundary", "processorIndex")


def frameToCSV(df):
    """CSV text of ``df`` with a stable row order (processor dirs come from glob)."""
    df = df.reset_index(drop=True)
    keys = [k for k in _SORT_KEYS if k in df.columns]
    df = df.sort_values(keys, na_position="first", kind="stable") if keys else df
    return df.to_csv(index=False)


def readCSV(text):
    return pandas.read_csv(io.StringIO(text))


def main():
    """Regenerate the PyFoam goldens. Only meaningful with the PyFoam implementation."""
    import PyFoam  # noqa: F401 -- fail loudly rather than freeze the foamlib output as "PyFoam"

    shutil.rmtree(PYFOAM_DIR, ignore_errors=True)
    os.makedirs(os.path.join(PYFOAM_DIR, "write"))
    os.makedirs(os.path.join(PYFOAM_DIR, "frames"))
    for name, fn in WRITE_SCENARIOS.items():
        out = os.path.join(PYFOAM_DIR, "write", name)
        try:
            fn(OFObjectHome(), out)
        except Exception as e:  # the old behaviour is part of the record
            shutil.rmtree(out, ignore_errors=True)
            with open(out + ".error", "w") as f:
                f.write(f"{type(e).__name__}\n")
    for name, fn in FRAME_SCENARIOS.items():
        base = os.path.join(PYFOAM_DIR, "frames", name)
        try:
            text = frameToCSV(fn(OFObjectHome()))
        except Exception as e:
            with open(base + ".error", "w") as f:
                f.write(f"{type(e).__name__}\n")
        else:
            with open(base + ".csv", "w") as f:
                f.write(text)


if __name__ == "__main__":
    main()
