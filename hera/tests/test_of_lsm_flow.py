"""
Tests for OFLSMToolkit OpenFOAM Lagrangian particle-file parsing.

What is tested
--------------
1. Module-level import smoke test (catches import regressions like the
   broken HERAMETADATA import that used to prevent the toolkit from loading).
2. _extractFile with a multi-row vector field (globalSigmaPositions format).
3. _extractFile with a compact uniform scalar field (origProcId "N{v}" format).
4. _extractFile with a multi-row scalar field (origId format).
5. _readRecord that assembles a full per-particle DataFrame from synthetic
   lagrangian files placed in the expected directory layout.

All tests are skipped when the Hermes workflow package is not importable,
because hera.simulations.openFoam requires hermes at import time.
An extra marker skips the full-OpenFOAM-run test when blockMesh is absent.
"""

import os
import shutil
import sys
import tempfile
import textwrap

import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock

# ---------------------------------------------------------------------------
# Make the Hermes submodule importable when it lives next to the repo root
# ---------------------------------------------------------------------------
_HERMES_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "Hermes")
)
if os.path.isdir(_HERMES_ROOT) and _HERMES_ROOT not in sys.path:
    sys.path.insert(0, _HERMES_ROOT)

hermes = pytest.importorskip("hermes")

from hera.simulations.openFoam.lagrangian.LSM.toolkit import OFLSMToolkit  # noqa: E402

# ---------------------------------------------------------------------------
# Skip marker for tests that need actual OpenFOAM binaries
# ---------------------------------------------------------------------------
_blockMesh_missing = shutil.which("blockMesh") is None
requires_openfoam = pytest.mark.skipif(
    _blockMesh_missing,
    reason="OpenFOAM not installed (blockMesh not in PATH)",
)


# ---------------------------------------------------------------------------
# Helpers for building synthetic OpenFOAM lagrangian files
# ---------------------------------------------------------------------------

_FOAM_HEADER_TMPL = textwrap.dedent("""\
    /*--------------------------------*- C++ -*----------------------------------*\\
      =========                 |
      \\\\      /  F ield         | OpenFOAM: test
       \\\\    /   O peration     |
        \\\\  /    A nd           |
         \\\\/     M anipulation  |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       {cls};
        location    "{loc}";
        object      {obj};
    }}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    {count}
    (
""")

_FOAM_FOOTER = textwrap.dedent("""\
    )

    // ************************************************************************* //

""")


def _write_vector_field(path, rows, cls="vectorField"):
    """Write a multi-row OF vector field file."""
    header = _FOAM_HEADER_TMPL.format(
        cls=cls,
        loc="test",
        obj=os.path.basename(path),
        count=len(rows),
    )
    with open(path, "w") as f:
        f.write(header)
        for x, y, z in rows:
            f.write(f"({x} {y} {z})\n")
        f.write(_FOAM_FOOTER)


def _write_scalar_field(path, values, cls="labelField"):
    """Write a multi-row OF scalar field file (one value per line)."""
    header = _FOAM_HEADER_TMPL.format(
        cls=cls,
        loc="test",
        obj=os.path.basename(path),
        count=len(values),
    )
    with open(path, "w") as f:
        f.write(header)
        for v in values:
            f.write(f"{v}\n")
        f.write(_FOAM_FOOTER)


def _write_uniform_scalar_field(path, repeat, value, cls="labelField"):
    """Write a compact OF uniform scalar field: ``N{v}``.

    Compact files have ONE blank line after ``// * * *``, then the
    compact data expression on line 18 (0-indexed line 17), which is
    where _extractFile's fallback parser reads from.
    """
    header = textwrap.dedent(f"""\
        /*--------------------------------*- C++ -*----------------------------------*\\
          =========                 |
          \\\\      /  F ield         | OpenFOAM: test
           \\\\    /   O peration     |
            \\\\  /    A nd           |
             \\\\/     M anipulation  |
        \\*---------------------------------------------------------------------------*/
        FoamFile
        {{
            version     2.0;
            format      ascii;
            class       {cls};
            location    "test";
            object      {os.path.basename(path)};
        }}
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        {repeat}{{{value}}}

        // ************************************************************************* //

        """)
    with open(path, "w") as f:
        f.write(header)


# ---------------------------------------------------------------------------
# Minimal mock that lets us call unbound OFLSMToolkit methods without MongoDB
# ---------------------------------------------------------------------------

class _FakeToolkit:
    """Minimal stand-in with just enough attributes for _extractFile / _readRecord."""
    _cloudName = "kinematicCloud"
    logger = MagicMock()

    _extractFile = OFLSMToolkit._extractFile
    _readRecord = OFLSMToolkit._readRecord


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestOFLSMImport:
    """Import-level sanity checks — catches broken import chains."""

    def test_toolkit_class_is_importable(self):
        assert OFLSMToolkit is not None

    def test_toolkit_has_extract_file(self):
        assert callable(getattr(OFLSMToolkit, "_extractFile", None))

    def test_toolkit_has_read_record(self):
        assert callable(getattr(OFLSMToolkit, "_readRecord", None))

    def test_toolkit_has_load_data(self):
        assert callable(getattr(OFLSMToolkit, "loadData", None))


class TestExtractFileVector:
    """_extractFile parses multi-row vector fields correctly."""

    def setup_method(self):
        self.tmp = tempfile.mkdtemp()
        self.tk = _FakeToolkit()

    def teardown_method(self):
        shutil.rmtree(self.tmp)

    def test_returns_dataframe(self):
        path = os.path.join(self.tmp, "globalSigmaPositions")
        _write_vector_field(path, [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)])
        df = self.tk._extractFile(path, ["x", "y", "z"])
        assert isinstance(df, pd.DataFrame)

    def test_row_count(self):
        rows = [(float(i), float(i + 1), float(i + 2)) for i in range(5)]
        path = os.path.join(self.tmp, "globalSigmaPositions")
        _write_vector_field(path, rows)
        df = self.tk._extractFile(path, ["x", "y", "z"])
        assert len(df) == 5

    def test_column_names(self):
        path = os.path.join(self.tmp, "globalSigmaPositions")
        _write_vector_field(path, [(1.0, 2.0, 3.0)])
        df = self.tk._extractFile(path, ["x", "y", "z"])
        assert list(df.columns) == ["x", "y", "z"]

    def test_values_correct(self):
        rows = [(10.5, 20.5, 30.5), (11.0, 21.0, 31.0), (12.0, 22.0, 32.0)]
        path = os.path.join(self.tmp, "globalSigmaPositions")
        _write_vector_field(path, rows)
        df = self.tk._extractFile(path, ["x", "y", "z"])
        expected = pd.DataFrame(rows, columns=["x", "y", "z"])
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True),
            expected.astype(float).reset_index(drop=True),
        )


class TestExtractFileScalar:
    """_extractFile parses scalar fields — both multi-row and compact uniform."""

    def setup_method(self):
        self.tmp = tempfile.mkdtemp()
        self.tk = _FakeToolkit()

    def teardown_method(self):
        shutil.rmtree(self.tmp)

    def test_multirow_scalar_row_count(self):
        path = os.path.join(self.tmp, "origId")
        _write_scalar_field(path, list(range(7)))
        df = self.tk._extractFile(path, ["id"], vector=False)
        assert len(df) == 7

    def test_multirow_scalar_values(self):
        path = os.path.join(self.tmp, "origId")
        _write_scalar_field(path, [0, 1, 2, 3])
        df = self.tk._extractFile(path, ["id"], vector=False)
        np.testing.assert_array_equal(df["id"].values, [0.0, 1.0, 2.0, 3.0])

    def test_compact_uniform_scalar(self):
        path = os.path.join(self.tmp, "origProcId")
        _write_uniform_scalar_field(path, repeat=5, value=0)
        df = self.tk._extractFile(path, ["procId"], vector=False)
        assert len(df) == 5
        assert (df["procId"] == 0.0).all()


class TestReadRecord:
    """_readRecord assembles a complete per-particle DataFrame."""

    def setup_method(self):
        self.tmp = tempfile.mkdtemp()
        self.tk = _FakeToolkit()
        time_name = "1"
        cloud_dir = os.path.join(self.tmp, time_name, "lagrangian", "kinematicCloud")
        os.makedirs(cloud_dir)

        self._n = 4
        self._positions = [(float(i), float(i + 10), float(i + 20)) for i in range(self._n)]
        self._global_pos = [(float(i + 100), float(i + 200), float(i + 300)) for i in range(self._n)]

        _write_vector_field(
            os.path.join(cloud_dir, "globalSigmaPositions"), self._positions
        )
        _write_scalar_field(
            os.path.join(cloud_dir, "origId"), list(range(self._n)), cls="labelField"
        )
        _write_uniform_scalar_field(
            os.path.join(cloud_dir, "origProcId"), repeat=self._n, value=0
        )
        _write_vector_field(
            os.path.join(cloud_dir, "globalPositions"), self._global_pos
        )

    def teardown_method(self):
        shutil.rmtree(self.tmp)

    def test_returns_dataframe(self):
        df = self.tk._readRecord("1", self.tmp)
        assert isinstance(df, pd.DataFrame)

    def test_row_count(self):
        df = self.tk._readRecord("1", self.tmp)
        assert len(df) == self._n

    def test_required_columns_present(self):
        df = self.tk._readRecord("1", self.tmp)
        for col in ["x", "y", "z", "id", "procId", "globalX", "globalY", "globalZ", "time"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_time_column_value(self):
        df = self.tk._readRecord("1", self.tmp)
        assert (df["time"] == 1.0).all()

    def test_global_id_computed(self):
        df = self.tk._readRecord("1", self.tmp)
        assert "globalID" in df.columns
        expected_ids = 1_000_000_000 * df["procId"] + df["id"]
        pd.testing.assert_series_equal(df["globalID"], expected_ids, check_names=False)

    def test_positions_match(self):
        df = self.tk._readRecord("1", self.tmp)
        expected_x = [p[0] for p in self._positions]
        np.testing.assert_allclose(df["x"].values, expected_x)

    def test_global_positions_match(self):
        df = self.tk._readRecord("1", self.tmp)
        expected_gx = [p[0] for p in self._global_pos]
        np.testing.assert_allclose(df["globalX"].values, expected_gx)


class TestOpenFOAMIntegration:
    """Full OpenFOAM run test — skipped unless blockMesh is in PATH."""

    @requires_openfoam
    def test_blockMesh_runs(self):
        # If this point is reached, OpenFOAM is installed; confirm binaries work
        import subprocess
        result = subprocess.run(["blockMesh", "-help"], capture_output=True, text=True)
        # blockMesh returns non-zero for --help on older OF; just check it runs
        assert result.returncode in (0, 1)
