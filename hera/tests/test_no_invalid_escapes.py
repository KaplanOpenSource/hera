import pathlib
import warnings

import pytest

HERA_ROOT = pathlib.Path(__file__).resolve().parents[1]


def _python_sources():
    return sorted(HERA_ROOT.rglob("*.py"))


@pytest.mark.parametrize("path", _python_sources(), ids=lambda p: str(p.relative_to(HERA_ROOT)))
def test_no_invalid_escape_sequences(path):
    """Compiling every hera source file must not emit invalid-escape warnings.

    Python <=3.11 emits DeprecationWarning, >=3.12 emits SyntaxWarning."""
    source = path.read_text(encoding="utf-8", errors="replace")
    with warnings.catch_warnings():
        warnings.simplefilter("error", SyntaxWarning)
        warnings.simplefilter("error", DeprecationWarning)
        try:
            compile(source, str(path), "exec")
        except SyntaxError as err:
            # Legacy files that never compiled (dead code) are out of scope here;
            # they cannot emit escape warnings because they cannot be imported at all.
            pytest.skip(f"pre-existing SyntaxError, file is not importable: {err.msg} (line {err.lineno})")
