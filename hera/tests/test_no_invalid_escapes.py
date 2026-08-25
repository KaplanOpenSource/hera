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

    # Two compiles, deliberately.  With SyntaxWarning escalated to an error,
    # Python surfaces an invalid escape AS a SyntaxError -- so a single
    # compile inside the filter cannot tell "this file does not parse" from
    # "this file has the very defect we are looking for", and the handler for
    # the former used to swallow the latter into a skip.  The failure this
    # test exists to detect could therefore never be reported.
    #
    # Pass 1 decides parseability, with no filter active.
    try:
        compile(source, str(path), "exec")
    except SyntaxError as err:
        pytest.skip(
            f"pre-existing SyntaxError, file is not importable: "
            f"{err.msg} (line {err.lineno})"
        )

    # Pass 2 looks for escape warnings in a file that is known to parse, so
    # anything raised here is the real thing.
    with warnings.catch_warnings():
        warnings.simplefilter("error", SyntaxWarning)
        warnings.simplefilter("error", DeprecationWarning)
        compile(source, str(path), "exec")
