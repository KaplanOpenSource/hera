"""
Regression guard for issue #935 — CLI startup latency.

`hera-riskassessment` and `hera-experiment` used to take ~5s and ~2.6s just to
print --help, because their import chains pulled pint (via
hera.utils.unitHandler) and argos eagerly. Both are now deferred to call time.

Each check runs in a fresh interpreter: the expensive modules must be absent
from sys.modules after importing what the CLI entry point imports.
"""
import subprocess
import sys

import pytest

# (entry point import, modules that must NOT be loaded as a side effect)
CASES = [
    ("hera.riskassessment.CLI", ("pint",)),
    ("hera.measurements.experiment.CLI", ("pint", "argos")),
]


@pytest.mark.parametrize("module,forbidden", CASES, ids=[c[0] for c in CASES])
def test_cli_import_does_not_pull_heavy_dependencies(module, forbidden):
    # Imported modules may print banners of their own, so the answer is tagged
    # and picked out of stdout rather than read as the whole of it.
    script = (
        f"import importlib, sys; importlib.import_module({module!r}); "
        f"print('LOADED:' + ','.join(sorted(m for m in {forbidden!r} if m in sys.modules)))"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, timeout=120
    )
    assert result.returncode == 0, result.stderr
    tagged = [ln for ln in result.stdout.splitlines() if ln.startswith("LOADED:")]
    assert len(tagged) == 1, f"no result line in stdout: {result.stdout!r}"
    loaded = tagged[0][len("LOADED:"):]
    assert not loaded, f"{module} eagerly imported: {loaded}"


def test_riskassessment_lazy_names_still_resolve():
    """The lazy __init__ must expose exactly the names the eager one did."""
    import hera.riskassessment as ra

    for name in ra.__all__:
        assert getattr(ra, name) is not None
    assert set(ra.__all__) <= set(dir(ra))
    with pytest.raises(AttributeError):
        ra.thisDoesNotExist
