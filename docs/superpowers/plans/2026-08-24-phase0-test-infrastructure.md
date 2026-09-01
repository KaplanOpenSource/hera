# Phase 0 — תשתית טסטים: תוכנית ביצוע

> **לעובדים אג'נטיים:** REQUIRED SUB-SKILL: השתמשו ב-`superpowers:subagent-driven-development` (מומלץ) או `superpowers:executing-plans` לביצוע צעד-אחר-צעד. הצעדים מסומנים ב-`- [ ]` למעקב.

**מטרה:** להקים שכבת טסטים הרמטית (`hera/tests/unit/`) שרצה בלי MongoDB, בלי S3 ובלי רשת, יחד עם מדידת coverage ו-gate ב-CI — בלי לשנות שורת קוד ייצור אחת.

**ארכיטקטורה:** כל הבידוד נשען על seam אחד — `hera.datalayer.document.connectToDatabase` — שמוחלף בזמן טסט בגרסה שמפנה את mongoengine ל-`mongomock`. ה-seam מותקן ב-`hera/tests/unit/conftest.py` בסדר קפדני, כי `hera/datalayer/__init__.py` בונה collections בזמן import. שלושה guards (רשת, ספריית בית, matplotlib) הופכים את הבידוד מהצהרה לאכיפה.

**Tech Stack:** pytest · pytest-cov · coverage.py · mongomock 4.3.0 · mongoengine · pymongo

**Spec:** [`docs/superpowers/specs/2026-08-24-unit-test-expansion-design.md`](../specs/2026-08-24-unit-test-expansion-design.md)

## Global Constraints

אלה חלים על **כל** משימה בתוכנית. הפרה שלהם היא כשל של המשימה, גם אם הטסטים עוברים.

1. **אפס שינוי בקוד ייצור.** אין לגעת באף קובץ תחת `hera/` פרט ל-`hera/tests/`. מותר לשנות: `hera/tests/**`, `.coveragerc`, `coverage_floor.txt`, `requirements-test.txt`, `Makefile`, `.github/workflows/ci.yml`.
2. **`requirements.txt` לא נגעת.** תלויות טסט נכנסות ל-`requirements-test.txt` חדש בלבד, כדי שהתקנות ייצור לא ישתנו.
3. **באגים מדווחים, לא מתוקנים.** כל פגם שנמצא הופך ל-issue (משימה 9). אין תיקוני קוד ייצור ב-Phase 0.
4. **`--require-data` פעיל רק ב-CI.** ברירת המחדל היא `skip` כמו היום. הרצת מפתח מקומית עם `TEST_HERA` חסר או חלקי **ממשיכה לדלג ולא נכשלת**. ההפעלה היא דרך `HERA_REQUIRE_TEST_DATA=1` (נקבע ב-CI) או `--require-data` (opt-in ידני).
5. **ענף ייעודי, בלי PR.** כל העבודה בענף `tests/phase0-infra`. בסוף — `git push -u origin tests/phase0-infra` בלבד. **אין לפתוח Pull Request** ללא אישור מפורש.
6. **תקציב זמן לשכבת ה-unit: מתחת ל-60 שניות.** נמדד במשימה 7 ונאכף ב-CI.
7. **הכל בעברית בתיעוד, קוד וזהויות באנגלית.**

---

## מפת קבצים

| קובץ | אחריות | משימה |
|---|---|---|
| `requirements-test.txt` | תלויות טסט בלבד: `pytest-cov`, `mongomock` | 1 |
| `.coveragerc` | הגדרות coverage: branch, source, omit, exclude_lines | 1 |
| `coverage_floor.txt` | רצפת ה-gate, מספר בודד בשליטת גרסה | 1, 7 |
| `hera/tests/unit/__init__.py` | הפיכת התיקייה לחבילה, לייבוא אבסולוטי דטרמיניסטי | 2 |
| `hera/tests/unit/_stubs.py` | `sys.modules` stubs לחבילות שאינן ניתנות להתקנה ב-CI | 2 |
| `hera/tests/unit/_seam.py` | בידוד HOME, ה-seam של mongomock, איפוס DB | 3, 5 |
| `hera/tests/unit/conftest.py` | bootstrap בסדר הנכון + שלושת ה-guards + fixtures | 3, 4, 5 |
| `hera/tests/unit/_factories.py` | בוני DataFrame / GeoDataFrame / xarray סינתטיים | 5 |
| `hera/tests/unit/test_infra_*.py` | טסטים שמאמתים את התשתית עצמה | 2-6 |
| `hera/tests/conftest.py` | תוספת בלבד: `--require-data` + helper. שאר הקובץ לא נגע | 6 |
| `Makefile` | יעדים חדשים: `test-unit`, `coverage-unit`, `coverage` | 7 |
| `.github/workflows/ci.yml` | job חדש `unit` לפני `test`; `HERA_REQUIRE_TEST_DATA` ב-`test` | 8 |

---

## הרקע שכל מבצע חייב להכיר

שלוש עובדות שנמדדו בקוד ובלעדיהן הצעדים למטה נראים שרירותיים:

**א. אי אפשר לבנות toolkit בלי DB.** `abstractToolkit` יורש מ-`Project` (`hera/toolkit.py:33`), ו-`Project.__init__` (`hera/datalayer/project.py:142`) בונה ארבעה collections, **קורא** מה-DB (`getConfig`), **כותב** ל-DB (`setConfig`) ועושה `os.makedirs` תחת `~/.hera/`.

**ב. `import hera.datalayer` נכשל בלי קובץ config.** `hera/datalayer/__init__.py:7-10` בונה collections ברמת המודול. בלי רשומה למשתמש הנוכחי ב-`~/.pyhera/config.json` נזרק `KeyError: user <name> not found`. `import hera` לבד תקין — `hera/__init__.py` עצל לחלוטין. **לכן בידוד ה-HOME וכתיבת ה-config חייבים לקרות לפני ה-import הראשון של `hera.datalayer`.**

**ג. `mongoengine.connect()` עצל.** הוא לא פותח חיבור, ולכן קובץ config עם IP מדומה מספיק כדי שה-import יעבור בלי Mongo רץ.

---

## Task 1: ענף, תלויות טסט, והגדרות coverage

**Files:**
- Create: `requirements-test.txt`
- Create: `.coveragerc`
- Create: `coverage_floor.txt`

**Interfaces:**
- Consumes: אין.
- Produces: `requirements-test.txt` (מותקן ע"י CI ב-משימה 8); `.coveragerc` (נקרא אוטומטית ע"י coverage.py); `coverage_floor.txt` (מכיל מספר שלם, נקרא ב-משימות 7 ו-8).

- [ ] **Step 1: יצירת הענף**

```bash
git checkout master
git pull
git checkout -b tests/phase0-infra
```

- [ ] **Step 2: כתיבת `requirements-test.txt`**

```
# Test-only dependencies.
#
# Deliberately NOT in requirements.txt: production installs must stay
# byte-identical.  CI installs both files; developers install this one
# only when they want to run the unit layer or measure coverage.
#
# mongomock 4.3.0 is the version the seam in hera/tests/unit/_seam.py was
# verified against.  Note that mongoengine removed the `mongomock://` URI
# scheme, so the seam uses `mongo_client_class=` instead — that keyword
# works on both the pinned (0.29.1) and the currently installed (0.29.3)
# mongoengine, which is why it was chosen.
pytest-cov==5.0.0
mongomock==4.3.0
```

- [ ] **Step 3: כתיבת `.coveragerc`**

```ini
[run]
branch = True
source = hera
parallel = True
omit =
    hera/tests/*
    hera/doc/*
    hera/bin/*
    */__pycache__/*

[report]
show_missing = True
skip_covered = False
exclude_lines =
    pragma: no cover
    if TYPE_CHECKING:
    raise NotImplementedError
    if __name__ == .__main__.:

[paths]
source =
    hera
    */hera
```

`parallel = True` גורם לכל הרצה לכתוב `.coverage.<host>.<pid>` נפרד, ש-`coverage combine` ממזג. זה מה שמאפשר למדוד את שכבת ה-unit ואת האינטגרציה בשתי הרצות ולקבל מספר אחד.

- [ ] **Step 4: יצירת `coverage_floor.txt` עם ערך זמני**

```bash
echo "0" > coverage_floor.txt
```

הערך האמיתי נקבע במשימה 7, אחרי מדידה. `0` כאן הוא placeholder מכוון ולא ניחוש של יעד.

- [ ] **Step 5: התקנת התלויות ואימות**

```bash
pip install -r requirements-test.txt
python -c "import mongomock, pytest_cov; print('mongomock', mongomock.__version__)"
coverage --version
```

Expected: `mongomock 4.3.0` ומספר גרסה של coverage.

- [ ] **Step 6: Commit**

```bash
git add requirements-test.txt .coveragerc coverage_floor.txt
git commit -m "test: add test-only dependencies and coverage configuration

requirements.txt is deliberately left untouched so production installs do
not change.  coverage_floor.txt holds a placeholder until the baseline is
measured."
```

---

## Task 2: מודול ה-stubs

**Files:**
- Create: `hera/tests/unit/__init__.py`
- Create: `hera/tests/unit/_stubs.py`
- Test: `hera/tests/unit/test_infra_stubs.py`

**Interfaces:**
- Consumes: אין.
- Produces: `_stubs.install() -> None` — אידמפוטנטי, רושם רשומות ב-`sys.modules` לחבילות שאינן מותקנות. נקרא מ-`conftest.py` (משימה 3) לפני איסוף הטסטים.

**רקע:** 19 מודולים ב-hera מייבאים ברמת המודול חבילות שלא מותקנות ב-CI. הטכניקה קיימת ומוכחת ב-`hera/tests/test_toolkit_coverage.py:26-60`. **הקובץ ההוא לא נגע** — אנחנו יוצרים מודול נפרד שניתן לייבוא. הכפילות הזמנית מכוונת: איחוד ידרוש שינוי בטסט עובד, וזה מחוץ להיקף Phase 0.

**מה לא מסטבים:** `torch`. `hera/simulations/machineLearningDeepLearning/torch/modelContainer.py` יוצר קלאסים שיורשים מ-`torch.nn.Module`, וירושה מ-`MagicMock` זורקת `TypeError`. זה נושא של אצווה 9, לא של Phase 0.

- [ ] **Step 1: יצירת החבילה ואימות שהיא נאספת**

```bash
mkdir -p hera/tests/unit
printf '' > hera/tests/unit/__init__.py
python -m pytest hera/tests/unit -q
```

Expected: `no tests ran` — התיקייה נאספת ואין בה עדיין טסטים. **לא** `ERROR: file or directory not found`.

`__init__.py` נדרש כדי ש-`from hera.tests.unit import _seam` יהיה ייבוא אבסולוטי דטרמיניסטי; `hera/tests/__init__.py` כבר קיים (ריק), כך ש-`hera.tests.unit` הופכת לחבילה תקינה.

- [ ] **Step 2: כתיבת `_stubs.py`**

```python
"""``sys.modules`` stubs for packages that cannot be installed in CI.

Nineteen hera modules import PyFoam, paraview, FreeCAD, hermes, argos or
evtk at module level.  Registering placeholder modules before those imports
run lets the surrounding pure-Python logic be tested without the external
binaries.

Deliberately NOT stubbed: ``torch``.  modelContainer.py subclasses
``torch.nn.Module``, and subclassing a MagicMock raises TypeError.  Batch 9
addresses that separately.
"""
import sys
import types
from unittest.mock import MagicMock

# Namespace packages: must expose __path__ so that `import a.b` works.
_NAMESPACE_PACKAGES = (
    "PyFoam",
    "PyFoam.RunDictionary",
    "PyFoam.Basics",
    "paraview",
    "paraview.simple",
)

# Leaf modules whose attributes are only ever called, never subclassed.
_MOCK_MODULES = (
    "PyFoam.RunDictionary.ParsedParameterFile",
    "PyFoam.RunDictionary.BoundaryDict",
    "PyFoam.Basics.DataStructures",
    "evtk",
    "evtk.hl",
    "dask.distributed",
)

# Stubbed only when genuinely absent, so a real local install wins.
_STUB_IF_MISSING = ("hermes", "argos", "FreeCAD")


def _namespace_module(name):
    module = types.ModuleType(name)
    module.__path__ = []
    module.__package__ = name
    sys.modules[name] = module
    return module


def install():
    """Register every stub. Idempotent — safe to call more than once."""
    for name in _NAMESPACE_PACKAGES:
        if name not in sys.modules:
            _namespace_module(name)

    for name in _MOCK_MODULES:
        if name not in sys.modules:
            sys.modules[name] = MagicMock()

    for name in _STUB_IF_MISSING:
        if name in sys.modules:
            continue
        try:
            __import__(name)
        except Exception:
            sys.modules[name] = MagicMock()

    # dask.dataframe must be imported before analysis.addDatesColumns runs.
    # Order matters here; see test_toolkit_coverage.py's original note.
    import dask.dataframe  # noqa: F401
```

- [ ] **Step 3: כתיבת הטסטים**

`hera/tests/unit/test_infra_stubs.py`:

```python
"""The stub layer must let hera's externally-dependent modules import."""
import sys

import pytest


@pytest.mark.unit
def test_stub_modules_are_registered():
    """install() is called by conftest at collection time."""
    for name in ("PyFoam", "paraview", "evtk", "PyFoam.Basics.DataStructures"):
        assert name in sys.modules, f"{name} was not stubbed"


@pytest.mark.unit
def test_namespace_stub_supports_submodule_import():
    """A namespace stub must carry __path__ or `import a.b` fails."""
    assert sys.modules["PyFoam"].__path__ == []


@pytest.mark.unit
def test_install_is_idempotent():
    from hera.tests.unit import _stubs

    before = sys.modules["PyFoam"]
    _stubs.install()
    assert sys.modules["PyFoam"] is before, "install() replaced an existing stub"


@pytest.mark.unit
def test_openfoam_toolkit_imports_under_stubs():
    """The whole point: a module that needs PyFoam+evtk must import."""
    from hera.simulations.openFoam.toolkit import OFToolkit

    assert OFToolkit is not None
```

- [ ] **Step 4: הרצה — הטסטים חייבים להיכשל**

```bash
python -m pytest hera/tests/unit/test_infra_stubs.py -q
```

Expected: FAIL — `install()` עוד לא נקרא מאף מקום (`conftest.py` נוצר במשימה 3). זה נכון ומכוון: הכשל מאמת שהטסטים באמת בודקים משהו.

- [ ] **Step 5: Commit**

```bash
git add hera/tests/unit/__init__.py hera/tests/unit/_stubs.py hera/tests/unit/test_infra_stubs.py
git commit -m "test: add sys.modules stub layer for CI-unavailable packages

Extracted as an importable module from the proven inline block in
test_toolkit_coverage.py, which is left untouched.  Tests fail until
conftest.py wires install() in (next task)."
```

---

## Task 3: ה-seam של mongomck ו-bootstrap ה-conftest

**Files:**
- Create: `hera/tests/unit/_seam.py`
- Create: `hera/tests/unit/conftest.py`
- Test: `hera/tests/unit/test_infra_seam.py`

**Interfaces:**
- Consumes: `_stubs.install()` מ-משימה 2.
- Produces:
  - `_seam.UNIT_DB_NAME: str` = `"heraUnitTest"`
  - `_seam.UNIT_DB_ALIAS: str` = `"heraUnitTest-alias"`
  - `_seam.isolate_home() -> str` — מחזיר את נתיב ה-HOME הזמני. **חייב להיקרא לפני כל `import hera.datalayer`.**
  - `_seam.install() -> None` — מתקין את ה-seam ורושם מחדש את החיבור תחת שם המשתמש הנוכחי.
  - `_seam.reset() -> None` — מוחק את ה-DB בזיכרון. בשימוש ב-משימה 5.

**סדר ההקמה — לא ניתן לשינוי.** הוא נגזר מעובדה ב' ברקע למעלה, ואומת מקצה-לקצה בסביבה נקייה בלי Mongo ובלי config:

```
1. isolate_home()   — HOME זמני + ~/.pyhera/config.json עם IP מדומה
2. _stubs.install() — לפני ייבוא מודולים שתלויים בחבילות חסרות
3. import hera.datalayer.document  — גוף החבילה רץ כאן, ומצליח בזכות (1)
4. החלפת connectToDatabase        — ה-seam
5. createDBConnection(user)       — רישום מחדש מול mongomock
6. rebind של ה-singletons          — ראה הערה למטה
```

**מדוע שלב 6 קיים:** `hera/datalayer/__init__.py:7-10` קושר `Measurements`/`Simulations`/`Cache`/`All` בזמן import — לפני שה-seam קיים. אומת ש-`hera/measurements/meteorology/highfreqdata/analysis/abstractcalculator.py` שורות 172, 190 ו-221 משתמש ב-`datalayer.Cache` **ישירות**. בלי ה-rebind, טסטים באזור הזה יעקפו את mongomock בשקט.

- [ ] **Step 1: כתיבת הטסט הנכשל**

`hera/tests/unit/test_infra_seam.py`:

```python
"""The seam must route the entire Project/toolkit stack at mongomock."""
import pytest


@pytest.mark.unit
def test_project_constructs_without_mongodb(unit_project):
    assert unit_project.projectName == "UNIT_TEST_PROJECT"


@pytest.mark.unit
def test_document_roundtrip(unit_project):
    unit_project.addMeasurementsDocument(
        resource="/synthetic/a.parquet",
        dataFormat="parquet",
        type="Probe",
        desc=dict(station="A"),
    )
    found = unit_project.getMeasurementsDocuments(type="Probe", station="A")
    assert len(found) == 1
    assert found[0].resource == "/synthetic/a.parquet"


@pytest.mark.unit
def test_config_roundtrip(unit_project):
    unit_project.setConfig(alpha=1)
    assert unit_project.getConfig()["alpha"] == 1


@pytest.mark.unit
def test_real_toolkit_is_usable(unit_toolkit_factory):
    """A real toolkit, built through toolkitHome, backed by mongomock."""
    from hera import toolkitHome

    toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
    toolkit.addDataSource("SRC", "/rel/path.tif", "GeoTIFF", version=[0, 0, 1])
    names = [d.desc.get("datasourceName") for d in toolkit.getDataSourceDocuments("SRC")]
    assert names == ["SRC"]


@pytest.mark.unit
def test_module_level_singletons_are_rebound():
    """abstractcalculator.py uses datalayer.Cache directly; it must be live."""
    import hera.datalayer as dl
    from hera.datalayer.document import getDBObject
    import getpass

    assert dl.Cache._metadataCol is getDBObject("Cache", getpass.getuser())
```

- [ ] **Step 2: הרצה — כשל מצופה**

```bash
python -m pytest hera/tests/unit/test_infra_seam.py -q
```

Expected: FAIL — `fixture 'unit_project' not found`.

- [ ] **Step 3: כתיבת `_seam.py`**

```python
"""Test-only seam that points hera's datalayer at an in-memory MongoDB.

Why a seam and not mocks: ``abstractToolkit`` inherits from ``Project``
(hera/toolkit.py:33), and ``Project.__init__`` reads from the DB, writes to
the DB and creates a directory under ``~/.hera``.  No toolkit can be built
without a database.  Rather than mock every call site, we substitute the one
function that creates the connection, so the real datalayer code — query
building, desc filtering, version resolution, document tagging — runs for
real against an in-memory store.

The substituted function is ``hera.datalayer.document.connectToDatabase``.
``createDBConnection`` calls it by module-global name, so rebinding the
attribute is enough.
"""
import getpass
import json
import os
import pathlib
import tempfile

UNIT_DB_NAME = "heraUnitTest"
UNIT_DB_ALIAS = f"{UNIT_DB_NAME}-alias"

_UNIT_MONGO_CONFIG = dict(
    dbIP="127.0.0.1",
    dbName=UNIT_DB_NAME,
    username="unit",
    password="unit",
)


def isolate_home():
    """Redirect HOME to a temp dir and write a placeholder pyhera config.

    MUST be called before the first ``import hera.datalayer``:
    hera/datalayer/__init__.py:7-10 builds collection singletons at import
    time, and ``getDBObject`` raises
    ``KeyError: user <name> not found`` when ~/.pyhera/config.json has no
    entry for the current user.

    The placeholder dbIP is never dialled — ``mongoengine.connect()`` is
    lazy, and install() replaces the connection before any query runs.

    Returns the temp HOME path.
    """
    home = tempfile.mkdtemp(prefix="hera_unit_home_")
    os.environ["HOME"] = home
    pyhera = pathlib.Path(home, ".pyhera")
    pyhera.mkdir(parents=True, exist_ok=True)
    with open(pyhera / "config.json", "w", encoding="utf-8") as handle:
        json.dump({getpass.getuser(): dict(_UNIT_MONGO_CONFIG)}, handle)
    return home


def install():
    """Replace connectToDatabase with a mongomock variant and re-register."""
    import mongomock
    from mongoengine import connect, disconnect

    import hera.datalayer.document as hdoc

    def _mongomock_connect(mongoConfig, alias=None):
        if isinstance(mongoConfig, str):
            mongoConfig = hdoc.parseConnectionString(mongoConfig)
        alias = "%s-alias" % mongoConfig["dbName"] if alias is None else alias
        disconnect(alias)
        return connect(
            alias=alias,
            db=mongoConfig["dbName"],
            mongo_client_class=mongomock.MongoClient,
            uuidRepresentation="standard",
        )

    hdoc.connectToDatabase = _mongomock_connect
    hdoc.createDBConnection(
        connectionName=getpass.getuser(),
        mongoConfig=dict(_UNIT_MONGO_CONFIG),
    )
    _rebind_module_level_collections()


def _rebind_module_level_collections():
    """Rebuild the singletons bound at import time, before the seam existed.

    hera/datalayer/__init__.py:7-10 creates Measurements/Simulations/Cache/All
    eagerly.  abstractcalculator.py:172,190,221 uses ``datalayer.Cache``
    directly, so a stale binding would silently bypass mongomock.
    """
    import hera.datalayer as dl

    dl.Measurements = dl.Measurements_Collection()
    dl.Simulations = dl.Simulations_Collection()
    dl.Cache = dl.Cache_Collection()
    dl.All = dl.AbstractCollection()


def reset():
    """Drop the in-memory database. Called between tests."""
    from mongoengine.connection import get_connection

    get_connection(UNIT_DB_ALIAS).drop_database(UNIT_DB_NAME)
```

- [ ] **Step 4: כתיבת `conftest.py`**

`hera/tests/unit/conftest.py`:

```python
"""Hermetic unit layer: no MongoDB, no S3 data, no network.

The module-level block below runs at collection time, before any test module
is imported.  Its order is load-bearing — see hera/tests/unit/_seam.py.
"""
# --- bootstrap: nothing from hera may be imported above this block ---------
from hera.tests.unit import _seam, _stubs   # safe: hera/__init__.py is lazy

_UNIT_HOME = _seam.isolate_home()
_stubs.install()
_seam.install()
# --------------------------------------------------------------------------

import getpass
import tempfile

import pytest

UNIT_PROJECT_NAME = "UNIT_TEST_PROJECT"


@pytest.fixture()
def unit_files_directory(tmp_path):
    """A per-test files directory, so nothing is written to the real home."""
    directory = tmp_path / "heraFiles"
    directory.mkdir()
    return str(directory)


@pytest.fixture()
def unit_project(unit_files_directory):
    """A real hera Project backed by the in-memory database."""
    from hera.datalayer import Project

    return Project(projectName=UNIT_PROJECT_NAME, filesDirectory=unit_files_directory)


@pytest.fixture()
def unit_toolkit_factory(unit_files_directory):
    """Build any real toolkit, backed by the in-memory database.

    Usage:
        toolkit = unit_toolkit_factory(toolkitHome.GIS_RASTER_TOPOGRAPHY)
    """
    from hera import toolkitHome

    def _build(toolkitName, projectName=UNIT_PROJECT_NAME):
        return toolkitHome.getToolkit(
            toolkitName,
            projectName=projectName,
            filesDirectory=unit_files_directory,
        )

    return _build
```

**הערה על `from hera.tests.unit import ...` בראש הקובץ:** אומת ש-`hera/__init__.py` עצל לחלוטין (הוא חושף `datalayer` דרך `__getattr__` בלבד) ו-`hera/tests/__init__.py` ריק. לכן שרשרת הייבוא הזו **לא** מפעילה את `hera/datalayer/__init__.py`, ובידוד ה-HOME קורה לפניו כנדרש.

- [ ] **Step 5: הרצה — הטסטים חייבים לעבור**

```bash
python -m pytest hera/tests/unit -q
```

Expected: כל הטסטים מ-משימה 2 ומ-משימה 3 עוברים (9 טסטים).

- [ ] **Step 6: אימות שזה עובד בלי Mongo ובלי config — הבדיקה שקובעת**

```bash
docker stop hera-mongo 2>/dev/null || true
env -u TEST_HERA HOME=/nonexistent-home python -m pytest hera/tests/unit -q
```

Expected: PASS. **אם זה נכשל, ה-bootstrap שבור** — השכבה לא הרמטית, ואין להמשיך למשימה 4.

- [ ] **Step 7: Commit**

```bash
git add hera/tests/unit/_seam.py hera/tests/unit/conftest.py hera/tests/unit/test_infra_seam.py
git commit -m "test: add mongomock seam and hermetic unit-layer bootstrap

Substitutes hera.datalayer.document.connectToDatabase so the real datalayer
runs against an in-memory store: no production code is modified.  Bootstrap
order is load-bearing because hera/datalayer/__init__.py builds collection
singletons at import time and raises KeyError without a pyhera config.

Also rebinds the module-level Measurements/Simulations/Cache/All singletons,
which are created before the seam exists and are used directly by
abstractcalculator.py:172,190,221."
```

---

## Task 4: שלושת ה-guards

**Files:**
- Modify: `hera/tests/unit/conftest.py` (הוספת fixtures, בלי שינוי הקיים)
- Test: `hera/tests/unit/test_infra_guards.py`

**Interfaces:**
- Consumes: ה-bootstrap מ-משימה 3.
- Produces: שלושה autouse fixtures — `_block_network`, `_matplotlib_agg`, `_no_home_writes`. אין להם API; הם נאכפים על כל טסט בתיקייה.

**מדוע זה נחוץ ולא קוסמטי:** נמדד שהרצת `pytest -m unit` על הסוויטה הקיימת לוקחת 66.76 שניות, ושמתוכן 60 הן `serverSelectionTimeoutMS` של pymongo — כי fixture ברמת session מנסה להתחבר ל-Mongo כבוי. אותו טסט לבד רץ ב-0.09 שניות. **בלי חסימת רשת, דרדור לתלות בשירות מתגלה כהמתנה שקטה של 30 שניות ולא ככשל.**

- [ ] **Step 1: כתיבת הטסט הנכשל**

`hera/tests/unit/test_infra_guards.py`:

```python
"""The unit layer must be unable to reach the network or the real home."""
import socket

import pytest


@pytest.mark.unit
def test_network_is_blocked():
    with pytest.raises(RuntimeError, match="network access"):
        socket.socket(socket.AF_INET, socket.SOCK_STREAM)


@pytest.mark.unit
def test_matplotlib_backend_is_agg():
    import matplotlib

    assert matplotlib.get_backend().lower() == "agg"


@pytest.mark.unit
def test_home_is_the_isolated_one():
    import os

    assert os.environ["HOME"].startswith(
        __import__("tempfile").gettempdir()
    ), "HOME was not isolated; a test could pollute the real ~/.hera"
```

- [ ] **Step 2: הרצה — כשל מצופה**

```bash
python -m pytest hera/tests/unit/test_infra_guards.py -q
```

Expected: FAIL על `test_network_is_blocked` — `socket.socket` עוד לא חסום (הוא יחזיר socket אמיתי ולא יזרוק).

- [ ] **Step 3: הוספת ה-guards ל-`conftest.py`**

הוסף בסוף `hera/tests/unit/conftest.py`:

```python
# ---------------------------------------------------------------------------
# Guards: turn the hermetic promise into something enforced
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _block_network(monkeypatch):
    """Fail loudly and instantly on any socket creation.

    Without this, a unit test that reaches for MongoDB waits out pymongo's
    30-second serverSelectionTimeoutMS and reports as merely slow.  That is
    exactly what costs the existing suite 60 of its 66 seconds.
    """
    import socket

    def _refuse(*args, **kwargs):
        raise RuntimeError(
            "network access is not permitted in the unit layer; "
            "use the unit_project / unit_toolkit_factory fixtures instead"
        )

    # mongomock never opens a socket, so nothing legitimate needs the real one.
    monkeypatch.setattr(socket, "socket", _refuse)


@pytest.fixture(autouse=True)
def _matplotlib_agg():
    """Force a headless backend and close figures, to stop state leaking."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture(autouse=True)
def _no_home_writes():
    """Assert that no test creates a .hera directory in the isolated home."""
    import pathlib

    hera_dir = pathlib.Path(_UNIT_HOME, ".hera")
    existed = hera_dir.exists()
    yield
    if not existed and hera_dir.exists():
        raise AssertionError(
            f"a test created {hera_dir}; pass filesDirectory (see the "
            "unit_files_directory fixture) instead of relying on the default"
        )
```

- [ ] **Step 4: הרצה — הטסטים חייבים לעבור**

```bash
python -m pytest hera/tests/unit -q
```

Expected: כל 12 הטסטים עוברים. **אם `test_infra_seam.py` נשבר כאן, זה אומר ש-mongomock פותח socket** — עצור ודווח, כי זה מפיל את הנחת היסוד של ה-seam.

- [ ] **Step 5: Commit**

```bash
git add hera/tests/unit/conftest.py hera/tests/unit/test_infra_guards.py
git commit -m "test: enforce the unit layer's isolation with three guards

Network access raises immediately instead of waiting out pymongo's 30s
server-selection timeout — the measured cause of 60 of the existing unit
run's 66 seconds.  Also pins matplotlib to Agg and asserts no test writes
a .hera directory."
```

---

## Task 5: איפוס DB בין טסטים ו-factories לנתונים סינתטיים

**Files:**
- Modify: `hera/tests/unit/conftest.py` (הוספת autouse fixture)
- Create: `hera/tests/unit/_factories.py`
- Test: `hera/tests/unit/test_infra_isolation.py`

**Interfaces:**
- Consumes: `_seam.reset()` מ-משימה 3.
- Produces:
  - `_factories.timeseries(n=10, freq="1min", columns=("u", "v", "w"), seed=0) -> pandas.DataFrame` — אינדקס `DatetimeIndex`, ערכים דטרמיניסטיים.
  - `_factories.elevation_grid(nx=3, ny=2, base=100.0, step=50.0) -> pandas.DataFrame` — עמודות `X`, `Y`, `Elevation`.
  - `_factories.points_geodataframe(points, crs=4326) -> geopandas.GeoDataFrame`.

**מדוע factories ולא נתוני S3:** הקלט הסינתטי הוא מה שנותן שליטה על ה-edge cases שנתוני S3 במקרה לא מכילים — NaN, שורה בודדת, טור קבוע, קלט ריק. זו הסיבה שדרגה A היא 83% מהקוד.

- [ ] **Step 1: כתיבת הטסט הנכשל**

`hera/tests/unit/test_infra_isolation.py`:

```python
"""Each test must start from an empty database."""
import pytest


@pytest.mark.unit
def test_writes_a_document(unit_project):
    unit_project.addMeasurementsDocument(
        resource="/leak-check", dataFormat="string", type="LeakCheck", desc={}
    )
    assert len(unit_project.getMeasurementsDocuments(type="LeakCheck")) == 1


@pytest.mark.unit
def test_previous_documents_are_gone(unit_project):
    """Depends on test order only in the sense that it must ALWAYS see zero."""
    assert unit_project.getMeasurementsDocuments(type="LeakCheck") == []


@pytest.mark.unit
def test_timeseries_factory_is_deterministic():
    from hera.tests.unit import _factories

    first = _factories.timeseries(n=5, seed=7)
    second = _factories.timeseries(n=5, seed=7)
    assert first.equals(second)
    assert list(first.columns) == ["u", "v", "w"]
    assert len(first) == 5


@pytest.mark.unit
def test_elevation_grid_shape():
    from hera.tests.unit import _factories

    grid = _factories.elevation_grid(nx=3, ny=2)
    assert len(grid) == 6
    assert set(grid.columns) == {"X", "Y", "Elevation"}
```

- [ ] **Step 2: הרצה — כשל מצופה**

```bash
python -m pytest hera/tests/unit/test_infra_isolation.py -q
```

Expected: FAIL — `test_previous_documents_are_gone` רואה את המסמך מהטסט הקודם, ו-`_factories` לא קיים.

- [ ] **Step 3: הוספת ה-reset ל-`conftest.py`**

הוסף ל-`hera/tests/unit/conftest.py`:

```python
@pytest.fixture(autouse=True)
def _reset_unit_database():
    """Drop the in-memory database after every test.

    Test order must never matter.  A test that only passes because an
    earlier test left a document behind is a bug, not a feature.
    """
    yield
    _seam.reset()
```

- [ ] **Step 4: כתיבת `_factories.py`**

```python
"""Deterministic synthetic inputs for the unit layer.

Synthetic input is not a compromise here — it is what gives control over
the edge cases the S3 data set happens not to contain: NaN, a single row,
a constant column, an empty frame.
"""
import numpy as np
import pandas as pd


def timeseries(n=10, freq="1min", columns=("u", "v", "w"), seed=0):
    """A DatetimeIndex frame with reproducible pseudo-random values."""
    rng = np.random.default_rng(seed)
    index = pd.date_range("2020-01-01 00:00:00", periods=n, freq=freq)
    return pd.DataFrame(
        {name: rng.standard_normal(n) for name in columns},
        index=index,
    )


def elevation_grid(nx=3, ny=2, base=100.0, step=50.0):
    """A tabular elevation grid with columns X, Y, Elevation."""
    xs = np.arange(1.0, nx + 1.0)
    ys = np.arange(10.0, 10.0 * (ny + 1), 10.0)
    rows = []
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            rows.append((x, y, base + step * i + 10.0 * j))
    return pd.DataFrame(rows, columns=["X", "Y", "Elevation"])


def points_geodataframe(points, crs=4326):
    """A GeoDataFrame of point geometries from (lon, lat) tuples."""
    import geopandas as gpd
    from shapely.geometry import Point

    return gpd.GeoDataFrame(
        {"id": list(range(len(points)))},
        geometry=[Point(x, y) for x, y in points],
        crs=crs,
    )
```

- [ ] **Step 5: הרצה — הטסטים חייבים לעבור, בשני סדרים**

```bash
python -m pytest hera/tests/unit -q
python -m pytest hera/tests/unit -q -p no:randomly 2>/dev/null || true
python -m pytest hera/tests/unit/test_infra_isolation.py::test_previous_documents_are_gone \
                hera/tests/unit/test_infra_isolation.py::test_writes_a_document -q
```

Expected: כל ההרצות עוברות. ההרצה השלישית הופכת את הסדר במפורש — היא מאמתת שהבידוד אמיתי ולא תוצר של סדר האיסוף.

- [ ] **Step 6: Commit**

```bash
git add hera/tests/unit/conftest.py hera/tests/unit/_factories.py hera/tests/unit/test_infra_isolation.py
git commit -m "test: reset the in-memory DB between tests, add data factories

Test order must never matter.  Factories give deterministic control over
the edge cases the S3 data set does not contain."
```

---

## Task 6: `--require-data` — כשל קשה ב-CI, skip מקומי

**Files:**
- Modify: `hera/tests/conftest.py` — **תוספות בלבד**: `parser.addoption` נוסף בתוך `pytest_addoption` הקיים (שורות 159-166), פונקציית helper חדשה, החלפת ארבע קריאות `pytest.skip` בקריאה ל-helper, והוספת `request` לחתימת ארבעה fixtures.
- Test: `hera/tests/unit/test_infra_require_data.py`

**Interfaces:**
- Consumes: אין.
- Produces: `hera.tests.conftest._missing_test_data(request, reason) -> NoReturn` — מדלג כברירת מחדל, נכשל כשההפעלה דלוקה.

**האילוץ המנחה (Global Constraint 4):** ברירת המחדל היא בדיוק ההתנהגות הקיימת. מפתח שמריץ מקומית בלי `TEST_HERA` מקבל skip. רק `HERA_REQUIRE_TEST_DATA=1` (שנקבע ב-CI במשימה 8) או `--require-data` הופכים את זה לכשל.

**המוטיבציה:** ארבע נקודות ב-`hera/tests/conftest.py` (שורות 226, 235, 254, 276) מדלגות כשחסרים `TEST_HERA`, `data_config.json`, ה-result set או `test_repository.json`. אם ההורדה מ-S3 נכשלת חלקית, 76 טסטי האינטגרציה מדולגים בשקט וה-CI מדווח הצלחה.

- [ ] **Step 1: כתיבת הטסט הנכשל**

`hera/tests/unit/test_infra_require_data.py`:

```python
"""The missing-data gate must skip locally and fail only when armed."""
import pytest


class _FakeConfig:
    def __init__(self, require_data):
        self._require_data = require_data

    def getoption(self, name):
        if name == "--require-data":
            return self._require_data
        raise KeyError(name)


class _FakeRequest:
    def __init__(self, require_data):
        self.config = _FakeConfig(require_data)


@pytest.mark.unit
def test_skips_by_default(monkeypatch):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.delenv("HERA_REQUIRE_TEST_DATA", raising=False)
    with pytest.raises(pytest.skip.Exception):
        _missing_test_data(_FakeRequest(False), "no data")


@pytest.mark.unit
def test_fails_when_cli_flag_is_set(monkeypatch):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.delenv("HERA_REQUIRE_TEST_DATA", raising=False)
    with pytest.raises(pytest.fail.Exception):
        _missing_test_data(_FakeRequest(True), "no data")


@pytest.mark.unit
@pytest.mark.parametrize("value", ["1", "true", "True", "yes"])
def test_fails_when_env_var_is_set(monkeypatch, value):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.setenv("HERA_REQUIRE_TEST_DATA", value)
    with pytest.raises(pytest.fail.Exception):
        _missing_test_data(_FakeRequest(False), "no data")


@pytest.mark.unit
@pytest.mark.parametrize("value", ["", "0", "false", "no"])
def test_still_skips_for_falsey_env_values(monkeypatch, value):
    from hera.tests.conftest import _missing_test_data

    monkeypatch.setenv("HERA_REQUIRE_TEST_DATA", value)
    with pytest.raises(pytest.skip.Exception):
        _missing_test_data(_FakeRequest(False), "no data")
```

- [ ] **Step 2: הרצה — כשל מצופה**

```bash
python -m pytest hera/tests/unit/test_infra_require_data.py -q
```

Expected: FAIL — `ImportError: cannot import name '_missing_test_data'`.

- [ ] **Step 3: הוספת האופציה ב-`hera/tests/conftest.py`**

בתוך `pytest_addoption` הקיים (שורות 159-166), **אחרי** קריאת ה-`--result-set` הקיימת, הוסף:

```python
    parser.addoption(
        "--require-data",
        action="store_true",
        default=False,
        help="Fail instead of skipping when the TEST_HERA data set is "
             "incomplete.  Intended for CI; local runs keep skipping.",
    )
```

- [ ] **Step 4: הוספת ה-helper**

מיד אחרי `pytest_addoption`, הוסף:

```python
_REQUIRE_DATA_TRUTHY = frozenset({"1", "true", "yes"})


def _missing_test_data(request, reason):
    """Skip locally; fail when the missing-data gate is armed.

    A local developer run with an absent or partial TEST_HERA must keep
    skipping — that is the historical behaviour and it stays the default.
    CI sets HERA_REQUIRE_TEST_DATA=1 so that a partial S3 download fails
    the build instead of reporting green with 76 silent skips.
    """
    armed = request.config.getoption("--require-data") or (
        os.environ.get("HERA_REQUIRE_TEST_DATA", "").strip().lower()
        in _REQUIRE_DATA_TRUTHY
    )
    if armed:
        pytest.fail(
            f"{reason} — the missing-data gate is armed "
            "(HERA_REQUIRE_TEST_DATA / --require-data)",
            pytrace=False,
        )
    pytest.skip(reason)
```

- [ ] **Step 5: חיבור ארבע נקודות הדילוג**

ארבעה שינויים נקודתיים. בכל אחד: הוסף `request` לחתימה, והחלף את `pytest.skip(...)` ב-`_missing_test_data(request, ...)`.

```python
# ~line 221
@pytest.fixture(scope="session")
def test_hera_root(request):
    root = os.environ.get("TEST_HERA", os.path.expanduser("~/hera_unittest_data"))
    root = Path(root)
    if not root.is_dir():
        _missing_test_data(request, f"TEST_HERA directory not found: {root}")
    return root


# ~line 231
@pytest.fixture(scope="session")
def data_config(request, test_hera_root):
    cfg_path = test_hera_root / "data_config.json"
    if not cfg_path.is_file():
        _missing_test_data(request, f"data_config.json not found at {cfg_path}")
    with open(cfg_path, encoding="utf-8") as fh:
        return json.load(fh)


# ~line 250
@pytest.fixture(scope="session")
def expected_dir(request, test_hera_root, result_set):
    d = test_hera_root / "expected" / result_set
    if not d.is_dir():
        _missing_test_data(request, f"Expected result set directory not found: {d}")
    return d


# ~line 263
@pytest.fixture(scope="session")
def hera_test_project(request, test_hera_root):
    # ... unchanged body up to:
    repo_json_path = test_hera_root / "test_repository.json"
    if not repo_json_path.is_file():
        _missing_test_data(request, f"test_repository.json not found at {repo_json_path}")
    # ... rest of the body unchanged
```

- [ ] **Step 6: הרצה — טסטי היחידה חייבים לעבור**

```bash
python -m pytest hera/tests/unit/test_infra_require_data.py -q
```

Expected: 11 טסטים עוברים.

- [ ] **Step 7: אימות מקצה-לקצה של שני המצבים**

```bash
echo "--- local, missing data: must SKIP, exit 0 ---"
TEST_HERA=/nonexistent python -m pytest hera/tests/test_topography.py -q; echo "exit=$?"

echo "--- CI-armed, missing data: must FAIL, exit non-zero ---"
TEST_HERA=/nonexistent HERA_REQUIRE_TEST_DATA=1 python -m pytest hera/tests/test_topography.py -q; echo "exit=$?"
```

Expected: הראשונה — skips ו-`exit=0`. השנייה — failures ו-`exit=1`. **זו ההוכחה לאילוץ 4.**

- [ ] **Step 8: אימות שהסוויטה הקיימת לא נשברה**

```bash
make mongo-up
TEST_HERA=$HOME/hera_unittest_data python -m pytest hera/tests -m "not notebook" -q
```

Expected: אותה תוצאה כמו לפני השינוי. **אם נוספו כשלים, השינוי לא היה additive** — חזור ל-Step 5.

- [ ] **Step 9: Commit**

```bash
git add hera/tests/conftest.py hera/tests/unit/test_infra_require_data.py
git commit -m "test: fail on missing test data only when explicitly armed

Four fixtures skip when TEST_HERA is absent or partial, so a failed S3
download makes CI green with 76 silent skips.  The gate is off by default:
local runs keep skipping.  CI arms it with HERA_REQUIRE_TEST_DATA=1."
```

---

## Task 7: מדידת בסיס, קביעת הרצפה, ויעדי Makefile

**Files:**
- Modify: `coverage_floor.txt`
- Modify: `Makefile`

**Interfaces:**
- Consumes: `.coveragerc` ו-`requirements-test.txt` מ-משימה 1; שכבת ה-unit ממשימות 2-5.
- Produces: יעדים `test-unit`, `coverage-unit`, `coverage`; ו-`coverage_floor.txt` עם המספר הנמדד.

- [ ] **Step 1: מדידת זמן שכבת ה-unit מול תקציב 60 השניות**

```bash
time python -m pytest hera/tests/unit -m unit -q
```

Expected: מתחת ל-60 שניות. אם לא — עצור ודווח איזה טסט אשם (`--durations=10`) לפני שממשיכים.

- [ ] **Step 2: מדידת coverage של שכבת ה-unit לבד**

```bash
rm -f .coverage .coverage.*
python -m pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
coverage combine
coverage report | tail -5
```

- [ ] **Step 3: מדידת ה-coverage המשולב**

```bash
make mongo-up
rm -f .coverage .coverage.*
python -m pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
TEST_HERA=$HOME/hera_unittest_data python -m pytest hera/tests -m "not notebook and not unit" -q --cov=hera --cov-report=
coverage combine
coverage report | tail -5
```

רשום את ה-`TOTAL`. **זו הרצפה.** היא לא נקבעת מראש כי יעד שנקבע לפני מדידה הוא ניחוש.

- [ ] **Step 4: כתיבת הרצפה הנמדדת**

```bash
# החלף <MEASURED> במספר השלם מ-Step 3, מעוגל למטה
echo "<MEASURED>" > coverage_floor.txt
cat coverage_floor.txt
```

- [ ] **Step 5: הוספת היעדים ל-`Makefile`**

הוסף `test-unit coverage-unit coverage` לרשימת `.PHONY` הקיימת, ואת היעדים אחרי `test-notebooks`:

```make
test-unit:
	PYTHONPATH=.$${PYTHONPATH:+:$$PYTHONPATH} pytest hera/tests/unit -m unit -q

coverage-unit:
	rm -f .coverage .coverage.*
	PYTHONPATH=.$${PYTHONPATH:+:$$PYTHONPATH} pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
	coverage combine
	coverage report

coverage:
	rm -f .coverage .coverage.*
	PYTHONPATH=.$${PYTHONPATH:+:$$PYTHONPATH} pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
	PYTHONPATH=.$${PYTHONPATH:+:$$PYTHONPATH} TEST_HERA=$(TEST_HERA) pytest hera/tests -m "not notebook and not unit" -q --cov=hera --cov-report=
	coverage combine
	coverage report
	coverage html -d cache/coverage_html
	@echo "floor is $$(cat coverage_floor.txt)% — see coverage_floor.txt"
```

- [ ] **Step 6: הוספת שורות עזרה**

ליד השורה הקיימת `make test` (שורה ~118) הוסף:

```make
	@echo "    make test-unit     Run the hermetic unit layer (no MongoDB, no S3)"
	@echo "    make coverage      Measure combined coverage (unit + integration)"
```

- [ ] **Step 7: אימות היעדים**

```bash
make test-unit
make coverage-unit
```

Expected: שניהם מסתיימים בהצלחה ומדפיסים דוח coverage.

- [ ] **Step 8: Commit**

```bash
git add coverage_floor.txt Makefile
git commit -m "test: measure the coverage baseline and add make targets

The floor comes from an actual combined measurement, not a chosen target.
Measurement is combined across both layers on purpose: gating on the unit
layer alone would reward mocking every boundary."
```

---

## Task 8: CI — שלב unit לפני האינטגרציה

**Files:**
- Modify: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: `requirements-test.txt`, `coverage_floor.txt`, `hera/tests/unit/`.
- Produces: job בשם `unit`; ה-job הקיים `test` מקבל `needs: unit` ו-`HERA_REQUIRE_TEST_DATA: "1"`.

- [ ] **Step 1: הוספת ה-job החדש**

הוסף לפני ה-job `test` הקיים, תחת `jobs:`:

```yaml
  unit:
    name: Unit tests (hermetic — no services)
    runs-on: ubuntu-latest

    steps:
      - name: Checkout hera
        uses: actions/checkout@v4

      - name: Set up Python 3.11
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'
          cache: 'pip'
          cache-dependency-path: |
            requirements.txt
            requirements-test.txt

      - name: Install system dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y libgdal-dev gdal-bin

      - name: Install Python dependencies
        run: |
          python -m pip install --upgrade pip setuptools wheel
          pip install -r requirements.txt
          pip install -r requirements-test.txt
          pip install -e . --no-deps

      # No MongoDB service and no S3 test data on purpose: this job proves
      # the unit layer is hermetic.  If it needs either, the seam is broken.
      - name: Run the unit layer
        run: pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
        env:
          MPLBACKEND: Agg
          JUPYTER_PLATFORM_DIRS: "1"

      - name: Report unit coverage
        run: |
          coverage combine
          coverage report

      - name: Upload coverage data
        uses: actions/upload-artifact@v4
        with:
          name: coverage-unit
          path: .coverage
          include-hidden-files: true
```

- [ ] **Step 2: חיווט ה-job הקיים `test`**

ב-job `test`: הוסף `needs: unit` מיד אחרי `runs-on`, והוסף לבלוק ה-`env:` הקיים:

```yaml
      # Arm the missing-data gate: a partial S3 download must fail the build
      # rather than skip 76 tests and report green.
      HERA_REQUIRE_TEST_DATA: "1"
```

- [ ] **Step 3: הוספת בדיקת ה-gate ל-job `test`**

הוסף אחרי הצעד `Run full test suite` הקיים (ושנה אותו כך שיאסוף coverage):

```yaml
      - name: Download unit coverage
        uses: actions/download-artifact@v4
        with:
          name: coverage-unit

      - name: Enforce the coverage floor on the combined measurement
        run: |
          coverage combine
          FLOOR="$(cat coverage_floor.txt)"
          echo "Enforcing combined coverage floor: ${FLOOR}%"
          coverage report --fail-under="${FLOOR}"
```

ושנה את הצעד `Run full test suite` להוסיף `--cov=hera --cov-report=` לפקודת ה-pytest, כך שה-coverage של האינטגרציה נאסף ומתמזג.

- [ ] **Step 4: אימות תחביר ה-YAML**

```bash
python -c "import yaml,sys; d=yaml.safe_load(open('.github/workflows/ci.yml')); print('jobs:', list(d['jobs'])); print('test needs:', d['jobs']['test'].get('needs'))"
```

Expected: `unit` מופיע ברשימת ה-jobs, ו-`test needs: unit`.

- [ ] **Step 5: הדמיית שלב ה-unit של CI מקומית**

```bash
docker stop hera-mongo 2>/dev/null || true
env -u TEST_HERA MPLBACKEND=Agg python -m pytest hera/tests/unit -m unit -q --cov=hera --cov-report=
coverage combine && coverage report | tail -3
```

Expected: PASS בלי Mongo ובלי `TEST_HERA`. זו אותה בדיקה כמו Step 6 במשימה 3, עכשיו עם coverage.

- [ ] **Step 6: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: run the hermetic unit layer before the integration suite

The unit job has no MongoDB service and no S3 data by design: if it needs
either, the seam is broken.  A logic failure now surfaces in under a minute
instead of after service setup and an S3 download.

The integration job arms HERA_REQUIRE_TEST_DATA and enforces the coverage
floor on the combined measurement of both layers."
```

---

## Task 9: דיווח הפגמים כ-issues

**Files:** אין. אין תיקוני קוד ב-Phase 0 (Global Constraint 3).

**Interfaces:**
- Consumes: הממצאים בסעיף 6 של ה-spec.
- Produces: ארבעה issues ב-GitHub. מספריהם נרשמים ב-spec כהפניה.

- [ ] **Step 1: אישור לפני פרסום**

יצירת issues היא פעולה חוצה-חוץ ופומבית. **בקש אישור מפורש לפני ההרצה**, גם אם דיווח באגים אושר באופן כללי. הצג את ארבע הכותרות והגוף, וחכה לאישור.

- [ ] **Step 2: אימות שהאותנטיקציה קיימת**

```bash
gh auth status
gh repo view KaplanOpenSource/hera --json nameWithOwner
```

- [ ] **Step 3: יצירת ה-issues**

```bash
gh issue create --repo KaplanOpenSource/hera \
  --title "dictToMongoQuery flattens list values: every list-valued query operator silently returns nothing" \
  --body 'hera/utils/query.py:dictToMongoQuery flattens a list value into indexed keys:

    station__in=["A","C"]  ->  {"desc__station__in__0": "A", "desc__station__in__1": "C"}

The query then looks for a document where `desc.station.in.0 == "A"`, which is
meaningless. The result is zero documents and no error.

Verified that mongomock and MongoDB both support `$in` on a nested field
(a direct pymongo query returned 2 of 3 documents), so this is in hera, not
in the storage engine.

Real-world impact: hera/measurements/meteorology/highfreqdata/analysis/abstractcalculator.py
lines 172 and 190 pass `params__all=self._AllCalculatedParams`, and
_AllCalculatedParams is a list (line 58 initialises it to [], line 160 extends
it). The cache lookup therefore never matches, so the high-frequency
turbulence statistics are recomputed on every call instead of being read back
from the Cache collection.

Found while writing the test-infrastructure plan; reported rather than fixed
so that Phase 0 changes no production code.'

gh issue create --repo KaplanOpenSource/hera \
  --title "toAzimuthAngle does not normalise its input and returns angles outside [0, 360)" \
  --body 'hera/utils/angle.py:4

    toAzimuthAngle = lambda mathematical_angle: (90 - mathematical_angle) if ((90 - mathematical_angle) >= 0) else (450 - mathematical_angle)

There is no `% 360` before the branch, so any input outside the first cycle
produces an invalid azimuth:

    toAzimuthAngle(800)   ->  -350
    toAzimuthAngle(-400)  ->   490

An azimuth must be in [0, 360).'

gh issue create --repo KaplanOpenSource/hera \
  --title "toMathematicalAngle and toMeteorologicalAngle are the same object, undocumented" \
  --body 'hera/utils/angle.py:2-3 binds both names to one lambda; `toMeteorologicalAngle is toMathematicalAngle` is True.

Mathematically this is sound — `(270 - x) % 360` is an involution, so the
conversion is identical in both directions. But it is undocumented, and it
means the two names provide no protection against converting in the wrong
direction, which is the specific risk CLAUDE.md calls out when it requires
explicit conversion.

Request: a docstring stating the involution, and a test that pins it, so that
a future "fix" to one of the two names fails loudly.'

gh issue create --repo KaplanOpenSource/hera \
  --title "getToolkitDocument swallows every exception and cannot signal failure" \
  --body 'hera/toolkit.py:getToolkitDocument contains three consecutive
`except Exception: pass` blocks and then returns None.

A connection error, a malformed query and a genuinely absent document are
therefore indistinguishable to the caller: all three return None. This also
makes the function untestable in the sense that matters, because no input can
make it report a problem.

Request: distinguish "not found" from "lookup failed".'
```

- [ ] **Step 4: רישום מספרי ה-issues ב-spec**

הוסף לסוף סעיף 6 של `docs/superpowers/specs/2026-08-24-unit-test-expansion-design.md`:

```markdown
**מספרי ה-issues:** 6.1 → #___ · 6.2 → #___ · 6.3 → #___ · 6.6 → #___
```

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/specs/2026-08-24-unit-test-expansion-design.md
git commit -m "docs: record the issue numbers for the defects found in Phase 0"
```

---

## Task 10: דחיפת הענף

**Files:** אין.

- [ ] **Step 1: אימות סופי של הסוויטה כולה**

```bash
make mongo-up
env -u TEST_HERA python -m pytest hera/tests/unit -m unit -q
TEST_HERA=$HOME/hera_unittest_data python -m pytest hera/tests -m "not notebook" -q
```

Expected: שכבת ה-unit עוברת בלי `TEST_HERA`; הסוויטה המלאה עוברת כמו לפני התוכנית.

- [ ] **Step 2: אימות שקוד הייצור לא נגע**

```bash
git diff --name-only master...HEAD | grep '^hera/' | grep -v '^hera/tests/' || echo "OK: no production file touched"
git diff --name-only master...HEAD
```

Expected: השורה `OK: no production file touched` מודפסת. **אם מופיע קובץ תחת `hera/` שאינו תחת `hera/tests/`, הופר Global Constraint 1** — עצור ודווח.

- [ ] **Step 3: דחיפה — בלי PR**

```bash
git push -u origin tests/phase0-infra
```

**אין להריץ `gh pr create`.** פתיחת PR דורשת אישור מפורש ונפרד (Global Constraint 5).

- [ ] **Step 4: דיווח**

דווח: שם הענף, מספר הטסטים החדשים, זמן שכבת ה-unit, הרצפה הנמדדת ב-`coverage_floor.txt`, ומספרי ארבעת ה-issues.

---

## מה Phase 0 לא עושה

- **לא כותב טסטים לפונקציות ייצור.** זו אצווה 1 והלאה. Phase 0 מספק את הכלים, וטסטי התשתית הם מה שמאמת שהכלים עובדים.
- **לא מתקן את הבאגים שנמצאו.** דווחו כ-issues (משימה 9).
- **לא מיישר גרסאות.** סעיף 6.7 ב-spec הציע לתקן הצמדות ב-`requirements.txt`; זה נפסל כאן כי הוא משפיע על התקנות ייצור. התלויות החדשות מבודדות ב-`requirements-test.txt`, וה-seam בנוי על `mongo_client_class` דווקא כדי להיות עמיד לשתי גרסאות ה-mongoengine.
- **לא מסטב `torch`.** ירושה מ-`MagicMock` נכשלת; זה נושא של אצווה 9.
- **לא מאחד את ה-stubs עם `test_toolkit_coverage.py`.** האיחוד ידרוש שינוי בטסט עובד.
