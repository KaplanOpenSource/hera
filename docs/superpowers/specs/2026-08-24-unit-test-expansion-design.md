# הרחבת סוויטת הטסטים של Hera — מסמך עיצוב

**תאריך:** 2026-08-24
**סטטוס:** לסקירה
**ענף:** `issue971` (התוכנית עצמה תבוצע בענפים נפרדים לפי אצווה)

---

## 1. תקציר מנהלים

המטרה: כיסוי טסטים מקיף ל-1,333 הפונקציות הציבוריות של Hera, עם **יעד coverage מדיד ו-gate ב-CI** שמונע נסיגה.

החקירה שקדמה למסמך הזה שינתה את הגישה המומלצת בשתי נקודות מהותיות:

1. **נמצא seam אחד** שדרכו כל מחסנית ה-`Project`/toolkit רצה בזיכרון בלי MongoDB — `hera.datalayer.document.connectToDatabase`. אומת מקצה-לקצה, **בלי שינוי אף שורת קוד ייצור**. זה מזיז 96 פונקציות מ"דורש שירות חי" ל"טסט הרמטי".
2. **מדידה מדויקת של התלות בנתונים** מראה שרק 134 פונקציות (10%) זקוקות לקבצי נתונים אמיתיים מ-S3 — לא 576 כפי שהוערך בהתחלה.

בנוסף, החקירה חשפה **שני באגי ייצור, שני פגמי תשתית טסטים, וסחיפת גרסאות** — כולם מפורטים בסעיף 6. הם הראיה הטובה ביותר שהמאמץ הזה מייצר ערך ולא רק מספר.

---

## 2. מצב קיים — מדידה, לא הערכה

| מדד | ערך |
|---|---|
| מודולים ב-`hera/` (ללא tests/docs) | 184 |
| קלאסים | 223 |
| פונקציות/מתודות ציבוריות | **1,333** |
| פונקציות פרטיות / dunder | 256 / 191 |
| קבצי טסט | 21 |
| פונקציות `def test_` | 407 |
| שורות קוד טסט | ~7,150 |
| items שנאספים ע"י pytest | 1,397 |
| מדידת coverage | **לא קיימת** — `pytest-cov` לא מותקן, אין `.coveragerc`, אין gate |

**הערה על 1,397 מול 407:** הפער נובע משני קבצים שעושים `parametrize` בהיקף רחב — `test_notebooks.py` (**741** items) ו-`test_no_invalid_escapes.py` (264 items, אחד לכל קובץ מקור). יחד 1,005 מתוך 1,397. מכיוון ש-`test_notebooks.py` מסומן `notebook` ומוחרג מ-CI ומ-`make test`, **הסוויטה האפקטיבית ב-CI היא 656 items**, ומספר טסטי הלוגיקה האמיתי הוא ~407. חשוב לדעת את זה כדי לא לטעות ולחשוב שהסוויטה גדולה פי שלושה ממה שהיא.

---

## 3. המחסום הארכיטקטוני המרכזי

`abstractToolkit` **יורש מ-**`Project` (`hera/toolkit.py:33`). ו-`Project.__init__` (`hera/datalayer/project.py:142`) מבצע, בכל בנייה:

1. בניית ארבעה אובייקטי collection — כל אחד קורא ל-`getDBObject(...)`.
2. **קריאה מה-DB** — `self.getConfig()` כדי לאתר `filesDirectory`.
3. **כתיבה ל-DB** — `self.setConfig(filesDirectory=...)` אם אין ערך שמור.
4. **תופעת לוואי בדיסק** — `os.makedirs` תחת `~/.hera/<projectName>`.

המשמעות: **אי אפשר לבנות שום toolkit בלי DB חי**, וכל בנייה מזהמת גם את ה-DB וגם את ספריית הבית של המשתמש. זו הסיבה ש-`conftest.py` נושא מכונת ניקוי מורכבת (`_no_trace_guard`, `_TEST_PROJECT_PREFIXES`).

הפתרון הקיים בקוד לעקיפת זה הוא העברת `None` כ-`self` למתודות לא-קשורות:

```python
# hera/tests/test_toolkit_coverage.py
TopographyToolkit.convertPointsCRS(None, [(34.78, 32.08)], WSG84, ITM)
topographyAnalysis(None).calculateStastics(elevation_df)
```

זה עובד, אבל **רק למתודות שלא נוגעות ב-`self` בכלל**. הוא לא מרחיב ל-96 הפונקציות שכן משתמשות ב-datalayer, והוא נשבר בשקט ברגע שמתודה מתחילה לגעת ב-`self`.

---

## 4. ה-seam שאומת

`createDBConnection` קורא ל-`connectToDatabase` **כשם גלובלי במודול**. לכן החלפת הפונקציה הזו בזמן טסט מפנה את כל המחסנית ל-`mongomock`:

```python
# hera/tests/unit/conftest.py — קוד טסט בלבד
def _mongomock_connect(mongoConfig, alias=None):
    alias = '%s-alias' % mongoConfig['dbName'] if alias is None else alias
    disconnect(alias)
    return connect(alias=alias, db=mongoConfig['dbName'],
                   mongo_client_class=mongomock.MongoClient,
                   uuidRepresentation='standard')

hera.datalayer.document.connectToDatabase = _mongomock_connect
```

**מדוע `mongo_client_class` ולא `mongomock://`:** הסכימה `mongomock://` הוסרה ב-mongoengine והגרסה המותקנת (0.29.3) זורקת עליה שגיאה מפורשת. `mongo_client_class` נתמך גם בגרסה שמוצמדת ב-`requirements.txt` (0.29.1) וגם בזו שמותקנת בפועל — כלומר הוא הבחירה היחידה שעמידה לסחיפת הגרסאות שתוארה בסעיף 6.5.

### מה אומת בפועל

| שלב | תוצאה |
|---|---|
| `createDBConnection` דרך mongomock | ✅ |
| בניית `Project(...)` | ✅ |
| `addMeasurementsDocument` + שאילתה חוזרת | ✅ |
| `setConfig` / `getConfig` | ✅ |
| `toolkitHome.getToolkit(GIS_RASTER_TOPOGRAPHY)` — toolkit אמיתי עם analysis | ✅ |
| `addDataSource` + `setDataSourceDefaultVersion` + `getDataSourceDocuments` | ✅ |

### משטח השאילתות שאומת מול mongomock

שוויון שדה שטוח · שדה מקונן (`nested__k`) · טווח מספרי (`height__gte`) · regex (`__contains`) · `getProjectList` (distinct) · `getDocumentsAsDict` · אוסף simulations · אוסף cache כולל `getData()` · מחיקה וספירה אחריה · `getConfig` — **כולם עוברים**.

**זה חשוב מכיוון שהטסטים בדרגה B מפעילים את קוד ה-datalayer האמיתי** — בניית שאילתות, סינון `desc`, פתרון גרסאות, תיוג מסמכים. הם לא בודקים mock; הם בודקים את Hera, מול מנוע אחסון בזיכרון.

---

## 5. שלוש דרגות הטסטים

הסיווג הוא **לפי מה שהפונקציה עושה**, לא לפי נוחות:

| דרגה | קריטריון | מספר | קלט | שירותים נדרשים |
|---|---|---|---|---|
| **A — חישוב טהור** | מקבלת נתון ומחשבת | **1,103** (83%) | DataFrame/מערך סינתטי | אין |
| **B — מתווכת ע"י Mongo** | עבודתה היא שאילתה/כתיבה של מסמכים | **96** (7%) | `mongomock` | אין |
| **C — קבצי נתונים אמיתיים** | קוראת/כותבת parquet, HGT, TIF, netcdf | **134** (10%) | נתוני S3 | `TEST_HERA` |

**שכבת ה-`unit` = A + B** (1,199 פונקציות, 90%). שתיהן הרמטיות, בלי שירות חי, ורצות בכל סביבה.
**שכבת ה-`integration` = C**, וגם הרצה חוזרת של B מול MongoDB אמיתי כדי לאמת שההתנהגות בזיכרון תואמת.

זו התשובה לשאלה "אבל אם אין את המידע?" — Mongo ו-S3 חיוניים, אבל לקבוצה מוגדרת וקטנה בהרבה ממה שנראה בהתחלה. לפונקציה שתפקידה לקרוא HGT — צריך HGT אמיתי. לפונקציה שמקבלת DataFrame ומחשבת סטטיסטיקות טורבולנציה — DataFrame סינתטי הוא הקלט הכי **כנה**, כי הוא נותן שליטה על ה-edge cases שנתוני S3 במקרה לא מכילים (NaN, שורה בודדת, טור קבוע, קלט ריק).

---

## 6. באגים ופגמים שנמצאו בחקירה

כולם אומתו בהרצה. הם נכנסים לתוכנית כתיקונים, לא כהערות שוליים.

### 6.1 `dictToMongoQuery` משטח ערכי-רשימה — כל אופרטור רשימה מוחזר ריק בשקט

```
station__in=["A","C"]  →  {'desc__station__in__0': 'A', 'desc__station__in__1': 'C'}
```

השאילתה מחפשת מסמך שבו `desc.station.in.0 == "A"` — חסר משמעות. **התוצאה: 0 מסמכים, בלי שגיאה.** אומת ש-mongomock עצמו תומך ב-`$in` (החזיר 2 מתוך 3), כלומר זה באג ב-`hera/utils/query.py` ולא מגבלת mock. חמור במיוחד מפני שגרסאות נשמרות כרשימות (`version=[0,0,1]`), כך שהגבול בין "עובד במקרה" ל"נכשל בשקט" דק.

**זו גם ההוכחה לכך ש-golden tests לא מתאימים ללב הקוד** — snapshot היה מקפיא את ההתנהגות השגויה הזו כ"נכונה".

### 6.2 `toAzimuthAngle` לא מנרמל קלט — מחזיר זוויות לא חוקיות

```python
toAzimuthAngle(800)   # -> -350
toAzimuthAngle(-400)  # ->  490
```

אזימוט חייב להיות ב-`[0, 360)`. הביטוי (`hera/utils/angle.py:4`) בודק `>= 0` פעם אחת בלי `% 360` מקדים, כך שכל קלט מחוץ למחזור הראשון מייצר פלט לא חוקי.

### 6.3 `toMathematicalAngle` הוא **אותה פונקציה** כמו `toMeteorologicalAngle`

```python
# hera/utils/angle.py:2-3
toMeteorologicalAngle = lambda mathematical_angle: (270 - mathematical_angle) % 360
toMathematicalAngle   = toMeteorologicalAngle   # אימות בהרצה: `met is math_` -> True
```

מתמטית זה תקף — `(270-x) % 360` היא אינוולוציה, ולכן ההמרה זהה בשני הכיוונים. אבל זה לא מתועד, ומשמעותו שהשמות לא נותנים שום הגנה מפני שימוש בכיוון הלא-נכון. `CLAUDE.md` דורש המרה מפורשת דווקא בגלל הסיכון הזה. **נדרש טסט שמקבע את תכונת האינוולוציה** כדי שאם מישהו "יתקן" אחת מהן — זה ייכשל.

### 6.4 `pytest.skip` על נתונים חסרים → CI ירוק כשאין נתונים

`conftest.py` שורות 226, 235, 254, 276 מדלגות כשחסרים `TEST_HERA`, `data_config.json`, ה-result set או `test_repository.json`. אם ההורדה מ-S3 נכשלת חלקית, **76 טסטי האינטגרציה מדולגים בשקט וה-CI מדווח הצלחה**. זה כרגע מסוכן יותר מכל פער כיסוי.

### 6.5 טסטים מסומנים `unit` נגררים ל-Mongo — 60 שניות של timeout

הרצת `pytest -m unit` ארכה 66.76 שניות. הפילוח מראה שהעלות כמעט כולה שני פריטים:

```
30.43s setup     test_toolkit_dynamic_loading_unit.py::test_dynamic_class_loading_via_getdata
30.03s teardown  test_workflow_execution.py::test_custom_target_task
```

הרצת אותו טסט **לבד** לוקחת **0.09 שניות**. 30 השניות הן `serverSelectionTimeoutMS` הדיפולטיבי של pymongo, כי fixture ברמת session מה-conftest הראשי מנסה להתחבר ל-Mongo כבוי. כלומר טסט שמסומן `unit` תלוי בפועל בשירות חי, וזה מתגלה כהמתנה שקטה ולא ככשל.

**המסקנה לתכנון:** העלות האמיתית של טסט unit היא ~0.005–0.09 שניות, ו-`import hera` עולה 0.04 שניות. שכבת unit של ~1,200 טסטים היא ריאלית **מתחת לדקה**, ברגע שהבידוד נאכף.

### 6.6 `getToolkitDocument` בולע כל שגיאה

`hera/toolkit.py` — שלושה בלוקים של `except Exception: pass` בזה אחר זה, ואז `return None`. הפונקציה לא מסוגלת לאותת על כשל: תקלת חיבור, שאילתה שגויה ומסמך שלא קיים כולם מחזירים `None`. צריך טסט שמקבע את ההתנהגות הרצויה, ותיקון שמבדיל בין "לא נמצא" ל"נכשל".

### 6.7 סחיפת גרסאות בין `requirements.txt` לסביבה בפועל

| חבילה | מוצמד ב-requirements | מותקן בפועל |
|---|---|---|
| pytest | 7.2.0 | 9.1.1 |
| mongoengine | 0.29.1 | 0.29.3 |
| pymongo | 4.11.2 | 4.17.0 |

ה-CI מתקין מ-`requirements.txt` על Python 3.11; הסביבה המקומית היא Python 3.12 עם גרסאות חדשות יותר. **טסט יכול לעבור מקומית ולהיכשל ב-CI, או ההפך.** זה גם מה שהופך את `mongo_client_class` (סעיף 4) לבחירה היחידה שעמידה לשתי הגרסאות.

---

## 7. ארכיטקטורת הסוויטה

```
hera/tests/
├── conftest.py              # קיים — session Project מול Mongo אמיתי (integration)
├── unit/
│   ├── conftest.py          # ה-seam של mongomock + guards; אין Mongo, אין S3, אין רשת
│   ├── _stubs.py            # sys.modules stubs ל-7 החבילות החסרות
│   ├── _factories.py        # בוני DataFrame / GeoDataFrame / xarray סינתטיים
│   ├── _unitproject.py      # פקטורי Project/toolkit על גבי mongomock (דרגה B)
│   ├── test_<area>_*.py     # דרגות A ו-B
│   └── golden/              # snapshots — לשוליים בלבד
└── test_*.py                # הקיימים — נשארים integration, לא נוגעים בהם
```

### `hera/tests/unit/conftest.py` — שלושה מנגנוני אכיפה

1. **ה-seam** — החלפת `connectToDatabase` ב-mongomock (סעיף 4), ואיפוס ה-DB בין טסטים.
2. **חסימת רשת** — autouse fixture שמחליף `socket.socket` בכפול שזורק. בלי זה שכבת ה-unit מתדרדרת בשקט לתלות בשירות, בדיוק כמו בסעיף 6.5. **הכשל חייב להיות מיידי ורועש, לא 30 שניות של המתנה.**
3. **חסימת כתיבה לספריית הבית** — `filesDirectory` מוזרם מ-`tmp_path`, ו-fixture נכשל אם `~/.hera` השתנה במהלך הטסט.

### 7 ה-stubs הנדרשים

`PyFoam`, `paraview`, `FreeCAD`, `hermes`, `argos`, `evtk`, `torch` — בסך הכול 19 מודולים בלבד תלויים בהם. הטכניקה קיימת ומוכחת ב-`test_toolkit_coverage.py:26-60`; היא עוברת ל-`_stubs.py` ומיושמת דרך hook שרץ לפני איסוף הטסטים.

---

## 8. מדידת coverage וה-gate

### כלים

`pytest-cov` נכנס ל-`requirements.txt` עם הצמדה מדויקת. `.coveragerc`: `branch = True`, `omit` על `hera/tests/*`, `hera/doc/*`, `hera/bin/*`, ו-`exclude_lines` ל-`if TYPE_CHECKING:` ול-`raise NotImplementedError`.

### מדידה משולבת — לא לפי שכבה

הרצת ה-unit וה-integration מייצרות כל אחת `.coverage.*`, ואחריהן `coverage combine`. ה-gate נבדק **רק על האיחוד**.

**זה לא פרט טכני אלא החלטת תמריצים.** אם ה-gate נמדד רק על שכבת ה-unit, הדרך הזולה להעלות את המספר היא למקק כל גבול — וזה בדיוק מה שהופך coverage למדד ריק. במדידה משולבת אין תגמול על mocking.

### ה-gate

`--cov-fail-under` נקרא מ-`coverage_floor.txt` בשליטת גרסה. כל PR שמעלה כיסוי מעדכן את הקובץ; הרצפה עולה מונוטונית. **אין יעד אחוזי מוצהר במסמך הזה** — הרצפה ההתחלתית נקבעת מהמדידה בפועל באצווה 0, כי יעד שנקבע לפני מדידה הוא ניחוש.

### שני שלבי CI

| שלב | תוכן | תלויות | תקציב זמן |
|---|---|---|---|
| `unit` | `pytest hera/tests/unit -m unit` | אין — בלי service, בלי S3 | **< 60 שניות** |
| `integration` | הסוויטה הקיימת + דרגה C | MongoDB service + נתוני S3 | כמו היום |

שלב ה-unit רץ **ראשון**. כשל בלוגיקה טהורה מתגלה בפחות מדקה, ולא אחרי 15 דקות של הקמת service והורדה מ-S3.

---

## 9. תוכנית האצוות

כל אצווה = ענף + PR עם עליית coverage מוכחת. המספרים מדודים, לא מוערכים.

| # | אזור | A טהור | B mongomock | C נתונים | סה"כ |
|---|---|---|---|---|---|
| **0** | תשתית — ראה פירוט למטה | — | — | — | — |
| 1 | `utils` (כולל `data`, `logging`, `rag`) | 100 | 9 | 12 | 121 |
| 2 | `toolkit.py` + `datalayer` | 77 | 45 | 24 | 146 |
| 3 | `simulations/gaussian` | 114 | 0 | 0 | 114 |
| 4 | `evaporation`+`deposition`+`hydrodynamics`+`windProfile`+`utils`+`analysis` | 91 | 0 | 3 | 94 |
| 5 | `riskassessment` + `presentation` | 109 | 2 | 3 | 114 |
| 6 | `measurements/meteorology` | 136 | 2 | 12 | 150 |
| 7 | `measurements/GIS` | 78 | 8 | 11 | 97 |
| 8 | `measurements/experiment` | 77 | 4 | 5 | 86 |
| 9 | `openFoam`+`LSM`+`mlDL`+`WRF`+`hermesWorkflow`+CLIs | 321 | 26 | 64 | 411 |
| | **סה"כ** | **1,103** | **96** | **134** | **1,333** |

### אצווה 0 — תשתית (חובה, חוסמת את כל השאר)

1. תלויות חדשות ב-`requirements.txt`, מוצמדות: `pytest-cov` ו-`mongomock` (אומת מול 4.3.0). שתיהן נדרשות בשלב ה-unit ב-CI.
2. `.coveragerc` + `coverage_floor.txt` + מדידת בסיס.
3. `hera/tests/unit/` עם ה-seam, שלושת ה-guards, ה-stubs וה-factories.
4. יישור גרסאות (סעיף 6.7) — ראה החלטה מעודכנת בתוכנית הביצוע.
5. `--require-data` / `HERA_REQUIRE_TEST_DATA=1` — הופך את ארבעת ה-`pytest.skip` לכשל קשה, ונכנס ל-CI (סעיף 6.4).
6. תיקון גרירת ה-Mongo מטסטי unit (סעיף 6.5).
7. שני שלבי CI, עם `coverage combine`.

### נימוק לסדר

אצוות 1-5 (589 פונקציות) הן לוגיקה טהורה שבה assertion התנהגותי זול ואמיתי, והן מכילות את הקוד שכל השאר מסתמך עליו — כולל שני הבאגים מסעיפים 6.1 ו-6.2. אצווה 9 היא הגדולה ביותר (411) והכי פחות שווה לפונקציה: 321 מ-A שלה הן עטיפות של OpenFOAM ובניית דיקשנרים. **היא אחרונה בכוונה** — אם ההתקדמות נתקעת שם, תשע האצוות שלפניה עומדות בזכות עצמן.

---

## 10. קונבנציות כתיבת טסטים

### assertions מדורגים

**בלב הקוד — assertion נגזר מהמפרט, לא מהקוד.** `toMeteorologicalAngle(90) == 180` כי מטאורולוגי 0°=צפון ומתמטי 0°=מזרח; המרות `ureg` מול יחסי המרה ידועים; `BriggsRural` מול טבלאות ה-σ המפורסמות; `convertCRS` מול קואורדינטות ITM מתועדות.

**בשוליים — golden.** פלטי VTK, דיקשנרים של OpenFOAM, גיאומטריות מורכבות: קובץ snapshot מתחת ל-`unit/golden/`, שנוצר פעם אחת ונסקר ידנית לפני commit.

### ציפייה מוצהרת: התהליך הזה ימצא באגים

החקירה מצאה שניים לפני שנכתב טסט אחד. **טסט ראשון שמתנגש עם התיעוד מדווח כ-issue — הוא לא "מיושר" בשקט לפי הקוד.** ההיסט הזה הוא הנקודה שבה coverage מייצר ערך במקום מספר.

### כללי היגיינה — נגזרים מהפגמים בסעיף 6

| כלל | מדוע |
|---|---|
| `sorted()` על כל `glob` | סדר קבצים לא דטרמיניסטי |
| `pytest.approx` בכל השוואת float | סחיפה בין גרסאות numpy |
| `Agg` נכפה, איפוס figures | דליפת state של matplotlib |
| `tmp_path` בכל כתיבה לדיסק | אף טסט לא כותב ל-`~/.hera` |
| כל edge case מפורש: ריק, NaN, שורה בודדת, טור קבוע | בדיוק מה שנתוני S3 לא מכילים |
| `-p no:randomly` לא בשימוש; סדר הטסטים לא משנה | טסט שתלוי בסדר הוא באג |

### מבנה ושמות

`Test<Subject>` לפי הקונבנציה ב-`CLAUDE.md`; קובץ אחד לכל אזור-מודול; מסמן `@pytest.mark.unit` על כל דרגה A ו-B; `@pytest.mark.integration` על C.

---

## 11. סיכונים ומיטיגציה

| סיכון | חומרה | מיטיגציה |
|---|---|---|
| **פער התנהגות mongomock ↔ Mongo אמיתי** | גבוהה | 12 סוגי שאילתה אומתו מראש (סעיף 4). דרגה B רצה **גם** מול Mongo אמיתי בשלב ה-integration; פער מתגלה כשוני בין שני השלבים ולא מתגלה מאוחר בייצור. |
| flakiness מיובא (matplotlib, glob, float) | בינונית | כללי ההיגיינה בסעיף 10, נאכפים בסקירת PR |
| שכבת ה-unit מתדרדרת לתלות בשירות | בינונית | חסימת socket — כשל מיידי ורועש. זה מה שהיה חוסך את 60 השניות בסעיף 6.5 |
| זמן CI גדל | נמוכה | נמדד: 0.005-0.09s לטסט. תקציב 60 שניות נבדק בכל PR |
| golden tests מקפיאים באגים | בינונית | golden מוגבל לשוליים; הלב מקבל assertions התנהגותיים. סעיף 6.1 הוא המקרה הממחיש |
| סחיפת גרסאות CI ↔ מקומי | בינונית | יישור ב-אצווה 0; ה-seam בנוי על `mongo_client_class` שעמיד לשתי הגרסאות |

---

## 12. Definition of Done

**לכל אצווה:**
1. שלב ה-unit ב-CI ירוק, מתחת ל-60 שניות.
2. שלב ה-integration ירוק, בלי skips חדשים.
3. `coverage_floor.txt` עלה, וה-PR מציג את הדלתא.
4. כל התנגשות בין קוד לתיעוד דווחה כ-issue ולא יושרה בשקט.

**למאמץ כולו:** לכל אחת מ-1,333 הפונקציות הציבוריות יש כיסוי מדוד; ה-gate מונע נסיגה; ושכבת ה-unit רצה בכל סביבה בלי MongoDB, בלי S3 ובלי רשת.

---

## 13. מחוץ להיקף

- הסוויטות של `ui/client` ו-`ui/server` — יש להן שלבי CI משלהן ותהליך נפרד (`ui/client/TEST_UI.md`).
- `Hermes/` ו-`pyargos` — ריפוזיטוריים חיצוניים.
- טסטי notebooks — נשארים תחת המסמן `notebook`, מחוץ ל-`make test`.
- ריפקטורינג של `Project.__init__` להסרת תופעות הלוואי. **מוצדק בפני עצמו** (סעיף 3), אבל הוא שינוי שובר, וה-seam בסעיף 4 מייתר את הצורך לבצע אותו כתנאי מקדים למאמץ הטסטים.
