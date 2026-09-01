# אצווה 1 — `hera/utils`: תוכנית ביצוע

**ענף:** `tests/batch1-utils` (מסתעף מ-`tests/phase0-infra`, שמכיל את התשתית)
**Spec:** [`../specs/2026-08-24-unit-test-expansion-design.md`](../specs/2026-08-24-unit-test-expansion-design.md)
**קודם:** [`2026-08-24-phase0-test-infrastructure.md`](2026-08-24-phase0-test-infrastructure.md)

**מטרה:** כיסוי התנהגותי ל-14 המודולים הניתנים-לייבוא ב-`hera/utils`, בשכבת ה-unit ההרמטית, בלי לשנות שורת קוד ייצור.

---

## Global Constraints

זהים ל-Phase 0, וחלים על כל משימה:

1. **אפס שינוי בקוד ייצור.** מותר לשנות רק `hera/tests/**` ו-`coverage_floor.txt`.
2. **באגים מדווחים, לא מתוקנים.** כל באג מקבל טסט `xfail(strict=True)` שמאשר את ההתנהגות **הנכונה** — כך שכשמישהו יתקן, הטסט ייכשל וידרוש להסיר את המסמן. הממצאים מרוכזים ל-**Issue אחד** בסוף כל האצוות, לפי החלטת המשתמש.
3. **ענף ייעודי, בלי PR.** בסוף — push בלבד.
4. **שכבת ה-unit נשארת מתחת ל-60 שניות.**
5. **assertion נגזר מהמפרט, לא מהקוד.** אם הקוד והתיעוד מתנגשים — הטסט מאמת את התיעוד ומסומן xfail.

---

## היקף

### בהיקף — 14 מודולים, 912 statements

| מודול | stmts | לא מכוסה | כיסוי נוכחי |
|---|---|---|---|
| `jsonutils.py` | 224 | 172 | 21% |
| `unitHandler.py` | 173 | 95 | 38% |
| `latex.py` | 116 | 116 | 0% |
| `matplotlibCountour.py` | 65 | 62 | 3% |
| `dataframeutils.py` | 54 | 54 | 0% |
| `logging/helpers.py` | 55 | 14 | 70% |
| `zipUtils.py` | 45 | 38 | 9% |
| `query.py` | 40 | 18 | 53% |
| `filter_immediate.py` | 38 | 34 | 6% |
| `statistics.py` | 29 | 6 | 73% |
| `slurm.py` | 28 | 26 | 6% |
| `lazy.py` | 22 | 7 | 62% |
| `utils/__init__.py` | 20 | 4 | 81% |
| `angle.py` | 3 | 0 | 100%\* |
| **סה"כ** | **912** | **646** | **23%** |

\* `angle.py` מדווח 100% כי שלוש הלמבדות מוגדרות בזמן import. אף טסט לא בודק את **ההתנהגות** שלהן — הדגמה חיה לכך ש-coverage לבד הוא מדד חלקי.

### מחוץ להיקף, עם נימוק

| מודול | stmts | מדוע נדחה |
|---|---|---|
| `utils/rag/*` | 410 | `ModuleNotFoundError: httpx`. `chromadb`, `sentence_transformers`, `llama_index` גם חסרים, ו**אף אחד מהם לא ב-`requirements.txt`**. לא ניתן לייבא בהתקנה תקנית. |
| `utils/data/CLI.py` | 480 | dispatch של argparse, כמעט כולו I/O. אצווה נפרדת. |
| `utils/data/toolkit*.py` | 449 | דרגות B ו-C — mongomock ונתוני S3. נכנס לאצווה 2 עם `datalayer`. |
| `SALibUtils.py` | 63 | `ModuleNotFoundError: SALib`, שאינו ב-`requirements.txt`. |
| `freeCAD.py` | 20 | דורש FreeCAD. |

---

## באגים שאומתו לפני כתיבת טסט

כולם אומתו בהרצה. כל אחד מקבל טסט `xfail(strict=True)`.

| # | מודול | הממצא |
|---|---|---|
| B1 | `angle.py:4` | `toAzimuthAngle` לא מנרמל: `az(800) = -350`, `az(-400) = 490`. אזימוט חייב להיות ב-`[0,360)`. |
| B2 | `angle.py:2-3` | `toMeteorologicalAngle is toMathematicalAngle` — אותו אובייקט. תקף מתמטית (אינוולוציה), לא מתועד. |
| B3 | `query.py` | `dictToMongoQuery` משטח רשימות למפתחות ממוספרים. השיטוח **מכוון** לשדות-רשימה (התיעוד אומר זאת), אבל הוא מחיל את אותו כלל על סיומות אופרטור של mongoengine שמקבלות רשימה — `__in`, `__nin`, `__all` — ואז הסיומת נבלעת בנתיב והשאילתה מחזירה אפס בשקט. |
| B4 | `dataframeutils.py` | `return pandas.Dataframe({})` — שגיאת כתיב (`DataFrame`). מסלול "קלט ריק" המתועד זורק `AttributeError`. |
| B5 | `zipUtils.py` | `raise(f"Type Error: ...")` — זריקת מחרוזת. התוצאה `TypeError: exceptions must derive from BaseException`, וההודעה שנכתבה לא מגיעה למשתמש. |
| B6 | `zipUtils.py` | `list_json_files_in_zip` על חבר `.json` בינארי מדפיס ואז נופל ב-`UnboundLocalError` על `content`. |
| B7 | `slurm.py` | `#SBATCH -mem={memoryInGB}G` — הדגל הוא `--mem=`. |
| B8 | `unitHandler.py` | ההערה מצהירה "celsius and K are NOT overridden — pint versions are kept", והקוד מיד אחריה כן דורס אותם ב-unum. `type(celsius).__module__ == "unum"`. |

---

## מבנה הטסטים

קובץ אחד למודול, תחת `hera/tests/unit/`:

```
test_utils_angle.py            test_utils_query.py
test_utils_filter.py           test_utils_lazy.py
test_utils_statistics.py       test_utils_dataframeutils.py
test_utils_slurm.py            test_utils_zip.py
test_utils_latex.py            test_utils_units.py
test_utils_contour.py          test_utils_jsonutils.py
test_utils_package.py          test_utils_logging.py
```

כל טסט מסומן `@pytest.mark.unit`. כל מודול מקבל מחלקה `Test<Subject>` לפי `CLAUDE.md`.

**קונבנציית באגים:**

```python
@pytest.mark.unit
@pytest.mark.xfail(strict=True, reason="B1: toAzimuthAngle does not normalise its input; see the consolidated findings issue")
def test_azimuth_normalises_out_of_range_input():
    """An azimuth must be in [0, 360) for any real input."""
    assert 0 <= toAzimuthAngle(800) < 360
```

---

## המשימות

כל משימה: כתיבת הטסטים → הרצה → commit. הטסטים נכתבים מול המפרט שבדוקסטרינג, ונקראים לפני כתיבה.

- [ ] **Task 1 — `angle.py` + `query.py`** (58 stmts). מקבע את האינוולוציה, את שיטוח הרשימות המכוון, את בריחת `__type`, את `prefixExclude`, ואת כל ענפי הטיפוס ב-`andClause`. שני xfail: B1, B3.
- [ ] **Task 2 — `filter_immediate.py` + `lazy.py` + `statistics.py`** (89 stmts). עשר ה-prepositions, גבולות `outsideInterval`, סמנטיקת `inplace`, פרוקסי הייבוא העצל, שלוש הנורמליזציות ו-`ValueError`.
- [ ] **Task 3 — `dataframeutils.py` + `slurm.py`** (82 stmts). שלושת פורמטי הקלט, `longFormat`, `changeDotToUnderscore`, שגיאות המפתח; ולידציה של `slurm`, ספירת השורות, תוכן ה-SBATCH. שני xfail: B4, B7.
- [ ] **Task 4 — `zipUtils.py` + `latex.py`** (161 stmts). זיפ מקובץ/תיקייה/dict, שגיאות טיפוס, קריאת JSON מזיפ; זיהוי עברית, שמות יוניקוד, המרת bibitem. שני xfail: B5, B6.
- [ ] **Task 5 — `unitHandler.py` + `matplotlibCountour.py`** (238 stmts). `dunam = 1000 m²`, `mmH2O`, `tonumber`/`tounit` בכל ארבעת השילובים pint/unum, המרות טמפרטורה; standardize_polygon ו-toGeopandas. xfail אחד: B8.
- [ ] **Task 6 — `jsonutils.py`** (224 stmts, המודול הגדול). `ConfigurationToJSON` / `JSONToConfiguration` הלוך-חזור, `compareJSONS`, `loadJSON`, `setJSONPath`, `JSONVariations`, `stripConfigurationUnits`.
- [ ] **Task 7 — `utils/__init__.py` + `logging/helpers.py`** (75 stmts). שכל שם ב-`__all__` נפתר דרך ה-`__getattr__` העצל, ש-`_`-names זורקים, ש-`__dir__` שלם; ושמות ה-loggers.
- [ ] **Task 8 — מדידה והעלאת הרצפה.** מדידה משולבת, עדכון `coverage_floor.txt`, אימות תקציב 60 השניות.
- [ ] **Task 9 — push.** אימות שאין שינוי בקוד ייצור, `git push -u origin tests/batch1-utils`. **בלי PR.**

---

## Definition of Done

1. שכבת ה-unit ירוקה ומתחת ל-60 שניות.
2. כיסוי `hera/utils` עלה מ-23% משמעותית; המספר המדויק נקבע במדידה ולא נבטח מראש.
3. `coverage_floor.txt` עלה, וה-commit מציג את הדלתא.
4. אף קובץ ייצור לא נגע.
5. כל באג מסומן `xfail(strict=True)` עם הפניה, וממתין ל-Issue המרוכז.
