# ממצאים מהרחבת סוויטת הטסטים — יומן מתמשך

**מטרת הקובץ:** לרכז את כל הפגמים שנמצאו במהלך אצוות הטסטים, כדי שבסוף כל האצוות ייפתח **Issue אחד מקיף** ולא עשרות נפרדים (החלטת המשתמש).

**כלל:** אין תיקוני קוד ייצור באצוות האלה. כל באג מקובע בטסט `xfail(strict=True)` שמאשר את ההתנהגות **הנכונה** — כשמישהו יתקן, הטסט ייכשל וידרוש להסיר את המסמן ולהפוך אותו ל-assertion אמיתי.

כל ממצא כאן אומת בהרצה בפועל, לא בקריאת קוד.

---

## Phase 0 — תשתית

### P1. `dictToMongoQuery` משטח סיומות אופרטור עם רשימה
**קובץ:** `hera/utils/query.py` · **מקובע ב:** `test_utils_query.py::test_list_valued_operator_suffix_survives`

שיטוח רשימות למפתחות ממוספרים הוא **מכוון ומתועד** לשדות שמאוחסנים כרשימה (`version=[0,0,1]` → `version__0`, `version__1`, `version__2`). הבאג הוא שאותו כלל מוחל גם על סיומות אופרטור של mongoengine שמקבלות רשימה לגיטימית:

```
station__in=["A","C"]  ->  {'desc__station__in__0': 'A', 'desc__station__in__1': 'C'}
```

השאילתה מחפשת מסמך שבו `desc.station.in.0 == "A"` — חסר משמעות. **אפס תוצאות, בלי שגיאה.** אומת ש-mongomock ו-MongoDB תומכים ב-`$in` על שדה מקונן (שאילתת pymongo ישירה החזירה 2 מתוך 3), כלומר הבאג ב-hera ולא במנוע האחסון.

**השפעה אמיתית:** `hera/measurements/meteorology/highfreqdata/analysis/abstractcalculator.py:172,190` מעבירים `params__all=self._AllCalculatedParams`, שהוא רשימה (שורה 58 מאתחלת ל-`[]`, שורה 160 עושה `extend`). חיפוש ה-cache מעולם לא מוצא התאמה, ולכן סטטיסטיקות הטורבולנציה בתדר גבוה מחושבות מחדש בכל קריאה במקום להיקרא מאוסף ה-Cache.

### P2. `toAzimuthAngle` לא מנרמל קלט
**קובץ:** `hera/utils/angle.py:4` · **מקובע ב:** `test_utils_angle.py::test_normalises_input_outside_one_cycle`

אין `% 360` לפני התנאי: `toAzimuthAngle(800) = -350`, `toAzimuthAngle(-400) = 490`. אזימוט חייב להיות ב-`[0,360)`.

### P3. שני שמות, פונקציה אחת
**קובץ:** `hera/utils/angle.py:2-3` · **מקובע ב:** `test_utils_angle.py::TestConversionAliasing`

`toMeteorologicalAngle is toMathematicalAngle` הוא `True`. תקף מתמטית — `(270-x) % 360` היא אינוולוציה — אבל לא מתועד, והשמות לא מספקים הגנה מפני המרה בכיוון הלא-נכון, שהוא הסיכון ש-`CLAUDE.md` מצביע עליו כשהוא דורש המרה מפורשת. **בקשה:** docstring שמצהיר על האינוולוציה. הטסט מקבע אותה כך שתיקון של אחד השמות ייכשל ברעש.

### P4. `getToolkitDocument` בולע כל שגיאה
**קובץ:** `hera/toolkit.py:358,367,374`

שלושה `except Exception: pass` רצופים ואז `return None`. תקלת חיבור, שאילתה שגויה ומסמך חסר בלתי-מבחינים מבחוץ. (בקובץ 12 מופעים של `except Exception:` בסך הכול.) **בקשה:** להבדיל בין "לא נמצא" ל"החיפוש נכשל".

### P5. שישה קבצים תחת `hera/` אינם Python תקין
**מקובע ב:** `.coveragerc` (omit, עם נימוק)

לא עוברים `ast.parse`, כלומר אי אפשר לייבא אותם, ו-`coverage report` קורס עליהם (`Couldn't parse ... as Python source`):

```
hera/simulations/gaussian/source.py:12                            invalid syntax
hera/simulations/utils/interpolation/interpolations.py:137        tabs/spaces
hera/simulations/openFoam/toberewritten/netcdf2of.py:11           unexpected indent
hera/simulations/openFoam/toberewritten/xarrayDataset2OF.py:285   bad dedent
hera/simulations/openFoam/eulerian/NavierStokes.old/postprocess/gui/ml.py:1161
hera/simulations/openFoam/eulerian/NavierStokes.old/postprocess/gui/plotter.py:1161
```

הם מועמדים למחיקה — קוד מת מוגמר, בתיקיות ששמן `toberewritten` ו-`NavierStokes.old`.

### P6. שני כשלים קיימים ב-`TestAutoCache`
**קובץ:** `hera/datalayer/autocache.py` + `hera/tests/test_datalayer.py`

```
ImportError: cannot import name 'cacheFunction' from partially initialized
            module 'hera.datalayer.autocache' (most likely due to a circular import)
TypeError:   tuple indices must be integers or slices, not str
```

**אומת על worktree נקי של `origin/master` (4d44428f), בשתי גרסאות pytest (7.2.0 ו-9.0.3) — זהים בדיוק.** לא נגרמו על ידי עבודת הטסטים. משמעות מעשית: שלב אכיפת ה-coverage ב-CI יושב אחרי פקודת ה-pytest של האינטגרציה, ולכן **הרצפה לא תיאכף עד שאלה יתוקנו**, וה-job `test` אדום גם עכשיו על master.

### P7. `pytest.skip` על נתונים חסרים → CI ירוק
**קובץ:** `hera/tests/conftest.py:261,270,289,311` · **תוקן בתשתית**

ארבע נקודות דילגו כשחסרים `TEST_HERA`, `data_config.json`, ה-result set או `test_repository.json`. הורדה חלקית מ-S3 הייתה מדלגת על עשרות טסטים בשקט ומדווחת הצלחה. נפתר ב-`--require-data` / `HERA_REQUIRE_TEST_DATA=1`, חמוש ב-CI בלבד; הרצת מפתח מקומית ממשיכה לדלג.

### P8. סחיפת גרסאות בין `requirements.txt` לסביבה
| חבילה | מוצמד | מותקן בפועל |
|---|---|---|
| pytest | 7.2.0 | 9.0.3 |
| mongoengine | 0.29.1 | 0.29.3 |
| pymongo | 4.11.2 | 4.17.0 |

ה-CI מתקין מ-`requirements.txt` על Python 3.11; הסביבה המקומית Python 3.12. שכבת ה-unit אומתה על **שתי** גרסאות ה-pytest כדי לסגור את הסיכון. ה-seam בנוי על `mongo_client_class` דווקא כי `mongomock://` הוסר ב-mongoengine 0.29.3 בעוד `mongo_client_class` נתמך בשתיהן.

---

## אצווה 1 — `hera/utils`

### B9. `bibItem.convert` — ענף ה-math mode הוא קוד מת
**קובץ:** `hera/utils/latex.py` · **מקובע ב:** `test_utils_latex.py::test_math_mode_is_preserved`

הענף בודק `word[indx] == 'DOLLAR SIGN'` — משווה **תו בודד** למחרוזת של 11 תווים, ולכן לא יכול לפעול. `$` נופל לענף הליטרלי והביטוי מתקלקל. **הדוגמה מהדוקסטרינג של המחלקה עצמה:**

```
קלט:  האם $4\frac{m}{s}$ נכון
פלט:  $$4\$\L{frac}{\L{m}}{\L{s}}$
```

### B10. `bibtexFile.items` מחזיר פריט ריק ראשון תמיד
**קובץ:** `hera/utils/latex.py:_parseBibItems` · **מקובע ב:** `test_utils_latex.py::test_item_count_matches_the_bibitem_count`

הלולאה מוסיפה את המצבר בכל שורת `\bibitem`, כולל בפעם הראשונה כשהוא עוד `[]`. שני פריטים אמיתיים → `len(items) == 3`.

### B11. קובץ bib ריק זורק `IndexError`
**קובץ:** `hera/utils/latex.py:__init__` · **מקובע ב:** `test_utils_latex.py::test_empty_file_reports_clearly`

`self._first_last_Lines = (thefileStr[0], thefileStr[-1])` על קובץ ריק → `IndexError: list index out of range` במקום שגיאה שמסבירה שהביבליוגרפיה ריקה.

### B5. `zip_items` זורק מחרוזת
**קובץ:** `hera/utils/zipUtils.py` · **מקובע ב:** `test_utils_zip.py::test_a_path_that_is_neither_file_nor_directory_names_itself`

`raise(f"Type Error: The item {item} is not a file or a directory")` — זריקת מחרוזת. Python מדווח `TypeError: exceptions must derive from BaseException`, וההודעה שמזהה את הנתיב הבעייתי **לא מגיעה למשתמש**.

### B6. `list_json_files_in_zip` נופל על חבר בינארי
**קובץ:** `hera/utils/zipUtils.py` · **מקובע ב:** `test_utils_zip.py::test_an_undecodable_json_member_is_skipped`

מטפל ב-`UnicodeDecodeError` בהדפסה `"Binary file, skipped reading content."` ואז ממשיך לשורה שמשתמשת ב-`content` — שלא הוגדר. `UnboundLocalError`. הכוונה המוצהרת היא לדלג.

### B7. דגל הזיכרון של Slurm שגוי
**קובץ:** `hera/utils/slurm.py` · **מקובע ב:** `test_utils_slurm.py::test_memory_directive_uses_slurms_option_name`

`#SBATCH -mem={memoryInGB}G` — הדגל של Slurm הוא `--mem=`. מקף אחד גורם ל-sbatch לדחות את הסקריפט.

### B8. הערה שסותרת את הקוד שמתחתיה
**קובץ:** `hera/utils/unitHandler.py` · **מקובע ב:** `test_utils_units.py::test_constants_are_the_pint_units_the_comment_promises`

ההערה: `"celsius and K are NOT overridden — pint versions are kept"`. הקוד מיד אחריה דורס את שניהם ב-unum. `type(celsius).__module__ == "unum"`.

### B12. `tounit` מייצר Quantity ב-registry הלא-נכון — **החמור באצווה**
**קובץ:** `hera/utils/unitHandler.py` · **מקובע ב:** `test_utils_units.py::test_result_belongs_to_heras_registry`, `::test_custom_units_are_reachable_through_tounit`

`tounit` — פונקציה ציבורית ב-`__all__` — בונה `Quantity(x, unit)` מהייבוא הישיר של pint (`from pint import Quantity`) ולא מ-`ureg.Quantity`. התוצאה שייכת ל-application registry הדיפולטיבי של pint. שתי השלכות, שתיהן אומתו:

```
tounit(5,"m") + 5*ureg.m   -> ValueError: Cannot operate with Quantity
                                and Quantity of different registries
tounit(1,"dunam")          -> UndefinedUnitError: 'dunam' is not defined
```

`dunam` היא היחידה המותאמת שכל הפואנטה של ה-registry של hera היא להגדיר, והיא בדיוק זו שה-helper הציבורי לא מצליח לייצר. `mmH2O` עובד רק במקרה — ל-pint יש `millimeter_H2O` מובנה משלו.

### B13. `returnStandardize` לא יכול לתקנן לבד
**קובץ:** `hera/utils/jsonutils.py` · **מקובע ב:** `test_utils_jsonutils.py::test_decoder_can_standardize_on_its_own`

מתועד כ-"If true, return the units in MKS". בפועל המפענח קורא את שדה `units`, שמכיל יחידות בסיס **רק אם המקודד קיבל `standardize=True`**. בקשת MKS על מסמך שלא קודד כך מחזירה את היחידות המקוריות בשקט:

```
enc = ConfigurationToJSON({"v": 5*ureg.km}, splitUnits=True)      # standardize=False
JSONToConfiguration(enc, returnStandardize=True)  ->  5 kilometer   # ביקשנו MKS
```

בנוסף, שדה בשם `units` מכיל מחרוזת כמות שלמה (`'5 kilometer'`) ולא יחידה — שם מטעה.

### B14. אנוטציה שגויה ב-`with_logger`
**קובץ:** `hera/utils/logging/helpers.py:82` · **מקובע ב:** `test_utils_package.py::test_the_return_annotation_matches_the_arity`

מסומן `-> (str, dict)` ומחזיר `(logger_name, logger_dict, 'loggers')` — tuple של 3. פריקה לפי האנוטציה זורקת `ValueError: too many values to unpack`. הגרסה הראשונה של הטסט נכשלה בדיוק כך.

### B15. `dataframeutils` — מסלול הקלט הריק שבור
**קובץ:** `hera/utils/dataframeutils.py` · **מקובע ב:** `test_utils_dataframeutils.py::test_empty_input_returns_an_empty_frame`

`return pandas.Dataframe({})` — `DataFrame` באות D גדולה. הקוד מתעד ומתעד ב-log `"Empty input, returning {}"` ואז זורק `AttributeError: module 'pandas' has no attribute 'Dataframe'`.

### B16. מודולים עם תלויות לא-מוצהרות שאינם ניתנים לייבוא
**קבצים:** `hera/utils/SALibUtils.py`, `hera/utils/rag/*`

```
hera.utils.SALibUtils   -> ModuleNotFoundError: No module named 'SALib'
hera.utils.rag.indexer  -> ModuleNotFoundError: No module named 'httpx'
hera.utils.rag.search   -> ModuleNotFoundError: No module named 'httpx'
```

`SALib`, `httpx`, `chromadb`, `sentence_transformers` ו-`llama_index` **אינם מוצהרים ב-`requirements.txt`**. זה לא פער כיסוי אלא 473 statements של קוד שבור בהתקנה תקנית.
