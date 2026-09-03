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

**השפעה אמיתית — תוקן באצווה 6 אחרי בדיקה בהרצה.** הניסוח המקורי כאן טען ש"חיפוש ה-cache מעולם לא מוצא התאמה". **זה היה חזק מדי.** הבדיקה בפועל מראה ש-mongoengine מזהה את `all` כאופרטור גם כשאינו בסוף המפתח, ומוציא אותו:

```
params__all=["u","v"]  ->  desc__params__all__0:"u", desc__params__all__1:"v"
                       ->  {"desc.params.0": {"$all": ["u"]},
                            "desc.params.1": {"$all": ["v"]}}
```

כלומר `$all` — שמוגדר כהכלה **בלתי-תלוית-סדר** — מתדרדר ל**התאמת קידומת לפי מיקום**. מסמך שמאוחסן עם `params=["u","v","w"]`:

| שאילתה | תוצאה |
|---|---|
| `params__all=["u"]` | ✅ מתאים |
| `params__all=["u","v"]` | ✅ מתאים |
| `params__all=["u","v","w"]` | ✅ מתאים |
| `params__all=["v"]` | ❌ **מחטיא** — למרות ש-`v` קיים |
| `params__all=["w"]` | ❌ **מחטיא** |
| `params__all=["v","u"]` | ❌ **מחטיא** — רק הסדר שונה |

ההשלכה על `abstractcalculator.py:172,190`: `_AllCalculatedParams` נבנה ב-`extend()` לפי סדר הבקשות, ולכן ה-cache **פוגע רק כשרשימת הפרמטרים המבוקשת היא קידומת של המאוחסנת באותו סדר**. התנהגות חלקית ותלוית-סדר, לא כישלון מוחלט — וההחטאה שקטה לגמרי.

**מקובע ב:** `test_meteorology_cache_query.py` (שני הכיוונים: מה שכן מתאים ומה שמחטיא).

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

---

## אצווה 2 — `toolkit.py` + `datalayer`

### B18. `getDataSourceList` מחזיר כפילויות
**קובץ:** `hera/toolkit.py:186` · **מקובע ב:** `test_toolkit_datasources.py::test_list_reports_each_name_once`

מתועד כ-"Return a list of data source **names** for this toolkit", אבל מחזיר רשומה אחת לכל **מסמך**. מקור עם שלוש גרסאות מופיע שלוש פעמים:

```
['SRC', 'LSRC', 'MULTI', 'MULTI', 'MULTI', 'MULTI']
```

האסימטריה מול `getDataSourceMap`, שמתועד במפורש כ-"all data sources **and their versions**", מלמדת שהרשימה אמורה להיות שמות ייחודיים.

### B19. עטיפות הוספת המסמכים משנות את המילון של הקורא — **החמור באצווה**
**קובץ:** `hera/toolkit.py:146,158,170` · **מקובע ב:** `test_toolkit_documents.py::TestCallerDictIsNotMutated`

שלוש העטיפות קוראות `desc.setdefault(TOOLKIT_TOOLKITNAME_FIELD, self.toolkitName)` על המילון שהתקבל מהקורא. שתי השלכות, שתיהן אומתו בהרצה:

```python
metadata = {"station": "A"}
topography.addCacheDocument(resource="/topo", ..., desc=metadata)
# metadata == {'station': 'A', 'toolkit': 'TopographyToolkit'}   <-- שונה בלי בקשה

landcover.addCacheDocument(resource="/land", ..., desc=metadata)
# setdefault לא דורס, ולכן המסמך של landcover נשמר עם 'TopographyToolkit'
```

**זה שיבוש נתונים חוצה-toolkits:** המסמך של LandCover נעלם מהשאילתות של LandCover ומופיע בשאילתות של Topography. התיוג הזה הוא בדיוק המנגנון שמפריד בין שני toolkits שחולקים פרויקט, ולכן הפגם פוגע בהפרדה עצמה. תרחיש הטריגר תמים לחלוטין — שימוש חוזר במילון מטא-דאטה.

### B20. מפת הפורמטים ממופתחת על נתיבי מודול פנימיים של צד-שלישי
**קובץ:** `hera/datalayer/datahandler.py:typeDatatypeMap` · **מקובע ב:** `test_datalayer_datahandler.py::test_a_dataframe_is_stored_as_parquet`, `::test_a_series_is_stored_as_json`

המפה ממופתחת על נתיבים כמו `'pandas.core.frame.DataFrame'`. pandas 3 מדווח `'pandas.DataFrame'`, שהמפה לא מכירה, ולכן היא נופלת למפתח `"object"` → **כל DataFrame נשמר כ-pickle במקום parquet**, בלי שגיאה. זה נוגד את המדיניות המוצהרת ב-`CLAUDE.md` ("Tabular data → parquet", "Prefer parquet over CSV"), ו-pickle גרוע יותר: לא בר-החלפה, שביר בין גרסאות, ובעייתי אבטחתית.

```
DataFrame  computed='pandas.DataFrame'                        map key='pandas.core.frame.DataFrame'
Series     computed='pandas.Series'                           map key='pandas.core.series.Series'
dask DF    computed='dask.dataframe.dask_expr._collection.DataFrame'
                                                              map key='dask_expr._collection.DataFrame'
```

`requirements.txt` מצמיד `pandas==2.2.3`, שבו המפתח **כן** מתאים, ולכן **ה-CI אינו חשוף** וזו הסביבה המקומית בלבד. זהו מופע קונקרטי של סחיפת הגרסאות ב-P8. הפגם השורשי אינו "המפה שגויה" אלא ש**היא נשענת על פרטי מימוש פנימיים של תלויות, כך ששדרוג מוריד בשקט את פורמט האחסון**. שני הטסטים נושאים `xfail` עם `condition` שמחושב בזמן איסוף מ-pandas המותקן, כדי שהציפייה תהיה נכונה בשתי הסביבות במקום מתנדנדת באחת.

---

## חשדות שנבדקו ונדחו

לתיעוד, כדי שלא ייבדקו שוב:

- **`setConfig` עם `raise` חשוף** — נראה כאילו מסלול השגיאה של הפרויקט הדיפולטיבי מסתיים ב-`raise` בלי חריגה פעילה, מה שהיה מייצר `RuntimeError` אטום במקום ההודעה המפורטת. **לא נכון:** `project.py:247` הוא `raise ValueError(err)`, וההודעה כן מוצגת. אומת בהרצה.
- **קונסטנטת פורמט בלי handler** — כל 18 הקונסטנטות ב-`datatypes` נבדקו, ולכולן יש `DataHandler_*` עם `getData` ו-`saveData`. אין פגם.
- **סמנטיקת המונים** — הנחתי post-increment וטסטים נכשלו. **הקוד צודק:** הקריאה הראשונה מגדירה את המונה ל-0 בלי להוסיף, ורק הקריאות הבאות מוסיפות. מתועד ב-docstring, והטסטים מקבעים זאת במפורש כי זו התנהגות מפתיעה.

### B21. `_get_data_toolkit` מצהיר על `projectName` ומתעלם ממנו
**קובץ:** `hera/toolkit.py:697` · **מקובע ב:** `test_toolkit_home.py::test_the_helper_ignores_the_project_it_is_given`, `::test_defaults_are_per_project`

```python
def _get_data_toolkit(self, projectName: str = None):
    from hera.utils.data.toolkit import dataToolkit
    return dataToolkit()          # projectName לא בשימוש
```

`dataToolkit.__init__` מקבל `connectionName` בלבד ועובד תמיד על הפרויקט הדיפולטיבי. לכן `setDefaultRepository` / `getDefaultRepository`, שמתועדים כ-"Persist default repository name **for a project**", אינם מוגבלים לפרויקט בכלל — ההגדרה גלובלית.

### B22. `setDefaultRepository` תמיד זורק — הפיצ'ר מת
**קובץ:** `hera/toolkit.py:1164` · **מקובע ב:** `test_toolkit_home.py::test_round_trip`

מכיוון ש-B21 מפנה את הכתיבה ל-`defaultProject`, ו-`Project` שומר עליו כקריא-בלבד (`project.py:757`, קוד ייצור ללא תנאי), כל קריאה נכשלת:

```
RuntimeError: project defaultProject is read-only.
```

`registerToolkit` מציב `dtk._allowWritingToDefaultProject = True` לפני הכתיבה (`toolkit.py:1271`); `setDefaultRepository` לא. **הפיצ'ר של ברירת מחדל למאגר אינו פעיל כלל.**

### B23. חיפוש פורמט מת ב-`setDefaultRepository`
**קובץ:** `hera/toolkit.py:1164` · **מקובע ב:** `test_toolkit_home.py::test_the_json_format_lookup_is_dead`

```python
dfmt = getattr(_dt, "JSON", None) or getattr(_dt, "json", None) or getattr(_dt, "TEXT", None)
```

אף אחד משלושת השמות לא קיים ב-`datatypes` — הקונסטנטות הן `JSON_DICT`, `JSON_PANDAS`, `JSON_GEOPANDAS`. הבלוק כולו מת ו-`dataFormat` לעולם לא נקבע. מדווח כקוד מת ולא כפלט שגוי: המידע יושב ב-`desc`, ולכן הפורמט אינו נושא משמעות כאן. הטסט נכשל אם אחד השמות **כן** יתווסף בעתיד, כדי שהחיפוש ייבדק מחדש.

---

## אצווה 3 — `simulations/gaussian`

### B24. `getSigma` לא מוודא את המרחק
**קובץ:** `hera/simulations/gaussian/Sigma.py` · **מקובע ב:** `test_gaussian_sigma.py::test_an_upwind_distance_is_rejected`

```
x = -100,   stability D  ->  sigmaY = -8.04 m      (סטיית תקן שלילית)
x = -20000, stability D  ->  sigmaZ = nan          (+ numpy RuntimeWarning)
```

הסיבה ל-NaN: `(1 + b·x)` הופך שלילי, ומועלה בחזקה שברית (`C = -0.5`). אין ולידציה, ואין שגיאה — קלט שגוי מתפשט לשדה ריכוזים.

### B25. תחילית `r` על המרכאות הסוגרות
**קבצים:** `hera/simulations/gaussian/Meteorology.py:197,216,246,265` · `gasCloud.py:648` · **מקובע ב:** `test_gaussian_meteorology.py::test_compiling_the_module_emits_no_warning`

```python
def getAirDensity(self, height):
    """
        ... \rho_{air} ... [g\cdot cm^{-3}]
    r"""        # <-- ה-r על המרכאות הסוגרות
```

הכוונה הייתה `r"""` בפתיחה. בפועל ה-`r` הופך לתו האחרון בתוך דוקסטרינג שאינו raw, ולכן ה-backslashes של ה-LaTeX הם escape sequences לא חוקיים. שתי תוצאות: `SyntaxWarning` בכל ייבוא תחת Python 3.12+, ו-`r` תועה בסוף כל דוקסטרינג מושפע.

### B26. `test_no_invalid_escapes.py` אינו מסוגל להיכשל על escape לא-חוקי
**קובץ:** `hera/tests/test_no_invalid_escapes.py:27` · **מקובע ב:** `test_gaussian_meteorology.py::TestRawStringPrefixOnTheWrongQuotes`

זהו הממצא המעניין באצווה. הטסט — עם ~990 מקרים מפורמטרים — מקמפל כל קובץ עם `simplefilter("error", SyntaxWarning)`. אבל **פייתון מציג SyntaxWarning שהוסלם לשגיאה כ-`SyntaxError` בזמן קימפול**, וה-handler שנועד לקבצי legacy שאינם מתקמפלים בולע אותו:

```python
except SyntaxError as err:
    pytest.skip(f"pre-existing SyntaxError, file is not importable: ...")
```

התוצאה, כפי שנצפתה בהרצה:

```
SKIPPED [1] test_no_invalid_escapes.py:27: pre-existing SyntaxError,
            file is not importable: invalid escape sequence '\c' (line 201)
```

**ההודעה שגויה** — `Meteorology.py` כן ניתן לייבוא (אומת). וחמור מכך: **הכשל שהטסט קיים כדי לגלות יכול רק להפוך ל-skip.** ככל שבעיית ה-escape חמורה יותר, כך גדל הסיכוי שתסווג כ"קובץ legacy מחוץ להיקף".

**התיקון המוצע:** לקמפל פעמיים — פעם אחת בלי מסנן, כדי לקבוע אם הקובץ בכלל מתפרסר, ורק אם כן — פעם שנייה עם המסנן, כדי לבדוק escapes. שני המצבים בלתי-מבחינים כרגע. **לא תיקנתי** משום שזהו טסט עובד קיים, וההנחיה היא לדווח ולא לתקן.

---

## חשדות שנבדקו ונדחו — אצווה 3

- **מקדמי Briggs שגויים** — 36 השוואות (6 מחלקות יציבות × 6 מרחקים) מול נוסחאות Briggs (1973) המפורסמות: **תואם עד `rel=1e-9`**. הטבלאות נבדקו גם רשומה-רשומה. אין פגם.
- **בניית המקור הווירטואלי** — חשדתי שאולי לא עקבית. **אומת:** עם `sigma0=10 m`, `getSigma(x=0)` מחזיר בדיוק 10.0 מ'. זו התכונה שמצדיקה את כל טריק הזזת הראשית, והיא מדויקת.
- **`getWindVelocity` וחספוס** — הנחתי ששטח מחוספס **מאט** את הרוח בגובה, וה-assertion נכשל. **הקוד צודק:** עם `u10` קבוע, `ln(z/z0)/ln(10/z0)` **גדל** עם z0 עבור z>10, ולכן חספוס גדול יותר משמעו גרדיינט תלול יותר ורוח חזקה יותר מעל. z0=1m נותן 10.0 m/s ב-100 מ' מול 6.67 עבור z0=1cm.
- **`getTKE` מתעלם מהגובה** — נראה כמו באג, אבל מתועד במפורש: "Currently implemented only neutral conditions". מקובע כטסט עובר כדי שהוספת תלות בגובה תהיה שינוי מכוון.

### B27. `MeshUtils` — שם עמודת sigmaY מוגדר לעמודת ה-X
**קובץ:** `hera/simulations/gaussian/MeshUtils.py:55` · **מקובע ב:** `test_gaussian_mesh.py::TestSigmaColumnNames`

```python
self.sigmaXName = "sigmaXCorrected"
self.sigmaYName = "sigmaXCorrected"      # <-- אותו שם
```

`_gaussianToMesh` קורא `sigy = gaussian[self.sigmaYName]`, ולכן הפרישה הרוחבית קוראת תמיד את הסיגמה האורכית. **אומת מספרית:** קלט `sigmaYCorrected = 50` מייצר σ_y משוחזר של **5.000** — נתוני ה-Y מתעלמים בשקט והענן נכפה איזוטרופי. במודל פיזור אטמוספרי σ_x ו-σ_y שונים בתכלית, ולכן זו שגיאה פיזיקלית ולא רק אי-דיוק.

הגלעין עצמו תקין: טסט עובר מראה שכשמכוונים `sigmaYName` לעמודה שלו, σ_y המשוחזר עוקב אחריה בדיוק. **רק ברירת המחדל בבנאי חוסמת ענן אליפטי.**

### B28. `_defineCoordinates` מבלבל תווית עם מיקום
**קובץ:** `hera/simulations/gaussian/MeshUtils.py:76-84` · **מקובע ב:** `test_gaussian_mesh.py::TestDataFrameIndexing`

```python
maxXID = gaussians[gaussians['x']==gaussians['x'].max()].index[0]   # תווית
maxX   = gaussians.iloc[maxXID]['x'] + ...                          # מיקום
```

`.index[0]` מחזיר תווית ו-`.iloc` מצפה למיקום. כל DataFrame שהאינדקס שלו אינו RangeIndex מאופס זורק `IndexError: single positional indexer is out-of-bounds`. זה מכסה **כל** פלט של `groupby`, `filter` או `concat` — כלומר תרחישי השימוש הרגילים.

### B29. `GaussianIntegrationToMesh` מצמצם את הקלט של ההורה
**קובץ:** `hera/simulations/gaussian/MeshUtils.py:201` מול `:157` · **מקובע ב:** `test_gaussian_mesh.py::test_the_subclass_accepts_a_plain_quantity_like_its_parent`

```python
# GaussianToMesh:            ret = fullX*fullY*gaussian[QName] * (1*ureg.kg).m_as(unumToPint(QUnits))
# GaussianIntegrationToMesh: ret = fullX*fullY*unumToPint(gaussian[QName]).m_as(unumToPint(QUnits))
```

הבסיס מכפיל במספר וממיר גורם יחידות; תת-המחלקה מפעילה `unumToPint` על **הערך עצמו**, ולכן float רגיל הופך ל-dimensionless ו-`.m_as(kg)` זורק. אומת:

| | `Q = 1.0` | `Q = 1.0 * ureg.kg` |
|---|---|---|
| `GaussianToMesh` | ✅ אינטגרל 1.0 | ✅ אינטגרל 1.0 |
| `GaussianIntegrationToMesh` | ❌ `DimensionalityError` | ✅ אינטגרל 1.0 |

הפרת Liskov: קוד שנכתב מול מחלקת הבסיס נשבר בהחלפה שהדוקסטרינג של תת-המחלקה עצמה מזמין ("This is more accurate than estimating the gaussian").

### B30. שתי נקודות כניסה, שתי ברירות מחדל שונות לפרופיל הרוח
**קבצים:** `hera/simulations/gaussian/toolkit.py:62` מול `Meteorology.py:343` · **מקובע ב:** `test_gaussian_toolkit.py::test_the_two_defaults_differ_between_toolkit_and_factory`

`gaussianToolkit.getMeteorologyFromU10` מוגדר `verticalProfileType="powerLaw"`, ו-`MeteorologyFactory.getMeteorologyFromU10` מוגדר `"log"`. אותה קריאה, אותם ארגומנטים — שני פרופילי רוח שונים לגמרי, בהתאם לנקודת הכניסה. לא קריסה, אבל מלכודת שקטה.

---

## אצווה 4 — `evaporation` · `deposition` · `hydrodynamics` · `windProfile` · `simulations/utils`

### B31. `evaporation` ו-`deposition` אינם ניתנים לייבוא — ייבוא יחסי ברמה שגויה
**קבצים:** `evaporation/models.py:3` · `evaporation/monaghan.py:5` · `deposition/models.py:2` · **מקובע ב:** `test_simulations_module_imports.py::TestTheBrokenRelativeImport`

```python
from ..utils import tonumber, tounit     # שתי נקודות = hera.simulations.utils
```

`hera/simulations/utils/__init__.py` **ריק**; `tonumber` ו-`tounit` יושבים ב-`hera.utils`, שלוש רמות למעלה. התוצאה:

```
ImportError: cannot import name 'tonumber' from 'hera.simulations.utils'
```

**שתי חבילות פיזיקה שלמות מתות בכל התקנה.** באותו עץ, `gaussian/DropletCloud.py:5` עושה `from ...utils import tounit,tonumber` — שלוש נקודות — ועובד. גם שורה 4 ב-`evaporation/models.py` עצמו עושה `from hera.utils.unitHandler import ...` נכון. כלומר התיקון הוא נקודה אחת בשלושה קבצים.

זהו הפגם הבודד עם ההשפעה הרחבה ביותר שנמצא עד כה: הוא הפיל כמחצית מהיקף האצווה המתוכנן.

### B32. `calculateR` מערבב נורמליזציות — מחזיר `(N-1)/N` מהמתאם
**קובץ:** `hera/simulations/analysis/errorCalculation.py` · **מקובע ב:** `test_simulations_errorcalculation.py::test_correlation_is_one`, `::test_the_shortfall_is_exactly_the_sample_size_factor`

```python
R = ((model-model.mean())*(measure-measure.mean())).mean() / (model.std()*measure.std())
```

המונה מחלק ב-**N** (`.mean()`); המכנה משתמש ב-`.std()`, שברירת המחדל ב-pandas היא סטיית תקן **מדגמית** (ddof=1, חלוקה ב-N−1). כל תוצאה קטנה בפקטור `(N-1)/N` בדיוק, ואומת על חמישה גדלי מדגם:

| N | `calculateR` למודל מושלם | שגיאה |
|---|---|---|
| 3 | 0.6667 | **33%** |
| 5 | 0.8000 | 20% |
| 10 | 0.9000 | 10% |
| 50 | 0.9800 | 2% |
| 200 | 0.9950 | 0.5% |

מודל מושלם **לעולם לא יקבל 1.0**, וההטיה נעלמת רק אסימפטוטית. במדדי הערכת מודלים (Chang & Hanna) שמדווחים על מדגמים קטנים זו הטיה מהותית. התיקון: `std(ddof=0)`.

### B33. `skin_friction` תמיד זורק — `numpy.log()` בלי ארגומנט
**קובץ:** `hera/simulations/hydrodynamics/nearWallFlow.py:174` ו-`:298` · **מקובע ב:** `test_simulations_nearwallflow.py::TestSkinFriction`

```python
Lambda = 2*numpy.log()
```

`TypeError: log() takes from 1 to 2 positional arguments but 0 were given`. הבאג קיים **בשתי המחלקות** — `channelFlow` ו-`couetteFlow`. שתי המתודות לא יכלו לרוץ מעולם.

### B34. `ReynoldsUm` מחזיר מספר ריינולדס עם יחידות
**קובץ:** `hera/simulations/hydrodynamics/nearWallFlow.py:59,63,97` · **מקובע ב:** `test_simulations_nearwallflow.py::test_a_reynolds_number_is_dimensionless`

```python
@property
def height(self):
    return self._channelHeight.m_as(ureg.m)      # מפשיט יחידות
```

`hydraulicHeight = 2*height` יוצא חסר-ממד, ולכן `ustar*hydraulicHeight/nu` = `(m/s × 1)/(m²/s)` = **1/מטר**. מספר ריינולדס חייב להיות חסר-ממד. `Re_tau` באותה מחלקה נמנע מזה כי הוא משתמש ב-`_channelHeight` (ה-Quantity) ישירות — כלומר שתי גישות סותרות באותו קובץ.

### B35. `ReynoldsUm` מתעלם מהפרמטר שעל שמו הוא קרוי
**קובץ:** `hera/simulations/hydrodynamics/nearWallFlow.py:97` · **מקובע ב:** `test_simulations_nearwallflow.py::test_the_bulk_velocity_changes_the_answer`

הדוקסטרינג: "Return the reynolds based on the **Um** and hydraulic diameter". הגוף מחשב `ustar*hydraulicHeight/nu` ואינו קורא את `Um` כלל. זהו מספר ריינולדס של חיכוך, לא של הזרימה הממוצעת, ושינוי `Um` אינו משנה דבר.

---

## אצווה 5 — `riskassessment` + `presentation`

### B36. ברירת המחדל של יחידות ל-pandas אינה ניתנת להשגה
**קובץ:** `hera/riskassessment/agents/effects/Calculator.py` · **מקובע ב:** `test_risk_calculator.py::TestPandasInput::test_the_documented_default_units_apply`

```python
inUnits = concentrationField.attrs[field] if hasattr(concentrationField, "attrs") else 1*(ureg.mg/ureg.m**3)
```

הבדיקה `hasattr(..., "attrs")` נכתבה כשרק ל-xarray היה `.attrs`. **pandas מגרסה 1.0 מספק `.attrs` גם הוא**, ולכן DataFrame נכנס למסלול ה-xarray ומבצע `df.attrs[None]` (כי `field=None` בברירת המחדל של Haber) → `KeyError: None`. ברירת המחדל המתועדת — "The default units of pandas are mg/m**3" — לא ניתנת להשגה.

### B37. מסלול ה-pandas מחזיר NaN בלבד
**קובץ:** `hera/riskassessment/agents/effects/Calculator.py` · **מקובע ב:** `test_risk_calculator.py::TestPandasInput::test_a_constant_exposure_accumulates_linearly`

```python
dt_min = concentrationField.reset_index()[time].diff().apply(lambda x: x.seconds)/60.
return (concentrationField[:-1].fillna(0)*dt_min[1:]).cumsum()*...
```

`concentrationField[:-1]` שומר `DatetimeIndex`; `dt_min[1:]` נושא את ה-`RangeIndex` שנוצר מ-`reset_index()`. pandas מיישר לפי אינדקס, אין חפיפה, והתוצאה **כולה NaN** (וגם באורך N−1 במקום N). גם עם `inUnits` מפורש. כלומר **`CalculatorHaber.calculate` אינו מסוגל להחזיר מספר עבור pandas** — למרות שהדוקסטרינג מפרט את המסלול ואף מבטיח "For pandas, we do not assume that the time is equispaced".

מסלול ה-xarray, לעומת זאת, **נכון לחלוטין**: C=10 mg/m³ בצעדי דקה נותן 10,20,30,40,50, ו-ten Berge עם n=1 זהה ל-Haber.

### B38. `InjuryLevelThreshold` נכשל באופן אטום על סף מספרי
**קובץ:** `hera/riskassessment/agents/effects/InjuryLevel.py:341` · **מקובע ב:** `test_risk_injurylevel.py::test_a_plain_number_is_accepted_like_its_sibling_takes_tl50`

`ureg(parameters["threshold"])` — `parse_expression` של pint מקבל מחרוזות בלבד:

| קלט | תוצאה |
|---|---|
| `'50 mg/m**3'` | ✅ |
| `'50'` | `DimensionalityError` (הודעה סבירה) |
| `50.0` | `AttributeError: 'float' object has no attribute 'replace'` |

המחלקה האחות `InjuryLevelLognormal10DoseResponse` **כן** מקבלת `TL_50=100.0` דרך `tounit`, והדוקסטרינג כאן אומר רק "Must include ``threshold``". יש גם שורה מתה: `thr = ureg(parameters["threshold"])` — מחושבת, לא בשימוש, ומחושבת שוב בשורה הבאה.

### B39. B12 עם קורבן קונקרטי — `getPercent` זורק על עומס רעילות מספרי
**קובץ:** `hera/riskassessment/agents/effects/InjuryLevel.py:351` · **מקובע ב:** `test_risk_injurylevel.py::test_a_plain_numeric_toxic_load_is_evaluated`

```python
ret = numpy.array([1 if tounit(x,self.units) > threshold else 0 for x in ToxicLoad])
```

ה-`threshold` נבנה ב-`ureg(...)` של hera; `tounit(x, ...)` על מספר רגיל מייצר Quantity ב-registry הדיפולטיבי של pint (**זהו B12**). ההשוואה חוצה registries:

```
getPercent(10.0)                    -> ValueError: Cannot operate with Quantity
                                        and Quantity of different registries
getPercent(10*ureg.mg/ureg.m**3)    -> 0   (עובד)
```

**זו ההוכחה ש-B12 אינו אי-עקביות רדומה אלא באג פעיל:** צורת הקריאה הרגילה — עומס רעילות כמספר — אינה ניתנת להערכה כלל. תיקון B12 מתקן גם את זה.

### מה שנמצא תקין
מודל התגובה הלוגנורמלי מדויק: `getPercent(TL_50) = 0.5` בדיוק, `getToxicLoad` הופכית מדויקת בשבע נקודות, מונוטוני וחסום ב-[0,1], ו-`TL/10` ו-`10·TL` נותנים 0.02275 ו-0.97725 — בדיוק ±2σ עבור σ=0.5 בבסיס 10. גם חוקי Haber ו-ten Berge מדויקים במסלול ה-xarray.

---

## אצווה 6 — `measurements/meteorology`

### B40. `radiosonde.py` יורש ממחלקה שהוסרה
**קובץ:** `hera/measurements/meteorology/radiosonde.py` · **מקובע ב:** `test_meteorology_module_imports.py::TestRadiosonde`

```python
class DataLayer(datalayer.ProjectMultiDBPublic):
```

`ProjectMultiDBPublic` **אינו קיים** ב-`hera.datalayer` — שריד לריפקטור שלא הושלם. המודול זורק `AttributeError` בייבוא ומת בכל התקנה.

### B41. נתיב מוחלט של מפתח מסוים, נקרא בזמן ייבוא
**קובץ:** `hera/measurements/meteorology/highfreqdata/__main__.py` · **מקובע ב:** `test_meteorology_module_imports.py::TestHighFreqMain`

```
FileNotFoundError: /home/ilay/hera_unittest_data/measurements/meteorology/
                   highfreqdata/slicedYamim_sonic.parquet
```

הקובץ קורא נתיב מוחלט **בזמן ייבוא**, והנתיב מצביע לספריית הבית של מפתח מסוים. `CLAUDE.md` אוסר נתיבים מוחלטים במפורש. המודול אינו ניתן לייבוא בשום מכונה אחרת.

### לא פגם: `GFS.py`
נכשל מקומית על `ModuleNotFoundError: sklearn`, אבל `scikit-learn==1.6.1` **כן** מוצהר ב-`requirements.txt:266`. זו הסביבה המקומית, לא הקוד. הטסט משתמש ב-`importorskip` כדי לא לדווח על כך כבאג.

### מה שנמצא תקין
סיווג העונות ב-`addDatesColumns` נכון לכל 12 החודשים ותואם את `seasonsdict` שהמודול עצמו מגדיר (DJF/MAM/JJA/SON), כולל מקרה הגלישה שבו דצמבר וינואר חולקים עונה. גם קידוד ה-HHMM (`06:30 → 630`) נכון לאורך היממה, והפונקציה אינה משנה את ה-DataFrame של הקורא.

---

## אצווה 7 — `measurements/GIS`

### B42. ייבוא מודול מריץ דמו וכותב 8 מגה-בייט לתיקייה הנוכחית
**קובץ:** `hera/measurements/GIS/raster/hill2stl.py` · **מקובע במקור ב:** `test_gis_utils.py::TestImportSideEffects`

בתחתית הקובץ, תחת ההערה `# Run the function`, יושב קוד דמו ברמת המודול:

```python
minx=-2; maxx=2; miny=-3; maxy=3
filename='test1.stl'
generate_solid_stl(function, x_range=(minx,maxx), y_range=(miny,maxy), resolution=100, filename=filename)
```

כל `import hera.measurements.GIS.raster.hill2stl` מדפיס ל-stdout וכותב **8.1MB** בשם `test1.stl` לתיקיית העבודה. כל כלי שסורק את עץ החבילה — בניית תיעוד, אינדוקס ב-IDE, איסוף טסטים — משאיר את הקובץ אחריו.

**נסיבה מקלה שאומתה:** `import hera.measurements.GIS` לבדו **אינו** מפעיל את זה; רק ייבוא ישיר של המודול.

**עדכון:** במקביל, מאמץ ניקוי-קוד-מת נפרד ב-master מחק את `hill2stl.py` כולו (אופיין כ"אפס הפניות" לפני שהטסטים האלה נכתבו). לכן `TestImportSideEffects` הוסר מ-`test_gis_utils.py` בזמן פתיחת ה-PR — המחיקה מ-master מתקבלת, אין יותר קובץ לכסות.

### הבדל בין שתי נקודות כניסה להמרת CRS
`GIS/utils.py:convertCRS` מתועד ומחזיר **`list of shapely.geometry.Point`**, בעוד `TopographyToolkit.convertPointsCRS` מחזיר **GeoDataFrame**. שתי דרכים להמיר קואורדינטות עם טיפוסי החזרה שונים. מקובע כטסט עובר ולא כבאג — שתיהן מתועדות — אבל ראוי לאחד.

### מה שנמצא תקין
ההמרה עצמה מדויקת: תל אביב (34.78°E, 32.08°N) נופלת בתחומי רשת ITM הצפויים, ההמרה הלוך-חזור משחזרת עד `1e-6` מעלות, מעלת קו רוחב אחת נותנת 111 ק"מ בתחום הצפוי, והכיווניות נשמרת בשני הצירים. הקונסטנטות `WSG84=4326` ו-`ITM=2039` נכונות.

### B43. שתי טבלאות חספוס סותרות לאותם קודי IGBP — **הממצא המשמעותי מדעית ביותר**
**קובץ:** `hera/measurements/GIS/raster/landcover.py:333` מול `:504` · **מקובע ב:** `test_gis_landcover_roughness.py::TestTheSecondTable`

`_handleType1` (שורה 504) מצטט מאמר אמיתי — Floors et al. 2021, *WES* 6, 1379, טבלה a2 — ומחזיר ערכי אורך חספוס פיזיקליים. `getRoughnessAtPoint` (שורה 333) מחזיק **טבלה שנייה משלו** לאותם קודים בדיוק, מסומנת בקוד `# Example values`:

| קוד IGBP | Floors et al. (מפורסם) | "Example values" | יחס |
|---|---|---|---|
| 0 מים | 0.0001 | 0.01 | **×100** |
| 1 יער מחטני | 1.0 | 0.1 | ÷10 |
| 13 עירוני | 0.8 | 1.1 | ×1.4 |
| 15 שלג וקרח | 0.001 | 1.3 | **×1300** |

הטבלה השנייה היא **רמפה אריתמטית** — הפרשים קבועים של 0.05 או 0.1 מקוד 2 עד 16 — ולא טבלה פיזיקלית. היא גם **הופכת את הסדר**: לפיה יער חלק יותר משלג.

אורך החספוס האווירודינמי מזין ישירות את פרופיל הרוח, ומשם את כל חישוב הפיזור. קורא שמגיע לענף `type_name="IGBP"` מקבל **סימולציה שגויה מהותית** — שגיאה של פי 1300 על שלג וקרח אינה אי-דיוק אלא מספר אחר לגמרי.

שני המסלולים נגישים; מה שמתקבל תלוי בארגומנט `type_name`.

---

## אצווה 8 — `measurements/experiment`

### B44. `Parser_TOA5.parse` הוא stub שמחזיר `None` בשקט
**קובץ:** `hera/measurements/experiment/parsers.py:359` · **מקובע ב:** `test_experiment_parsers.py::TestTOA5Stub`

```python
def parse(self, file):
    """..."""
    pass
```

גוף ריק, ולכן כל קריאה מחזירה `None` — גם על נתיב שאינו קיים. לא מנתח ולא מסרב. הדוקסטרינג מתעד את הפרמטר ואינו מבטיח דבר על ההחזרה. **קיים parser אמיתי ל-TOA5** ב-`meteorology/highfreqdata/parsers/TOA5.ASCIIParser`, ולכן ה-stub אינו המימוש היחיד אלא כפילות ריקה.

### המודול כולו יתום
`hera/measurements/experiment/parsers.py` (41 פונקציות) **אינו מיובא מאף מקום** בעץ הקוד — נבדק בסריקה של כל קבצי ה-`.py`. קוד ה-dispatch שהיה בונה נתיבי `Parser_{name}` נמצא ב-`lowfreqdata/toolkit.py:116,188` **כהערה מוערת**. ה-parsers החיים יושבים ב-`meteorology/highfreqdata/parsers/`. מועמד למחיקה, אבל ראוי לוודא שאין צרכן חיצוני.

### מה שנמצא תקין
קונבנציית ה-factory ש-`CLAUDE.md` דורשת עובדת: `pydoc.locate("...Calculator%s")` פותר את `CalculatorHaber`, `CalculatorTenBerge` ו-`CalculatorMaxConcentration`, מחזיר `None` (ולא זורק) לשם לא מוכר — שזה בדיוק החוזה ש-`Injury.py:57` מסתמך עליו — והמחלקה שנפתרה ניתנת לבנייה ולשימוש.

---

## אצווה 9 — `openFoam` · `LSM` · `mlDL` · `WRF`

### B45. `prepareParams` משנה את התבנית שהיא מקבלת
**קובץ:** `hera/simulations/LSM/template.py:335` · **מקובע ב:** `test_lsm_template_params.py::TestTheTemplateIsNotMutated`

```python
params = template_desc.get('params', {})   # מחזיר את המילון עצמו
params.update(paramsToPrepare)             # וכותב לתוכו
```

הפרמטרים של הקורא **נכתבים לתוך המילון של התבנית**. אומת:

```
template = {"params": {"a": 1}, "units": {}}
prepareParams(template, {"b": 2})   ->  template["params"] == {"a": 1, "b": 2}
prepareParams(template, {"c": 3})   ->  הריצה השנייה רואה גם את b
```

תבנית LSM נועדה לשימוש חוזר בין סימולציות. אחרי הריצה הראשונה היא נושאת את הפרמטרים שלה, וכל ריצה נוספת מקבלת אותם בשקט. אותו דפוס כמו B19 ב-`abstractToolkit`.

### הודעות ל-stdout בזמן ייבוא
`hera/simulations/WRF/wrfDatalayer.py:5,12` — שני `print` חשופים בזמן ייבוא במקום `warnings.warn`. גרסה מתונה של B42: לא כותב קבצים, אבל מזהם את הפלט ואינו ניתן לסינון או לניתוב דרך מנגנון ה-logging.

---

## תיקון עצמי: הנימוק שרשמתי ב-Phase 0 לגבי `torch` היה שגוי

**מה שרשמתי** (ב-`_stubs.py`, בתוכנית Phase 0 ובסיכומים): "modelContainer.py subclasses `torch.nn.Module`, and subclassing a MagicMock raises TypeError".

**מה שנכון:** `modelContainer.py` **אינו** יורש מ-`torch.nn.Module`. הוא מייבא `LightningModuleHera` — מחלקה רגילה — ומשתמש ב-`isinstance` בלבד; האזכור של `nn.Module` נמצא בדוקסטרינג.

**הסיבה האמיתית ש-stub ל-torch נכשל,** שנבדקה בהרצה:

```
ModuleNotFoundError: No module named 'torch.utils'; 'torch' is not a package
```

MagicMock ב-`sys.modules` חסר `__path__`, ולכן ייבוא תת-מודול נכשל. הפתרון הוא טיפול namespace-package כמו ש-PyFoam מקבל — בר-ביצוע, אבל היה מחוץ להיקף. `torch` מוצהר ב-`requirements.txt`, ולכן ב-CI קיים הדבר האמיתי.

הנימוק תוקן ב-`_stubs.py` והמנגנון האמיתי מקובע בטסט.

---

## אצווה 10 — `evaporation` · `deposition` (הקוד שהתיקון של B31 שחרר)

**הערה מקדימה:** שני המודולים האלה **לא היו ניתנים לייבוא** עד תיקון B31. כלומר אף אחד מהם לא רץ מעולם, וששת הממצאים הבאים הם התוצאה הישירה של פתיחתם לבדיקה.

### B46. `diameter.setter` קרוי `ustar` — ודורס את ה-property `ustar`
**קובץ:** `hera/simulations/deposition/models.py:44-46` · **מקובע ב:** `test_deposition_models.py::TestTangledProperties`

```python
@property
def diameter(self): return self._diameter

@diameter.setter
def ustar(self,newdiameter):        # <-- השם ustar
    self._diameter=newdiameter
```

`@diameter.setter` מייצר property חדש, וההשמה שלו לשם `ustar` **דורסת את ה-property `ustar`** שהוגדר למעלה. אומת:

| פעולה | בפועל |
|---|---|
| `obj.ustar` | מחזיר `_diameter` |
| `obj.ustar = 5` | כותב ל-`_diameter` |
| `obj.diameter = 7` | `AttributeError` — אין setter |

`depositionRate_Petroff` עושה `ustar = self.ustar`, כלומר **מחשב שיקוע יבש עם קוטר החלקיק במקום מהירות החיכוך** — אורך במקום מהירות. העברת `ustar` לבנאי אינה משפיעה על התוצאה בכלל.

### B47. `heatFlux.setter` כותב ל-`_ustar`
**קובץ:** `hera/simulations/deposition/models.py:37-38`

```python
@heatFlux.setter
def heatFlux(self,newheatFlux):
    self._ustar=newheatFlux        # <-- השדה הלא נכון
```

`obj.heatFlux = 42` משאיר את `_heatFlux` על 0.1 ו**דורס את `_ustar`**.

### B48. קצב השיקוע אינו תלוי בחספוס פני השטח
**מקובע ב:** `test_deposition_models.py::test_the_rate_responds_to_the_surface_roughness`

`zrough=0.001` ו-`zrough=0.5` נותנים **תוצאה זהה ביט-אחר-ביט** (`0.0000079953`), גם ב-`ustar` ריאלי. הפרמטר `surface` — הסיבה שהמחלקה מקבלת אותו בכלל — אינו משפיע על התשובה.

### B49. קצב השיקוע אינו תלוי בצפיפות החלקיק
**מקובע ב:** `test_deposition_models.py::test_the_rate_responds_to_particle_density`

`500` ו-`5000 kg/m³` נותנים תוצאה זהה ביט-אחר-ביט, למרות ששיקוע גרביטציוני תלוי בצפיפות ישירות. אומת ב-`ustar` ריאלי כדי לשלול תלות ב-B46.

### B50. `evaporationModels` אינו ניתן לבנייה
**קובץ:** `hera/simulations/evaporation/models.py:71` · **מקובע ב:** `test_evaporation_models.py::test_the_class_can_be_constructed`

```python
self._agent = RiskToolkit.getAgent(agent)
```

`getAgent` היא **מתודת מופע** בחתימה `(self, nameOrDesc, version=None)`. הקריאה על המחלקה קושרת את ה-agent ל-`self` ומשאירה את `nameOrDesc` חסר:

```
TypeError: RiskToolkit.getAgent() missing 1 required positional argument: 'nameOrDesc'
```

**כל בנייה נכשלת.** הטסטים בונים את האובייקט ישירות כדי לבדוק את הקורלציות, שהן טהורות ונכונות.

### B51. `flux` אינו יכול לרוץ — יחידות מופשטות ואז נדרשות
**קובץ:** `hera/simulations/evaporation/models.py:122` · **מקובע ב:** `test_evaporation_models.py::TestEvaporativeFlux`

`flux_US` עושה `temperature = tonumber(temperature, ureg.K)` — מספר חשוף — ואז מעביר אותו ל-`agent.physicalproperties.vaporPressure`, שמתחיל ב-`unumToPint(temperature).m_as(ureg.kelvin)`. מספר חשוף הוא dimensionless:

```
DimensionalityError: Cannot convert from 'dimensionless' to 'kelvin'
```

לכל קלט. חישוב שטף האידוי — הפונקציה המרכזית של המודול — אינו עובד.

### מה שנמצא תקין
הקורלציות עצמן מדויקות. FSG ו-EPA תואמים את צורתן המפורסמת, ומדד החזקה משוחזר מ-log-log ל-`1.75` בדיוק. צמיגות האוויר `1.8205e-2·√(T/293)` נותנת את ערך הייחוס ב-293K במדויק, וכפל T ב-2 מכפיל אותה ב-`√2`. Reynolds ליניארי במהירות ובאורך, ו-Schmidt שווה ל-`ν/D` המחושב עצמאית.

---

## אצווה 11 — `openFoam` (חלקי)

### B52. `pandasToFoamFormat` אינו ניתן לקריאה בשום צורה
**קובץ:** `hera/simulations/openFoam/preprocessOFObjects/OFObjectHome.py:238` · **מקובע ב:** `test_openfoam_objects.py::TestPandasToFoamFormat`

```python
@staticmethod
def pandasToFoamFormat(self,data):
    ...
    D = data if self.componentNames is None else data[self.componentNames]
```

**שבור פעמיים.** הדקורטור `@staticmethod` על חתימה שמתחילה ב-`self`, ולכן:

```
home.pandasToFoamFormat(df)              -> TypeError: missing 1 required positional argument: 'data'
OFObjectHome.pandasToFoamFormat(home,df) -> AttributeError: 'OFObjectHome' object has no attribute 'componentNames'
```

והשגיאה השנייה מגלה שהמתודה יושבת **על המחלקה הלא נכונה**: `componentNames` הוא שדה של `OFField`, לא של `OFObjectHome`. המתודה מעולם לא רצה.

### מה שנמצא תקין
וקטור הממדים של OpenFOAM — `[kg m s K mol A cd]` — נכון לכל שבעת המקומות, ונבדק מול חמש כמויות פיזיקליות מוכרות (מהירות `[0 1 -1 ...]`, לחץ `[1 -1 -2 ...]`, צפיפות `[1 -3 0 ...]`, TKE, טמפרטורה). רגיסטר הגדרות השדות נטען נכון, `U` הוא וקטור ו-`p` סקלר, `overwrite` עובד, וההגדרות לא משותפות בין מופעים.

---

## אצווה 12 — `riskassessment` (protectionpolicy, analysis/riskAreas, riskToolkit)

### B53. `ProtectionPolicy.addActions` שבור עם נתיב קובץ — **נפתר במקביל**
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:130` · **מקובע במקור ב:** `test_riskassessment_protectionpolicy.py::TestActionDispatch::test_addactions_from_a_json_string_loads_the_file`

```python
if os.JSONpath.exists(jsonStrOrFile):
```

ל-`os` אין תכונה `JSONpath` (הכוונה הייתה `os.path`). כל קריאה ל-`addActions` עם מחרוזת נתיב קובץ נכשלה עם `AttributeError` — הענף הזה מעולם לא רץ.

**עדכון:** זה היה ממצא אמיתי ומדויק בזמן שתועד. commit `14508043` (PR #1010) תיקן את זה בחזרה ל-`os.path.exists`, במקביל ובאופן בלתי תלוי במאמץ הרחבת הטסטים הזה. הטסט כבר לא xfail — עובד נכון.

### B54. בניית `ActionIndoor` עם `alpha=` — הדוגמה בדוקסטרינג של המחלקה עצמה — שבורה
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:355` · **מקובע ב:** `test_riskassessment_protectionpolicy.py::TestActionIndoor::test_constructing_with_alpha_directly`

```python
elif "alpha" in kwargs:
    self._alpha = (kwargs["alpha"]).ureg(1/ureg.h)
```

ל-`pint.Quantity` אין מתודה `.ureg`. הביטוי הראשון בדוקסטרינג של `ProtectionPolicy` עצמה — `ProtectionPolicy.indoor(alpha=0.2,...)` — נכשל בפועל; רק `turnover=` עובד.

### B55–B57. `RiskToolkit.loadData` — שלושה כשלים ב-`saveMode`
**קובץ:** `hera/riskassessment/riskToolkit.py:182-207` · **מקובע ב:** `test_risk_toolkit.py::TestLoadData`

```python
if agentDoc is None:
    self.addDataSource(...)                    # B55: רץ תמיד, בלי בדיקת saveMode
elif saveMode == TOOLKIT_SAVEMODE_FILEANDDB:
    raise ValueError(...)
else:
    agentDoc.resource = json.dumps(agentDescription)
    agentDoc.desc['version'] = version           # B57: רק version מתעדכן
    agentDoc.save()
```

- **B55:** `TOOLKIT_SAVEMODE_NOSAVE` — "*Just load the data from file and return the datafile*" — נבדק רק בענף `elif`. בטעינה ראשונה של agent חדש (`agentDoc is None`) הקוד קורא ל-`addDataSource` **בלי תלות ב-saveMode בכלל**, כך ש-NOSAVE שומר בפועל.
- **B56:** ה-lookup לפני ה-replace הוא `getDataSourceDocument(datasourceName=name, version=version)` עם ה-**version החדש**, לא הישן. Replace עם version שונה לא מוצא את הרשומה הקיימת, נופל לענף ה-`if agentDoc is None`, ומשאיר **כפילות** בשם אותו agent במקום להחליף.
- **B57:** גם כשה-version תואם וההחלפה בפועל קורית, רק `desc['version']` מתעדכן; שאר `desc` (כמו `effectParameters`) נשאר ישן בזמן ש-`resource` (ה-JSON הגולמי) מוחלף במלואו — metadata ו-payload מתפצלים.

### מה שנמצא תקין ב-ProtectionPolicy
`abstractAction.getAction` מתעל שם פעולה (`"indoor"`/`"masks"`) ל-`ActionIndoor`/`ActionMasks` נכון, שרשור `.indoor().masks()` מחזיר `self` כראוי, ואימות `ValueError` בהיעדר `alpha`/`turnover`/`protectionFactor` עובד. חילוק ריכוז ע"י מסכה (`ActionMasks.compute`) נכון ומדויק.

### B58. **הוסר** — נבדק מול pandas שגוי, לא קיים תחת הגרסה הנעולה
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:390,481` · **מקובע במקור ב:** `test_riskassessment_protectionpolicy.py::TestActionIndoor::test_compute_with_no_time_window_at_all`

```python
abegin = data[self.policy.datetimename].to_series()[0] if abegin is None else abegin
```

אינדקס positional (`[0]`) על `Series` עם `DatetimeIndex` תועד כלא-נתמך (זורק `KeyError`) — אבל זה אומת מול pandas 3.0.2 שהתקין מקומית, לא מול הגרסה שנעולה בפועל ב-`requirements.txt` (pandas==2.2.3). תחת 2.2.3, ה-fallback ל-position עדיין עובד (רק אזהרת `FutureWarning`). קריאה ל-`.indoor(turnover=...)` בלי `begin`/`end`/`enter`/`stay` כלל עובדת נכון בפועל. הטסט כבר לא xfail — בודק את ההתנהגות התקינה.

### B59. חישוב ה-indoor הוא no-op מוחלט — התוצאה תמיד אפס
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:374` · **מקובע ב:** `test_riskassessment_protectionpolicy.py::TestActionIndoor::test_compute_builds_a_low_pass_filtered_field_towards_the_outdoor_value`

```python
Cin[curstep].values = ((Cin[prevstep] + alphanum*dt*Cout[curstep])/(1+alphanum*dt)).values
```

`Cin[curstep]` (אינדוקס לפי dict) מחזיר **עותק**, לא view, בגרסת ה-xarray המותקנת. ההשמה ל-`.values` על העותק נזרקת — `Cin` נשאר אפסים מוחלטים בלי קשר ל-alpha, dt או לריכוז החיצוני. זהו החישוב המרכזי של כל מודל ה-indoor, ומעולם לא באמת רץ.

### B60. `getRiskAreaAlgorithm` — dispatcher מת לחלוטין — **נפתר במקביל**
**קובץ:** `hera/riskassessment/analysis/riskAreas.py:19` · **מקובע במקור ב:** `test_risk_riskareas.py::TestGetRiskAreaAlgorithm`

```python
estimatorCLS = pydoc.locate("pyriskassessment.datalayer.riskAreas.riskAreaAlgorithm_%s" % algorithmName.title())
```

החבילה `pyriskassessment` אינה קיימת בפרויקט או בתלויותיו — המודול הזה יושב ב-`hera.riskassessment.analysis.riskAreas`. `pydoc.locate` החזיר `None` תמיד, ולכן ה-dispatcher מעולם לא הצליח, אפילו לא עבור `"Sweep"` — האלגוריתם היחיד שממומש בקובץ.

**עדכון:** זה היה ממצא אמיתי ומדויק בזמן שתועד. commit `193249e6` ("fix: repair stale module path and nonexistent method call in risk analysis") תיקן את הנתיב ל-`hera.riskassessment.analysis.riskAreas` האמיתי, במקביל ובאופן בלתי תלוי במאמץ הרחבת הטסטים הזה. `getRiskAreaAlgorithm("Sweep")` עובד נכון, וההודעה לשם לא-קיים מציגה נכון "Sweep" כזמין. הטסטים כבר לא xfail.

### מה שנמצא תקין ב-riskAreas/riskToolkit
תכונות ה-setter/getter של `riskAreaAlgorithm_Sweep` (`dxdy`, `workerCount`, `parallel`) עובדות וממירות טיפוסים נכון; `outlayers` הוא read-only כמתועד. `RiskToolkit.getAgent` מבחין נכון בין שם (str), תיאור (dict) וקלט לא-תקין; `loadAgent` מדביק name/version לתיאור לפני העברה ל-`loadData`. `loadData` מטפל נכון בקלט dict, JSON string, ומקרה "לא קובץ/dict" (raises).

`riskAreaAlgorithm_Sweep.calculate()` ו-`RiskToolkit.analysis.getRiskAreas` לא נבדקו — הראשון דורש multiprocessing + geopandas מלא עם demog/isopleths, השני דורש סימולציית LSM חיה ו-agent דוז-רספונס מלא; שניהם מועמדים לטסט אינטגרציה, לא הרמטי.

---

## אצווה 13 — `openFoam` preprocessing, המשך (`OFObject`, `OFList`, `preprocessOFObjects/utils.py`)

### B61. `OFList` הוא קוד מת — אף קורא לא בונה אותו, ואי אפשר לכתוב איתו כלום — **נפתר במחיקה**
**קובץ:** `hera/simulations/openFoam/preprocessOFObjects/OFList.py:55` (נמחק) · **מקובע במקור ב:** `test_openfoam_objects_base.py::TestOFListIsUnusable`

```python
fileStrContent = self.getHeader()
```

זו הייתה השורה הראשונה של `_writeNew` (וגם `_updateExisting` קוראת ל-`_writeNew`). אף מחלקה בהיררכיה של `OFList` (`OFList`, `OFObject`) הגדירה `getHeader` — היא מוגדרת רק על `OFField`, מחלקת-אח לא קשורה. שורה אחת אחרי זה, `self.columnNames` (בניגוד ל-`componentNames` שכן מוגדר) גם אף פעם לא נקבע. כל קריאה, גם בענף הסקלרי, קרסה ב-`AttributeError` לפני שההסתעפות סקלר/וקטור בכלל הגיעה להרצה. חיפוש בכל עץ `hera/simulations/openFoam` לא מצא אף מקום שבונה `OFList(...)`.

**עדכון:** מאמץ ניקוי-קוד-מת נפרד ב-master מחק את `OFList.py` כולו, מאותה סיבה בדיוק (לא מופנה, שבור פנימית) — במקביל ובאופן בלתי תלוי במאמץ הרחבת הטסטים הזה. `TestOFListIsUnusable` הוסר מ-`test_openfoam_objects_base.py` בזמן מיזוג master.

### מה שנמצא תקין
`OFObject` — מחלקת הבסיס הטהורה של `OFField`/`OFList` — עובדת נכון על כל 10 המתודות הציבוריות שלה: וקטור הממדים (`getDimensions`/`dimensionsStr`/`dimensionsList`, כולל השלמת ברירת מחדל לאפס), שמות הרכיבים לפי סוג שדה (סקלר/וקטור/טנזור), גישה ל-`internalField`/`boundaryField`/`processors`/`processorItems` לפי פרוססור, אימות `fieldType` לא תקין, וכתיבה לדיסק (`writeToCase`) בפיצול single/multi-processor כולל ה-escaping של `"proc.*"`.

`ParsedParameterFileToDataFrame` (ב-`preprocessOFObjects/utils.py`) — הפונקציה שמתרגמת שדה OpenFOAM מפורק לטבלת pandas — נבדקה מול אובייקט מזויף שמדמה את הצורה של PyFoam (`obj['key'].val`), בגלל שה-PyFoam המסוטב הוא MagicMock ולא היה מפעיל את ההסתעפויות האמיתיות. כל ההסתעפויות עובדות נכון: אינדקס פרוססור לשדה הפנימי, שורה נפרדת לכל boundary patch עם `value`, שורת NaN ל-patch בלי `value` (כמו `zeroGradient`), דילוג מלא על patch עם `value` ריק, סינון patches של פרוססור (`proc*`) כברירת מחדל, וכיבוד `patchNameList` מפורש. `extractFieldFile` עוטף כשל בפריסה כ-`ValueError` (נבדק ע"י monkeypatch ל-`ParsedParameterFile` עצמה, כי הסטאב לא נכשל על נתיב לא קיים).

---

## אצווה 14 — `riskassessment/agents` (Agents.py, Calculator, Injury, InjuryLevel, thresholdGeoDataFrame)

### B62. `Injury.calculateToxicLoads` שולח `field` למקום הלא נכון עבור Haber
**קובץ:** `hera/riskassessment/agents/effects/Injury.py:282` · **מקובע ב:** `test_risk_injury_pandas_path.py::TestInjuryCalculateToxicLoadsOnPandas::test_toxic_loads_via_haber_raises_on_the_positional_mismatch`

```python
return self.calculator.calculate(concentrationField, field, breathingRate=breathingRate, time=time)
```

`field` מועבר כפרמטר positional שני. זה תואם את `CalculatorTenBerge`/`CalculatorMaxConcentration` (ששם הפרמטר השני שלהן הוא `field`), אבל **לא** את `CalculatorHaber`, שהחתימה שלה היא `(self, concentrationField, breathingRate=..., time=..., field=None, ...)`. כל `Injury` שמוגדר עם מחשבון Haber קורס ב-`TypeError: got multiple values for argument 'breathingRate'` בפעם הראשונה שהוא בשימוש.

### B63. אין דרך להעביר `inUnits` דרך ה-API של `Injury` — קלט pandas קורס תמיד
**קובץ:** `hera/riskassessment/agents/effects/Injury.py:259` · **מקובע ב:** `test_risk_injury_pandas_path.py::TestInjuryCalculateToxicLoadsOnPandas`

שלושת המחשבונים בודקים `hasattr(concentrationField,"attrs")` כדי להבחין בין xarray ל-pandas — אבל ל-`pandas.DataFrame` יש `.attrs` (ריק) כבר מגרסה 1.0, כך שהבדיקה לא מבחינה יותר בין השניים, ו-`df.attrs[field]` מתפוצץ ב-`KeyError`. קריאה ישירה למחשבון עוקפת את זה עם `inUnits=` מפורש, אבל ל-`calculateToxicLoads`/`calculatePointWiseFractionInjured` (החתימות של שתיהן) **אין פרמטר `inUnits` בכלל** — כך שדרך ה-API של `Injury` אין שום דרך לעקוף את זה, בלי קשר לאיזה מחשבון מוגדר.

### B64. הענף הפנדס של Haber מכפיל DataFrame ב-Series בלי `axis=0` — התוצאה NaN מוחלטת
**קובץ:** `hera/riskassessment/agents/effects/Calculator.py:112` · **מקובע ב:** `test_risk_injury_pandas_path.py::TestCalculatorHaberOnPandas`

```python
return (concentrationField[:-1].fillna(0)*dt_min[1:]).cumsum()*breathingRatio * CunitConversion
```

`dt_min[1:]` הוא `Series` עם אינדקס שלם (0,1,2...) בעוד ה-DataFrame מאונדקס בזמן. הכפלת DataFrame ב-Series **בלי `axis=0`** מיישרת לפי **עמודות**, לא לפי שורות — ומכיוון שאין התאמה בין תוויות העמודות ("P1") לתוויות ה-Series (מספרים שלמים), התוצאה היא DataFrame עם עמודות מדומות נוספות ו-NaN מוחלט. `CalculatorTenBerge` נמנעת בדיוק מהמלכודת הזו באותה שורה, ע"י המרה ל-`.values.reshape(...)` לפני הכפל.

### B65. `thresholdGeoDataFrame.project` — קוד שלא ניתן להרצה בפייתון הנוכחי — **נפתר במקביל**
**קובץ:** `hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:77` · **מקובע במקור ב:** `test_risk_thresholdgeodataframe.py::TestProjectAngleDispatch`

```python
if isinstance(meteorological_angle,collections.Iterable):
```

`collections.Iterable` הוסר מ-Python 3.10 (הועבר ל-`collections.abc.Iterable` עוד ב-3.3, וה-alias הוסר לגמרי ב-3.10). כל קריאה ל-`project` קרסה ב-`AttributeError` לפני שהיא בכלל הגיעה לטפל בפרמטרים.

**עדכון:** זה היה ממצא אמיתי ומדויק בזמן שתועד. commit `db405cec` ("fix: repair Python 3.10+ and typo crashes in risk assessment modules") תיקן ל-`collections.abc.Iterable`, במקביל ובאופן בלתי תלוי במאמץ הרחבת הטסטים הזה. ה-dispatch (זווית בודדת מול רשימת זוויות, שתי מוסכמות הזווית) נבדק כעת מול `_project` מדומה (monkeypatch) — ה-pipeline המלא דורש שכבה דמוגרפית geopandas אמיתית ולכן נותר מחוץ לתחום.

### B97. `shiftLocationAndAngle` מאבד את תת-המחלקה בשקט
**קובץ:** `hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:34` · **מקובע ב:** `test_risk_thresholdgeodataframe.py::TestShiftLocationAndAngleLosesTheSubclass`

```python
def shiftLocationAndAngle(self,loc,meteorological_angle=None,mathematical_angle=None,geometry="ThresholdPolygon"):
    ret = self.copy()
    ret[geometry] = self._shiftPolygons(...)
    ret = ret.set_geometry(geometry)
    return ret
```

`thresholdGeoDataFrame` (תת-מחלקה של `geopandas.GeoDataFrame`) לא מגדירה `_constructor` — האידיום הסטנדרטי לשימור זהות תת-המחלקה דרך פעולות pandas/geopandas. תחת הגרסה הנעולה (`geopandas==1.0.1`), `self.copy()` **כבר** מאבד את הטיפוס ומחזיר `GeoDataFrame` רגיל (אומת ישירות: `type(frame.copy())` הוא `GeoDataFrame`, לא `thresholdGeoDataFrame`). התוצאה: כל שרשור של מתודה שמוגדרת רק על `thresholdGeoDataFrame` (כמו `shiftLocationAndAngle` נוסף, או `project`) על התוצאה נכשל ב-`AttributeError`.

### B66. `getVolatility` מחזיר יחידות שגויות — g²/(mol·cm³) במקום g/cm³
**קובץ:** `hera/riskassessment/agents/Agents.py:264` · **מקובע ב:** `test_risk_agents.py::TestPhysicalPropetiesCorrelations`

```python
MW = self.getMolecularWeight().to(ureg.g/ureg.mol)
...
return 1.585287951807229e-5*MW*V1/(temperature+273.16)*ureg.g/ureg.cm**3
```

`MW` הוא `Quantity` ביחידות g/mol, ומוכפל ישירות בביטוי הסופי בלי שה-`mol⁻¹` מבוטל אף פעם — ואז מוכפל שוב ב-`ureg.g/ureg.cm**3`. הממד שמתקבל בפועל הוא `g²/(mol·cm³)`, לא ריכוז מסה `g/cm³` כמתועד. כל שימוש ב-`getVolatility` נותן ערך במימד שגוי.

### מה שנמצא תקין
`Agent`/`PhysicalPropeties` — 10 השיטות הציבוריות שנותרו — עובדות נכון: `fullDescription`, `effectproperties`, `toJSON` (כולל nesting של effects), וכל ה-getters/setters של `PhysicalPropeties` (המרות יחידות למול g/mol ו-cm/s, ברירות מחדל). קורלציית הצפיפות (`getDensity`) ולחץ האדים (`vaporPressure`, נוסחת Antoine) מדויקות מול הנוסחה המתועדת. `CalculatorTenBerge.toJSON`, `InjuryLevelThreshold.toJSON`, ו-`InjuryLevelExponential` המלאה (constructor, `getPercent = 1-exp(-k·D)`, מונוטוניות, `toJSON`) — כולם תקינים. `Injury.toJSON`, האזהרה המתועדת ב-`calculate` (deprecated), ו-`calculatePointWiseFractionInjured`'s דחיית קלט לא-pandas — כולם עובדים כצפוי. `thresholdGeoDataFrame.shiftLocationAndAngle` מסובב ומזיז פוליגונים נכון תחת שתי מוסכמות הזווית, ודוחה קריאה בלי זווית כלל.

---

## אצווה 15 — `OFField.py` (חמש הפונקציות שנותרו מאצווה 11)

### B67. `setFieldFromDataFrame` דורש עמודת `boundary` גם כשאין boundaries בכלל
**קובץ:** `hera/simulations/openFoam/preprocessOFObjects/OFField.py:436` · **מקובע ב:** `test_openfoam_field.py::TestSetFieldFromDataFrame::test_a_dataframe_with_no_boundary_column_currently_raises`

```python
for boundaryName, boundaryData in processordataFrame.query("region=='boundaryField'").groupby("boundary"):
```

`.groupby("boundary")` נכשל ב-`KeyError` אם העמודה `boundary` לא קיימת בכלל ב-DataFrame — וזה מדויק המצב עבור שדה בלי boundary patches, בדיוק כמו ש-`getDataFrame()` (מאותה מחלקה) מחזיר עבור שדה כזה (ראו B61-adjacent ב-`ParsedParameterFileToDataFrame` מאצווה 13). כלומר round-trip `field.setFieldFromDataFrame(field.getDataFrame())` נכשל בוודאות עבור כל שדה בלי boundaries.

### B68. הענף המקבילי לא יכול לאתחל שדה חדש מאפס
**קובץ:** `hera/simulations/openFoam/preprocessOFObjects/OFField.py:432` · **מקובע ב:** `test_openfoam_field.py::TestSetFieldFromDataFrame::test_a_fresh_single_processor_field_currently_raises_on_parallel_data`

```python
for boundaryName, boundaryData in self.data[processorName]['boundaryField'].items():
```

לפני שהוא כותב boundary חדש, `_processorToPyFoam` קורא boundary **קיים** מ-`self.data[processorName]`. שדה שנבנה טרי (למשל דרך `OFObjectHome.getEmptyField`) מאותחל כ-single-processor כברירת מחדל, כך של-`self.data` אין מפתחות `processorN` בכלל. האכלה שלו ב-DataFrame מקבילי — הדרך הטבעית להפוך שדה למקבילי — קורסת ב-`KeyError`; זה עובד רק אם השדה כבר אותחל כמקבילי מראש עם `initialize(noOfProc=...)` תואם.

### הערה טכנית: `preprocessOFObjects/__init__.py` מסתיר את מודול `OFField`
`from .OFField import OFField` ב-`__init__.py` מחליף את התכונה `OFField` על ה-package (שבד"כ מפנה למודול) במחלקה עצמה. `import hera...preprocessOFObjects.OFField as m` מחזיר את **המחלקה**, לא את המודול — כדי לעשות monkeypatch ל-`ParsedParameterFile`/`BoundaryDict` האמיתיים חובה לגשת ישירות ל-`sys.modules["hera.simulations.openFoam.preprocessOFObjects.OFField"]`.

### מה שנמצא תקין
`getDataFrame` (single ו-multi processor, כולל parsing מספר הפרוססור מהמפתח) עובד נכון על מבנה מדומה שמדמה את הפלט של `ParsedParameterFileToDataFrame`. ה-dispatch של `setFieldFromDataFrame` (single לפי היעדר עמודת `processor`, מקבילי לפי נוכחותה) עובד נכון כשה-boundary וה-processors הקיימים תואמים. `readFromCase`/`readBoundariesFromCase` עוטפים `FileNotFoundError` כ-`ValueError("Field not found")` (אומת ב-monkeypatch על PyFoam האמיתי, כי הסטאב לא נכשל על נתיב לא קיים). `addProcBoundary`, `readBoundariesFromCase`, ו-`readFromCase` על תיקייה ריקה רצים בלי לקרוס, אבל לא ניתן לאמת את הערכים שהם כותבים בפועל: PyFoam מסוטב כ-`MagicMock` גורף בלי אחסון אמיתי, כך ש-`mock["key"]=value` לא נשמר וכל קריאה ל-`WriteParameterFile(...)` מחזירה את אותו אובייקט mock משותף — התקרה הטכנית הזו תועדה כבר באצווה 11 לגבי `addBoundaryField`.

---

## אצווה 16 — `measurements/meteorology` (highfreqdata: abstractcalculator, meandatacalculator, turbulencestatistics, parsers)

### B69. שמירה כ-HDF5 שבורה — `to_HDF` במקום `to_hdf`
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/abstractcalculator.py:263` · **מקובע ב:** `test_meteorology_savedata.py::TestSaveDataHandlerHDFIsBroken`

```python
data.to_HDF(path, key)
```

ל-pandas יש רק `to_hdf` (אותיות קטנות). כל שמירה ל-HDF5 קורסת ב-`AttributeError`.

### B70. בניית `MeanDataCalculator` ישירות מ-DataFrame — קוד מתועד שלא רץ
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/meandatacalculator.py:79` · **מקובע ב:** `test_meteorology_meandatacalculator.py::TestConstructionRejectsPlainDataFrames`

```python
if type(TurbCalcOrData) == pandas:
    ...
elif type(TurbCalcOrData) == dask:
```

הבדיקה משווה את ה-**טיפוס** של הארגומנט למודול `pandas`/`dask` עצמם — השוואה שלעולם לא מתקיימת (אין אובייקט שה-`type()` שלו שווה למודול). התיעוד והטיפוסים המוצהרים אומרים במפורש ש-`TurbCalcOrData`/`AverageCalcOrData` יכולים להיות `pandas.DataFrame` או `dask.DataFrame` ישירות — אבל זה קורס תמיד ב-`ValueError`. רק מופע `singlePointTurbulenceStatistics`/`AveragingCalculator` (או `None` ל-`AverageCalcOrData`) עובד בפועל.

### B71. `anisotropyEigs` — `_eig` חסרה `self`
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/meandatacalculator.py:548` · **מקובע ב:** `test_meteorology_meandatacalculator.py::TestStabilityAndAnisotropy`

```python
def _eig(series):
    ...
# בקריאה:
eig_data = self.MeanData.apply(self._eig, axis=1)
```

`_eig` מוגדרת בלי `self` כפרמטר ראשון, אבל נגישה כ-`self._eig` — bound method שמעביר את `self` אוטומטית. כשה-DataFrame.apply קורא לה עם השורה, מתקבלות שתי ארגומנטים לפונקציה שמקבלת רק אחד: `TypeError`. `anisotropyCats` קורסת באותה צורה כי היא קוראת ל-`anisotropyEigs` קודם.

**עדכון (batch29):** B70 חוזרת גם ב-`thresholds()` ו-`filterDates()` (לא רק בבנאי): שני אלה, ב-mode `inplace=False` (ברירת המחדל), בונים `MeanDataCalculator` חדש מ-`filter_obj.data` שהוא `pandas.DataFrame` רגיל — אותה תבנית השוואה שבורה, אז זה תמיד קורס ב-`ValueError`. רק `inplace=True` עובד בפועל. מקובע ב-`TestThresholdsAndFilterDatesInplace`.

### B72. `InMemoryRawData.append` — מתודת pandas שהוסרה
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/turbulencestatistics.py:1478` · **מקובע ב:** `test_meteorology_turbulencestatistics.py::TestInMemoryRawData`

```python
ret = super(InMemoryRawData, self).append(other, ...)
```

`pandas.DataFrame.append` הוסרה ב-pandas 2.0 (הסביבה כאן על 3.0.2). כל קריאה קורסת ב-`AttributeError`.

### B73. `zoL_Sonic` — מנגנון מניעת כפילויות בודק רשימה שאף פעם לא מתמלאת
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/turbulencestatistics.py:968` · **מקובע ב:** `test_meteorology_turbulencestatistics.py::TestMoninObukhovStability`

```python
while 'zoL_Sonic%s' % i in self._TemporaryData.columns:
    if ['zoL_Sonic%s' % i, {'zmd': zmd}] in self._AllCalculatedParams:
        return self
    i += 1
```

הבדיקה קוראת מ-`self._AllCalculatedParams` — אבל בכל המחלקה נכתב רק ל-`self._CalculatedParams` (רשימה **אחרת**, בשם דומה). `_AllCalculatedParams` נשאר `[]` לתמיד, כך שהבדיקה לעולם לא מתקיימת: קריאה חזרה ל-`zoL_Sonic` עם אותו `zmd` יוצרת עמודה כפולה (`zoL_Sonic2`, `zoL_Sonic3`,...) במקום no-op.

### מודולים חסומים בתלות חבילה — תועד ולא נבדק
- **`GFS.py`**: מייבא `sklearn.metrics` וגם `osgeo.gdal` — שתיהן חסרות בסביבה. תועד כבר ב-batch 6 (`TestGFS`), לא הצדיק בדיקה נוספת.
- **`radiosonde.py`**: כבר מתועד כ-B40 (`ProjectMultiDBPublic` לא קיים) — מודול מת, שום פונקציה בו לא ניתנת לבדיקה.
- **`CampbellBinary.py`** (20 פונקציות) ו-**`TOA5.py`**/**`HighFreqToolKit.parseData`/`.loadData`** (5 פונקציות): דורשים קובצי TOB1/TOA5 בינאריים/ASCII אמיתיים עם מבנה מדויק (metadata rows, multi-device columns) — מאמץ בנייה משמעותי, נדחה למעבר עתידי ולא חסם את שאר האצווה.

### מה שנמצא תקין
כל שרשרת החישוב הטהורה ב-`singlePointTurbulenceStatistics` — `sigmaH`, `sigmaHOverUstar/sigmaWOverUstar`, `wind_dir_std`, `sigmaHOverWindSpeed/sigmaWOverWindSpeed`, `w3/w4/w3OverSigmaW3`, `wTKE`, `uStarOverWindSpeed`, `Rvw/Ruw` (חסומים ב-±1), `zOverL_Sonic`, `Lminus1_masked_Sonic` (מסנן נכון לפי סף wT/Ustar), ו-`StabilityMOLength_Sonic` — כולם נכונים על נתוני רוח מלאכותיים. כך גם השרשרת המקבילה ב-`MeanDataCalculator` (horizontalSpeed, sigma, Ustar, sigmaH/UOverUstar, alignedStress עם שימור trace, absWOverSigmaW, effectivez, zOverL, StabilityMOLength, Rvw/Ruw). `AveragingCalculator.getData` מחזיר resample נכון. פונקציות הנרמול (`normalise_sonic_df`, `normalise_trh_df`, `detect_device_type`) עובדות נכון על שינוי שם עמודות, חילוץ metadata, הטלה ל-float, ו-DatetimeIndex — כולל fallback עדין כשאין אינדקס תקין. `InMemoryRawData.to_hdf`/`read_hdf` עובדים נכון (round-trip data+attrs) כשמותקן `pytables`.

---

## אצווה 17 — `measurements/GIS` (חלקי: stlFactory, TilesToolkit, LandCoverToolkit)

**עדכון:** `test_gis_hill2stl.py` (הכיסוי הייעודי ל-4 הפונקציות של `hill2stl.py`) הוסר בזמן מיזוג master — מאמץ ניקוי-קוד-מת נפרד מחק את `hill2stl.py` כולו (ראה גם ההערה באצווה 7/B42). המחיקה מתקבלת, כמו שם.

### B74. `vectorToSTL` — סניף dask שקורא לשם שלא מיובא בקובץ
**קובץ:** `hera/measurements/GIS/utils.py:389` · **מקובע ב:** `test_gis_utils_stlfactory.py::TestVectorToSTLDispatch`

```python
elif isinstance(gpandas, pandas.DataFrame) or isinstance(gpandas, dask.dataframe.DataFrame):
```

`dask` אף פעם לא מיובא במודול הזה. כל עוד הקלט הוא `geopandas.GeoDataFrame` או `pandas.DataFrame` רגיל, ה-`or` הקצר-מעגל מונע מהביטוי השני להיבדק — אבל כל קלט שהוא לא אחד מהשניים (כולל `dask.dataframe.DataFrame` אמיתי, המתועד כנתמך במקביל ל-`rasterizePandas`) קורס ב-`NameError: name 'dask' is not defined` במקום דחייה נקייה.

### הערת קוד: `roughnesslength2sandgrainroughness` מוגדרת פעמיים באותה מחלקה
**קובץ:** `hera/measurements/GIS/raster/landcover.py:561,717` — שתי ההגדרות `@staticmethod` בתוך `LandCoverToolkit` מחשבות את אותה נוסחה (`z0*30`, לפי Desmond et al. 2017 eq. 5). ההגדרה השנייה (717) מחליפה את הראשונה בסמנטיקת גוף המחלקה הרגילה — ההגדרה הראשונה היא קוד מת, אבל ההתנהגות בפועל לא נפגעת (שתיהן מחשבות את אותו הדבר). לא נפתח כממצא B נפרד כי אין כשל התנהגותי.

### מה שנמצא תקין
`stlFactory.heightColumnsNames` (getter/setter) ו-`rasterizeGeopandas` — אינטרפולציה לינארית נכונה מקווי-מפלס (LineString) לרשת רגולרית, מכבדת שם עמודת גובה מותאם. `TilesToolkit` — מתמטיקת Slippy Map סטנדרטית: `tileScaleAtLatLonZoom` נחצה בכל רמת zoom ומצטמצם לקוטבים, `deg2tile`/`tile2deg` הם היפוכים מקורבים, `doctype` ו-`setDefaultTileServer` (נשמר ב-config הפרויקט) עובדים נכון.

### נדחה לעבר עתידי — לא נבדק
`vector/buildings/{analysis,presentation,toolkit}.py` (18), `vector/{demography,topography,toolkit}.py` (16), ו-`raster/tiles.py`'s `getImageFromCorners`/`listImages`/`presentation.plot` — דורשים shapefiles/DEM אמיתיים, שרת tile חי, או נתוני תמונה אמיתיים לבדיקה משמעותית. `GIS/CLI.py` (6 פקודות argparse) — עטיפות דקות מעל הטולקיטים שלמעלה, לא נבדק בנפרד.

---

## אצווה 18 — `measurements/experiment` (experimentHome, experimentAnalysis, dataEngine)

### B75. `getTurbulenceStatistics` — משתנים שלא קיימים בשום מקום
**קובץ:** `hera/measurements/experiment/analysis.py:74` · **מקובע ב:** `test_experiment_analysis_dataengine.py::TestGetTurbulenceStatisticsIsBroken`

```python
def getTurbulenceStatistics(self,sonicData,samplingWindow,height=1):
    ...
    analysis = highfreqtk.analysis.singlePointTurbulenceStatistics(sonicData=kaijoData,
                                                                   ...
                                                                   height=kaijoHeight,
                                                                   ...)
```

הפרמטרים הם `sonicData`/`height`, אבל הגוף קורא ל-`kaijoData`/`kaijoHeight` — שמות שלא מוגדרים בשום מקום בפונקציה. כל קריאה קורסת ב-`NameError` לפני שהיא בכלל מגיעה לטולקיט שהיא מנסה לקרוא.

### B76. `DBconfiguration` — property בלי `return`, בשתי המחלקות
**קובץ:** `hera/measurements/experiment/dataEngine.py:52,236` · **מקובע ב:** `test_experiment_analysis_dataengine.py::TestDaskDataEngineDB`, `TestPandasDataEngineDBConstructionValidation`

```python
@property
def DBconfiguration(self):
    self._dbConfiguration
```

אין `return`. גם ב-`pandasDataEngineDB` וגם ב-`daskDataEngineDB` (אותו קוד בדיוק, פעמיים). ה-property מחזירה `None` תמיד, בלי קשר לתצורה בפועל.

### B77. `dbConnect` — `OperationFailure` שלא מיובא מסווה כשל חיבור אמיתי
**קובץ:** `hera/measurements/experiment/dataEngine.py:66` · **מקובע ב:** `test_experiment_analysis_dataengine.py::TestPandasDataEngineDBConnect`

```python
try:
    self._mongo_client.server_info()
except OperationFailure as e:
    return e
except Exception as e:
    return e
```

`OperationFailure` אף פעם לא מיובא בקובץ. כשל חיבור אמיתי (כל חריגה שאינה `OperationFailure` בעצמה) גורם ל-Python לנסות להתאים אותה ל-`OperationFailure` — הערכת השם הזה עצמה זורקת `NameError`, שנתפס במקום ע"י ה-`except Exception` **החיצוני** ומוחזר בשקט **במקום** החריגה האמיתית. כל כשל חיבור אמיתי מוסתר מאחורי `NameError` על שם לא קיים.

**אזהרת בדיקה:** בניית `pandasDataEngineDB` עם config תקין קוראת בפועל ל-`dbConnect()`, שפותח socket אמיתי דרך `pymongo.MongoClient` ומחכה ל-timeout מובנה של pymongo (~30 שניות) לפני שהחריגה נתפסת — לא נחסם ע"י ה-guard ההרמטי הרגיל. כל הטסטים כאן עוקפים את זה ע"י monkeypatch ל-`pymongo.MongoClient` או ע"י `__new__` (לא קוראים ל-constructor האמיתי).

### מה שנמצא תקין
`experimentHome` — `experimentMap`/`experimentsTable` (שכבת alias מעל `getExperimentsMap`/`getDataSourceTable`), `experimentDataType()` (מחזירה `None` כברירת מחדל — לא באג, רק hook לתאימות לאחור), `keys()`, `__getitem__`, ו-`getExperiment` (שגיאה ברורה כשהניסוי לא נטען). `daskDataEngineDB.__init__` (בדיקת מפתח `'DB'`, **בלי** לפתוח חיבור רשת בבנייה — בניגוד ל-`pandasDataEngineDB`) ו-`connectionString` (בניית URI נכונה). `dataEngineFactory.getDataEngine` — dispatch נכון ל-DASKDB ודחיית dataType לא ידוע.

### נדחה לעבר עתידי — לא נבדק
`experimentSetupWithData` וכל מה שתלוי בו (`defaultTrialSet`, `trialsOfDefaultTrialSet`, `get_devices_image_coordinates`) — דורשים קובץ zip אמיתי של הקמת ניסוי מ-argos. `getDeviceTypeTransmissionFrequencyOfTrial`/`getDeviceTypePlannedMessageCount`/`addMetadata` ב-`experimentAnalysis`, ו-`getDataFromTrial`/`getData`/`getDeviceList`/`getDeviceTable` בשלוש מחלקות ה-data engine — דורשים trialSet/cache אמיתיים. `parsers.py` (27 פונקציות, `Parser_CampbellBinary`/`CampbellBinaryInterface`) — עותק כפול של אותו parser בינארי שנדחה באצווה 16 (B44: `parsers.py` המקורי מתועד כ-orphaned). `presentation.py` (14 פונקציות) — matplotlib בעיקרו, ערך מוגבל לבדיקה. `CLI.py`'s `create_experiment` — אורכסטרציה מלאה של יצירת project+repository, מתאים יותר לטסט אינטגרציה.

---

## אצווה 19 — `simulations/gaussian` (gasCloud.py)

25 הפונקציות שנדחו מאצווה 3 (`abstractGasCloud`, `instantaneousReleaseGasCloud`, `continuousReleaseGasCloud`, `Continuous`), נבנו על גבי `Sigma`/`Meteorology` שאומתו כבר מול Briggs 1973.

### B78. `continuousReleaseGasCloud` — כל ארבע השיטות הציבוריות שלו קורסות
**קובץ:** `hera/simulations/gaussian/gasCloud.py:606,619,625,632` · **מקובע ב:** `test_gaussian_gascloud.py::TestContinuousReleaseGasCloudIsUnusable`

```python
class continuousReleaseGasCloud(abstractGasCloud):
    def getConcentration_cont(self, ...):
        C_noQ = self.getDosage_inst_noQ(...)   # מוגדר רק על instantaneousReleaseGasCloud!
```

`continuousReleaseGasCloud` יורש **ישירות** מ-`abstractGasCloud` — **לא** מ-`instantaneousReleaseGasCloud`. אבל כל ארבע השיטות הציבוריות שלו (`getConcentration_cont`, `getConcentration_cont_NoERF`, `getDosage_cont_NoERF`, `getDosage_cont_doubleNoERF`) קוראות ל-`self.getDosage_inst_noQ`/`self.getDosage_inst_NoERF_noQ` — מתודות שמוגדרות **רק** על `instantaneousReleaseGasCloud`, מחלקת-אח לא קשורה. כל קריאה קורסת ב-`AttributeError`. **המחלקה כולה — הענף היחיד שמטפל בשחרור רציף (Continuous release) — לא שמישה בכלל.** `createGasCloud` בכל זאת מפנה קלט עם יחידות mass/time אליה (בהתאם לתיעוד), כך שהבאג פוגע בכל שימוש שנעשה בדרך המתועדת.

### מה שנמצא תקין
`abstractGasCloud.createGasCloud` — dispatch נכון בין instantaneous/continuous לפי יחידות `sourceQ`, ודחיית קלט לא מוכר. אימות `wind_profile_type`. `getDF_xarray` — Depletion Factor תמיד בין 0 ל-1 ומונוטוני-יורד במרחק הרוח (דעיכה, לא צמיחה). `trapezoidal_integration` — אינטגרציה מצטברת נכונה. כל השרשרת של `instantaneousReleaseGasCloud`: יחס Q ישיר בדיוק בין הגרסאות עם/בלי Q (`getConcentration_inst`/`_noQ`, `getDosage_inst`/`_noQ`/`_NoERF`), מיסוך חלון זמן מדויק (`_bounded`), DF שרק מקטין (לא מגדיל) מנה, המרות רדיולוגיה (Bq, TIAC) עם יחידות פלט נכונות, ו-`get_TIAC_for_dist` שמחזיר את נקודת הרשת הקרובה ביותר עם TIAC מונוטוני-יורד במרחק. `Continuous.__init__` (הקרנל של יהודה) בונה קרנל דעיכה מעריכי בגודל הנכון שמתכנס לכ-90% שחרור עד `timetofinish` כמתועד.

### נדחה לעבר עתידי — לא נבדק
`FallingNonEvaporatingDroplets.py` (217 statements) ו-`DropletCloud.py` (83 statements) — נדחו כבר באצווה 3 ונותרו נדחים; היקף הזמן שנותר בסבב הזה הופנה ל-gasCloud.py, שהניב את הממצא המשמעותי ביותר (B78).

---

## אצווה 20 — `measurements/GIS` המשך (buildings/analysis, buildings/toolkit, vector/toolkit, vector/demography)

המשך על מה שנדחה באצווה 17 — עכשיו עם fixtures סינתטיים קטנים (geopandas אמיתי, בלי shapefiles/DEM חיצוניים).

### B79. `Blocks.getHc` — הבדיקה בודקת מתודה, לא ערך
**קובץ:** `hera/measurements/GIS/vector/buildings/analysis.py` (המחלקה `Blocks`) · **מקובע ב:** `test_gis_buildings_analysis.py::TestLambdaPipeline::test_getHc_currently_returns_the_area_weighted_height_never_none`

```python
def getHc(self):
    if self._LambdaP is not None:
        return self._hc
    else:
        return None
```

`self._LambdaP` (בלי סוגריים) הוא **המתודה עצמה** — bound method — שלעולם לא `None`. הבדיקה לא יכולה אף פעם להגיע ל-`else`, כך ש-`getHc()` תמיד מחזירה `self._hc` (0 כברירת מחדל) בלי קשר אם `_LambdaP()` בכלל רץ. **לטנטי, לא פוגע כרגע בפועל** — הקורא היחיד הקיים (`Lambda()`) תמיד קורא ל-`_LambdaP()` לפני `getHc()` — אבל החוזה המתועד ("מחזיר None אם lambda_P לא הוגדר") לא יכול להתקיים לעולם.

### B80. `cutRegionFromSource` — המרת CRS מבוטלת בשקט
**קובץ:** `hera/measurements/GIS/vector/toolkit.py:134-141` · **מקובע ב:** `test_gis_vector_toolkit.py::TestCutRegionFromSource`

```python
dct = dict(bbox=regionWithCRS) if isBounds else dict(mask=regionWithCRS)

if regionWithCRS.crs is None:
    regionWithCRS.crs = datasourceDocument.desc['desc']['crs']      # מוטציה במקום — dct "רואה" את זה
elif regionWithCRS.crs.to_epsg() != datasourceDocument.desc['desc']['crs']:
    regionWithCRS = regionWithCRS.to_crs(datasourceDocument.desc['desc']['crs'])  # אובייקט חדש — dct לא רואה את זה!
```

`dct` נבנה **לפני** בדיקת ה-CRS, ומחזיק רפרנס לאובייקט `regionWithCRS` המקורי. הענף "אין CRS" עובד במקרה כי הוא ממוטט את האובייקט הקיים (`data.crs = ...` היא הקצאה במקום, לא reprojection). אבל הענף החשוב באמת — **אי-התאמת CRS** — קורא ל-`.to_crs(...)` שמחזיר **אובייקט חדש**, וההקצאה מחדש ל-`regionWithCRS` לא נוגעת ב-`dct` שכבר נבנה. `getData()` מקבל את הצורה **בקואורדינטות המקוריות, לא הממוירות**.

### B81. `DemographyToolkit.shapes` — property שקורא לתכונה שלא קיימת
**קובץ:** `hera/measurements/GIS/vector/demography.py:36` · **מקובע ב:** `test_gis_demography_toolkit.py::TestShapesIsBroken`

```python
@property
def shapes(self):
    return self._shapes
```

`_shapes` לא מוצהר בשום מקום במחלקה — לא כתכונת מחלקה, לא ב-`__init__`. כל גישה ל-`.shapes` זורקת `AttributeError`.

### מה שנמצא תקין
`Blocks` — כל שרשרת חישוב ה-lambda (`_LambdaP`, `_LambdaF`, `_A_f`, `initBuildingsBlock`, `_BuildIndexList`, `iterBlocks`) עובדת נכון על נתונים סינתטיים: אימות פרמטרי חלוקה (`size`/`npxy`/`width`/`height`, כולל דחיית שילובים סותרים), חלוקת דומיין לרשת נכונה, lambdaP/lambdaF מדויקים מול חישוב ידני, hc כממוצע משוקלל-שטח, החרגת בניינים לפי FTYPE. `analysis.ConvexPolygons` מקבצת נכון בניינים קרובים ומפרידה רחוקים, ומחזירה ממוין לפי שטח יורד. `BuildingsToolkit.get_buildings_height`/`filter_buildings_in_area` (סטטיות) — נכונות על השלמת גובה מ-levels, ברירת מחדל "Unnamed", וסינון bounding-box על פני GeoJSON features. `VectorToolkit.geopandasToGeoJson` (ייצוא נכון + דחיית קלט לא-GeoDataFrame) ו-`cutRegionFromSource`'s dispatch ל-`bbox`/`mask` לפי `isBounds` (מלבד B80). `DemographyToolkit` — בנייה (כולל BuildingsToolkit מקונן), `populationTypes`, `setDefaultDirectory`/`FilesDirectory`, `loadData` (רישום datasource נכון), ו-`projectPolygonOnPopulation` (אזהרת deprecation + האצלה מלאה לפרמטרים).

### נדחה עדיין — לא נבדק
`vector/topography.py` (9 פונקציות) — דורש DEM אמיתי או mongomock מורכב יותר; `BuildingsToolkit.getBuildingsFromRectangle`/`getBuildingHeightFromRasterTopographyToolkit`/`buildingsGeopandasToSTLRasterTopography` — דורשים datasource רשום אמיתי או FreeCAD מלא; `analysis.LambdaFromBuildingData`/`calculatePopulationInPolygon` — pipeline מלא עם cache, דורש datasource אמיתי; שתי מחלקות ה-`presentation` (buildings + demography, ~20 פונקציות) — matplotlib בעיקרו.

---

## אצווה 21 — `simulations/openFoam/toolkit.py` (החלקים הטהורים, בלי solver אמיתי)

### B83. `getTimeList` — הענף single-processor מחזיר תמיד רשימה ריקה
**קובץ:** `hera/simulations/openFoam/toolkit.py` (בתוך `getTimeList`) · **מקובע ב:** `test_openfoam_toolkit_pure.py::TestGetTimeListSingleProcessorIsBroken`

```python
timeList = sorted([float(x) for x in os.listdir(case) if (
        os.path.isdir(x) and   # <- x הוא שם יחסי, לא os.path.join(case, x)
        ...
```

הבדיקה `os.path.isdir(x)` נבדקת יחסית ל-**תיקיית העבודה הנוכחית**, לא ל-`case`. הענף המקביל (multi-processor) כמה שורות מטה עושה את זה נכון (`os.path.isdir(os.path.join(case, processorList[0], x))`). התוצאה: לכל case שלא רץ מבפנים תיקיית ה-case עצמה — כלומר כמעט תמיד — הרשימה חוזרת ריקה, גם כשיש תתי-תיקיות זמן אמיתיות.

### מה שנמצא תקין
`processorList` (glob על `processor*`), `read_points_file` (פרסור נכון של קובץ נקודות OpenFOAM, כולל עצירה בסוגר הסוגר), ו-`getMeshExtent` (bounding box נכון מנקודות, `FileNotFoundError` כשהקובץ חסר).

---

## אצווה 22 — `openFoam/CLI.py` (חלק טהור) + `eulerian/abstractEulerianSolver.py`

### B84-B86. `absractEulerianSolver_toolkitExtension` — שלושה NameError ברצף — **נפתרו במקביל**
**קובץ:** `hera/simulations/openFoam/eulerian/abstractEulerianSolver.py` · **מקובע במקור ב:** `test_openfoam_eulerian_solver.py`

- **B84**: `flowType` — ענף ה-`else` הפנה ל-`SIMULATIONTYPE_COMPRESSIBLE`, שם שלא קיים בשום מקום בקובץ (הייבוא בראש הקובץ הוא `FLOWTYPE_COMPRESSIBLE`). כל קריאה עם `incompressible=False` קרסה.
- **B85**: `blockMesh_setBoundFromFile` — הפרמטר נקרא `eulerianWorkFlow`, אבל הגוף בדק `isinstance(eulerianWF, workflow_Eulerian)` — `eulerianWF` לא היה קיים בכלל כפרמטר. קרס תמיד לפני שנוגע בקלט.
- **B86**: `blockMesh_setDomainHeight` — העתק-הדבק של `blockMesh_setBoundFromFile` בלי התאמה: אותו `eulerianWF` לא מוגדר, וגם הפנה ל-`fileName`/`dx`/`dy` — אף אחד מהם לא פרמטר של המתודה הזו (`eulerianWorkFlow, Z, dz`). `Z`, הקלט האמיתי היחיד שלה, אף פעם לא היה בשימוש.

**עדכון:** שלושתם היו ממצאים אמיתיים ומדויקים בזמן שתועדו. commit `276d9a93` ("fix: repair OpenFOAM Eulerian solver name errors and method mismatches") תיקן את כולם — `flowType` משתמש ב-`FLOWTYPE_COMPRESSIBLE`, ושתי מתודות ה-`blockMesh_*` בודקות נכון את `eulerianWorkFlow` — במקביל ובאופן בלתי תלוי במאמץ הרחבת הטסטים הזה. תוך כדי אימות מול venv נעול התגלה גם פער נפרד בשכבת ה-stub (`hermes.workflow` כ-`MagicMock` גולמי גרם ל-`isinstance` לזרוק `TypeError` במקום לעבוד) — תוקן ב-`_stubs.py`. הטסטים כעת בודקים גם את הנתיב התקין (workflow אמיתי) וגם את הדחייה הנקייה (workflow לא תקין).

### מה שנמצא תקין
`Foam_parser_FieldDescription` (CLI) — כותב JSON תקין לתיאור שדות, עם ברירת מחדל ל-`exampleField` כשלא סופקו שדות.

### נדחה
שאר `openFoam/CLI.py` (~25 פקודות) — כולן דורשות טולקיט OpenFOAM אמיתי מחובר (mongomock + hermes workflow), לא רק לוגיקה טהורה.

---

## אצווה 23 — `openFoam/lagrangian/abstractLagrangianSolver.py` (מקורות חלקיקים + ריכוז)

בלי ממצא חדש — כל 17 הפונקציות שכוסו עבדו נכון: `sourcesTypeList`, שש שיטות `makeSource_*` (Point/Circle/Sphere/Cylinder/Rectangle/Cube — כולן בתוך התחום הגיאומטרי המתועד), `writeParticlePositionFile` (דחיית סוג לא ידוע, כתיבת קובץ קואורדינטות OpenFOAM תקין), ו-`analysis.calcConcentrationPointWise` (חלוקה לתאים, סכימת מסה, C=mass/dV נכון).

---

## אצווה 24 — `datalayer/datahandler.py` (המשך פורמטים)

### B87. `DataHandler_JSON_geopandas.saveData` — קורא ל-API הלא נכון, קורס תמיד
**קובץ:** `hera/datalayer/datahandler.py:550` · **מקובע ב:** `test_datalayer_datahandler_formats.py::TestJsonGeopandasHandlerIsBroken`

```python
resource.to_json(fileName,**kwargs)
```

`geopandas.GeoDataFrame.to_json` הפרמטר הפוזיציוני הראשון שלה הוא `na` (ערכים אפשריים: `'null'/'drop'/'keep'`) — **לא נתיב קובץ**, והמתודה **מחזירה מחרוזת JSON**, לא כותבת לדיסק בכלל. `fileName` (נתיב) נופל ל-`na=`, ותמיד נכשל בבדיקת תקינות. `saveData` קורס בכל קריאה. `getData` (טעינה) עצמאית ותקינה — נבדקה מול קובץ שנכתב ידנית.

### מה שנמצא תקין
`DataHandler_time`, `DataHandler_dict` (no-op לשמירה, מעביר ישר בטעינה), `DataHandler_csv_pandas`, `DataHandler_JSON_pandas` (DataFrame ו-Series כולל דגל `pandasSeries`), `DataHandler_netcdf_xarray` (round-trip נתונים ו-attrs — עם הערה ש-`JSONToConfiguration` הופך כל מחרוזת שניתנת לפרסור כיחידת pint ל-`Quantity`, ללא קשר לשם המפתח — לא באג בקובץ הזה), `DataHandler_geopandas` (round-trip GPKG כולל CRS), ו-`DataHandler_image` (round-trip מערך RGB דרך matplotlib).

### נדחה
`DataHandler_geotiff` (דורש raster GDAL אמיתי), `DataHandler_HDF` (`pytables` לא מותקן), `DataHandler_zarr_xarray` (`zarr` לא מותקן).

---

## אצווה 25 — `datalayer/project.py` (השלמת שכבת השאילתות/שמירה)

### B88. `Project.getProjectList` — decorator שגוי, קריאה תקינה קורסת
**קובץ:** `hera/datalayer/project.py:1051-1052` · **מקובע ב:** `test_datalayer_project_remaining.py::TestClassmethodGetProjectListIsBroken`

```python
@staticmethod
def getProjectList(cls,user=None):
    ...
    return list(set(AbstractCollection(connectionName=user).getProjectList()))
```

מתודה שמוגדרת `@staticmethod` אך מקבלת פרמטר ראשון בשם `cls` שאף פעם לא נעשה בו שימוש בגוף הפונקציה. קריאה טבעית `Project.getProjectList()` קורסת ב-`TypeError: missing 1 required positional argument: 'cls'`. הקריאה "עובדת" רק אם המשתמש "מבריח" ערך שרירותי לתוך `cls` (למשל `Project.getProjectList(None)`), מה שמצביע על כך שזה היה אמור להיות `@classmethod` (או סתם השמטת הפרמטר), ואף אחד לא קרא לזה נכון עד כה. פונקציית המודול המקבילה `getProjectList(connectionName=None)` תקינה ולא מושפעת.

### מה שנמצא תקין
`Project.FilesDirectory`/`filesDirectory` (מתאימים זה לזה), `updateProjectNameOnDoc` (תמיד מחתים/דורס את שם הפרויקט הנוכחי), משפחת `save*Data` (Measurement/Cache/Simulation — רושמים לקולקציה הנכונה, כולל override של `type`), `getDocumentByID`, `getMetadata`, `getAllDocuments`, שלישיית get/delete `*DocumentsAsDict`/`delete*Documents` לשלוש הקולקציות, ופונקציית המודול `getProjectList`. הערה חשובה: `getCounterAndAdd` מחזירה את הערך **אחרי** ההוספה (לא לפני, כפי שה-docstring עלול לרמז) — לא באג, אבל התנהגות שקל לפרש לא נכון.

### הערת בדיקה
השוואת `QuerySet` ריק (מ-mongoengine) ל-`[]` ישירות (`queryset == []`) מחזירה `False` למרות ש-`repr` שלו מציג `[]` — יש לעטוף ב-`list(...)` לפני ההשוואה. לא באג בקוד הפרויקט, מלכודת נפוצה בבדיקות בלבד.

---

## אצווה 26 — `measurements/experiment/parsers.py` (המשך: לוגיקת הפרסור בפועל)

### B89. `Parser_CampbellBinary.getData` — שכפול מדויק של B82 בקובץ נפרד
**קובץ:** `hera/measurements/experiment/parsers.py:347` · **מקובע ב:** `test_experiment_parsers_logic.py::TestParserCampbellBinaryGetDataIsBroken`

קובץ זה מכיל עותק עצמאי ומלא (לא import) של `CampbellBinaryInterface` ו-`Parser_CampbellBinary`/`getData` מתוך `measurements/meteorology/highfreqdata/parsers/CampbellBinary.py`, כולל אותו באג בדיוק שתועד שם כ-B82: `columnsIndexes` מחושב כאילו הוא מאנדקס לתוך הרשומה הגולמית (כולל שני שדות timestamp מובילים), בעוד `line` (מ-`_getDataFromStream`, `retval[2:]`) כבר חתך אותם. התוצאה: עמודת נתונים ראשונה חסרה בכל שורה, ו-pandas מסרב לבנות את ה-DataFrame. מכיוון שזה עותק עצמאי, תיקון אחד לא יתקן את השני.

### B90. `Parser_OldStyleMetaDataParquet._getLists` — כיוון שגוי ב-`DataFrame.from_dict`
**קובץ:** `hera/measurements/experiment/parsers.py:111` · **מקובע ב:** `test_experiment_parsers_logic.py::TestParserOldStyleMetaDataParquetGetListsIsBroken`

```python
stationsData = pandas.DataFrame.from_dict(descriptionData['stationLocations'])\
    .reset_index()\
    .rename(columns={"index": "station_name"})
```

`from_dict` בברירת המחדל (`orient="columns"`) הופך כל שם תחנה למפתח **עמודה**, ואת שמות התכונות שלה (`lat`/`lon`/`station_code`) ל-**אינדקס שורות** — ההפך הגמור ממה שהקוד בהמשך מצפה לו (`stndata['lat'].item()` אחרי `.query("station_name==@stn")`). יש להשתמש ב-`orient="index"`. כתוצאה מכך, כל קמפיין עם יותר משדה בודד לכל תחנה קורס ב-`KeyError: 'lat'` בתוך `_getLists`/`getExperimentDict`/`parse`.

### מה שנמצא תקין
`CampbellBinaryInterface` (הקורא הבינארי הנמוך-רמה, זהה למקור התקין) בקובץ זה — עדיין תקין. `Parser_TOA5` כבר תועד קודם כ-stub ריק (B44, ב-`test_experiment_parsers.py`).

### הערה
מודול שלם זה (`experiment/parsers.py`) מתועד כ"יתום" — שום קוד לא מייבא אותו בפרודקשן (ראה `TestParserModuleIsOrphaned` הקיים). שני הבאגים החדשים כאן משפיעים רק אם/כאשר מישהו יחבר את המודול הזה בעתיד.

---

## אצווה 26 (המשך) — `utils/data/toolkit_repository.py` ו-`utils/data/repositoryExport.py`

### B91. `mergeDocumentsIntoRepository` — הדוח מדווח שם פריט לפני הפתרון של התנגשות שמות
**קובץ:** `hera/utils/data/repositoryExport.py:141,146` · **מקובע ב:** `test_utils_data_repository_export.py::TestMergeReportNameIsBrokenOnCollision`

```python
else:
    report["added"].append(item_name)   # <-- לפני הרזולוציה

sectionDict = toolkitDict.setdefault(section, {})
item_name = _uniqueItemName(sectionDict, item_name, entry["contentHash"])   # <-- item_name משתנה כאן
sectionDict[item_name] = entry
```

כאשר שני מסמכים שונים בתוכן חולקים אותו `item_name` (למשל אותו `sourceId`), `_uniqueItemName` פותר את ההתנגשות ומשנה את המפתח שבפועל נשמר במילון — אבל `report["added"]`/`report["overwritten"]` כבר נרשם עם השם *לפני* השינוי. התוצאה: הדוח המוחזר למשתמש עלול להצביע על שם פריט ששייך בפועל למסמך אחר לגמרי (זה שכבר היה קיים תחת אותו שם), בעוד המסמך שבאמת נוסף שוכן תחת שם מופרד (`<name>_<hash>`) שהדוח אף פעם לא מזכיר.

### מה שנמצא תקין
`ToolkitRepository` (register/getToolkitDocument/getToolkitTable, כולל overwrite ואי-overwrite) — כל ההתנהגות תואמת את התיעוד. שאר `repositoryExport.py`: `documentContentHash` (שתי אסטרטגיות, כולל שגיאות), `documentToRepositoryItem` (מיפוי `_cls`, בחירת שם פריט), `deduplicateRepository` (איחוד לפי contentHash, אי-שינוי הקלט).

---

## אצווה 26 (המשך 2) — `measurements/meteorology/lowfreqdata/presentationLayer.py`

בלי ממצא חדש. כוסו: `presenation` (חיווט datalayer/analysis/seasonalPlots/dailyPlots), `Plots.__init__` (מילוני ברירת מחדל לעיצוב), `_getCountourDict`/`_getContourfDict` (מפרטי contour/contourf, כולל שימור מפת הצבעים), ו-`_getcmap` (under/over color, ואי-שינוי `matplotlib.colormaps["jet"]` הגלובלי — `copy.copy` אכן מגן על המפה המשותפת). המתודות הכבדות יותר (`SeasonalPlots`/`DailyPlots` — ציורי contour בפועל) נותרו לא מכוסות, דורשות נתוני רוח אמיתיים בצורה ספציפית.

### הערה (לא באג)
`presenation.__init__` ו-`Plots.__init__` מדפיסים debug prints עם אימוג'ים (`print("📥 ... called")`) ישירות ל-stdout בקוד פרודקשן — לא משפיע על נכונות, אבל שווה ניקוי מתישהו.

---

## אצווה 27 — `simulations/LSM/toolkit.py` (LSMToolkit)

### B92. `getTemplatesTable` — קורס תמיד, כי `loadData` אף פעם לא שומר את ה-`desc` של הקובץ עצמו
**קובץ:** `hera/simulations/LSM/toolkit.py:153` (התוצאה) · השורש ב-`hera/simulations/LSM/toolkit.py:210-214` (`loadData`) · **מקובע ב:** `test_simulations_lsm_toolkit.py::TestGetTemplatesTableIsBroken`

```python
templateName = desc['name']   # desc נטען מקובץ ה-JSON
...
self.addDataSource(dataSourceName=templateName,
                   resource=fileNameOrData,
                   dataFormat=datalayer.datatypes.STRING,
                   version=version,
                   **kwargs)   # <-- ה-desc שנטען מהקובץ אף פעם לא מועבר!
```

`loadData` קורא את קובץ ה-JSON של התבנית לתוך משתנה מקומי `desc`, אבל משתמש בו רק כדי לשלוף `desc['name']`. שאר תוכן ה-JSON (למשל `params`) **לעולם לא** מגיע למסמך שנשמר ב-DB — רק ה-`**kwargs` שהקורא עצמו מעביר ל-`loadData` נשמרים. בהמשך, `getTemplatesTable()` מניחה בעיוורון ש-`desc.pop('params')` יעבוד על כל מסמך שמוחזר מה-DB — וכך קורסת ב-`KeyError: 'params'` על כל תבנית שנרשמה אי-פעם דרך הנתיב הרגיל (ברירת המחדל `saveMode=TOOLKIT_SAVEMODE_FILEANDDB`).

### מה שנמצא תקין
בנאי `LSMToolkit` וה-properties שלו (`to_xarray`/`to_database`/`forceKeep` עם ולידציית boolean, `analysis`, `singleSimulation`), מחלקת `analysis` הפנימית (`coordinateHandler`/`datalayer`), `loadData` עם `TOOLKIT_SAVEMODE_NOSAVE` (לא נוגע ב-DB כלל), `getTemplates`/`getTemplateByName` (מוצאים תבנית שנרשמה, מחזירים `None`/רשימה ריקה כשלא קיימת), ודחיית רישום כפול לאותו שם תבנית.

### נדחה
`getSimulations`/`getSimulationsList` (דורשים מסמכי `SingleSimulation` אמיתיים) ו-`prepareSlurmLSMExecution` (דורש סביבת Slurm/כתיבת קבצי הרצה מורכבת) — נותרו לא מכוסים.

---

## אצווה 28 — `simulations/utils/interpolations.py` (spatialInterpolate)

### B93. `windprofile` — הוגדר בתוך מחלקה בלי `self`
**קובץ:** `hera/simulations/utils/interpolations.py:101` · **מקובע ב:** `test_simulations_utils_interpolations.py::TestWindprofileIsBroken`

```python
class spatialInterpolate():
    ...
    def windprofile(z, uref=3, href=24, he=24, lambdap=0.3, lambdaf=0.3, beta=0.3):
```

מוגדרת בתוך גוף המחלקה בלי פרמטר `self`. קריאה טבעית על מופע (`spatialInterpolate().windprofile(24)`) קושרת את המופע עצמו ל-`z`, וכל שאר הפרמטרים זזים מקום אחד ימינה (`24` הופך ל-`uref`) — התוצאה קריסה עמוקה בתוך גוף הפונקציה עם `TypeError` שלא מרמז על הסיבה האמיתית. פועלת רק אם קוראים לה ישירות דרך המחלקה (`spatialInterpolate.windprofile(24)`), עוקפים למעשה את כל הסיבה להיותה בתוך מחלקה.

### B94. **נשלל** — `interpPandas` נבדק מול pandas שגוי, לא קיים תחת הגרסה הנעולה
**קובץ:** `hera/simulations/utils/interpolations.py:87` · **מקובע במקור ב:** `test_simulations_utils_interpolations.py::TestInterpPandas`

```python
points["interpulation"][i] = self.interp(point=point, stations=stations, topography=topography, ...)
```

תועד כ-assignment משורשר שתחת copy-on-write של pandas **לעולם לא** מעדכן את ה-DataFrame המקורי — אבל זה אומת מול pandas 3.0.2 שהתקין מקומית (שם copy-on-write חובה), לא מול הגרסה הנעולה בפועל ב-`requirements.txt` (pandas==2.2.3). תחת 2.2.3, copy-on-write הוא opt-in ולא ברירת המחדל — ה-assignment המשורשר עדיין עובד (רק `FutureWarning`/`SettingWithCopyWarning`). `interpPandas` **כן** ממלאת נכון את עמודת `interpulation` היום. הבאג הזה עשוי לחזור אם/כש-pandas ישודרג מעבר ל-2.x עם copy-on-write כברירת מחדל — שווה לדעת, לא שווה לקבע כממצא חי.

### B95. `interpArray` — מעביר `topography=` שלא קיים ב-`interpPandas`, קורס תמיד
**קובץ:** `hera/simulations/utils/interpolations.py:98` · **מקובע ב:** `test_simulations_utils_interpolations.py::TestInterpArrayIsBroken`

```python
def interpArray(self, points, stations, columnNames={...}, dx=20, dy=20, C=1000, D=5, Hsl=100, b=150):
    newPoints = pandas.DataFrame({...})
    return self.interpPandas(points=newPoints, stations=stations, topography=newPoints["topography"], ...)
```

`interpPandas` לא מקבל פרמטר `topography` בכלל (יש רק `columnNames` עם מיפוי `"topography"`) — כל קריאה ל-`interpArray` קורסת ב-`TypeError` עוד לפני שהיא מגיעה ל-B94 שהייתה אמורה לרשת ממנו.

### מה שנמצא תקין
`interp` (התאמה מדויקת, שקלול לפי מרחק הפוך, ערכים סקלריים ווקטוריים, התאמת טופוגרפיה), `checkInterpulation` (הפרשי אחוזים בין נתונים מדודים למשוקללים — עובד נכון עבור תחנות בפורמט ה-vector המתועד).

### נדחה
`canopyWindProfile.urbanLogExponentProfile` (פונקציה מונוליתית יחידה, דורשת geopandas גיאומטרי מלא + רשת תלת-ממדית מורכבת) — נותרה לא מכוסה.

---

## אצווה 28 (המשך) — `simulations/gaussian/FallingNonEvaporatingDroplets.py`

### B96. הבנאי שבור לחלוטין — כל יצירת מופע קורסת
**קובץ:** `hera/simulations/gaussian/FallingNonEvaporatingDroplets.py:219` · **מקובע ב:** `test_gaussian_falling_droplets.py::TestConstructorIsBroken`

```python
self._meteorology = MeteorologyFactory().getMeteorology(meteorologyName,**met_kwargs)
```

`MeteorologyFactory` (ב-`Meteorology.py`) **לא מגדירה מתודה בשם `getMeteorology` בכלל** — יש לה רק `getMeteorologyFromU10` ו-`getMeteorologyFromURefHeight`, ששתיהן דורשות פרמטרים שונים לגמרי (לא שם מחרוזת יחיד). כתוצאה מכך, **כל** קריאה ל-`FallingNonEvaporatingDroplets(...)` קורסת ב-`AttributeError` עוד לפני שהבנאי מסיים — המחלקה כולה לא שמישה דרך ה-API הציבורי שלה. נבדק גם ישירות: `hasattr(MeteorologyFactory(), "getMeteorology")` הוא `False`.

בגלל זה, הטסטים בדקו את החלקים הטהורים (property setters/getters, ארבע נוסחאות מקדם הגרר, `hit_ground`) דרך `__new__` תוך עקיפת הבנאי השבור — כדי לוודא שהם *יעבדו נכון* ברגע שמישהו יתקן את B96.

### מה שנמצא תקין (נבדק תוך עקיפת הבנאי השבור)
`position`/`cloudSigma` (ולידציית 3-tuple, גישה ל-x/y/z), `particleMass` (נגזר מ-`rho_p`+`particleDiameter`, `None` עד ששניהם מוגדרים), `beta`/`SpreadFactor`/`g` (קבועים), `N` (מסה כוללת חלקי מסת חלקיק), ארבע פונקציות `_DragCoefficient_*` (Ik/Kelbaliyev15/Fan/Haugen — כל הענפים לפי Re, מונוטוניות, רוויה בערכי Re גבוהים), `correctionCloud_None`, ו-`hit_ground` (אירוע עצירה טרמינלי).

### נדחה
`getTerminalVelocity`, `_VTFunc`, `correctionCloud_Plume`/`_Puff`, `solveToTime`, `_fallingParticle` — כולן תלויות ב-`self.meteorology`, שאי אפשר להשיג כלל כל עוד B96 לא תוקן (אי אפשר גם "לזייף" אותו בקלות בלי לדעת את הממשק המדויק שה-methods הללו מצפים לו).

---

## אצווה 28 (המשך 2) — `utils/data/CLI.py` (הפונקציות הטהורות)

בלי ממצא חדש. כוסו `_parse_query_value` (פרסור מחרוזת query ל-int/float/bool/None/מחרוזת מצוטטת) ו-`load_project_name` (קריאת `projectName` מ-`caseConfiguration.json`, כולל שגיאות קובץ חסר/שדה חסר/שדה ריק). שאר ~26 פקודות ה-CLI (`project_list`, `toolkit_register` וכו') נותרו לא מכוסות — כולן נוגעות ב-Project/toolkit אמיתיים דרך `_setup()` העצל, ודורשות harness גדול יותר מבוסס DB.

---

## אצווה 29 — `datalayer/datahandler.py` (המשך: tif, numpy_dict_array, Class)

בלי ממצא חדש. כוסו: `guessHandler` (פונקציית המודול, מאצילה נכון ל-`datatypes.guessHandler`), `DataHandler_tif.getData` (נבדק מול GeoTIFF אמיתי שנוצר עם `rasterio` — הספרייה כן מותקנת תחת ה-venv הנעול, בניגוד ל-`osgeo`/`tables`/`zarr` שדרושות ל-geotiff/HDF/zarr ונותרות חסומות), `DataHandler_tif.saveData` (NotImplementedError מתועד), `DataHandler_numpy_dict_array` (round-trip `.npz` תקין), ו-`DataHandler_Class` — טעינת מחלקה דינמית: גם נתיב `sys.path` רגיל וגם ה-fallback לטעינה ישירה מ-`resource` (כש-`sys.path` לא מוצא את המודול), כלל המיזוג הנכון בין `desc.parameters` ל-`kwargs` (הראשון גובר), `instantiate=False` (מחזיר את המחלקה עצמה), חסר `classpath` (`ValueError`), ו-`saveData` (`NotImplementedError`).

### הערת בדיקה
כל הטסטים באצווה זו הורצו ואומתו ישירות מול ה-venv הנעול-לפי-`requirements.txt` (`/home/ilay-falach/hera-pinned-venv`), לא מול `heraenv` (שהתגלה סוטה משמעותית — ראה התיקונים בסבב הקודם). ממשיכים כך מכאן ואילך למניעת עוד ממצאים שגויים בגלל גרסאות.

---

## אצווה 29 (המשך) — `simulations/LSM/singleSimulation.py`

### B98. `getConcentrationAtPoint` — קורס תמיד עם נקודה בודדת, בדיוק כפי שמתועד
**קובץ:** `hera/simulations/LSM/singleSimulation.py:129` · **מקובע ב:** `test_lsm_single_simulation.py::TestGetConcentrationAtPointIsBroken`

```python
return self.getConcentration(...)['C'].interp(x=x, y=y, datetime=datetime).values[0]
```

כאשר `x`, `y`, ו-`datetime` הם כולם סקלרים בודדים (בדיוק כפי שהדוקסטרינג מתאר — "x: int... y: int... הערך המבוקש") — `.interp()` מצמצם את כל הממדים המתאימים ומחזיר `DataArray` **אפס-ממדי**. `.values` על מערך כזה הוא מערך numpy 0-ממדי (צורה `()`), ואינדוקס שלו עם `[0]` זורק `IndexError: too many indices for array: array is 0-dimensional`. במילים אחרות: השימוש הכי בסיסי והכי מתועד של הפונקציה — לקבל את הריכוז בנקודה אחת — קורס תמיד.

### מה שנמצא תקין
בנאי `SingleSimulation` — כל שלושת הנתיבים (Dataset ישיר, מסמך עם `getData()`, מחרוזת נתיב — הראשון והשני נבדקו). `params`/`version` (קריאה מ-`_document['desc']`). `getDosage` — גם עם ציר זמן `datetime64` וגם עם ציר זמן מספרי (שניות), חישוב `dt` נכון בשני המקרים, קנה-מידה נכון לפי `Q`, ויחידות `Q`/`C` לפי הבקשה. `getConcentration` — נגזרת ה-Dosage לאורך זמן (ממד `datetime` קטן באחד, כמצופה מ-`diff`), חישוב `C = dDosage/dt` נכון.

---

## אצווה 29 (המשך 2) — `simulations/LSM/template.py`

### B99. `getSimulationByID` — קורא למתודה שלא קיימת על עצמו
**קובץ:** `hera/simulations/LSM/template.py:469` · **מקובע ב:** `test_lsm_template.py::TestGetSimulationByIDIsBroken`

```python
def getSimulationByID(self,id):
    return SingleSimulation(self.getDocumentByID(id))
```

`getDocumentByID` היא מתודה של `Project`/הטולקיט (`self.Toolkit.getDocumentByID`), לא של `LSMTemplate` עצמו — `self.getDocumentByID` לא קיימת. כל קריאה ל-`getSimulationByID` קורסת ב-`AttributeError`.

### הערה: `LSMToolkit.loadData` (B92, כבר מתועד באצווה 27) חוסם בדיקה ישירה של תבניות שנטענות דרך ה-API הרגיל
כדי לבודד את הלוגיקה של `LSMTemplate` עצמה מה-B92 הקיים, הטסטים כאן בונים את המסמך (`_document`) ידנית במקום דרך `toolkit.loadData(...)` — כי B92 גורם ל-`params`/`modelFolder` לזרוק `KeyError` על כל תבנית שנטענת דרך הנתיב הרגיל.

### מה שנמצא תקין
כל ה-properties (`Toolkit`, `doctype_simulation`, `dirPath`, `params`, `version`, `templateName`, `modelFolder` — גם יחסי וגם מוחלט). `prepareParams` (סטטית) — ללא `template_desc`, עם המרת יחידות (`pint`, כולל `duration` שמטופל כדקות גם בלי הגדרה מפורשת), וקאסט ל-`int` עבור `TopoXn`/`TopoYn`. `getSimulationByName`/`getSimulations` — לוגיקת השאילתה נכונה (נבדקה מול `SingleSimulation` מדומה, כי בניית אחת אמיתית דורשת קובצי netcdf אמיתיים): `getSimulations()` בלי ארגומנטים מסנן לפי ה-`params` **של התבנית עצמה** (לא "הכל"), וניתן לדרוס עם קלט מפורש.

### נדחה
`run()`/`_toNetcdf()` — דורשים בינארי LSM אמיתי, קבצי תבנית קלט אמיתיים, והרצת subprocess. `_getSimulationsList`/`getSimulationsTable` — נותרו לא מכוסים.

---

## אצווה 29 (המשך 3) — `datalayer/autocache.py`, `utils/logging/helpers.py`, `gaussian/DropletCloud.py`

בלי ממצא חדש ב-autocache/logging (decorator `cacheFunction` מקצה-לקצה, `clearFunctionCache`, `add_FileHandler`/`add_formatter`).

### B100. `FixedPointClippedDropletCloud` — סדר אתחול הפוך, קורס לפני שמגיע ל-B96
**קובץ:** `hera/simulations/gaussian/DropletCloud.py:148` · **מקובע ב:** `test_gaussian_dropletcloud.py`

`__init__` קורא ל-`super().__init__()` (שמפעיל דרך פולימורפיזם את ה-`_initDropletPosition` **של המחלקה הזו עצמה**, שמשתמש ב-`self.clippedDiameter`) **לפני** שההגדרה `self.clippedDiameter = clippedDiameter` בכלל רצה. כתוצאה, `clippedDiameter` הוא `None` בזמן השימוש, וקורס ב-`TypeError` (numpy `log` על `None`) — עוד לפני שמגיעים לקריסה של B96 (שאר מחלקות הקובץ, שבונות ישירות `FallingNonEvaporatingDroplets`, כן מגיעות ל-B96 בצורה נקייה).

### הערה
כל מחלקות `DropletCloud.py` בונות `FallingNonEvaporatingDroplets` — B96 (מאצווה 28: הבנאי תמיד קורס) מתפשט אליהן כולן.

---

## אצווה 29 (המשך 4) — hermes/pyargos אמיתיים מקושרים ל-venv הנעול + `openFoam/toolkit.py`

מ-נקודה זו, ה-venv הנעול (`hera-pinned-venv`) מקושר קבוע (`.pth`) ל-hermes/pyargos האמיתיים שכבר קיימים ב-`$HOME` (לא stub, לא clone) — פותח גישה לכל openFoam/LSM lagrangian/hermesWorkflowToolkit לבדיקה אמיתית. torch **לא** ב-`requirements.txt` בכלל (ההערה הישנה ב-`_stubs.py` הייתה שגויה) — נשאר לא מותקן, machineLearningDeepLearning עדיין חסום.

### B101. `clearVTKPipelineCache` — קורס בגלל שורת לוג debug שמוערכת תמיד
**קובץ:** `hera/simulations/openFoam/toolkit.py:935` · **מקובע ב:** `test_openfoam_toolkit_vtk_cache.py`

```python
logger.debug(f"Deleting resource {doc['desc']['filterName']} : {doc['desc']['pipeline']['filters']} ")
```

ה-f-string נבנה תמיד (גם אם רמת הלוג לא debug), וקורא ל-`doc['desc']['pipeline']['filters']`. `pipeline` אינו מפתח מתועד/נדרש בסכמה של `vtk_filter` (התיעוד של `getVTKPipelineCacheDocuments` מזכיר רק `regularMesh`/`filterName`/`simulation`) — כל מסמך cache בלי `pipeline` (למשל אחד שנוצר ידנית תואם לסכמה המתועדת) גורם ל-`clearVTKPipelineCache` לקרוס ב-`KeyError`.

### מה שנמצא תקין
`OFToolkit.__init__` (analysis/presentation/OFObjectHome), `getVTKPipelineCacheDocuments` (סינון לפי regularMesh/workflowName/groupName), `getVTKPipelineCacheTable`.

### B102. `splitWorkflowName` — `"".join` במקום `"_".join`, מאבד underscore בשם הבסיס
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py:628` · **מקובע ב:** `test_hermes_workflow_toolkit_names.py`

```python
return "".join(split_name[:-1]), split_name[-1]
```

לא ההיפוך של `getworkFlowName` (שבונה `f"{baseName}_{formatted_number}"`) — כל שם בסיס שמכיל underscore משלו (כמו `"my_flow"`) מאבד אותו: `splitWorkflowName("my_flow_0007")` מחזיר `("myflow", "0007")` במקום `("my_flow", "0007")`.

### B103. `OFLSMToolkit.__init__` — קורא לקבוע `toolkitHome` שלא קיים, לא ניתן לבנייה
**קובץ:** `hera/simulations/openFoam/lagrangian/LSM/toolkit.py:97` · **מקובע ב:** `test_openfoam_lagrangian_lsm_toolkit.py`

```python
self._topography = toolkitHome.getToolkit(toolkitName=toolkitHome.GIS_TOPOGRAPHY,projectName=projectName)
```

`toolkitHome.GIS_TOPOGRAPHY` לא קיים (הקבועים האמיתיים הם `GIS_RASTER_TOPOGRAPHY`/`GIS_VECTOR_TOPOGRAPHY`). כל בנייה של `OFLSMToolkit` דרך `toolkitHome.getToolkit` קורסת ב-`AttributeError` — כל המחלקה לא שמישה. ה-properties (`doctype`, `casePath`, `cloudName`, `parallelCase`) נבדקו בעקיפת הבנאי (`__new__`) ועובדים נכון.

### B104. `meteorology/analysis.py` — `addDatesColumns`/`calcHourlyDist` מוגדרות בלי `self`
**קובץ:** `hera/measurements/meteorology/analysis.py:31,75` · **מקובע ב:** `test_meteorology_analysis_missing_self.py`

אותו דפוס באג כמו B93/B100: שתי המתודות מוגדרות בתוך מחלקת `analysis` בלי פרמטר `self` (`def addDatesColumns(data, ...)`, `def calcHourlyDist(data, Field, ...)`). קריאה טבעית על מופע (`analysis(dl).addDatesColumns(df)`) קושרת את המופע עצמו ל-`data`, וקורסת ב-`AttributeError` (`'analysis' object has no attribute 'assign'`/`'dropna'`). פועלות רק כשקוראים להן ישירות דרך המחלקה. ה-`lowfreqdata/analysis.py` המקביל מגדיר את אותן מתודות נכון (עם `self`).

### B105. `AbstractCalculator` — כל נתיב שנוגע ב-DB הוא קוד מת (`datalayer.Cache` לא קיים)
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/abstractcalculator.py:172,190,221` · **מקובע ב:** `test_meteorology_abstractcalculator.py::TestSaveToDbIsUnreachable`

```python
from hera.datalayer import collection as datalayer
...
datalayer.Cache.getDocuments(...)   # ב-_compute_from_db_and_save / _compute_from_db_and_not_save
...
datalayer.Cache.addDocument(**doc)  # ב-_save_to_db
```

המודול `hera.datalayer.collection` מגדיר רק את המחלקה `Cache_Collection` — אין בו singleton בשם `Cache`. ה-singleton האמיתי חי ב-`hera.datalayer.Cache` (רמת ה-package, ייבוא אחד למעלה). לכן `_compute_from_db_and_save`, `_compute_from_db_and_not_save`, ו-`_save_to_db` (וכפועל יוצא גם `_compute_not_from_db_and_save`, שקוראת ל-`_save_to_db`) קורסים תמיד ב-`AttributeError: module 'hera.datalayer.collection' has no attribute 'Cache'`. מכל ארבעת מצבי ה-`compute()` (`not_from_db_and_not_save` / `from_db_and_save` / `from_db_and_not_save` / `not_from_db_and_save`), רק ברירת המחדל — זו שלא נוגעת ב-DB בכלל — עובדת בפועל.

### B106. `AbstractCalculator._saveProperties` — dict mutable משותף בין כל המופעים
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/abstractcalculator.py:28` · **מקובע ב:** `test_meteorology_abstractcalculator.py::TestSaveSharedAcrossInstancesIsBroken`

```python
class AbstractCalculator(object):
    _saveProperties = {'dataFormat': None}   # class attribute
    def __init__(self, rawData, metadata):
        ...  # לא מעתיק _saveProperties למופע
```

`_saveProperties` מוגדר כ-attribute של המחלקה (dict אחד משותף), ו-`__init__` אף פעם לא יוצר עותק ברמת המופע. `set_saveProperties()` מבצע `self._saveProperties['dataFormat'] = dataFormat` — מאחר ואין attribute מקביל ברמת ה-instance, זו מוטציה על אותו ה-dict המשותף. כתוצאה מכך, קריאה ל-`set_saveProperties` על מופע אחד "דולפת" לכל שאר המופעים — קיימים וגם עתידיים — עד שמישהו יאתחל dict חדש במפורש. גילוי אגבי: הבדיקות ל-B105 נתקלו בזה כי בדיקה קודמת בקובץ קבעה `dataFormat`, וזה "דלף" למופע חדש בבדיקה הבאה.

### B107. `CampbellBinaryInterface._newfloatConvert` — ה-sentinel של NaN לא מחזיר את מה שהוא שומר ב-cache
**קובץ:** `hera/measurements/meteorology/highfreqdata/parsers/CampbellBinary.py:503-512` · **מקובע ב:** `test_meteorology_campbellbinary.py::TestNewFloatConvertNanSentinelIsInconsistent`

```python
def _newfloatConvert(self, key):
    try:
        return self._lut[key]
    except:
        if key == 65183:
            self._lut[key] = float('nan')
            return                      # <-- מחזיר None, לא את ה-nan שנשמר
        val = self._floatConvert(int(key % 256), key / 256)
        self._lut[key] = val
        return val
```

`key == 65183` הוא ה-sentinel של Campbell ל"אין נתונים" בפורמט FP2. הענף שומר `float('nan')` ב-`_lut[key]`, אבל ה-`return` בלי ערך מחזיר בפועל `None` — לא את ה-NaN שנשמר. קריאה שנייה עם אותו `key` כן מחזירה `nan`, כי היא פוגעת ב-`return self._lut[key]` שבתחילת הפונקציה (ה-cache) ומדלגת על הענף השבור לגמרי. התוצאה: אותו קלט מחזיר תשובה שונה בהתאם לכך אם זו הפעם הראשונה שהוא נקרא או לא — `None` בפעם הראשונה, `nan` בכל פעם אחרי זה.

### B108. `thresholdGeoDataFrame._calculatePopulationInPolygon` — `not` על Series רב-איברי, קורס תמיד
**קובץ:** `hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:125` · **מקובע ב:** `test_risk_thresholdgeodataframe.py::TestCalculatePopulationInPolygonIsBroken`

```python
res_intersect_poly = demography.loc[not demography["geometry"].intersection(poly).is_empty]
```

`demography["geometry"].intersection(poly).is_empty` הוא Series בוליאני (per-row). הפעלת `not` על Series כזה — גם אם יש בו רק שורה אחת — מעלה תמיד `ValueError: The truth value of a Series is ambiguous` (בניגוד לסקלר numpy, ל-pandas Series אין קיצור דרך ל-Series באורך 1). המתודה קורסת ב-100% מהמקרים, לא משנה כמה שורות יש ב-demography. השלילה האלמנטרית שהקוד כנראה התכוון אליה היא `~` (או `.apply(...)`), לא `not`. **חיזוק (batch31):** למימוש המקביל ב-`GIS/vector/demography.py::calculatePopulationInPolygon` יש בדיוק את אותה שורה בצורה **נכונה** — `demography.loc[demography["geometry"].intersection(poly).is_empty == False]` — כך שהצורה התקינה כבר קיימת בקוד הבסיס, וזה מחזק שמדובר בטעות מקומית ולא בכוונה. מכיוון ש-`project`→`_project`→`_calculatePopulationInPolygon` היא שרשרת הקריאות המלאה, זה שובר כל נתיב אמיתי של חישוב אוכלוסייה בפוליגון סף.

### B109. `riskassessment/CLI.py createRepository` — נתיב לא מחובר, ה-repository תמיד ריק
**קובץ:** `hera/riskassessment/CLI.py:28-29` · **מקובע ב:** `test_riskassessment_cli.py::TestCreateRepositoryIsBroken`

```python
for file_path in os.listdir(arguments.path):
    agent_json_description = _check_if_agent(file_path)   # <-- לא os.path.join(arguments.path, file_path)
```

`os.listdir(arguments.path)` מחזיר שמות קבצים חשופים (בלי הנתיב המלא), אבל `_check_if_agent` מקבל אותם ישירות בלי `os.path.join` עם `arguments.path`. בפועל `_check_if_agent` פותח אותם יחסית ל-cwd הנוכחי של התהליך, לא יחסית לתיקייה שנסרקת. אלא אם התהליך רץ במקרה כשה-cwd שלו זהה בדיוק ל-`arguments.path`, כל חיפוש קובץ נכשל בשקט (נתפס ב-`except`, מחזיר `False`), ו-`repo['RiskAssessment']['DataSource']` יוצא ריק תמיד — גם כשהתיקייה מלאה בקבצי agent תקינים. התיקון הברור הוא `os.path.join(arguments.path, file_path)`.

### B110. `getWindProfile` — `U_star` מחושב מאותו `z` שהוא מעריך, אז הפרופיל שטוח לגמרי
**קובץ:** `hera/simulations/windProfile/toolkit.py:64,76` · **מקובע ב:** `test_windprofile_toolkit.py::TestGetWindProfileHasNoShear`

```python
U_star = (wind_speed * KARMAN) / np.log(z / z0)
...
U_z = (U_star / KARMAN) * np.log(z / z0)
```

שני הלוגריתמים מצטמצמים אלגברית: `U_z == wind_speed` עבור **כל** גובה. הפונקציה שכל תפקידה לבנות פרופיל רוח מחזירה בפועל את מהירות הרוח של התחנה בכל גובה — בלי שכבת גזירה (shear) בכלל, ובלי שום תלות ב-`z0`: שינוי של פי 50 באורך הגֵרוּיוּת מחזיר פרופיל זהה נומרית. כדי שחוק הלוג יַסֶה משהו, `U_star` חייב להיגזר מגובה מדידה קבוע (גובה התחנה) ולא מה-`z` של הלולאה. אותה בעיה חוזרת בענף האורבני עבור `z > hc`.

### B111. `getWindProfile` — שורת `z=0` סינגולרית, ושני הענפים לא מסכימים מה יש בה
**קובץ:** `hera/simulations/windProfile/toolkit.py:63` · **מקובע ב:** `test_windprofile_toolkit.py::TestGetWindProfileGroundRow`

לולאת הגבהים היא `numpy.arange(0, height + dz, dz)`, כלומר היא תמיד מעריכה את חוק הלוג ב-`z=0` — נקודה סינגולרית (`log(0/z0) = -inf`). הענף הלא-אורבני מחזיר `NaN` ל-`u`/`v`/`U_z` באותה שורה, והענף האורבני מחזיר `-0.0` לאותו קלט. כל פרופיל שמוחזר מתחיל בשורת זבל, ושתי הענפים חולקים חוסר-עקביות לגבי מה בדיוק היא מכילה.

### B112. `_getStationsInRegion` — `wind_stations.json` בנתיב יחסי ל-CWD
**קובץ:** `hera/simulations/windProfile/toolkit.py:118` · **מקובע ב:** `test_windprofile_toolkit.py::TestGetStationsInRegionCannotFindItsDataFile`

```python
with open('wind_stations.json', 'r') as json_file:
```

הקובץ אכן נשלח עם המודול, אבל הוא נפתח כנתיב יחסי חשוף — כלומר יחסית ל-cwd של התהליך. מכל תיקייה שאינה תיקיית המודול, המתודה קורסת ב-`FileNotFoundError`. ה"אחות" שלה `getSpaceTime` קוראת את אותו קובץ נכון, דרך `os.path.dirname(os.path.abspath(__file__))`. אותה משפחת באגים כמו B109.

### B113. `plotOrigin` — `plt.scatter` במקום `ax.scatter`, מתעלם מה-`ax` שהועבר
**קובץ:** `hera/measurements/experiment/presentation.py` · **מקובע ב:** `test_experiment_presentation.py::TestPlotOriginIgnoresItsAxes`

המתודה מקבלת `ax`, מחזירה אותו — אבל מציירת עם `plt.scatter(...)` הגלובלי. כשהקורא מעביר צירים שאינם ה"נוכחיים", הסמן נוחת על צירים אחרים לגמרי וה-`ax` שהוחזר נשאר ריק. אומת: 0 collections על ה-`ax` שהועבר, 1 על הצירים הנוכחיים.

### B114. `generateLatexTable` — קורא ל-`plotDevices` עם פרמטרים שלא קיימים בחתימה
**קובץ:** `hera/measurements/experiment/presentation.py` · **מקובע ב:** `test_experiment_presentation.py::TestGenerateLatexTableIsDead`

```python
fig , _ = self.plotDevices(trialSetName=..., trialName=..., device=device_name, display=False)
```

החתימה האמיתית היא `plotDevices(self, trialSetName, trialName, deviceType, mapName, ax=None, plotkwargs=None)` — אין `device` ואין `display`, ו-`mapName` הוא פרמטר חובה. כל קריאה קורסת ב-`TypeError: plotDevices() got an unexpected keyword argument 'device'` לפני שנוצר משהו. כל נתיב הייצוא ל-LaTeX הוא קוד מת.

### B115. `getRecordByIndex` — המרת שניות דרך מספר ימים float מאבדת ננושנייה
**קובץ:** `hera/measurements/experiment/parsers.py` **וגם** `hera/measurements/meteorology/highfreqdata/parsers/CampbellBinary.py` · **מקובע ב:** `test_experiment_parsers_more.py::TestRecordTimestampsLoseTheirLastNanosecond`

```python
time = pandas.Timestamp(1990, 1, 1) + pandas.Timedelta(days=lastSec / 86400.0, milliseconds=lastmili)
```

החלוקה ב-float לא עושה round-trip: רשומה ששדה ה-SECONDS שלה הוא 11 מפוענחת ל-`1990-01-01 00:00:10.999999999`. אומת: 77 מתוך 600 ערכי השניות הראשונים נפגעים (בערך רשומה אחת מכל 8). מכיוון ש-`getRecordIndexByTime` משווה חותמות זמן בשוויון מדויק, רשומה כזו לא ניתנת למצוא לפי השנייה שהיא באמת מסומנת בה. `pandas.Timedelta(seconds=lastSec)` מדויק.

### B116. `_getDataFromStream` — לולאת ההמרה מתחילה שדה אחד מאוחר מדי
**קובץ:** שני הקבצים מ-B115 · **מקובע ב:** `test_experiment_parsers_more.py::TestFirstDataColumnIsNeverConverted`

```python
for i in range(3, len(retval)):   # אבל retval[0],retval[1] הם SECONDS/NANOSECONDS
...
return retval[0], retval[1] / 1000000, retval[2:]
```

שדה הנתונים הראשון הוא `retval[2]`, והלולאה מדלגת עליו. אומת: רשומה עם שני שדות FP2 **זהים בבתים** מפוענחת ל-`[31, 7936.12109375]` — הראשון חוזר כ-`uint16` גולמי, השני מומר. אותו דבר עם `ASCII(n)`: העמודה הראשונה חוזרת כ-`bytes` גולמי.

### B117. `_newfloatConvert` — `key / 256` במקום `key // 256`
**קובץ:** שני הקבצים מ-B115 · **מקובע ב:** `test_experiment_parsers_more.py::TestFp2LowByteIsDividedInsteadOfShifted`

```python
val = self._floatConvert(int(key % 256), key / 256)
```

הפרמטר השני של `_floatConvert` מתועד כ-low byte, אבל `/` היא חלוקה אמיתית, כך שמה שמועבר הוא `lowbyte + hbyte/256` — ה-high byte "דולף" בחזרה וכל קריאת FP2 סוחבת היסט מדומה. אומת: `_newfloatConvert(0x001F)` מחזיר `7936.12109375` במקום `7936.0` שהבתים האלה מקודדים.

### B118. `_newfloatConvert` — B107 משוכפל בעותק השני של הפרסר
**קובץ:** `hera/measurements/experiment/parsers.py` · **מקובע ב:** `test_experiment_parsers_more.py::TestNanSentinelDoesNotReturnItsOwnCachedValue`

בדיוק אותו פגם שתועד כ-B107 (ענף ה-sentinel של NaN שומר ב-cache אבל לא מחזיר), בעותק השני של הקוד. **הערה רוחבית:** `experiment/parsers.py` ו-`meteorology/highfreqdata/parsers/CampbellBinary.py` הם שני עותקים עצמאיים של אותו פרסר בינארי — ההבדל היחיד ביניהם הוא `parse`. לכן B115, B116 ו-B117 קיימים פעמיים כל אחד, וכל תיקון צריך להיעשות בשני המקומות (או, עדיף, שהכפילות תבוטל).

### B119. `_VTFunc` — מספר ריינולדס בלי צפיפות האוויר
**קובץ:** `hera/simulations/gaussian/FallingNonEvaporatingDroplets.py:331` · **מקובע ב:** `test_gaussian_falling_droplets.py::TestTerminalVelocity`

```python
Re = (uVt*self.particleDiameter/nu_air).magnitude   # nu_air היא צמיגות דינמית
```

`Re = ρVd/µ`, אבל כאן חסרה הצפיפות, כך שהביטוי אינו חסר-מידות בכלל — ו-`.magnitude` פשוט מחזיר בשקט את הערך המספרי של ביטוי מורכב, כלומר פי `10/ρ_air ≈ 8.3` מ-Re האמיתי. שורה 422 ב-`_fallingParticle` מבצעת את אותו חישוב **נכון** (`rho_air*Uabs*self.particleDiameter/nu_air`), מה שמראה את הכוונה. השגיאה מתבטלת במשטר Stokes (שם `Cd·Re = 24` קבוע ו-Re נושר) אבל לא מעליו: טיפה של 3 מ"מ מקבלת 2.831 מ'/ש' במקום 8.156 מ'/ש' שמאזן הכוחות נותן (טיפת גשם של 3 מ"מ נמדדת בערך 8 מ'/ש').

### B120. `correctionCloud_Plume`/`_Puff` — `numpy.mean` על רשימת Quantity, קורס תמיד
**קובץ:** `hera/simulations/gaussian/FallingNonEvaporatingDroplets.py:342,346` · **מקובע ב:** `test_gaussian_falling_droplets.py::TestCloudSigmaCorrection`

```python
Ubar = numpy.mean([self.meteorology.getWindVelocity(z) for z in numpy.arange(0, ...)])
```

זו רשימת פייתון של אובייקטי `Quantity` של pint; numpy לא יכול לבנות מהם מערך, אז שתי המתודות קורסות ב-`ValueError: setting an array element with a sequence.` לכל קלט. הכוונה הייתה `numpy.mean(ureg.Quantity([...], "m/s"))`. מכיוון ש-`correctionCloudFunc` ברירת המחדל היא `"Plume"`, גם `solveToTime` בלתי שמיש בנתיב הזה.

### B121. `numpy.max(z, 0)` — ה-0 מתפרש כ-`axis`, לא כרצפה (פגם רדום)
**קובץ:** `hera/simulations/gaussian/FallingNonEvaporatingDroplets.py:402`

הכוונה הברורה היא "אל תרד מתחת לקרקע", אבל הפרמטר השני של `numpy.max` הוא `axis` ולא רצפה, כך שהחיתוך הוא no-op: `numpy.max(-5.0, 0) == -5.0`. הצורה הנכונה היא `numpy.maximum(z, 0)`. הפגם רדום כרגע כי השורה בלתי-נגישה כל עוד B12 עומד.

### B122. `Injury.__str__` ו-`InjuryLevel.__str__` — בלי `return`
**קובץ:** `hera/riskassessment/agents/effects/Injury.py`, `InjuryLevel.py:170` · **מקובע ב:** `test_risk_injury_effects.py::TestSerialisation`

שתי המתודות מחשבות `json.dumps(self.toJSON(), indent=4)` וזורקות את התוצאה, אז הן מחזירות `None`. כל `str()`/`print()`/f-string על injury או injury level קורס ב-`TypeError: __str__ returned non-string (type NoneType)`. אומת ישירות.

### B123. הרמז "Injuries found" — חיתוך `x[6:]` עיוור
**קובץ:** `hera/riskassessment/agents/effects/Injury.py` · **מקובע ב:** `test_risk_injury_effects.py::TestInjuryFactoryValidation`

```python
[x[6:] for x in dir(module) if x.startswith("Injury")]
```

החיתוך העיוור של 6 תווים מפרסם גם את `Factory` (מ-`InjuryFactory`, שאינו injury) וגם שם ריק (מהמחלקה המופשטת `Injury`). ההודעה בפועל: `Injuries found: ,Exponential,Factory,Lognormal10,Threshold`. אומת ישירות.

### B124. `calculator` ריק → `IndexError` במקום ה-`ValueError` המתועד
**קובץ:** `hera/riskassessment/agents/effects/Injury.py` · **מקובע ב:** `test_risk_injury_effects.py`

```python
try:
    calculatorTypeAndParams = cfgJSON["calculator"]
    calcType,calcParam = [x for x in calculatorTypeAndParams.items()][0]
except KeyError:
    raise ValueError("Calculator not defined")
```

עם `{"calculator": {}}` הגישה לאיבר 0 ברשימה ריקה מעלה `IndexError`, ששומר ה-`except KeyError` לא יכול לתפוס. בדיוק המקרה שהקוד מתכוון לדווח עליו בצורה ידידותית הופך לשגיאה פנימית שלא ניתן להבדיל אותה מבאג. אומת ישירות.

### B125. `units = ureg(units)` — מייצר `Quantity` ולא `pint.Unit`
**קובץ:** `hera/riskassessment/agents/effects/Injury.py` · **מקובע ב:** `test_risk_injury_effects.py::TestDescriptorUnits`

`ureg("mg/m**3")` מחזיר `Quantity(1.0, 'milligram / meter ** 3')`, וה-Quantity הזה מועבר לכל level בתור `units`. לכן `InjuryLevel.units` (מתועד "pint.Unit or None") מקבל טיפוס שגוי, `toJSON()` מסדר אותו כ-`"1.0 milligram / meter ** 3"`, ו-`hera.utils.tounit` (שמקבל רק `Unit` או `str`) נופל לענף ה-fallback שלו. הצורה הנכונה היא `ureg.Unit(units)`. אומת ישירות.

### B126. `hera-GIS` כולו מת — `_setup` מייבא את קבועי ה-CRS מהמודול הלא-נכון
**קובץ:** `hera/measurements/GIS/CLI.py:14` · **מקובע ב:** `test_gis_cli.py::TestSetup`, `TestEveryCommandIsDeadOnArrival`

```python
from hera.utils import WSG84 as _w, ITM as _i
```

`WSG84`/`ITM` חיים ב-`hera/measurements/GIS/utils.py` ואינם מיוצאים מ-`hera.utils` (ה-`__getattr__` העצל שלו מחפש רק ב-unitHandler/jsonutils/query/matplotlibCountour/angle/zipUtils). כל ה-handlers עטופים ב-`@_lazy_setup`, כך שאף פקודת `hera-GIS` לא יכולה לרוץ בכלל: `ImportError: cannot import name 'WSG84' from 'hera.utils'`. מחמיר: `_setup` מסמן `_initialized = True` **לפני** הייבואים, כך שניסיון שני חוזר בשקט עם `WSG84`/`ITM` שנשארו `None` — שגיאה רועשת מתדרדרת לנתיב שקט עם ערכים שגויים.

### B127. `get_landocver` — קורא את `arguments.windDirectionis` (שגיאת כתיב)
**קובץ:** `hera/measurements/GIS/CLI.py:213` · **מקובע ב:** `test_gis_cli.py::TestGetLandcover`

```python
windDirection = None if arguments.windDirection is None else float(arguments.windDirectionis)
```

השם `windDirectionis` לא קיים, כך ש-`hera-GIS landcover getLandcover --windDirection 270` (חובה כש-`isBuilding=True`) קורס ב-`AttributeError` לפני שה-toolkit נקרא בכלל.

### B128. `buildings_raster_toSTL` — מחשב ארבעה ארגומנטים ולא מעביר אף אחד
**קובץ:** `hera/measurements/GIS/CLI.py:168-177` · **מקובע ב:** `test_gis_cli.py::TestBuildingsRasterToSTL`

הפונקציה גוזרת `dxdy`, `inputCRS`, `outputCRS` ו-`dataSourceName` — ואז קוראת `getBuildingsFromRectangle(minx, miny, maxx, maxy, withElevation=True)` בלבד. `--dataSourceName`, `--inputCRS`, `--outputCRS` ו-`--dxdy` נבלעים בשקט, ותמיד נעשה שימוש ב-datasource ברירת המחדל ב-WGS84.

### B129. `hera-experiment create` תמיד נכשל — סדר הקריאות הפוך
**קובץ:** `hera/measurements/experiment/CLI.py:159-160,287-289` · **מקובע ב:** `test_experiment_cli.py::TestCreateExperiment`

`create_experiment` קורא ל-`_create_repository` **לפני** `_make_runtimeExperimentData`, אבל `_create_repository` מסיים בכתיבת `runtimeExperimentData/Datasources_Configurations.json` — תיקייה שרק `_make_runtimeExperimentData` יוצר. הפקודה תמיד מתה ב-`FileNotFoundError`, ומשאירה מאחור `code/`, `data/` ואת ה-repository JSON, בלי שאף פעם נקרא `repository_load`.

### B130. `createNodeRedDeviceMap.sh` נכתב בשם קובץ חשוף ונוחת ב-CWD
**קובץ:** `hera/measurements/experiment/CLI.py:166` · **מקובע ב:** `test_experiment_cli.py`

אותה משפחה כמו B109/B112: הסקריפט נכתב עם שם חשוף במקום `os.path.join(experiment_path, ...)`, כך שעם `--path` הוא נוחת בתיקיית העבודה של התהליך ולא בתיקיית הניסוי. (הפגם הזה באמת לכלך את שורש הריפו בהרצת בדיקה ראשונה, לפני שהבדיקה קובעה עם `chdir`.)

### B131. `workflow_delete` — כותב `completeDelete.py` ומדפיס `completeRemove.py`
**קובץ:** `hera/simulations/CLI.py:120,123` · **מקובע ב:** `test_simulations_cli.py::TestWorkflowDelete`

הסקריפט נכתב ל-`completeDelete.py` אבל ההודעה למשתמש היא `"In order to remove all directories from disk type: python completeRemove.py"`. מי שיעקוב אחרי ההוראה יקבל "No such file". אומת בקוד: שורה 120 פותחת `completeDelete.py`, שורה 123 מדפיסה `completeRemove.py`.

### B132. `workflow_list` — קורא `arguments.object` שלא קיים, ובלי `return` אחריו
**קובץ:** `hera/simulations/CLI.py:436` · **מקובע ב:** `test_simulations_cli.py::TestWorkflowList`

ענף ה"לא נמצא כלום" מדפיס `arguments.object`, אבל ה-parser של `list workflows` מגדיר רק `group`, `projectName`, `nodes`, `parameters`. לכן `hera-workflows list workflows nosuchgroup` קורס ב-`AttributeError`. ה-attribute החסר גם מסתיר פגם שני: לענף אין `return`, כך שגם אם היה `object`, הביצוע היה ממשיך ל-`simDocument[0]` ומעלה `IndexError` על רשימה ריקה.

### B133. `create_workflow_variations` — מדווח על קובץ variations חסר וממשיך לפתוח אותו
**קובץ:** `hera/simulations/CLI.py:291-298` · **מקובע ב:** `test_simulations_cli.py`

`logger.error(f"{path} doesn't point to anything.")` ואחריו — בלי `return` — `open(variation_json_path)`. `hera-workflows variation wf nosuchfile.json` קורס ב-`FileNotFoundError` אחרי שכבר דיווח על השגיאה בצורה מסודרת.

### B134. `create_workflow_variations` — מדווח "workflow not found" וממשיך לאנדקס רשימה ריקה
**קובץ:** `hera/simulations/CLI.py:283-285` · **מקובע ב:** `test_simulations_cli.py`

אותה תבנית כמו B133: `logger.error(...)` בלי `return`, ומיד `workflow_doc_list[0]` → `IndexError`.

### B135. `workflow_compare --file` שבור בדיוק בפורמט ברירת המחדל
**קובץ:** `hera/simulations/CLI.py:592-593,615` · **מקובע ב:** `test_simulations_cli.py::TestWorkflowCompare`

עבור `--format pandas` (ברירת המחדל של ה-argparse) `output` נשאר ה-DataFrame הגולמי, ואז `outputFile.write(output)` מעלה `TypeError: write() argument must be str` — ומשאיר קובץ ריק מאחור. כלומר `hera-workflows compare wf_0 wf_1 --file out` נכשל תמיד אלא אם המשתמש מציין פורמט אחר במפורש.

### B136. `getHermesWorkflowFromJSON` — `pydoc.locate("hermes.workflow")` מחזיר מודול, לא מחלקה
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py` · **מקובע ב:** `test_hermes_workflow_toolkit_db.py::TestGetHermesWorkflowFromJSON`

נתיב ה-fallback (workflow בלי solver) מאתר `hermes.workflow` — שהוא ה**מודול** `hermes/workflow/__init__.py`; המחלקה עצמה היא `hermes.workflow.workflow.workflow`. לכן טעינת workflow ללא solver קורסת ב-`TypeError: 'module' object is not callable`.

### B137. מסנן ה-`component` מקבל prefix כפול, אז אף template לא נמצא
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py` · **מקובע ב:** `test_hermes_workflow_toolkit_db.py::TestComponentFilteredTemplateLookups`

`getHermesFlowTemplate`, `getHermesNodeTemplate` ו-`listHermesNodesTemplates` מעבירים `desc__component=...`, אבל `collection.getDocuments` מזרים כל kwarg נוסף דרך `dictToMongoQuery(..., prefix="desc")` — והתוצאה היא השדה `desc__desc__component`. לכן החיפוש תמיד מחזיר `None`, בזמן ש-`getDataSourceDocumentsList(component="Flow")` כן מוצא את אותו מסמך. (ל-`listHermesNodesTemplates` יש פגם שני שמוסתר מאחורי הראשון: היא קוראת `doc.desc['desc']`/`doc.desc['datasourceName']`, ו-`addDataSource` לא יוצר את ה-`desc` המקונן הזה.)

### B138. `getHermesWorkflowFromDB` עובד רק אם מעבירים לו רשימה
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py` · **מקובע ב:** `test_hermes_workflow_toolkit_db.py::TestGetHermesWorkflowFromDB`

עבור מזהה שאינו רשימה, `getWorkflowListDocumentFromDB` מחזיר את ה-`QuerySet` של mongoengine כמו שהוא; אז `getHemresWorkflowFromDocument` בודק `isinstance(documentList, list)` (שקר עבור QuerySet), עוטף את כל ה-QuerySet ברשימה בת איבר אחד וקורא `.desc` עליו. `tk.getHermesWorkflowFromDB("flow_0000")` קורס ב-`AttributeError: 'QuerySet' object has no attribute 'desc'`, בזמן ש-`tk.getHermesWorkflowFromDB(["flow_0000"])` עובד.

### B139. `longFormat=True` בלתי שמיש מכל נקודות הכניסה להשוואה
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py` + `hera/utils/dataframeutils.py:137` · **מקובע ב:** `test_hermes_workflow_toolkit_db.py::TestCompareWorkflowObj`

`compareWorkflowObj` מקודד `changeDotToUnderscore=True` בקוד, ואז `dataframeutils` מפעיל `.replace(".", "_")` על `ret.T.columns` — שבפורמט הארוך הם אינדקס שורה שלם (int). התוצאה: `AttributeError: 'int' object has no attribute 'replace'`. שורש הבעיה ב-`dataframeutils`, אבל רק הקומבינציה הזו מגיעה אליה, ו-`test_utils_dataframeutils.py` לא בדק את שני הדגלים יחד.

### B140. כל השוואה מבוססת-DB מתה: `resource=` במקום `Resource_path=`
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py:907,937` · **מקובע ב:** `test_hermes_workflow_toolkit_db.py::TestCompareWorkflows`, `TestCompareWorkflowInGroupAndWorkflowTable`

```python
workflow(..., name=..., resource=simulationDoc['resource'])
```

הבנאי של `hermes.workflow` מקבל `Resource_path` (`workflow.py:91`), לא `resource` — ושאר נקודות הקריאה באותו קובץ עצמו (שורות 279, 721, 740) מעבירות נכון `Resource_path=`. לכן `compareWorkflowInGroup`, `workflowTable` ו-`compareWorkflows` על קבוצה קורסים ב-`TypeError: workflow.__init__() got an unexpected keyword argument 'resource'`. רק הנתיב של קבוצה ריקה (שלא בונה כלום) ונתיב הקבצים של `compareWorkflows` עובדים. אומת בקוד.

### B141. `findAvailableName` קורא מונה אחד ו-`addWorkflowToGroup` מקדם מונה אחר
**קובץ:** `hera/simulations/hermesWorkflowToolkit.py:591` מול `:718` · **מקובע ב:** `test_hermes_workflow_toolkit_db.py::TestFindAvailableName`

```python
newID = self.getCounterAndAdd(f"simulations_{simulationGroup}")   # findAvailableName
groupID = self.getCounterAndAdd(groupName)                        # addWorkflowToGroup
```

שני מונים נפרדים לחלוטין, שלעולם לא מסונכרנים. אחרי הוספת שני workflows לקבוצה `flow` (שקיבלו `flow_0000`, `flow_0001`), `findAvailableName("flow")` מחזיר `(0, "flow_0000")` — שם שכבר תפוס.

### B142. `InputForModelsCreator.render` — שומר הבטיחות עצמו קורס ב-`UnboundLocalError`
**קובץ:** `hera/simulations/utils/inputForModelsCreation.py` · **מקובע ב:** `test_simulations_inputformodels.py::TestRenderGuardIsBroken`

```python
def render(self, savePath=None):
    if self._templateName is None or self._paramsMap is None:
        print("templateName and paramsMap are not set yet")
    else:
        ...
        renderedTemplate = template.render(self._paramsMap)
        ...
    return renderedTemplate      # <-- מוגדר רק בתוך ה-else
```

בנתיב השומר, `renderedTemplate` מעולם לא נקשר, ולכן ה-`return` מעלה `UnboundLocalError`. כלומר הענף שנועד להתמודד עם "עוד לא הוגדר" גרוע יותר מאי-בדיקה בכלל: המשתמש מקבל גם את ההודעה וגם קריסה. אומת ישירות.

### B143. `PlotFields.plotFields` — לא יכול לצייר שדה בודד
**קובץ:** `hera/presentation/basicplots.py` · **מקובע ב:** `test_presentation_basicplots.py::TestPlotFieldsSingleFieldIsBroken`

```python
fig, axes = plt.subplots(1, plotsNum, figsize=[20*plotsNum, 10])
for i in range(plotsNum):
    self.plotField(..., ax=axes[i], ...)
```

עבור `plotsNum == 1` הפונקציה `plt.subplots(1, 1)` מחזירה `Axes` בודד ולא מערך, ואז `axes[i]` מעלה `TypeError: 'Axes' object is not subscriptable`. רק שני שדות ומעלה עובדים. אומת ישירות.

### B144. `PlotFields.plotFields` — מקבל `axes` וזורק אותו מיד
**קובץ:** `hera/presentation/basicplots.py` · **מקובע ב:** `test_presentation_basicplots.py::TestPlotFieldsIgnoresItsAxes`

החתימה כוללת `axes=None`, אבל השורה הראשונה בגוף הפונקציה היא `fig, axes = plt.subplots(...)` — כלומר הפרמטר נדרס לפני שנעשה בו שימוש. קורא שמעביר צירים מוכנים מקבל בשקט פיגורה חדשה לגמרי, והצירים שלו נשארים ריקים. אותה משפחה כמו B113 בשכבת ההצגה של הניסוי. אומת ישירות.

### B145. `abstractWorkflow.controlDict` — קורא `'ControlDict'` בעוד שהצומת נקרא `'controlDict'`
**קובץ:** `hera/simulations/openFoam/OFWorkflow.py:88` מול `:46` · **מקובע ב:** `test_openfoam_ofworkflow.py`

`_requiredNodeList` (שורה 46) דורש `controlDict` באות קטנה — וזו גם הצורה בכל תבניות ה-Flow תחת `hera/doc` — אבל ה-property מחזירה `self['ControlDict']` באות גדולה. ה-`__getitem__` של hermes מחזיר `None` לצומת לא מוכר במקום לזרוק, כך שה-accessor מחזיר `None` בשקט על כל workflow שרק עתה עבר את הוודיקציה של המחלקה עצמה. (זה עובד רק על תבניות הפיזור הלגראנז'יאניות, שבמקרה משתמשות בצורה עם האות הגדולה.)

### B146. `abstractWorkflow.fvScheme` — קורא שם צומת שלא קיים באף מקום
**קובץ:** `hera/simulations/openFoam/OFWorkflow.py:99` · **מקובע ב:** `test_openfoam_ofworkflow.py`

ה-property מחזירה `self['fvScheme']` (יחיד), בעוד שהצומת הוא `fvSchemes` ב-`_requiredNodeList` ובכל תבנית בריפו; חיפוש `"fvScheme"` על כל קבצי ה-`.py`/`.json` לא מחזיר כלום. לכן ה-property לא יכולה להחזיר שום דבר מלבד `None`.

### B147. `abstractWorkflow.__setitem__` — קריאה ל-super עם שם משובש ומְמוּנְגָּל
**קובץ:** `hera/simulations/openFoam/OFWorkflow.py:158` · **מקובע ב:** `test_openfoam_ofworkflow.py`

```python
super().__setitem(key=key,value=value)      # חסרים שני קווים תחתונים בסוף
```

חסרים ה-`__` הסופיים, ובתוך גוף מחלקה שם שמתחיל בשני קווים תחתונים עובר name mangling ל-`_abstractWorkflow__setitem` — שלא קיים על שום דבר. לכן **כל** השמה של צומת מעלה `AttributeError`. מחמיר: שתי תופעות הלוואי (`self.buildAllRun[key] = ...` ו-`self.fileWriter[key] = ...`) מתבצעות **לפני** הכשל, כך שחריגה שנתפסת משאירה רשומת הרצה ורשומת file-writer לצומת שמעולם לא נוסף.

### B148. `set_blockMesh_blockHeight` — מזיז את הפאה התחתונה במקום העליונה
**קובץ:** `hera/simulations/openFoam/OFWorkflow.py` · **מקובע ב:** `test_openfoam_ofworkflow.py`

```python
for i in range(len(verticsList[4:])):
    verticsList[i][2] = Z
```

הפרוסה `verticsList[4:]` משמשת רק לאורך שלה (4), ואז האינדוקס נעשה ב-`i` עצמו — כך שקודקודים 0-3 (הבסיס) עולים ל-`Z` והפאה העליונה נשארת בגובה הקודם, כלומר הבלוק יוצא הפוך. האינדקס הנכון הוא `i+4`. אומת: עם גבולות `z ∈ [0,30]` וקריאה `set_blockMesh_blockHeight(Z=100, dz=5)` מתקבל `[100,100,100,100,30,30,30,30]`.

### B149. `set_blockMesh_blockHeight` — מספר תאי ה-z שנגזר הוא תמיד אפס
**קובץ:** `hera/simulations/openFoam/OFWorkflow.py` · **מקובע ב:** `test_openfoam_ofworkflow.py`

בהמשך לאותה מתודה, `minZ = verticsList[0][2]` קורא את קודקוד 0 — אחד מארבעת אלה שזה עתה נדרסו ל-`Z` (B148) — כך ש-`minZ == Z` ולכן `cellCount[2] = (Z-minZ)/dz` יוצא **תמיד** `0.0`. מספר תאים אפס הופך את הרשת לבלתי-ניתנת-לבנייה, וגם הטיפוס הוא float בניגוד ל-int ש-`set_blockMesh_blockBoundaries` כותב.

### B150. `workflowGroupID` — ה-setter פוסל בדיוק את הטיפוס שה-getter שלו מחזיר
**קובץ:** `hera/simulations/openFoam/OFWorkflow.py` · **מקובע ב:** `test_openfoam_ofworkflow.py`

`addWorkflowToGroup` שומר את `groupID` כ-int (מונה הקבוצה), וה-getter מחזיר אותו כמו שהוא — אבל ה-setter דורש `isinstance(value, str)` ואחרת מעלה `ValueError("Group name must be a string")` (הודעה שהועתקה מ-setter של *שם* הקבוצה). לכן `wf.workflowGroupID = wf.workflowGroupID` קורס, וכתיבה של id בחזרה מחייבת מחרוזת — מה שגורם לטיפוס המאוחסן לסתור כל מסמך אחר באותה קבוצה.

### B151. `dataToolkit.addRepository` — קובע סיומת לפי חיפוש תת-מחרוזת בכל הנתיב
**קובץ:** `hera/utils/data/toolkit.py` · **מקובע ב:** `test_utils_data_toolkit.py::TestAddRepositoryExtensionHeuristic`

הבדיקה היא `"json" not in repositoryPath` — חיפוש תת-מחרוזת על **כל** הנתיב ולא על הסיומת. מאגר שיושב בתוך תיקייה שיש בשמה "json" נרשם בלי סיומת, עם resource שלא קיים. (זה תפס את מערך הבדיקות עצמו: `tmp_path` של pytest נקרא על שם הבדיקה, אז בדיקה עם "json" בשם משחזרת את זה בטעות.)

### B152. `addRepository` — משאיר את פרויקט ברירת המחדל פתוח לכתיבה אחרי כשל
**קובץ:** `hera/utils/data/toolkit.py` · **מקובע ב:** `test_utils_data_toolkit.py::TestAddRepositoryLeavesTheDefaultProjectWritable`

המתודה מדליקה `_allowWritingToDefaultProject = True`, קוראת ל-`addDataSource`, ואז מכבה — בלי `try/finally`. אם `addDataSource` זורק (שם כפול בלי `overwrite`), הדגל נשאר `True` ופרויקט ברירת המחדל, שאמור להיות לקריאה בלבד, מקבל כתיבות שרירותיות לכל אורך חיי המופע.

### B153. סעיף `Cache` ב-repository JSON לא יכול להיטען בכלל
**קובץ:** `hera/utils/data/toolkit.py:222` מול `:266` · **מקובע ב:** `test_utils_data_toolkit.py::TestCacheSectionIsUnreachable`

```python
Measurements=lambda toolkit, itemName, docTypeDict, overwrite, basedir: ...
Cache=lambda toolkit, itemName, itemDesc,      overwrite, basedir: ...   # שם פרמטר שונה
...
handler(..., docTypeDict=docTypeDict)                                    # תמיד בשם הזה
```

ה-lambda של `Cache` קורא לפרמטר השלישי `itemDesc` בעוד ש-Measurements/Simulations קוראים לו `docTypeDict`, והמפעיל תמיד מעביר `docTypeDict=`. לכן כל סעיף `Cache` מעלה `TypeError`, שנבלע ב-`except Exception: logger.error` שעוטף את הקריאה. התוצאה: הטעינה נראית מוצלחת ואפס מסמכי cache נוצרים. אומת בקוד.

### B154. `_makeItemPathAbsolute` — `bool("False")` הוא `True`
**קובץ:** `hera/utils/data/toolkit.py:485` · **מקובע ב:** `test_utils_data_toolkit.py::TestMakeItemPathAbsoluteMisreadsTheStringFlag`

```python
isRelativePath = bool(theItem.get("isRelativePath", True))
```

ב-repository JSON הדגל נשמר כמחרוזת `"True"`/`"False"` (ה"אח" `_handle_DataSource` מניח בדיוק את זה), ו-`bool("False")` הוא `True` — כך שresource שסומן במפורש כלא-יחסי עדיין מחובר ל-`basedir`. אומת בקוד.

### B155. `_DocumentHandler` — ה-`ValueError` הידידותי הוא קוד מת
**קובץ:** `hera/utils/data/toolkit.py` · **מקובע ב:** `test_utils_data_toolkit.py::TestDocumentHandlerRejectsUnknownDocumentTypes`

הקוד עושה `getattr(toolkit, f"get{documentType}Documents")` ורק **אחר כך** `if retrieveFunc is None: raise ValueError(...)`. `getattr` בלי ברירת מחדל זורק `AttributeError` קודם, כך שה-`ValueError` שמפרט את סוגי המסמכים החוקיים לא יכול לרוץ אף פעם.

### B156. `_handle_Function` — מפעיל את המתודות על ה-dataToolkit במקום על ה-toolkit היעד
**קובץ:** `hera/utils/data/toolkit.py:426` · **מקובע ב:** `test_utils_data_toolkit.py::TestHandleFunctionTargetsTheWrongObject`

הקוד עושה `getattr(self, itemName)` — כאשר `self` הוא ה-dataToolkit שיושב על פרויקט ברירת המחדל (לקריאה בלבד) — במקום `getattr(toolkit, itemName)`. סעיף `Function` יושב **תחת מפתח של toolkit**, כלומר הקריאות מיועדות ל-toolkit ההוא בפרויקט היעד; הפרמטר `toolkit` לא בשימוש בכלל. כל קריאה שנוגעת בקונפיגורציה מתה, כי פרויקט ברירת המחדל אוסר קונפיגורציה. הדוקסטרינג עצמו מתעד את ההתנהגות השגויה ("method name on ``self``"). אומת בקוד.

### B157. `_handle_DataSource` — `pop` מוקדם מדי, ואז קריסה שמסתירה את הסיבה האמיתית
**קובץ:** `hera/utils/data/toolkit.py` · **מקובע ב:** `test_utils_data_toolkit.py::TestHandleDataSourcePopsResourceFilePathTooEarly`

```python
with open(theItem.pop("resourceFilePath")) as f:
```

ה-`pop` מתבצע בזמן חישוב הארגומנטים, כך שאם הקריאה נכשלת ה-`except` מתעד את הסיבה האמיתית — אבל הביצוע ממשיך ל-`addDataSource(**theItem)` **בלי** `resource` וגם בלי `resourceFilePath`. מה שהמשתמש רואה בסוף הוא `TypeError: addDataSource() missing 1 required positional argument: 'resource'` במקום `FileNotFoundError`.

### B158. `_resolveDocumentsForExport` — בדיקת `is None` מול API שזורק
**קובץ:** `hera/utils/data/toolkit.py` · **מקובע ב:** `test_utils_data_toolkit.py::TestResolveDocumentsForExportRejectsUnknownIds`

`doc = proj.getDocumentByID(d)` ואחריו `if doc is None: raise ValueError("Document id not found...")`. אבל `getDocumentByID` עובר דרך `QuerySet.get()` של mongoengine, שזורק `DoesNotExist` במקום להחזיר `None` — כך שהשומר הוא קוד מת וההודעה הידידותית לא נראית אף פעם.

### B159. `_DocumentHandler` — מעביר את דגל הבקרה `isRelativePath` כשדה נתונים
**קובץ:** `hera/utils/data/toolkit.py` · **מקובע ב:** `test_utils_data_toolkit.py::TestDocumentHandlerForwardsTheControlFlagAsData`

ה-handler **קורא** את `isRelativePath` מתוך ה-item (דרך `_makeItemPathAbsolute`) אבל אף פעם לא מוציא אותו ב-`pop`, ואז קורא `add<Type>Document(**theItem)` — שמקבל רק `resource`/`dataFormat`/`type`/`desc`. לכן כל רשומת Measurements/Simulations/Cache שמצהירה `isRelativePath` מתה ב-`TypeError`, בשקט, דרך הלואדר. בנוסף שני ה-handlers **לא מסכימים איפה הדגל יושב**: `_handle_DataSource` קורא אותו מהעוטף שליד `"item"`, ו-`_DocumentHandler` מתוך `"item"` עצמו.

### B160. `_process_row` — קורא `.x`/`.y` מרשימה שמוחזרת מ-`convertCRS`
**קובץ:** `hera/measurements/experiment/experiment.py:458` מול `hera/measurements/GIS/utils.py:92` · **מקובע ב:** `test_experiment_experiment.py::TestProcessRow`

```python
return pd.Series([pp.x[0], pp.y[0]])          # experiment.py:458
...
return list(gdf.to_crs(outputCRS).geometry)   # utils.py:92 — רשימה של Points
```

`convertCRS` מחזירה **רשימה** של נקודות shapely, ולרשימה אין `.x`. הצורה הנכונה היא `pp[0].x`. זה הורג את כל ענף ה-`outputCRS=ITM` של `get_devices_image_coordinates`. אומת בקוד.

### B161. `get_devices_image_coordinates` — קורא `Latitude`/`Longitude` בעוד ש-argos כותב באות קטנה
**קובץ:** `hera/measurements/experiment/experiment.py` · **מקובע ב:** `test_experiment_experiment.py::TestGetDevicesImageCoordinatesColumnNames`

argos פורש מיקום של מכשיר בניסוי לעמודות `latitude`/`longitude` באות קטנה (`fillContained.spread_attributes`), אבל hera קוראת `Latitude`/`Longitude`. לכן כל ניסוי שבו המכשירים ממוקמים על המפה — כלומר כל ייצוא לדוגמה של argos — נופל ב-`KeyError: 'Latitude'`. המתודה עובדת רק אם למכשירים יש במקרה מאפייני טקסט חופשי בשמות `Latitude`/`Longitude` בדיוק.

### B162. `TrialWithdata.getData(withMetadata=True)` — ממזג על עמודה שלא קיימת בטבלה הזו
**קובץ:** `hera/measurements/experiment/experiment.py` · **מקובע ב:** `test_experiment_experiment.py::TestTrialGetDataWithMetadata`

המיזוג נעשה עם `right_on="entityName"`, אבל `argos.Trial.entitiesTable` קורא לעמודה `deviceItemName`; `entityName` קיימת רק ב-`Experiment.entitiesTable`/`EntityType.entitiesTable`. לכן `trial.getData(deviceType=..., withMetadata=True)` נופל ב-`KeyError: 'entityName'`.

### B163. שני initializers של hera פחות הגנתיים מ-argos, ולכן לא טוענים ייצוא אמיתי
**קובץ:** `hera/measurements/experiment/experiment.py:500,590` מול `argos/experimentSetup/dataObjects.py:792` · **מקובע ב:** `test_experiment_experiment.py::TestTrialSetInitialisationIsNotDefensive`

```python
for trial in self._metadata["trials"]:        # hera:500
for entity in self._metadata["entities"]:     # hera:590
for trial in self._metadata.get('trials', []) # argos:792 — הגנתי
```

ייצוא v3.0.0 אמיתי **משמיט** את `trials` עבור סוג trial שאין לו אף אחד — וייצוא הדוגמה של argos עצמו עושה בדיוק את זה (אומת: `example_exp/exp_simple`, סוג "New Trial Type 1", בלי מפתח `trials`). hera זורקת `KeyError: 'trials'` מתוך הבנאי ולא מצליחה לטעון אותו בכלל. אומת בקוד, זה מול זה.

### B164. `EntityTypeWithData.getDataTrial` — מדליק `perDevice` בלי להעביר `deviceName`
**קובץ:** `hera/measurements/experiment/experiment.py` · **מקובע ב:** `test_experiment_experiment.py::TestEntityTypeGetDataTrialPerDevice`

המתודה מעבירה `perDevice=StoreDataPerDevice` אבל אף פעם לא `deviceName`, בעוד שהענף הפר-מכשירי ב-`parquetDataEngineHera.getData` נפתח ב-`assert deviceName, "If perDeivce=True then deviceName should be defined!"`. כלומר לכל סוג ישות שהנתונים שלה **כן** נשמרים פר-מכשיר — המקרה היחיד שבו הדגל נדלק — הקריאה נופלת ב-`AssertionError`. `EntityWithData.getData` כן מעביר `deviceName`, מה שמראה שאותו דגל עובד שם.

### B165. `parquetDataEngineHera.getDataFromTrial` — קוד מת, ועוד שני פגמים מאחוריו
**קובץ:** `hera/measurements/experiment/dataEngine.py` · **מקובע ב:** `test_experiment_dataengine_parquet.py::TestGetDataFromTrial`

המתודה קוראת `self.experimentObj.experimentSetup.trialSet[trialSet][trialName]`, אבל ל-`experimentSetupWithData` אין attribute בשם `experimentSetup` — הוא **עצמו** ה-setup וחושף `.trialSet` ישירות. לכן `AttributeError` בהוראה הראשונה שמשתמשת בארגומנט. מאחורי זה, על אותו נתיב, יש עוד שניים: ברירת המחדל `trialSet=None` משימה את **כל המילון** של ה-trial sets במקום שם/מפתח, וענף ה-`withMetadata` קורא `entitiesTable()` למרות שזה property (`TypeError`). שלושת אלה מועתקים גם ל-`pandasDataEngineDB.getDataFromTrial` ול-`daskDataEngineDB.getDataFromTrial`.

### B166. `openFoam/CLI.py` — הגדרות כפולות הופכות פקודות רשת ל-no-op שקט
**קובץ:** `hera/simulations/openFoam/CLI.py:259,263,267` מול `:640,654,708` · **מקובע ב:** `test_openfoam_cli.py::TestPlaceholderHandlers`

`foam_mesh_blockMesh`, `foam_mesh_setDomainHeight` ו-`IC_hydrostaticPressure` מוגדרות **פעמיים** בקובץ. ההגדרות המאוחרות מסתירות את המוקדמות, ושתיים מהן הן `pass` ריק — כך שהפקודה יוצאת בקוד 0 בלי לעשות כלום. (ההגדרה המוקדמת של `IC_hydrostaticPressure` בשורה 267 גם מקבלת פרמטר בשם `argumets` בשגיאת כתיב, אבל היא בלתי נגישה בכל מקרה.) אומת בקוד.

### B167. `if 'projectName' not in arguments` — בדיקה שלעולם לא מתקיימת
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestFoamCreateEmptyProjectNameFallbackIsDead`

argparse מגדיר `--projectName` עם `default=None`, כך שה-attribute **תמיד** קיים ונפילת החזרה ל-`caseConfiguration.json` שמאחורי הבדיקה הזו היא קוד מת; `projectName=None` מגיע ישר ל-`getToolkit`. משפיע על `Foam_createEmpty`, `stochasticLagrangian_dispersionFlow_list` ו-`IC_hydrostaticPressure` החי. האידיום הנכון (`'projectName' in arguments and arguments.projectName is not None`) מופיע ב-`stochasticLagrangian_dispersion_create` באותו קובץ.

### B168. `caseConfiguration.json` פגום נמחק ונדרס בשקט
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestFoamCreateEmptyClobbersABrokenConfiguration`

`Foam_createEmpty`/`IC_hydrostaticPressure` עוטפים את `loadJSON` ב-`except:` עירום, מדווחים "not found", ואז **דורסים** את הקובץ ב-`{"projectName": null}`. קונפיגורציה תקינה שנשמרה עם שגיאת תחביר קטנה נמחקת במקום שהמשתמש יתבקש לתקן אותה.

### B169. `dispersionFlow create` בולע `FileExistsError`
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestDispersionFlowCreateSwallowsFileExists`

`logger.error` ואז נפילה מסוף הפונקציה — הפקודה יוצאת בקוד 0 בלי שנוצר שום דבר. אותה משפחה כמו B133/B134.

### B170. `postProcess toParquet` — מחשב `cache` ולא מעביר אותו
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestPostProcessToParquetDropsTheCacheFlag`

ה-handler קובע `cache = False` (ענף התיקייה) או `True` (ענף ה-DB), ואז קורא ל-`getCaseResults(...)` **בלי** `cache` — שברירת המחדל שלו `True`. כלומר ענף התיקייה תמיד עובד בניגוד למה שחושב.

### B171. `UnboundLocalError` על שם dispersion לא מוכר
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestPostProcessUnknownDispersionName`

ב-`toParquet` וגם ב-`toVTK`: `outputName` מוקצה רק בתוך `if os.path.isdir(...)`, כך ששם שאינו workflow שמור וגם לא תיקייה מת ב-`UnboundLocalError` באתר השימוש במקום בהודעת "not found" המיועדת.

### B172. `foam_solver_simulations_list --format pandas --file` — `TypeError`
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestFoamSolverSimulationsListPandasToFile`

`output = res` (ה-DataFrame) ואז כתיבה לקובץ טקסט. בדיוק אותה צורה כמו B135 ב-`workflow_compare`.

### B173. `dispersionFlow list` — פורמט ברירת המחדל קורס תמיד
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestDispersionFlowListPandasFormat`

ענף `--format pandas` מדפיס כל קבוצה אבל **לא מקצה** `output`, והביצוע ממשיך ל-`print(output)` → `UnboundLocalError`.

### B174. `dispersionFlow list` — קורא למתודה שלא קיימת
**קובץ:** `hera/simulations/openFoam/CLI.py:330` · **מקובע ב:** `test_openfoam_cli.py::TestDispersionFlowListCallsAMethodThatDoesNotExist`

```python
res = tk.workflowCompare(workflowsType="stochasticLagrangianSolver")
```

חיפוש `workflowCompare` בכל `hera/` מחזיר **בדיוק את אתר הקריאה הזה** ואף הגדרה. המתודה האמיתית היא `OFToolkit.compareWorkflows`. אומת בקוד.

### B175. `dispersion create --overwrite` לא יכול להשתמש בשדה זרימה מה-DB
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestDispersionCreateOverwriteIgnoresTheDatabase`

`if len(doc) == 0 or arguments.overwrite:` גורם ל-`--overwrite` לזרוק גם את חיפוש שדה הזרימה ב-DB ולא רק את תיקיית ה-case, כך שיצירה מחדש נכשלת ב-"myFlow not found!".

### B176. `--outputDirectory` נבלע בשקט ב-`postProcess toVTK`
**קובץ:** `hera/simulations/openFoam/CLI.py:569` · **מקובע ב:** `test_openfoam_cli.py::TestPostProcessToVTKIgnoresTheOutputDirectoryFlag`

ה-handler קורא `arguments.outputdir` בעוד שה-parser מגדיר dest בשם `outputDirectory`. הפלט תמיד נוחת תחת תיקיית העבודה של התהליך. אומת בקוד.

### B177. `createReleaseRateFile` לא יכול לרוץ בכלל
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestSourceMakeEscapedMassFile`

`stochasticLagrangian_source_makeEscapedMassFile` קורא `toolkitHome.OF_LSM`, שאינו קבוע של `ToolkitHome` (רק מפתח פנימי ב-`_toolkits`) → `AttributeError`. בנוסף הוא מעביר את שם הפרויקט **פוזיציונלית** ל-`getToolkit(toolkitName, filesDirectory, ...)`, כך ש-`"tmpProject"` הופך לתיקיית הקבצים; ותת-הפרסר שלו לא מגדיר אף אחד משבעת ה-attributes שה-handler קורא.

### B178. חוסר ב-hermes הופך ל-`NameError`
**קובץ:** `hera/simulations/openFoam/CLI.py` · **מקובע ב:** `test_openfoam_cli.py::TestBuildExecuteWithoutHermes`

שני ה-handlers של buildExecute תופסים `ImportError` עם `warnings.warn` בלבד, ואז קוראים לשם שהייבוא הכושל היה אמור לקשור → `NameError: handler_buildExecute`. השגיאה שהמשתמש רואה לא מרמזת על התלות החסרה.

### B179. `toolkit_register` — משתמש בשמות שלא הוגדרו, והכשל נבלע
**קובץ:** `hera/utils/data/CLI.py:681-682` · **מקובע ב:** `test_utils_data_cli.py::TestToolkitRegisterIsBroken`

```python
th.registerToolkit(
    toolkit_name = toolkit_name,      # לא מוגדר — המשתנה המקומי הוא `name`
    toolkit_path = toolkit_path,      # לא מוגדר — המשתנה המקומי הוא `cls_path`
    ...
```

אף אחד משני השמות לא קיים בסקופ ואין גלובלים כאלה במודול. ה-`NameError` נבלע ב-`except Exception` של הפונקציה, שמדפיס `[ERROR] name 'toolkit_name' is not defined` ומחזיר 0 — כך שהפקודה **נכשלת בשקט** וכל פרסור ה-`--version`/`--repository` שמעליה הוא קוד מת. אומת בקוד.

### B180. `hera-project project version update` — קורא מתודה שלא קיימת על `Project`
**קובץ:** `hera/utils/data/CLI.py:487-488` · **מקובע ב:** `test_utils_data_cli.py::TestUpdateDatasourceDefaultVersion`

ה-handler בונה `hera.datalayer.Project` וקורא `proj.setDataSourceDefaultVersion(...)` — אבל המתודה מוגדרת על `hera.toolkit.abstractToolkit`, לא על `Project`, ול-`Project` אין `__getattr__` חלופי. אומת: `hasattr(Project, 'setDataSourceDefaultVersion')` הוא `False`. כלומר **אי אפשר לקבוע גרסת ברירת מחדל מה-CLI בכלל** — וזה סותר ישירות את ההוראה ב-`CLAUDE.md` לקרוא ל-`setDataSourceDefaultVersion` אחרי רישום datasource מגורסן.

### B181. `add_toolkit` — מדווח על JSON תקין ככשל פרסור
**קובץ:** `hera/utils/data/CLI.py:396-404` · **מקובע ב:** `test_utils_data_cli.py::TestAddToolkit`

ה-`raise ValueError("Params must be a JSON object.")` המפורש יושב בתוך אותו `try` שתופס שגיאות פרסור, וה-`except Exception` העירום זורק הכול מחדש בתור `Invalid JSON passed to --params: ...`. לכן `--params '[1,2]'` אומר למשתמש שה-JSON שלו שגוי כשהוא נפרס מצוין; ההודעה שהמחבר עצמו כתב בלתי-נגישה (היא שורדת רק כ-`__cause__`).

### B182. `hera-project ... --default` — דגל בלי `store_true`, ולכן בולע ערך
**קובץ:** `hera/bin/hera-project:118`

```python
parser_display.add_argument('--default', dest="default", default=False, help='If to show Default Versions only.')
```

חסר `action="store_true"`, כך שהארגומנט **צורך ערך**: `--default false` נותן את המחרוזת האמיתית-לוגית `"false"`, ואי אפשר לכבות את הדגל אחרי שהועבר. לא מקובע בבדיקה כי `hera/bin/hera-project` הוא סקריפט בלי סיומת שאינו ניתן לייבוא. אומת בקוד.

### B183. `_getUrbanRoughnessFromLandCover` — כותב **קוד סיווג** של כיסוי קרקע לתוך `z0`
**קובץ:** `hera/measurements/GIS/raster/landcover.py:673` · **מקובע ב:** `test_gis_landcover_more.py::TestUrbanRoughnessOutsideTheBuildings`

תא שלא חותך אף פוליגון lambda מקבל `landcover['z0'].values[i,j] = self.getLandCoverAtPoint(...)` — כלומר את קוד IGBP (0–16) — בעוד שכל ההשמות ה"אחיות" כותבות `lambdas['zz0']`. התוצאה: תאים כפריים בתוך תחום אורבני מקבלים z0 של 12–14 **מטר** במקום ~0.15 מ'. אומת: `z0` יצא זהה בדיוק ל-`landcover`.

### B184. טבלת z0 **שלישית** שלא מסכימה עם השתיים האחרות
**קובץ:** `hera/measurements/GIS/raster/landcover.py:407-412` · **מקובע ב:** `test_gis_landcover_more.py::TestGetRoughnessFromLandcover`

ה-lambda הפנימי של `getRoughnessFromLandcover`, שנמצא בנתיב ה"ברירת מחדל: IGBP" המתועד, לא מסכים לא עם `_handleType1` (Floors et al. 2021 המפורסם) ולא עם הרמפה ב-`getRoughnessAtPoint` שכבר קובעה כ-B43: מים 0.01 מול 0.0001 (פי 100), croplands 0.55 מול 0.15, barren 0.0001 מול 0.01. רק סיווג 15 מתלכד.

### B185. `getLandCover` — לא מטפל ב-CRS של הראסטר, בניגוד ל-`getLandCoverAtPoint`
**קובץ:** `hera/measurements/GIS/raster/landcover.py` · **מקובע ב:** `test_gis_landcover_more.py::TestGetLandCoverIgnoresTheRasterCRS`

הפונקציה קוראת `src.transform` אבל אף פעם לא משווה את `src.crs` ל-CRS הקלט. אומת מול ראסטר ב-EPSG:2039: `getLandCoverAtPoint` מחזיר 13 נכון, ואילו `getLandCover` על אותו תיבה מת ב-`IndexError: index 6568 is out of bounds for axis 0 with size 20`. ראסטר מוקרן שהאינדקסים שלו במקרה קטנים יחזיר בשקט סיווגים שגויים במקום לקרוס.

### B186. `presentation.plotLambdas` — מנרמל את סרגל הצבע על `z0` במקום על השדה המבוקש
**קובץ:** `hera/measurements/GIS/raster/landcover.py` · **מקובע ב:** `test_gis_landcover_more.py::TestPlotLambdasNormalisesOnZ0`

המלבנים נצבעים לפי `Normalize(landcover[field])`, אבל הסרגל עצמו הוא `Normalize(landcover.z0...)` ומתויג `f"{field} Value"`. אומת: עם `lambdaF` בטווח [0.05,0.4] ו-`z0` בטווח [1,5], ה-ylim של הסרגל יוצא (1.0, 5.0). בלי `z0` בכלל — כלומר המקרה הלא-אורבני — הקריאה מתה ב-`AttributeError: 'DataArray' object has no attribute 'z0'`.

### B187. `np.vectorize` בלי `otypes` — קוטע כל z0 למספר שלם
**קובץ:** `hera/measurements/GIS/raster/landcover.py:417` · **מקובע ב:** `test_gis_landcover_more.py::TestRoughnessDtypeFromTheHandler`

numpy גוזר את ה-dtype מהתא (0,0) בלבד. המילון של `_handleType1` מערבב טיפוסים — חמשת סיווגי היער הם ה-`int` `1`, כל השאר float — כך שרשת שהתא הראשון בה הוא יער מקבלת `int64` ומאפסת בשקט את כל השאר. אומת: `[[1,10],[10,10],[10,10]]` נותן `z0 == [[1,0],[0,0],[0,0]]`; עשב 0.12 → 0, אורבני 0.8 → 0. **z0 אפס הופך את פרופיל הרוח הלוגריתמי לסינגולרי.** התיקון הוא `1 → 1.0` בטבלה או `otypes=[float]`.

### B188. שמות אריחי SRTM — בלי padding לרוחב ובלי אותיות חצי-כדור
**קובץ:** `hera/measurements/GIS/raster/topography.py:124,178,188` · **מקובע ב:** `test_gis_raster_topography_more.py::TestSRTMTileNaming`

```python
filename = 'N'+str(int(lat))+'E'+str(int(long)).zfill(3)+'.hgt'
```

אומת (השם נקרא בחזרה מתוך ה-`FileNotFoundError`): רוחב 5.5 → `N5E007.hgt` (צריך `N05E007.hgt`); רוחב ‎−1.5 → `N-1E007.hgt` (צריך `S02E007.hgt`, ו-`int` קוטע לכיוון אפס במקום לרצף למטה); אורך ‎−2.5 → `N40E-02.hgt` (צריך `N40W003.hgt`). ישראל (29–33N/34–36E) היא הסיבה שזה מעולם לא נשך.

### B189. `getPointListElevation` — מקבל `inputCRS` ולא קורא אותו בכלל
**קובץ:** `hera/measurements/GIS/raster/topography.py` · **מקובע ב:** `test_gis_raster_topography_more.py::TestGetPointListElevationIgnoresInputCRS`

אין המרה, אין ולידציה, אין אזהרה. אומת: המרת (35.15625, 32.90625) ל-ITM והעברה עם `inputCRS=ITM` גורמת לקריאת מטרים כמעלות → `FileNotFoundError: N753355E212021.hgt not found`.

### B190. ההודעה "datasource לא נמצא" בלתי-נגישה
**קובץ:** `hera/measurements/GIS/raster/topography.py:128,197` · **מקובע ב:** `test_gis_raster_topography_more.py::TestMissingDefaultSRTM`

שתי מתודות הגובה עושות `dataSourceName = self.getConfig()['defaultSRTM']` ורק **אחר כך** `if dataSourceName is None: raise ValueError(...)`. פרויקט שלא הגדיר ברירת מחדל פשוט אין לו את המפתח, כך שהמנוי העירום זורק `KeyError('defaultSRTM')` קודם. השומר נורה רק אם הקונפיג שומר `None` במפורש. `.get('defaultSRTM')` היה מתקן. אומת בקוד.

### B191. `TopographyToolkit.create_xarray` — `minx`/`miny` הם (lat, lon) ל-WSG84 ו-(x, y) לכל השאר
**קובץ:** `hera/measurements/GIS/raster/topography.py` · **מקובע ב:** `test_gis_raster_topography_more.py::TestCreateXarrayArgumentOrder`

ענף ה-WSG84 קורא `self.convertPointsCRS(points=[[miny, minx]], ...)`, כלומר `minx`/`maxx` הם **רוחב** — ההפך משמם, ההפך מהחוזה (x,y)=(lon,lat) של `convertPointsCRS` באותה מחלקה עצמה, וההפך מענף ה-`else` של אותה מתודה (`Point(minx, miny)`, שם `minx` אכן מזרח). הדוקסטרינג אומר רק "bounding box coordinates". אומת: `create_xarray(minx=35.15, miny=32.90, ...)` מחזיר רשת ב-lon≈32 / lat≈35 — מדינה אחרת — בלי שום שגיאה. `getElevation` ו-`createElevationSTL` יורשים את המלכודת.

**הערה נוספת:** `roughnesslength2sandgrainroughness` מוגדרת **פעמיים** ב-`LandCoverToolkit` (שורות 561 ו-717) והשנייה מסתירה את הראשונה. שתיהן מחשבות `z0 * 30`, כך שההתנהגות זהה — קובע בבדיקה ברמת המקור כדי שעריכה עתידית של עותק אחד בלבד תיתפס. אומת בקוד.

## אצווה 31 — `openFoam/toolkit.py` ו-`lagrangian/LSM/toolkit.py` (B192-B211)

כל הפגמים למטה מקובעים ב-`test_openfoam_toolkit_more.py` / `test_openfoam_lagrangian_lsm_toolkit_more.py`, כל אחד עם `xfail(strict=True)` ובדיקת אפיון עוברת.

### `openFoam/toolkit.py`

### B192. `getMesh` — `-case <dir>` נדחף כאיבר **אחד** ברשימת argv
`f"-case {caseDirectory}"` נכנס כאלמנט בודד (ו-`""` כשה-case הוא ה-cwd), כך ש-`foamJob` מקבל ארגומנט אחד בלתי-פרסבילי. אומת: `['foamJob', '-case /tmp/...', '-wait', ...]`.

### B193. `getMeshFromName`/`getMeshExtentFromName` — מתעדים פרמטרים ולא מעבירים אותם
שניהם מתעדים `readParallel` ו-`time` אבל קוראים `self.getMesh(doc.getData())` פוזיציונלית, כך ששניהם נבלעים בשקט (ל-`getMeshExtent` אין בכלל פרמטרים כאלה).

### B194. `prepareSlurmWorkflowExecution` — `UnboundLocalError` על הטיפוס המתועד
`workflow` נקשר רק בענף `isinstance(baseConfiguration, str)`, כך ש-workflow כ-**dict** (הטיפוס המתועד) מת ב-`UnboundLocalError: workflow`. ענף הטיפוס הלא-חוקי רק עושה `logger.error` וממשיך ליפול (`baseConfiguration=12345` → `TypeError: 'int' object is not iterable` מתוך `JSONVariations`).

### B195. `prepareSlurmWorkflowExecution` — בודק `jsonVariations` כ-dict בעוד ש-`JSONVariations` רוצה list
ה-dict שמתקבל מעלה `AttributeError: 'str' object has no attribute 'items'`, בזמן שצורת ה-list שכן עובדת היא בדיוק זו שהקוד מדווח עליה כשגיאה.

### B196. `Analysis.getFiltersDocuments` — `is None` מול API שמחזיר רשימה
`getWorkflowDocumentFromDB` מחזיר list, כך שסימולציה לא מוכרת נותנת `TypeError: list indices must be integers` ולא ה-`ValueError` המיועד.

### B197. `Presentation.toUnstructuredVTK` — שומר ה-overwrite בודק שם קובץ שגוי
הבדיקה היא `os.path.exists("<base>_0")` בעוד ש-evtk כותב `<base>_0.vtu` — כך שגם ה-`FileExistsError` וגם המחיקה ב-`overwrite=True` לא יכולים לפעול על פלט אמיתי.

### B198. `Presentation.toStructuredVTK` — קורא למתודה שקיימת רק ב-Analysis של ה-LSM
`self.analysis.calcConcentrationTimeStepFullMesh` לא קיימת על ה-Analysis של OFToolkit → `AttributeError` לכל קלט.

### B199. `Presentation.loadLagrangianDataParallel` — מפיל את צעד הזמן הראשון
הרשומות נבנות מ-`timeList[1:]` גם כאשר `times` הועבר במפורש; `times=["100"]` נותן `TypeError: Must supply at least one delayed object` מ-dask.

### B200. אותה מתודה — צעדי זמן מספריים קורסים
`times` עובר `numpy.atleast_1d` ואז ישר ל-`os.path.join` → `TypeError: join() argument must be str…`. רק מחרוזות עובדות.

### `lagrangian/LSM/toolkit.py`

### B201. `_extractFile`/`_readRecord` — משתמשים ב-`self.logger` שלא קיים
**קובץ:** שורות 137, 237, 241. לא `abstractToolkit` ולא `Project` מגדירים `logger` (הם בונים loggers מקומיים דרך `get_classMethod_logger`), כך ש**כל** קריאה ל-`_readRecord` נגמרת ב-`AttributeError: 'OFLSMToolkit' object has no attribute 'logger'` — עם נתונים או בלעדיהם. מאותה סיבה גם הפרסר החלופי הכתוב-ביד של OpenFOAM ב-`_extractFile` בלתי-נגיש. (שורה 137 גם קוראת `self.logger.execute(...)` — שאינה מתודה של logger בכלל.) אומת בקוד.

### B202. `loadData` — שומר ה-saveMode קורס במקום לדווח
ההודעה היא `",".join([...])` על רשימה שמכילה את `TOOLKIT_SAVEMODE_NOSAVE`, שהוא `None` → `TypeError: sequence item 0: expected str instance` במקום ה-`ValueError` המיועד. אותה רשימה גם **משמיטה** את `TOOLKIT_SAVEMODE_ONLYFILE` (וכופלת את `FILEANDDB`), כך שהמצב המתועד הזה נדחה למרות שהגוף מטה מטפל בו.

### B203. `loadData` — ענף המעבד הבודד מסנן על שם חשוף
`os.path.isdir(x)` על השם בלי הנתיב (אותו שורש כמו B83 שכבר תועד) → אין צעדי זמן → `TypeError: Must supply at least one delayed object`.

### B204. `loadData` — ענף ה-DB: QuerySet לעולם אינו `None`
`doc = self.getSimulationsDocuments(...)` הוא QuerySet, כך ש-`FILEANDDB` **תמיד** מעלה `FileExistsError("Data already in the DB")` (ו-`addSimulationsDocument` הוא קוד מת — אף מסמך `LSMRuns` לא נכתב), ואילו `FILEANDDB_REPLACE` נופל ב-`AttributeError: 'QuerySet' object has no attribute 'save'`. ה-parquet נכתב לפני זה, כך שהקובץ כן קיים.

### B205. `loadData` עם `ONLYFILE_REPLACE` — `UnboundLocalError` על `doc`
הקובץ נכתב ואז `return doc` קורס: `doc` נקשר רק בענפי ה-DB וב-`else` של NOSAVE.

### B206. `makeSource` — כותרת מקובעת שמתעלמת מ-`fileName`
ה-header הוא מחרוזת קשיחה שבה `object` הוא תמיד `kinematicCloudPositions`, גם כשה-`fileName` שונה.

### B207. `createRootCaseMeshLink` — `NameError` על `lastTS`
**קובץ:** שורה 717 משתמשת ב-`lastTS` שאף פעם לא הוקצה (מופע יחיד בקובץ, בלי השמה) → `NameError` לכל root case שיש בו תיקיית processor. בנוסף צעד היעד מקובע ל-`"3600"`. אומת בקוד.

### B208. `calcConcentrationTimeStepFullMesh` — השמת attribute על `DataArray`
`fulldata.filterType = "C"` — xarray מסרב להשמת attribute → `AttributeError` בכל קריאה. זה מפיל גם את `calcConcentrationFieldFullMesh`, שרושמת את מסמך ה-cache ואת ה-resource `Concentrations*.nc` **לפני** הכשל, ומשאירה מסמך שמצביע על קבצים שמעולם לא נכתבו — וגם את שני המימושים של `toStructuredVTK`. רדום מאחורי זה: `transpose("yI","xI","zI","time ")` עם רווח נגרר בשם המימד.

### B209. `calcDocumentConcentrationPointWise` — מתעד `:return: Document` ואין `return`
**קובץ:** שורה 946: `ret = doc` ואז סוף הפונקציה → תמיד `None`, למרות שה-parquet ומסמך ה-cache כן נוצרים כמו שצריך. אומת בקוד.

### B210. `getConcentrationField(returnFirst=False)` — מתעד רשימה, מחזיר `docList[0]`
מתועד שיחזיר רשימה שאולי ריקה, אבל מחזיר `docList[0].getData(usePandas=True)` → `IndexError` כשלא נמצא כלום.

### B211. `_extractVelocityField` — לא יכול לקרוא קובץ שדה מציאותי
קורא את כל קובץ ה-U עם `pandas.read_csv(names=['u','v','w'])` וחותך ל-`nCells` רק **אחרי כן**, כך שכל שורת boundary עם יותר משלושה טוקנים (למשל `value uniform (0 0 0);`) מעלה `pandas.errors.ParserError: Expected 3 fields`. התיקון הוא `nrows=nCells`.

## אצווה 32 — `experiment/analysis.py` ו-`experiment/presentation.py` (B212-B226)

מקובע ב-`test_experiment_analysis.py` / `test_experiment_presentation_plots.py`, כל אחד עם `xfail(strict=True)` ובדיקת אפיון עוברת. שני הקבצים בונים **ניסוי argos v3.0.0 אמיתי** ב-`tmp_path` (בלי שום stub). התוצאה המרכזית: מתוך 13 המתודות שנבדקו, **רק אחת** (`plotDevicesOnImage`, ובפריסת מכשירים אחת ספציפית) רצה בכלל — כל השאר מתות בהוראה הראשונה או השנייה.

### B212. `entitiesTable` נקראת כפונקציה למרות שהיא property (2 מקומות)
**קובץ:** `hera/measurements/experiment/analysis.py:52,280` · `argos/experimentSetup/dataObjects.py:208`

`getDeviceLocations` ו-`addMetadata` כותבות `...[trialName].entitiesTable()`, אבל ב-argos זה `@property` שמחזיר DataFrame → `TypeError: 'DataFrame' object is not callable`. אומת בקוד, זה מול זה.

### B213. `getDeviceLocations` — מסננת על עמודה שלא קיימת בטבלה של trial
**קובץ:** `hera/measurements/experiment/analysis.py:52`

השאילתה היא `entityType==@entityTypeName`, אבל טבלת ה-trial קוראת לעמודה `deviceTypeName`; `entityType` קיימת רק בטבלאות של `Experiment`/`EntityType` → `pandas.errors.UndefinedVariableError`.

### B214. `addMetadata` — ממזגת על `entityName` במקום `deviceItemName`
אותו חוסר-התאמה בדיוק כמו B162, מועתק לתוך `analysis.py`.

### B215. `getDeviceTypePlannedMessageCount` — קוראת ל-`getOptimalFrequencyHz` שלא קיימת בשום מקום
**קובץ:** `hera/measurements/experiment/analysis.py:249,252`

חיפוש בכל `hera/` מחזיר **רק את שני אתרי הקריאה** ואף הגדרה → `AttributeError`. זה הורג גם את `normalize=True`, שהוא **ברירת המחדל** של `getDeviceTypeTransmissionFrequencyOfTrial`. אומת בקוד.

### B216. `TrialSet.trials` הוא dict, לא DataFrame
ענף החישוב-מחדש עושה `.trials.query("trialName == @trialName")` → `AttributeError: 'dict' object has no attribute 'query'`. גם ל-`trialsTable` (ה-DataFrame האמיתי) אין עמודת `trialName`. מכיוון ש-`recalculate` ברירת המחדל היא `True`, שום דבר לא יכול למלא את ה-cache. מאחורי זה: `freq.set_index('timestamp')` שהתוצאה שלו נזרקת, ו-`resample("1min")` שמקבע את הבין במקום להשתמש ב-`samplingWindow`.

### B217. `addTrialProperties` — חיסור זמן שחרור בלי תנאי
שולף `ReleaseStart` עם `.get(...,None)` ואז מחשב `tmp.timestamp - releaseTime` → `TypeError` לכל trial בלי שחרור, ואיתו נאבדות גם עמודות ה-`fromStart`.

### B218. `plotImage` — מפתחות extent שגויים
קורא `metadata['xleft'/'xright'/'ybottom'/'ytop']`, בעוד ש-`getImageMetadata` של argos מחזיר `left`/`right`/`lower`/`upper` → `KeyError: 'xleft'`, אחרי שהתמונה כבר פוענחה.

### B219. `plotMap` — `self.trialSet` שיושב על ה-datalayer
`AttributeError`. בנוסף: `ax`/`plot_kwargs` נקראים למרות שאף אחד מהם אינו פרמטר, `deviceType`/`toolkitDataSource` לא מוגדרים, ואין `return`.

### B220. `_plotEntityLocationScatter`/`_plotEntityLocationNames` — `datalayer.experimentSetup` לא קיים
שתיהן קוראות `self.datalayer.experimentSetup`, ול-`experimentSetupWithData` אין attribute כזה — **הוא עצמו** ה-setup (אותה טעות בדיוק כמו B165) → `AttributeError`. מאחורי זה: `entitiesTable(status)`, `self._entityMarkers` חסר, ו-`FLOOR_PLATFORM`/`FLOOR_CONCOURSE` לא מוגדרים.

### B221. `plotDevices` — `UnboundLocalError` על שם פרמטר שגוי
`plot_kwargs = plot_kwargs or {}` קורא משתנה מקומי לפני השמה — הפרמטר נקרא `plotkwargs`. אם מעבירים `ax`, הוא נופל במקום זה על `self.trialSet`, ואחר כך על `_process_row` (B160) ועל `row.stationName`.

### B222. `plotDevicesOnImage` — mutable default שנכתב אליו
**קובץ:** `hera/measurements/experiment/presentation.py:395`

`scatterkwargs.setdefault("s",50)`/`("c","r")` על ברירת מחדל `scatterkwargs={}` כותב **לתוך ה-dict המשותף** לתמיד, כך שאף קורא מאוחר יותר לא יכול לקבל scatter בלי עיצוב על ידי השמטת הארגומנט. אומת בקוד.

### B223. `plotDevicesOnImage` — קורא עמודה אופציונלית בלי תנאי
`row.containedIn` נקרא תמיד, אבל argos יוצר את העמודה הזו רק כשמכשיר מצהיר על הכלה → `AttributeError` לניסוי הרגיל שממוקם על מפה, למרות שלמתודה **כבר יש** ענף חלופי בשבילה.

### B224. `plotDevicesOnImage` — אגרגציה על כל ה-DataFrame
`devices_df.max().longitude` / `.min().latitude` מבצעים אגרגציה על **כל** העמודות; NaN בעמודת המחרוזות `containedIn` (trial מעורב) → `TypeError: '>=' not supported between instances of 'str' and 'float'`. יחד עם B223 זה משאיר בדיוק פריסה אחת שעובדת.

### B225. `plotDevicesOnImage` — הצירים מוחלפים
מפזר `x=row.latitude, y=row.longitude`, בעוד ש-`plotImage` שם קווי אורך על x (`left`/`right`) ורוחב על y (`lower`/`upper`). היסט התווית מחמיר את זה בכך שהוא מזיז קו אורך לפי מרווח הרוחב. אומת: ההיסטים יוצאים `(32.0, 34.0)` במקום `(34.0, 32.0)`.

### B226. `plotDeviceTypeFunctionality` — מת כתוצאה מ-B215/B216
קורא למתודת התדירות עם `normalize=True` ובלי `recalculate=False`, כך שהוא תמיד נכנס לענף המת של B216 ואז פוגע ב-B215. מקובע כאן כתוצאה, לא כפגם עצמאי.

## אצווה 32 — GIS vector/raster (B227-B242)

מקובע ב-`test_gis_vector_topography.py`, `test_gis_demography_more.py`, `test_gis_buildings_toolkit_more.py`, `test_gis_utils_more.py`, `test_gis_raster_topography_statistics.py`, `test_gis_buildings_analysis_more.py`.

### B227. `vector/topography.py::cutRegionFromSource` — קורא ל-super עם חתימה אחרת לגמרי
**קובץ:** `hera/measurements/GIS/vector/topography.py:77` מול `vector/toolkit.py:98`

```python
# הילד:
super().cutRegionFromSource(shapeDataOrName=..., datasourceName=..., isBounds=..., crs=...)
# ההורה:
def cutRegionFromSource(self, datasourceDocument, shape, isBounds=False, inputCRS=WGS84)
```

אף אחד מארבעת שמות הפרמטרים לא תואם → כל קריאה נופלת ב-`TypeError: ... unexpected keyword argument 'shapeDataOrName'`. זה מפיל איתו גם את `regionToSTL`. אומת בקוד, זה מול זה.

### B228. `geoPandasToSTL` — מעביר את הליטרל במקום את הפרמטר
**קובץ:** `hera/measurements/GIS/vector/topography.py:100,114`

הפרמטר `solidName="Topography"` מוגדר בחתימה, אבל השורה המחזירה כותבת `solidName="Topography"` כליטרל. `geoPandasToSTL(c, solidName="MyHill")` מייצר `solid Topography`. אומת בקוד.

### B229. `toDEM` — `getDatasourceData` בשגיאת כתיב
**קובץ:** `hera/measurements/GIS/vector/topography.py:180`

המתודה האמיתית היא `getDataSourceData` (S גדולה). `toDEM("CONTOURS")` → `AttributeError`, וזה גם הופך את נתיב ה-fallback של geoJSON-string לקוד מת. אומת בקוד.

### B230. `toDEM` — אינדקסים מוחלפים בין הרשת לגובה
**קובץ:** `hera/measurements/GIS/vector/topography.py:213-215`

`Nx=grid_x.shape[0]`, `Ny=shape[1]`, ואז `height[j, i]` כאשר `i∈Nx` ו-`j∈Ny` — המנויים הפוכים. קווי גובה של 800 מ' × 400 מ' ב-`dxdy=200` נותנים `IndexError`; רשתות מרובעות פשוט משתטחות בטרנספוזיציה בשקט.

### B231. `analysis.addHeight` — מייבא את המודול ומצפה למחלקה
**קובץ:** `hera/measurements/GIS/vector/topography.py:14,262`

`from hera.simulations.utils import coordinateHandler` קושר את **המודול** (כי `hera/simulations/utils/__init__.py` ריק), בעוד שהמחלקה נמצאת בתוכו. כל קריאה → `AttributeError: module ... has no attribute 'regularizeTimeSteps'`. אומת בקוד.

### B232. `analysis.addHeight` — `toolkit` נדרס, ולכן **כל** מצב שמירה נכשל
**קובץ:** `hera/measurements/GIS/vector/topography.py:13,282`

שורה 13 עושה `from hera.measurements.GIS.vector import toolkit`, כלומר השם `toolkit` מצביע על מודול ה-toolkit הווקטורי — ואז שורה 282 קוראת `toolkit.TOOLKIT_SAVEMODE_FILEANDDB`. רשימת המצבים נבנית **לפני** בדיקת החברות, כך שגם `NOSAVE` נכשל; `addHeight` לא יכול להחזיר בשום מצב. אומת בקוד. (מקובע בבדיקה מפורמטרת על חמשת מצבי השמירה.)

### B233. `demography.createNewArea` — אותה דריסת `toolkit`, אחרי שהקובץ כבר נכתב
**קובץ:** `hera/measurements/GIS/vector/demography.py:9-10,305`

`saveMode=FILEANDDB` → `AttributeError: ... no attribute 'TOOLKIT_DATASOURCE_NAME'`, **אחרי** שהקובץ נכתב לדיסק — כלומר נשאר קובץ יתום בלי מסמך.

### B234. `createNewArea` — `FILEANDDB_REPLACE` לא נוגע ב-DB בכלל
**קובץ:** `hera/measurements/GIS/vector/demography.py:303`

ענף ה-DB מוגן ב-`== TOOLKIT_SAVEMODE_FILEANDDB` בלבד, כך ש-`FILEANDDB_REPLACE` כותב את הקובץ ולא מעדכן מסמך — בניגוד לדוקסטרינג שלו.

### B235. `createNewArea` — `ONLYFILE` מתעד שגיאה על קובץ קיים ולא בודק כלום
**קובץ:** `hera/measurements/GIS/vector/demography.py:298-301`

הדוקסטרינג אומר "raise exception if file exists", אבל אין שום `os.path.exists`, כך שהוא דורס בשקט ואי אפשר להבדיל אותו מ-`ONLYFILE_REPLACE`.

### B236. `getBuildingHeightFromRasterTopographyToolkit` — מקבל `topographyDataSource` ולא מעביר אותו
**קובץ:** `hera/measurements/GIS/vector/buildings/toolkit.py:47-65`

הפרמטר מתועד ומקובל, אבל לא מועבר ל-`getPointListElevation`. עם שני ראסטרים רשומים, בקשה למי-שאינו-ברירת-המחדל מחזירה את ערכי ברירת המחדל, ושם שלא קיים מתקבל בשקט.

### B237. `buildingsGeopandasToSTLRasterTopography` — בניינים בשטח שטוח יוצאים גבוהים מדי
**קובץ:** `hera/measurements/GIS/vector/buildings/toolkit.py:129,150`

`LengthFwd = wallsheight + nonFlatTopographyShift` מחושב ללא תנאי, אבל ענף `flatTerrain` ממקם את המסגרת ב-`referenceTopography` **בלי** להנמיך אותה בהיסט — כך שבניין של 10 מ' עם `nonFlatTopographyShift=10` יוצא עם גג ב-20 מ'.

### B238. ייבוא FreeCAD בבלוק אחד — השומר של המודול נעקף
**קובץ:** `hera/measurements/GIS/vector/buildings/toolkit.py:8-13`

`try: import FreeCAD; import Part; import Mesh` בבלוק אחד: התקנה חלקית משאירה את `FreeCAD` קשור אבל לא את `Part`/`Mesh`, כך שההודעה שהקוד עצמו כתב ("FreeCAD not found. Install before using…") נעקפת והמתודה מתה מאוחר יותר על `NameError: name 'Part' is not defined`.

### B239. `stlFactory.rasterizePandas` — צירי הקואורדינטות וגריד הגבהים בגדלים שונים
**קובץ:** `hera/measurements/GIS/utils.py:414-417`

הצירים באים מ-`numpy.mgrid[xmin:xmax:dxdy]` (מספר נקודות מעוגל **למעלה**) בעוד שגריד הגבהים נבנה לפי `Nx=int((xmax-xmin)/dxdy)` (**נקטע**). מרווח של 100 מ' ב-`dxdy=30` נותן `x.shape==(4,4)` מול `height.shape==(3,3)`, ואז `rasterToSTL` פולט בשקט גוף קטן מדי (אין אף קודקוד מעבר ל-x=60 למרות שהציר מגיע ל-90).

### B240. `rasterizePandas` — גם כשהגדלים מסתדרים, הקואורדינטות שגויות
**קובץ:** `hera/measurements/GIS/utils.py:414-417`

`mgrid` עם **צעד** מחריג את נקודת הקצה, בעוד ש-`regularizeTimeSteps` מבצע אינטרפולציה על `mgrid[...:complex(n)]` שכן **כולל** אותה. גובה שנדגם ב-x=100 מתויג כ-x=50 — שגיאה של `dxdy` שלם בכל תא מלבד הראשון.

### B241. `topographyAnalysis.calculateStastics` — שמות משתנים לא תואמים למפיק היחיד שלהם
**קובץ:** `hera/measurements/GIS/raster/topography.py:561-573`

קורא `['X']`/`['Y']`/`['Elevation']`, אבל `getElevation` — המפיק היחיד של גבהים ב-toolkit — מחזיר `lat`/`lon` ו-coord בשם `elevation`. הזנה של האחד לשני נותנת `KeyError: 'X'`. שינוי שם שלושת המשתנים גורם לאותו dataset לעבוד, כלומר הפגם הוא בדיוק בגבול השמות. גם הדוגמה בדוקסטרינג בלתי שמישה: היא קוראת ל-`tk.getDomainElevation()` שלא קיים בשום מקום.

### B242. `LambdaFromBuildingData` — ה-cache הוא write-only בגלל `numpy.ndarray` במתאר
**קובץ:** `hera/measurements/GIS/vector/buildings/analysis.py:139-147`

מתאר ה-cache נושא `"bounds": data.total_bounds`, שהוא `numpy.ndarray`. הכתיבה מצליחה, אבל ה**חיפוש** משווה את המערך בהקשר בוליאני → `ValueError: The truth value of an array with more than one element is ambiguous`. הקריאה הראשונה (אוסף ריק) עובדת, ומאז **כל** קריאה עם `saveCache=True` נכשלת — כולל `overwrite=True`, כי החיפוש קודם לבדיקת הדריסה. גם ענף "החזר נתונים שנמצאו" וגם ענף "עדכן רשומה ישנה", וגם עצת ה-`FileNotFoundError` שבדוקסטרינג — כולם קוד מת. אותו מפתח עם `bounds` כרשימה רגילה נשאל בלי בעיה.

## אצווה 32 — זנב הפיזיקה/utils (B243-B258)

מקובע ב-`test_simulations_nearwallflow_more.py`, `test_evaporation_models_more.py`, `test_simulations_canopywindprofile.py`, `test_wrf_datalayer.py`, `test_deposition_models_more.py`, `test_windprofile_toolkit_more.py`, `test_utils_freecad.py`.

### B243. `nearWallFlow.Cplus` — הענף הגס-לחלוטין מחזיר קבוע של חוק אחר
**קובץ:** `hera/simulations/hydrodynamics/nearWallFlow.py`

הענף מחזיר `8` חשוף — קבוע Nikuradse של חוק `ln(y/k)` — במקום `8 − ln(k⁺)/κ`, כשהחוק כאן נמצא ביחידות צמיגות. C⁺ מקפץ מ-‎−2.47 ל-‎+8 סביב k⁺=70, כך ש-`channelFlow(Ra=2e-3).get_Ucenter_from_Ustar(0.3)` נותן 9.93 מ'/ש' לעומת 9.03 מ'/ש' בקיר חלק — כלומר **חספוס הקיר מוסיף תנע**.

**חשוב — יש כאן חוסר הסכמה בתוך מערך הבדיקות עצמו:** בדיקה קיימת, `test_simulations_nearwallflow.py::test_the_fully_rough_limit_is_eight`, מתעדת את הרמה הזו כהתנהגות **מכוונת**. שתי הבדיקות עוברות יחד (הישנה כאפיון, החדשה כ-xfail), אבל צריך הכרעה אנושית מי מהן צודקת. אני מציג את שתי הטענות ולא מכריע.

### B244. `canopyWindProfile.calcU` — גובה החופה עצמו נופל לענף השגוי
**קובץ:** `hera/simulations/utils/canopyWindProfile.py`

`if z < hc` ו-`elif z > hc` שניהם ממשיים, כך ש-`z == hc` נופל ל-`else` — שדורס את z ב-2000 מ'. בגובה החופה מתקבל 11.10 מ'/ש' במקום 2.44 מ'/ש' ששני הענפים בנויים לתת (פי 4.55), ערך זהה לזה של z=3000 מ'.

### B245. `urbanLogExponentProfile` — כתיבת האינטרפולציה היא chained assignment
**קובץ:** `hera/simulations/utils/canopyWindProfile.py`

`data['Ux'].iloc[i] = …` — כבר עכשיו מפיק `FutureWarning: ChainedAssignmentError`, ותחת `mode.copy_on_write` (ברירת המחדל בפנדס 3) 60 מתוך 80 תאים נשארים NaN: השלב כולו הופך ל-no-op שקט.

### B246. `wrfDatalayer.find_i` — `UnboundLocalError` על רוב תחום WRF
**קובץ:** `hera/simulations/WRF/wrfDatalayer.py`

`delta`/`request_delta` נקשרים רק בענף המתאים אבל מוחזרים ללא תנאי, כך שבקשה שנופלת ברווח שבין שתי עמודות מוטות (כלומר רוב תחום WRF אמיתי) מעלה `UnboundLocalError`.

### B247. `wrfDatalayer.find_i` — `isel(dim=i+1)` בלי בדיקת גבול
השורה/עמודה האחרונה של **כל** תחום אינה ניתנת לכתובת: בקשה בעמודה האחרונה מעלה `IndexError: index 4 is out of bounds`.

### B248. `wrfDatalayer` — `import geopandas` מסומן כהערה אבל geopandas בשימוש
**קובץ:** `hera/simulations/WRF/wrfDatalayer.py:8` מול `:84,94,150`

שורה 8 היא `#import geopandas`, ו-`getPandas` משתמש בו בשלושה מקומות → `NameError: name 'geopandas' is not defined`. גם ה-GeoDataFrame המסיים מת, כך שהפונקציה לא יכולה להחזיר כלום. אומת בקוד.

### B249. `wrfDatalayer` — כשל `import wrf` מודפס בלבד
ה-`ImportError` רק עושה `print`, ואז `getPandas` מת באמצע על `NameError: name 'wrf' is not defined` — בלי לרמז על התלות החסרה. (`python-wrf` גם אינו ב-`requirements.txt`.)

### B250. `getSpatialWind` — מעביר רזולוציה למקום של כיוון רוח
**קובץ:** `hera/simulations/windProfile/toolkit.py`

קורא `getRoughnessFromLandcover(xarray, dxdy)`, אבל הפרמטר השני שם הוא `windMeteorologicalDirection`. מרווח הרשת נצרך ככיוון רוח (30°) ו-`resolution` נשאר `None`.

### B251. `_getWindSpeedDirection` — `max()` על רשימה ריקה
`max(datetime_objects)` על רשימה שמתמלאת רק עבור תחנות שימושיות → `ValueError: max() iterable argument is empty` כשכל התחנות מנותקות או מחוץ לטווח.

### B252. `_getWindSpeedDirection` — `ws`/`wd` נקשרים רק בתוך לולאת הערוצים
תחנה עם גשם בלבד מעלה `UnboundLocalError`; וכשחסר ערוץ אחד בלבד, הערך של **התחנה הקודמת** נעשה בו שימוש בשקט.

### B253. `freeCAD.getObjFileBoundaries` — כשל ייבוא שמתדרדר ל-`NameError`
**קובץ:** `hera/utils/freeCAD.py`

כשל הייבוא מומתן ל-`warnings.warn` והשמות נשארים לא מוגנים → `NameError: name 'Mesh' is not defined`, שלא מזכיר לא את FreeCAD ולא את הפתרון. (ה-conftest מסטב את FreeCAD אבל לא את `Mesh`, וזה בדיוק המצב של התקנה חלקית.)

### B254. `depositionModels.__init__` — `[0]` על תוצאת חיפוש לא מסוננת
`p.getCacheDocuments(type="surface", surface=surface)[0]` בלי בדיקה → שם משטח שלא נרשם מעלה `IndexError: list index out of range` במקום הודעה מובנת.

### B255. `depositionRate_Petroff` — המרת יחידות כפולה של קוטר החלקיק
`dpm = 0.000001*dp` ממיר שוב קוטר שהבנאי כבר שם במטרים (כל הקבועים השכנים הם SI), כך שחלקיק של 1 מיקרון מדומה כ-1 **פיקומטר**: `vds/u*` יוצא `8.5e4` עבור גודל שהוא **יעילות** שיקוע (≤ 1, בסדר גודל 1e-5 מצופה).

### B256. `evaporationModels.agent` setter — קורא מתודת מופע על המחלקה
`RiskToolkit.getAgent(newAgent)` על המחלקה → `TypeError: missing 1 required positional argument: 'nameOrDesc'`; אי אפשר להחליף agent אחרי הבנייה. בנוסף: ה-getter של `agent` מוגדר **פעמיים** והראשון מת.

### B257. `depositionRate_Petroff` — קבוע בולצמן עם ספרות מוחלפות
**קובץ:** `hera/simulations/deposition/models.py:100`

```python
kb = 1.83E-23     # הנכון: 1.380649e-23
```

הספרות הוחלפו (1.83 מול 1.38), כך שהדיפוזיביות הבראונית גדולה ב-**32.6%** בכל מקום שבו מודל Petroff רץ. אומת בקוד.

### B258. `_getWindSpeedDirection` — לולאת הניסיונות מאמתת אחרי ההשמה
`data` מוקצה לפני המנוי שמאמת אותו, כך שמטען שגיאה שורד 20 ניסיונות, עובר את `if data:`, ומפוענח מחוץ ל-`try` → `KeyError: 'data'` על טוקן IMS שנדחה.

## אצווה 32 — `ToolkitHome`, datahandler, וזנב ה-risk (B259-B266)

### B259. `ToolkitHome.auto_register_and_get` — קורא ל-`registerToolkit` בחתימה שגויה
**קובץ:** `hera/toolkit.py`

קורא `registerToolkit(toolkitclass=…, datasource_name=…, repositoryName=…)`, בעוד שהחתימה היא `(toolkit_name, toolkit_path, params, version, overwrite, **kwargs)` → `TypeError: registerToolkit() missing 2 required positional arguments`. המתודה לא יכולה להחזיר בשום מצב. (זהו אחד מהיעדים שכמה אצוות דחו כ"לוגיקת רג'יסטרי מורכבת" — מסתבר שהוא פשוט מת.)

### B260. `ToolkitHome.import_toolkits_from_json` — אותו חוסר-התאמה בחתימה
אי אפשר לרשום שום toolkit מ-JSON; אותו `TypeError` בדיוק.

### B261. `AbstractCollection.addDocumentFromJSON` — no-op שקט לכל JSON שנושא `_id`
**קובץ:** `hera/datalayer/collection.py`

`mongoengine.Document.from_json` מוגדר `(json_data, created=False, **kwargs)`, ולפי הדוקסטרינג של mongoengine עצמו `created=False` אומר "אם יש id, הנח שכבר נשמר". ה-`.save()` שאחריו מנפיק update עם change set ריק — שום דבר לא נכתב, שום דבר לא נכשל, אין שגיאה. מכיוון ש-`asDict(with_id=True)`, `to_json()` ו-`Project.exportProject` **כולם** פולטים `_id`, ייבוא חוזר של מסמכים שיוצאו מפיל את כולם בשקט. התיקון: `from_json(json_data, created=True)`. אומת: `inspect.signature(mongoengine.Document.from_json)` → `(json_data, created=False, **kwargs)`.

**זה גם סוגר שאלה פתוחה:** הדוקסטרינג של `test_datalayer_project_export.py` תלה את זה ב"התנהגות לא עקבית של mongomock". זו לא mongomock — זו התנהגות מתועדת של mongoengine. ההערה שם עודכנה בהתאם.

### B262. `DataHandler_zarr_xarray.getData` — קורא `.attrs` מהנתיב במקום מה-dataset
**קובץ:** `hera/datalayer/datahandler.py`

`resource.attrs = JSONToConfiguration(resource.attrs)` כאשר `resource` הוא מחרוזת הנתיב ולא ה-dataset שנטען → `AttributeError: 'str' object has no attribute 'attrs'`. צריך להיות `df.attrs`, כמו ב-`DataHandler_netcdf_xarray.getData`. הקורא של zarr לא יכול להחזיר כלום.

### B263. `riskAreaAlgorithm_Sweep` — לא יכול לרוץ בכלל
**קובץ:** `hera/riskassessment/analysis/riskAreas.py`

`_doCalculation` קורא `effectIsopleths.datalayer(demog, releaseLoc, mathematical_angle=…)`, אבל ל-`thresholdGeoDataFrame` אין `datalayer` — המתודה עם בדיוק החתימה הזו נקראת `project`. כל נקודת שחרור מעלה `AttributeError`, ולכן `calculate()` לעולם לא מסתיים ו**האלגוריתם היחיד שמומש לאזורי סיכון בלתי שמיש**.

### B264. `InjuryLevel.__str__` — בלי `return`
**קובץ:** `hera/riskassessment/agents/effects/InjuryLevel.py:172`

`json.dumps(self.toJSON(),indent=4)` והתוצאה נזרקת → `TypeError: __str__ returned non-string`. מתודה נפרדת מ-B122 (`Injury.__str__`), ואף תת-מחלקה לא דורסת אותה, כך שכל ארבע מחלקות הרמות בלתי-ניתנות להדפסה. אומת בקוד.

### B265. `plotCasualtiesProjection` — שתי דרגות חומרה ומעלה לא מציירות כלום
**קובץ:** `hera/riskassessment/presentation/casualtiesFigs.py`

`query("severity in %s" % numpy.atleast_1d(plumSeverity))` מייצר `"severity in ['Severe' 'Light']"` — ה-`str()` של numpy משמיט פסיקים, ואיחוי מחרוזות צמודות של פייתון מקפל את זה ל-`['SevereLight']`, שלא תואם כלום. דרגה אחת עובדת במקרה. התיקון: להשתמש ב-`list(...)` או במשתני `@`.

### B266. `_findBoundingBox` — נשבר בדיוק במקרה של פגיעה בדרגה אחת
**קובץ:** `hera/riskassessment/analysis/riskAreas.py`

`set_index("datetime").loc[maxdatetime]` מחזיר **Series** כשרק שורה אחת נושאת את חותמת הזמן האחרונה, והשורה הבאה קוראת `.set_geometry` עליו → `AttributeError: 'Series' object has no attribute 'set_geometry'`. שורה אחת לצעד הזמן האחרון היא המצב **הרגיל** לפגיעה חד-רמתית. התיקון: `.loc[[maxdatetime]]`.

## אצווה 32 — LSM ו-gaussian (B267-B286)

### B267. `_LoggingShim.get_logger` — שובר בדיוק את מה שהוא נועד לשמר
**קובץ:** `hera/simulations/LSM/CLI.py` · ה-shim מעביר ל-`get_logger(instance, name=None)` של hera, שמתייחס לארגומנט הראשון כ-**מופע**. אומת: `CLI.logging.get_logger("hera.bin.hera_lsm.load_template").name` → `"builtins.str"`. ה-shim קיים כדי שקריאות לפי שם ימשיכו לעבוד — והוא שובר בדיוק אותן.

### B268. `setup_template` — השומר "symlink קיים ומצביע למקום אחר" הוא קוד מת
`os.unlink(dst)` רץ קודם באותה פונקציה. אומת: symlink קיים לתיקיית קוד אחרת מופנה מחדש בשקט, בלי חריגה.

### B269. `setup_template` — גם השומר "קיים ואינו symlink" מת מאותה סיבה
המקרה היחיד ש-`os.unlink` לא מסתדר איתו הוא תיקייה, ואז מתקבל `IsADirectoryError` חשוף במקום ההודעה המאובחנת.

### B270. `setup_template` — התקנה ראשונה נקייה מדווחת שתי שגיאות
הרצה ראשונה מדפיסה שני `logger.error("file … not found")` על היעדר צפוי לחלוטין של הקישורים שהיא עצמה עומדת ליצור.

### B271. `_toNetcdf` — ממיר דוסאז' בכיוון ההפוך מהמתועד
**קובץ:** `hera/simulations/LSM/template.py` · הדוקסטרינג אומר "dosage are converted to s/m³ instead of min/m³", אבל הקוד מכפיל ב-`(1*ureg.s/ureg.m**3).m_as(ureg.minute/ureg.m**3)` = 1/60 — כלומר פקטור s→min. אומת: שורת OUTD עם Dosage 6.0 חוזרת כ-0.1.

### B272. `_getSimulationsList` — מכפלה קרטזית במקום זיווג
**קובץ:** `hera/simulations/LSM/template.py`

`product(desc_df_list, params_df_list)` במקום לזווג. אומת: שתי סימולציות → 4 פריימים, כלומר ה-id/שם של הרצה אחת מוצמד לפרמטרים של אחרת. `getSimulationsTable` יורש את זה (4 שורות, כל id פעמיים), ו-`LSMToolkit.getSimulationsList` ו-`getTemplatesTable` משחזרים את אותה בנייה.

### B273. `run(topography=..., stations=None)` — `AttributeError` על `stations.columns`
התנאי הוא `if (topography is not None or stations is not None)` ואז נקרא `stations.columns` → קורס כשרק טופוגרפיה הועברה.

### B274. `run` — `os.chdir(saveDir)` בלי `try/finally`
כל כשל אחרי השורה הזו משאיר את **כל המפרש** בתוך תיקיית הסימולציה.

### B275. `run(stations=...)` — `KeyError: 'station'` בנתיב המתועד
מחשב `stations["station"] = x+y`, ואז בונה מחדש את הפריים מ-xarray ממודגם שמפיל את העמודה, ומצרף בחזרה `onlyStations` שנבנה **לפני** שהעמודה נוצרה → `KeyError` בשורה הבאה. עובד רק אם הקורא במקרה סיפק עמודת `station`, מה שלא מתועד בשום מקום.

### B276. אותו נתיב — צירוף many-to-many שמשכפל כל רשומה
`onlyStations` לא עובר דה-דופליקציה והצירוף על `(x,y)` הוא many-to-many, כך שכל רשומה ממודגמת נכתבת פעם אחת לכל תצפית מקורית. אומת: 3 דגימות של עשר דקות → 5 צעדים מאונטרפולטים שנכתבים ×3 = 15 רשומות; הפותר קורא סדרה מגומגמת.

### B277. `getSimulations` — `_document` בלי בדיקת None
**קובץ:** `hera/simulations/LSM/toolkit.py` · `getTemplateByName(unitsTemplateVersion)._document` → `AttributeError: 'NoneType' object has no attribute '_document'` בכל פרויקט שאין בו תבנית בשם `"v4-general"` (ברירת המחדל).

### B278. `getSimulations` — `except:` עירום שמפיל מסמכים בשקט
עוטף `SingleSimulation(doc)` ב-`except:` שמדפיס אזהרה ומשמיט את המסמך: הקוראים מקבלים רשימה קצרה בשקט, האבחון עובר ליד ה-logger, וגם `KeyboardInterrupt`/`SystemExit` נתפסים.

### B279. `getSimulationsList(wideFormat=True)` — מדווח "לא נמצא" על פרויקט עם שתיים
הכפילויות מ-B272 גורמות ל-`pivot` לזרוק `ValueError("Index contains duplicate entries")`, וה-`except ValueError` (שנועד לדווח על ריקנות) הופך את זה ל-`FileNotFoundError("No simulations.old found")`.

### B280. `prepareSlurmLSMExecution` — נקודת הכניסה ל-Slurm בלתי שמישה לחלוטין
**קובץ:** `hera/simulations/LSM/toolkit.py:366` מול `template.py:330`

```python
LSMTemplate.prepareParams(desc=None, paramsToPrepare=jsonConfig)   # החתימה: (template_desc, paramsToPrepare)
```

כל קריאה → `TypeError: prepareParams() got an unexpected keyword argument 'desc'`, **אחרי** שתיקיית הסקריפטים ו-`stations.parquet` כבר נכתבו. אומת בקוד.

### B281. `prepareSlurmLSMExecution` — `logger.error` בלי `return`, ×2
גם ל-`baseParameters` שאינו dict וגם ל-`jsonVariations` בלתי שמיש. אומת: `baseParameters="notADictionary"` → `ValueError: dictionary update sequence element #0 has length 1; 2 is required` מתוך `JSONVariations`. מקרה ה-dict גרוע יותר: ההודעה אומרת "path or dict" בעוד שבדיקת ה-`isinstance` דוחה dict.

### B282. `Continuous.calc` — הקרנל בצורה שגויה, לא יכול להחזיר
**קובץ:** `hera/simulations/gaussian/gasCloud.py`

בונה קרנל `(kernelsize, nx, ny, nz)`, אבל `rolling(...).reduce()` של xarray מעביר ל-reducer צורה `(nt, nx, ny, nz, kernelsize)` עם `axis=-1`; `_convolve` מאנדקס את הקרנל על ציר 0 לפי `data.shape[0]` (מספר צעדי הזמן) ומחתך 3 מתוך 4 מימדים. אומת: `Continuous(dt=1*ureg.min, kernelsize=3).calc(ones((5,2,2,2)))` → `ValueError: operands could not be broadcast together with shapes (5,2,2,2,3) (2,2,2,2)`.

### B283. `CirclePositionClippedDropletsCloud` — הרוחב הצידי נלקח מ-`position[0]`
**קובץ:** `hera/simulations/gaussian/DropletCloud.py`

המיקומים נבנים כ-`(position[0]+r·cos, position[0]+r·sin, position[2])`. אומת: `position=(100 m, 250 m, 10 m)` שם את הטבעת סביב (100, 100) — כ-150 מ' מהמרכז המבוקש. `LinePositionDropletsCloud` עושה את זה נכון.

### B284. אותה מחלקה — ענף `radius is None` מת
`tounit(radius, ureg.m)` בשורה שמעליו זורק `TypeError: Invalid magnitude for Quantity: None`.

### B285. אותה מחלקה — `linspace` עם שני הקצוות, ולכן 0° מקבל מנה כפולה
הזוויות הן `numpy.linspace(0, 2*pi, circlepositions)` כולל קצוות, כך ש-0 ו-2π מתלכדים: n מיקומים תופסים n−1 מקומות, זה של 0° מקבל מנה כפולה, והמרווח הוא 2π/(n−1). צריך `endpoint=False`.

### B286. `hera-LSM` — ה-CLI לא יכול לשגר שום תת-פקודה
**קובץ:** `hera/bin/hera-LSM:41`

```python
parsed = parser.parse_args()()
```

הסוגריים הכפולים **קוראים** ל-`Namespace` שהוחזר → `TypeError: 'Namespace' object is not callable`. אומת בקוד. לא מקובע בבדיקה כי הקובץ הוא סקריפט בלי סיומת שאינו ניתן לייבוא (אותה מגבלה כמו B182).

## אצווה 32 — meteorology (B296-B308) ו-side effect גלובלי (B309)

מקובע ב-`test_meteorology_turbulencestatistics_more.py`, `test_meteorology_highfreq_analysislayer.py`, `test_meteorology_lowfreq_analysis_more.py`, `test_meteorology_lowfreq_presentation_more.py`.

### B296. `Plots._scatterdict` — `size=` במקום `s=` ל-seaborn
**קובץ:** `hera/measurements/meteorology/lowfreqdata/presentationLayer.py:111` · `size` הוא שם של **משתנה סמנטי** ב-seaborn, לא שטח הסמן (`s`). הסמנים מצוירים בברירת המחדל 18, והמפתח `size` התקוע גם מסתיר `s` מפורש של הקורא. אומת: `get_sizes()` → `[18.0, …]`.

### B297. `dateLinePlot` — דורס את תווית ה-y שהוגדרה
מחיל את `set_ylabel` מ-`_plotfieldaxfuncdict` ואז דורס ללא תנאי בשם העמודה החשוף: `WS` הוא `'Wind Speed [m/s]'` בגרף פיזור/קונטור אבל `'WS'` בגרף קו. גם `set_ylabel` של הקורא נזרק כך.

### B298. `_labelsdict['levels']` — מוקפא מברירת המחדל ולא נבנה מחדש
כל `contour_values` שהקורא מספק גורם ל-`ax.clabel` לזרוק `ValueError: Specified levels … don't match available levels`. `withLabels` ברירת המחדל `True`, כך ש-`plotProbContourf(df,"WS",contour_values=…)` נכשל מיד.

### B299. `plotProbContourf` — ה-clamp של "רמות עולות" לא מוחל על תוויות
בענף `y_normalized` ה-clamp מוחל על רמות הקונטור אבל לא על רמות ה-clabel (שנבנות מחדש מ-`M_hist.max()` הלא-מוגבל). היסטוגרמה שהמקסימום שלה נופל מתחת ל-`under_value` נותנת רמות תוויות **יורדות** → `clabel` זורק.

### B300. `plotProbContourf_bySeason` — ה-`ax` המתועד בלתי שמיש
מאנדקס `ax.shape`/`ax[i,j]` (צריך מערך 2×2) ומזין את אותו אובייקט ל-`plt.sca` (צריך Axes בודד). מערך → `AttributeError: 'numpy.ndarray' object has no attribute 'figure'`; Axes בודד → `AttributeError` על `.shape`.

### B301. `singlePointTurbulenceStatistics`/`AveragingCalculator` — `inmemory` מתועד ולא נקרא
לא מגיע לא למזהה ולא למחשבן, כך ש-`inmemory=True` ו-`False` מייצרים מטא-דאטה זהה.

### B302. `AveragingCalculator` — `isMissingData` חסר מהמזהה
בניגוד ל"אח" שלו, הוא לא נכנס למזהה — כך שהרצה עם חורים בנתונים בלתי-נבדלת מהרצה שלמה, והן **מתנגשות** בשאילתת ה-Cache ש-`AbstractCalculator` בונה מאותו מטא-דאטה.

### B303. `StrucFun` — `raise` על מחרוזת
**קובץ:** `hera/measurements/meteorology/highfreqdata/analysis/turbulencestatistics.py:1260`

```python
raise("mode must be either MeanDir or 3dMeanDir")
```

הקורא מקבל `TypeError: exceptions must derive from BaseException`; ההודעה האבחנתית לא מגיעה אליו לעולם. אומת בקוד.

### B304. `StrucFun`/`ThirdStrucFun` — השומר בודק שם שלא נוצר
השומר הוא `if "u_mag" not in self.TemporaryData.columns` אבל הכתיבה/הרישום הם ל-`"u_mag" + title_additions`. עם `title_additions="X"` השומר בודק שם שלא קיים, וכל קריאה חוזרת מוסיפה עוד `["u_magX", {}]` ל-`_CalculatedParams` (2 רשומות, עמודה אחת), ו-`_compute` מקרין לפי הרשימה הזו.

### B305. משפחת פונקציות המבנה — dask בלבד, למרות שמתועד pandas או dask
`StrucFunDir` → `AttributeError: repartition`; שדה המיצוע של `StrucFun` → `AttributeError: compute`.

### B306. `SinglePointStatisticsSpark.fluctuations` — לא קובע `self.data`
הבסיס כן קובע. לכן `getData()` מחזיר `None` גם אחרי שהכול חושב — מה שהדוקסטרינג של המתודה עצמה מגדיר כ"לא בוצעו חישובים".

### B307. אותו override — `wind_dir_bar` נכתב ולא נרשם
נכנס ל-`_TemporaryData` אבל לא ל-`_CalculatedParams`, כך ש-`compute()` מפיל בשקט את כיוון הרוח הממוצע (`_param_names == ['u_bar','v_bar','w_bar','T_bar']`). המימוש בבסיס כן רושם אותו.

### B308. `InMemoryRawData.to_hdf` — נתיב בלי סיומת נכתב ל-CWD
שני הנתיבים נגזרים מ-`path_or_buf.rpartition('.')[0]`, שהוא ריק לנתיב בלי סיומת: כתיבה ל-`<dir>/rawdata` שמה את ה-`.hdf` ואת ה-`.json` ב**תיקיית העבודה** ולא ב-`<dir>`.

### B309. `presentation/basicplots.py` — `seaborn.set()` בזמן ייבוא משנה צבעים גלובלית
**קובץ:** `hera/presentation/basicplots.py:4`

```python
import seaborn as sns; sns.set()
```

הקריאה מתבצעת ב**זמן ייבוא המודול**, ו-`seaborn.set()` מחדש את קיצורי הצבע החד-אותיים של matplotlib ברמת התהליך. אומת ניסויית: אותה קריאה בדיוק ל-`ax.scatter(..., c="r")` ב-`experiment/presentation.py` מחזירה `(1.0, 0.0, 0.0, 1.0)` כשהמודול לא יובא, ו-`(0.769, 0.306, 0.322, 1.0)` אחרי שיובא. כלומר **ייבוא של מודול hera אחד משנה את מה ש"אדום" אומר בכל שאר האפליקציה** — כולל בקוד שלא יודע ש-seaborn קיים. זה גם מה שגרם לבדיקה שעברה בבידוד להיכשל בהרצה מלאה; הבדיקה תוקנה להשוות מול הצבע ש-matplotlib מפענח בפועל (`to_rgba("r")`) במקום מול RGBA ספרותי. הכשל הזה אינו רק אסתטי: הוא הופך פלט גרפי לתלוי בסדר הייבוא.
