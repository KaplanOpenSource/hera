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
**קובץ:** `hera/measurements/GIS/raster/hill2stl.py` · **מקובע ב:** `test_gis_utils.py::TestImportSideEffects`

בתחתית הקובץ, תחת ההערה `# Run the function`, יושב קוד דמו ברמת המודול:

```python
minx=-2; maxx=2; miny=-3; maxy=3
filename='test1.stl'
generate_solid_stl(function, x_range=(minx,maxx), y_range=(miny,maxy), resolution=100, filename=filename)
```

כל `import hera.measurements.GIS.raster.hill2stl` מדפיס ל-stdout וכותב **8.1MB** בשם `test1.stl` לתיקיית העבודה. כל כלי שסורק את עץ החבילה — בניית תיעוד, אינדוקס ב-IDE, איסוף טסטים — משאיר את הקובץ אחריו.

**נסיבה מקלה שאומתה:** `import hera.measurements.GIS` לבדו **אינו** מפעיל את זה; רק ייבוא ישיר של המודול.

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

### B53. `ProtectionPolicy.addActions` שבור עם נתיב קובץ
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:130` · **מקובע ב:** `test_riskassessment_protectionpolicy.py::TestActionDispatch::test_addactions_from_a_json_string_is_broken`

```python
if os.JSONpath.exists(jsonStrOrFile):
```

ל-`os` אין תכונה `JSONpath` (הכוונה הייתה `os.path`). כל קריאה ל-`addActions` עם מחרוזת נתיב קובץ נכשלת עם `AttributeError` — הענף הזה מעולם לא רץ.

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

### B58. `ProtectionPolicy` קורס בשימוש הפשוט ביותר — בלי חלון זמן
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:390,481` · **מקובע ב:** `test_riskassessment_protectionpolicy.py::TestActionIndoor::test_compute_with_no_time_window_at_all`

```python
abegin = data[self.policy.datetimename].to_series()[0] if abegin is None else abegin
```

אינדקס positional (`[0]`) על `Series` עם `DatetimeIndex` — לא נתמך יותר ב-pandas המותקן (זורק `KeyError` במקום fallback ל-position). קריאה ל-`.indoor(turnover=...)` **בלי** `begin`/`end`/`enter`/`stay` כלל — השימוש הבסיסי ביותר — קורסת תמיד.

### B59. חישוב ה-indoor הוא no-op מוחלט — התוצאה תמיד אפס
**קובץ:** `hera/riskassessment/protectionpolicy/ProtectionPolicy.py:374` · **מקובע ב:** `test_riskassessment_protectionpolicy.py::TestActionIndoor::test_compute_builds_a_low_pass_filtered_field_towards_the_outdoor_value`

```python
Cin[curstep].values = ((Cin[prevstep] + alphanum*dt*Cout[curstep])/(1+alphanum*dt)).values
```

`Cin[curstep]` (אינדוקס לפי dict) מחזיר **עותק**, לא view, בגרסת ה-xarray המותקנת. ההשמה ל-`.values` על העותק נזרקת — `Cin` נשאר אפסים מוחלטים בלי קשר ל-alpha, dt או לריכוז החיצוני. זהו החישוב המרכזי של כל מודל ה-indoor, ומעולם לא באמת רץ.

### B60. `getRiskAreaAlgorithm` — dispatcher מת לחלוטין
**קובץ:** `hera/riskassessment/analysis/riskAreas.py:19` · **מקובע ב:** `test_risk_riskareas.py::TestGetRiskAreaAlgorithm`

```python
estimatorCLS = pydoc.locate("pyriskassessment.datalayer.riskAreas.riskAreaAlgorithm_%s" % algorithmName.title())
```

החבילה `pyriskassessment` אינה קיימת בפרויקט או בתלויותיו — המודול הזה יושב ב-`hera.riskassessment.analysis.riskAreas`. `pydoc.locate` מחזיר `None` תמיד, ולכן ה-dispatcher **לעולם לא מצליח**, אפילו לא עבור `"Sweep"` — האלגוריתם היחיד שממומש בקובץ. הדרך היחידה להשתמש ב-`riskAreaAlgorithm_Sweep` היא import + בנייה ישירה.

### מה שנמצא תקין ב-riskAreas/riskToolkit
תכונות ה-setter/getter של `riskAreaAlgorithm_Sweep` (`dxdy`, `workerCount`, `parallel`) עובדות וממירות טיפוסים נכון; `outlayers` הוא read-only כמתועד. `RiskToolkit.getAgent` מבחין נכון בין שם (str), תיאור (dict) וקלט לא-תקין; `loadAgent` מדביק name/version לתיאור לפני העברה ל-`loadData`. `loadData` מטפל נכון בקלט dict, JSON string, ומקרה "לא קובץ/dict" (raises).

`riskAreaAlgorithm_Sweep.calculate()` ו-`RiskToolkit.analysis.getRiskAreas` לא נבדקו — הראשון דורש multiprocessing + geopandas מלא עם demog/isopleths, השני דורש סימולציית LSM חיה ו-agent דוז-רספונס מלא; שניהם מועמדים לטסט אינטגרציה, לא הרמטי.

---

## אצווה 13 — `openFoam` preprocessing, המשך (`OFObject`, `OFList`, `preprocessOFObjects/utils.py`)

### B61. `OFList` הוא קוד מת — אף קורא לא בונה אותו, ואי אפשר לכתוב איתו כלום
**קובץ:** `hera/simulations/openFoam/preprocessOFObjects/OFList.py:55` · **מקובע ב:** `test_openfoam_objects_base.py::TestOFListIsUnusable`

```python
fileStrContent = self.getHeader()
```

זו השורה הראשונה של `_writeNew` (וגם `_updateExisting` קוראת ל-`_writeNew`). אף מחלקה בהיררכיה של `OFList` (`OFList`, `OFObject`) מגדירה `getHeader` — היא מוגדרת רק על `OFField`, מחלקת-אח לא קשורה. שורה אחת אחרי זה, `self.columnNames` (בניגוד ל-`componentNames` שכן מוגדר) גם אף פעם לא נקבע. כל קריאה, גם בענף הסקלרי, קורסת ב-`AttributeError` לפני שההסתעפות סקלר/וקטור בכלל מגיעה להרצה. חיפוש בכל עץ `hera/simulations/openFoam` לא מוצא אף מקום שבונה `OFList(...)` — נראה שהמחלקה מעולם לא שולבה בזרימת העבודה בפועל.

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

### B65. `thresholdGeoDataFrame.project` — קוד שלא ניתן להרצה בפייתון הנוכחי
**קובץ:** `hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:77` · **מקובע ב:** `test_risk_thresholdgeodataframe.py::TestProjectIsUnusableOnThisPython`

```python
if isinstance(meteorological_angle,collections.Iterable):
```

`collections.Iterable` הוסר מ-Python 3.10 (הועבר ל-`collections.abc.Iterable` עוד ב-3.3, וה-alias הוסר לגמרי ב-3.10). כל קריאה ל-`project` — כולל עם `demographic=None` — קורסת ב-`AttributeError` **לפני** שהיא בכלל מגיעה לטפל בפרמטרים.

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

## אצווה 17 — `measurements/GIS` (חלקי: hill2stl, stlFactory, TilesToolkit, LandCoverToolkit)

### B74. `vectorToSTL` — סניף dask שקורא לשם שלא מיובא בקובץ
**קובץ:** `hera/measurements/GIS/utils.py:389` · **מקובע ב:** `test_gis_utils_stlfactory.py::TestVectorToSTLDispatch`

```python
elif isinstance(gpandas, pandas.DataFrame) or isinstance(gpandas, dask.dataframe.DataFrame):
```

`dask` אף פעם לא מיובא במודול הזה. כל עוד הקלט הוא `geopandas.GeoDataFrame` או `pandas.DataFrame` רגיל, ה-`or` הקצר-מעגל מונע מהביטוי השני להיבדק — אבל כל קלט שהוא לא אחד מהשניים (כולל `dask.dataframe.DataFrame` אמיתי, המתועד כנתמך במקביל ל-`rasterizePandas`) קורס ב-`NameError: name 'dask' is not defined` במקום דחייה נקייה.

### הערת קוד: `roughnesslength2sandgrainroughness` מוגדרת פעמיים באותה מחלקה
**קובץ:** `hera/measurements/GIS/raster/landcover.py:561,717` — שתי ההגדרות `@staticmethod` בתוך `LandCoverToolkit` מחשבות את אותה נוסחה (`z0*30`, לפי Desmond et al. 2017 eq. 5). ההגדרה השנייה (717) מחליפה את הראשונה בסמנטיקת גוף המחלקה הרגילה — ההגדרה הראשונה היא קוד מת, אבל ההתנהגות בפועל לא נפגעת (שתיהן מחשבות את אותו הדבר). לא נפתח כממצא B נפרד כי אין כשל התנהגותי.

### מה שנמצא תקין
`hill2stl.py` (4 פונקציות): `compute_normal` מחזיר וקטור יחידה מאונך לשני הצלעות; `write_triangle`/`generate_solid_stl` מייצרים STL ASCII תקין (facet/endfacet מאוזנים, בסיס שטוח מתחת למשטח). `stlFactory.heightColumnsNames` (getter/setter) ו-`rasterizeGeopandas` — אינטרפולציה לינארית נכונה מקווי-מפלס (LineString) לרשת רגולרית, מכבדת שם עמודת גובה מותאם. `TilesToolkit` — מתמטיקת Slippy Map סטנדרטית: `tileScaleAtLatLonZoom` נחצה בכל רמת zoom ומצטמצם לקוטבים, `deg2tile`/`tile2deg` הם היפוכים מקורבים, `doctype` ו-`setDefaultTileServer` (נשמר ב-config הפרויקט) עובדים נכון.

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

### B84-B86. `absractEulerianSolver_toolkitExtension` — שלושה NameError ברצף
**קובץ:** `hera/simulations/openFoam/eulerian/abstractEulerianSolver.py` · **מקובע ב:** `test_openfoam_eulerian_solver.py`

- **B84**: `flowType` — ענף ה-`else` מפנה ל-`SIMULATIONTYPE_COMPRESSIBLE`, שם שלא קיים בשום מקום בקובץ (הייבוא בראש הקובץ הוא `FLOWTYPE_COMPRESSIBLE`). כל קריאה עם `incompressible=False` קורסת.
- **B85**: `blockMesh_setBoundFromFile` — הפרמטר נקרא `eulerianWorkFlow`, אבל הגוף בודק `isinstance(eulerianWF, workflow_Eulerian)` — `eulerianWF` לא קיים בכלל כפרמטר. קורס תמיד לפני שנוגע בקלט.
- **B86**: `blockMesh_setDomainHeight` — העתק-הדבק של `blockMesh_setBoundFromFile` בלי התאמה: אותו `eulerianWF` לא מוגדר, וגם מפנה ל-`fileName`/`dx`/`dy` — אף אחד מהם לא פרמטר של המתודה הזו (`eulerianWorkFlow, Z, dz`). `Z`, הקלט האמיתי היחיד שלה, אף פעם לא בשימוש.

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
