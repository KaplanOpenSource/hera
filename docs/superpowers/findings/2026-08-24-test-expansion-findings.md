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
