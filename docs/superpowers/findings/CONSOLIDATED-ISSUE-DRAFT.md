# טיוטת ה-Issue המרוכז

> **סטטוס: טיוטה בלבד. לא נפתח ולא פורסם.** ממתין לאישור.
> מקור מלא: [`2026-08-24-test-expansion-findings.md`](2026-08-24-test-expansion-findings.md)

---

**כותרת מוצעת:** `78 defects surfaced by the unit-test expansion (Phase 0 + batches 1-19)`

---

## תקציר

הרחבת סוויטת הטסטים (Phase 0 ותשע-עשרה אצוות) הוסיפה **1,370 טסטים הרמטיים** שרצים בכ-6-10 שניות בלי MongoDB, בלי נתוני S3 ובלי רשת, והעלתה את רצפת ה-coverage מ-20% ל-**37%** (כיסוי שורות+branch משולב עלה מ-22% ל-38%; כיסוי **פונקציות** — מדד מחמיר יותר, פונקציה נחשבת מכוסה רק אם שורה בתוך ה-body שלה רצה בפועל — עלה מ-27% ל-**46.5%**, 665 מתוך 1,429 פונקציות ציבוריות). במהלכה נמצאו **78 פגמים** (B1–B78, בתוספת 8 ממצאי תשתית/תהליך מסומנים P1–P8), כולם אומתו בהרצה ולא בקריאת קוד, וכל אחד מקובע בטסט `xfail(strict=True)` — כך שתיקון יגרום לטסט להיכשל ולדרוש שיוסר המסמן.

**לא בוצע שום שינוי בקוד ייצור**, מלבד ענף `fix/test-blockers` אחד שתיקן חוסמי-הרצה טכניים (import מעגלי, imports ברמה לא נכונה, side effects בזמן import) — לא לוגיקה עסקית.

הדפוס החוזר, וההצדקה העיקרית למאמץ: **הפגמים החמורים ביותר אינם חישובים שגויים אלא קוד שמעולם לא רץ** — מתודות שקוראות לשם שלא קיים, מחלקות שיורשות מהאב הלא נכון, בדיקות טיפוס שלעולם לא מתקיימות. מדד coverage לבדו לא היה מוצא אף אחד מהם; הדרך היחידה לגלות אותם היא לקרוא לפונקציה בפועל.

---

## 1. קוד שאינו ניתן לייבוא או שאינו רץ

| # | מה | היכן |
|---|---|---|
| **B31** | `evaporation` ו-`deposition` — `from ..utils import` מפנה לחבילה ריקה במקום ל-`hera.utils`. **שתי חבילות פיזיקה מתות.** התיקון: נקודה אחת בשלושה קבצים | `evaporation/models.py:3`, `monaghan.py:5`, `deposition/models.py:2` |
| **B40** | `radiosonde.py` יורש מ-`datalayer.ProjectMultiDBPublic` — מחלקה שהוסרה | `measurements/meteorology/radiosonde.py` |
| **B41** | נתיב מוחלט לספריית בית של מפתח, נקרא בזמן ייבוא | `highfreqdata/__main__.py` |
| **B16** | `rag/*` ו-`SALibUtils` — `httpx`, `chromadb`, `sentence_transformers`, `llama_index`, `SALib` לא מוצהרים ב-`requirements.txt` (473 statements) | `utils/rag/*`, `utils/SALibUtils.py` |
| **B22** | `setDefaultRepository` תמיד זורק `RuntimeError` — **הפיצ'ר לא עבד מעולם** | `toolkit.py:1164` |
| **B33** | `skin_friction` — `numpy.log()` בלי ארגומנט, בשתי המחלקות | `nearWallFlow.py:174,298` |
| **B9** | ענף ה-math-mode ב-LaTeX משווה תו בודד למחרוזת בת 11 תווים; הדוגמה מהדוקסטרינג עצמו יוצאת שבורה | `utils/latex.py` |
| **B23** | חיפוש פורמט שמנסה שלושה שמות שאינם קיימים | `toolkit.py:1164` |
| **B44** | `Parser_TOA5.parse` — גוף ריק, מחזיר `None` בשקט | `experiment/parsers.py:359` |
| **B50** | `evaporationModels` לא ניתן לבנייה — `RiskToolkit.getAgent(agent)` נקרא על המחלקה, אבל זו מתודת מופע; כל בנייה נכשלת ב-`TypeError` | `simulations/evaporation/models.py:71` |
| **B51** | `flux` לא יכול לרוץ — טמפרטורה מומרת ל-מספר חשוף (`dimensionless`) ואז מועברת ל-`vaporPressure` שדורש יחידות | `simulations/evaporation/models.py:122` |
| **B52** | `pandasToFoamFormat` לא ניתן לקריאה בשום צורה — `@staticmethod` על חתימה שמתחילה ב-`self`, וגם יושבת על מחלקה לא נכונה | `preprocessOFObjects/OFObjectHome.py:238` |
| **B53** | `ProtectionPolicy.addActions` עם נתיב קובץ — `os.JSONpath.exists` (אין ל-`os` תכונה כזו) | `protectionpolicy/ProtectionPolicy.py:130` |
| **B60** | `getRiskAreaAlgorithm` — dispatcher מת לחלוטין: מפנה לחבילה `pyriskassessment` שלא קיימת בפרויקט בכלל, אפילו לא לאלגוריתם היחיד שהוא בעצמו מיישם | `riskassessment/analysis/riskAreas.py:19` |
| **B61** | `OFList` לא יכול לכתוב שום דבר בשום ענף — קורא ל-`self.getHeader()` שלא מוגדר בשום מקום בהיררכיה שלו; שום קוד בפרויקט לא בונה `OFList` | `preprocessOFObjects/OFList.py:55` |
| **B65** | `thresholdGeoDataFrame.project` — קוד שלא ניתן להרצה בפייתון הנוכחי: `collections.Iterable` הוסר מ-Python 3.10 | `agents/effects/thresholdGeoDataFrame.py:77` |
| **B69** | שמירה כ-HDF5 שבורה — `to_HDF` במקום `to_hdf` | `highfreqdata/analysis/abstractcalculator.py:263` |
| **B70** | בניית `MeanDataCalculator` ישירות מ-`pandas.DataFrame`/`dask.DataFrame` — קוד מתועד שלא רץ: הבדיקה משווה טיפוס למודול עצמו, שלעולם לא מתקיים | `highfreqdata/analysis/meandatacalculator.py:79` |
| **B71** | `anisotropyEigs` — `_eig` חסרה `self` בחתימה אבל נגישה כ-bound method; `TypeError` בכל קריאה | `highfreqdata/analysis/meandatacalculator.py:548` |
| **B72** | `InMemoryRawData.append` קוראת ל-`pandas.DataFrame.append` שהוסרה ב-pandas 2.0 | `highfreqdata/analysis/turbulencestatistics.py:1478` |
| **B74** | `vectorToSTL` — סניף dask שמפנה לשם `dask` שלא מיובא בקובץ בכלל | `measurements/GIS/utils.py:389` |
| **B75** | `experimentAnalysis.getTurbulenceStatistics` — קוראת ל-`kaijoData`/`kaijoHeight`, שמות שלא קיימים בשום מקום בפונקציה | `experiment/analysis.py:74` |
| **B78** | `continuousReleaseGasCloud` — **כל ארבע השיטות הציבוריות שלו קורסות.** יורש מ-`abstractGasCloud`, לא מ-`instantaneousReleaseGasCloud`, אבל קורא למתודות שמוגדרות רק על האח — ענף השחרור-הרציף של מודל הפיזור הגאוסי כולו לא שמיש | `simulations/gaussian/gasCloud.py:606,619,625,632` |
| **P5** | שישה קבצים תחת `hera/` אינם Python תקין | `gaussian/source.py` ועוד 5 |

## 2. שגיאות מדעיות — משפיעות על תוצאות סימולציה

| # | מה |
|---|---|
| **B43** | **שתי טבלאות חספוס סותרות לאותם קודי IGBP באותו קובץ.** אחת מצטטת Floors et al. 2021; השנייה מסומנת `# Example values` והיא רמפה אריתמטית. פער של **×100 למים ו-×1300 לשלג וקרח**, וסדר הפוך. `z0` מזין את פרופיל הרוח ומשם את כל הפיזור |
| **B27** | `MeshUtils.sigmaYName` מוגדר לעמודת ה-X — הענן נכפה איזוטרופי. קלט `sigmaY=50` מייצר σ_y של **5.000** |
| **B32** | `calculateR` מערבב נורמליזציות (מונה ב-N, מכנה ב-N−1). כל מתאם קטן ב-`(N-1)/N`; מודל מושלם מקבל 0.667 ב-N=3 ולעולם לא 1.0 |
| **B34** | `ReynoldsUm` מחזיר `80000/מטר` — מספר ריינולדס עם ממדים |
| **B24** | `getSigma` מחזיר סטיית תקן שלילית ל-x שלילי, ו-NaN ל-x שלילי מאוד |
| **B37** | מסלול ה-pandas ב-`CalculatorHaber` מחזיר **NaN בלבד** (אי-התאמת אינדקסים) |
| **B13** | `returnStandardize` מתועד כ-MKS ומחזיר יחידות מקוריות בשקט |
| **B49** | קצב השיקוע (deposition rate) אינו תלוי בצפיפות החלקיק בכלל, למרות ששיקוע גרביטציוני תלוי בצפיפות ישירות | `simulations/deposition/models.py` |
| **B59** | **מודל הישארות-במבנה (indoor dilution) מחשב אפס תמיד**, בלי קשר לפרמטרים — ההשמה לתוך `Cin[...]` נופלת על עותק של xarray, לא view | `protectionpolicy/ProtectionPolicy.py:374` |
| **B64** | מסלול ה-pandas ב-`CalculatorHaber` (המחשבון של חוק Haber) מכפיל DataFrame ב-Series בלי `axis=0` — התוצאה NaN מוחלטת עם עמודות מדומות | `agents/effects/Calculator.py:112` |
| **B66** | `getVolatility` מחזיר יחידות שגויות — g²/(mol·cm³) במקום ריכוז מסה g/cm³, כי `mol⁻¹` של המשקל המולקולרי אף פעם לא מבוטל | `riskassessment/Agents.py:264` |

## 3. שיבוש נתונים חוצה-הקשרים

| # | מה |
|---|---|
| **B19** | עטיפות הוספת המסמכים עושות `desc.setdefault` על מילון הקורא. מילון מטא-דאטה משותף לשני toolkits → **המסמך של השני נשמר בשם הראשון**, נעלם משאילתותיו ומופיע באחר |
| **B45** | `prepareParams` כותבת לתוך `template_desc['params']`. תבנית LSM לשימוש חוזר צוברת פרמטרים מריצות קודמות |
| **B12** | `tounit` בונה `Quantity` ב-registry הדיפולטיבי של pint ולא ב-`ureg` של hera. `tounit(1,"dunam")` זורק — היחידה היחידה שה-registry קיים בשבילה |
| **B39** | הקורבן הקונקרטי של B12: `InjuryLevelThreshold.getPercent` **זורק על כל עומס רעילות מספרי** |
| **P1** | `dictToMongoQuery` משטח רשימות, ולכן `$all` מתדרדר להתאמת **קידומת לפי מיקום**. ה-cache של הטורבולנציה פוגע רק כשהבקשה היא קידומת באותו סדר |
| **B56** | `RiskToolkit.loadData` — `TOOLKIT_SAVEMODE_FILEANDDB_REPLACE` עם version שונה מהקיים לא מוצא את הרשומה הישנה (מחפש לפי ה-version **החדש**) ו**יוצר כפילות** במקום להחליף |
| **B57** | גם כשההחלפה בפועל קורית (version תואם), רק `desc['version']` מתעדכן — שאר המטא-דאטה נשאר ישן בעוד ה-`resource` (JSON גולמי) מוחלף במלואו; metadata ו-payload מתפצלים |

## 4. תשתית ובנייה

| # | מה |
|---|---|
| **B26** | **`test_no_invalid_escapes.py` אינו מסוגל להיכשל על escape לא-חוקי** — הוא מסלים `SyntaxWarning` לשגיאה, פייתון מציג אותה כ-`SyntaxError`, וה-handler של קבצי legacy בולע אותה ל-`skip`. ~990 מקרים שלא יכולים לתפוס את מה שהם קיימים בשבילו |
| **P6** | שני כשלים ב-`TestAutoCache` **קיימים על `origin/master` נקי** (אומת בשתי גרסאות pytest). כל עוד הם קיימים, שלב אכיפת ה-coverage ב-CI לא מגיע לרוץ |
| **P7** | ~~`pytest.skip` על נתונים חסרים → CI ירוק~~ **נפתר בתשתית** ב-`--require-data` |
| **B42** | ~~ייבוא `hill2stl` מריץ דמו וכותב **8.1MB** לתיקייה הנוכחית~~ **נפתר**: הקובץ נמחק כקוד מת במאמץ ניקוי נפרד ב-master |
| **B20** | מפת הפורמטים ממופתחת על נתיבי מודול פנימיים של pandas; שדרוג מוריד DataFrame ל-pickle במקום parquet (ה-CI מוגן ע"י ההצמדה) |
| **P8** | סחיפת גרסאות בין `requirements.txt` לסביבה |

## 5. חוזי API והתנהגות מתועדת שאינה מתקיימת

**B18** כפילויות ב-`getDataSourceList` · **B21** `_get_data_toolkit` מתעלם מ-`projectName` · **B29** `GaussianIntegrationToMesh` מצמצם את קלט ההורה (הפרת Liskov) · **B30** שתי ברירות מחדל לפרופיל הרוח · **B36** ברירת מחדל של יחידות ל-pandas לא ניתנת להשגה · **B38** `InjuryLevelThreshold` נכשל אטומית על סף מספרי · **B14** אנוטציה `-> (str, dict)` על פונקציה שמחזירה tuple של 3 · **B15** `pandas.Dataframe` באות קטנה · **B5** זריקת מחרוזת · **B6** `UnboundLocalError` · **B28** תווית מוזנת ל-`.iloc` · **B7** `-mem=` במקום `--mem=` · **B35** פרמטר מתועד שאינו נקרא · **B25** תחילית `r` על המרכאות הסוגרות · **B8/B3/B10/B11** ועוד.

**מאצוות 12–19:** **B54** בניית `ActionIndoor` עם `alpha=` — הדוגמה הראשונה בדוקסטרינג של המחלקה עצמה — שוברת (`.ureg()` אינה מתודה של `Quantity`) · **B55** `TOOLKIT_SAVEMODE_NOSAVE` עדיין שומר בפועל בטעינה ראשונה של agent חדש · **B58** `ProtectionPolicy.compute` קורס בשימוש הבסיסי ביותר — בלי חלון זמן כלל (positional indexing על Series עם DatetimeIndex) · **B62** `Injury.calculateToxicLoads` מעביר `field` למקום positional שמתאים ל-TenBerge/MaxConcentration אבל לא ל-Haber · **B63** ל-`Injury.calculateToxicLoads`/`calculatePointWiseFractionInjured` אין פרמטר `inUnits`, ומאז ש-pandas.DataFrame קיבל `.attrs` (מגרסה 1.0) הבדיקה שנועדה להבחין מ-xarray כבר לא עובדת — כל קלט pandas קורס · **B67** `setFieldFromDataFrame` דורש עמודת `boundary` גם כשאין boundaries בכלל — round-trip עם `getDataFrame()` נכשל בוודאות · **B68** הענף המקבילי של `setFieldFromDataFrame` לא יכול לאתחל שדה חדש מאפס (קורא boundary קיים לפני כתיבה) · **B73** `zoL_Sonic` — מנגנון מניעת כפילויות בודק רשימה (`_AllCalculatedParams`) שאף פעם לא מתמלאת · **B76** `DBconfiguration` — property בלי `return`, בשתי מחלקות ה-data engine, מחזירה `None` תמיד · **B77** `dbConnect` — `except OperationFailure` על שם שלא מיובא, מסווה כל כשל חיבור אמיתי מאחורי `NameError` שנתפס בשקט.

---

## תיקונים מהירים ובעלי תשואה גבוהה

1. **B31** — נקודה אחת בשלושה קבצים מחזירה שתי חבילות פיזיקה לחיים
2. **B12** — `ureg.Quantity` במקום `Quantity` מתקן גם את B39
3. **B27** — תו אחד: `sigmaYName = "sigmaYCorrected"`
4. **B32** — `std(ddof=0)`
5. **B33** — הארגומנט החסר ל-`numpy.log()`
6. **B43** — למחוק את טבלת ה-`# Example values` ולהפנות ל-`_handleType1`
7. **B59** — לתקן את ה-indexing ל-`.loc[]`/`isel` שמחזיר view, לא dict-getitem שמחזיר עותק — מתקן את מודל ה-indoor dilution כולו
8. **B71** — להוסיף `self` לחתימת `_eig`
9. **B78** — להחליף את מחלקת האב של `continuousReleaseGasCloud` ל-`instantaneousReleaseGasCloud` (או להעביר את המתודות החולקות לבסיס המשותף)
10. **B69/B72** — `to_hdf`/הסרת הקריאה ל-`DataFrame.append` (או מעבר ל-`pandas.concat`)

---

## שני תיקונים לעצמי

לשקיפות, שתי טענות שרשמתי במהלך העבודה והתבררו כשגויות:

1. **P1** — טענתי ש"חיפוש ה-cache מעולם לא מוצא התאמה". **חזק מדי**: `$all` מתדרדר להתאמת קידומת לפי מיקום, כלומר עובד חלקית ותלוי-סדר.
2. **torch** — רשמתי שהוא לא מסוטב כי `modelContainer.py` יורש מ-`torch.nn.Module`. **לא נכון** — הוא לא יורש. הסיבה האמיתית: MagicMock חסר `__path__`, ולכן `import torch.utils` נכשל.

שתיהן תוקנו בקוד, בטסטים ובמסמך הממצאים.

---

## מה עדיין לא נבדק, ולמה (לא נכשל — לא נוגע)

כדי שהכיסוי לא ייקרא כ"הושלם" בטעות: כ-764 פונקציות ציבוריות (מ-1,429) נותרו לא מכוסות, מרוכזות ברובן באזורים שדורשים תשתית שלא הייתה זמינה בסביבה הרמטית:

- **`simulations/openFoam`** (222 חסרות) — פותר OpenFOAM אמיתי, VTK/paraview לרינדור, hermes ל-workflow.
- **`machineLearningDeepLearning`** (57) — `torch` מסוטב כ-MagicMock בודד בלי `__path__` לתת-מודולים; דורש namespace-package stub ייעודי.
- **`simulations/LSM`** ו-**`hermesWorkflowToolkit.py`** (51+28) — דורשים hermes/luigi אמיתיים לאורכסטרציית workflow.
- **`utils/rag`** (18) — תלויות RAG חסרות (ראו B16).
- קבצים ספציפיים שנדחו במפורש עם נימוק בכל אצווה: `CampbellBinary.py` (פרסר בינארי, שני עותקים), `experimentSetupWithData` (דורש zip אמיתי מ-argos), `buildings/demography/topography` ב-GIS (דורשים shapefiles/DEM אמיתיים), `FallingNonEvaporatingDroplets.py`/`DropletCloud.py` (נדחו מאצווה 3).

