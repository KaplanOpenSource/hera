# סיכום מהיר לפגישה — מה עשינו ב-21 הענפים

מסמך רפרנס מהיר, לא רשמי. המקורות המלאים: `CONSOLIDATED-ISSUE-DRAFT.md` (כל הבאגים), `2026-08-25-how-the-new-test-layer-works.md` (הסבר טכני מלא, פורסם גם ב-[Issue #1061](https://github.com/KaplanOpenSource/hera/issues/1061)).

---

## תקציר ב-30 שניות

בנינו שכבת טסטים חדשה שרצה תוך שניות בלי תלות בתשתית (Mongo/S3), כיסינו איתה כמעט חצי מהפונקציות בקוד שמעולם לא נבדקו, ומצאנו **78 באגים אמיתיים** — קוד שקורס בפועל, לא שגיאות תיאורטיות. שום דבר לא תוקן בייצור (חוץ מחוסמי-הרצה טכניים). הכל בענפים נפרדים, מוכן למיזוג, אבל **עדיין לא ממוזג**.

---

## מה זה mongomock, ולמה זה שונה מהטסטים הקודמים

**הטסטים הקודמים (713 טסטי אינטגרציה):** רצים מול **MongoDB אמיתי** שצריך להרים (`make mongo-up`), ולפעמים גם נתוני S3 אמיתיים (כמה מאות MB). הם בודקים "מלמעלה" — קוראים ל-Project/Toolkit ומשווים תוצאה סופית ל-baseline. **הבעיה:** כדי להריץ אותם צריך תשתית מלאה, ולכן הם רצים רק כשמישהו טורח — לא על כל שינוי קטן. וגם: הם לא מגיעים לכל פונקציה בנפרד, רק לזרימות שלמות.

**mongomock:** ספרייה שמדמה את הפרוטוקול של MongoDB אבל שומרת הכל **בזיכרון** (RAM), לא על דיסק ולא בשרת. **קריטי להבין:** זה **לא** מוק של הקוד של hera. שכבת ה-datalayer, ה-Project, ה-Toolkit — **כל הקוד האמיתי של hera רץ במלואו**, בדיוק כמו בייצור. רק הרגע שבו הוא היה צריך "לדבר" עם שרת Mongo אמיתי — שם הוא מדבר במקום עם עותק בזיכרון. זו לא סימולציה של ההתנהגות, זו ההתנהגות האמיתית, רק בלי השרת.

**למה זה "נכון" (תקף כשיטת בדיקה):** כי אנחנו לא בודקים "האם Mongo עובד" (זה כבר ידוע/מובטח), אלא "האם **הקוד שלנו** מתנהג נכון מול Mongo". mongomock מספיק דומה כדי לבדוק את זה, ומהיר מספיק כדי לרוץ על כל commit. זה לא תחליף לטסטי האינטגרציה — זה **שכבה נוספת**, כי יש דברים ש-mongomock לא תופס (למשל ביצועים אמיתיים, אינדקסים, replica sets).

**התוצאה המעשית:** טסט שלקח 30+ שניות (המתנה לחיבור אמיתי) רץ עכשיו במילישניות. זה מה שמאפשר להריץ 1,370 טסטים ב-6-10 שניות.

---

## הקשר לגרסת פייתון/pandas — למה זה משמעותי

חלק ניכר מ-78 הבאגים **לא** קשור ללוגיקה עסקית שגויה — הם קוד שנכתב לפני שנים, מול גרסאות ישנות של Python/pandas, ומעולם לא רץ מחדש מול הגרסה המותקנת היום (**Python 3.12, pandas 3.0**). דוגמאות קונקרטיות שמצאנו:

- **`collections.Iterable`** הוסר מ-Python **3.10** (הועבר ל-`collections.abc.Iterable` עוד ב-3.3, אבל ה-alias הישן נשמר עד 3.9). קוד שכתוב `isinstance(x, collections.Iterable)` **קורס תמיד** על 3.10+.
- **`pandas.DataFrame.append`** הוסרה לגמרי ב-pandas **2.0**. קוד ישן שקורא לה קורס.
- אינדוקס positional על `Series` עם `DatetimeIndex` (כמו `series[0]`) — פעם fallback-ל-מיקום, בגרסאות pandas חדשות זורק שגיאה.
- `pandas.DataFrame` קיבלה `.attrs` (מטא-דאטה) רק מגרסה 1.0 — קוד ישן שהשתמש ב-`hasattr(x,"attrs")` כדי להבחין בין pandas ל-xarray פשוט הפסיק לעבוד, כי גם ל-pandas יש את זה עכשיו.

**המשמעות:** הקוד **מעולם לא היה תואם** לגרסאות המודרניות האלה, אבל אף אחד לא ידע כי **אף טסט לא קרא לפונקציות האלה**. זו לא "תקלה חדשה שנוצרה" — זו תקלה ישנה שהתגלתה רק עכשיו כי לראשונה מישהו (הטסטים החדשים) בפועל קרא לקוד. זה בדיוק הטיעון המרכזי לכל הפרויקט: **coverage נמוך = חורים עיוורים שמצטברים ככל שהתלויות מתקדמות בלי שהקוד מתעדכן איתן.**

---

## מפת 21 הענפים — מה כל אחד עשה

| ענף | תוכן | ממצאים עיקריים |
|---|---|---|
| `tests/phase0-infra` | תשתית הבסיס: seam ה-mongomock, שכבת ה-stubs, conftest, CI, Makefile | — |
| `tests/batch1-utils` | `hera/utils` (angle, query, filter, lazy, statistics, slurm, zip, latex, units, contour, json, package) | B9 (LaTeX math-mode) |
| `tests/batch2-datalayer` | `datalayer` + `toolkit.py` הליבה | B18, B19, B21, B22, B23 |
| `tests/batch3-gaussian` | `Sigma.py`, `Meteorology.py`, `MeshUtils.py` (אומת מול Briggs 1973) | B24, B27, B29, B30 |
| `tests/batch4-physics-utils` | evaporation, deposition, hydrodynamics, windProfile, simulations/utils | B33, B34 |
| `tests/batch5-riskassessment` | Calculator, InjuryLevel (הבסיס) | B36, B37, B38, B39 |
| `tests/batch6-meteorology` | lowFreq analysis, module imports | B40, B41, תיקון-עצמי ל-P1 |
| `tests/batch7-gis` | GIS utils, landcover roughness | B43 (טבלאות חספוס סותרות) |
| `tests/batch8-experiment` | experiment parsers | B44 |
| `tests/batch9-simulations` | openFoam/LSM/mlDL/WRF imports, LSM template params | — |
| `fix/test-blockers` | **היחיד עם שינוי בקוד ייצור** — תיקון import מעגלי, imports ברמה לא נכונה, side-effects בזמן import | — |
| `tests/batch10-pure-physics` | deposition/evaporation models (שוחררו ע"י תיקון ה-blockers) | B49, B50, B51 |
| `tests/batch11-openfoam` | OFObjectHome (dimension vector, field registry) | B52 |
| `tests/batch12-risk-pure` | ProtectionPolicy, riskAreas, riskToolkit | B53-B60 (כולל B59: מודל indoor תמיד אפס, B60: dispatcher מת) |
| `tests/batch13-openfoam-deep` | OFObject, OFList, utils.py (parser) | B61 |
| `tests/batch14-risk-agents` | Agents, Calculator, Injury, InjuryLevel | B62-B66 |
| `tests/batch15-openfoam-field` | OFField (5 מתודות אחרונות) | B67, B68 |
| `tests/batch16-meteorology` | highfreq calculators (abstractcalculator, meandatacalculator, turbulencestatistics, parsers) | B69-B73 |
| `tests/batch17-gis` | hill2stl, stlFactory, TilesToolkit, landcover | B74 |
| `tests/batch18-experiment` | experimentHome, experimentAnalysis, dataEngine | B75-B77 |
| `tests/batch19-gaussian-gascloud` | gasCloud.py (מודל הפיזור הגאוסי המלא) | **B78** — כל מחלקת השחרור-הרציף לא עובדת |

**כולם מוערמים בשרשרת לינארית אחת** (כל ענף מכיל את כל מה שלפניו) — `tests/batch19-gaussian-gascloud` הוא הקצה עם הכל.

---

## התוכנית הכללית — לאן זה הולך

1. **מיזוג** — 21 ענפים, אפס PRs פתוחים כרגע. אפשר PR אחד גדול (הקצה מכיל הכל) או שרשרת PRs לפי הסדר למעלה.
2. **החלטה על הבאגים** — טיוטת Issue מרוכזת עם כל 78 הממצאים קיימת (`CONSOLIDATED-ISSUE-DRAFT.md`) אבל לא פורסמה. הבאגים החמורים ביותר לתיקון: **B78** (שחרור רציף בפיזור גאוסי) ו-**B59** (מודל indoor dilution).
3. **מה שנשאר לא מכוסה** (764 מתוך 1,429 פונקציות) — מרוכז באזורים שדורשים תשתית אמיתית שלא זמינה כרגע: פותר OpenFOAM, VTK/paraview, torch מלא, hermes/luigi, GIS shapefiles. לא המשך-עבודה טבעי בלי השקעה נפרדת.
4. **המנגנון שנשאר לצמיתות:** מרגע שממזגים, כל PR עתידי ל-master **חייב** לעבור את שני ה-jobs ב-CI (unit הרמטי → integration עם Mongo אמיתי) ולא להוריד את הכיסוי המשולב מתחת ל-37% — כך שהחורים העיוורים האלה לא יחזרו.
