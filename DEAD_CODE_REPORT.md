# דוח קוד מת וקוד מיותר — hera

דוח מרוכז מסריקה מקיפה של כל חבילת `hera/` (datalayer, measurements/GIS, measurements/meteorology, measurements/experiment, simulations, riskassessment, utils, toolkit.py, bin/CLI). כל ממצא אומת ב-grep על כל הריפו (כולל בדיקות, תיעוד, JSON) לפני שסומן כ"מת" — כדי לשלול שימוש דינמי (pydoc.locate, getattr, factories, __all__).

**המלצה כללית:** לפני מחיקה בפועל — להריץ טסטים קיימים, ולוודא עם בעלי הקוד (בעיקר ב-simulations/riskassessment שיש הרבה API ציבורי) שאין שימוש חיצוני (notebooks/סקריפטים מחוץ לריפו).

---

## 🔴 עדיפות דחופה — קוד שבור (לא רק מת, אלא יקרוס אם ירוץ)

אלו הממצאים החמורים ביותר: קוד שנראה "בשימוש" (מגיע אליו קריאה אמיתית) אבל מכיל שגיאה שתגרום לקריסה. חלקם חוסמים למעשה טולקיטים שלמים.

1. **`hera/simulations/gaussian/gasCloud.py`** — כל טולקיט הגאוסיאן שבור בפועל:
   - שורה 393-397: `instantaneousReleaseGasCloud` משתמש בברירת מחדל `dz=1*m, dt=1*min` — `m`/`min` הם שמות לא מוגדרים (הכוונה ל-`ureg.m`/`ureg.min`). ברירות מחדל של פרמטרים מוערכות **בזמן import**, כך שה-`NameError` קורה כבר בטעינת המודול. אומת בפועל (`ast.parse`+הרצה עם stubs).
   - כתוצאה מכך `hera/simulations/gaussian/toolkit.py` (עושה `from hera.simulations.gaussian.gasCloud import abstractGasCloud`) **לא ניתן לטעינה** — כלומר `gaussianToolkit` הרשום ב-`toolkitHome` שבור כרגע. סביר שזה מוסתר כי `pydoc.locate` בולע שגיאות import בשקט.
   - שורות 325-390: הגדרה כפולה וישנה יותר של `instantaneousReleaseGasCloud` — נדרסת מיידית ע"י ההגדרה השנייה (מתה לחלוטין, unreachable).
   - **תיקון:** להחליף `m`→`ureg.m`, `min`→`ureg.min` (שורה 397), למחוק את ההגדרה הראשונה הכפולה (325-390).

2. **`hera/simulations/gaussian/source.py:12`** — שגיאת תחביר ממש (`self.Q = #need to see...` בלי ערך) — `SyntaxError`, המודול לא ניתן לטעינה בכלל. אף אחד לא מייבא אותו כרגע, אז זה לא שובר משהו חי, אבל זה קוד לא תקין ב-100%.

3. **`hera/simulations/utils/interpolation/` (כל התיקייה)** — `interpolations.py:137` מכיל `TabError` (טאבים ורווחים מעורבבים) — לא ניתן לטעינה. זהו עותק כפול ושבור של `hera/simulations/utils/interpolations.py` (שכן בשימוש אמיתי ע"י `windProfile/toolkit.py`). **מחיקה מלאה של התיקייה הכפולה מומלצת** (גם דוח קודם `CODE_REVIEW_REPORT.md:188` ציין זאת).

4. **`hera/riskassessment/agents/effects/thresholdGeoDataFrame.py:78,90`** — `isinstance(x, collections.Iterable)` — `collections.Iterable` הוסר החל מ-Python 3.10 (עבר ל-`collections.abc.Iterable`). כל קריאה ל-`project()` עם רשימת זוויות תקרוס. **תיקון קל**: `collections.abc.Iterable`.

5. **`hera/riskassessment/protectionpolicy/ProtectionPolicy.py:128`** — `os.JSONpath.exists(...)` — טעות הקלדה (`os` אין לו `JSONpath`), צ"ל `os.path.exists`. כרגע לא נפגע כי `addActions` נקרא תמיד עם `dict`, אבל ה-API המתועד (string/קובץ JSON) שבור לגמרי.

6. **`hera/riskassessment/riskToolkit.py:274` (`analysis.getRiskAreas`)** — קורא ל-`calculateRaw(...)` שלא קיים בשום מקום ב-`Injury`/`InjuryLevel`. יקרוס תמיד אם ייקרא.

7. **`hera/riskassessment/analysis/riskAreas.py:21-24` (`getRiskAreaAlgorithm`)** — משתמש ב-`pydoc.locate("pyriskassessment...")` — שם חבילה ישן מלפני השינוי ל-`hera` — תמיד יחזיר `None` ויזרוק `ValueError`.

8. **`hera/simulations/openFoam/eulerian/abstractEulerianSolver.py`** — צביר שלם של מתודות ש"ילכו" ל-`NameError`/`AttributeError` ברגע שייקראו (`flowType` שורה 30 מפנה ל-`SIMULATIONTYPE_COMPRESSIBLE` הלא קיים; `blockMesh_setBoundFromFile`/`blockMesh_setDomainHeight` שורות 88/111/116 מפנים לפרמטרים לא קיימים; קוראים ל-`set_blockMesh_boundaries` שלא קיים בכלל ב-`OFWorkflow.py`). API "ציבורי" שמעולם לא רץ בפועל.

9. **`hera/simulations/openFoam/OFWorkflow.py:74-83`** — `@workflowGroup.setter` (צ"ל `@workflowGroupID.setter`) — דורס בטעות את ה-property `workflowGroupID` המקורי (הופך אותו ל-unreachable).

10. **`hera/simulations/openFoam/postProcess/pvOpenFOAMBase.py`** — שיטת `write_parquet` מוגדרת פעמיים זהה (שורות 631-634 ו-637-640) — ההגדרה הראשונה dead-on-arrival. שתיהן קוראות לפונקציה `writeNonRegularCase` שלא קיימת (צ"ל `self.writeCase`) — יקרוס תמיד. גם `write_netcdf`, `to_pandas`, `to_xarray`, `to_dataFrame`, `to_dataArray` (כולם `@deprecated`) קוראים למתודות/פרמטרים לא קיימים — שבורים לגמרי.

11. **`hera/simulations/machineLearningDeepLearning/torch/modelContainer.py:211-222, 266-277` (`fit`/`only_validate`)** — `ckpt_path = None`, ואז `elif os.path.exists(ckpt_path):` נבדק כש-`ckpt_path` עדיין `None` → `TypeError`. זה קוד **חי** (לא רק "מת") עם באג אמיתי — worth fixing, not deleting.

---

## 1️⃣ פונקציות/מחלקות/מתודות/משתנים לא בשימוש (Unused Definitions)

### datalayer / utils / toolkit.py / bin
| קובץ:שורה | מה | למה מת |
|---|---|---|
| `datalayer/collection.py:171` | `AbstractCollection.addDocumentFromJSON` | grep — hit יחיד (ההגדרה) |
| `utils/query.py:1` | `andClause(excludeFields=[], **kwargs)` | 0 קריאות; גם באג: mutable default argument |
| `utils/zipUtils.py:9` | `add_directory_to_zip` | 0 קריאות; אפילו `zip_items` באותו קובץ לא משתמש בו — כפילות לוגיקה |
| `utils/unitHandler.py:304` | `extractUnumUnitsFromPint` | 0 קריאות; כפילות מלאה של גוף `pintToUnum` |
| `utils/unitHandler.py:357` | `unumToBaseUnits` | 0 קריאות, לא ב-`__all__` |
| `toolkit.py:128` | `abstractToolkit.classLoggerName` (property) | 0 קריאות |
| `bin/hera-OF-postProcess.old` | קובץ שלם | מסונן בכוונה ע"י `setup.py` (סיומת `.old`) — קוד מת בפועל בעץ |
| `bin/jupyter-lab-server` | קובץ CLI | **תקול-רישום**: תבנית ה-glob ב-`setup.py:19` (`hera-*`) לא תופסת אותו כי אין לו קידומת `hera-` → אף פעם לא מותקן, למרות שמתועד |

### GIS (raster + vector)
- `raster/hill2stl.py` — **קובץ שלם** מת: מבצע `generate_solid_stl(...)` בזמן import (side-effect מסוכן — כותב `test1.stl` לדיסק אם מישהו יעשה import).
- `raster/landcover.py:561-585` — `roughnesslength2sandgrainroughness` מוגדר **פעמיים** (גם ב-717) — ההגדרה הראשונה unreachable.
- `raster/landcover.py:488-522` — `_handleType1` — 0 קריאות בכלל.
- `raster/tiles.py` — `tileScaleAtLatLonZoom`, `listImages`, `setDefaultTileServer` — 0 קריאות, אין `test_tiles.py`.
- `vector/buildings/toolkit.py:202-238,241-292` — `get_buildings_height`, `filter_buildings_in_area` (סטטיות) — 0 קריאות.
- `vector/toolkit.py:9,37-57` — `TOOLKIT_VECTOR_REGIONNAME`, `geopandasToGeoJson` — 0 קריאות.
- `vector/topography.py:99-113,158-217` — `geoPandasToSTL`, `toDEM` — 0 קריאות.
- `vector/demography.py:113-119` — `projectPolygonOnPopulation` — shim מיושן (`DeprecationWarning`), 0 קריאות — התיעוד עצמו אומר להסיר.

### meteorology (highfreqdata)
- `toolkit.py:26-27` — `DOCTYPE_STATIONS`, `DOCTYPE_MEASUREMENTS` — 0 קריאות.
- `analysis/meandatacalculator.py:197-198` — `_UV_to_SpdDir` — 0 קריאות.
- `analysis/turbulencestatistics.py:1352-1425` — מחלקה שלמה `SinglePointStatisticsSpark` — 0 ייבוא/שימוש בשום מקום.
- `analysis/turbulencestatistics.py:1428-1549` — `InMemoryRawData.append/read_hdf/to_hdf` — 0 קריאות (המחלקה היורשת `InMemoryAvgData` בשימוש, אך לא דרך מתודות אלה).
- `parsers/CampbellBinary.py` — `instrument`, `station`, `firstTime`, `lastTime` (properties), `getRecordByTime` — כולם 0 קריאות; הקוד החי קורא ישירות ל-`headers[0].split(...)` במקום דרך ה-properties.
- `analysis/abstractcalculator.py:62-65` — `JoinMethod` property — 0 קריאות.

### experiment
- `dataEngine.py` — `pandasDataEngineDB`, `daskDataEngineDB` (מחלקות שלמות) — נוצרות רק דרך factory שתמיד מקבל את ברירת המחדל `PARQUETHERA`; שתיהן גם שבורות פנימית (`self.experimentObj` לא מוגדר אף פעם, `dask_mongo` לא מיובא).
- `analysis.py:29-52,232-254` — `getDeviceLocations`, `getDeviceTypePlannedMessageCount` (קורא ל-`getOptimalFrequencyHz` שלא קיים) — 0 בדיקות אמיתיות.
- `presentation.py:219-248,449-497,648-700` — `plotMap`, `plotDevices`, `generateLatexTable` — משתמשים בפרמטרים/attributes לא קיימים (`self.trialSet` צ"ל `self.datalayer.trialSet` וכו') — שבורים וגם לא נקראים.
- `parsers.py:356-372` — `Parser_TOA5` — 0 קריאות, גוף ריק (`pass`).
- `experiment.py:221-225,382-399` — `experimentDataType` (תמיד `None`), `_initAnalysisAndPresentation`, `trialsOfDefaultTrialSet` — 0 קריאות.

### riskassessment
- שלוש מחלקות factory-eligible שמעולם לא נקראות דרך JSON/config: **`CalculatorHaber`**, **`InjuryExponential`**, **`InjuryLevelExponential`** (אומת: 0 מופעים של המחרוזות `"Haber"`/`"Exponential"` בכל JSON/config בריפו).
- `Agents.py:60-80` — `Agent.fullDescription`, `Agent.effectproperties` — 0 קריאות.
- `Injury.py:190-195,285-357` — `_postCalculatePointWise` (גוף `pass`), `calculatePointWiseFractionInjured` — 0 קריאות.
- `ProtectionPolicy.py:205-213` — `hdfkey` (ברמת המדיניות, לא ברמת הפעולה) — 0 קריאות.
- `riskToolkit.py:127-152` — `listAgentsNames`, `loadAgent` — 0 קריאות, אין wiring ב-CLI.
- `analysis/riskAreas.py` — כל המודול (`getRiskAreaAlgorithm`, `riskAreaAlgorithm_Sweep`) — 0 קריאות ושבור (ראה עדיפות דחופה #7).

### simulations (הרשימה הגדולה ביותר — ראו גם "עדיפות דחופה" למקרים השבורים)
- `analysis/errorCalculation.py` — מחלקה שלמה `errorCalculation` (7 מתודות) — 0 שימוש בכל הריפו.
- `gaussian/DropletCloud.py`, `gaussian/MeshUtils.py` — קבצים שלמים, 0 ייבוא בשום מקום.
- `hydrodynamics/nearWallFlow.py` — קובץ/חבילה שלמה, לא רשומה ב-`toolkitHome` בכלל.
- `gaussian/gasCloud.py:668-716` — מחלקת `Continuous` ("Yehuda's Code") — 0 שימוש.
- `gaussian/Meteorology.py:317-322` — `MeteorologyProfile` — stub ריק (`pass`), לא רשום ב-factory.
- `windProfile/toolkit.py:115` — `_getStationsInRegion` — 0 קריאות, גם עם באג נתיב יחסי.
- `CLI.py:480,525` — `workflowNodes_list`, `workflowNodes_listParameters` — לא מחוברים ל-`hera-workflows` CLI.
- `WRF/wrfDatalayer.py:16` — `wrfDatalayer` — מיוצא מ-`__init__.py` אך לא נרשם ב-`toolkit.py`, לא נוצר בשום מקום. **דורש אימות ידני**.
- `LSM/hermesWorkflowToolkit.py` — **קובץ שלם (783 שורות)**: `workflowToolkit` הוא fork ישן ומיושן של `hera/simulations/hermesWorkflowToolkit.py` שאיש לא יורש ממנו (`LSMToolkit` יורש ישירות מ-`abstractToolkit`).
- `deposition/models.py`, `evaporation/models.py`, `evaporation/monaghan.py` — מחלקות `depositionModels`, `evaporationModels`, `MonaghanConstantConditions` — 0 שימוש חיצוני; `monaghan.py` גם מייבא `pyriskassessment` שלא מותקן כלל.
- `machineLearningDeepLearning/dataanalisys/` — **חבילה שלמה** (`stat()`, מחלקת `ml`) — 0 ייבוא בשום מקום.
- `machineLearningDeepLearning/toolkit.py:323-415` — מחלקת `analysis` (`sensitivityAnalysis_morris`) — אף פעם לא מאותחלת ב-`__init__`.
- `LSM/toolkit.py:33-34`, `LSM/template.py:26-28` — קבועים `TRUE`/`FALSE`, `STABILITY_*` — 0 קריאות (הקוד החי משתמש ב-string literals ישירות).
- `openFoam/toolkit.py` — שורת ארוכה של מתודות ללא קוראים: `template_add` (497-520, גוף `pass`), `clearVTKPipelineCache`, `getVTKPipelineCacheTable`, `getHermesWorkflow_Flow`, `getMeshFromName`, `getMeshExtentFromName`, `xarrayToSetFieldsDictDomain` (התיעוד עצמו אומר "Not debugged").
- `openFoam/postProcess/VTKpipelineExecutionContext.py` — **מחלקה שלמה** לא בשימוש (0 הפניות מלבד ההגדרה עצמה); כפילות מיותרת ל-`VTKPipeline.py`.
- `openFoam/postProcess/VTKPipeline.py:201-619` — מחלקת `registeredVTKPipeLine` (~420 שורות) ו-`registerPipeline` — 0 קריאות; `toolkit.py:814` בונה רק `VTKPipeLine` הבסיסי בלי לקרוא ל-`registerPipeline`.
- `openFoam/preprocessOFObjects/OFList.py` — **מחלקה שלמה**, מיובאת ב-`__init__.py` אך לעולם לא נוצרת; גם שבורה פנימית (`self.columnNames`/`getHeader()` לא קיימים).
- `openFoam/preprocessOFObjects/OFObject.py:22-23,59-81` — `REGION_INTERNSALFIELD` (יש גם טעות הקלדה בשם), `internalField()`, `.processors`, `.processorItems`, `.dimensionsStr` — 0 שימוש.
- `openFoam/lagrangian/abstractLagrangianSolver.py` — `getOriginalFlowFieldMesh`, `getDispersionDocument`, `getDispersionFlowDocument` — 0 קריאות.
- `openFoam/lagrangian/LSM/toolkit.py:699-730` — `createRootCaseMeshLink` — 0 קריאות.
- `eulerian/simpleFoam.py:4` — `simpleFoam_toolkitExtension` — לא מוגדר/מחובר ב-`OFToolkit.__init__` (בניגוד ל-`buoyantReactingFoam_toolkitExtension` שכן מחובר).
- `openFoam/toberewritten/netcdf2of.py`, `xarrayDataset2OF.py` — **שגיאות תחביר ממש** (`SyntaxError`/`TabError`), לא ניתנים לטעינה כלל. `wrf2of.py` — לא שבור תחבירית אך עם נתיבים מוחלטים קשיחים (`/data5/NOBACKUP/...`) ואף אחד לא מייבא אותו. (`toberewritten/utils.py` דווקא **כן בשימוש** ע"י `LSM/toolkit.py` — לא למחוק).

---

## 2️⃣ קוד שהפך להערה (Commented-Out Code) — לא נחוץ יותר

בלוקים משמעותיים (מעל כמה שורות) שמומלץ למחוק:

- `hera/simulations/openFoam/postProcess/pvOpenFOAMBase.py:643-718` — מתודת `writeRegularCase` שלמה בהערה (~76 שורות).
- `hera/simulations/openFoam/preprocessOFObjects/utils.py:52-270` — בלוק ענק (~219 שורות): `emptyParallelField`, `extractFile`, `extractBoundaryField` — גם לא שלם תחבירית אם יבוטל ההערה.
- `hera/simulations/openFoam/OFWorkflow.py:445-935` — בלוק ענק (~490 שורות) של פונקציות handler ישנות ל-geometry/blockMesh.
- `hera/simulations/openFoam/CLI.py:772-981` — עותק ישן שלם (~210 שורות) של `stochasticLagrangian_dispersion_create`, כפול לגרסה החיה מעליו.
- `hera/simulations/openFoam/lagrangian/abstractLagrangianSolver.py:1879-2056` — מימוש ישן (~178 שורות) של `_extractFile`/`_readRecord`, הוחלף בפונקציות חדשות.
- `hera/simulations/openFoam/lagrangian/LSM/toolkit.py:1190-1324` — כולל **שאריות conflict-marker של git merge שלא נוקו** (`# >> >> >> > b93a26...`) — עדות לresolve מבולגן שלא טופל.
- `hera/simulations/CLI.py:434-446` — ענף שלם בהערה שסותר את הלוגיקה החיה מתחתיו.
- `hera/simulations/machineLearningDeepLearning/dataanalisys/ml.py` — כמה בלוקים גדולים (שורות ~38-51, 161-171, 454-479, 487-521) של קוד ניסיוני מבוטל.
- `hera/measurements/GIS/utils.py:456-587` — **חמש** פונקציות שלמות בהערה: `ITMtolatlong`, `PolygonDataFrameIntersection`, `getBoundaries`, `makePolygonFromEndPoints`, `ConvexPolygons` (האחרונה כפילות ישנה של גרסה חיה ב-`buildings/analysis.py:33`).
- `hera/measurements/experiment/presentation.py:703-784` — שתי מתודות שלמות בהערה: `plotNDIRFrequencyDistribution`, `plotMessageFrequencyDistribution`.
- `hera/measurements/meteorology/highfreqdata/analysis/meandatacalculator.py:452-472` — מתודת `zoL` שלמה בהערה, הוחלפה ב-`zOverL`.
- `hera/riskassessment/agents/effects/InjuryLevel.py:9-12,276-279` — ייבואי `unum` ישנים (הפרויקט עבר ל-`pint`) + מימוש חלופי עם "bug" מתועד.
- `hera/datalayer/collection.py:247-249`, `hera/datalayer/document/__init__.py:22` — שאריות קטנות של קוד ישן.

---

## 3️⃣ לוגיקה תנאית שלעולם לא תתקיים (Unreachable Code)

- `hera/simulations/gaussian/gasCloud.py:325-390` — הגדרת מחלקה כפולה, נדרסת מיידית (ראו עדיפות דחופה #1).
- `hera/simulations/openFoam/CLI.py` — שלוש פונקציות מוגדרות **פעמיים**; Python שומר רק את ההגדרה האחרונה: `foam_mesh_blockMesh` (259 מול 640), `foam_mesh_setDomainHeight` (263 מול 654), `IC_hydrostaticPressure` (267 מול 708) — ההגדרות הראשונות מתות לחלוטין.
- `hera/measurements/GIS/vector/buildings/analysis.py:201-202` — `_LambdaF = 0`/`_LambdaP = 0` כתכונת מחלקה, נדרסות ע"י הגדרות `def` מאוחרות יותר של אותם שמות (385, 494) — הערכים הראשוניים בלתי-נצפים לעולם.
- `hera/measurements/GIS/vector/topography.py:50-95` — `cutRegionFromSource` קורא ל-`super()...` עם keyword arguments שלא תואמים לחתימת ה-parent (יזרוק `TypeError` בכל קריאה); הענף השני קורא ל-`getRegionData` שלא קיים בכלל בריפו.
- `hera/simulations/hydrodynamics/nearWallFlow.py:174,298` — `Lambda = 2*numpy.log()` — קריאה בלי ארגומנטים, תמיד `TypeError`.
- `hera/simulations/machineLearningDeepLearning/dataanalisys/ml.py:117,154,252` — שלושה בלוקים גדולים מאחורי תנאים קבועים כמו `if 4==5:`/`if 5==6:` — לעולם לא ירוצו.
- `hera/simulations/LSM/template.py:449,458` — `getSimulationByID`/`getSimulationByName` קוראים למתודות שלא קיימות על `LSMTemplate` (`getDocumentByID`/`getSimulationsDocuments` — הכוונה ל-`self.Toolkit....`) — `AttributeError` מובטח.
- `hera/simulations/openFoam/eulerian/abstractEulerianSolver.py:30` — property `flowType` מפנה לקבוע לא קיים.
- ראו גם סעיף "עדיפות דחופה" למקרים החמורים ביותר בקטגוריה זו.

---

## 4️⃣ קבצי בדיקות שבודקים רכיבים שנמחקו/לא בשימוש

**ממצא טוב:** לא נמצאו imports שבורים ממש בקבצי הבדיקות (`hera/tests/`) — לכל import שנבדק יש מחלקה/פונקציה תואמת בקוד המקור.

הערה אחת: `hera/tests/test_datalayer.py:959` — `@pytest.mark.xfail(reason="Project.load references _iter_pickled_docs which is not implemented")` — הסיבה **לא מדויקת יותר**: `_iter_pickled_docs` כן מומש ב-`project.py:342`. ה-xfail עצמו לא "שבור" אך התיעוד שלו מיושן — כדאי לעדכן/להסיר את ה-xfail אם הבדיקה עוברת כעת.

**פער כיסוי משמעותי (לא "שבור" אלא "לא נבדק בכלל"):** כל תיקיות ה-`simulations/openFoam/{postProcess,preprocessOFObjects,toberewritten}` — **אין להן שום קובץ טסט** בריפו. זה מסביר מדוע כל-כך הרבה מהבאגים שנמצאו שם (שיטות שקוראות לפונקציות לא קיימות) מעולם לא התגלו.

**ממצא נוסף (לא "בדיקה שבורה" אלא באג מקור סמוי):** ב-`hera/datalayer/autocache.py` — אם המודול מיובא ישירות לפני שה-lazy loading הפנימי של `hera/__init__.py` (`_load_deferred()`) הופעל, השורה `from hera import Project` בתוך `autocache.py` גורמת לקריאה חוזרת ל-`_load_deferred()` שמנסה לייבא מחדש את `autocache` בזמן שהוא עדיין `partially initialized` → `ImportError: cannot import name 'cacheFunction'...`. בפועל זה לא מתפוצץ היום כי שלושת הטסטים הרלוונטיים (`test_datalayer.py:898,917,931`) מסומנים `@requires_mongo` ומדולגים בסביבה בלי MongoDB חי — אבל זה יקרוס עם MongoDB אמיתי. דורש תיקון ב-`hera/datalayer/autocache.py` (למשל import מקומי/lazy במקום import ברמת המודול).

---

## 5️⃣ משתני סביבה / הגדרות קונפיג שלא נקראים בשום מקום

| משתנה | נקרא ב- | מוגדר איפשהו? |
|---|---|---|
| `HERA_REPO_ROOT` | `toolkit.py:771`, `experiment.py:141` | ✅ כן (`set_hera_environment.sh`, `activate_hera.sh`) |
| `TEST_HERA` | `conftest.py:223` | ✅ כן (`.github/workflows/ci.yml:22`) |
| `HERA_FULL_LOGGING_TESTS` | `logging/test_helpers.py:16` | ⚠️ מתועד בלבד, לא מוגדר ב-CI — טסטים המותנים בו אף פעם לא רצים ב-CI |
| `RESULT_SET` | `conftest.py:245` | ⚠️ מתועד בלבד, לא מוגדר ב-CI |
| `GDF_TOL_AREA` | `conftest.py:417` | ⚠️ מתועד בלבד, לא מוגדר ב-CI |
| **`PYARGOS_PATH`** | `tests/dynamic_loading_tests_pack/test_experiment_cli_shortcuts.py:58` | ❌ **לא מוגדר בשום מקום** בריפו (לא ב-CI, לא ב-Dockerfile, לא ב-docs, לא ב-.env) — הכי "מת" מבין כולם |

קבועים לא נקראים:
- `hera/measurements/meteorology/highfreqdata/toolkit.py:26-27` — `DOCTYPE_STATIONS`, `DOCTYPE_MEASUREMENTS`.
- `hera/measurements/GIS/vector/toolkit.py:9` — `TOOLKIT_VECTOR_REGIONNAME`.
- `hera/simulations/LSM/toolkit.py:33-34`, `hera/simulations/LSM/template.py:26-28` — קבועי `TRUE`/`FALSE`/`STABILITY_*`.
- `hera/simulations/openFoam/preprocessOFObjects/OFObject.py:22-23` — `REGION_INTERNSALFIELD`/`REGION_BOUNDARYFIELD` (יש גם טעות כתיב בשם הראשון).
- `hera/simulations/openFoam/postProcess/pvOpenFOAMBase.py:12,31` — ייבואים לא בשימוש: `CASETYPE_RECONSTRUCTED`, `hera_logging`.

---

## סיכום המלצות לפי סדר עדיפות

1. **לתקן מיד** (קוד חי שבור): `gasCloud.py` (m/min), `thresholdGeoDataFrame.py` (collections.Iterable), `ProtectionPolicy.py` (os.JSONpath), `modelContainer.py` (ckpt_path), `abstractEulerianSolver.py`.
2. **למחוק תיקיות/קבצים שלמים כפולים ושבורים**: `simulations/utils/interpolation/`, `simulations/LSM/hermesWorkflowToolkit.py`, `openFoam/postProcess/VTKpipelineExecutionContext.py`, `openFoam/preprocessOFObjects/OFList.py`, `openFoam/toberewritten/{netcdf2of,xarrayDataset2OF}.py`, `GIS/raster/hill2stl.py`.
3. **למחוק בלוקי הערות ישנים** (סעיף 2 למעלה) — ניקוי קל וללא סיכון.
4. **לבדוק ידנית לפני מחיקה** (API ציבורי שאולי משמש קוד חיצוני): `wrfDatalayer`, מתודות ה-presentation/analysis הלא-נבדקות, `PYARGOS_PATH`.
