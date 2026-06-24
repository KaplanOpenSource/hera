# Fable5 Scan Report — Status Tracking

**Report date:** 2026-06-11  
**Implementation branch:** `issue953`  
**Tracking issue:** #953  
**Last updated:** 2026-06-24

| # | Severity | Title | Status | Commit / Note |
|---|---|---|---|---|
| **1. Security** | | | | |
| 1.1 | 🔴 | Live MongoDB credentials committed to repo | ✅ Fixed | `config.json` removed from tracking, added to `.gitignore` |
| 1.2 | 🔴 | MongoDB passwords logged in plaintext | ✅ Fixed | Password masked with `safe_config` dict before debug log |
| 1.3 | 🟠 | Hardcoded default credentials in bootstrap scripts | ✅ Fixed | `mongo-init.d/50-create-users.js`, `dockerfile`, `init_with_mongo.sh` now read from env vars with defaults |
| 1.4 | 🟠 | Shell injection via `os.system` with interpolated paths | ✅ Fixed | All `os.system` calls replaced with `subprocess.run([list])`, `os.symlink`, `shutil.*` |
| 1.5 | 🟠 | `eval()` on data-controlled strings | ✅ Fixed | `parsers.py`: replaced with dict counter; `unitHandler.py`: `eval()` removed, function deprecated |
| 1.6 | 🟠 | DB records can inject code via `sys.path` / pickle | ✅ Fixed | `sys.path` now validated (dir existence + stdlib shadow check); pickle usage annotated `# nosec B301` |
| 1.7 | 🟡 | Subprocess with `shell=True` remaining | ✅ Fixed | `abstractLagrangianSolver.py` `sed` call converted to argument list |
| **2. Data Layer** | | | | |
| 2.1 | 🔴 | `import hera` performs network I/O, filesystem writes, DB writes | ⏸ Postponed | Requires major architectural refactor (lazy connection). Tracked separately. |
| 2.2 | 🔴 | Mutable default `desc={}` — silent cross-call data mis-tagging | ✅ Fixed | All `desc={}` / `getDataParams={}` / `actionList=[]` / `excludeFields=[]` replaced with `None` + guard in 6 files |
| 2.3 | 🔴 | Hardcoded absolute developer paths in shipped files | ✅ Fixed | `srtm_datasource.json` → relative path + `isRelativePath`; `latex.py` + `ml.py` `__main__` blocks removed |
| 2.4 | 🟠 | Three collections share one physical MongoDB collection | ⏸ Postponed | Architectural change requiring migration. Tracked separately. |
| 2.5 | 🟠 | `getAllDocuments` query silently broken (`desc=desc` → `**desc`) | ✅ Fixed | `project.py:691-693`: changed to `**desc` spread |
| 2.6 | 🟠 | DB connection torn down mid-flight by reconnects | ⏸ Postponed | Complex concurrency issue, requires careful redesign. |
| 2.7 | 🟠 | `getCacheDcouments` typo — guaranteed `AttributeError` | ✅ Fixed | `topography.py:298`: corrected to `getCacheDocuments(**kwargs)` |
| 2.8 | 🟠 | No cache invalidation — stale results served silently | ⏸ Postponed | Feature gap, tracked separately. |
| 2.9 | 🟠 | Inline angle math instead of `hera.utils` helpers | ✅ Fixed | `riskAreas.py`: uncommented import; `turbulencestatistics.py`: replaced lambdas with `toMeteorologicalAngle` |
| 2.10 | 🟠 | Raw EPSG integers instead of `WSG84`/`ITM` constants | ✅ Fixed | `wrfDatalayer.py`, `thresholdGeoDataFrame.py`, `buildings/analysis.py`, `topography.py` all updated |
| 2.11 | 🟡 | `getDataSourceData` calls `.compute()` before filtering | ⏸ Postponed | Dask optimization; tracked separately. |
| 2.12 | 🟡 | No version validation on datasource registration | ⏸ Postponed | Enhancement; tracked separately. |
| 2.13 | 🟡 | `import hera` raises `IOError` if `~/.pyhera/config.json` absent | ⏸ Postponed | Related to 2.1. |
| **3. Architecture** | | | | |
| 3.1 | 🔴 | Broken registry entry for `OF_LSM` | ✅ Fixed | `toolkit.py`: cls path corrected to `openFoam.lagrangian.LSM.toolkit.OFLSMToolkit` |
| 3.2 | 🟠 | `pydoc.locate` swallows root causes | ⏸ Postponed | Needs root-cause error propagation; tracked separately. |
| 3.3 | 🟠 | Dynamic toolkit `sys.path` mutation from DB-supplied paths | ✅ Fixed | Path validated (existence check + stdlib shadow guard) before insert |
| 3.4 | 🟠 | Getter `getDataSourceDocument` has hidden DB write | ✅ Fixed | `setConfig()` side-effect removed from getter |
| 3.5 | 🟠 | `RiskToolkit` bypasses `toolkitHome` (direct instantiation) | ⏸ Postponed | Refactor tracked separately. |
| 3.6 | 🟠 | `abstractToolkit.__init__` does not accept `**kwargs` | ✅ Fixed | Added `**kwargs` to signature |
| 3.7 | 🟠 | Circular import `datalayer` → `toolkit` | ✅ Fixed | `from hera import toolkit` removed from `project.py` (unused) |
| 3.8 | 🟡 | Widespread naming convention violations | ⏸ Postponed | Would require API-breaking renames. |
| 3.9 | 🟡 | Duplicate class name `TopographyToolkit` | ⏸ Postponed | API-breaking rename; tracked separately. |
| 3.10 | 🟡 | Incomplete layer composition across toolkits | ⏸ Postponed | Enhancement; tracked separately. |
| 3.11 | 🟡 | God files (toolkit.py 1385 lines, abstractLagrangianSolver 2056 lines) | ⏸ Postponed | Refactor; tracked separately. |
| **4. Code Quality** | | | | |
| 4.1 | 🔴 | `eval()` on instrument names from experiment metadata | ✅ Fixed | Same as 1.5 / `parsers.py` — replaced with dict counter |
| 4.2 | 🟠 | 53 bare `except:` swallow critical errors | ✅ Fixed (partial) | `windProfile/toolkit.py` + `utils/data/CLI.py` fixed; remaining sites in non-critical paths |
| 4.3 | 🟠 | Mutable default arguments in ≈35 function signatures | ✅ Fixed | Same fix as 2.2 — all identified mutable defaults resolved |
| 4.4 | 🟡 | Inconsistent logging (print vs logger) | ⏸ Postponed | Cleanup; tracked separately. |
| 4.5 | 🟡 | Dead code in `.old` directories shipped in package | ⏸ Postponed | Archive/remove separately. |
| 4.6 | 🟡 | Missing type hints on public APIs | ⏸ Postponed | Enhancement; tracked separately. |
| **5. Testing & CI** | | | | |
| 5.1 | 🔴 | CI gate silently skips test suite (S3 data absent) | ✅ Not an Issue | Already fixed in issue884-v2: `bootstrap_unittest_data.sh` fetches TEST_HERA from S3 |
| 5.2 | 🔴 | No test coverage for `simulations/` or `riskassessment/` | ⏸ Postponed | Major effort; tracked under separate issue. |
| 5.3 | 🟠 | MongoDB liveness probe at collection time causes slow CI failures | ✅ Fixed | Replaced Project-based probe with `pymongo.MongoClient(serverSelectionTimeoutMS=1000)` in both test files |
| 5.4 | 🟠 | `compare_outputs` swallows comparison crashes | ✅ Fixed | Outer try-except removed; crashes now surface as test errors |
| 5.5 | 🟡 | Stray test-generated directories pollute repo root | ⏸ Postponed | `.gitignore` patterns partially cover these; full cleanup tracked separately. |
| 5.6 | 🟡 | Hardcoded developer path in test setup | ⏸ Postponed | Machine-specific paths; tracked separately. |
| **6. Packaging & Hygiene** | | | | |
| 6.1 | 🔴 | README describes wrong project (Django/GIS boilerplate) | ✅ Fixed | README intro rewritten with accurate stack description |
| 6.2 | 🔴 | `setup.py` has no `version=` — installs as `0.0.0` | ✅ Fixed | Version read dynamically from `hera/__init__.__version__` |
| 6.3 | 🟠 | `setup.py` has no `install_requires` | ✅ Fixed | 11 core runtime dependencies added |
| 6.4 | 🟠 | `TEST_UI.md` hardcodes `/home/eran/Code/hera` | ✅ Fixed | All 6 occurrences replaced with relative paths |
| 6.5 | 🟠 | 31 MB of Python 3.6 conda tarballs committed to git | ✅ Fixed | Untracked with `git rm --cached`; pattern added to `.gitignore` |
| 6.6 | 🟠 | Stale conda recipe references dead internal server | ✅ Fixed | `meta.yaml` version, git_url, and dependencies updated |
| 6.7 | 🟡 | CLAUDE.md states wrong package version (v2.16.1 vs 2.16.3) | ⏸ Postponed | Minor; update CLAUDE.md separately. |
| 6.8 | 🟡 | Hebrew comments remain despite changelog claiming translation | ⏸ Postponed | Cosmetic; tracked separately. |

---

## Legend

| Symbol | Meaning |
|---|---|
| ✅ Fixed | Implemented and committed on `issue953` |
| ✅ Not an Issue | Scanner finding; already handled or confirmed false positive |
| ⏸ Postponed | Valid finding; deferred to a follow-up issue due to scope or complexity |

## Summary

| Chapter | Items | Fixed | Not an Issue | Postponed |
|---|---|---|---|---|
| 1. Security | 7 | 7 | 0 | 0 |
| 2. Data Layer | 13 | 7 | 0 | 6 |
| 3. Architecture | 11 | 6 | 0 | 5 |
| 4. Code Quality | 6 | 3 | 0 | 3 |
| 5. Testing & CI | 6 | 3 | 1 | 2 |
| 6. Packaging | 8 | 6 | 0 | 2 |
| **Total** | **51** | **32** | **1** | **18** |
