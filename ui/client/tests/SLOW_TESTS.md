# Slow UI client tests - map

Measured 2026-09-15 on the `ui-projname` branch. 48 files, 424 tests.

All "4-core" numbers come from `taskset -c 0-3 npx vitest run`, which matches a
GitHub `ubuntu-latest` runner. "Solo" numbers come from running that one file
alone on an idle 20-core machine.

## Headline

| Where the time goes (4 cores) | Seconds |
|---|---|
| import | 86.7 |
| tests | 45.5 |
| environment (jsdom setup) | 31.1 |
| **wall clock** | **42.8** |

Those first three are summed across 4 workers, so they overlap.

Read that table before optimising single tests. Import is the biggest block, and
no change to an individual test touches it. The 30 slow tests below are 23.7s of
the 45.5s spent in `tests`, so roughly a quarter of the wall clock.

## The 30 tests over 500ms

Grouped by cause. Cause is measured, not guessed, except where marked.

### Cause A - jsdom parses the emotion stylesheet on the first query (18 tests)

This is the single biggest per-test cost, and it is almost always the **first
rendering test in a file**.

Probe on `AgentConfigEditor`, a 395-node tree:

```
render                          = 167.0ms
getByRole(spinbutton, {name})   = 504.1ms   <- first call
getByRole(spinbutton, {name})   =   9.3ms   <- same call again
getByRole(spinbutton) no name   =  13.6ms
getByRole(spinbutton) hidden    =   1.2ms
getByLabelText                  =   7.0ms
querySelector('input')          =   0.5ms
getComputedStyle on all 395     = 834.8ms
```

The first call costs 55x the second. MUI/emotion injects 222 separate
stylesheets into the document. Any `*ByRole` query filters out invisible
elements, which calls `getComputedStyle` on every candidate, which makes jsdom
parse and cascade all 222 sheets. After that it is cached, so every later query
in the same file is cheap.

The tax is paid once per test file, because isolation gives each file a fresh
document.

Tests hit by this, 4-core ms (solo ms):

| ms | solo | file | test |
|---|---|---|---|
| 1110 | 360 | TwoUiSync | two UIs viewing the same document both see the same data |
| 1062 | 483 | agentConfig | renders Ten Berge Coefficient field |
| 1005 | 358 | agentConfig | prevents adding duplicate effect names |
| 950 | 321 | addDocument | opens dialog when clicked |
| 882 | 311 | TwoUiSync | a second UI sees the new value after its project updates |
| 795 | 313 | addDocument | creates a regular document with name and resource |
| 784 | 329 | addDocument | creates an agent document with effects resource |
| 773 | - | deleteProject | opens confirmation dialog on click |
| 730 | - | deleteDocumentButton | calls execPython to delete document on confirm |
| 712 | 515 | agentConfig | calls setAgentResource when Ten Berge Coefficient changes |
| 686 | - | WorkflowFlowNode | updates only the type when the field is typed into |
| 666 | 236 | agentConfig | renders empty state with add controls |
| 648 | - | runWorkflow | starts the run and marks the button running |
| 598 | - | deleteDescField | removes the field from desc while keeping hidden desc fields |
| 574 | 273 | addDocument | shows resource field with auto-generated path when notebook selected |
| 563 | 257 | addDocument | updates the store after successful creation |
| 538 | 313 | addDocument | creates a notebook document with resource derived from name |
| 529 | 170 | TwoUiSync | clicking the refresh button in the tree updates the details panel |

The declining pattern inside a file is the giveaway. `deleteDocumentButton`
runs 730, 298, 166, 50, 188. `runWorkflow` runs 648, then nothing above 170.

### Cause B - real 500ms timer, not faked (2 tests) - FIXED

| ms | solo | file | test |
|---|---|---|---|
| 507 | 506 | workflowRunPoller | keeps polling while running, then stops on done |
| 522 | 505 | workflowRunPoller | writes partial chunks to the store while still running |

These were the only tests whose solo time equalled their 4-core time. They were
not CPU bound, they were waiting.

`WorkflowRunPoller.tsx:8` sets `POLL_MS = 500` and line 62 does
`setTimeout(poll, POLL_MS)`. Both tests needed a second poll, so both sat
through a real 500ms. The sibling test that needs only one poll ran in 18ms.

Fixed by calling `vi.useFakeTimers()` in those two tests and stepping the clock
with `await vi.advanceTimersByTimeAsync(POLL_MS)`, plus `vi.useRealTimers()` in
`afterEach`. The other three tests in the file keep real timers, because they
finish on the first poll and never wait. `POLL_MS` is now exported so the test
advances by the real interval instead of a copied number.

The whole file went from 1054ms to 24ms. Checked that it is not passing
vacuously: advancing by `POLL_MS - 1` makes both tests fail.

### Cause C - large or repeated MUI renders (10 tests)

No single dominant factor found. These render big trees, sometimes several
times in one test, and pay Cause A on top.

| ms | solo | file | test |
|---|---|---|---|
| 1965 | 783 | WorkflowContextMenu | references another node output via two fly-out submenus |
| 1233 | 421 | projectViewSettings | toggles "Always save workflow before running" through the store |
| 1226 | 429 | addProject | creates the project, then scans it for notebooks in a separate call |
| 1116 | 341 | itemTypeObject | turns a scalar field into an empty substructure when set to object |
| 746 | - | fillProjectName | fills a desc field renamed to projectname |
| 668 | - | centralRepoFolder | opens settings dialog and re-fetches with new folder |
| 570 | - | detailsViewRepo | does not call execPython for temp repos |
| 567 | 216 | projectViewSettings | reflects the current store value when opened |
| 524 | - | repoTreeAddButton | calls execPython with repo data on confirm |
| 506 | - | projectActions | opens a popover with the Actions title and both actions |

`WorkflowContextMenu` is the worst single test. It opens two nested MUI
Autocomplete fly-outs and does four `findBy*` queries, each re-running the
accessibility scan against a menu that keeps growing.

## Contention multiplier

Everything except Cause B runs about **3x slower on 4 cores** than solo. That is
not a test problem, it is 4 workers on 4 cores. It is why the same suite reads
25s locally and 78s on CI.

## What was tried and rejected

`configure({ defaultHidden: true })` makes `*ByRole` skip the visibility filter,
which should sidestep Cause A. All 424 tests still passed, but wall clock went
42.8s -> 43.4s, and `tests` only dropped 45.5s -> 42.9s. Not worth the semantic
change. MUI calls `getComputedStyle` during render anyway, so the stylesheet
parse happens with or without the query.

## Open leads, in order of expected payoff

1. **`isolate: false`.** Attacks import (86.7s) and Cause A at once, since files
   sharing a worker share both the module registry and the parsed styles.
   Locally: import 52.6s -> 9.8s, wall 25.4s -> 10.2s. Blocked: 17-29 tests fail,
   the count varying by run. 11 test files mock `../src/io/fetchPython` with
   their own `vi.fn()`, and with a shared registry whichever file loads first
   wins, so the others see 0 calls. Needs a shared mock helper first.
2. ~~**Fake timers in `workflowRunPoller.test.tsx`.**~~ Done, see Cause B.
   1054ms -> 24ms.
3. **`// @vitest-environment node`** on the 17 logic-only `.ts` test files.
   Measured on those files alone: 2.36s -> 0.61s, environment cost 6.15s -> 1ms.
   No risk, they touch no DOM.
4. **Check the runner's core count.** `maxWorkers: 4` in `vite.config.ts:26` is
   hardcoded. If CI gives 2 vCPU, 4 workers oversubscribe it.
