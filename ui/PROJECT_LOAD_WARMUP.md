# Handling half-loaded project data on startup (design plan)

Status: **not built** — plan for review.

## Problem

Issue #1011: the UI can come up blank / half-working when the server's
`toolkitHome` is still warming up (Mongo connecting). The startup fetch in
`FetchProjects.tsx` bundles two Mongo queries into one `/exec` call:

```python
toolkitDocs = toolkitHome.getToolkitDocuments()          # query 1
docs = All.getDocumentsAsDict('<proj>', with_id=True)    # query 2
```

While the DB is cold this returns **half data**:

- **(b)** toolkit docs missing `desc` / `desc.classpath` — and `parseToolkits`
  crashes reading `desc.classpath`.
- **(a)** documents whose `desc.toolkit` points at a toolkit that isn't in the
  loaded set.

We do NOT want to split the fetch (no extra server round-trips / friction). We
fuse what's usable, mark the data as "not warm", and retry.

## Why not the obvious fix

A `setTimeout(retry, 5000)` placed inside the async `fetchProjectData` is racy:
it has **no owner and no cleanup**. `fetchProjectData` is called from the URL
effect, again by React StrictMode's double-invoke, and again by its own retry —
so timers stack up and nothing can cancel a stale one when the project changes.

The fix is to give the timer a single owner with cleanup: a `useEffect`.

## Plan: separate read → analyze → effect

Three clearly separated pieces.

### 1. Reader (I/O only, no logic)

`fetchProjectDataRaw(projectName)` — issues the same single bundled `/exec` and
returns the **raw** response (`{ toolkitDocs, project }`) or `null` on
problem/network error. No filtering, no store writes, no timers.

The raw result is held in **local state of `FetchProjects`** (`useState`), keyed
so a stale project's raw data can be ignored:

```ts
const [raw, setRaw] = useState<{ projectName: string, data: any } | null>(null);
```

Everything *derived* still lands in the zustand store — local state holds only
the transient raw payload, so there's no second source of truth for good data.

### 2. Analyzer (pure function, exported, unit-tested)

`analyzeProjectData(raw) => { toolkits, project, ok }`

- **(b)** keep only toolkit docs with a `toolkit` name and `desc.classpath`;
  drop the rest. If any dropped → not `ok`.
- **(a)** if any document's `desc.toolkit` matches no loaded toolkit → not `ok`.
  We don't chase the missing toolkit; its absence is just the signal.
- `ok = !dropped && !referencesMissingToolkit`.

Pure and synchronous → tested directly with no mocks or fake timers.

### 3. Effect (owns store writes + the single retry timer)

A `useEffect` in `FetchProjects` keyed on `raw`:

```ts
useEffect(() => {
  if (!raw || raw.projectName !== currProjectName) return;   // ignore stale
  const { toolkits, project, ok } = analyzeProjectData(raw.data);
  setToolkits(toolkits);
  setCurrentProject(project);                                 // fuse: show what we have
  if (ok) return;
  const timer = setTimeout(() => { fetchProjectDataRaw(currProjectName).then(setRaw); }, 5000);
  return () => clearTimeout(timer);                           // cleanup kills the race
}, [raw, currProjectName]);
```

The cleanup cancels a pending retry whenever `raw` changes, the project changes,
or the component unmounts. **One owner, one live timer** — no stacking.

## Flow

```
URL effect ──► fetchProjectDataRaw ──► setRaw(local state)
                                          │
                                   analyze effect
                                    ├─ store: toolkits + project (fused)
                                    └─ if !ok: setTimeout(retry, 5s) + cleanup
                                          │ (retry) ──► fetchProjectDataRaw ──► setRaw ─┐
                                          └──────────────────────────────────────────┘
                                    loop self-terminates once ok
```

## Files touched

- `src/io/FetchProjects.tsx` — split `fetchProjectData` into
  `fetchProjectDataRaw` (reader) + `analyzeProjectData` (pure); add raw state
  and the analyze/retry effect. `parseToolkits` folds into the analyzer.
- `tests/FetchProjects.test.ts` — reader test (raw passthrough), analyzer tests
  (drop incomplete docs, missing-reference → not ok, warm → ok), effect/retry
  tests with fake timers (retry once then stop; cleanup on project switch).

No server changes. `fetchProjectsNames` and `fetchProjectDetails` (auto-reload)
are untouched.

## Caveat to decide

An **unregistered** toolkit — a document referencing a toolkit that legitimately
has no toolkit doc (supported via `ToolkitObj.unregistered`) — reads as "not
warm" under (a) and would retry every 5s forever even when Mongo is warm.
Options:
- treat a missing reference as not-warm **only when zero toolkits loaded**
  (strong "cold" signal), or
- cap the retries (e.g. stop after N attempts).

Recommendation: start with the "zero toolkits loaded" guard for (a); keep (b)
(dropped incomplete docs) as an always-valid not-warm signal.

## Effort

Small–medium. One file restructured, tests updated. Pure analyzer makes the
logic cheap to test; the effect isolates the only timing-sensitive part.
