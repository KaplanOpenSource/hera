# Server startup: fast serve + fewer exec calls (design plan)

Status: **not built** — plans for review.

The server now warms `hera` at import time (`from hera import toolkitHome` in
`server.py`). That fixed the cold-start import races, but the warm import is
slow and **blocks the whole server from starting** — the browser gets nothing
until hera finishes loading.

Two independent ideas below. They can ship separately.

---

## Plan 1 — Serve the page immediately, gate Python until ready

### Goal
Load `index.html` and the UI instantly. Python calls (`/exec`) that arrive
before hera finished warming return a clear "not ready" answer; the client
shows a spinner and retries.

### Why it works
Serving static files (`index.html`, assets) needs **no hera import**. Only
`/exec` and the workflow endpoints do. So we can boot the web server first and
warm hera in the background.

### Server changes (`ui/server/server.py`)
1. Move the warm-up off the import line into a **background thread** started at
   app startup. A module-level flag tracks state:
   ```python
   _ready = False
   def _warm():
       global _ready
       from hera import toolkitHome  # noqa: F401
       _ready = True
   threading.Thread(target=_warm, daemon=True).start()
   ```
2. Add a readiness endpoint:
   ```python
   @app.get("/ready")
   def ready(): return {"ready": _ready}
   ```
3. In `/exec`, short-circuit while not ready:
   ```python
   if not _ready:
       return ExecResponse(problem=Problem(error="WARMING_UP", traceback=""))
   ```
   Use a distinct, machine-readable error string (`WARMING_UP`) so the client
   can tell "still starting" apart from a real error.

### Client changes
- `fetchPython` (`src/io/fetchPython.ts`): when the response problem is
  `WARMING_UP`, don't show an error toast — return a special "not ready" marker
  and let the caller retry.
- A small startup gate component: poll `/ready` (or react to `WARMING_UP`),
  show a full-page **spinner** ("Starting Hera…") until ready, then let the
  normal fetches run.

### Tradeoffs
- The spinner replaces the current blank/half-loaded window with an explicit
  "starting" state — better UX, and it composes with issue #1011's retry.
- One new flag + endpoint; `/exec` gains a cheap guard. No change to how code
  runs once warm.

---

## Plan 2 — Unify the concurrent startup exec calls

### Today
At startup the client fires **three** separate `/exec` calls at once, from two
components:
- `fetchProjectsNames` — project names
- `fetchProjectData` — toolkit docs + project (already 2 commands in 1 call)
- `ServerConstantReader.readAllConstants` — datatypes

Three round-trips, three cold-start hits.

### Option A — Manual single loader
One startup function calls `fetchPython(namesCmd, toolkitsCmd, projectCmd,
datatypesCmd)`. `fetchPython`/`assembleCode` **already** support many commands
in one call, so this is a one-request change.

- **Pro:** one round-trip instead of three.
- **Con:** couples independent fetches. `assembleCode` concatenates all
  commands into one `exec`; if **any** line throws, the whole call fails and
  nothing loads. Startup becomes all-or-nothing.

### Option B — Auto-batch layer (debounce)
A thin queue in `fetchPython`: collect commands issued within the same tick,
flush them as one `/exec`, split the result back to each caller. Callers stay
independent; batching is transparent.

- **Pro:** fewer round-trips **and** components stay decoupled.
- **Con:** more machinery; still shares the all-or-nothing `exec` failure unless
  the server runs each command in its own try/except and returns per-command
  results.

### Does it break the architecture?
No. `/exec` is **stateless** — it just runs code and returns `result`. Batching
is purely a client-side concern (how many commands per call). The only real
risk is **error isolation**: one bundled `exec` fails as a unit. If we care
about partial success, the server's `exec_code` would need to run commands
separately and collect per-command errors — a bigger change.

### Recommendation
Start with **Option A** (trivial, big win on round-trips), but only after Plan 1
— once the server gates on readiness, a single failed bundle is far less likely
(no more cold-start races). Revisit Option B only if error isolation matters.

---

## Effort
- Plan 1: small–medium (background thread + flag + endpoint + client spinner/retry).
- Plan 2 Option A: small. Option B: medium.
```
