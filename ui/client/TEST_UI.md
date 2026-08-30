# UI Validation Checklist

Run these steps **in order** from the repo root to validate the UI client.

## 1. TypeScript type-checking

```bash
cd ui/client && npx tsc --noEmit
```

## 2. Unit tests

```bash
cd ui/client && npm run test
```

## 3. Production build

```bash
cd ui/client && npx vite build --outDir /tmp/hera-vite-build --emptyOutDir
```

## 4. Server tests

```bash
pytest ui/server/tests
```

Needs only `fastapi`, `pydantic`, `argcomplete` and `pytest`
(`pip install -r ui/server/requirements.txt`) — the tests install a fake `hera`
module, so no database is involved. In Docker: `bash ui/server/scripts/run_tests_docker.sh`.

## 5. Integration tests

```bash
cd ui/client && npm run test:integ
```

Builds the test image and runs the suite against a real server and MongoDB inside
the container. To run them directly instead, start MongoDB on port 27018 with the
`hera`/`heracles` user and run:

```bash
cd ui/client && npx vitest run tests/integration -c vitest.integ.config.ts --no-file-parallelism --retry 1
```

`tests/integration/globalSetup.ts` starts `ui/server/server.py` itself on a free
port and waits for `/ready`.

---

CI runs all five steps on every pull request — see the `ui-client`, `ui-server`
and `ui-integration` jobs in `.github/workflows/ci.yml`.
