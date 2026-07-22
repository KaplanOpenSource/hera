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
cd ui/client && npm run build
```

## 4. Clean build artifacts

**CRITICAL: Always run from repo root.** The build creates new hash-named files and modifies buildNumber.ts. Both must be reverted.

```bash
rm -f ui/client/bundle/assets/index-*.js ui/client/bundle/assets/index-*.css && git checkout -- ui/client/bundle/ ui/client/src/buildNumber.ts
```

## 5. Verify clean state

```bash
git status ui/client/bundle/ ui/client/src/buildNumber.ts
```

Must show "nothing to commit, working tree clean" with no untracked files.
