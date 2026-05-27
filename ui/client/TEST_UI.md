# UI Validation Checklist

Run these steps **in order** from the repo root (`/home/eran/Code/hera`) to validate the UI client.

## 1. TypeScript type-checking

```bash
cd /home/eran/Code/hera/ui/client && npx tsc --noEmit
```

## 2. Unit tests

```bash
cd /home/eran/Code/hera/ui/client && npm run test
```

## 3. Production build

```bash
cd /home/eran/Code/hera/ui/client && npm run build
```

## 4. Clean build artifacts

**CRITICAL: Always run from repo root.** The build creates new hash-named files and modifies commitId.ts. Both must be reverted.

```bash
cd /home/eran/Code/hera && rm -f ui/client/bundle/assets/index-*.js ui/client/bundle/assets/index-*.css && git checkout -- ui/client/bundle/ ui/client/src/commitId.ts
```

## 5. Verify clean state

```bash
cd /home/eran/Code/hera && git status ui/client/bundle/ ui/client/src/commitId.ts
```

Must show "nothing to commit, working tree clean" with no untracked files.
