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
