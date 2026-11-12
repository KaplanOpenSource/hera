/// <reference types="vitest/config" />
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { fileURLToPath } from 'url';
import { dirname } from 'path';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

export default defineConfig({
  root: __dirname,
  plugins: [react()],
  server: {
    port: 5173,
    strictPort: true
  },
  build: {
    outDir: 'bundle',
    chunkSizeWarningLimit: 2000,
    rollupOptions: {
      // Keep lazy imports (e.g. geotiff decoders) in the one bundle file.
      output: { inlineDynamicImports: true },
    },
  },
  test: {
    environment: 'jsdom',
    setupFiles: ['./tests/setup.ts'],
    // Cap worker processes: the default is one per CPU (20 here), each a full
    // Node+jsdom fork, which exhausts memory on a loaded machine and hangs it.
    // Two still runs the suite in well under a minute and leaves the machine
    // usable while it runs.
    maxWorkers: 2,
    minWorkers: 1,
  }
});


