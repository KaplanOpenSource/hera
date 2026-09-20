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
  },
  test: {
    environment: 'jsdom',
    setupFiles: ['./tests/setup.ts'],
    // Cap worker processes: the default is one per CPU (20 here), each a full
    // Node+jsdom fork, which exhausts memory on a loaded machine and hangs it.
    maxWorkers: 4,
    minWorkers: 1,
  }
});


