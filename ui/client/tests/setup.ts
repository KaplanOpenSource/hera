import { vi } from 'vitest';

global.ResizeObserver = class ResizeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
};

// Suppress verbose logging noise during tests
vi.spyOn(console, 'log').mockImplementation(() => {});
vi.spyOn(console, 'trace').mockImplementation(() => {});
