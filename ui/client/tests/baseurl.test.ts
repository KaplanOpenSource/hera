/// <reference types="vite/client" />
import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';

const originalLocation = window.location;
const originalImportMeta = import.meta.env;

function mockLocation(hostname: string, port: string) {
  Object.defineProperty(window, 'location', {
    value: { ...originalLocation, hostname, port },
    writable: true,
    configurable: true,
  });
}

function restoreLocation() {
  Object.defineProperty(window, 'location', {
    value: originalLocation,
    writable: true,
    configurable: true,
  });
}

async function loadBaseurl() {
  vi.resetModules();
  const mod = await import('../src/shared/baseurl');
  return mod.BASEURL;
}

describe('BASEURL', () => {
  afterEach(() => {
    restoreLocation();
    import.meta.env.DEV = originalImportMeta.DEV;
    delete import.meta.env.VITE_API_PORT;
  });

  it('returns empty string when served from localhost in production', async () => {
    mockLocation('localhost', '8000');
    import.meta.env.DEV = false;
    const baseurl = await loadBaseurl();
    expect(baseurl).toBe('');
  });

  it('returns empty string for any port on localhost in production', async () => {
    mockLocation('localhost', '9000');
    import.meta.env.DEV = false;
    const baseurl = await loadBaseurl();
    expect(baseurl).toBe('');
  });

  it('returns full URL in dev mode even on localhost', async () => {
    mockLocation('localhost', '5173');
    import.meta.env.DEV = true;
    const baseurl = await loadBaseurl();
    expect(baseurl).toBe('http://localhost:8000');
  });

  it('uses VITE_API_PORT in dev mode', async () => {
    mockLocation('localhost', '5173');
    import.meta.env.DEV = true;
    import.meta.env.VITE_API_PORT = '9000';
    const baseurl = await loadBaseurl();
    expect(baseurl).toBe('http://localhost:9000');
  });

  it('returns full URL for non-localhost hostname', async () => {
    mockLocation('10.0.0.5', '8000');
    import.meta.env.DEV = false;
    const baseurl = await loadBaseurl();
    expect(baseurl).toBe('http://localhost:8000');
  });

  it('returns empty string for 127.0.0.1 in production', async () => {
    mockLocation('127.0.0.1', '8000');
    import.meta.env.DEV = false;
    const baseurl = await loadBaseurl();
    expect(baseurl).toBe('');
  });
});
