import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, cleanup, waitFor } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

const mockFetchInit = vi.fn();
const mockNavigate = vi.fn();
vi.mock('../src/io/fetchInit', () => ({
  fetchInit: (...args: any[]) => mockFetchInit(...args),
}));
vi.mock('react-router-dom', async () => {
  const actual = await vi.importActual<any>('react-router-dom');
  return { ...actual, useNavigate: () => mockNavigate };
});

const { ServerReadyGate } = await import('../src/components/ServerReadyGate');

const renderGate = () => {
  return render(
    <MemoryRouter>
      <ServerReadyGate><div>app mounted</div></ServerReadyGate>
    </MemoryRouter>
  );
};

beforeEach(() => {
  vi.clearAllMocks();
  window.history.replaceState({}, '', '/');
});

afterEach(() => {
  cleanup();
});

describe('ServerReadyGate', () => {
  it('shows the spinner while the server is warming', async () => {
    mockFetchInit.mockResolvedValue({ ready: false });

    renderGate();

    expect(await screen.findByText('Starting Hera…')).toBeTruthy();
    expect(screen.queryByText('app mounted')).toBeNull();
  });

  it('mounts the app once the server is ready', async () => {
    mockFetchInit.mockResolvedValue({ ready: true, project: { name: 'Alpha', documents: [] } });
    window.history.replaceState({}, '', '/Alpha');

    renderGate();

    expect(await screen.findByText('app mounted')).toBeTruthy();
    expect(mockNavigate).not.toHaveBeenCalled();
  });

  it('passes the URL project name, decoded', async () => {
    mockFetchInit.mockResolvedValue({ ready: true, project: { name: 'My Project', documents: [] } });
    window.history.replaceState({}, '', '/My%20Project');

    renderGate();

    await screen.findByText('app mounted');
    expect(mockFetchInit).toHaveBeenCalledWith('My Project');
  });

  it('corrects the URL when another project was loaded', async () => {
    mockFetchInit.mockResolvedValue({ ready: true, project: { name: 'Alpha', documents: [] } });
    window.history.replaceState({}, '', '/NoSuchProject');

    renderGate();

    await screen.findByText('app mounted');
    expect(mockNavigate).toHaveBeenCalledWith('/Alpha', { replace: true });
  });

  it('keeps polling when the server does not answer', async () => {
    mockFetchInit.mockRejectedValueOnce(new Error('down'));
    mockFetchInit.mockResolvedValue({ ready: true, project: { name: 'Alpha', documents: [] } });

    renderGate();

    // The retry waits out the poll interval, so allow more than the default.
    expect(await screen.findByText('app mounted', {}, { timeout: 3000 })).toBeTruthy();
    expect(mockFetchInit.mock.calls.length).toBeGreaterThan(1);
  });
});
