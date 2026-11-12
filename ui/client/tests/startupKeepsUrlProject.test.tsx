import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, waitFor } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));
vi.mock('../src/io/snackbar', () => ({
  pushRunning: () => 'k',
  pushError: () => 'k',
  pushInfo: () => 'k',
  dismiss: () => {},
}));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const App = (await import('../src/App')).default;
const { useProjectStore, NO_PROJECT } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

const CHOSEN = { name: 'AgentTestProject', documents: [] };
const WANTED = { name: 'DocTestProject', documents: [] };
const NAMES = [{ name: 'AgentTestProject' }, { name: 'DocTestProject' }];

const json = (data: any) => {
  return {
    ok: true,
    status: 200,
    headers: { get: () => 'application/json' },
    json: () => Promise.resolve(data),
  };
};

// The server answers /init for the project it was asked for, and picks the first
// project when asked for none - which is what happens when the gate cannot see
// the URL.
const answerAsServer = () => {
  mockFetch.mockImplementation(async (url: string) => {
    const asked = new URL(url, 'http://test').searchParams.get('project');
    if (url.includes('/init')) {
      let project = CHOSEN;
      if (asked === 'DocTestProject') {
        project = WANTED;
      }
      return json({ ready: true, projects: NAMES, datatypes: {}, toolkitDocs: [], project, server: {} });
    }
    if (url.includes('/project/toolkits')) {
      return json({ toolkitDocs: [] });
    }
    if (url.includes('/project/names')) {
      return json({ projects: NAMES });
    }
    if (url.includes('projectName=DocTestProject')) {
      return json({ project: WANTED });
    }
    return json({ project: CHOSEN });
  });
};

beforeEach(() => {
  vi.clearAllMocks();
  answerAsServer();
  window.history.pushState({}, '', '/');
  useProjectListStore.setState({ projectNames: [] });
  useProjectStore.setState({ currProjectName: NO_PROJECT, currProject: null, editedDocs: {} });
});

afterEach(() => {
  cleanup();
});

describe('startup with a project in the router URL', () => {
  it('keeps the project the router URL names, not the one the server picked', async () => {
    render(
      <MemoryRouter initialEntries={['/DocTestProject']}>
        <App />
      </MemoryRouter>
    );

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('DocTestProject');
    }, { timeout: 5000 });
  });

  it('asks the server for the project the router URL names', async () => {
    render(
      <MemoryRouter initialEntries={['/DocTestProject']}>
        <App />
      </MemoryRouter>
    );

    await waitFor(() => {
      const initCalls = mockFetch.mock.calls.map(c => String(c[0])).filter(u => u.includes('/init'));
      expect(initCalls.some(u => u.includes('project=DocTestProject'))).toBe(true);
    }, { timeout: 5000 });
  });
});
