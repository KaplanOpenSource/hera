import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, waitFor, act } from '@testing-library/react';
import { BrowserRouter, Route, Routes, useParams } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { mockApiCall } from './helpers/mockApiCall';

let mockGetJson: ReturnType<typeof vi.fn>;

const { FetchProjects } = await import('../src/io/FetchProjects');
const { useProjectStore, NO_PROJECT } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

const CHOSEN = { name: 'AgentTestProject', documents: [{ _id: { $oid: 'a' } }] };
const WANTED = { name: 'DocTestProject', documents: [{ _id: { $oid: 'd' } }] };
const NAMES = [{ name: 'AgentTestProject' }, { name: 'DocTestProject' }];

const Page = () => {
  const { projectName } = useParams<{ projectName: string }>();
  return <FetchProjects urlProjectName={projectName} />;
};

const renderAt = (path: string) => {
  window.history.pushState({}, '', path);
  render(
    <BrowserRouter>
      <Routes>
        <Route path="/:projectName/:docId" element={<Page />} />
        <Route path="/:projectName" element={<Page />} />
        <Route path="/" element={<Page />} />
      </Routes>
    </BrowserRouter>
  );
};

// What fetchInit does when the server picked a project the URL did not ask for.
const startupChose = (project: typeof CHOSEN) => {
  useProjectListStore.getState().setProjectNames(NAMES);
  useProjectStore.getState().selectProject(project.name);
  useProjectStore.getState().setCurrentProject(project);
};

beforeEach(() => {
  vi.clearAllMocks();
  const call = mockApiCall();
  mockGetJson = call.getJson;
  useProjectListStore.setState({ projectNames: [] });
  useProjectStore.setState({ currProjectName: NO_PROJECT, currProject: null, editedDocs: {} });
  mockGetJson.mockImplementation(async (path: string, params: any) => {
    if (path === '/project/toolkits') {
      return { toolkitDocs: [] };
    }
    if (path === '/project/names') {
      return { projects: NAMES };
    }
    if (params?.projectName === 'DocTestProject') {
      return { project: WANTED };
    }
    return { project: CHOSEN };
  });
});

afterEach(() => {
  cleanup();
});

describe('the URL decides which project is shown', () => {
  it('switches to the URL project when startup chose another one', async () => {
    renderAt('/DocTestProject');
    await act(async () => { startupChose(CHOSEN); });

    await waitFor(() => {
      expect(useProjectStore.getState().currProjectName).toBe('DocTestProject');
      expect(useProjectStore.getState().currProject?.name).toBe('DocTestProject');
    });
    expect(window.location.pathname).toBe('/DocTestProject');
  });

  it('switches when the URL also names a document', async () => {
    renderAt('/DocTestProject/d');
    await act(async () => { startupChose(CHOSEN); });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('DocTestProject');
    });
    // The document part of the URL must survive, or the opened document is lost.
    expect(window.location.pathname).toBe('/DocTestProject/d');
  });

  it('keeps the startup project when the URL names none', async () => {
    renderAt('/');
    await act(async () => { startupChose(CHOSEN); });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('AgentTestProject');
      expect(window.location.pathname).toBe('/AgentTestProject');
    });
  });
});
