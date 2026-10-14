import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { BrowserRouter, Route, Routes, useParams } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

// The dialog also lists the registered repositories, which still uses /exec.
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: vi.fn().mockResolvedValue({ data: null }),
  fetchPythonClean: vi.fn().mockResolvedValue({ data: null }),
  WARMING_UP: 'WARMING_UP',
}));

const mockPostJson = vi.fn();
const mockGetJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: (...args: any[]) => mockGetJson(...args),
}));

const { AddProjectButton } = await import('../src/components/header/AddProjectButton');
const { FetchProjects } = await import('../src/io/FetchProjects');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

const NEW_PROJECT = { name: 'NewOne', documents: [{ _id: { $oid: 'cfg' } }] };
const OLD_PROJECT = { name: 'Alpha', documents: [{ _id: { $oid: 'old' } }] };

// Stands in for Dashboard: it reads the project name from the URL, and React
// Router remounts it when the matched route changes (/ -> /:projectName), which
// is what hid a freshly created project.
const Page = () => {
  const { projectName } = useParams<{ projectName: string }>();
  return (<>
    <AddProjectButton />
    <FetchProjects urlProjectName={projectName} />
  </>);
};

// A real history, because FetchProjects reads window.location to fix the URL.
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

beforeEach(() => {
  vi.clearAllMocks();
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }] });
  useProjectStore.setState({ currProjectName: 'Alpha', currProject: OLD_PROJECT, editedDocs: {} });
  mockPostJson.mockResolvedValue({
    projects: [{ name: 'Alpha' }, { name: 'NewOne' }],
    project: NEW_PROJECT,
    addedNotebooks: 0,
  });
  mockGetJson.mockImplementation(async (path: string, params: any) => {
    if (path === '/project/toolkits') {
      return { toolkitDocs: [] };
    }
    if (path === '/project/names') {
      return { projects: [{ name: 'Alpha' }, { name: 'NewOne' }] };
    }
    if (params?.projectName === 'NewOne') {
      return { project: NEW_PROJECT };
    }
    return { project: OLD_PROJECT };
  });
});

afterEach(() => {
  cleanup();
});

const createProject = async () => {
  const addWrapper = await screen.findByLabelText('Add project');
  await act(async () => { fireEvent.click(within(addWrapper).getByRole('button')); });
  const dialog = await screen.findByRole('dialog');
  fireEvent.change(within(dialog).getByRole('textbox', { name: /project name/i }), {
    target: { value: 'NewOne' },
  });
  await act(async () => {
    fireEvent.click(within(dialog).getByRole('button', { name: /add project/i }));
  });
};

describe('creating a project then showing it', () => {
  it('shows the new project when the page remounts, coming from a document URL', async () => {
    renderAt('/Alpha/old');

    await createProject();

    await waitFor(() => {
      expect(useProjectListStore.getState().projectNames.map(p => p.name)).toContain('NewOne');
      expect(useProjectStore.getState().currProjectName).toBe('NewOne');
      expect(useProjectStore.getState().currProject?.name).toBe('NewOne');
    });
  });

  it('shows the new project when there was no project in the URL', async () => {
    useProjectStore.setState({ currProjectName: 'Alpha', currProject: OLD_PROJECT, editedDocs: {} });
    renderAt('/');

    await createProject();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('NewOne');
    });
  });

  it('costs one call: the project is shown without being loaded again', async () => {
    renderAt('/Alpha');

    await createProject();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('NewOne');
    });
    expect(mockPostJson).toHaveBeenCalledTimes(1);
    const asked = mockGetJson.mock.calls.map(c => `${c[0]} ${c[1]?.projectName ?? ''}`.trim());
    expect(asked).not.toContain('/project/documents NewOne');
  });

  it('shows the new project when another project was already open', async () => {
    renderAt('/Alpha');

    await createProject();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('NewOne');
    });
  });
});
