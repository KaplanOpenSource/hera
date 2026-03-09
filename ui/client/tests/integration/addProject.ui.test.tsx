/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter, Routes, Route, useParams } from 'react-router-dom';
import { startDockerEnv, type DockerEnv } from './dockerSetup';

const TEST_SERVER_PORT = 8002;

// Point BASEURL to the test server
vi.mock('../../src/shared/baseurl', () => ({
  BASEURL: `http://localhost:${TEST_SERVER_PORT}`,
}));

// Mock components that aren't needed for this test
vi.mock('../../src/components/header/DeleteProjectButton', () => ({
  DeleteProjectButton: () => null,
}));

let env: DockerEnv;

// Lazy-import components after mocks are set up
const loadComponents = async () => {
  const { ProjectChooser } = await import('../../src/components/header/ProjectChooser');
  const { FetchProjects } = await import('../../src/io/FetchProjects');
  const { useProjectStore } = await import('../../src/stores/useProjectStore');
  return { ProjectChooser, FetchProjects, useProjectStore };
};

describe('Add Project UI integration', () => {
  let ProjectChooser: Awaited<ReturnType<typeof loadComponents>>['ProjectChooser'];
  let FetchProjects: Awaited<ReturnType<typeof loadComponents>>['FetchProjects'];
  let useProjectStore: Awaited<ReturnType<typeof loadComponents>>['useProjectStore'];

  beforeAll(async () => {
    env = await startDockerEnv({
      network: 'hera-test-ui-net',
      mongoContainer: 'hera-test-ui-mongo',
      serverContainer: 'hera-test-ui-server',
      serverPort: TEST_SERVER_PORT,
      dbName: 'hera_test_ui',
    });

    const components = await loadComponents();
    ProjectChooser = components.ProjectChooser;
    FetchProjects = components.FetchProjects;
    useProjectStore = components.useProjectStore;
  }, 30000);

  beforeEach(() => {
    useProjectStore.setState({
      projectNames: [],
      currProjectName: '* NONE *',
      currProject: null,
      toolkits: [],
    });
  });

  afterEach(() => {
    cleanup();
  });

  afterAll(() => {
    cleanup();
    env.cleanup();
  }, 15000);

  it('clicking Add Project creates the project and it appears in the list', async () => {
    const Page = () => {
      const { projectName } = useParams<{ projectName: string }>();
      return (
        <>
          <FetchProjects urlProjectName={projectName} />
          <ProjectChooser />
        </>
      );
    };
    const App = () => (
      <MemoryRouter initialEntries={['/']}>
        <Routes>
          <Route path="/:projectName" element={<Page />} />
          <Route path="/" element={<Page />} />
        </Routes>
      </MemoryRouter>
    );

    render(<App />);

    // Wait for the add button to be available
    const addWrapper = await screen.findByLabelText('Add project');
    const addButton = within(addWrapper).getByRole('button');

    // Click the "+" button to open the dialog
    await act(async () => {
      fireEvent.click(addButton);
    });

    // Wait for the Dialog to appear and find the Project Name input
    const dialog = await screen.findByRole('dialog', {}, { timeout: 3000 });
    const nameInput = within(dialog).getByRole('textbox', { name: /project name/i });
    fireEvent.change(nameInput, { target: { value: 'UITestProject' } });

    // Uncheck "Load Repositories" to speed things up
    const loadReposCheckbox = within(dialog).getByRole('checkbox', { name: /load repositories/i });
    if ((loadReposCheckbox as HTMLInputElement).checked) {
      fireEvent.click(loadReposCheckbox);
    }

    // Click "Add Project"
    await act(async () => {
      fireEvent.click(within(dialog).getByText('Add Project'));
    });

    // Wait for the project to appear in the store's project list
    await waitFor(() => {
      const names = useProjectStore.getState().projectNames;
      expect(names.some(p => p.name === 'UITestProject')).toBe(true);
    }, { timeout: 15000 });
  }, 30000);

  it('project persists after simulated refresh (store reset + refetch)', async () => {
    // Reset the store to simulate a fresh page load
    useProjectStore.setState({
      projectNames: [],
      currProjectName: '* NONE *',
      currProject: null,
      toolkits: [],
    });

    const Page = () => {
      const { projectName } = useParams<{ projectName: string }>();
      return (
        <>
          <FetchProjects urlProjectName={projectName} />
          <ProjectChooser />
        </>
      );
    };
    const App = () => (
      <MemoryRouter initialEntries={['/']}>
        <Routes>
          <Route path="/:projectName" element={<Page />} />
          <Route path="/" element={<Page />} />
        </Routes>
      </MemoryRouter>
    );

    render(<App />);

    // FetchProjects should refetch project names on mount.
    // The project created in the previous test should still be there.
    await waitFor(() => {
      const names = useProjectStore.getState().projectNames;
      expect(names.some(p => p.name === 'UITestProject')).toBe(true);
    }, { timeout: 10000 });
  }, 15000);
}, 90000);
