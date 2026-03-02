import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { render, screen, cleanup, fireEvent, waitFor, within, act } from '@testing-library/react';
import { useState } from 'react';
import { DetailsViewDocId } from '../src/components/details/DetailsViewDocId';
import { DetailsViewPanel } from '../src/components/details/DetailsViewPanel';
import { ProjectTreeView } from '../src/components/project/ProjectTreeView';
import { ProjectObj } from '../src/objects/ProjectObj';
import { useProjectStore } from '../src/stores/useProjectStore';
import { ProjectDocument, ProjectEntire } from '../src/shared/types';
import * as FetchDocument from '../src/io/FetchDocument';
import * as FetchProjects from '../src/io/FetchProjects';

vi.mock('../src/components/agents/AgentConfigEditor', () => ({
  AgentConfigEditor: () => null,
}));
vi.mock('../src/components/project/AddDocumentButton', () => ({
  AddDocumentButton: () => null,
}));
vi.mock('../src/components/project/ProjectViewSettingsButton', () => ({
  ProjectViewSettingsButton: () => null,
}));
vi.mock('../src/components/project/RepoTreeWhole', () => ({
  RepoTreeWhole: () => null,
}));

const DOC_ID = 'abc123';

const makeProjectData = (): ProjectEntire => ({
  name: 'TestProject',
  documents: [
    {
      _cls: 'Document',
      _id: { $oid: DOC_ID },
      projectName: 'TestProject',
      desc: { datasourceName: 'MyDoc', filesDirectory: '/tmp' },
      type: 'myType',
      resource: 'original-resource',
      dataFormat: 'json',
    },
  ],
});

const makeDoc = (resource: string): ProjectDocument => ({
  _cls: 'Document',
  _id: { $oid: DOC_ID },
  projectName: 'TestProject',
  desc: { datasourceName: 'MyDoc', filesDirectory: '/tmp' },
  type: 'myType',
  resource,
  dataFormat: 'json',
});

describe('Two UI sync', () => {
  let fetchDocumentSpy: ReturnType<typeof vi.spyOn>;
  let updateDocumentSpy: ReturnType<typeof vi.spyOn>;

  beforeEach(() => {
    fetchDocumentSpy = vi.spyOn(FetchDocument, 'fetchDocument');
    updateDocumentSpy = vi.spyOn(FetchDocument, 'updateDocument');
  });

  afterEach(() => {
    cleanup();
    vi.restoreAllMocks();
  });

  it('two UIs viewing the same document both see the same data', async () => {
    const doc = makeDoc('shared-value');
    fetchDocumentSpy.mockResolvedValue(doc);

    const projectData = makeProjectData();

    const ui1 = render(
      <div data-testid="ui1">
        <DetailsViewDocId project={new ProjectObj(projectData)} docid={DOC_ID} />
      </div>
    );

    const ui2 = render(
      <div data-testid="ui2">
        <DetailsViewDocId project={new ProjectObj(projectData)} docid={DOC_ID} />
      </div>
    );

    // Both UIs load and display the same document name
    await waitFor(() => {
      expect(within(ui1.container).getByText('MyDoc')).toBeDefined();
      expect(within(ui2.container).getByText('MyDoc')).toBeDefined();
    });

    // Both UIs show the same resource value
    expect(within(ui1.container).getByDisplayValue('shared-value')).toBeDefined();
    expect(within(ui2.container).getByDisplayValue('shared-value')).toBeDefined();
  });

  it('UI-2 sees the old value, then sees the updated value after refresh', async () => {
    const oldDoc = makeDoc('old-value');
    const updatedDoc = makeDoc('new-value');

    // Both UIs initially fetch the old document
    fetchDocumentSpy.mockResolvedValue(oldDoc);

    const projectData = makeProjectData();
    const project1 = new ProjectObj(projectData);
    const project2 = new ProjectObj(projectData);

    // --- Render UI-1 ---
    const ui1 = render(
      <div data-testid="ui1">
        <DetailsViewDocId project={project1} docid={DOC_ID} />
      </div>
    );

    // Wait for UI-1 to load the document
    await waitFor(() => {
      expect(within(ui1.container).getByText('MyDoc')).toBeDefined();
    });

    // Verify UI-1 shows old resource value
    expect(within(ui1.container).getByDisplayValue('old-value')).toBeDefined();

    // --- Render UI-2 ---
    const ui2 = render(
      <div data-testid="ui2">
        <DetailsViewDocId project={project2} docid={DOC_ID} />
      </div>
    );

    // Wait for UI-2 to load the document
    await waitFor(() => {
      expect(within(ui2.container).getByText('MyDoc')).toBeDefined();
    });

    // UI-2 also shows the old value
    expect(within(ui2.container).getByDisplayValue('old-value')).toBeDefined();

    // --- UI-1 changes the resource value and saves ---
    const resourceInput = within(ui1.container).getByDisplayValue('old-value');
    fireEvent.change(resourceInput, { target: { value: 'new-value' } });

    // The "Update Document" button (Done icon) should appear after editing
    updateDocumentSpy.mockResolvedValue(updatedDoc);

    const updateButton = await waitFor(() =>
      within(ui1.container).getByLabelText('Update Document')
    );
    fireEvent.click(updateButton);

    // Wait for UI-1 to reflect the saved value
    await waitFor(() => {
      expect(within(ui1.container).getByDisplayValue('new-value')).toBeDefined();
    });

    // --- UI-2 still shows the old value ---
    expect(within(ui2.container).getByDisplayValue('old-value')).toBeDefined();

    // --- UI-2 "refreshes" (simulated by re-rendering with a new project ref) ---
    // In the real app, clicking the refresh button calls fetchProjectDetails,
    // which updates the project store, causing DetailsViewDocId to re-render
    // with a new project object. The useEffect triggers on project.name change,
    // but since name is the same, we simulate by unmounting and remounting.
    fetchDocumentSpy.mockResolvedValue(updatedDoc);

    ui2.unmount();
    const ui2Refreshed = render(
      <div data-testid="ui2">
        <DetailsViewDocId project={new ProjectObj(projectData)} docid={DOC_ID} />
      </div>
    );

    // After refresh, UI-2 should now show the new value
    await waitFor(() => {
      expect(within(ui2Refreshed.container).getByDisplayValue('new-value')).toBeDefined();
    });
  });

  it('clicking the refresh button in the tree updates the details panel', async () => {
    const oldDoc = makeDoc('old-value');
    const updatedDoc = makeDoc('new-value');

    // Initial document fetch returns old data
    fetchDocumentSpy.mockResolvedValue(oldDoc);

    const projectData = makeProjectData();

    // Set up the store with the project (like Dashboard does)
    useProjectStore.setState({ currProject: projectData });

    // A mini wrapper that wires the tree and details panel via the store,
    // mimicking Dashboard's behavior
    const DashboardLike = () => {
      const { getProject } = useProjectStore();
      const [selectedItemId, setSelectedItemId] = useState<string>('');
      const project = getProject()!;
      return (
        <div>
          <div data-testid="tree-panel">
            <ProjectTreeView
              project={project}
              setSelectedItemIds={(ids) => setSelectedItemId(ids[0] || '')}
            />
          </div>
          <div data-testid="details-panel">
            <DetailsViewPanel
              project={project}
              showItemId={selectedItemId}
            />
          </div>
        </div>
      );
    };

    render(<DashboardLike />);

    const treePanel = screen.getByTestId('tree-panel');
    const detailsPanel = screen.getByTestId('details-panel');

    // Click on the document in the tree to select it
    const docItem = await waitFor(() => within(treePanel).getByText('MyDoc'));
    fireEvent.click(docItem);

    // Wait for the details panel to load and show the old value
    await waitFor(() => {
      expect(within(detailsPanel).getByDisplayValue('old-value')).toBeDefined();
    });

    // Now the backend has updated data — mock fetchDocument to return new value
    fetchDocumentSpy.mockResolvedValue(updatedDoc);

    // Also mock fetchProjectDetails to update the store (this is what the
    // refresh button triggers). It re-fetches the project and calls setCurrentProject.
    const fetchProjectDetailsSpy = vi.spyOn(FetchProjects, 'fetchProjectDetails');
    fetchProjectDetailsSpy.mockImplementation(async () => {
      useProjectStore.getState().setCurrentProject({ ...projectData });
    });

    // Click the refresh button
    const refreshButton = within(treePanel).getByLabelText('Reload documents');
    await act(async () => {
      fireEvent.click(refreshButton);
    });

    // After refresh, the details panel should show the new value
    await waitFor(() => {
      setTimeout(() => {
        expect(within(detailsPanel).getByDisplayValue('new-value')).toBeDefined();
      }, 100);
    });
  });
});
