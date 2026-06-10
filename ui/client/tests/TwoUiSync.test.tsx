import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { render, screen, cleanup, fireEvent, waitFor, within, act } from '@testing-library/react';
import { useState } from 'react';
import { MemoryRouter } from 'react-router-dom';
import { DetailsViewDocId } from '../src/components/details/DetailsViewDocId';
import { DetailsViewPanel } from '../src/components/details/DetailsViewPanel';
import { ProjectTreeView } from '../src/components/project/ProjectTreeView';
import { ProjectObj } from '../src/objects/ProjectObj';
import { useProjectStore } from '../src/stores/useProjectStore';
import { useServerConstants } from '../src/stores/useServerConstants';
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
vi.mock('../src/io/snackbar', () => ({
  pushRunning: () => 'mock-key',
  pushError: () => 'mock-key',
  dismiss: () => {},
}));

const DOC_ID = 'abc123';

// The document data now lives in the project (loaded centrally into the store),
// so tests build a project carrying the resource value they expect to see.
const makeProjectData = (resource: string): ProjectEntire => ({
  name: 'TestProject',
  documents: [
    {
      _cls: 'Document',
      _id: { $oid: DOC_ID },
      projectName: 'TestProject',
      desc: { datasourceName: 'MyDoc', filesDirectory: '/tmp' },
      type: 'myType',
      resource,
      dataFormat: 'json',
    },
  ],
});

const makeDoc = (resource: string): ProjectDocument => makeProjectData(resource).documents[0];

describe('Two UI sync', () => {
  let updateDocumentSpy: ReturnType<typeof vi.spyOn>;

  beforeEach(() => {
    updateDocumentSpy = vi.spyOn(FetchDocument, 'updateDocument');
    useServerConstants.setState({ dataTypes: { JSON: 'json' } });
  });

  afterEach(() => {
    cleanup();
    vi.restoreAllMocks();
  });

  it('two UIs viewing the same document both see the same data', async () => {
    const projectData = makeProjectData('shared-value');

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

  it('a second UI sees the new value after its project updates', async () => {
    const ui1 = render(
      <div data-testid="ui1">
        <DetailsViewDocId project={new ProjectObj(makeProjectData('old-value'))} docid={DOC_ID} />
      </div>
    );
    await waitFor(() => {
      expect(within(ui1.container).getByDisplayValue('old-value')).toBeDefined();
    });

    const ui2 = render(
      <div data-testid="ui2">
        <DetailsViewDocId project={new ProjectObj(makeProjectData('old-value'))} docid={DOC_ID} />
      </div>
    );
    await waitFor(() => {
      expect(within(ui2.container).getByDisplayValue('old-value')).toBeDefined();
    });

    // --- UI-1 changes the resource value and saves ---
    const resourceInput = within(ui1.container).getByDisplayValue('old-value');
    fireEvent.change(resourceInput, { target: { value: 'new-value' } });

    updateDocumentSpy.mockResolvedValue(makeDoc('new-value'));
    // The save pulls the project back into the store; no-op it here.
    vi.spyOn(FetchProjects, 'fetchProjectDetails').mockResolvedValue(undefined);

    const updateWrapper = await waitFor(() =>
      within(ui1.container).getByLabelText('Update Document')
    );
    fireEvent.click(within(updateWrapper).getByRole('button'));
    await waitFor(() => expect(updateDocumentSpy).toHaveBeenCalled());

    // UI-2 still shows the old value (its project hasn't been reloaded yet)
    expect(within(ui2.container).getByDisplayValue('old-value')).toBeDefined();

    // UI-2's project updates (as the auto-reload would do) → it shows the new value
    ui2.rerender(
      <div data-testid="ui2">
        <DetailsViewDocId project={new ProjectObj(makeProjectData('new-value'))} docid={DOC_ID} />
      </div>
    );
    await waitFor(() => {
      expect(within(ui2.container).getByDisplayValue('new-value')).toBeDefined();
    });
  });

  it('clicking the refresh button in the tree updates the details panel', async () => {
    // Set up the store with the project (like Dashboard does)
    useProjectStore.setState({ currProject: makeProjectData('old-value') });

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
              onSelectItem={(id) => setSelectedItemId(id ?? '')}
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

    render(<MemoryRouter><DashboardLike /></MemoryRouter>);

    const treePanel = screen.getByTestId('tree-panel');
    const detailsPanel = screen.getByTestId('details-panel');

    // Click on the document in the tree to select it
    const docItem = await waitFor(() => within(treePanel).getByText('MyDoc'));
    fireEvent.click(docItem);

    // Wait for the details panel to show the old value (read from the store)
    await waitFor(() => {
      expect(within(detailsPanel).getByDisplayValue('old-value')).toBeDefined();
    });

    // The refresh button calls fetchProjectDetails; simulate it updating the store
    // with the new backend value.
    const fetchProjectDetailsSpy = vi.spyOn(FetchProjects, 'fetchProjectDetails');
    fetchProjectDetailsSpy.mockImplementation(async () => {
      useProjectStore.getState().setCurrentProject(makeProjectData('new-value'));
    });

    // Click the refresh button (the actual <button> inside the labelled wrapper)
    const refreshWrapper = within(treePanel).getByLabelText('Reload documents');
    const refreshButton = within(refreshWrapper).getByRole('button');
    await act(async () => {
      fireEvent.click(refreshButton);
    });

    // After refresh, the details panel should show the new value
    await waitFor(() => {
      expect(within(detailsPanel).getByDisplayValue('new-value')).toBeDefined();
    });
  });
});
