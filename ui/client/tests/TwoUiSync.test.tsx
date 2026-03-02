import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { render, screen, cleanup, fireEvent, waitFor, within } from '@testing-library/react';
import { DetailsViewDocId } from '../src/components/details/DetailsViewDocId';
import { ProjectObj } from '../src/objects/ProjectObj';
import { ProjectDocument, ProjectEntire } from '../src/shared/types';
import * as FetchDocument from '../src/io/FetchDocument';

vi.mock('../src/components/agents/AgentConfigEditor', () => ({
  AgentConfigEditor: () => null,
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
});
