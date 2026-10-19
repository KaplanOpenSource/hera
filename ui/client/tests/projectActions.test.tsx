import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

const { DeleteSelectedButton } = await import('../src/components/project/DeleteSelectedButton');
const { ProjectActionsButton } = await import('../src/components/project/ProjectActionsButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');
const { idDocId, idRepoId } = await import('../src/shared/idDocId');

const configDoc = {
  _id: { $oid: 'cfg1' },
  desc: { datasourceName: 'config', filesDirectory: '/tmp/testproject' },
  type: 'TestProject__config__',
  _cls: 'Cache',
  resource: '',
  dataFormat: 'string',
  projectName: 'TestProject',
};
const doc = (oid: string) => ({
  _id: { $oid: oid },
  desc: { datasourceName: oid },
  type: 'T',
  _cls: 'Measurements',
  resource: '',
  dataFormat: 'string',
  projectName: 'TestProject',
});

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents: [configDoc, doc('doc1'), doc('doc2')] },
  });
  useProjectListStore.setState({ projectNames: [{ name: 'TestProject' }] });
});

afterEach(() => {
  cleanup();
});

const confirmDelete = async () => {
  const dialog = await screen.findByRole('dialog');
  await act(async () => {
    fireEvent.click(within(dialog).getByRole('button', { name: /^yes$/i }));
  });
};

// The Actions trigger is an icon button (no text label), so open it via its icon.
const openActions = async () => {
  await act(async () => {
    fireEvent.click(screen.getByTestId('MoreVertIcon').closest('button')!);
  });
};

describe('DeleteSelectedButton', () => {
  it('is disabled when no documents are selected', () => {
    render(<DeleteSelectedButton selectedIds={[idRepoId('repoA'), 'split_field=value']} />);
    const btn = screen.getByRole('button', { name: /delete selected/i }) as HTMLButtonElement;
    expect(btn.disabled).toBe(true);
  });

  it('shows the count of selected documents', () => {
    render(<DeleteSelectedButton selectedIds={[idDocId('doc1'), idDocId('doc2')]} />);
    expect(screen.getByRole('button', { name: /delete selected \(2\)/i })).toBeTruthy();
  });

  it('deletes only the selected documents, skipping repos/splits/config', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [configDoc] } });
    const onDeleted = vi.fn();
    render(
      <DeleteSelectedButton
        selectedIds={[idDocId('doc1'), idDocId('doc2'), idRepoId('repoA'), 'split_field=value', idDocId('cfg1')]}
        onDeleted={onDeleted}
      />
    );

    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /delete selected/i }));
    });
    await confirmDelete();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostJson.mock.calls[0][0]).toBe('/document/delete');
    expect(mockPostJson.mock.calls[0][1]).toEqual({ docIds: ['doc1', 'doc2'] });

    // The answer carries the project, so nothing is reloaded after it.
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);
    expect(onDeleted).toHaveBeenCalledTimes(1);
  });

  it('does nothing when the confirmation is cancelled', async () => {
    const onDeleted = vi.fn();
    render(<DeleteSelectedButton selectedIds={[idDocId('doc1')]} onDeleted={onDeleted} />);

    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /delete selected/i }));
    });
    const dialog = await screen.findByRole('dialog');
    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /^no$/i }));
    });

    expect(mockPostJson).not.toHaveBeenCalled();
    expect(onDeleted).not.toHaveBeenCalled();
  });

  it('does not refresh or notify when the call fails', async () => {
    mockPostJson.mockResolvedValueOnce(null);
    const onDeleted = vi.fn();
    render(<DeleteSelectedButton selectedIds={[idDocId('doc1')]} onDeleted={onDeleted} />);

    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /delete selected/i }));
    });
    await confirmDelete();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(onDeleted).not.toHaveBeenCalled();
  });
});

describe('ProjectActionsButton', () => {
  it('opens a popover with the Actions title and both actions', async () => {
    render(<ProjectActionsButton selectedIds={[idDocId('doc1')]} />);

    await openActions();

    expect(screen.getByText('Actions')).toBeTruthy();
    expect(screen.getByRole('button', { name: /add document/i })).toBeTruthy();
    expect(screen.getByRole('button', { name: /delete selected/i })).toBeTruthy();
  });

  it('clears the selection after a successful bulk delete', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [configDoc] } });
    const onSelectDocument = vi.fn();
    render(<ProjectActionsButton selectedIds={[idDocId('doc1')]} onSelectDocument={onSelectDocument} />);

    await openActions();
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /delete selected/i }));
    });
    await confirmDelete();

    await waitFor(() => {
      expect(onSelectDocument).toHaveBeenCalledWith(undefined);
    });
  });
});
