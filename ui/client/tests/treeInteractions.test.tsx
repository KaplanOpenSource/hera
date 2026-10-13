import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const mockFetchProjectDetails = vi.fn();
vi.mock('../src/io/FetchProjects', () => ({
  fetchProjectDetails: (...args: any[]) => mockFetchProjectDetails(...args),
}));

vi.mock('../src/components/project/RepoTreeWhole', () => ({
  RepoTreeWhole: () => null,
}));

const { ProjectTreeView } = await import('../src/components/project/ProjectTreeView');
const { ProjectObj } = await import('../src/objects/ProjectObj');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { idDocId } = await import('../src/shared/idDocId');
const { SplitTree } = await import('../src/utils/splitTree');
const { collectSubtreeKeys, findNodeByKey } = await import('../src/utils/collectSubtreeKeys');
const { useViewSettingsStore } = await import('../src/stores/useViewSettingsStore');

const configDoc = {
  _cls: 'Cache',
  _id: { $oid: 'cfg1' },
  projectName: 'TestProject',
  desc: { datasourceName: 'config', filesDirectory: '/tmp/test' },
  type: 'TestProject__config__',
  resource: '',
  dataFormat: 'string',
};

const doc = (oid: string, toolkit: string) => ({
  _cls: 'Measurements',
  _id: { $oid: oid },
  projectName: 'TestProject',
  desc: { datasourceName: oid, toolkit },
  type: 'T',
  resource: '',
  dataFormat: 'string',
});

const documents = [
  configDoc,
  doc('doc1', 'toolkitA'),
  doc('doc2', 'toolkitA'),
  doc('doc3', 'toolkitB'),
  doc('doc4', 'toolkitB'),
];

// Opens a collapsed branch by clicking its chevron.
const expand = (name: string) => {
  const chevron = row(name).querySelector('[class*="iconContainer"]') as HTMLElement;
  fireEvent.click(chevron);
};

const renderTree = (onSelectItem = vi.fn()) => {
  const project = new ProjectObj({ name: 'TestProject', documents } as any);
  const result = render(
    <MemoryRouter>
      <ProjectTreeView project={project} onSelectItem={onSelectItem} />
    </MemoryRouter>
  );
  expand('toolkitA');
  expand('toolkitB');
  return { ...result, onSelectItem };
};

const row = (name: string) => {
  return screen.getByText(name).closest('[role="treeitem"]') as HTMLElement;
};

const isSelected = (name: string) => {
  return row(name).getAttribute('aria-selected') === 'true';
};

const click = (name: string) => {
  fireEvent.click(screen.getByText(name));
};

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents } as any,
  });
});

afterEach(() => {
  cleanup();
});

describe('tree selection', () => {
  it('a single click selects the row without opening it', () => {
    const { onSelectItem } = renderTree();
    onSelectItem.mockClear();

    click('doc1');

    expect(isSelected('doc1')).toBe(true);
    expect(onSelectItem).not.toHaveBeenCalled();
  });

  it('a double click opens the clicked document', () => {
    const { onSelectItem } = renderTree();
    onSelectItem.mockClear();

    click('doc1');
    fireEvent.doubleClick(screen.getByText('doc1'));

    expect(onSelectItem).toHaveBeenCalledWith(idDocId('doc1'));
  });

  it('a double click on a row that is already selected still opens it', () => {
    const { onSelectItem } = renderTree();
    click('doc1');
    onSelectItem.mockClear();

    click('doc1');
    fireEvent.doubleClick(screen.getByText('doc1'));

    expect(onSelectItem).toHaveBeenCalledWith(idDocId('doc1'));
  });

  it('clicking a branch selects every document under it', () => {
    renderTree();

    click('toolkitA');

    expect(isSelected('toolkitA')).toBe(true);
    expect(isSelected('doc1')).toBe(true);
    expect(isSelected('doc2')).toBe(true);
    expect(isSelected('doc3')).toBe(false);
  });

  it('clicking the panel background clears the selection', () => {
    const { container } = renderTree();
    click('doc1');
    expect(isSelected('doc1')).toBe(true);

    fireEvent.click(container.firstChild as HTMLElement);

    expect(isSelected('doc1')).toBe(false);
  });
});

describe('tree context menu', () => {
  const openMenu = async (name: string) => {
    await act(async () => {
      fireEvent.contextMenu(screen.getByText(name));
    });
    return screen.getByRole('menu');
  };

  it('right click on an unselected row selects it and opens the menu', async () => {
    renderTree();

    const menu = await openMenu('doc1');

    expect(isSelected('doc1')).toBe(true);
    expect(within(menu).getByText(/delete 1 document/i)).toBeTruthy();
    expect(within(menu).getByText(/duplicate 1 document/i)).toBeTruthy();
    expect(within(menu).getByText(/cancel/i)).toBeTruthy();
  });

  it('right click on a branch acts on all its documents', async () => {
    renderTree();

    const menu = await openMenu('toolkitA');

    expect(within(menu).getByText(/delete 2 documents/i)).toBeTruthy();
  });

  it('right click keeps an existing multi selection', async () => {
    renderTree();
    click('doc1');
    fireEvent.click(screen.getByText('doc3'), { ctrlKey: true });

    const menu = await openMenu('doc1');

    expect(within(menu).getByText(/delete 2 documents/i)).toBeTruthy();
  });

  it('duplicate copies the selected documents', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: {} });
    renderTree();

    const menu = await openMenu('doc1');
    await act(async () => {
      fireEvent.click(within(menu).getByText(/duplicate 1 document/i));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    const code = mockFetchPython.mock.calls[0][0].code as string;
    expect(code).toContain("duplicateDocumentByID('doc1')");
    expect(mockFetchProjectDetails).toHaveBeenCalledWith('TestProject');
  });

  it('the delete question names a single document', async () => {
    renderTree();

    const menu = await openMenu('doc1');
    await act(async () => {
      fireEvent.click(within(menu).getByText(/delete 1 document/i));
    });

    const dialog = await screen.findByRole('dialog');
    expect(within(dialog).getByText('Delete doc1?')).toBeTruthy();
  });

  it('the delete question counts many documents', async () => {
    renderTree();

    const menu = await openMenu('toolkitA');
    await act(async () => {
      fireEvent.click(within(menu).getByText(/delete 2 documents/i));
    });

    const dialog = await screen.findByRole('dialog');
    expect(within(dialog).getByText('Delete 2 documents?')).toBeTruthy();
  });

  it('delete asks first, then deletes the selected documents', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: {} });
    renderTree();

    const menu = await openMenu('doc1');
    await act(async () => {
      fireEvent.click(within(menu).getByText(/delete 1 document/i));
    });
    const dialog = await screen.findByRole('dialog');
    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /^yes$/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockFetchPython.mock.calls[0][0].code).toContain("All.deleteDocumentByID('doc1')");
  });
});

describe('collectSubtreeKeys', () => {
  it('returns the branch and the documents below it', () => {
    const { viewSettings } = useViewSettingsStore.getState();
    const project = new ProjectObj({ name: 'TestProject', documents } as any);
    const tree = new SplitTree(project.documents, viewSettings.maxDepth, viewSettings);

    const node = findNodeByKey(tree.nodes, 'split_/toolkit=toolkitA');

    expect(node).toBeTruthy();
    expect(collectSubtreeKeys(node!)).toEqual([
      'split_/toolkit=toolkitA',
      idDocId('doc1'),
      idDocId('doc2'),
    ]);
  });

  it('returns undefined for a key that is not in the tree', () => {
    const { viewSettings } = useViewSettingsStore.getState();
    const project = new ProjectObj({ name: 'TestProject', documents } as any);
    const tree = new SplitTree(project.documents, viewSettings.maxDepth, viewSettings);

    expect(findNodeByKey(tree.nodes, 'split_/toolkit=nope')).toBeUndefined();
  });
});
