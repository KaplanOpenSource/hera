import { describe, it, expect, beforeEach, afterEach, vi } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockGetJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  getJson: (...args: any[]) => mockGetJson(...args),
  postJson: vi.fn(),
}));

import { CentralRepoFolder } from '../src/components/repo/CentralRepoFolder';

const renderInTree = ({
  expanded = false,
}: {
  expanded?: boolean,
} = {}) =>
  render(
    <SimpleTreeView defaultExpandedItems={expanded ? ['central-repo-folder'] : []}>
      <CentralRepoFolder />
    </SimpleTreeView>,
  );

const clickTreeItem = async () => {
  const treeItem = screen.getByRole('tree').querySelector(':scope > [role="treeitem"]')!;
  await act(async () => { fireEvent.click(treeItem); });
};

beforeEach(() => {
  vi.clearAllMocks();
  localStorage.clear();
});

afterEach(() => { cleanup(); });

describe('CentralRepoFolder', () => {
  it('renders default folder path', () => {
    renderInTree();
    expect(screen.getByText('~/hera/repositories/')).toBeTruthy();
  });

  it('uses folder from localStorage when set', () => {
    localStorage.setItem('hera-central-repo-folder', '/custom/path/');
    renderInTree();
    expect(screen.getByText('/custom/path/')).toBeTruthy();
  });

  it('shows "Loading..." when expanded but not yet fetched', () => {
    renderInTree({ expanded: true });
    expect(screen.getByText('Loading...')).toBeTruthy();
  });
});

describe('CentralRepoFolder — the call it makes', () => {
  it('asks the folder entrypoint for the names only on expand', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [] });
    renderInTree();
    await clickTreeItem();

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledTimes(1);
    });
    expect(mockGetJson.mock.calls[0][0]).toBe('/repository/jsonsInFolder');
    expect(mockGetJson.mock.calls[0][1]).toEqual({
      folder: '~/hera/repositories/',
      withContents: 'false',
    });
  });

  it('sends the folder the user set', async () => {
    localStorage.setItem('hera-central-repo-folder', '/my/repos/');
    mockGetJson.mockResolvedValueOnce({ files: [] });
    renderInTree();
    await clickTreeItem();

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledTimes(1);
    });
    expect(mockGetJson.mock.calls[0][1].folder).toBe('/my/repos/');
  });
});

describe('CentralRepoFolder — response handling', () => {
  it('renders returned file paths as tree items', async () => {
    mockGetJson.mockResolvedValueOnce({ files: ['/repos/a.json', '/repos/b.json'] });
    renderInTree({ expanded: true });
    await clickTreeItem();

    await waitFor(() => {
      expect(screen.getByText('/repos/a.json')).toBeTruthy();
      expect(screen.getByText('/repos/b.json')).toBeTruthy();
    });
  });

  it('shows "No JSON files found" when server returns empty list', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [] });
    renderInTree({ expanded: true });
    await clickTreeItem();

    await waitFor(() => {
      expect(screen.getByText('No JSON files found')).toBeTruthy();
    });
  });

  it('shows "No JSON files found" when data is null', async () => {
    mockGetJson.mockResolvedValueOnce(null);
    renderInTree({ expanded: true });
    await clickTreeItem();

    await waitFor(() => {
      expect(screen.getByText('No JSON files found')).toBeTruthy();
    });
  });

  it('does not re-fetch on second click', async () => {
    mockGetJson.mockResolvedValueOnce({ files: ['/repos/a.json'] });
    renderInTree({ expanded: true });
    await clickTreeItem();
    await waitFor(() => { expect(screen.getByText('/repos/a.json')).toBeTruthy(); });

    await clickTreeItem();
    expect(mockGetJson).toHaveBeenCalledTimes(1);
  });
});

describe('CentralRepoFolder — change folder dialog', () => {
  it('opens settings dialog and re-fetches with new folder', async () => {
    mockGetJson.mockResolvedValue({ files: [] });
    renderInTree();

    const settingsWrapper = screen.getByLabelText('Change repository folder');
    await act(async () => {
      fireEvent.click(within(settingsWrapper).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    const input = within(dialog).getByRole('textbox', { name: /folder path/i });
    fireEvent.change(input, { target: { value: '/new/folder/' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /save/i }));
    });

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalled();
      const lastCall = mockGetJson.mock.calls[mockGetJson.mock.calls.length - 1];
      expect(lastCall[1].folder).toBe('/new/folder/');
    });

    expect(localStorage.getItem('hera-central-repo-folder')).toBe('/new/folder/');
    expect(screen.getByText('/new/folder/')).toBeTruthy();
  });
});
