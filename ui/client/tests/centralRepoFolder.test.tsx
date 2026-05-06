import { describe, it, expect, beforeEach, afterEach, vi } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
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

describe('CentralRepoFolder — fetchFiles Python code', () => {
  it('sends isRepoJson filter code on expand', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: { jsonFiles: [] } });
    renderInTree();
    await clickTreeItem();

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });

    const call = mockFetchPython.mock.calls[0][0];
    expect(call.code).toContain('def isRepoJson(path)');
    expect(call.code).toContain('SECTIONS');
    expect(call.code).toContain("jsonFiles = sorted");
    expect(call.results).toEqual(['jsonFiles']);
  });

  it('interpolates the folder path into the Python code', async () => {
    localStorage.setItem('hera-central-repo-folder', '/my/repos/');
    mockFetchPython.mockResolvedValueOnce({ data: { jsonFiles: [] } });
    renderInTree();
    await clickTreeItem();

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockFetchPython.mock.calls[0][0].code).toContain("/my/repos/");
  });

  it('excludes caseConfiguration.json in the Python code', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: { jsonFiles: [] } });
    renderInTree();
    await clickTreeItem();

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockFetchPython.mock.calls[0][0].code).toContain("caseConfiguration.json");
  });
});

describe('CentralRepoFolder — response handling', () => {
  it('renders returned file paths as tree items', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { jsonFiles: ['/repos/a.json', '/repos/b.json'] },
    });
    renderInTree({ expanded: true });
    await clickTreeItem();

    await waitFor(() => {
      expect(screen.getByText('/repos/a.json')).toBeTruthy();
      expect(screen.getByText('/repos/b.json')).toBeTruthy();
    });
  });

  it('shows "No JSON files found" when server returns empty list', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: { jsonFiles: [] } });
    renderInTree({ expanded: true });
    await clickTreeItem();

    await waitFor(() => {
      expect(screen.getByText('No JSON files found')).toBeTruthy();
    });
  });

  it('shows "No JSON files found" when data is null', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: null });
    renderInTree({ expanded: true });
    await clickTreeItem();

    await waitFor(() => {
      expect(screen.getByText('No JSON files found')).toBeTruthy();
    });
  });

  it('does not re-fetch on second click', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { jsonFiles: ['/repos/a.json'] },
    });
    renderInTree({ expanded: true });
    await clickTreeItem();
    await waitFor(() => { expect(screen.getByText('/repos/a.json')).toBeTruthy(); });

    await clickTreeItem();
    expect(mockFetchPython).toHaveBeenCalledTimes(1);
  });
});

describe('CentralRepoFolder — change folder dialog', () => {
  it('opens settings dialog and re-fetches with new folder', async () => {
    mockFetchPython.mockResolvedValue({ data: { jsonFiles: [] } });
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
      expect(mockFetchPython).toHaveBeenCalled();
      const lastCall = mockFetchPython.mock.calls[mockFetchPython.mock.calls.length - 1][0];
      expect(lastCall.code).toContain('/new/folder/');
    });

    expect(localStorage.getItem('hera-central-repo-folder')).toBe('/new/folder/');
    expect(screen.getByText('/new/folder/')).toBeTruthy();
  });
});
