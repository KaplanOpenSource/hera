import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const mockFetchProjectDetails = vi.fn();
vi.mock('../src/io/FetchProjects', () => ({
  fetchProjectDetails: (...args: any[]) => mockFetchProjectDetails(...args),
}));

// These two rows talk to python on their own, which is not what we test here.
vi.mock('../src/components/repo/CentralRepoFolder', () => ({
  CentralRepoFolder: () => null,
}));
vi.mock('../src/components/repo/RegisteredRepositories', () => ({
  RegisteredRepositories: () => null,
}));

const { RepoTreeWhole } = await import('../src/components/project/RepoTreeWhole');
const { useProjectStore } = await import('../src/stores/useProjectStore');

// The one repository row that RepoTreeWhole starts with.
const REPO_PATH = 'hera/doc/jupyter/Developer/Documentation_Repository.json';

const renderRepos = (onOpenItem = vi.fn()) => {
  const result = render(
    <RepoTreeWhole
      selectedIds={[]}
      onSelectedItemsChange={() => {}}
      expandedItems={[]}
      onExpandedItemsChange={() => {}}
      onOpenItem={onOpenItem}
    />
  );
  return { ...result, onOpenItem };
};

// The row shows the path in two parts, so look for its last segment.
const repoRow = () => {
  return screen.getByText('/Documentation_Repository.json');
};

const answer = async (button: RegExp) => {
  const dialog = await screen.findByRole('dialog');
  await act(async () => {
    fireEvent.click(within(dialog).getByRole('button', { name: button }));
  });
  return dialog;
};

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({ currProjectName: 'TestProject' });
});

afterEach(() => {
  cleanup();
});

describe('repositories tree', () => {
  it('a double click offers to load the repository into the project', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: {} });
    renderRepos();

    await act(async () => {
      fireEvent.doubleClick(repoRow());
    });
    await answer(/^yes$/i);

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockFetchPython.mock.calls[0][0].code).toContain(`repositoryName='${REPO_PATH}'`);
    expect(mockFetchProjectDetails).toHaveBeenCalledWith('TestProject');
  });

  it('answering no loads nothing', async () => {
    renderRepos();

    await act(async () => {
      fireEvent.doubleClick(repoRow());
    });
    await answer(/^no$/i);

    expect(mockFetchPython).not.toHaveBeenCalled();
  });

  it('a right click offers the same load, and a cancel', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: {} });
    renderRepos();

    await act(async () => {
      fireEvent.contextMenu(repoRow());
    });
    const menu = screen.getByRole('menu');
    expect(within(menu).getByText(/cancel/i)).toBeTruthy();
    await act(async () => {
      fireEvent.click(within(menu).getByText(/into project/i));
    });
    await answer(/^yes$/i);

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockFetchPython.mock.calls[0][0].code).toContain(`repositoryName='${REPO_PATH}'`);
  });

  it('cancel in the menu does nothing', async () => {
    renderRepos();

    await act(async () => {
      fireEvent.contextMenu(repoRow());
    });
    await act(async () => {
      fireEvent.click(within(screen.getByRole('menu')).getByText(/cancel/i));
    });

    expect(mockFetchPython).not.toHaveBeenCalled();
    expect(screen.queryByRole('dialog')).toBeNull();
  });
});
