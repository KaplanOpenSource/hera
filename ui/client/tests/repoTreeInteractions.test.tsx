import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

// These two rows talk to python on their own, which is not what we test here.
vi.mock('../src/components/repo/CentralRepoFolder', () => ({
  CentralRepoFolder: () => null,
}));
vi.mock('../src/components/repo/RegisteredRepositories', () => ({
  RegisteredRepositories: () => null,
}));

const mockPushError = vi.fn();
const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushError: (...args: any[]) => mockPushError(...args),
  pushInfo: (...args: any[]) => mockPushInfo(...args),
  pushRunning: vi.fn(),
  dismiss: vi.fn(),
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
  it('a double click opens the repository, so its contents show', async () => {
    const { onOpenItem } = renderRepos();

    await act(async () => {
      fireEvent.click(repoRow());
      fireEvent.doubleClick(repoRow());
    });

    expect(onOpenItem).toHaveBeenCalledWith(`__repo_*_${REPO_PATH}`);
    expect(screen.queryByRole('dialog')).toBeNull();
    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('answering no loads nothing', async () => {
    renderRepos();

    await act(async () => {
      fireEvent.contextMenu(repoRow());
    });
    await act(async () => {
      fireEvent.click(within(screen.getByRole('menu')).getByText(/into project/i));
    });
    await answer(/^no$/i);

    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('a right click offers the same load, and a cancel', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] }, added: 2, problems: [] });
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
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostJson.mock.calls[0][0]).toBe('/repository/loadIntoProject');
    // The row carries its path, so nothing has to guess a name apart from a file.
    expect(mockPostJson.mock.calls[0][1]).toEqual({
      projectName: 'TestProject',
      path: REPO_PATH,
      overwrite: false,
    });
  });

  it('cancel in the menu does nothing', async () => {
    renderRepos();

    await act(async () => {
      fireEvent.contextMenu(repoRow());
    });
    await act(async () => {
      fireEvent.click(within(screen.getByRole('menu')).getByText(/cancel/i));
    });

    expect(mockPostJson).not.toHaveBeenCalled();
    expect(screen.queryByRole('dialog')).toBeNull();
  });

  it('loading a repository is one call, and stores the project it answers', async () => {
    const loaded = { name: 'TestProject', documents: [] };
    mockPostJson.mockResolvedValueOnce({ project: loaded, added: 2, problems: [] });
    renderRepos();

    await act(async () => {
      fireEvent.contextMenu(repoRow());
    });
    await act(async () => {
      fireEvent.click(within(screen.getByRole('menu')).getByText(/into project/i));
    });
    await answer(/^yes$/i);

    await waitFor(() => {
      expect(useProjectStore.getState().currProject).toEqual(loaded);
    });
    expect(mockPostJson).toHaveBeenCalledTimes(1);
  });
});
