import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPushError = vi.fn();
const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushError: (...args: any[]) => mockPushError(...args),
  pushInfo: (...args: any[]) => mockPushInfo(...args),
  pushRunning: vi.fn(),
  dismiss: vi.fn(),
}));

const { UpdateRepositoriesButton } = await import('../src/components/repo/UpdateRepositoriesButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { mockApiCall } = await import('./helpers/mockApiCall');

afterEach(() => { cleanup(); });

let mockPostJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ postJson: mockPostJson } = mockApiCall());
  useProjectStore.setState({ currProjectName: 'TestProject', currProject: undefined });
});

const clickUpdate = async () => {
  const wrapper = screen.getByLabelText('Update repositories into project');
  await act(async () => { fireEvent.click(within(wrapper).getByRole('button')); });
};

describe('UpdateRepositoriesButton', () => {
  it('sends no source, which means every registered repository', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [] }, added: 1, problems: [],
    });

    render(<UpdateRepositoriesButton />);
    await clickUpdate();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostJson.mock.calls[0][0]).toBe('/repository/loadIntoProject');
    expect(mockPostJson.mock.calls[0][1]).toEqual({
      projectName: 'TestProject',
      overwrite: false,
    });
  });

  it('stores the project the answer carries, and makes no second call', async () => {
    const updated = {
      name: 'TestProject',
      documents: [
        { _id: { $oid: 'new1' }, desc: {}, type: 'T', _cls: 'M', resource: '', dataFormat: 'string', projectName: 'TestProject' },
      ],
    };
    mockPostJson.mockResolvedValueOnce({ project: updated, added: 1, problems: [] });

    render(<UpdateRepositoriesButton />);
    await clickUpdate();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject).toEqual(updated);
    });
    expect(mockPostJson).toHaveBeenCalledTimes(1);
  });

  it('stores nothing when the call failed', async () => {
    mockPostJson.mockResolvedValueOnce(null);

    render(<UpdateRepositoriesButton />);
    await clickUpdate();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProject).toBeUndefined();
  });

  it('says nothing was added, and why, when hera skipped everything', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [] },
      added: 0,
      problems: ['GIS_Tiles: ImportError: Cannot locate class: x.y.TilesToolkit'],
    });

    render(<UpdateRepositoriesButton />);
    await clickUpdate();

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledTimes(1);
    });
    expect(mockPushError.mock.calls[0][0]).toContain('nothing was added');
    expect(mockPushError.mock.calls[0][0]).toContain('GIS_Tiles');
  });

  it('invents no reason when hera reported none', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [] }, added: 0, problems: [],
    });

    render(<UpdateRepositoriesButton />);
    await clickUpdate();

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledTimes(1);
    });
    expect(mockPushError.mock.calls[0][0]).toContain('hera reported nothing');
  });

  it('a load that added something but skipped some of it says both', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [] },
      added: 3,
      problems: ['experiment: ValueError: no such toolkit'],
    });

    render(<UpdateRepositoriesButton />);
    await clickUpdate();

    await waitFor(() => {
      expect(mockPushInfo).toHaveBeenCalledTimes(1);
    });
    expect(mockPushInfo.mock.calls[0][0]).toContain('added 3, 1 skipped');
    expect(mockPushError).not.toHaveBeenCalled();
  });
});
