import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { mockApiCall } from './helpers/mockApiCall';

const { ExportProjectButton } = await import('../src/components/header/ExportProjectButton');
const { useProjectStore, DEFAULT_PROJECT } = await import('../src/stores/useProjectStore');

let mockPostBytes: ReturnType<typeof vi.fn>;
let createObjectURL: ReturnType<typeof vi.fn>;
let revokeObjectURL: ReturnType<typeof vi.fn>;
let clickedDownloads: string[];

beforeEach(() => {
  vi.clearAllMocks();
  ({ postBytes: mockPostBytes } = mockApiCall());
  createObjectURL = vi.fn().mockReturnValue('blob:zip');
  revokeObjectURL = vi.fn();
  URL.createObjectURL = createObjectURL as any;
  URL.revokeObjectURL = revokeObjectURL as any;
  clickedDownloads = [];
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(function (this: HTMLAnchorElement) {
    clickedDownloads.push(this.download);
  });
  useProjectStore.setState({
    currProjectName: 'Alpha',
    currProject: { name: 'Alpha', documents: [] },
  });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const clickExport = async () => {
  const wrapper = screen.getByLabelText('Export project');
  await act(async () => { fireEvent.click(within(wrapper).getByRole('button')); });
};

describe('ExportProjectButton', () => {
  it('asks the server for the zip and saves it as <project>.zip', async () => {
    mockPostBytes.mockResolvedValueOnce(new Uint8Array([1, 2, 3]).buffer);
    const onExported = vi.fn();
    render(<ExportProjectButton onExported={onExported} />);

    await clickExport();

    await waitFor(() => {
      expect(clickedDownloads).toEqual(['Alpha.zip']);
    });
    expect(mockPostBytes).toHaveBeenCalledWith('/project/export',
      { projectName: 'Alpha' }, { label: 'export project Alpha' });
    expect((createObjectURL.mock.calls[0][0] as Blob).size).toBe(3);
    expect(onExported).toHaveBeenCalled();
  });

  it('saves nothing when the export failed', async () => {
    mockPostBytes.mockResolvedValueOnce(null);
    const onExported = vi.fn();
    render(<ExportProjectButton onExported={onExported} />);

    await clickExport();

    await waitFor(() => {
      expect(mockPostBytes).toHaveBeenCalledTimes(1);
    });
    expect(clickedDownloads).toEqual([]);
    expect(onExported).not.toHaveBeenCalled();
  });

  it('is disabled for the default project', () => {
    useProjectStore.setState({ currProjectName: DEFAULT_PROJECT });
    render(<ExportProjectButton />);

    const btn = within(screen.getByLabelText('Export project')).getByRole('button') as HTMLButtonElement;
    expect(btn.disabled).toBe(true);
  });
});
