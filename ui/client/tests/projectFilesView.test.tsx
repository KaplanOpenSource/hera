import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, cleanup, fireEvent, waitFor } from '@testing-library/react';
import { ProjectObj } from '../src/objects/ProjectObj';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { ProjectFilesView } = await import('../src/components/details/ProjectFilesView');
const { mockApiCall } = await import('./helpers/mockApiCall');

let getJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson } = mockApiCall());
});

afterEach(() => {
  cleanup();
});

const makeProject = (filesDirectory?: string): ProjectObj => {
  let desc = {};
  if (filesDirectory) {
    desc = { filesDirectory };
  }
  return new ProjectObj({
    name: 'P',
    documents: [{
      _cls: 'Config',
      _id: { $oid: 'config-id' },
      projectName: 'P',
      desc,
      type: 'P__config__',
      resource: '',
      dataFormat: 'json',
    }],
  });
};

describe('ProjectFilesView', () => {
  it('shows the project folder contents', async () => {
    getJson.mockResolvedValueOnce({ files: [{ name: 'run', isFolder: true }, { name: 'a.txt', isFolder: false }] });

    render(<ProjectFilesView project={makeProject('/data/p')} />);

    expect(await screen.findByText('run')).toBeDefined();
    expect(screen.getByText('a.txt')).toBeDefined();
    expect(getJson).toHaveBeenCalledWith('/project/files', { projectName: 'P', path: '' }, { label: 'project files' });
  });

  it('reads a subfolder only when it is opened', async () => {
    getJson.mockResolvedValueOnce({ files: [{ name: 'run', isFolder: true }] });
    getJson.mockResolvedValueOnce({ files: [{ name: 'out.nc', isFolder: false }] });

    render(<ProjectFilesView project={makeProject('/data/p')} />);
    fireEvent.click(await screen.findByText('run'));

    expect(await screen.findByText('out.nc')).toBeDefined();
    expect(getJson).toHaveBeenCalledTimes(2);
    expect(getJson.mock.calls[1][1]).toEqual({ projectName: 'P', path: 'run' });
  });

  it('says so when the folder is empty or not found', async () => {
    getJson.mockResolvedValueOnce({ files: [] });

    render(<ProjectFilesView project={makeProject('/data/gone')} />);

    expect(await screen.findByText('The folder is empty or not found.')).toBeDefined();
  });

  it('asks nothing when the project has no files directory', async () => {
    render(<ProjectFilesView project={makeProject()} />);

    expect(screen.getByText('This project has no files directory.')).toBeDefined();
    await waitFor(() => expect(getJson).not.toHaveBeenCalled());
  });
});
