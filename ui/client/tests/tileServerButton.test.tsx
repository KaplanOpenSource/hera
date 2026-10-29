import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, fireEvent, screen, waitFor } from '@testing-library/react';
import { mockApiCall } from './helpers/mockApiCall';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { TileServerButton } = await import('../src/components/details/maps/TileServerButton');
const { TILE_SERVERS } = await import('../src/components/details/maps/tileServers');
const { useProjectStore } = await import('../src/stores/useProjectStore');

const NEW_URL = 'https://tile.opentopomap.org/{z}/{x}/{y}.png';

const configDoc = {
  _id: { $oid: 'cfg' },
  _cls: 'Metadata.Cache',
  projectName: 'P',
  desc: { other: 1 },
  type: 'P__config__',
  resource: '',
  dataFormat: 'string',
};

let mockPostJson: ReturnType<typeof vi.fn>;

const setDocs = (documents: object[]) => {
  useProjectStore.setState({
    currProjectName: 'P',
    currProject: { name: 'P', documents: documents as any },
    editedDocs: {},
  });
};

const openAndSave = async () => {
  render(<TileServerButton currentUrl="https://old/{z}/{x}/{y}.png" />);
  fireEvent.click(screen.getByRole('button'));
  const input = await screen.findByRole('combobox');
  expect((input as HTMLInputElement).value).toBe('https://old/{z}/{x}/{y}.png');
  fireEvent.change(input, { target: { value: NEW_URL } });
  fireEvent.click(screen.getByText('Save'));
  await waitFor(() => expect(mockPostJson).toHaveBeenCalled());
};

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.clearAllMocks();
  mockPostJson = mockApiCall().postJson;
  mockPostJson.mockResolvedValue({ project: { name: 'P', documents: [] }, document: {}, addedDocId: 'x' });
});

describe('TileServerButton', () => {
  it('updates defaultTileServer on the existing config document', async () => {
    setDocs([configDoc]);

    await openAndSave();

    const [path, body] = mockPostJson.mock.calls[0];
    expect(path).toBe('/document/update');
    expect(body.docId).toBe('cfg');
    expect(body.fields.desc).toEqual({ other: 1, defaultTileServer: NEW_URL });
  });

  it('shows the server list as soon as the dialog opens', async () => {
    setDocs([configDoc]);
    render(<TileServerButton currentUrl="https://old/{z}/{x}/{y}.png" />);
    fireEvent.click(screen.getByRole('button'));

    expect(await screen.findByText('OpenTopoMap')).toBeTruthy();
  });

  it('saves a server picked from the list', async () => {
    setDocs([configDoc]);
    render(<TileServerButton currentUrl="https://old/{z}/{x}/{y}.png" />);
    fireEvent.click(screen.getByRole('button'));
    fireEvent.mouseDown(await screen.findByRole('combobox'));
    fireEvent.click(await screen.findByText('OpenTopoMap'));
    fireEvent.click(screen.getByText('Save'));
    await waitFor(() => expect(mockPostJson).toHaveBeenCalled());

    expect(mockPostJson.mock.calls[0][1].fields.desc.defaultTileServer).toBe(TILE_SERVERS[2].url);
  });

  it('creates the config document when the project has none', async () => {
    setDocs([]);

    await openAndSave();

    const [path, body] = mockPostJson.mock.calls[0];
    expect(path).toBe('/document/create');
    expect(body).toMatchObject({
      projectName: 'P',
      collection: 'Cache',
      type: 'P__config__',
      desc: { defaultTileServer: NEW_URL },
    });
  });
});
