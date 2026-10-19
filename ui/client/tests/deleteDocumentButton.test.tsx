import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { mockApiCall } from './helpers/mockApiCall';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { DeleteDocumentButton } = await import('../src/components/details/DeleteDocumentButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => { cleanup(); });

let mockPostJson: ReturnType<typeof vi.fn>;

const document = {
  _id: { $oid: 'doc1' },
  desc: { datasourceName: 'MyDoc' },
  type: 'TestType',
  _cls: 'Metadata.Cache',
  resource: 'res',
  dataFormat: 'string',
  projectName: 'TestProject',
} as any;

beforeEach(() => {
  vi.clearAllMocks();
  mockPostJson = mockApiCall().postJson;
  useProjectStore.setState({
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents: [document] },
    editedDocs: {},
  });
});

const renderButton = () => render(
  <DeleteDocumentButton
    document={document}
    displayName="MyDoc"
    isConfig={false}
  />
);

const clickDeleteAndConfirm = async () => {
  const wrapper = screen.getByLabelText('Delete Document');
  await act(async () => { fireEvent.click(within(wrapper).getByRole('button')); });
  const yesBtn = await screen.findByRole('button', { name: /yes/i });
  await act(async () => { fireEvent.click(yesBtn); });
};

describe('DeleteDocumentButton', () => {
  it('calls document.delete with the one id on confirm', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] } });
    renderButton();

    await clickDeleteAndConfirm();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostJson.mock.calls[0][0]).toBe('/document/delete');
    expect(mockPostJson.mock.calls[0][1]).toEqual({ docIds: ['doc1'] });
  });

  it('stores the project the answer carries, with no second call', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] } });
    renderButton();

    await clickDeleteAndConfirm();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toEqual([]);
    });
    expect(mockPostJson).toHaveBeenCalledTimes(1);
  });

  it('leaves the project alone when the delete fails', async () => {
    mockPostJson.mockResolvedValueOnce(null);
    renderButton();

    await clickDeleteAndConfirm();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProject?.documents).toEqual([document]);
  });

  it('is disabled for the project config document', async () => {
    render(
      <DeleteDocumentButton
        document={document}
        displayName="TestProject config"
        isConfig={true}
      />
    );

    const wrapper = screen.getByLabelText('Project Config is deleted only with project');
    const button = within(wrapper).getByRole('button') as HTMLButtonElement;
    expect(button.disabled).toBe(true);
  });

  it('handles delayed response', async () => {
    let resolveCall!: (v: any) => void;
    mockPostJson.mockImplementationOnce(() => new Promise(r => { resolveCall = r; }));
    renderButton();

    await clickDeleteAndConfirm();

    expect(mockPostJson).toHaveBeenCalledTimes(1);
    expect(useProjectStore.getState().currProject?.documents).toEqual([document]);

    await act(async () => {
      resolveCall({ project: { name: 'TestProject', documents: [] } });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toEqual([]);
    });
  });
});
