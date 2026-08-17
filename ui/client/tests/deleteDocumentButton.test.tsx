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

const { DeleteDocumentButton } = await import('../src/components/details/DeleteDocumentButton');

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.clearAllMocks();
  mockFetchProjectDetails.mockResolvedValue(undefined);
});

const document = {
  _id: { $oid: 'doc1' },
  desc: { datasourceName: 'MyDoc' },
  type: 'TestType',
  _cls: 'Metadata.Cache',
  resource: 'res',
  dataFormat: 'string',
  projectName: 'TestProject',
} as any;

const renderButton = () => render(
  <DeleteDocumentButton
    document={document}
    projectName="TestProject"
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
  it('calls execPython to delete document on confirm', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: {}, problem: undefined });
    renderButton();

    await clickDeleteAndConfirm();

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("All.deleteDocumentByID('doc1')");
    });
  });

  it('calls fetchProjectDetails after successful delete', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: {}, problem: undefined });
    renderButton();

    await clickDeleteAndConfirm();

    await waitFor(() => {
      expect(mockFetchProjectDetails).toHaveBeenCalledWith('TestProject');
    });
  });

  it('does not call fetchProjectDetails when delete fails', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: undefined, problem: 'error' });
    renderButton();

    await clickDeleteAndConfirm();

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockFetchProjectDetails).not.toHaveBeenCalled();
  });

  it('is disabled for the project config document', async () => {
    render(
      <DeleteDocumentButton
        document={document}
        projectName="TestProject"
        displayName="TestProject config"
        isConfig={true}
      />
    );

    const wrapper = screen.getByLabelText('Project Config is deleted only with project');
    const button = within(wrapper).getByRole('button') as HTMLButtonElement;
    expect(button.disabled).toBe(true);
  });

  it('handles delayed response', async () => {
    let resolveExec!: (v: any) => void;
    mockFetchPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );
    renderButton();

    await clickDeleteAndConfirm();

    // execPython called but not yet resolved
    expect(mockFetchPython).toHaveBeenCalledTimes(1);
    expect(mockFetchProjectDetails).not.toHaveBeenCalled();

    // Resolve after delay
    await act(async () => {
      resolveExec({ data: {}, problem: undefined });
    });

    await waitFor(() => {
      expect(mockFetchProjectDetails).toHaveBeenCalledWith('TestProject');
    });
  });
});
