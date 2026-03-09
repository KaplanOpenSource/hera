import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const { useServerConstants } = await import('../src/stores/useServerConstants');

beforeEach(() => {
  vi.clearAllMocks();
  useServerConstants.setState({ dataTypes: {} });
});

describe('readAllConstants', () => {
  it('sets dataTypes in store on success', async () => {
    const datatypes = { STRING: 'string', INTEGER: 'integer' };
    mockFetchPython.mockResolvedValueOnce({
      data: { datatypes },
      problem: undefined,
    });

    await useServerConstants.getState().readAllConstants();

    expect(useServerConstants.getState().dataTypes).toEqual(datatypes);
  });

  it('does not update store when data is falsy', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: undefined,
      problem: 'error',
    });

    await useServerConstants.getState().readAllConstants();

    expect(useServerConstants.getState().dataTypes).toEqual({});
  });

  it('sends correct Python code', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: { datatypes: {} }, problem: undefined });

    await useServerConstants.getState().readAllConstants();

    const code = mockFetchPython.mock.calls[0][0].code;
    expect(code).toContain('from hera import datalayer');
    expect(code).toContain('datatypes');
  });

  it('store not updated while response is in flight', async () => {
    const datatypes = { FLOAT: 'float' };
    let resolve!: (v: any) => void;
    mockFetchPython.mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    const promise = useServerConstants.getState().readAllConstants();
    expect(useServerConstants.getState().dataTypes).toEqual({});

    resolve({ data: { datatypes }, problem: undefined });
    await promise;

    expect(useServerConstants.getState().dataTypes).toEqual(datatypes);
  });
});
