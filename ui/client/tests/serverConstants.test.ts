import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockExecPython = vi.fn();
vi.mock('../src/io/execPython', () => ({
  execPython: (...args: any[]) => mockExecPython(...args),
}));

const { useServerConstants } = await import('../src/stores/useServerConstants');

beforeEach(() => {
  vi.clearAllMocks();
  useServerConstants.setState({ dataTypes: {} });
});

describe('readAllConstants', () => {
  it('sets dataTypes in store on success', async () => {
    const datatypes = { STRING: 'string', INTEGER: 'integer' };
    mockExecPython.mockResolvedValueOnce({
      data: { datatypes },
      problem: undefined,
    });

    await useServerConstants.getState().readAllConstants();

    expect(useServerConstants.getState().dataTypes).toEqual(datatypes);
  });

  it('does not update store when data is falsy', async () => {
    mockExecPython.mockResolvedValueOnce({
      data: undefined,
      problem: 'error',
    });

    await useServerConstants.getState().readAllConstants();

    expect(useServerConstants.getState().dataTypes).toEqual({});
  });

  it('sends correct Python code', async () => {
    mockExecPython.mockResolvedValueOnce({ data: { datatypes: {} }, problem: undefined });

    await useServerConstants.getState().readAllConstants();

    const code = mockExecPython.mock.calls[0][0];
    expect(code).toContain('from hera import datalayer');
    expect(code).toContain('datatypes');
  });

  it('works with delayed response', async () => {
    const datatypes = { FLOAT: 'float' };
    mockExecPython.mockImplementationOnce(() =>
      new Promise(resolve => setTimeout(() => resolve({
        data: { datatypes },
        problem: undefined,
      }), 50))
    );

    await useServerConstants.getState().readAllConstants();

    expect(useServerConstants.getState().dataTypes).toEqual(datatypes);
  });
});
