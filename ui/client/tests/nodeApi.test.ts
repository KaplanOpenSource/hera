import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';

const { nodeApi } = await import('../src/api/api');

let getJson: ReturnType<typeof vi.fn>;
let postJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson, postJson } = mockApiCall());
});

describe('nodeApi.catalog', () => {
  it('asks for the catalog and returns it', async () => {
    const catalog = [{ type: 'general.Parameters', parameters: [] }];
    getJson.mockResolvedValueOnce(catalog);

    expect(await nodeApi.catalog()).toEqual(catalog);
    expect(getJson).toHaveBeenCalledWith('/node-catalog', {}, { label: 'node catalog' });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await nodeApi.catalog()).toBeNull();
  });
});

describe('nodeApi.validateParams', () => {
  it('sends the type and the values, and returns the check', async () => {
    postJson.mockResolvedValueOnce({ ok: false, message: 'Agent is missing' });

    const answer = await nodeApi.validateParams({ type: 'RiskAssessment.x', params: { Agent: '' } });

    expect(answer).toEqual({ ok: false, message: 'Agent is missing' });
    expect(postJson).toHaveBeenCalledWith('/node-params/validate',
      { type: 'RiskAssessment.x', params: { Agent: '' } },
      { label: 'node params validation' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await nodeApi.validateParams({ type: 'x', params: {} })).toBeNull();
  });
});
