import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';

const { workflowApi } = await import('../src/api/api');

let getJson: ReturnType<typeof vi.fn>;
let postJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson, postJson } = mockApiCall());
});

describe('workflowApi.start', () => {
  const doc = { _id: { $oid: 'w1' }, desc: { workflowName: 'W' } } as any;

  it('sends the project and the document quietly, and returns the token', async () => {
    postJson.mockResolvedValueOnce({ token: 't1' });

    expect(await workflowApi.start('A', doc)).toEqual({ token: 't1' });
    expect(postJson).toHaveBeenCalledWith('/start_workflow', { projectName: 'A', doc },
      { label: 'run workflow', silent: true });
  });

  it('returns busy as the server said it', async () => {
    postJson.mockResolvedValueOnce({ status: 'busy' });

    expect(await workflowApi.start('A', doc)).toEqual({ status: 'busy' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await workflowApi.start('A', doc)).toBeNull();
  });
});

describe('workflowApi.poll', () => {
  it('puts the token in the path and asks quietly', async () => {
    const status = { status: 'running', error: '', chunks: [{ name: '__preamble__', text: 'go' }] };
    getJson.mockResolvedValueOnce(status);

    expect(await workflowApi.poll('t1')).toEqual(status);
    expect(getJson).toHaveBeenCalledWith('/workflow_status/t1', {}, { label: 'poll workflow', silent: true });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await workflowApi.poll('t1')).toBeNull();
  });
});

describe('workflowApi.runs', () => {
  it('asks quietly for the saved runs of one workflow and answers the list', async () => {
    const runs = [{ id: 'r1', token: 't1', startTime: 's', endTime: 'e', status: 'done' }];
    getJson.mockResolvedValueOnce({ runs });

    expect(await workflowApi.runs('A', 'W')).toEqual(runs);
    expect(getJson).toHaveBeenCalledWith('/workflow_runs', { projectName: 'A', workflowName: 'W' },
      { label: 'list workflow runs', silent: true });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await workflowApi.runs('A', 'W')).toBeNull();
  });
});

describe('workflowApi.run', () => {
  it('asks for one saved run by id', async () => {
    const run = { id: 'r1', token: 't1', startTime: 's', endTime: 'e', status: 'done', error: '', chunks: [] };
    getJson.mockResolvedValueOnce(run);

    expect(await workflowApi.run('A', 'W', 'r1')).toEqual(run);
    expect(getJson).toHaveBeenCalledWith('/workflow_run', { projectName: 'A', workflowName: 'W', runId: 'r1' },
      { label: 'read workflow run' });
  });
});
