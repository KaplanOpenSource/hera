import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act, within } from '@testing-library/react';
import { mockApiCall } from './helpers/mockApiCall';

const { WorkflowOutputPanel } = await import('../src/components/workflow/log/WorkflowOutputPanel');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { shownRun, useWorkflowRunStore, WorkflowRunStatus } = await import('../src/stores/useWorkflowRunStore');

let getJson: ReturnType<typeof vi.fn>;

const summary = (id: string, token: string, startTime: string) => {
  return { id, token, startTime, endTime: startTime, status: 'done' };
};

const savedRun = (id: string, token: string, text: string) => {
  return { ...summary(id, token, 's'), error: '', chunks: [{ name: '__between__', text }], doc: { _id: { $oid: 'd' }, desc: { workflow: {} } } };
};

// Answers the list call with runs, and a run call with the matching saved run.
const serve = (runs: ReturnType<typeof summary>[], byId: { [id: string]: ReturnType<typeof savedRun> }) => {
  getJson.mockImplementation(async (path: string, params: { [key: string]: string }) => {
    if (path === '/workflow_runs') {
      return { runs };
    }
    return byId[params.runId] ?? null;
  });
};

const openSelect = () => {
  fireEvent.click(screen.getByRole('button', { name: 'Past runs' }));
  return within(screen.getByRole('menu'));
};

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson } = mockApiCall());
  useProjectStore.setState({ currProjectName: 'P' });
  useWorkflowRunStore.setState({ runs: {}, oldRuns: {} });
});

afterEach(() => {
  cleanup();
});

describe('WorkflowRunSelect', () => {
  it('lists the latest run and the older saved runs, but not the latest twice', async () => {
    useWorkflowRunStore.setState({ runs: { w: { token: 't2', status: WorkflowRunStatus.Done, error: '', chunks: [] } } });
    serve([summary('r2', 't2', '2026-02-01 10:00:00'), summary('r1', 't1', '2026-01-01 10:00:00')], {});
    render(<WorkflowOutputPanel workflowName="w" />);
    await waitFor(() => { expect(getJson).toHaveBeenCalled(); });

    const list = openSelect();
    expect(list.getByText('Latest run - done')).toBeTruthy();
    expect(list.getByText('2026-01-01 10:00:00 - done')).toBeTruthy();
    expect(list.queryByText('2026-02-01 10:00:00 - done')).toBeNull();
  });

  it('shows an old run when picked, and the latest again when picked back', async () => {
    useWorkflowRunStore.setState({
      runs: { w: { token: 't2', status: WorkflowRunStatus.Done, error: '', chunks: [{ name: '__between__', text: 'latest output' }] } },
    });
    serve([summary('r1', 't1', '2026-01-01 10:00:00')], { r1: savedRun('r1', 't1', 'old output') });
    render(<WorkflowOutputPanel workflowName="w" />);
    expect(screen.getByText('latest output')).toBeTruthy();
    await waitFor(() => { expect(getJson).toHaveBeenCalled(); });

    fireEvent.click(openSelect().getByText('2026-01-01 10:00:00 - done'));
    await waitFor(() => { expect(screen.getByText('old output')).toBeTruthy(); });
    expect(shownRun(useWorkflowRunStore.getState(), 'w')?.token).toBe('t1');

    fireEvent.click(openSelect().getByText('Latest run - done'));
    await waitFor(() => { expect(screen.getByText('latest output')).toBeTruthy(); });
  });

  it('shows the newest saved run when this session has no run', async () => {
    serve([summary('r2', 't2', '2026-02-01 10:00:00'), summary('r1', 't1', '2026-01-01 10:00:00')], {
      r2: savedRun('r2', 't2', 'newest saved output'),
    });
    render(<WorkflowOutputPanel workflowName="w" />);

    await waitFor(() => { expect(screen.getByText('newest saved output')).toBeTruthy(); });
    // It stands in for the latest run, so the canvas keeps today's graph.
    expect(useWorkflowRunStore.getState().oldRuns.w.doc).toBeUndefined();
  });

  it('keeps the saved graph of a run picked from the list', async () => {
    useWorkflowRunStore.setState({ runs: { w: { token: 't2', status: WorkflowRunStatus.Done, error: '', chunks: [] } } });
    serve([summary('r1', 't1', '2026-01-01 10:00:00')], { r1: savedRun('r1', 't1', 'old output') });
    render(<WorkflowOutputPanel workflowName="w" />);
    await waitFor(() => { expect(getJson).toHaveBeenCalled(); });

    fireEvent.click(openSelect().getByText('2026-01-01 10:00:00 - done'));
    await waitFor(() => { expect(useWorkflowRunStore.getState().oldRuns.w?.doc).toBeTruthy(); });
    expect(useWorkflowRunStore.getState().oldRuns.w.startTime).toBe('s');
  });

  it('drops the picked old run when a new run starts', () => {
    useWorkflowRunStore.setState({
      oldRuns: { w: { id: 'r1', token: 't1', startTime: 's', status: WorkflowRunStatus.Done, error: '', chunks: [] } },
    });

    act(() => { useWorkflowRunStore.getState().startRun('w', 't3'); });

    expect(shownRun(useWorkflowRunStore.getState(), 'w')?.token).toBe('t3');
  });

  it('hides the run list until the menu button is clicked', async () => {
    serve([summary('r1', 't1', '2026-01-01 10:00:00')], { r1: savedRun('r1', 't1', 'old output') });
    render(<WorkflowOutputPanel workflowName="w" />);
    await waitFor(() => { expect(screen.getByText('old output')).toBeTruthy(); });

    expect(screen.queryByRole('menu')).toBeNull();
    expect(openSelect().getByText('2026-01-01 10:00:00 - done')).toBeTruthy();
  });

  it('shows the number of past runs on the menu button', async () => {
    useWorkflowRunStore.setState({ runs: { w: { token: 't3', status: WorkflowRunStatus.Done, error: '', chunks: [] } } });
    serve([
      summary('r3', 't3', '2026-03-01 10:00:00'),
      summary('r2', 't2', '2026-02-01 10:00:00'),
      summary('r1', 't1', '2026-01-01 10:00:00'),
    ], {});
    render(<WorkflowOutputPanel workflowName="w" />);

    await waitFor(() => { expect(within(screen.getByTestId('past-runs-badge')).getByText('2')).toBeTruthy(); });
  });

  it('says in the tooltip that the button holds the past runs', async () => {
    serve([summary('r2', 't2', '2026-02-01 10:00:00'), summary('r1', 't1', '2026-01-01 10:00:00')], {
      r2: savedRun('r2', 't2', 'newest saved output'),
    });
    render(<WorkflowOutputPanel workflowName="w" />);
    await waitFor(() => { expect(screen.getByText('newest saved output')).toBeTruthy(); });

    fireEvent.mouseOver(screen.getByRole('button', { name: 'Past runs' }));
    await waitFor(() => {
      expect(screen.getByText('Past runs: 2 saved runs. Click to see the output of an earlier run.')).toBeTruthy();
    });
  });

  it('asks nothing when no project is open', () => {
    useProjectStore.setState({ currProjectName: '* NONE *' });
    render(<WorkflowOutputPanel workflowName="w" />);

    expect(getJson).not.toHaveBeenCalled();
    expect(screen.getByText('No run yet.')).toBeTruthy();
  });
});
