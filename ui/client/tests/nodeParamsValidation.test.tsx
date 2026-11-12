import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { act, cleanup, render, screen } from '@testing-library/react';

// The node renders ReactFlow Handles, which need the ReactFlow store.
vi.mock('@xyflow/react', () => ({
  Handle: () => null,
  NodeResizer: () => null,
  Position: { Left: 'left', Right: 'right', Top: 'top', Bottom: 'bottom' },
  useUpdateNodeInternals: () => () => {},
}));

import { mockApiCall } from './helpers/mockApiCall';
import { WorkflowFlowNode } from '../src/components/workflow/WorkflowFlowNode';

let postJson: ReturnType<typeof vi.fn>;

const catalog = [{ type: 'RiskAssessment.calculateThresholds', parameters: [] }];

const nodeElement = (node: any) => {
  return (
    <WorkflowFlowNode
      // The message shows in the editor, which opens on hover.
      data={{ name: 'node1', node, catalog, onRename: vi.fn(), onChange: vi.fn(), onDelete: vi.fn(), expanded: true }}
      selected={false}
      {...({} as any)}
    />
  );
};

const renderNode = (node: any) => {
  return render(nodeElement(node));
};

// Past the hook's settle time, so the pending check actually goes out.
const settle = async () => {
  await act(async () => {
    vi.advanceTimersByTime(600);
  });
};

beforeEach(() => {
  vi.useFakeTimers();
  vi.clearAllMocks();
  ({ postJson } = mockApiCall());
  postJson.mockResolvedValue({ ok: true, message: '' });
});

afterEach(() => {
  cleanup();
  vi.useRealTimers();
});

describe('node parameter validation', () => {
  it('shows the message Hermes returned for a bad value', async () => {
    postJson.mockResolvedValue({
      ok: false,
      message: "Agent 'NotAnAgent' doesn't exists, choose one of: H2S",
    });
    renderNode({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'NotAnAgent' } },
    });
    await settle();
    expect(screen.getByText(/doesn't exists/)).toBeTruthy();
  });

  it('shows nothing when the values are fine', async () => {
    renderNode({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H2S' } },
    });
    await settle();
    expect(screen.queryByText(/doesn't exists/)).toBeNull();
  });

  it('sends the node type and its parameters', async () => {
    renderNode({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H2S', Calculator: 'AEGL10min' } },
    });
    await settle();
    expect(postJson).toHaveBeenCalledWith('/node-params/validate', {
      type: 'RiskAssessment.calculateThresholds',
      params: { Agent: 'H2S', Calculator: 'AEGL10min' },
    }, expect.objectContaining({ label: expect.any(String) }));
  });

  it('does not ask before the settle time passes', () => {
    renderNode({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H2S' } },
    });
    expect(postJson).not.toHaveBeenCalled();
  });

  it('sends one check for a run of edits, with the last value', async () => {
    const { rerender } = renderNode({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H' } },
    });
    rerender(nodeElement({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H2' } },
    }));
    rerender(nodeElement({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H2S' } },
    }));
    expect(postJson).not.toHaveBeenCalled();
    await settle();
    expect(postJson).toHaveBeenCalledTimes(1);
    expect(postJson).toHaveBeenLastCalledWith('/node-params/validate', {
      type: 'RiskAssessment.calculateThresholds',
      params: { Agent: 'H2S' },
    }, expect.objectContaining({ label: expect.any(String) }));
  });

  it('does not ask for a node with no type', async () => {
    renderNode({ type: '', Execution: { input_parameters: { Agent: 'H2S' } } });
    await settle();
    expect(postJson).not.toHaveBeenCalled();
  });
});
