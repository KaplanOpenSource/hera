import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { act, cleanup, render, screen } from '@testing-library/react';

// The node renders ReactFlow Handles, which need the ReactFlow store.
vi.mock('@xyflow/react', () => ({
  Handle: () => null,
  NodeResizer: () => null,
  Position: { Left: 'left', Right: 'right', Top: 'top', Bottom: 'bottom' },
}));

const validateNodeParams = vi.fn();
vi.mock('../src/io/validateNodeParams', () => ({
  validateNodeParams: (args: any) => validateNodeParams(args),
}));

import { WorkflowFlowNode } from '../src/components/workflow/WorkflowFlowNode';

const catalog = [{ type: 'RiskAssessment.calculateThresholds', parameters: [] }];

const renderNode = (node: any) => {
  render(
    <WorkflowFlowNode
      data={{ name: 'node1', node, catalog, onRename: vi.fn(), onChange: vi.fn(), onDelete: vi.fn() }}
      selected={false}
      {...({} as any)}
    />,
  );
};

// Past the hook's settle time, so the pending check actually goes out.
const settle = async () => {
  await act(async () => {
    vi.advanceTimersByTime(600);
  });
};

beforeEach(() => {
  vi.useFakeTimers();
  validateNodeParams.mockReset();
  validateNodeParams.mockResolvedValue({ ok: true, message: '' });
});

afterEach(() => {
  cleanup();
  vi.useRealTimers();
});

describe('node parameter validation', () => {
  it('shows the message Hermes returned for a bad value', async () => {
    validateNodeParams.mockResolvedValue({
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
    expect(validateNodeParams).toHaveBeenCalledWith({
      type: 'RiskAssessment.calculateThresholds',
      params: { Agent: 'H2S', Calculator: 'AEGL10min' },
    });
  });

  it('does not ask before the settle time passes', () => {
    renderNode({
      type: 'RiskAssessment.calculateThresholds',
      Execution: { input_parameters: { Agent: 'H2S' } },
    });
    expect(validateNodeParams).not.toHaveBeenCalled();
  });

  it('does not ask for a node with no type', async () => {
    renderNode({ type: '', Execution: { input_parameters: { Agent: 'H2S' } } });
    await settle();
    expect(validateNodeParams).not.toHaveBeenCalled();
  });
});
