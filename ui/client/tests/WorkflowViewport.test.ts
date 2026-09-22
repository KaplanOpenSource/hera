import { describe, it, expect, vi } from 'vitest';
import { FIT_MAX_ZOOM, WorkflowViewport } from '../src/components/workflow/WorkflowViewport';

const fakeFlow = (node?: { position: { x: number, y: number }, measured?: { width?: number, height?: number } }) => {
  return {
    fitView: vi.fn(),
    getViewport: vi.fn(() => { return { x: 10, y: 20, zoom: 2 }; }),
    setViewport: vi.fn(),
    setCenter: vi.fn(),
    getNode: vi.fn(() => { return node; }),
  };
};

describe('WorkflowViewport', () => {
  it('fits the whole graph with the zoom cap', () => {
    const flow = fakeFlow();
    new WorkflowViewport(flow).fitAll();
    expect(flow.fitView).toHaveBeenCalledWith({ duration: 300, maxZoom: FIT_MAX_ZOOM });
  });

  it('centres on a node, keeping the zoom', () => {
    const flow = fakeFlow({ position: { x: 10, y: 20 }, measured: { width: 100, height: 40 } });
    new WorkflowViewport(flow).centerOnNode('A');
    expect(flow.setCenter).toHaveBeenCalledWith(60, 40, { zoom: 2, duration: 300 });
  });

  it('does nothing for a node the canvas does not know', () => {
    const flow = fakeFlow(undefined);
    new WorkflowViewport(flow).centerOnNode('gone');
    expect(flow.setCenter).not.toHaveBeenCalled();
  });

  it('scales the whole view by the ratio', () => {
    const flow = fakeFlow();
    new WorkflowViewport(flow).scaleZoom(0.5);
    expect(flow.setViewport).toHaveBeenCalledWith({ x: 5, y: 10, zoom: 1 });
  });
});
