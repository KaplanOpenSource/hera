import { describe, it, expect } from 'vitest';
import { Model, TabNode } from 'flexlayout-react';
import { LayoutModel } from '../src/components/layout/LayoutModel';

const canvasTabs = (model: Model): TabNode[] => {
  const tabs: TabNode[] = [];
  model.visitNodes((node) => {
    if (node.getType() === 'tab' && node.getId().startsWith('canvas:')) {
      tabs.push(node as TabNode);
    }
  });
  return tabs;
};

describe('LayoutModel.openOrFocusCanvasTab', () => {
  it('opens one canvas tab below the details panel', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');

    expect(canvasTabs(layout.model)).toHaveLength(1);
  });

  it('reuses the tab when the same workflow is opened again', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');

    expect(canvasTabs(layout.model)).toHaveLength(1);
  });

  it('puts every canvas in the same tabset, so the details panel keeps its size', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');
    const rowsAfterFirst = (layout.model.toJson().layout.children as any[]).length;

    layout.openOrFocusCanvasTab('doc2', 'Workflow2');
    layout.openOrFocusCanvasTab('doc3', 'Workflow3');

    const tabs = canvasTabs(layout.model);
    expect(tabs).toHaveLength(3);

    const tabsetIds = new Set(tabs.map(t => t.getParent()!.getId()));
    expect(tabsetIds.size).toBe(1);

    // No extra row was added, so nothing halved the details panel.
    expect((layout.model.toJson().layout.children as any[]).length).toBe(rowsAfterFirst);
  });
});
