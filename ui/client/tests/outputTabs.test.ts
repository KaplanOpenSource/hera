import { describe, it, expect } from 'vitest';
import { Model, TabNode } from 'flexlayout-react';
import { LayoutModel } from '../src/components/layout/LayoutModel';

const outputTabs = (model: Model): TabNode[] => {
  const tabs: TabNode[] = [];
  model.visitNodes((node) => {
    if (node.getType() === 'tab' && node.getId().startsWith('output:')) {
      tabs.push(node as TabNode);
    }
  });
  return tabs;
};

describe('LayoutModel.openOrFocusOutputTab', () => {
  it('opens one output tab', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusOutputTab('hello_1');

    const tabs = outputTabs(layout.model);
    expect(tabs).toHaveLength(1);
    expect(tabs[0].getConfig().workflowName).toBe('hello_1');
  });

  it('sits to the right of the canvas, in the canvas row', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');
    layout.openOrFocusOutputTab('Workflow1');

    const canvasTabset = layout.model.getNodeById('canvas:doc1')!.getParent()!;
    const outputTabset = outputTabs(layout.model)[0].getParent()!;
    expect(outputTabset.getId()).not.toBe(canvasTabset.getId());
    // Same row, so they sit side by side rather than one above the other.
    expect(outputTabset.getParent()!.getId()).toBe(canvasTabset.getParent()!.getId());

    const row = outputTabset.getParent()!;
    const ids = row.getChildren().map(c => c.getId());
    expect(ids.indexOf(outputTabset.getId())).toBeGreaterThan(ids.indexOf(canvasTabset.getId()));
  });

  it('reuses the tab when the same workflow runs again', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusOutputTab('hello_1');
    layout.openOrFocusOutputTab('hello_1');

    expect(outputTabs(layout.model)).toHaveLength(1);
  });

  it('puts every output in the same tabset, so the details panel keeps its size', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusOutputTab('w1');
    const rowsAfterFirst = (layout.model.toJson().layout.children as any[]).length;

    layout.openOrFocusOutputTab('w2');
    layout.openOrFocusOutputTab('w3');

    const tabs = outputTabs(layout.model);
    expect(tabs).toHaveLength(3);

    const tabsetIds = new Set(tabs.map(t => t.getParent()!.getId()));
    expect(tabsetIds.size).toBe(1);

    // No extra row was added, so nothing halved the details panel.
    expect((layout.model.toJson().layout.children as any[]).length).toBe(rowsAfterFirst);
  });

  it('opens the tab again after it was closed', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusOutputTab('w1');
    layout.model.doAction({ type: 'FlexLayout_DeleteTab', data: { node: 'output:w1' } } as any);
    expect(outputTabs(layout.model)).toHaveLength(0);

    layout.openOrFocusOutputTab('w1');
    expect(outputTabs(layout.model)).toHaveLength(1);
  });
});
