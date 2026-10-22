import { describe, it, expect } from 'vitest';
import { Model, TabNode } from 'flexlayout-react';
import { LayoutModel } from '../src/components/layout/LayoutModel';
import { ProjectObj } from '../src/objects/ProjectObj';

const project = new ProjectObj({ name: 'P', documents: [] } as any);

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

describe('LayoutModel.closeCanvasTab', () => {
  it('closes the canvas of one document and leaves the others', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');
    layout.openOrFocusCanvasTab('doc2', 'Workflow2');

    layout.closeCanvasTab('doc1');

    expect(canvasTabs(layout.model).map(t => t.getId())).toEqual(['canvas:doc2']);
  });

  it('does nothing when the document has no canvas open', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusCanvasTab('doc1', 'Workflow1');

    layout.closeCanvasTab('doc2');

    expect(canvasTabs(layout.model)).toHaveLength(1);
  });
});

describe('LayoutModel.activeDetailsTab', () => {
  it('is undefined when no details tab is open', () => {
    const layout = LayoutModel.create(true);

    expect(layout.activeDetailsTab()).toBeUndefined();
  });

  it('is the selected details tab', () => {
    const layout = LayoutModel.create(true);
    layout.openOrFocusDetailsTab('doc:doc1', project);
    layout.openOrFocusDetailsTab('doc:doc2', project);

    expect(layout.activeDetailsTab()?.getId()).toBe('details:doc:doc2');
  });
});
