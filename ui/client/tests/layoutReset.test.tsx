import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { createLayoutModel } from '../src/components/layout/ProjectLayout';
import { DashboardHeader } from '../src/components/header/DashboardHeader';

vi.mock('../src/components/header/PageTitle', () => ({
  PageTitle: () => <div>PageTitle</div>,
}));
vi.mock('../src/components/header/ProjectChooser', () => ({
  ProjectChooser: () => <div>ProjectChooser</div>,
}));
vi.mock('../src/components/header/StatusIndicators', () => ({
  StatusIndicators: () => null,
}));

// Pull the tabset definitions out of the model JSON for assertions.
const tabsetsOf = (model: ReturnType<typeof createLayoutModel>) =>
  (model.toJson().layout.children as any[]);

describe('createLayoutModel (panel layout reset)', () => {
  it('builds the original two-panel arrangement: tree (25%) left, details (75%) right', () => {
    const tabsets = tabsetsOf(createLayoutModel(true));
    expect(tabsets).toHaveLength(2);

    const tree = tabsets.find(t => t.id === 'tree-tabset');
    const details = tabsets.find(t => t.id === 'details-tabset');
    expect(tree.weight).toBe(25);
    expect(details.weight).toBe(75);
    expect(tree.children).toHaveLength(1);
    expect(tree.children[0].id).toBe('tree');
  });

  it('omits the tree panel when the sidebar is collapsed', () => {
    const tabsets = tabsetsOf(createLayoutModel(false));
    expect(tabsets).toHaveLength(1);
    expect(tabsets[0].id).toBe('details-tabset');
  });

  it('keeps the supplied open details tabs in the details panel', () => {
    const openTabs = [
      { type: 'tab', id: 'details:doc-a', name: 'Doc A', component: 'details' },
      { type: 'tab', id: 'details:doc-b', name: 'Doc B', component: 'details' },
    ];
    const tabsets = tabsetsOf(createLayoutModel(true, openTabs));
    const details = tabsets.find(t => t.id === 'details-tabset');
    expect(details.children.map((c: any) => c.id)).toEqual(['details:doc-a', 'details:doc-b']);
  });

  it('marks the active tab as selected when a valid index is given', () => {
    const openTabs = [
      { type: 'tab', id: 'details:doc-a', name: 'Doc A', component: 'details' },
      { type: 'tab', id: 'details:doc-b', name: 'Doc B', component: 'details' },
    ];
    const details = tabsetsOf(createLayoutModel(true, openTabs, 1))
      .find(t => t.id === 'details-tabset');
    expect(details.selected).toBe(1);
  });

  it('leaves selection unset when no tab is active (index -1)', () => {
    const details = tabsetsOf(createLayoutModel(true, [], -1))
      .find(t => t.id === 'details-tabset');
    expect(details.selected ?? -1).toBe(-1);
  });
});

describe('DashboardHeader reset button', () => {
  afterEach(() => {
    cleanup();
  });

  it('calls onResetLayout when the reset panel layout button is clicked', () => {
    const onResetLayout = vi.fn();
    render(
      <DashboardHeader
        treeCollapsed={false}
        setTreeCollapsed={() => {}}
        onResetLayout={onResetLayout}
      />
    );
    const button = screen.getByTestId('ViewQuiltIcon').closest('button');
    fireEvent.click(button!);
    expect(onResetLayout).toHaveBeenCalledTimes(1);
  });
});
