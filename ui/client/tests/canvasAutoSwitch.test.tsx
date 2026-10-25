import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { ProjectLayout } from '../src/components/layout/ProjectLayout';
import { ProjectObj } from '../src/objects/ProjectObj';
import { idDocId } from '../src/shared/idDocId';
import { ProjectDocument } from '../src/shared/types';

// The panels are not under test; the tree stub only has to report a selection.
vi.mock('../src/components/layout/LayoutPanel', () => ({
  LayoutPanel: ({ component, onSelectItem }: any) => {
    if (component !== 'tree') {
      return null;
    }
    return (
      <>
        <button data-testid="pick-a" onClick={() => onSelectItem(idDocId('a'))}>pick a</button>
        <button data-testid="pick-b" onClick={() => onSelectItem(idDocId('b'))}>pick b</button>
      </>
    );
  },
}));

// flexlayout only draws a tab's contents once it knows its size, which jsdom does
// not give it: report a fixed size and tell the observer about it at once.
const LAYOUT_RECT = { x: 0, y: 0, width: 1200, height: 800, top: 0, left: 0, right: 1200, bottom: 800, toJSON: () => {} } as DOMRect;

beforeAll(() => {
  vi.spyOn(Element.prototype, 'getBoundingClientRect').mockReturnValue(LAYOUT_RECT);
  globalThis.ResizeObserver = class {
    constructor(private callback: ResizeObserverCallback) {}
    observe(target: Element) {
      this.callback([{ target, contentRect: LAYOUT_RECT } as ResizeObserverEntry], this as unknown as ResizeObserver);
    }
    unobserve() {}
    disconnect() {}
  } as unknown as typeof ResizeObserver;
});

const workflowDoc = (oid: string, name: string): ProjectDocument => ({
  _cls: 'Metadata.Simulations',
  _id: { $oid: oid },
  projectName: 'TestProject',
  desc: { datasourceName: name, workflow: { nodeList: [], nodes: {} } },
  type: 'hermesWorkflow',
  resource: '',
  dataFormat: 'JSON_dict',
});

// A fresh ProjectObj over the same data, as an auto-reload produces.
const project = () => {
  return new ProjectObj({
    name: 'TestProject',
    documents: [workflowDoc('a', 'wfA'), workflowDoc('b', 'wfB')],
  });
};

// The tab strip's button for one canvas. flexlayout also keeps an offscreen copy
// of each button for dragging, so the stamps are skipped.
const canvasTabButton = (docName: string): Element => {
  const buttons = screen.getAllByText(`Canvas: ${docName}`)
    .map(label => label.closest('.flexlayout__tab_button'))
    .filter(button => button && !button.className.includes('stamp'));
  expect(buttons).toHaveLength(1);
  return buttons[0]!;
};

const isSelected = (button: Element): boolean => {
  return button.className.includes('--selected');
};

describe('two open workflow canvases', () => {
  afterEach(() => {
    cleanup();
  });

  // Issue 1119: the canvas of the active document was re-selected on every
  // render, so an auto-reload pulled the view off the canvas the user picked.
  it('keeps the canvas the user picked when the project reloads', () => {
    const { rerender } = render(
      <MemoryRouter>
        <ProjectLayout project={project()} treeCollapsed={false} resetSignal={0} />
      </MemoryRouter>
    );

    // Open both canvases, leaving wfA the active document.
    fireEvent.click(screen.getByTestId('pick-b'));
    fireEvent.click(screen.getByTestId('pick-a'));
    expect(isSelected(canvasTabButton('wfA'))).toBe(true);

    // The user looks at the other canvas.
    fireEvent.click(canvasTabButton('wfB'));
    expect(isSelected(canvasTabButton('wfB'))).toBe(true);

    // An auto-reload hands down a new project object holding the same documents.
    rerender(
      <MemoryRouter>
        <ProjectLayout project={project()} treeCollapsed={false} resetSignal={0} />
      </MemoryRouter>
    );

    expect(isSelected(canvasTabButton('wfB'))).toBe(true);
    expect(isSelected(canvasTabButton('wfA'))).toBe(false);
  });
});
