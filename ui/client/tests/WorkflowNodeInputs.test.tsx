import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, render } from '@testing-library/react';

// Handles need the ReactFlow store; render them as plain markers instead.
vi.mock('@xyflow/react', () => ({
  Handle: ({ type, id, children }: any) => <div data-handle-type={type} data-handle-id={id}>{children}</div>,
  Position: { Left: 'left', Right: 'right' },
}));

import { WorkflowNodeInputs } from '../src/components/workflow/WorkflowNodeInputs';

afterEach(() => cleanup());

const sourceHandleIds = (): string[] => {
  return Array.from(document.querySelectorAll('[data-handle-type="source"]'))
    .map(el => el.getAttribute('data-handle-id') as string);
};

const renderInputs = () => {
  render(
    <WorkflowNodeInputs
      nodeName="A"
      params={{ cmd: 'ls', group: { inner: 'x' } }}
      paramsDef={{}}
      expandedItems={['input_parameters', 'input_parameters/group']}
      onExpandedItemsChange={vi.fn()}
      onChangeParams={vi.fn()}
      onFieldContextMenu={vi.fn()}
      onFieldInlineEdit={vi.fn()}
    />,
  );
};

describe('WorkflowNodeInputs', () => {
  it('puts a source handle on a top-level parameter row', () => {
    renderInputs();
    expect(sourceHandleIds()).toContain('A:param:cmd');
  });

  it('puts no source handle on a nested row', () => {
    renderInputs();
    expect(sourceHandleIds()).toEqual(['A:param:cmd', 'A:param:group']);
  });
});
