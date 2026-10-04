import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, render } from '@testing-library/react';

// Handles need the ReactFlow store; render them as plain markers instead.
vi.mock('@xyflow/react', () => ({
  Handle: ({ type, id, children }: any) => <div data-handle-type={type} data-handle-id={id}>{children}</div>,
  Position: { Left: 'left', Right: 'right' },
}));

import { WorkflowNodeInputs } from '../src/components/workflow/WorkflowNodeInputs';

afterEach(() => cleanup());

const handleIds = (type: string): string[] => {
  return Array.from(document.querySelectorAll(`[data-handle-type="${type}"]`))
    .map(el => el.getAttribute('data-handle-id') as string);
};

const sourceHandleIds = (): string[] => {
  return handleIds('source');
};

const targetHandleIds = (): string[] => {
  return handleIds('target');
};

const renderInputs = () => {
  render(
    <WorkflowNodeInputs
      nodeName="A"
      params={{ cmd: 'ls', group: { inner: 'x' }, list: ['a'] }}
      paramsDef={{}}
      expandedItems={['input_parameters', 'input_parameters/group', 'input_parameters/list']}
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
    expect(sourceHandleIds()).toEqual(['A:param:cmd', 'A:param:group', 'A:param:list']);
  });

  it('puts a target handle on a key inside a dict, named by its path', () => {
    renderInputs();
    expect(targetHandleIds()).toContain('A:in:group.inner');
  });

  it('puts a target handle on a list element, named by its position', () => {
    renderInputs();
    expect(targetHandleIds()).toContain('A:in:list.0');
  });

  it('puts no target handle on a dict or list row, which holds no value', () => {
    renderInputs();
    expect(targetHandleIds()).toEqual(['A:in:cmd', 'A:in:group.inner', 'A:in:list.0']);
  });
});
