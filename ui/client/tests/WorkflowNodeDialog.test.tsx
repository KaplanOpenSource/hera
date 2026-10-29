import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';

// A handle outside a canvas node is an error, so each one leaves a marker to count.
vi.mock('@xyflow/react', () => ({
  Handle: () => <span data-testid="handle" />,
  Position: { Left: 'left', Right: 'right', Top: 'top', Bottom: 'bottom' },
}));

import { WorkflowNodeDialog } from '../src/components/workflow/WorkflowNodeDialog';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { NodeParameterSource } from '../src/shared/types';

afterEach(() => cleanup());

const catalog: NodeCatalogEntry[] = [
  {
    type: 'openFOAM.constant.g',
    parameters: [{ name: 'x', is_required: true, source: NodeParameterSource.JsonForm }],
    outputs: [{ name: 'result', source: NodeParameterSource.Python }],
  },
];

const renderDialog = () => {
  const onClose = vi.fn();
  const onRename = vi.fn();
  const onChange = vi.fn();
  render(
    <WorkflowNodeDialog
      name="node1"
      node={{ type: 'openFOAM.constant.g', Execution: { input_parameters: { x: 'aaa' } } }}
      catalog={catalog}
      onClose={onClose}
      onRename={onRename}
      onChange={onChange}
      onFieldContextMenu={vi.fn()}
      onFieldInlineEdit={vi.fn()}
    />,
  );
  return { onClose, onRename, onChange };
};

describe('WorkflowNodeDialog', () => {
  it('shows the node editor', () => {
    renderDialog();
    expect((screen.getByLabelText('node name') as HTMLInputElement).value).toBe('node1');
    expect((screen.getByLabelText('type') as HTMLInputElement).value).toBe('openFOAM.constant.g');
    expect(screen.getByDisplayValue('aaa')).toBeDefined();
  });

  it('draws no line handles', () => {
    renderDialog();
    expect(screen.queryAllByTestId('handle')).toHaveLength(0);
  });

  it('sends edits back', () => {
    const { onChange, onRename } = renderDialog();
    fireEvent.change(screen.getByLabelText('type'), { target: { value: 'custom.type' } });
    expect(onChange).toHaveBeenCalledWith(expect.objectContaining({ type: 'custom.type' }));

    const nameInput = screen.getByLabelText('node name');
    fireEvent.change(nameInput, { target: { value: 'renamed' } });
    fireEvent.blur(nameInput);
    expect(onRename).toHaveBeenCalledWith('renamed');
  });

  it('closes from its X', () => {
    const { onClose } = renderDialog();
    fireEvent.click(screen.getByLabelText('close'));
    expect(onClose).toHaveBeenCalled();
  });
});
