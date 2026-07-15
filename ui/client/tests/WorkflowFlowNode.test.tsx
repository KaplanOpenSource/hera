import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';

// The node renders ReactFlow Handles, which need the ReactFlow store. We only
// care about the editing UI here, so stub them out.
vi.mock('@xyflow/react', () => ({
  Handle: () => null,
  NodeResizer: () => null,
  Position: { Left: 'left', Right: 'right', Top: 'top', Bottom: 'bottom' },
}));

import { WorkflowFlowNode } from '../src/components/workflow/WorkflowFlowNode';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';

afterEach(() => cleanup());

const catalog: NodeCatalogEntry[] = [
  {
    type: 'openFOAM.constant.g',
    parameters: [
      { name: 'x', is_required: true, source: 'jsonForm' },
      { name: 'y', is_required: true, source: 'jsonForm' },
      { name: 'z', is_required: true, source: 'jsonForm' },
    ],
  },
  { type: 'general.JinjaTransform', parameters: [] },
];

const renderNode = (node: any = {}) => {
  const onRename = vi.fn();
  const onChange = vi.fn();
  const onDelete = vi.fn();
  render(
    <WorkflowFlowNode
      data={{ name: 'node1', node, catalog, onRename, onChange, onDelete }}
      selected={false}
      // The rest of NodeProps is unused by the component.
      {...({} as any)}
    />,
  );
  return { onRename, onChange, onDelete };
};

describe('WorkflowFlowNode', () => {
  it('shows the current name and type', () => {
    renderNode({ type: 'general.JinjaTransform' });
    expect((screen.getByLabelText('node name') as HTMLInputElement).value).toBe('node1');
    expect((screen.getByLabelText('type') as HTMLInputElement).value).toBe('general.JinjaTransform');
  });

  it('updates only the type when the field is typed into', () => {
    const { onChange } = renderNode({ type: '', Execution: { input_parameters: {} } });
    fireEvent.change(screen.getByLabelText('type'), { target: { value: 'custom.type' } });
    expect(onChange).toHaveBeenCalledWith(expect.objectContaining({ type: 'custom.type' }));
    const last = onChange.mock.calls.at(-1)![0];
    expect(last.Execution.input_parameters).toEqual({});
  });

  it('prefills the parameter names when a catalog type is picked', async () => {
    const { onChange } = renderNode({ type: '', Execution: { input_parameters: {} } });
    const input = screen.getByLabelText('type');
    fireEvent.change(input, { target: { value: 'openFOAM.constant.g' } });
    const option = await screen.findByRole('option', { name: /openFOAM\.constant\.g/ });
    fireEvent.click(option);
    const last = onChange.mock.calls.at(-1)![0];
    expect(last.type).toBe('openFOAM.constant.g');
    expect(last.Execution.input_parameters).toEqual({ x: '', y: '', z: '' });
  });

  it('keeps existing parameter values when re-picking a type', async () => {
    const { onChange } = renderNode({ type: '', Execution: { input_parameters: { x: 5 } } });
    const input = screen.getByLabelText('type');
    fireEvent.change(input, { target: { value: 'openFOAM.constant.g' } });
    const option = await screen.findByRole('option', { name: /openFOAM\.constant\.g/ });
    fireEvent.click(option);
    const last = onChange.mock.calls.at(-1)![0];
    expect(last.Execution.input_parameters).toEqual({ x: 5, y: '', z: '' });
  });

  it('renames the node on blur of the name field', () => {
    const { onRename } = renderNode({ type: 't' });
    const nameInput = screen.getByLabelText('node name');
    fireEvent.change(nameInput, { target: { value: 'renamed' } });
    fireEvent.blur(nameInput);
    expect(onRename).toHaveBeenCalledWith('renamed');
  });

  it('does not rename when the name is unchanged', () => {
    const { onRename } = renderNode({ type: 't' });
    fireEvent.blur(screen.getByLabelText('node name'));
    expect(onRename).not.toHaveBeenCalled();
  });

  it('warns when the type is not in the catalog', () => {
    renderNode({ type: 'made.up' });
    expect(screen.getByText('unknown type: made.up')).toBeDefined();
  });

  it('shows "required" under empty required fields', () => {
    // x and z empty, y filled -> two "required" helper texts.
    renderNode({ type: 'openFOAM.constant.g', Execution: { input_parameters: { x: '', y: 2, z: '' } } });
    expect(screen.getAllByText('required')).toHaveLength(2);
  });

  // Regression for #990: the auto-reload re-renders the editor with fresh
  // (equal) objects every few seconds; the focused input must not remount and
  // lose focus while the user is typing.
  it('keeps input focus across a data reload (re-render with fresh objects)', () => {
    const build = () => (
      <WorkflowFlowNode
        data={{
          name: 'node1',
          node: { type: 'general.JinjaTransform', Execution: { input_parameters: { x: 'aaa' } } },
          catalog,
          onRename: vi.fn(),
          onChange: vi.fn(),
          onDelete: vi.fn(),
          onFieldContextMenu: vi.fn(),
          onFieldInlineEdit: vi.fn(),
        }}
        selected={false}
        {...({} as any)}
      />
    );
    const { rerender } = render(build());
    const input = screen.getByDisplayValue('aaa') as HTMLInputElement;
    input.focus();
    expect(document.activeElement).toBe(input);

    // A fresh but equal document arrives (as the periodic reload delivers).
    rerender(build());

    expect(document.activeElement).toBe(input);
    expect(screen.getByDisplayValue('aaa')).toBe(input);
  });
});
