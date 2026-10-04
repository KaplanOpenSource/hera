import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowNodePinButton } from '../src/components/workflow/WorkflowNodePinButton';

afterEach(() => cleanup());

describe('WorkflowNodePinButton', () => {
  it('calls onToggle when clicked', () => {
    const onToggle = vi.fn();
    render(<WorkflowNodePinButton pinned={false} onToggle={onToggle} />);
    fireEvent.click(screen.getByLabelText('keep open'));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });

  it('stops the click from bubbling to the node behind it', () => {
    const onToggle = vi.fn();
    const onParentClick = vi.fn();
    render(
      <div onClick={onParentClick}>
        <WorkflowNodePinButton pinned={false} onToggle={onToggle} />
      </div>,
    );
    fireEvent.click(screen.getByLabelText('keep open'));
    expect(onToggle).toHaveBeenCalled();
    expect(onParentClick).not.toHaveBeenCalled();
  });

  it('fills the pin in once the node is pinned', () => {
    const { rerender } = render(<WorkflowNodePinButton pinned={false} onToggle={vi.fn()} />);
    expect(document.querySelector('[data-testid="PushPinOutlinedIcon"]')).not.toBeNull();
    rerender(<WorkflowNodePinButton pinned onToggle={vi.fn()} />);
    expect(document.querySelector('[data-testid="PushPinIcon"]')).not.toBeNull();
  });
});
