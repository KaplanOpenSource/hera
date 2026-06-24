import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowNodeDeleteButton } from '../src/components/workflow/WorkflowNodeDeleteButton';

afterEach(() => cleanup());

describe('WorkflowNodeDeleteButton', () => {
  it('calls onDelete when clicked', () => {
    const onDelete = vi.fn();
    render(<WorkflowNodeDeleteButton onDelete={onDelete} />);
    fireEvent.click(screen.getByRole('button'));
    expect(onDelete).toHaveBeenCalledTimes(1);
  });

  it('stops the click from bubbling to the node behind it', () => {
    const onDelete = vi.fn();
    const onParentClick = vi.fn();
    render(
      <div onClick={onParentClick}>
        <WorkflowNodeDeleteButton onDelete={onDelete} />
      </div>,
    );
    fireEvent.click(screen.getByRole('button'));
    expect(onDelete).toHaveBeenCalled();
    expect(onParentClick).not.toHaveBeenCalled();
  });
});
