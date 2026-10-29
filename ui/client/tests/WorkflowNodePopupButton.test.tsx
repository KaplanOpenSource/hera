import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowNodePopupButton } from '../src/components/workflow/WorkflowNodePopupButton';

afterEach(() => cleanup());

describe('WorkflowNodePopupButton', () => {
  it('calls onOpen when clicked', () => {
    const onOpen = vi.fn();
    render(<WorkflowNodePopupButton onOpen={onOpen} />);
    fireEvent.click(screen.getByLabelText('open in a pop-up'));
    expect(onOpen).toHaveBeenCalledTimes(1);
  });

  it('stops the click from bubbling to the node behind it', () => {
    const onOpen = vi.fn();
    const onParentClick = vi.fn();
    render(
      <div onClick={onParentClick}>
        <WorkflowNodePopupButton onOpen={onOpen} />
      </div>,
    );
    fireEvent.click(screen.getByLabelText('open in a pop-up'));
    expect(onOpen).toHaveBeenCalled();
    expect(onParentClick).not.toHaveBeenCalled();
  });
});
