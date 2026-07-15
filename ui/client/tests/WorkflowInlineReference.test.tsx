import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { useState } from 'react';
import { WorkflowInlineReference } from '../src/components/workflow/WorkflowInlineReference';

afterEach(() => cleanup());

// Renders the popup anchored to a real input so its native key listener (bound to
// the anchor element) has something to attach to.
const Harness = ({
  options,
  onPick,
  onClose,
}: {
  options: string[],
  onPick: (option: string) => void,
  onClose: () => void,
}) => {
  const [el, setEl] = useState<HTMLInputElement | null>(null);
  return (
    <>
      <input ref={setEl} data-testid="field" />
      <WorkflowInlineReference anchorEl={el} options={options} onPick={onPick} onClose={onClose} />
    </>
  );
};

describe('WorkflowInlineReference', () => {
  it('shows the options and picks the one clicked', () => {
    const onPick = vi.fn();
    render(<Harness options={['A', 'B', 'C']} onPick={onPick} onClose={vi.fn()} />);
    fireEvent.click(screen.getByText('B'));
    expect(onPick).toHaveBeenCalledWith('B');
  });

  it('renders nothing when there are no options', () => {
    render(<Harness options={[]} onPick={vi.fn()} onClose={vi.fn()} />);
    expect(screen.queryByRole('menuitem')).toBeNull();
  });

  it('navigates with arrows and picks the highlighted option on Enter', () => {
    const onPick = vi.fn();
    render(<Harness options={['A', 'B', 'C']} onPick={onPick} onClose={vi.fn()} />);
    const field = screen.getByTestId('field');
    fireEvent.keyDown(field, { key: 'ArrowDown' });
    fireEvent.keyDown(field, { key: 'ArrowDown' });
    fireEvent.keyDown(field, { key: 'Enter' });
    expect(onPick).toHaveBeenCalledWith('C');
  });

  it('wraps the highlight around the ends', () => {
    const onPick = vi.fn();
    render(<Harness options={['A', 'B']} onPick={onPick} onClose={vi.fn()} />);
    const field = screen.getByTestId('field');
    fireEvent.keyDown(field, { key: 'ArrowUp' });
    fireEvent.keyDown(field, { key: 'Enter' });
    expect(onPick).toHaveBeenCalledWith('B');
  });

  it('closes on Escape', () => {
    const onClose = vi.fn();
    render(<Harness options={['A']} onPick={vi.fn()} onClose={onClose} />);
    fireEvent.keyDown(screen.getByTestId('field'), { key: 'Escape' });
    expect(onClose).toHaveBeenCalled();
  });
});
