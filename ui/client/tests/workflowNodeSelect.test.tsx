import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowNodeSelect } from '../src/components/workflow/log/WorkflowNodeSelect';

const chunks = [
  { name: '__preamble__', text: 'starting up\n' },
  { name: 'ListFiles', text: 'ran ls\n' },
  { name: 'Empty', text: '   \n' },
  { name: '__between__', text: 'between tasks\n' },
  { name: 'finalnode_xx', text: 'done\n' },
];

afterEach(() => {
  cleanup();
});

const openMenu = () => {
  fireEvent.mouseDown(screen.getByRole('combobox'));
};

describe('WorkflowNodeSelect', () => {
  it('renders nothing when no chunk has content', () => {
    const { container } = render(
      <WorkflowNodeSelect chunks={[{ name: 'A', text: '  ' }]} onPick={() => {}} />,
    );

    expect(container.firstChild).toBeNull();
  });

  it('lists only the task chunks that render a card', () => {
    render(<WorkflowNodeSelect chunks={chunks} onPick={() => {}} />);
    openMenu();

    const labels = screen.getAllByRole('option').map(o => o.textContent);
    // The setup / between filler is left out, and "Empty" renders no card.
    expect(labels).toEqual(['ListFiles', 'finalnode_xx']);
  });

  it('shows the chunk currently in view', () => {
    render(<WorkflowNodeSelect chunks={chunks} currentIndex={1} onPick={() => {}} />);

    expect(screen.getByRole('combobox').textContent).toBe('ListFiles');
  });

  it('shows nothing selected when the current chunk has no card', () => {
    const { container } = render(<WorkflowNodeSelect chunks={chunks} currentIndex={2} onPick={() => {}} />);

    expect(container.querySelector('input')?.value).toBe('');
  });

  it('renders nothing when the run has only filler chunks', () => {
    const { container } = render(
      <WorkflowNodeSelect chunks={[{ name: '__preamble__', text: 'starting up' }]} onPick={() => {}} />,
    );

    expect(container.firstChild).toBeNull();
  });

  it('reports the picked chunk by its index in the full chunk list', () => {
    const onPick = vi.fn();
    render(<WorkflowNodeSelect chunks={chunks} onPick={onPick} />);
    openMenu();

    fireEvent.click(screen.getByRole('option', { name: 'finalnode_xx' }));

    // Index 4: the filler and empty chunks still count in the chunk list.
    expect(onPick).toHaveBeenCalledWith(4);
  });
});
