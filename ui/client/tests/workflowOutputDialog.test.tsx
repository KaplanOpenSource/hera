import { describe, it, expect, afterEach } from 'vitest';
import { render, screen, cleanup } from '@testing-library/react';
import { WorkflowOutputDialog } from '../src/components/workflow/log/WorkflowOutputDialog';

afterEach(() => {
  cleanup();
});

describe('WorkflowOutputDialog', () => {
  it('shows the partial log and a running hint while running', () => {
    render(
      <WorkflowOutputDialog
        open
        running
        output={'partial line while running'}
        error={null}
        workflowName="w"
        onClose={() => {}}
      />,
    );

    expect(screen.getByText('Running…')).toBeTruthy();
    expect(screen.getByText('partial line while running')).toBeTruthy();
  });

  it('shows the final log when done', () => {
    render(
      <WorkflowOutputDialog
        open
        running={false}
        output={'final output line'}
        error={null}
        workflowName="w"
        onClose={() => {}}
      />,
    );

    expect(screen.queryByText('Running…')).toBeNull();
    expect(screen.getByText('final output line')).toBeTruthy();
  });

  it('shows the error when the run failed', () => {
    render(
      <WorkflowOutputDialog
        open
        running={false}
        output={null}
        error={'it broke'}
        workflowName="w"
        onClose={() => {}}
      />,
    );

    expect(screen.getByText('it broke')).toBeTruthy();
  });
});
