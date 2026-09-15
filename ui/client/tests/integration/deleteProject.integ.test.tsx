/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { createProjectViaUI, openProjectActions, renderApp, resetStore } from './integHelpers';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock());
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());
vi.mock('../../src/io/snackbar', async () => (await import('./mockFactories')).createSnackbarMock());

import { useProjectStore } from '../../src/stores/useProjectStore';

describe('Delete Project UI integration', () => {
  beforeAll(async () => {
    await createProjectViaUI('ProjectToKeep');
    await createProjectViaUI('ProjectToDelete');
  }, 90000);

  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });
  afterAll(() => { cleanup(); });

  it('both projects appear in the list', async () => {
    renderApp('/ProjectToDelete');

    await waitFor(() => {
      const names = useProjectStore.getState().projectNames.map(p => p.name);
      expect(names).toContain('ProjectToKeep');
      expect(names).toContain('ProjectToDelete');
    }, { timeout: 15000 });
  }, 30000);

  it('deleting a project removes it from the list and switches to another', async () => {
    renderApp('/ProjectToDelete');

    await waitFor(() => {
      expect(useProjectStore.getState().currProjectName).toBe('ProjectToDelete');
    }, { timeout: 10000 });

    await openProjectActions();
    const deleteWrapper = await screen.findByLabelText('Delete project');
    await act(async () => {
      fireEvent.click(within(deleteWrapper).getByRole('button'));
    });

    const input = await screen.findByRole('textbox');
    fireEvent.change(input, { target: { value: 'ProjectToDelete' } });

    const yesBtn = screen.getByRole('button', { name: /yes/i });
    await act(async () => { fireEvent.click(yesBtn); });

    await waitFor(() => {
      const names = useProjectStore.getState().projectNames.map(p => p.name);
      expect(names).not.toContain('ProjectToDelete');
      expect(names).toContain('ProjectToKeep');
    }, { timeout: 15000 });

    await waitFor(() => {
      const curr = useProjectStore.getState().currProjectName;
      expect(curr).not.toBe('ProjectToDelete');
      expect(curr).not.toBe('* NONE *');
    }, { timeout: 5000 });
  }, 45000);

  it('deleted project does not reappear after simulated refresh', async () => {
    resetStore();
    renderApp('/ProjectToKeep');

    await waitFor(() => {
      const names = useProjectStore.getState().projectNames.map(p => p.name);
      expect(names).toContain('ProjectToKeep');
      expect(names).not.toContain('ProjectToDelete');
    }, { timeout: 15000 });
  }, 30000);
}, 180000);
