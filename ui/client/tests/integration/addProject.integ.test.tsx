/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { renderApp, resetStore } from './integHelpers';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock());
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());
vi.mock('../../src/io/snackbar', async () => (await import('./mockFactories')).createSnackbarMock());

import { useProjectStore } from '../../src/stores/useProjectStore';

describe('Add Project UI integration', () => {
  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });
  afterAll(() => { cleanup(); });

  it('clicking Add Project creates the project and it appears in the list', async () => {
    renderApp('/');

    const addWrapper = await screen.findByLabelText('Add project');
    await act(async () => { fireEvent.click(within(addWrapper).getByRole('button')); });

    const dialog = await screen.findByRole('dialog', {}, { timeout: 3000 });
    fireEvent.change(within(dialog).getByRole('textbox', { name: /project name/i }), {
      target: { value: 'UITestProject' },
    });

    const loadRepos = within(dialog).getByRole('checkbox', { name: /load repositories/i });
    if ((loadRepos as HTMLInputElement).checked) fireEvent.click(loadRepos);

    await act(async () => {
      fireEvent.click(within(dialog).getByText('Add Project'));
    });

    await waitFor(() => {
      const names = useProjectStore.getState().projectNames;
      expect(names.some(p => p.name === 'UITestProject')).toBe(true);
    }, { timeout: 15000 });
  }, 30000);

  it('project persists after simulated refresh', async () => {
    resetStore();
    renderApp('/');

    await waitFor(() => {
      const names = useProjectStore.getState().projectNames;
      expect(names.some(p => p.name === 'UITestProject')).toBe(true);
    }, { timeout: 10000 });
  }, 15000);
}, 90000);
