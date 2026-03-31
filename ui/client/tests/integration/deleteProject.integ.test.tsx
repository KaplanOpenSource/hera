/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { startDockerEnv, type DockerEnv } from './dockerSetup';
import { createProjectViaUI, renderApp, resetStore } from './integHelpers';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock(8004));
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());

import { useProjectStore } from '../../src/stores/useProjectStore';

let env: DockerEnv;

describe('Delete Project UI integration', () => {
  beforeAll(async () => {
    env = await startDockerEnv({
      network: 'hera-test-delproj-net',
      mongoContainer: 'hera-test-delproj-mongo',
      serverContainer: 'hera-test-delproj-server',
      serverPort: 8004,
      dbName: 'hera_test_delproj',
    });
    await createProjectViaUI('ProjectToKeep');
    await createProjectViaUI('ProjectToDelete');
  }, 90000);

  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });

  afterAll(() => {
    cleanup();
    env?.cleanup();
  }, 15000);

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

    const deleteWrapper = screen.getByLabelText('Delete project');
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
      expect(useProjectStore.getState().currProjectName).toBe('ProjectToKeep');
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
