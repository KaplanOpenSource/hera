/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { startDockerEnv, type DockerEnv } from './dockerSetup';
import { createProjectViaUI, renderApp, resetStore } from './integHelpers';
import type { ProjectDocument } from '../../src/shared/types';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock(8003));
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());

import { useProjectStore } from '../../src/stores/useProjectStore';

let env: DockerEnv;

describe('Documents UI integration', () => {
  beforeAll(async () => {
    env = await startDockerEnv({
      network: 'hera-test-docs-net',
      mongoContainer: 'hera-test-docs-mongo',
      serverContainer: 'hera-test-docs-server',
      serverPort: 8003,
      dbName: 'hera_test_docs',
    });
    await createProjectViaUI('DocTestProject');
  }, 60000);

  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });

  afterAll(() => {
    cleanup();
    env.cleanup();
  }, 15000);

  const loadProject = async () => {
    renderApp('/DocTestProject');
    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('DocTestProject');
    }, { timeout: 15000 });
  };

  it('add a regular document via the UI and verify it appears', async () => {
    await loadProject();

    const addWrapper = await screen.findByLabelText('Add Document');
    await act(async () => { fireEvent.click(within(addWrapper).getByRole('button')); });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), {
      target: { value: 'IntegDoc1' },
    });
    fireEvent.change(within(dialog).getByRole('textbox', { name: /resource/i }), {
      target: { value: 'some-resource' },
    });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const docs: ProjectDocument[] = useProjectStore.getState().currProject?.documents ?? [];
      expect(docs.some(d => d.desc.datasourceName === 'IntegDoc1')).toBe(true);
    }, { timeout: 15000 });
  }, 30000);

  it('add an agent document via the UI', async () => {
    await loadProject();

    const addWrapper = await screen.findByLabelText('Add Document');
    await act(async () => { fireEvent.click(within(addWrapper).getByRole('button')); });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), {
      target: { value: 'AgentDoc1' },
    });
    fireEvent.click(within(dialog).getByRole('checkbox', { name: /agent/i }));

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const docs: ProjectDocument[] = useProjectStore.getState().currProject?.documents ?? [];
      const agentDoc = docs.find(d => d.desc.datasourceName === 'AgentDoc1');
      expect(agentDoc).toBeTruthy();
      expect(agentDoc!.type).toBe('ToolkitDataSource');
      expect(typeof agentDoc!.resource).toBe('object');
    }, { timeout: 15000 });
  }, 30000);

  it('delete a document via the tree view UI', async () => {
    await loadProject();

    const docLabel = await screen.findByText('IntegDoc1', {}, { timeout: 10000 });

    const docRow = docLabel.closest('[class*="MuiStack-root"]')!;
    const deleteWrapper = within(docRow as HTMLElement).getByLabelText(/delete document/i);
    await act(async () => {
      fireEvent.click(within(deleteWrapper).getByRole('button'));
    });

    const yesBtn = await screen.findByRole('button', { name: /yes/i });
    await act(async () => { fireEvent.click(yesBtn); });

    await waitFor(() => {
      const docs: ProjectDocument[] = useProjectStore.getState().currProject?.documents ?? [];
      expect(docs.some(d => d.desc.datasourceName === 'IntegDoc1')).toBe(false);
    }, { timeout: 15000 });
  }, 30000);

  it('document persists after simulated refresh', async () => {
    resetStore();
    await loadProject();

    await waitFor(() => {
      const docs: ProjectDocument[] = useProjectStore.getState().currProject?.documents ?? [];
      expect(docs.some(d => d.desc.datasourceName === 'AgentDoc1')).toBe(true);
      expect(docs.some(d => d.desc.datasourceName === 'IntegDoc1')).toBe(false);
    }, { timeout: 15000 });
  }, 30000);
}, 180000);
