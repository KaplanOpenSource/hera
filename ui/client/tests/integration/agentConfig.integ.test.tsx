/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { startDockerEnv, type DockerEnv } from './dockerSetup';
import { createProjectViaUI, addDocumentViaUI, renderApp, resetStore } from './integHelpers';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock(8005));
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());

import { useProjectStore } from '../../src/stores/useProjectStore';

let env: DockerEnv;

describe('Agent config UI integration', () => {
  beforeAll(async () => {
    env = await startDockerEnv({
      network: 'hera-test-agent-net',
      mongoContainer: 'hera-test-agent-mongo',
      serverContainer: 'hera-test-agent-server',
      serverPort: 8005,
      dbName: 'hera_test_agent',
    });
    await createProjectViaUI('AgentTestProject');
    await addDocumentViaUI('AgentTestProject', 'TestAgent', { agent: true });
  }, 90000);

  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });

  afterAll(() => {
    cleanup();
    env?.cleanup();
  }, 15000);

  // Render the full app, wait for project to load, click the agent doc in the tree
  const openAgentDoc = async () => {
    renderApp('/AgentTestProject');

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('AgentTestProject');
    }, { timeout: 15000 });

    // Click the agent document name in the tree to select it
    const docLabel = await screen.findByText('TestAgent', {}, { timeout: 10000 });
    await act(async () => { fireEvent.click(docLabel); });

    // Wait for the detail panel to load the document
    await waitFor(() => {
      // The document name appears in the detail panel header (h6)
      const headings = screen.getAllByText('TestAgent');
      expect(headings.length).toBeGreaterThanOrEqual(2); // one in tree, one in detail panel
    }, { timeout: 10000 });
  };

  it('shows agent config editor with effects and physical properties', async () => {
    await openAgentDoc();

    await waitFor(() => {
      expect(screen.getByText('Effects')).toBeTruthy();
      expect(screen.getByText('Physical Properties')).toBeTruthy();
    }, { timeout: 5000 });
  }, 45000);

  it('can add a new effect via the UI', async () => {
    await openAgentDoc();

    await waitFor(() => {
      expect(screen.getByText('Effects')).toBeTruthy();
    }, { timeout: 5000 });

    const newEffectInput = screen.getByRole('textbox', { name: /new effect name/i });
    fireEvent.change(newEffectInput, { target: { value: 'Inhalation' } });

    const addBtn = screen.getByRole('button', { name: /add effect/i });
    await act(async () => { fireEvent.click(addBtn); });

    expect(screen.getByText('Inhalation')).toBeTruthy();

    await waitFor(() => {
      expect(screen.getByLabelText('Update Document')).toBeTruthy();
    });

    const updateWrapper = screen.getByLabelText('Update Document');
    await act(async () => {
      fireEvent.click(within(updateWrapper).getByRole('button'));
    });

    await waitFor(() => {
      expect(screen.queryByLabelText('Update Document')).toBeNull();
    }, { timeout: 10000 });
  }, 45000);

  it('effect persists after simulated refresh', async () => {
    await openAgentDoc();

    await waitFor(() => {
      expect(screen.getByText('Inhalation')).toBeTruthy();
    }, { timeout: 10000 });
  }, 45000);

  it('can change the Ten Berge Coefficient and save', async () => {
    await openAgentDoc();

    const tbField = await screen.findByRole('spinbutton', { name: /ten berge coefficient/i });
    fireEvent.change(tbField, { target: { value: '3.5' } });

    const updateWrapper = await screen.findByLabelText('Update Document');
    await act(async () => {
      fireEvent.click(within(updateWrapper).getByRole('button'));
    });

    await waitFor(() => {
      expect(screen.queryByLabelText('Update Document')).toBeNull();
    }, { timeout: 10000 });

    expect((screen.getByRole('spinbutton', { name: /ten berge coefficient/i }) as HTMLInputElement).value).toBe('3.5');
  }, 45000);

  it('Ten Berge Coefficient persists after simulated refresh', async () => {
    await openAgentDoc();

    await waitFor(() => {
      const field = screen.getByRole('spinbutton', { name: /ten berge coefficient/i }) as HTMLInputElement;
      expect(field.value).toBe('3.5');
    }, { timeout: 10000 });
  }, 45000);

  it('can edit molecular weight and save', async () => {
    await openAgentDoc();

    await waitFor(() => {
      expect(screen.getByText('Physical Properties')).toBeTruthy();
    }, { timeout: 5000 });

    const mwField = screen.getByRole('textbox', { name: /molecular weight/i });
    fireEvent.change(mwField, { target: { value: '100*g/mol' } });

    const updateWrapper = await screen.findByLabelText('Update Document');
    await act(async () => {
      fireEvent.click(within(updateWrapper).getByRole('button'));
    });

    await waitFor(() => {
      expect(screen.queryByLabelText('Update Document')).toBeNull();
    }, { timeout: 10000 });
  }, 45000);

  it('molecular weight persists after simulated refresh', async () => {
    await openAgentDoc();

    await waitFor(() => {
      const field = screen.getByRole('textbox', { name: /molecular weight/i }) as HTMLInputElement;
      expect(field.value).toBe('100*g/mol');
    }, { timeout: 10000 });
  }, 45000);
}, 360000);
