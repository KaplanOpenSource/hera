/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { createProjectViaUI, addDocumentViaUI, renderApp, resetStore } from './integHelpers';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock());
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());
vi.mock('../../src/io/snackbar', async () => (await import('./mockFactories')).createSnackbarMock());

import { useProjectStore } from '../../src/stores/useProjectStore';

describe('Agent config UI integration', () => {
  beforeAll(async () => {
    await createProjectViaUI('AgentTestProject');
    await addDocumentViaUI('AgentTestProject', 'TestAgent', { agent: true });
  }, 90000);

  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });
  afterAll(() => { cleanup(); });

  const openAgentDoc = async () => {
    renderApp('/AgentTestProject');

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('AgentTestProject');
    }, { timeout: 15000 });

    const docLabel = await screen.findByText('TestAgent', {}, { timeout: 10000 });
    await act(async () => { fireEvent.click(docLabel); });

    await waitFor(() => {
      const headings = screen.getAllByText('TestAgent');
      expect(headings.length).toBeGreaterThanOrEqual(2);
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
