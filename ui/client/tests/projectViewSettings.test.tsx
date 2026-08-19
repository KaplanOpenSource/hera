import { describe, it, expect, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { ProjectViewSettingsButton } from '../src/components/project/ProjectViewSettingsButton';
import { useViewSettingsStore } from '../src/stores/useViewSettingsStore';

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  useViewSettingsStore.getState().setViewSettings({ alwaysSaveBeforeRun: false });
});

describe('ProjectViewSettingsButton', () => {
  it('toggles "Always save workflow before running" through the store', () => {
    render(<ProjectViewSettingsButton />);

    // Open the settings dialog.
    fireEvent.click(screen.getByRole('button'));

    const checkbox = screen.getByRole('checkbox', { name: /auto save workflow before running/i });
    expect(checkbox).toHaveProperty('checked', false);

    fireEvent.click(checkbox);
    expect(useViewSettingsStore.getState().viewSettings.alwaysSaveBeforeRun).toBe(true);
    expect(checkbox).toHaveProperty('checked', true);
  });

  it('reflects the current store value when opened', () => {
    useViewSettingsStore.getState().setViewSettings({ alwaysSaveBeforeRun: true });
    render(<ProjectViewSettingsButton />);

    fireEvent.click(screen.getByRole('button'));

    expect(screen.getByRole('checkbox', { name: /auto save workflow before running/i }))
      .toHaveProperty('checked', true);
  });
});
