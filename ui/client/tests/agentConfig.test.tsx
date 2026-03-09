import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { AgentConfigEditor } from '../src/components/agents/AgentConfigEditor';
import { EffectsListEditor } from '../src/components/agents/EffectsListEditor';
import { makeDefaultCalculator, getCalculatorType } from '../src/components/agents/CalculatorEditor';
import { LevelParamsEditor, makeDefaultLevelParams, EFFECT_TO_LEVEL_TYPE } from '../src/components/agents/LevelParamsEditor';
import type { AgentConfig, AgentEffect } from '../src/shared/AgentConfig';

afterEach(() => {
  cleanup();
});

// --- EffectsListEditor ---

describe('EffectsListEditor', () => {
  const makeEffects = (): Record<string, AgentEffect> => ({});

  it('renders empty state with add controls', () => {
    const onChange = vi.fn();
    render(<EffectsListEditor effects={makeEffects()} onChange={onChange} />);
    expect(screen.getByText('Effects')).toBeTruthy();
    expect(screen.getByRole('textbox', { name: /new effect name/i })).toBeTruthy();
    expect((screen.getByRole('button', { name: /add effect/i }) as HTMLButtonElement).disabled).toBe(true);
  });

  it('add button is enabled when name is typed', () => {
    render(<EffectsListEditor effects={{}} onChange={vi.fn()} />);
    const input = screen.getByRole('textbox', { name: /new effect name/i });
    fireEvent.change(input, { target: { value: 'Inhalation' } });
    expect((screen.getByRole('button', { name: /add effect/i }) as HTMLButtonElement).disabled).toBe(false);
  });

  it('calls onChange with new effect on add', () => {
    const onChange = vi.fn();
    render(<EffectsListEditor effects={{}} onChange={onChange} />);
    fireEvent.change(screen.getByRole('textbox', { name: /new effect name/i }), {
      target: { value: 'Inhalation' },
    });
    fireEvent.click(screen.getByRole('button', { name: /add effect/i }));

    expect(onChange).toHaveBeenCalledTimes(1);
    const newEffects = onChange.mock.calls[0][0];
    expect(newEffects).toHaveProperty('Inhalation');
    expect(newEffects.Inhalation.type).toBe('Threshold');
  });

  it('prevents adding duplicate effect names', () => {
    const existing: Record<string, AgentEffect> = {
      Inhalation: {
        type: 'Threshold',
        calculator: { Haber: { breathingRate: 10 } },
        parameters: { type: 'Threshold', levels: [], parameters: {} },
      },
    };
    render(<EffectsListEditor effects={existing} onChange={vi.fn()} />);
    fireEvent.change(screen.getByRole('textbox', { name: /new effect name/i }), {
      target: { value: 'Inhalation' },
    });
    expect((screen.getByRole('button', { name: /add effect/i }) as HTMLButtonElement).disabled).toBe(true);
  });

  it('renders existing effects', () => {
    const effects: Record<string, AgentEffect> = {
      Skin: {
        type: 'Threshold',
        calculator: { Haber: { breathingRate: 10 } },
        parameters: { type: 'Threshold', levels: [], parameters: {} },
      },
    };
    render(<EffectsListEditor effects={effects} onChange={vi.fn()} />);
    expect(screen.getByText('Skin')).toBeTruthy();
  });
});

// --- makeDefaultCalculator / getCalculatorType ---

describe('Calculator utilities', () => {
  it('makeDefaultCalculator creates Haber with breathingRate 10', () => {
    const calc = makeDefaultCalculator('Haber');
    expect(calc).toEqual({ Haber: { breathingRate: 10 } });
  });

  it('makeDefaultCalculator creates TenBerge', () => {
    const calc = makeDefaultCalculator('TenBerge');
    expect(calc).toEqual({ TenBerge: { breathingRate: 10 } });
  });

  it('makeDefaultCalculator creates MaxConcentration with sampling', () => {
    const calc = makeDefaultCalculator('MaxConcentration');
    expect(calc).toEqual({ MaxConcentration: { sampling: '10min', breathingRate: 10 } });
  });

  it('getCalculatorType identifies Haber', () => {
    expect(getCalculatorType({ Haber: { breathingRate: 10 } })).toBe('Haber');
  });

  it('getCalculatorType identifies TenBerge', () => {
    expect(getCalculatorType({ TenBerge: { breathingRate: 5 } })).toBe('TenBerge');
  });

  it('getCalculatorType identifies MaxConcentration', () => {
    expect(getCalculatorType({ MaxConcentration: { sampling: '1h', breathingRate: 10 } })).toBe('MaxConcentration');
  });
});

// --- LevelParamsEditor defaults ---

describe('makeDefaultLevelParams', () => {
  it('Threshold has threshold string', () => {
    expect(makeDefaultLevelParams('Threshold')).toEqual({ threshold: '1*mg/m**3' });
  });

  it('Lognormal10DoseResponse has TL_50 and sigma', () => {
    expect(makeDefaultLevelParams('Lognormal10DoseResponse')).toEqual({ TL_50: 1, sigma: 0.5 });
  });

  it('Exponential has k', () => {
    expect(makeDefaultLevelParams('Exponential')).toEqual({ k: 1 });
  });
});

// --- EFFECT_TO_LEVEL_TYPE mapping ---

describe('EFFECT_TO_LEVEL_TYPE', () => {
  it('maps Threshold to Threshold', () => {
    expect(EFFECT_TO_LEVEL_TYPE.Threshold).toBe('Threshold');
  });

  it('maps Lognormal10 to Lognormal10DoseResponse', () => {
    expect(EFFECT_TO_LEVEL_TYPE.Lognormal10).toBe('Lognormal10DoseResponse');
  });

  it('maps Exponential to Exponential', () => {
    expect(EFFECT_TO_LEVEL_TYPE.Exponential).toBe('Exponential');
  });
});

// --- LevelParamsEditor rendering ---

describe('LevelParamsEditor', () => {
  it('renders threshold field for Threshold type', () => {
    render(
      <LevelParamsEditor
        levelType="Threshold"
        params={{ threshold: '50*mg/m**3' }}
        onChange={vi.fn()}
      />
    );
    expect(screen.getByRole('textbox', { name: /threshold/i })).toBeTruthy();
  });

  it('renders TL_50 and Sigma for Lognormal10DoseResponse', () => {
    render(
      <LevelParamsEditor
        levelType="Lognormal10DoseResponse"
        params={{ TL_50: 100, sigma: 0.3 }}
        onChange={vi.fn()}
      />
    );
    expect(screen.getByRole('spinbutton', { name: /tl_50/i })).toBeTruthy();
    expect(screen.getByRole('spinbutton', { name: /sigma/i })).toBeTruthy();
  });

  it('renders k field for Exponential', () => {
    render(
      <LevelParamsEditor
        levelType="Exponential"
        params={{ k: 0.5 }}
        onChange={vi.fn()}
      />
    );
    expect(screen.getByRole('textbox', { name: /k \(rate constant\)/i })).toBeTruthy();
  });

  it('calls onChange when threshold is edited', () => {
    const onChange = vi.fn();
    render(
      <LevelParamsEditor
        levelType="Threshold"
        params={{ threshold: '1*mg/m**3' }}
        onChange={onChange}
      />
    );
    fireEvent.change(screen.getByRole('textbox', { name: /threshold/i }), {
      target: { value: '200*mg/m**3' },
    });
    expect(onChange).toHaveBeenCalledWith({ threshold: '200*mg/m**3' });
  });
});

// --- AgentConfigEditor ---

describe('AgentConfigEditor', () => {
  const baseConfig: AgentConfig = {
    effects: {},
    effectParameters: { tenbergeCoefficient: 2.0 },
    physicalProperties: {},
  };

  it('renders effects section and physical properties', () => {
    render(<AgentConfigEditor agentResource={baseConfig} setAgentResource={vi.fn()} />);
    expect(screen.getByText('Effects')).toBeTruthy();
    expect(screen.getByText('Physical Properties')).toBeTruthy();
  });

  it('renders Ten Berge Coefficient field', () => {
    render(<AgentConfigEditor agentResource={baseConfig} setAgentResource={vi.fn()} />);
    const field = screen.getByRole('spinbutton', { name: /ten berge coefficient/i });
    expect(field).toBeTruthy();
    expect((field as HTMLInputElement).value).toBe('2');
  });

  it('calls setAgentResource when Ten Berge Coefficient changes', () => {
    const setAgent = vi.fn();
    render(<AgentConfigEditor agentResource={baseConfig} setAgentResource={setAgent} />);
    fireEvent.change(screen.getByRole('spinbutton', { name: /ten berge coefficient/i }), {
      target: { value: '3.5' },
    });
    expect(setAgent).toHaveBeenCalledTimes(1);
    expect(setAgent.mock.calls[0][0].effectParameters.tenbergeCoefficient).toBe(3.5);
  });
});
