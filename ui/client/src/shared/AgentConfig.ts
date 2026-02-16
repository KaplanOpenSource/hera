// =============================================================================
// Agent Configuration Types — hera/riskassessment/agents/
// =============================================================================

// -----------------------------------------------------------------------------
// Calculators
// -----------------------------------------------------------------------------

export interface CalculatorHaber {
  Haber: {
    /** Breathing rate in L/min. Default: 10 */
    breathingRate?: number;
  };
}

export interface CalculatorTenBerge {
  TenBerge: {
    /** Breathing rate in L/min. Default: 10 */
    breathingRate?: number;
  };
}

export interface CalculatorMaxConcentration {
  MaxConcentration: {
    /** Rolling-window size, e.g. "10min", "30min", "1h", "8h" */
    sampling: string;
    /** Breathing rate in L/min. Default: 10 */
    breathingRate?: number;
  };
}

export type Calculator =
  | CalculatorHaber
  | CalculatorTenBerge
  | CalculatorMaxConcentration;

// -----------------------------------------------------------------------------
// Injury Level Types
// -----------------------------------------------------------------------------

export interface InjuryLevelThresholdParams {
  /** Concentration threshold as a unit expression, e.g. "106*mg/m**3" */
  threshold: string;
}

export interface InjuryLevelLognormal10DoseResponseParams {
  /** Toxic load at which 50 % of the population is affected (mg/m³) */
  TL_50: number;
  /** Standard deviation of the log-normal distribution (base 10) */
  sigma: number;
}

export interface InjuryLevelExponentialParams {
  /** Rate constant for the exponential model (may be stringified) */
  k: number | string;
}

// -----------------------------------------------------------------------------
// Injury Levels
// -----------------------------------------------------------------------------

interface InjuryLevelsBase {
  /** Severity levels ordered highest → lowest */
  levels: string[];
}

export interface InjuryLevelsThreshold extends InjuryLevelsBase {
  type: "Threshold";
  parameters: { [key: string]: InjuryLevelThresholdParams };
}

export interface InjuryLevelsLognormal10 extends InjuryLevelsBase {
  type: "Lognormal10DoseResponse";
  parameters: { [key: string]: InjuryLevelLognormal10DoseResponseParams };
}

export interface InjuryLevelsExponential extends InjuryLevelsBase {
  type: "Exponential";
  parameters: { [key: string]: InjuryLevelExponentialParams };
}

// -----------------------------------------------------------------------------
// Agent Effects
// -----------------------------------------------------------------------------

interface AgentEffectBase {
  calculator: Calculator;
  units?: string;
}

export interface AgentEffectThreshold extends AgentEffectBase {
  type: "Threshold";
  parameters: InjuryLevelsThreshold;
}

export interface AgentEffectLognormal10 extends AgentEffectBase {
  type: "Lognormal10";
  parameters: InjuryLevelsLognormal10;
}

export interface AgentEffectExponential extends AgentEffectBase {
  type: "Exponential";
  parameters: InjuryLevelsExponential;
}

export type AgentEffect =
  | AgentEffectThreshold
  | AgentEffectLognormal10
  | AgentEffectExponential;

// -----------------------------------------------------------------------------
// Vapor Pressure
// -----------------------------------------------------------------------------

export interface VaporPressure {
  A: number;
  B: number;
  C: number;
  D?: number;
  E?: number;
  F?: number;
  /** Unit expression, e.g. "mmHg" */
  units?: string;
}

// -----------------------------------------------------------------------------
// Top-Level Agent Configuration
// -----------------------------------------------------------------------------

export interface AgentConfig {
  /** Global parameters shared across all effects */
  effectParameters?: {
    /** Ten Berge exponent n (used by CalculatorTenBerge) */
    tenbergeCoefficient?: number;
    [key: string]: unknown;
  };

  /** Dictionary of named effects */
  effects: { [key: string]: AgentEffect };

  /** Physical / chemical constants */
  physicalProperties?: {
    /** Molecular weight as unit expression, e.g. "36.46*g/mol" */
    molecularWeight?: string;
    /** Sorption coefficient as unit expression, e.g. "0.13333*cm/s" */
    sorptionCoefficient?: string;
    spreadFactor?: number;
    /** Antoine-like constants [a, b, c, d] */
    volatilityConstants?: [number, number, number, number];
    /** Linear density model [a, b, c] → density = a - b*(T - c) g/cm³ */
    densityConstants?: [number, number, number];
    molecularVolume?: number;
    /** Extended Antoine equation parameters */
    vaporPressure?: VaporPressure;
  };
}