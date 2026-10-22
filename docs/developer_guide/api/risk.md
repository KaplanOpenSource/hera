# Risk Assessment API

Toolkit for hazard modeling, agent effects, and protection policy evaluation.

## RiskToolkit

::: hera.riskassessment.riskToolkit.RiskToolkit
    options:
      show_source: true
      members_order: source

## Agent

::: hera.riskassessment.agents.Agents.Agent
    options:
      show_source: true
      members_order: source

## ProtectionPolicy

::: hera.riskassessment.protectionpolicy.ProtectionPolicy.ProtectionPolicy
    options:
      show_source: true
      members_order: source

## Effects

Injury models and calculators. Subclass names follow the `Calculator<Name>`, `InjuryLevel<Name>`, `Injury<Name>` convention because the factories resolve them with `pydoc.locate`.

### Calculator

::: hera.riskassessment.agents.effects.Calculator
    options:
      show_source: true
      members_order: source

### InjuryLevel

::: hera.riskassessment.agents.effects.InjuryLevel
    options:
      show_source: true
      members_order: source

### Injury

::: hera.riskassessment.agents.effects.Injury
    options:
      show_source: true
      members_order: source
