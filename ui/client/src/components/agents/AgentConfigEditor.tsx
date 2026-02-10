import AddIcon from "@mui/icons-material/Add";
import DeleteIcon from "@mui/icons-material/Delete";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ScienceIcon from "@mui/icons-material/Science";
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Box,
  Button,
  Chip,
  Divider,
  FormControl,
  IconButton,
  InputLabel,
  MenuItem,
  Paper,
  Select,
  Stack,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import React, { useCallback, useState } from "react";
import type { AgentConfig, AgentEffect, Calculator, VaporPressure } from "../../shared/AgentConfig";
import { SectionHeader } from "../../elements/SectionHeader";

// =============================================================================
// Helpers
// =============================================================================

type EffectType = AgentEffect["type"];
type CalculatorType = "Haber" | "TenBerge" | "MaxConcentration";
type InjuryLevelType = "Threshold" | "Lognormal10DoseResponse" | "Exponential";

const EFFECT_TYPES: EffectType[] = ["Threshold", "Lognormal10", "Exponential"];
const CALCULATOR_TYPES: CalculatorType[] = ["Haber", "TenBerge", "MaxConcentration"];

const EFFECT_TO_LEVEL_TYPE: { [key in EffectType]: InjuryLevelType } = {
  Threshold: "Threshold",
  Lognormal10: "Lognormal10DoseResponse",
  Exponential: "Exponential",
};

function getCalculatorType(calc: Calculator): CalculatorType {
  if ("Haber" in calc) return "Haber";
  if ("TenBerge" in calc) return "TenBerge";
  return "MaxConcentration";
}

function getCalculatorParams(calc: Calculator) {
  if ("Haber" in calc) return calc.Haber;
  if ("TenBerge" in calc) return calc.TenBerge;
  return calc.MaxConcentration;
}

function makeDefaultCalculator(type: CalculatorType): Calculator {
  switch (type) {
    case "Haber":
      return { Haber: { breathingRate: 10 } };
    case "TenBerge":
      return { TenBerge: { breathingRate: 10 } };
    case "MaxConcentration":
      return { MaxConcentration: { sampling: "10min", breathingRate: 10 } };
  }
}

function makeDefaultLevelParams(levelType: InjuryLevelType) {
  switch (levelType) {
    case "Threshold":
      return { threshold: "1*mg/m**3" };
    case "Lognormal10DoseResponse":
      return { TL_50: 1, sigma: 0.5 };
    case "Exponential":
      return { k: 1 };
  }
}

function makeDefaultEffect(type: EffectType): AgentEffect {
  const levelType = EFFECT_TO_LEVEL_TYPE[type];
  return {
    type,
    calculator: makeDefaultCalculator("Haber"),
    parameters: {
      type: levelType,
      levels: [],
      parameters: {},
    },
  } as AgentEffect;
}

// =============================================================================
// Sub-editors
// =============================================================================

// -- Calculator editor --------------------------------------------------------

const CalculatorEditor = ({
  calculator,
  onChange,
}: {
  calculator: Calculator;
  onChange: (calc: Calculator) => void;
}) => {
  const type = getCalculatorType(calculator);
  const params = getCalculatorParams(calculator);

  const handleTypeChange = (newType: CalculatorType) => {
    onChange(makeDefaultCalculator(newType));
  };

  const updateParam = (key: string, value: any) => {
    const newParams = { ...params, [key]: value };
    switch (type) {
      case "Haber":
        onChange({ Haber: newParams as any });
        break;
      case "TenBerge":
        onChange({ TenBerge: newParams as any });
        break;
      case "MaxConcentration":
        onChange({ MaxConcentration: newParams as any });
        break;
    }
  };

  return (
    <Stack spacing={2}>
      <FormControl size="small" fullWidth>
        <InputLabel>Calculator</InputLabel>
        <Select
          value={type}
          label="Calculator"
          onChange={(e) => handleTypeChange(e.target.value as CalculatorType)}
        >
          {CALCULATOR_TYPES.map((t) => (
            <MenuItem key={t} value={t}>
              {t}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <TextField
        label="Breathing Rate (L/min)"
        type="number"
        size="small"
        value={params.breathingRate ?? 10}
        onChange={(e) => updateParam("breathingRate", parseFloat(e.target.value) || 0)}
      />

      {type === "MaxConcentration" && "sampling" in params && (
        <TextField
          label="Sampling Window"
          size="small"
          helperText='e.g. "10min", "30min", "1h", "8h"'
          value={(params as any).sampling ?? ""}
          onChange={(e) => updateParam("sampling", e.target.value)}
        />
      )}
    </Stack>
  );
};

// -- Injury level params editor -----------------------------------------------

const LevelParamsEditor = ({
  levelType,
  params,
  onChange,
}: {
  levelType: InjuryLevelType;
  params: any;
  onChange: (params: any) => void;
}) => {
  const update = (key: string, value: any) => onChange({ ...params, [key]: value });

  switch (levelType) {
    case "Threshold":
      return (
        <TextField
          label="Threshold"
          size="small"
          fullWidth
          helperText='Unit expression, e.g. "106*mg/m**3"'
          value={params.threshold ?? ""}
          onChange={(e) => update("threshold", e.target.value)}
        />
      );
    case "Lognormal10DoseResponse":
      return (
        <Stack direction="row" spacing={2}>
          <TextField
            label="TL_50"
            type="number"
            size="small"
            fullWidth
            value={params.TL_50 ?? ""}
            onChange={(e) => update("TL_50", parseFloat(e.target.value) || 0)}
          />
          <TextField
            label="Sigma"
            type="number"
            size="small"
            fullWidth
            inputProps={{ step: 0.1 }}
            value={params.sigma ?? ""}
            onChange={(e) => update("sigma", parseFloat(e.target.value) || 0)}
          />
        </Stack>
      );
    case "Exponential":
      return (
        <TextField
          label="k (rate constant)"
          size="small"
          fullWidth
          value={params.k ?? ""}
          onChange={(e) => {
            const num = parseFloat(e.target.value);
            update("k", isNaN(num) ? e.target.value : num);
          }}
        />
      );
  }
};

// -- Injury levels editor -----------------------------------------------------

const InjuryLevelsEditor = ({
  levelType,
  levels,
  parameters,
  onChange,
}: {
  levelType: InjuryLevelType;
  levels: string[];
  parameters: { [key: string]: any };
  onChange: (levels: string[], parameters: { [key: string]: any }) => void;
}) => {
  const [newLevelName, setNewLevelName] = useState("");

  const addLevel = () => {
    const name = newLevelName.trim();
    if (!name || levels.includes(name)) return;
    onChange([...levels, name], { ...parameters, [name]: makeDefaultLevelParams(levelType) });
    setNewLevelName("");
  };

  const removeLevel = (name: string) => {
    const { [name]: _, ...rest } = parameters;
    onChange(
      levels.filter((l) => l !== name),
      rest
    );
  };

  const reorderLevel = (index: number, direction: -1 | 1) => {
    const target = index + direction;
    if (target < 0 || target >= levels.length) return;
    const next = [...levels];
    [next[index], next[target]] = [next[target], next[index]];
    onChange(next, parameters);
  };

  const updateLevelParams = (name: string, params: any) => {
    onChange(levels, { ...parameters, [name]: params });
  };

  return (
    <Stack spacing={1}>
      <Typography variant="caption" color="text.secondary">
        Levels are ordered highest severity first
      </Typography>

      {levels.map((name, i) => (
        <Paper key={name} variant="outlined" sx={{ p: 1.5 }}>
          <Stack spacing={1.5}>
            <Stack direction="row" alignItems="center" justifyContent="space-between">
              <Stack direction="row" alignItems="center" spacing={1}>
                <Chip label={name} size="small" color="primary" variant="outlined" />
                <Typography variant="caption" color="text.secondary">
                  #{i + 1}
                </Typography>
              </Stack>
              <Stack direction="row" spacing={0.5}>
                <Tooltip title="Move up (higher severity)">
                  <span>
                    <IconButton size="small" disabled={i === 0} onClick={() => reorderLevel(i, -1)}>
                      ▲
                    </IconButton>
                  </span>
                </Tooltip>
                <Tooltip title="Move down (lower severity)">
                  <span>
                    <IconButton
                      size="small"
                      disabled={i === levels.length - 1}
                      onClick={() => reorderLevel(i, 1)}
                    >
                      ▼
                    </IconButton>
                  </span>
                </Tooltip>
                <Tooltip title="Remove level">
                  <IconButton size="small" color="error" onClick={() => removeLevel(name)}>
                    <DeleteIcon fontSize="small" />
                  </IconButton>
                </Tooltip>
              </Stack>
            </Stack>
            <LevelParamsEditor
              levelType={levelType}
              params={parameters[name] ?? makeDefaultLevelParams(levelType)}
              onChange={(p) => updateLevelParams(name, p)}
            />
          </Stack>
        </Paper>
      ))}

      <Stack direction="row" spacing={1} alignItems="center">
        <TextField
          label="New level name"
          size="small"
          value={newLevelName}
          onChange={(e) => setNewLevelName(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && addLevel()}
          sx={{ flex: 1 }}
        />
        <Button
          variant="outlined"
          size="small"
          startIcon={<AddIcon />}
          onClick={addLevel}
          disabled={!newLevelName.trim() || levels.includes(newLevelName.trim())}
        >
          Add
        </Button>
      </Stack>
    </Stack>
  );
};

// -- Effect editor ------------------------------------------------------------

const EffectEditor = ({
  name,
  effect,
  onChange,
  onRename,
  onDelete,
}: {
  name: string;
  effect: AgentEffect;
  onChange: (effect: AgentEffect) => void;
  onRename: (oldName: string, newName: string) => void;
  onDelete: () => void;
}) => {
  const [editingName, setEditingName] = useState(false);
  const [nameValue, setNameValue] = useState(name);

  const commitRename = () => {
    const trimmed = nameValue.trim();
    if (trimmed && trimmed !== name) {
      onRename(name, trimmed);
    } else {
      setNameValue(name);
    }
    setEditingName(false);
  };

  const handleTypeChange = (newType: EffectType) => {
    const newLevelType = EFFECT_TO_LEVEL_TYPE[newType];
    onChange({
      ...effect,
      type: newType,
      parameters: {
        type: newLevelType,
        levels: [],
        parameters: {},
      },
    } as AgentEffect);
  };

  const updateCalculator = (calc: Calculator) => {
    onChange({ ...effect, calculator: calc } as AgentEffect);
  };

  const updateLevels = (levels: string[], parameters: { [key: string]: any }) => {
    onChange({
      ...effect,
      parameters: { ...effect.parameters, levels, parameters },
    } as AgentEffect);
  };

  return (
    <Accordion defaultExpanded={false} disableGutters>
      <AccordionSummary expandIcon={<ExpandMoreIcon />}>
        <Stack direction="row" alignItems="center" spacing={1} sx={{ flex: 1, mr: 1 }}>
          <ScienceIcon fontSize="small" color="action" />
          {editingName ? (
            <TextField
              size="small"
              variant="standard"
              value={nameValue}
              autoFocus
              onChange={(e) => setNameValue(e.target.value)}
              onBlur={commitRename}
              onKeyDown={(e) => {
                if (e.key === "Enter") commitRename();
                if (e.key === "Escape") {
                  setNameValue(name);
                  setEditingName(false);
                }
              }}
              onClick={(e) => e.stopPropagation()}
            />
          ) : (
            <Typography
              variant="subtitle2"
              sx={{ cursor: "pointer", "&:hover": { textDecoration: "underline" } }}
              onClick={(e) => {
                e.stopPropagation();
                setEditingName(true);
              }}
            >
              {name}
            </Typography>
          )}
          <Chip label={effect.type} size="small" />
          <Chip label={getCalculatorType(effect.calculator)} size="small" variant="outlined" />
        </Stack>
        <IconButton
          size="small"
          color="error"
          onClick={(e) => {
            e.stopPropagation();
            onDelete();
          }}
        >
          <DeleteIcon fontSize="small" />
        </IconButton>
      </AccordionSummary>
      <AccordionDetails>
        <Stack spacing={2}>
          <Stack direction="row" spacing={2}>
            <FormControl size="small" sx={{ minWidth: 180 }}>
              <InputLabel>Effect Type</InputLabel>
              <Select
                value={effect.type}
                label="Effect Type"
                onChange={(e) => handleTypeChange(e.target.value as EffectType)}
              >
                {EFFECT_TYPES.map((t) => (
                  <MenuItem key={t} value={t}>
                    {t}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
            <TextField
              label="Units"
              size="small"
              value={effect.units ?? ""}
              onChange={(e) => onChange({ ...effect, units: e.target.value || undefined } as AgentEffect)}
            />
          </Stack>

          <Divider />
          <SectionHeader>Calculator</SectionHeader>
          <CalculatorEditor calculator={effect.calculator} onChange={updateCalculator} />

          <Divider />
          <SectionHeader>Injury Levels</SectionHeader>
          <InjuryLevelsEditor
            levelType={EFFECT_TO_LEVEL_TYPE[effect.type]}
            levels={effect.parameters.levels}
            parameters={effect.parameters.parameters as any}
            onChange={updateLevels}
          />
        </Stack>
      </AccordionDetails>
    </Accordion>
  );
};

// -- Physical properties editor -----------------------------------------------

const PhysicalPropertiesEditor = ({
  properties,
  onChange,
}: {
  properties: NonNullable<AgentConfig["physicalProperties"]>;
  onChange: (props: AgentConfig["physicalProperties"]) => void;
}) => {
  const update = (key: string, value: any) => {
    onChange({ ...properties, [key]: value === "" ? undefined : value });
  };

  const updateVaporPressure = (key: keyof VaporPressure, value: any) => {
    onChange({
      ...properties,
      vaporPressure: { ...properties.vaporPressure, [key]: value } as VaporPressure,
    });
  };

  const updateArrayConst = (key: "volatilityConstants" | "densityConstants", index: number, value: number) => {
    const current = properties[key] ? [...properties[key]!] : key === "volatilityConstants" ? [0, 0, 0, 0] : [0, 0, 0];
    (current as number[])[index] = value;
    onChange({ ...properties, [key]: current });
  };

  return (
    <Stack spacing={2}>
      <Stack direction="row" spacing={2}>
        <TextField
          label="Molecular Weight"
          size="small"
          fullWidth
          helperText='e.g. "36.46*g/mol"'
          value={properties.molecularWeight ?? ""}
          onChange={(e) => update("molecularWeight", e.target.value)}
        />
        <TextField
          label="Sorption Coefficient"
          size="small"
          fullWidth
          helperText='e.g. "0.13333*cm/s"'
          value={properties.sorptionCoefficient ?? ""}
          onChange={(e) => update("sorptionCoefficient", e.target.value)}
        />
      </Stack>

      <Stack direction="row" spacing={2}>
        <TextField
          label="Spread Factor"
          size="small"
          type="number"
          value={properties.spreadFactor ?? ""}
          onChange={(e) => update("spreadFactor", e.target.value ? parseFloat(e.target.value) : "")}
        />
        <TextField
          label="Molecular Volume"
          size="small"
          type="number"
          value={properties.molecularVolume ?? ""}
          onChange={(e) => update("molecularVolume", e.target.value ? parseFloat(e.target.value) : "")}
        />
      </Stack>

      <SectionHeader>Volatility Constants [a, b, c, d]</SectionHeader>
      <Stack direction="row" spacing={1}>
        {["a", "b", "c", "d"].map((label, i) => (
          <TextField
            key={label}
            label={label}
            size="small"
            type="number"
            value={properties.volatilityConstants?.[i] ?? ""}
            onChange={(e) => updateArrayConst("volatilityConstants", i, parseFloat(e.target.value) || 0)}
          />
        ))}
      </Stack>

      <SectionHeader>Density Constants [a, b, c]</SectionHeader>
      <Stack direction="row" spacing={1}>
        {["a", "b", "c"].map((label, i) => (
          <TextField
            key={label}
            label={label}
            size="small"
            type="number"
            value={properties.densityConstants?.[i] ?? ""}
            onChange={(e) => updateArrayConst("densityConstants", i, parseFloat(e.target.value) || 0)}
          />
        ))}
      </Stack>

      <SectionHeader>Vapor Pressure (Extended Antoine)</SectionHeader>
      <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
        {(["A", "B", "C", "D", "E", "F"] as const).map((key) => (
          <TextField
            key={key}
            label={key}
            size="small"
            type="number"
            sx={{ width: 100 }}
            value={properties.vaporPressure?.[key] ?? ""}
            onChange={(e) => updateVaporPressure(key, e.target.value ? parseFloat(e.target.value) : undefined)}
          />
        ))}
        <TextField
          label="Units"
          size="small"
          sx={{ width: 140 }}
          value={properties.vaporPressure?.units ?? ""}
          onChange={(e) => updateVaporPressure("units", e.target.value || undefined)}
        />
      </Stack>
    </Stack>
  );
};

// =============================================================================
// Main Editor
// =============================================================================

export const AgentConfigEditor = ({
  agentResource,
  setAgentResource,
}: {
  agentResource: AgentConfig;
  setAgentResource: (agent: AgentConfig) => void;
}) => {
  const [newEffectName, setNewEffectName] = useState("");

  const update = useCallback(
    (patch: Partial<AgentConfig>) => {
      setAgentResource({ ...agentResource, ...patch });
    },
    [agentResource, setAgentResource]
  );

  // -- Effect CRUD --

  const addEffect = () => {
    const name = newEffectName.trim();
    if (!name || agentResource.effects[name]) return;
    update({
      effects: { ...agentResource.effects, [name]: makeDefaultEffect("Threshold") },
    });
    setNewEffectName("");
  };

  const updateEffect = (name: string, effect: AgentEffect) => {
    update({ effects: { ...agentResource.effects, [name]: effect } });
  };

  const deleteEffect = (name: string) => {
    const { [name]: _, ...rest } = agentResource.effects;
    update({ effects: rest });
  };

  const renameEffect = (oldName: string, newName: string) => {
    if (agentResource.effects[newName]) return;
    const entries = Object.entries(agentResource.effects).map(([k, v]) =>
      k === oldName ? [newName, v] : [k, v]
    );
    update({ effects: Object.fromEntries(entries) });
  };

  // -- Physical properties toggle --

  const hasPhysProps = !!agentResource.physicalProperties;

  return (
    <Box sx={{ maxWidth: 800 }}>
      {/* Effects */}
      <Stack direction="row" alignItems="center" justifyContent="space-between" sx={{ mb: 1 }}>
        <Typography variant="h6">Effects</Typography>
      </Stack>

      <Stack spacing={1} sx={{ mb: 2 }}>
        {Object.entries(agentResource.effects).map(([name, effect]) => (
          <EffectEditor
            key={name}
            name={name}
            effect={effect}
            onChange={(e) => updateEffect(name, e)}
            onRename={renameEffect}
            onDelete={() => deleteEffect(name)}
          />
        ))}
      </Stack>

      <Stack direction="row" spacing={1} sx={{ mb: 3 }}>
        <TextField
          label="New effect name"
          size="small"
          value={newEffectName}
          onChange={(e) => setNewEffectName(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && addEffect()}
          sx={{ flex: 1 }}
        />
        <Button
          variant="contained"
          startIcon={<AddIcon />}
          onClick={addEffect}
          disabled={!newEffectName.trim() || !!agentResource.effects[newEffectName.trim()]}
        >
          Add Effect
        </Button>
      </Stack>

      <Divider sx={{ mb: 2 }} />

      <TextField
        label="Ten Berge Coefficient"
        type="number"
        size="small"
        inputProps={{ step: 0.1 }}
        helperText="Global exponent n used by TenBerge calculators"
        value={agentResource.effectParameters?.tenbergeCoefficient ?? ""}
        onChange={(e) =>
          update({
            effectParameters: {
              ...agentResource.effectParameters,
              tenbergeCoefficient: e.target.value ? parseFloat(e.target.value) : undefined,
            },
          })
        }
      />

      <Divider sx={{ mb: 2 }} />

      {/* Physical properties */}
      <Stack direction="row" alignItems="center" justifyContent="space-between" sx={{ mb: 2 }}>
        <Typography variant="h6">Physical Properties</Typography>
        {!hasPhysProps ? (
          <Button size="small" startIcon={<AddIcon />} onClick={() => update({ physicalProperties: {} })}>
            Add
          </Button>
        ) : (
          <Button
            size="small"
            color="error"
            startIcon={<DeleteIcon />}
            onClick={() => update({ physicalProperties: undefined })}
          >
            Remove
          </Button>
        )}
      </Stack>

      {hasPhysProps && (
        <PhysicalPropertiesEditor
          properties={agentResource.physicalProperties!}
          onChange={(p) => update({ physicalProperties: p })}
        />
      )}
    </Box>
  );
};