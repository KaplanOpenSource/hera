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
  Select,
  Stack,
  TextField,
  Typography
} from "@mui/material";
import { useCallback, useState } from "react";
import { SectionHeader } from "../../elements/SectionHeader";
import type { AgentConfig, AgentEffect, Calculator } from "../../shared/AgentConfig";
import { CalculatorEditor, getCalculatorType, makeDefaultCalculator } from "./CalculatorEditor";
import { InjuryLevelsEditor } from "./InjuryLevelsEditor";
import {
  EFFECT_TO_LEVEL_TYPE,
  EFFECT_TYPES,
  EffectType
} from "./LevelParamsEditor";
import { PhysicalPropertiesEditor } from "./PhysicalPropertiesEditor";

// =============================================================================
// Helpers
// =============================================================================

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
          <CalculatorEditor
            calculator={effect.calculator}
            setCalculator={(calc: Calculator) => {
              onChange({ ...effect, calculator: calc } as AgentEffect);
            }}
          />

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
        slotProps={{ htmlInput: { step: 0.1 } }}
        helperText="Global exponent n used by TenBerge calculators"
        value={agentResource.effectParameters?.tenbergeCoefficient ?? ""}
        onChange={(e) => {
          const tenbergeCoefficient = e.target.value ? parseFloat(e.target.value) : undefined;
          setAgentResource({
            ...agentResource,
            effectParameters: { ...agentResource.effectParameters, tenbergeCoefficient },
          })
        }}
      />

      <Divider sx={{ mb: 2 }} />

      {/* Physical properties */}
      <Stack direction="row" alignItems="center" justifyContent="space-between" sx={{ mb: 2 }}>
        <Typography variant="h6">Physical Properties</Typography>
      </Stack>

      <PhysicalPropertiesEditor
        properties={agentResource.physicalProperties || {}}
        onChange={(p) => setAgentResource({ ...agentResource, physicalProperties: p })}
      />
    </Box>
  );
};