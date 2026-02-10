import AddIcon from "@mui/icons-material/Add";
import {
  Box,
  Button,
  Divider,
  Stack,
  TextField,
  Typography
} from "@mui/material";
import { useCallback, useState } from "react";
import type { AgentConfig, AgentEffect } from "../../shared/AgentConfig";
import { makeDefaultCalculator } from "./CalculatorEditor";
import { EffectEditor } from "./EffectEditor";
import {
  EFFECT_TO_LEVEL_TYPE,
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