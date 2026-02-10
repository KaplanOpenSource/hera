import { Add } from "@mui/icons-material";
import { Button, Stack, TextField, Typography } from "@mui/material";
import { useState } from "react";
import type { AgentEffect } from "../../shared/AgentConfig";
import { makeDefaultCalculator } from "./CalculatorEditor";
import { EffectEditor } from "./EffectEditor";
import { EFFECT_TO_LEVEL_TYPE, EffectType } from "./LevelParamsEditor";

function makeDefaultEffect(type: EffectType): AgentEffect {
  return {
    type,
    calculator: makeDefaultCalculator("Haber"),
    parameters: {
      type: EFFECT_TO_LEVEL_TYPE[type],
      levels: [],
      parameters: {},
    },
  } as AgentEffect;
}

export const EffectsListEditor = ({
  effects,
  onChange,
}: {
  effects: Record<string, AgentEffect>;
  onChange: (effects: Record<string, AgentEffect>) => void;
}) => {
  const [newEffectName, setNewEffectName] = useState("");

  const addEffect = () => {
    const name = newEffectName.trim();
    if (!name || effects[name]) {
      return;
    }
    onChange({ ...effects, [name]: makeDefaultEffect("Threshold") });
    setNewEffectName("");
  };

  const handleUpdate = (oldName: string) => (newName?: string, newEffect?: AgentEffect) => {
    if (!newName || !newEffect) {
      // delete
      const { [oldName]: _, ...rest } = effects;
      onChange(rest);
    } else if (newName !== oldName) {
      // rename (also applies the new effect value)
      if (effects[newName]) return;
      const entries = Object.entries(effects).map(([k, v]) =>
        k === oldName ? [newName, newEffect] : [k, v]
      );
      onChange(Object.fromEntries(entries));
    } else {
      // update
      onChange({ ...effects, [oldName]: newEffect });
    }
  };

  return (
    <>
      <Stack direction="row" alignItems="center" justifyContent="space-between" sx={{ mb: 1 }}>
        <Typography variant="h6">Effects</Typography>
      </Stack>

      <Stack spacing={1} sx={{ mb: 2 }}>
        {Object.entries(effects).map(([name, effect]) => (
          <EffectEditor
            key={name}
            name={name}
            effect={effect}
            onUpdate={handleUpdate(name)}
          />
        ))}
      </Stack>

      <Stack direction="row" spacing={1}>
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
          startIcon={<Add />}
          onClick={addEffect}
          disabled={!newEffectName.trim() || !!effects[newEffectName.trim()]}
        >
          Add Effect
        </Button>
      </Stack>
    </>
  );
};
