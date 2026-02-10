import {
  Delete,
  ExpandMore,
  Science,
} from "@mui/icons-material";
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Chip,
  Divider,
  FormControl,
  IconButton,
  InputLabel,
  MenuItem,
  Select,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { useState } from "react";
import { SectionHeader } from "../../elements/SectionHeader";
import { AgentEffect, Calculator } from "../../shared/AgentConfig";
import { CalculatorEditor, getCalculatorType } from "./CalculatorEditor";
import { InjuryLevelsEditor } from "./InjuryLevelsEditor";
import { EFFECT_TO_LEVEL_TYPE, EFFECT_TYPES, EffectType } from "./LevelParamsEditor";

export const EffectEditor = ({
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
      <AccordionSummary expandIcon={<ExpandMore />}>
        <Stack direction="row" alignItems="center" spacing={1} sx={{ flex: 1, mr: 1 }}>
          <Science fontSize="small" color="action" />
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
          <Delete fontSize="small" />
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