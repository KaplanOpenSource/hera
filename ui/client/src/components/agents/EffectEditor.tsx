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
  TextField
} from "@mui/material";
import { SectionHeader } from "../../elements/SectionHeader";
import { AgentEffect, Calculator } from "../../shared/AgentConfig";
import { RenameField } from "../../elements/RenameField";
import { CalculatorEditor, getCalculatorType } from "./CalculatorEditor";
import { InjuryLevelsEditor } from "./InjuryLevelsEditor";
import { EFFECT_TO_LEVEL_TYPE, EFFECT_TYPES, EffectType } from "./LevelParamsEditor";

export const EffectEditor = ({
  name,
  effect,
  onUpdate,
}: {
  name: string;
  effect: AgentEffect;
  onUpdate: (newName?: string, newEffect?: AgentEffect) => void;
}) => {
  const handleTypeChange = (newType: EffectType) => {
    const newLevelType = EFFECT_TO_LEVEL_TYPE[newType];
    onUpdate(name, {
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
    onUpdate(name, {
      ...effect,
      parameters: { ...effect.parameters, levels, parameters },
    } as AgentEffect);
  };

  return (
    <Accordion defaultExpanded={false} disableGutters>
      <AccordionSummary
        expandIcon={<ExpandMore />}
        sx={{
          '& .effect-delete': { visibility: 'hidden' },
          '&:hover .effect-delete': { visibility: 'visible' },
        }}
      >
        <Stack direction="row" alignItems="center" spacing={1} sx={{ flex: 1, mr: 1 }}>
          <Science fontSize="small" color="action" />
          <RenameField
            value={name}
            setValue={(newName) => {
              if (newName && newName !== name) {
                onUpdate(newName, effect);
              }
            }}
          />
          <Chip
            label={effect.type}
            size="small"
            sx={{ bgcolor: 'rgba(56, 189, 248, 0.15)', color: '#38bdf8' }}
          />
          <Chip
            label={getCalculatorType(effect.calculator)}
            size="small"
            variant="outlined"
            sx={{ color: '#2dd4bf', borderColor: '#2dd4bf' }}
          />
        </Stack>
        <IconButton
          className="effect-delete"
          component="div"
          size="small"
          onClick={(e) => {
            e.stopPropagation();
            onUpdate();
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
              onChange={(e) => onUpdate(name, { ...effect, units: e.target.value || undefined } as AgentEffect)}
            />
          </Stack>

          <Divider />
          <SectionHeader>Calculator</SectionHeader>
          <CalculatorEditor
            calculator={effect.calculator}
            setCalculator={(calc: Calculator) => {
              onUpdate(name, { ...effect, calculator: calc } as AgentEffect);
            }}
          />

          <Divider />
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
