import { TextField, Stack } from "@mui/material";
import { AgentEffect } from "../../shared/AgentConfig";

export type InjuryLevelType = "Threshold" | "Lognormal10DoseResponse" | "Exponential";

export type EffectType = AgentEffect["type"];

export const EFFECT_TYPES: EffectType[] = ["Threshold", "Lognormal10", "Exponential"];

export const EFFECT_TO_LEVEL_TYPE: { [key in EffectType]: InjuryLevelType } = {
  Threshold: "Threshold",
  Lognormal10: "Lognormal10DoseResponse",
  Exponential: "Exponential",
};

export function makeDefaultLevelParams(levelType: InjuryLevelType) {
  switch (levelType) {
    case "Threshold":
      return { threshold: "1*mg/m**3" };
    case "Lognormal10DoseResponse":
      return { TL_50: 1, sigma: 0.5 };
    case "Exponential":
      return { k: 1 };
  }
}

export const LevelParamsEditor = ({
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
