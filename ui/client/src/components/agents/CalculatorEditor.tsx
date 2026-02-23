import { FormControl, InputLabel, MenuItem, Select, Stack, TextField } from "@mui/material";
import { Calculator } from "../../shared/AgentConfig";

type CalculatorType = "Haber" | "TenBerge" | "MaxConcentration";

const CALCULATOR_TYPES: CalculatorType[] = ["Haber", "TenBerge", "MaxConcentration"];

export function getCalculatorType(calc: Calculator): CalculatorType {
  if ("Haber" in calc) return "Haber";
  if ("TenBerge" in calc) return "TenBerge";
  return "MaxConcentration";
}

function getCalculatorParams(calc: Calculator) {
  if ("Haber" in calc) return calc.Haber;
  if ("TenBerge" in calc) return calc.TenBerge;
  return calc.MaxConcentration;
}

export function makeDefaultCalculator(type: CalculatorType): Calculator {
  switch (type) {
    case "Haber":
      return { Haber: { breathingRate: 10 } };
    case "TenBerge":
      return { TenBerge: { breathingRate: 10 } };
    case "MaxConcentration":
      return { MaxConcentration: { sampling: "10min", breathingRate: 10 } };
  }
}

export const CalculatorEditor = ({
  calculator,
  setCalculator,
}: {
  calculator: Calculator;
  setCalculator: (newCalc: Calculator) => void;
}) => {
  const type = getCalculatorType(calculator);
  const params = getCalculatorParams(calculator);

  const handleTypeChange = (newType: CalculatorType) => {
    setCalculator(makeDefaultCalculator(newType));
  };

  const updateParam = (key: string, value: any) => {
    const newParams = { ...params, [key]: value };
    switch (type) {
      case "Haber":
        setCalculator({ Haber: newParams as any });
        break;
      case "TenBerge":
        setCalculator({ TenBerge: newParams as any });
        break;
      case "MaxConcentration":
        setCalculator({ MaxConcentration: newParams as any });
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

      {type === "MaxConcentration" && (
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
