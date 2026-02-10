import { Stack, TextField } from "@mui/material";
import { SectionHeader } from "../../elements/SectionHeader";
import { AgentConfig, VaporPressure } from "../../shared/AgentConfig";

export const PhysicalPropertiesEditor = ({
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
