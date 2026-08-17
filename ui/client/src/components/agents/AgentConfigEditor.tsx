import {
  Box,
  Divider,
  Stack,
  TextField,
  Typography
} from "@mui/material";
import type { AgentConfig } from "../../shared/AgentConfig";
import { EffectsListEditor } from "./EffectsListEditor";
import { PhysicalPropertiesEditor } from "./PhysicalPropertiesEditor";

export const AgentConfigEditor = ({
  agentResource,
  setAgentResource,
}: {
  agentResource: AgentConfig;
  setAgentResource: (agent: AgentConfig) => void;
}) => {
  return (
    <Box sx={{ maxWidth: 800 }}>
      <EffectsListEditor
        effects={agentResource.effects}
        onChange={(effects) => setAgentResource({ ...agentResource, effects })}
      />

      <Divider sx={{ my: 2 }} />

      <TextField
        label="Ten Berge Coefficient (exponent n)"
        type="number"
        size="small"
        slotProps={{ htmlInput: { step: 0.1 } }}
        helperText="Global exponent n used by Ten Berge toxic gas exposure calculators"
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
