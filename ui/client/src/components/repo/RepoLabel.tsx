import { Typography, Stack } from "@mui/material";
import { Repository } from "../../shared/types";

export const RepoLabel = ({ repo }: { repo: Repository }) => (
  <Stack direction="column" spacing={0.5} sx={{ py: 0.5 }}>
    <Stack direction='row' spacing={0.5}>
      <Typography
        variant="body2"
        fontWeight="bold"
      >
        {repo.datasourceName}
      </Typography>
      <Typography variant="caption" color="text.secondary">{repo.version.join('.')}</Typography>
    </Stack>
    <Typography variant="caption" color="text.secondary">Resource: {repo.resource}</Typography>
    {/* <Typography variant="caption" color="text.secondary">Data Format: {repo.dataFormat}</Typography> */}
    <Typography variant="caption" color="text.secondary">Toolkit: {repo.toolkit}</Typography>
  </Stack>
)
