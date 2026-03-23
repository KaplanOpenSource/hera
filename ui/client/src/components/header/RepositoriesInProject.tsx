import { Typography, Stack } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { Repository } from "../../shared/types";

const RepoLabel = ({ repo }: { repo: Repository }) => (
  <Stack direction="column" spacing={0.5} sx={{ py: 0.5 }}>
    <Stack direction='row' spacing={0.5}>
      <Typography variant="body2" fontWeight="bold">{repo.datasourceName}</Typography>
      <Typography variant="caption" color="text.secondary">{repo.version.join('.')}</Typography>
    </Stack>
    <Typography variant="caption" color="text.secondary">Resource: {repo.resource}</Typography>
    {/* <Typography variant="caption" color="text.secondary">Data Format: {repo.dataFormat}</Typography> */}
    <Typography variant="caption" color="text.secondary">Toolkit: {repo.toolkit}</Typography>
  </Stack>
)

export const RepositoriesInProject = ({
  repositories,
}: {
  repositories: Repository[],
}) => {
  if (repositories.length === 0) {
    return null;
  }
  return (
    <SimpleTreeView>
      {repositories.map(r => (
        <TreeItem key={r.datasourceName} itemId={r.datasourceName}
          label={<RepoLabel repo={r} />}
        />
      ))}
    </SimpleTreeView>
  )
}
