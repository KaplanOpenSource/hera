import { TreeItem } from "@mui/x-tree-view";
import { Stack } from "@mui/material";
import { UpdateRepositoriesButton } from "./UpdateRepositoriesButton";
import { useRegisteredRepositories } from "./useRegisteredRepositories";
import { DEFAULT_PROJECT } from "../../stores/useProjectStore";
import { idRepoId } from "../../shared/idDocId";

export const RegisteredRepositories = ({
  showUpdateButton = false,
}: {
  showUpdateButton?: boolean,
}) => {
  const { repositories } = useRegisteredRepositories();

  return (
    <TreeItem itemId="registered-repos" label={
      <Stack direction='row' spacing={1} alignItems={'center'}>
        <span>
          Repositories in <b>{DEFAULT_PROJECT}</b>
        </span>
        {showUpdateButton && <UpdateRepositoriesButton />}
      </Stack>
    }>
      {repositories.length === 0
        ? <TreeItem itemId="registered-repos/__empty" label="No registered repositories" />
        : repositories.map(r => (
          <TreeItem
            key={r.datasourceName}
            itemId={idRepoId(r.datasourceName)}
            label={r.resource}
          />
        ))
      }
    </TreeItem>
  )
}
