import { Delete } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { idRepoId } from "../../shared/idDocId";

export const RepoTreeOne = ({
  repoPath,
  setRepoPath,
}: {
  repoPath: string,
  setRepoPath: (v: string | undefined) => void,
}) => {
  return (
    <TreeItem
      key={idRepoId(repoPath)}
      itemId={idRepoId(repoPath)}
      label={(
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
          {repoPath}
          <ButtonTooltip
            title={'Remove repository'}
            onClick={() => setRepoPath(undefined)}
          >
            <Delete />
          </ButtonTooltip>
        </Stack>
      )}
    >
    </TreeItem>
  )
}