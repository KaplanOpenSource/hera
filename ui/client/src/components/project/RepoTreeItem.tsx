import { Delete } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { RenameField } from "../../elements/RenameField";
import { idRepoId } from "../../shared/idDocId";

export const RepoTreeItem = ({
  repoPath,
  onRename,
  onRemove,
  defaultEditing = false,
}: {
  repoPath: string,
  onRename: (newPath: string) => void,
  onRemove: () => void,
  defaultEditing?: boolean,
}) => {
  return (
    <TreeItem
      itemId={idRepoId(repoPath)}
      label={(
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
          <RenameField
            value={repoPath}
            setValue={onRename}
            defaultEditing={defaultEditing}
          />
          <ButtonTooltip
            title={'Remove repository from this list'}
            onClick={onRemove}
          >
            <Delete />
          </ButtonTooltip>
        </Stack>
      )}
    >
    </TreeItem>
  )
}
