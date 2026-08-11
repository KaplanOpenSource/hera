import { Delete } from "@mui/icons-material";
import { Box, Stack } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { RenamePathField } from "../../elements/RenamePathField";
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
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'} sx={{ minWidth: 0 }}>
          <RenamePathField
            value={repoPath}
            setValue={onRename}
            defaultEditing={defaultEditing}
          />
          <Box sx={{ flexShrink: 0 }}>
            <ButtonTooltip
              title={'Remove repository from this list'}
              onClick={onRemove}
            >
              <Delete />
            </ButtonTooltip>
          </Box>
        </Stack>
      )}
    >
    </TreeItem>
  )
}
