import { Delete, Edit } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
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
  const [editing, setEditing] = useState(defaultEditing);

  return (
    <TreeItem
      itemId={idRepoId(repoPath)}
      label={(
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
          {editing
            ? (
              <RenameField
                value={repoPath}
                setValue={newPath => {
                  onRename(newPath);
                  setEditing(false);
                }}
              />
            )
            : (<>
              {repoPath}
              <ButtonTooltip
                title={'Rename repository address'}
                onClick={() => setEditing(true)}
              >
                <Edit />
              </ButtonTooltip>
              <ButtonTooltip
                title={'Remove repository from this list'}
                onClick={onRemove}
              >
                <Delete />
              </ButtonTooltip>
            </>)}
        </Stack>
      )}
    >
    </TreeItem>
  )
}
