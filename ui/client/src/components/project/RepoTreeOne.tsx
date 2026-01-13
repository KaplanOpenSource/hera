import { Add, Delete } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";

export const RepoTreeOne = ({
  repoPath,
  setRepoPath,
}: {
  repoPath: string,
  setRepoPath: (v: string | undefined) => void,
}) => {
  return (
    <TreeItem
      key={'__repo_*_' + repoPath}
      itemId={'__repo_*_' + repoPath}
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