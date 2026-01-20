import { Add, Delete } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";
import { idRepoId } from "../../shared/idDocId";

export const RepoTreeWhole = ({ }) => {
  const [repositories, setRepositories] = useState<string[]>(['hera/doc/jupyter/Developer/Documentation_Repository.json']);
  const { confirmOpen, ConfirmDialog } = useConfirm();
  return (
    <TreeItem key={'*repos*'} itemId={'*repos*'}
      label={(
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
          <Typography>
            Repositories
          </Typography>
          <ButtonTooltip
            title={'Add repository'}
            onClick={async (e) => {
              const { confirmed, text } = await confirmOpen({
                title: `Add repository`,
                requireText: true,
                textLabel: 'Repository location',
              });
              if (confirmed && text) {
                setRepositories([...repositories, text]);
              }
            }}
          >
            <Add />
          </ButtonTooltip>
          {ConfirmDialog}
        </Stack>
      )}
    >
      {repositories.map(repoPath => (
        <TreeItem
          key={idRepoId(repoPath)}
          itemId={idRepoId(repoPath)}
          label={(
            <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
              {repoPath}
              <ButtonTooltip
                title={'Remove repository'}
                onClick={() => {
                  setRepositories(repositories.filter(x => x !== repoPath))
                }}
              >
                <Delete />
              </ButtonTooltip>
            </Stack>
          )}
        >
        </TreeItem>
      ))}
    </TreeItem>
  )
}