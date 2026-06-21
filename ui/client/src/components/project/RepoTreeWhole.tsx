import { Add, Delete, EditNote } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { CentralRepoFolder } from "../repo/CentralRepoFolder";
import { RegisteredRepositories } from "../repo/RegisteredRepositories";

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
          <ButtonTooltip
            title={'Edit temporary repository'}
            onClick={() => {
              let num = 1;
              while (repositories.includes(`${TEMP_REPO_NAME} ${num}`)) num++;
              const name = `${TEMP_REPO_NAME} ${num}`;
              setRepositories([...repositories, name]);
            }}
          >
            <EditNote />
          </ButtonTooltip>
        </Stack>
      )}
    >
      <CentralRepoFolder />
      <RegisteredRepositories showUpdateButton />
      {repositories.map(repoPath => (
        <TreeItem
          key={idRepoId(repoPath)}
          itemId={idRepoId(repoPath)}
          label={(
            <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
              {repoPath}
              <ButtonTooltip
                title={'Remove repository from this list'}
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