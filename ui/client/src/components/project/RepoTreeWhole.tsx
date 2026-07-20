import { Add, EditNote } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { CentralRepoFolder } from "../repo/CentralRepoFolder";
import { RegisteredRepositories } from "../repo/RegisteredRepositories";
import { RepoTreeItem } from "./RepoTreeItem";

const SAMPLE_REPO_PATH = 'hera/path/to/repo.json';

export const RepoTreeWhole = ({ }) => {
  const [repositories, setRepositories] = useState<string[]>(['hera/doc/jupyter/Developer/Documentation_Repository.json']);
  const [newRepo, setNewRepo] = useState<string | null>(null);

  const addSampleRepo = () => {
    let sample = SAMPLE_REPO_PATH;
    let num = 1;
    while (repositories.includes(sample)) {
      sample = SAMPLE_REPO_PATH.replace('.json', `-${num}.json`);
      num++;
    }
    setNewRepo(sample);
    setRepositories([...repositories, sample]);
  };

  return (
    <TreeItem key={'*repos*'} itemId={'*repos*'}
      label={(
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
          <Typography>
            Repositories
          </Typography>
          <ButtonTooltip
            title={'Add repository'}
            onClick={addSampleRepo}
          >
            <Add />
          </ButtonTooltip>
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
        <RepoTreeItem
          key={idRepoId(repoPath)}
          repoPath={repoPath}
          defaultEditing={repoPath === newRepo}
          onRename={newPath => setRepositories(repositories.map(x => x === repoPath ? newPath : x))}
          onRemove={() => setRepositories(repositories.filter(x => x !== repoPath))}
        />
      ))}
    </TreeItem>
  )
}