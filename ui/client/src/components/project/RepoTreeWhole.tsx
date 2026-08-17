import { Add, EditNote } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { SimpleTreeView } from "@mui/x-tree-view/SimpleTreeView";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { CentralRepoFolder } from "../repo/CentralRepoFolder";
import { RegisteredRepositories } from "../repo/RegisteredRepositories";
import { RepoTreeItem } from "./RepoTreeItem";
import { treeSelectionSx } from "./treeSelectionSx";

const SAMPLE_REPO_PATH = 'hera/path/to/repo.json';

// The repositories section: its own titled tree, separate from the documents tree,
// with no single wrapping root. Selection/expansion state is owned by the parent.
export const RepoTreeWhole = ({
  selectedIds,
  onSelectedItemsChange,
  expandedItems,
  onExpandedItemsChange,
}: {
  selectedIds: string[],
  onSelectedItemsChange: (event: React.SyntheticEvent | null, ids: string[]) => void,
  expandedItems: string[],
  onExpandedItemsChange: (event: React.SyntheticEvent | null, ids: string[]) => void,
}) => {
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

  const addTempRepo = () => {
    let num = 1;
    while (repositories.includes(`${TEMP_REPO_NAME} ${num}`)) num++;
    const name = `${TEMP_REPO_NAME} ${num}`;
    setRepositories([...repositories, name]);
  };

  return (
    <>
      <Stack direction="row" alignItems="center" spacing={0.5} sx={{ mt: 2, mb: 1 }}>
        <Typography
          variant="overline"
          sx={{ color: 'text.secondary', fontWeight: 600, letterSpacing: 1 }}
        >
          Repositories
        </Typography>
        <ButtonTooltip title="Add repository" onClick={addSampleRepo}>
          <Add />
        </ButtonTooltip>
        <ButtonTooltip title="Edit temporary repository" onClick={addTempRepo}>
          <EditNote />
        </ButtonTooltip>
      </Stack>
      <SimpleTreeView
        selectedItems={selectedIds}
        onSelectedItemsChange={onSelectedItemsChange}
        expandedItems={expandedItems}
        onExpandedItemsChange={onExpandedItemsChange}
        expansionTrigger="content"
        multiSelect
        sx={treeSelectionSx}
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
      </SimpleTreeView>
    </>
  );
};
