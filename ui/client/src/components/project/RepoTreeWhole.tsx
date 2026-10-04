import { Add, EditNote } from "@mui/icons-material";
import { Box, Stack, Typography } from "@mui/material";
import { SimpleTreeView } from "@mui/x-tree-view/SimpleTreeView";
import { useRef, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";
import { loadRepositoryIntoProject } from "../../io/loadRepositoryIntoProject";
import { idFromRepoDomId, idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { clickedRow } from "../../utils/clickedRow";
import { CentralRepoFolder } from "../repo/CentralRepoFolder";
import { RegisteredRepositories } from "../repo/RegisteredRepositories";
import { RepoContextMenu } from "./RepoContextMenu";
import { RepoTreeItem } from "./RepoTreeItem";
import { MenuPosition } from "./TreeContextMenu";
import { treeSelectionSx } from "./treeSelectionSx";

const SAMPLE_REPO_PATH = 'hera/path/to/repo.json';

// The repositories section: its own titled tree, separate from the documents tree,
// with no single wrapping root. Selection/expansion state is owned by the parent.
export const RepoTreeWhole = ({
  selectedIds,
  onSelectedItemsChange,
  expandedItems,
  onExpandedItemsChange,
  onOpenItem,
}: {
  selectedIds: string[],
  onSelectedItemsChange: (event: React.SyntheticEvent | null, ids: string[]) => void,
  expandedItems: string[],
  onExpandedItemsChange: (event: React.SyntheticEvent | null, ids: string[]) => void,
  // Opens a row that is not a repository (the central folder, a repository field).
  onOpenItem: (rawId: string | undefined) => void,
}) => {
  const [repositories, setRepositories] = useState<string[]>(['hera/doc/jupyter/Developer/Documentation_Repository.json']);
  const [newRepo, setNewRepo] = useState<string | null>(null);
  const [menuPosition, setMenuPosition] = useState<MenuPosition | null>(null);
  const [menuRepoName, setMenuRepoName] = useState<string | undefined>(undefined);
  // The row clicked last, so a double click knows what it opens.
  const lastClickedRef = useRef<string | undefined>(undefined);
  const { confirmOpen, ConfirmDialog } = useConfirm();

  // Asks, then pulls the repository contents into the project.
  const askLoadRepository = async (repoName: string) => {
    const { confirmed } = await confirmOpen({
      title: `Load repository "${repoName}" into project?`,
    });
    if (!confirmed) return;
    await loadRepositoryIntoProject(repoName);
  };

  const handleSelectedItemsChange = (event: React.SyntheticEvent | null, ids: string[]) => {
    lastClickedRef.current = clickedRow(ids, selectedIds);
    onSelectedItemsChange(event, ids);
  };

  // The repository of the row under the mouse. Its id is only in the DOM, so read it there.
  const repoNameOfEvent = (event: React.MouseEvent) => {
    const target = event.target as HTMLElement | null;
    const rowId = target?.closest('[role="treeitem"]')?.id;
    return rowId ? idFromRepoDomId(rowId) : undefined;
  };

  // A double click loads a repository, or opens any other row.
  const handleDoubleClick = async (event: React.MouseEvent) => {
    const repoName = repoNameOfEvent(event);
    if (!repoName) {
      onOpenItem(lastClickedRef.current);
      return;
    }
    await askLoadRepository(repoName);
  };

  const handleContextMenu = (event: React.MouseEvent) => {
    const repoName = repoNameOfEvent(event);
    if (!repoName) return;
    event.preventDefault();
    setMenuRepoName(repoName);
    setMenuPosition({ x: event.clientX, y: event.clientY });
  };

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
    <Box onDoubleClick={handleDoubleClick} onContextMenu={handleContextMenu}>
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
        onSelectedItemsChange={handleSelectedItemsChange}
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
      <RepoContextMenu
        position={menuPosition}
        repoName={menuRepoName}
        onClose={() => setMenuPosition(null)}
        onLoad={askLoadRepository}
      />
      {ConfirmDialog}
    </Box>
  );
};
