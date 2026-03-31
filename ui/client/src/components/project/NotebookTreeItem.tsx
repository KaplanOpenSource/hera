import { Description } from '@mui/icons-material';
import { TreeItem } from '@mui/x-tree-view';
import { idNotebookId } from '../../shared/idDocId';

export const NotebookTreeItem = ({
  projectName,
}: {
  projectName: string,
}) => {
  return (
    <TreeItem
      itemId={idNotebookId(projectName)}
      label="Notebook"
      slots={{ icon: () => <Description fontSize="small" /> }}
    />
  );
};
