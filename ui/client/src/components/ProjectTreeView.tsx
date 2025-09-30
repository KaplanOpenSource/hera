import { Paper } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { Project } from '../shared/types';
import { ProjectTreeItem } from './ProjectTreeItem';

export const ProjectTreeView = ({
  projects,
  onProjectSelect,
  onProjectExpand,
}: {
  projects: Project[];
  onProjectSelect: (project: Project) => void;
  onProjectExpand: (projectId: string) => void;
}) => {
  const handleItemSelection = (event: React.SyntheticEvent | null, itemId: string | null) => {
    if (itemId) {
      const selectedProject = projects.find((p) => p.id === itemId);
      if (selectedProject) {
        onProjectSelect(selectedProject);
      }
    }
  };

  const handleItemExpansion = (event: React.SyntheticEvent | null, itemIds: string[]) => {
    const lastExpanded = itemIds[0];
    if (lastExpanded) {
      onProjectExpand(lastExpanded);
    }
  };

  return (
    <Paper sx={{ p: 2 }}>
      <SimpleTreeView onSelectedItemsChange={handleItemSelection} onExpandedItemsChange={handleItemExpansion}>
        {projects.map((project) => (
          <ProjectTreeItem key={project.id} project={project} />
        ))}
      </SimpleTreeView>
    </Paper>
  );
};
