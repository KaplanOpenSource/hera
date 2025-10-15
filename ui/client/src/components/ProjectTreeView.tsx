import { Paper } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { Project, ProjectEntire } from '../shared/types';
import { ProjectTreeItem } from './ProjectTreeItem';
import { TreeItem } from '@mui/x-tree-view';

export const ProjectTreeView = ({
  project,
  // onProjectSelect,
  // onProjectExpand,
}: {
  project: ProjectEntire | null;
  // onProjectSelect: (project: Project) => void;
  // onProjectExpand: (projectId: string) => void;
}) => {
  // const handleItemSelection = (event: React.SyntheticEvent | null, itemId: string | null) => {
  //   if (itemId) {
  //     const selectedProject = projects.find((p) => p.id === itemId);
  //     if (selectedProject) {
  //       onProjectSelect(selectedProject);
  //     }
  //   }
  // };

  // const handleItemExpansion = (event: React.SyntheticEvent | null, itemIds: string[]) => {
  //   const lastExpanded = itemIds[0];
  //   if (lastExpanded) {
  //     onProjectExpand(lastExpanded);
  //   }
  // };

  return (
    <Paper sx={{ p: 2 }}>
      <SimpleTreeView
      // onSelectedItemsChange={handleItemSelection}
      // onExpandedItemsChange={handleItemExpansion}
      >
        {!project
          ? (
            <TreeItem key={`no-project`} itemId={`no-project`} label="No project loaded" disabled />
          ) : (
            <TreeItem key={`project-documents`} itemId={`project-documents`} label={`Project ${project.name}`}>
              {project.documents.map(d => (
                <TreeItem
                  key={`document${d.desc.docid}`} itemId={`document${d.desc.docid}`} label={`Document: ${d.desc.datasourceName}`}
                >

                </TreeItem>
              ))}
            </TreeItem>
          )
        }
        {/* {projects.map((project) => (
          <ProjectTreeItem key={project.id} project={project} />
        ))} */}
      </SimpleTreeView>
    </Paper>
  );
};
