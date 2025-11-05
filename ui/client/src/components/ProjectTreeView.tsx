import { Paper } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { ProjectEntire } from '../shared/types';
import { useProjectStore } from '../stores/useProjectStore';
import { ProjectDocumentItem } from './ProjectDocumentItem';

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

  const { toolkits } = useProjectStore();

  console.log(toolkits)
  console.log(project)

  const documentsForToolkit = (toolkitName: string) => {
    return project?.documents.filter(d => d?.desc?.toolkit === toolkitName) || [];
  }

  const documentsWithoutToolkit = () => {
    return project?.documents.filter(d => {
      const found = toolkits.find(({ toolkit }) => toolkit === d?.desc?.toolkit);
      return found === undefined;
    });
  }

  return (
    <Paper sx={{ p: 2 }}>
      <SimpleTreeView
      >
        {!project
          ? (
            <TreeItem key={`no-project`} itemId={`no-project`} label="No project loaded" disabled />
          ) : (
            <TreeItem key={`project-documents`} itemId={`project-documents`}
              label={`Project ${project.name}`}
            >
              <TreeItem key={'no-toolkit'} itemId={'no-toolkit'} label='Without Toolkit'>
                {documentsWithoutToolkit()?.map(d => {
                  return (
                    <ProjectDocumentItem
                      key={`proj${project.name}_doc${d._id.$oid}`} project={project} document={d}
                    >
                    </ProjectDocumentItem>
                  )
                })}
              </TreeItem>
              {toolkits.map(({ toolkit, cls, description }) => {
                const docs = documentsForToolkit(toolkit);
                return !docs.length ? null : (
                  <TreeItem key={toolkit} itemId={toolkit} label={toolkit}>
                    {/* <Typography>
                    {cls}
                  </Typography>
                  {desc ?? (
                    <Typography>
                      {desc}
                    </Typography>
                  )} */}
                    {docs.map(d => {
                      return (
                        <ProjectDocumentItem
                          key={`proj${project.name}_doc${d._id.$oid}`} project={project} document={d}
                        >
                        </ProjectDocumentItem>
                      )
                    })}
                  </TreeItem>
                )
              })}
            </TreeItem>
          )
        }
      </SimpleTreeView>
    </Paper>
  );
};
