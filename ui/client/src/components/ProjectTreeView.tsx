import { Paper } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { ProjectEntire } from '../shared/types';
import { useProjectStore } from '../stores/useProjectStore';
import { ProjectDocumentItem } from './ProjectDocumentItem';
import { ToolkitTreeItem } from './ToolkitTreeItem';

export const ProjectTreeView = ({
  project,
}: {
  project: ProjectEntire | null;
}) => {
  const { toolkits } = useProjectStore();

  console.log(toolkits)
  console.log(project)

  const documentsWithoutToolkit = () => {
    return project?.documents.filter(d => {
      const found = toolkits.find(({ toolkit }) => toolkit === d?.desc?.toolkit);
      return found === undefined;
    });
  }

  return (
    <Paper sx={{ p: 2 }}>
      <SimpleTreeView
        defaultExpandedItems={['project-documents', 'no-toolkit']}
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
              {toolkits.map(toolkit => (
                <ToolkitTreeItem
                  project={project}
                  toolkit={toolkit}
                />
              ))}
            </TreeItem>
          )
        }
      </SimpleTreeView>
    </Paper>
  );
};
