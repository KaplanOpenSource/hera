import { Paper, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { ProjectEntire } from '../shared/types';
import { useProjectStore } from '../stores/useProjectStore';
import { ToolkitTreeItem } from './ToolkitTreeItem';
import { useState } from 'react';

export const ProjectTreeView = ({
  project,
}: {
  project: ProjectEntire;
}) => {
  const { toolkits } = useProjectStore();
  const [showEmptyToolkits, setShowEmptyToolkits] = useState(false);

  console.log(toolkits)
  console.log(project)

  return (
    <Paper sx={{ p: 2 }}>
      <SimpleTreeView
        defaultExpandedItems={['project-documents', 'no-toolkit']}
      >
        <TreeItem key={`project-documents`} itemId={`project-documents`}
          label={(<>
            <Typography>
              Project {project.name}
            </Typography>
          </>)}
        >
          {toolkits.map(toolkit => (
            <ToolkitTreeItem
              project={project}
              toolkit={toolkit}
            />
          ))}
          <ToolkitTreeItem
            project={project}
            toolkit={undefined}
          />
        </TreeItem>
      </SimpleTreeView>
    </Paper>
  );
};
