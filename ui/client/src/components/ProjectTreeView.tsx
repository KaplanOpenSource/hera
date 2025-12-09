import { IconButton, Paper, Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { ProjectEntire } from '../shared/types';
import { useProjectStore } from '../stores/useProjectStore';
import { ToolkitTreeItem } from './ToolkitTreeItem';
import { useState } from 'react';
import { Visibility, VisibilityOff } from '@mui/icons-material';
import { ButtonTooltip } from '../elements/ButtonTooltip';

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
          label={(
            <Stack direction='row' spacing={1} justifyContent="start" alignItems='center'>
              <Typography>
                Project {project.name}
              </Typography>
              <ButtonTooltip
                title={showEmptyToolkits ? 'Showing empty toolkits, click to hide' : 'Hiding empty toolkits, click to show'}
                onClick={() => setShowEmptyToolkits(!showEmptyToolkits)}
              >
                {showEmptyToolkits ? <Visibility /> : <VisibilityOff />}
              </ButtonTooltip>
            </Stack>
          )}
        >
          {toolkits.map(toolkit => (
            <ToolkitTreeItem
              project={project}
              toolkit={toolkit}
              showEmpty={showEmptyToolkits}
            />
          ))}
          <ToolkitTreeItem
            project={project}
            toolkit={undefined}
            showEmpty={showEmptyToolkits}
          />
        </TreeItem>
      </SimpleTreeView>
    </Paper>
  );
};
