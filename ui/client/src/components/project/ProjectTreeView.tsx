import { Folder, Visibility, VisibilityOff } from '@mui/icons-material';
import { Box, Paper, Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import type { ProjectEntire } from '../../shared/types';
import { useProjectStore } from '../../stores/useProjectStore';
import { ToolkitTreeItem } from './ToolkitTreeItem';
import { RepoAddButton } from '../repo/RepoAddButton';

export const ProjectTreeView = ({
  project,
  setSelectedItemIds,
}: {
  project: ProjectEntire;
  setSelectedItemIds: (v: string[]) => void,
}) => {
  const { toolkits } = useProjectStore();
  const [showEmptyToolkits, setShowEmptyToolkits] = useState(false);

  console.log(toolkits)
  console.log(project)

  return (
    <Paper sx={{ p: 2 }}>
      <SimpleTreeView
        defaultExpandedItems={['project-documents', 'no-toolkit']}
        onSelectedItemsChange={(e, itemIds) => {
          setSelectedItemIds(itemIds ? [itemIds] : [])
        }}
        multiSelect={false}
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
              <RepoAddButton
              />
            </Stack>
          )}
        >
          <Tooltip title='Files directory where the project is located'>
            <Stack direction={'row'} spacing={1} alignItems={'center'} justifyContent={'start'} sx={{ marginLeft: 5, width: 'fit-content' }}>
              <Folder />
              <Typography>
                {project.documents.find(d => d.type === project.name + '__config__')?.desc.filesDirectory}
              </Typography>
            </Stack>
          </Tooltip>
          {toolkits.map(toolkit => (
            <ToolkitTreeItem
              key={toolkit.toolkit}
              project={project}
              toolkit={toolkit}
              showEmpty={showEmptyToolkits}
            />
          ))}
          <ToolkitTreeItem
            key={'no_toolkit'}
            project={project}
            toolkit={undefined}
            showEmpty={showEmptyToolkits}
          />
        </TreeItem>
      </SimpleTreeView>
    </Paper>
  );
};
