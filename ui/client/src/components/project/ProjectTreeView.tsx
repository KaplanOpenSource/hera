import { Folder, Visibility, VisibilityOff } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ProjectObj } from '../../objects/ProjectObj';
import { useProjectStore } from '../../stores/useProjectStore';
import { RepoTreeWhole } from './RepoTreeWhole';
import { ToolkitTreeItem } from './ToolkitTreeItem';

export const ProjectTreeView = ({
  project,
  setSelectedItemIds,
}: {
  project: ProjectObj;
  setSelectedItemIds: (v: string[]) => void,
}) => {
  const { toolkits } = useProjectStore();
  const [showEmptyToolkits, setShowEmptyToolkits] = useState(false);

  console.log(toolkits)
  console.log(project)

  return (
    <SimpleTreeView
      defaultExpandedItems={['project-documents', 'no-toolkit', '*repos*']}
      onSelectedItemsChange={(_e, itemIds) => {
        setSelectedItemIds(itemIds ? [itemIds] : [])
      }}
      expansionTrigger={'iconContainer'}
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
          </Stack>
        )}
      >
        <Tooltip title='Files directory where the project is located'>
          <Stack direction={'row'} spacing={1} alignItems={'center'} justifyContent={'start'} sx={{ marginLeft: 5, width: 'fit-content' }}>
            <Folder />
            <Typography>
              {project.configDocument?.data.desc.filesDirectory}
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
      <RepoTreeWhole
      />
    </SimpleTreeView>
  );
};
