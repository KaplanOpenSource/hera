import { Folder, Preview, Refresh, Visibility, VisibilityOff } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ProjectObj } from '../../objects/ProjectObj';
import { useProjectStore } from '../../stores/useProjectStore';
import { RepoTreeWhole } from './RepoTreeWhole';
import { ToolkitTreeItem } from './ToolkitTreeItem';
import { fetchProjectDetails } from '../../io/FetchProjects';

export const ProjectTreeView = ({
  project,
  setSelectedItemIds,
}: {
  project: ProjectObj;
  setSelectedItemIds: (v: string[]) => void,
}) => {
  const { toolkits } = useProjectStore();
  const [showEmptyToolkits, setShowEmptyToolkits] = useState(false);
  const [showDocumentPreview, setShowDocumentPreview] = useState(true);

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
          <Stack direction='row' justifyContent="start" alignItems='center'>
            <Typography marginRight={1}>
              Project {project.name}
            </Typography>
            <ButtonTooltip
              title={showEmptyToolkits ? 'Showing empty toolkits, click to hide' : 'Hiding empty toolkits, click to show'}
              onClick={() => setShowEmptyToolkits(!showEmptyToolkits)}
              color={showEmptyToolkits ? 'primary' : 'default'}
            >
              <Visibility />
            </ButtonTooltip>
            <ButtonTooltip
              title={showDocumentPreview ? 'Showing document preview, click to hide' : 'Hiding document preview, click to show'}
              onClick={() => setShowDocumentPreview(!showDocumentPreview)}
              color={showDocumentPreview ? 'primary' : 'default'}
            >
              <Preview />
            </ButtonTooltip>
            <ButtonTooltip
              title={'Reload documents'}
              onClick={() => fetchProjectDetails(project.name)}
            >
              <Refresh />
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
            showDocumentPreview={showDocumentPreview}
            minGroupSize={2}
            maxDepth={5}
          />
        ))}
        <ToolkitTreeItem
          key={'no_toolkit'}
          project={project}
          toolkit={undefined}
          showEmpty={showEmptyToolkits}
          showDocumentPreview={showDocumentPreview}
        />
      </TreeItem>
      <RepoTreeWhole
      />
    </SimpleTreeView>
  );
};
