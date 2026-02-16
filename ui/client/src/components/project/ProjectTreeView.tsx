import { Folder, Preview, Refresh, Settings, Visibility } from '@mui/icons-material';
import { Stack, TextField, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useDialog } from '../../elements/useDialog';
import { fetchProjectDetails } from '../../io/FetchProjects';
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
  const [showDocumentPreview, setShowDocumentPreview] = useState(true);
  const [minGroupSize, setMinGroupSize] = useState(2);
  const [maxDepth, setMaxDepth] = useState(5);
  const { DialogComponent, openDialog } = useDialog();

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
            <ButtonTooltip
              title={'Change tree settings'}
              onClick={async () => {
                const result = await openDialog({
                  title: 'Change tree settings',
                  initialValues: { maxDepth, minGroupSize },
                  maxWidth: 'sm',
                  render:
                    ({ values, setValues }) => (
                      <Stack direction={'column'} spacing={1}>
                        <TextField
                          label="Max Depth"
                          value={values.maxDepth}
                          type='number'
                          variant='outlined'
                          slotProps={{ htmlInput: { min: 0, max: 25 } }}
                          size='small'
                          onChange={(e) => setValues({ ...values, maxDepth: e.target.value })}
                          onClick={e => e.stopPropagation()}
                        />
                        <TextField
                          label="Min Group Size"
                          value={values.minGroupSize}
                          type='number'
                          variant='outlined'
                          slotProps={{ htmlInput: { min: 1 } }}
                          size='small'
                          onChange={(e) => setValues({ ...values, minGroupSize: e.target.value })}
                          onClick={e => e.stopPropagation()}
                        />
                      </Stack>
                    )
                });

                if (result.confirmed && result.values) {
                  setMaxDepth(result.values.maxDepth)
                  setMinGroupSize(result.values.minGroupSize)
                }
              }}
            >
              <Settings />
              {DialogComponent}
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
            minGroupSize={minGroupSize}
            maxDepth={maxDepth}
          />
        ))}
        <ToolkitTreeItem
          key={'no_toolkit'}
          project={project}
          toolkit={undefined}
          showEmpty={showEmptyToolkits}
          showDocumentPreview={showDocumentPreview}
          minGroupSize={minGroupSize}
          maxDepth={maxDepth}
        />
      </TreeItem>
      <RepoTreeWhole
      />
    </SimpleTreeView>
  );
};
