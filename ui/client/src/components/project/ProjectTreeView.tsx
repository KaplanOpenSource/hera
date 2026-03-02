import { Folder, Refresh } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { ProjectObj } from '../../objects/ProjectObj';
import { useProjectStore } from '../../stores/useProjectStore';
import { AddDocumentButton } from './AddDocumentButton';
import { DocumentSplitGroup } from './DocumentSplitGroup';
import { ProjectViewSettingsButton } from './ProjectViewSettingsButton';
import { RepoTreeWhole } from './RepoTreeWhole';
import { useViewSettingsStore } from '../../stores/useViewSettingsStore';

export type ViewSettingsType = {
  minGroupSize: number;
  maxDepth: number;
  maxBranches: number;
  showEmptyToolkits: boolean;
  showDocumentPreview: boolean;
};

export const ProjectTreeView = ({
  project,
  selectedItemsIds = [],
  setSelectedItemIds,
}: {
  project: ProjectObj;
  selectedItemsIds?: string[];
  setSelectedItemIds: (v: string[]) => void,
}) => {
  const { toolkits } = useProjectStore();
  const { viewSettings } = useViewSettingsStore();

  console.log(toolkits)
  console.log(project)

  return (
    <SimpleTreeView
      defaultExpandedItems={['project-documents', 'no-toolkit', '*repos*']}
      selectedItems={selectedItemsIds[0] ?? null}
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
              title={'Reload documents'}
              onClick={() => fetchProjectDetails(project.name)}
            >
              <Refresh />
            </ButtonTooltip>
            <ProjectViewSettingsButton
            />
            <AddDocumentButton
              toolkit={undefined}
            />
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
        <DocumentSplitGroup
          docs={project?.documents}
          project={project}
          depth={viewSettings.maxDepth}
        />
      </TreeItem>
      <RepoTreeWhole
      />
    </SimpleTreeView>
  );
};
