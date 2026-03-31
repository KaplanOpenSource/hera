import { Folder, Refresh } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useCallback, useEffect, useRef, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { ProjectObj } from '../../objects/ProjectObj';
import { idDocId, idFromDocId } from '../../shared/idDocId';
import { useProjectStore } from '../../stores/useProjectStore';
import { SplitTree } from '../../utils/splitTree';
import { AddDocumentButton } from './AddDocumentButton';
import { DocumentSplitGroup } from './DocumentSplitGroup';
import { ProjectViewSettingsButton } from './ProjectViewSettingsButton';
import { RegisteredRepositories } from '../repo/RegisteredRepositories';
import { NotebookTreeItem } from './NotebookTreeItem';
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
  const [expandedItems, setExpandedItems] = useState<string[]>(['project-documents', 'no-toolkit', '*repos*', 'registered-repos']);
  const splitTreeRef = useRef<SplitTree | null>(null);

  const getSplitTree = useCallback(() => {
    const currentProject = useProjectStore.getState().getProject();
    const currentSettings = useViewSettingsStore.getState().viewSettings;
    if (!currentProject) return null;
    if (!splitTreeRef.current || splitTreeRef.current.needsRebuild(currentProject.documents, currentSettings)) {
      splitTreeRef.current = new SplitTree(currentProject.documents, currentSettings.maxDepth, currentSettings);
    }
    return splitTreeRef.current;
  }, []);

  const expandToDocument = useCallback((docOid: string) => {
    const tree = getSplitTree();
    if (!tree) return;
    const ancestors = tree.findAncestorKeys(docOid);
    console.log('[focus] docOid:', docOid, 'ancestors:', ancestors);
    if (ancestors) {
      setExpandedItems(prev => [...new Set([...prev, 'project-documents', ...ancestors])]);
    }
  }, [getSplitTree]);

  const handleDocumentCreated = useCallback((docOid: string) => {
    expandToDocument(docOid);
    setSelectedItemIds([idDocId(docOid)]);
  }, [expandToDocument, setSelectedItemIds]);

  // Expand branches to the initially selected document (e.g. from URL)
  const hasExpandedInitial = useRef(false);
  useEffect(() => {
    if (hasExpandedInitial.current) return;
    const docOid = idFromDocId(selectedItemsIds[0]);
    if (docOid && project?.documentIds.has(docOid)) {
      console.log('[expand-initial] expanding to', docOid, 'documents:', project?.documentIds);
      expandToDocument(docOid);
      hasExpandedInitial.current = true;
    }
  }, [project, selectedItemsIds, expandToDocument]);

  console.log(toolkits)
  console.log(project)

  return (
    <SimpleTreeView
      expandedItems={expandedItems}
      onExpandedItemsChange={(_e, itemIds) => setExpandedItems(itemIds)}
      selectedItems={selectedItemsIds[0] ?? null}
      onSelectedItemsChange={(_e, itemIds) => {
        setSelectedItemIds(itemIds ? [itemIds] : [])
      }}
      expansionTrigger={'content'}
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
              onDocumentCreated={handleDocumentCreated}
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
          onDocumentDeleted={() => setSelectedItemIds([])}
        />
      </TreeItem>
      <NotebookTreeItem projectName={project.name} />
      <RegisteredRepositories showUpdateButton />
      <RepoTreeWhole
      />
    </SimpleTreeView>
  );
};
