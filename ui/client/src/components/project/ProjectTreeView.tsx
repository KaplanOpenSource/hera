import { ContentCopy, Folder, Refresh } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useCallback, useEffect, useRef, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idDocId, idFromDocId } from '../../shared/idDocId';
import { useProjectStore } from '../../stores/useProjectStore';
import { SplitTree } from '../../utils/splitTree';
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
  onSelectItem,
}: {
  project: ProjectObj;
  onSelectItem: (rawId: string | undefined) => void,
}) => {
  const { docId } = useParams<{ docId: string }>();
  const navigate = useNavigate();
  const { toolkits } = useProjectStore();
  const { viewSettings } = useViewSettingsStore();
  const [selectedId, setSelectedId] = useState<string | null>(docId ? idDocId(docId) : null);
  const [expandedItems, setExpandedItems] = useState<string[]>(['project-documents', 'no-toolkit', '*repos*']);
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

  // Select an item: highlight it, sync the URL, and tell the layout to open its tab.
  const selectItem = useCallback((rawId: string | undefined) => {
    setSelectedId(rawId ?? null);
    const oid = rawId ? idFromDocId(rawId) : undefined;
    const basePath = '/' + encodeURIComponent(project?.name ?? '');
    const newPath = oid ? `${basePath}/${oid}` : basePath;
    if (location.pathname !== newPath) {
      navigate(newPath, { replace: true });
    }
    onSelectItem(rawId);
  }, [project?.name, navigate, onSelectItem]);

  const handleDocumentCreated = useCallback((docOid: string) => {
    expandToDocument(docOid);
    selectItem(idDocId(docOid));
  }, [expandToDocument, selectItem]);

  // Sync the selection from the URL on mount and when the project changes.
  useEffect(() => {
    const rawId = (docId && project?.documentIds.has(docId)) ? idDocId(docId) : undefined;
    setSelectedId(rawId ?? null);
    onSelectItem(rawId);
  }, [project?.name]);

  // Expand branches to the initially selected document (e.g. from URL)
  const hasExpandedInitial = useRef(false);
  useEffect(() => {
    if (hasExpandedInitial.current) return;
    const docOid = selectedId ? idFromDocId(selectedId) : undefined;
    if (docOid && project?.documentIds.has(docOid)) {
      expandToDocument(docOid);
      hasExpandedInitial.current = true;
    }
  }, [project, selectedId, expandToDocument]);

  console.log(toolkits)
  console.log(project)

  return (
    <SimpleTreeView
      expandedItems={expandedItems}
      onExpandedItemsChange={(e, itemIds) => {
        const wasExpanded = expandedItems.includes(CENTRAL_REPO_FOLDER_ID);
        const willBeExpanded = itemIds.includes(CENTRAL_REPO_FOLDER_ID);
        if (wasExpanded !== willBeExpanded) {
          const target = (e as React.SyntheticEvent | null)?.target as HTMLElement | null;
          const fromChevron = !!target?.closest('[class*="iconContainer"]');
          if (!fromChevron) {
            const corrected = wasExpanded
              ? [...new Set([...itemIds, CENTRAL_REPO_FOLDER_ID])]
              : itemIds.filter(id => id !== CENTRAL_REPO_FOLDER_ID);
            setExpandedItems(corrected);
            return;
          }
        }
        setExpandedItems(itemIds);
      }}
      selectedItems={selectedId}
      onSelectedItemsChange={(_e, itemId) => selectItem(itemId ?? undefined)}
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
        <Stack direction={'row'} spacing={0} alignItems={'center'} justifyContent={'start'} sx={{ marginLeft: 5, width: 'fit-content' }}>
          <Folder sx={{ mr: 1 }} />
          <Tooltip title='Files directory where the project is located'>
            <Typography color={project.configDocument?.data.desc.filesDirectory ? 'text.primary' : 'text.secondary'}>
              {project.configDocument?.data.desc.filesDirectory || 'No directory'}
            </Typography>
          </Tooltip>
          {project.configDocument?.data.desc.filesDirectory && (
            <ButtonTooltip
              title='Copy path'
              onClick={() => navigator.clipboard.writeText(project.configDocument?.data.desc.filesDirectory ?? '')}
              sx={{ ml: 0.5 }}
            >
              <ContentCopy sx={{ fontSize: 16 }} />
            </ButtonTooltip>
          )}
        </Stack>
        <DocumentSplitGroup
          docs={project?.documents}
          project={project}
          depth={viewSettings.maxDepth}
          onDocumentDeleted={() => selectItem(undefined)}
        />
      </TreeItem>
      <RepoTreeWhole />
    </SimpleTreeView>
  );
};
