import { ContentCopy, Folder } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idDocId, idFromDocId } from '../../shared/idDocId';
import { useProjectStore } from '../../stores/useProjectStore';
import { documentSearchText } from '../../utils/documentSearch';
import { collectBranchKeys, SplitTree } from '../../utils/splitTree';
import { DocumentSplitGroup } from './DocumentSplitGroup';
import { ProjectActionsButton } from './ProjectActionsButton';
import { RepoTreeWhole } from './RepoTreeWhole';
import { TreeSearchBar } from './TreeSearchBar';
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
  const [selectedIds, setSelectedIds] = useState<string[]>(docId ? [idDocId(docId)] : []);
  const [expandedItems, setExpandedItems] = useState<string[]>(['project-documents', 'no-toolkit', '*repos*']);
  const [search, setSearch] = useState('');
  const splitTreeRef = useRef<SplitTree | null>(null);

  // Search index: one lowercased value-blob per document, rebuilt only when the project changes.
  const searchIndex = useMemo(
    () => project.documents.map(doc => ({ doc, text: documentSearchText(doc.data) })),
    [project],
  );

  const query = search.trim().toLowerCase();
  const filteredDocs = useMemo(
    () => query ? searchIndex.filter(e => e.text.includes(query)).map(e => e.doc) : searchIndex.map(e => e.doc),
    [searchIndex, query],
  );

  // While searching, expand every matching branch so results aren't hidden in collapsed groups.
  const searchExpandedKeys = useMemo(() => {
    if (!query) return [];
    const tree = new SplitTree(filteredDocs, viewSettings.maxDepth, viewSettings);
    return ['project-documents', ...collectBranchKeys(tree.nodes)];
  }, [query, filteredDocs, viewSettings]);

  const effectiveExpandedItems = query
    ? [...new Set([...expandedItems, ...searchExpandedKeys])]
    : expandedItems;

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

  // Sync the URL to the given item (the active one), or to the project root if none.
  const navigateToItem = useCallback((rawId: string | undefined) => {
    const oid = rawId ? idFromDocId(rawId) : undefined;
    const basePath = '/' + encodeURIComponent(project?.name ?? '');
    const newPath = oid ? `${basePath}/${oid}` : basePath;
    if (location.pathname !== newPath) {
      navigate(newPath, { replace: true });
    }
  }, [project?.name, navigate]);

  // Multi-select highlights rows, but only the most recently added item opens/focuses a tab.
  const handleSelectionChange = useCallback((event: React.SyntheticEvent | null, ids: string[]) => {
    // Clicking the expand/collapse chevron shouldn't change the selection.
    const target = (event?.target as HTMLElement | null);
    if (target?.closest('[class*="iconContainer"]')) return;

    const added = ids.filter(id => !selectedIds.includes(id));
    setSelectedIds(ids);
    const clicked = added[added.length - 1];
    if (clicked) {
      navigateToItem(clicked);
      onSelectItem(clicked);
    }
  }, [selectedIds, navigateToItem, onSelectItem]);

  // Make the given document the selection (highlight, URL, open tab), or clear it when none.
  const selectDocument = useCallback((docOid?: string) => {
    if (!docOid) {
      setSelectedIds([]);
      navigateToItem(undefined);
      return;
    }
    expandToDocument(docOid);
    const id = idDocId(docOid);
    setSelectedIds([id]);
    navigateToItem(id);
    onSelectItem(id);
  }, [expandToDocument, navigateToItem, onSelectItem]);

  // Sync the selection from the URL on mount and when the project changes.
  useEffect(() => {
    const rawId = (docId && project?.documentIds.has(docId)) ? idDocId(docId) : undefined;
    setSelectedIds(rawId ? [rawId] : []);
    onSelectItem(rawId);
  }, [project?.name]);

  // Expand branches to the initially selected document (e.g. from URL)
  const hasExpandedInitial = useRef(false);
  useEffect(() => {
    if (hasExpandedInitial.current) return;
    const docOid = selectedIds[0] ? idFromDocId(selectedIds[0]) : undefined;
    if (docOid && project?.documentIds.has(docOid)) {
      expandToDocument(docOid);
      hasExpandedInitial.current = true;
    }
  }, [project, selectedIds, expandToDocument]);

  console.log(toolkits)
  console.log(project)

  return (
    <>
    <TreeSearchBar value={search} onChange={setSearch} />
    <SimpleTreeView
      expandedItems={effectiveExpandedItems}
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
      selectedItems={selectedIds}
      onSelectedItemsChange={(e, itemIds) => handleSelectionChange(e, itemIds)}
      expansionTrigger={'content'}
      multiSelect
      sx={{
        // Highlighted (selected) rows: soft blue.
        '& .MuiTreeItem-content.Mui-selected': {
          backgroundColor: 'rgba(25, 118, 210, 0.24)',
        },
        '& .MuiTreeItem-content.Mui-selected:hover': {
          backgroundColor: 'rgba(25, 118, 210, 0.34)',
        },
        // Active (focused) row stands out with a stronger blue.
        '& .MuiTreeItem-content.Mui-selected.Mui-focused': {
          backgroundColor: 'rgba(25, 118, 210, 0.44)',
        },
        '& .MuiTreeItem-content.Mui-selected.Mui-focused:hover': {
          backgroundColor: 'rgba(25, 118, 210, 0.54)',
        },
      }}
    >
      <TreeItem key={`project-documents`} itemId={`project-documents`}
        label={(
          <Stack direction='row' justifyContent="start" alignItems='center'>
            <Typography marginRight={1}>
              Project {project.name}
            </Typography>
            <ProjectActionsButton
              selectedIds={selectedIds}
              onSelectDocument={selectDocument}
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
          docs={filteredDocs}
          project={project}
          depth={viewSettings.maxDepth}
          onDocumentDeleted={() => {
            setSelectedIds([]);
            navigateToItem(undefined);
          }}
        />
      </TreeItem>
      <RepoTreeWhole />
    </SimpleTreeView>
    </>
  );
};
