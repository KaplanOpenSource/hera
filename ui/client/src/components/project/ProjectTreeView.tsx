import { Box, Typography } from '@mui/material';
import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idDocId, idFromDocId } from '../../shared/idDocId';
import { projectPath } from '../../shared/projectPath';
import { useProjectStore } from '../../stores/useProjectStore';
import { documentMatchesQuery, documentSearchText, parseSearchQuery, unknownSearchFields } from '../../utils/documentSearch';
import { collectBranchKeys, SplitTree } from '../../utils/splitTree';
import { DocumentTreeWhole } from './DocumentTreeWhole';
import { RepoTreeWhole } from './RepoTreeWhole';
import { TreeSearchBar } from './TreeSearchBar';
import { useViewSettingsStore } from '../../stores/useViewSettingsStore';

export const ProjectTreeView = ({
  project,
  onSelectItem,
}: {
  project: ProjectObj;
  onSelectItem: (rawId: string | undefined) => void,
}) => {
  const { docId } = useParams<{ docId: string }>();
  const navigate = useNavigate();
  const { viewSettings } = useViewSettingsStore();
  const [selectedIds, setSelectedIds] = useState<string[]>(docId ? [idDocId(docId)] : []);
  const [expandedItems, setExpandedItems] = useState<string[]>(['project-documents', 'no-toolkit']);
  // Repositories are a separate tree with their own selection/expansion state.
  const [repoSelectedIds, setRepoSelectedIds] = useState<string[]>([]);
  const [repoExpandedItems, setRepoExpandedItems] = useState<string[]>([]);
  const [search, setSearch] = useState('');
  const splitTreeRef = useRef<SplitTree | null>(null);

  // Search index: one lowercased value-blob per document, rebuilt only when the project changes.
  const searchIndex = useMemo(
    () => project.documents.map(doc => ({ doc, text: documentSearchText(doc.data) })),
    [project],
  );

  const searchTerms = useMemo(() => parseSearchQuery(search), [search]);
  const isSearching = searchTerms.length > 0;
  const unknownFields = useMemo(() => unknownSearchFields(searchTerms), [searchTerms]);
  const filteredDocs = useMemo(
    () => isSearching
      ? searchIndex.filter(e => documentMatchesQuery(e.doc.data, e.text, searchTerms)).map(e => e.doc)
      : searchIndex.map(e => e.doc),
    [searchIndex, searchTerms, isSearching],
  );

  // The tree as currently shown, so branch clicks only take the visible documents.
  const displayTree = useMemo(
    () => new SplitTree(filteredDocs, viewSettings.maxDepth, viewSettings),
    [filteredDocs, viewSettings],
  );

  // While searching, expand every matching branch so results aren't hidden in collapsed groups.
  const searchExpandedKeys = useMemo(() => {
    if (!isSearching) return [];
    return ['project-documents', ...collectBranchKeys(displayTree.nodes)];
  }, [isSearching, displayTree]);

  const effectiveExpandedItems = isSearching
    ? [...new Set([...expandedItems, ...searchExpandedKeys])]
    : expandedItems;

  const searchWarning = unknownFields.length
    ? `Unknown field${unknownFields.length > 1 ? 's' : ''}: ${unknownFields.join(', ')}`
    : undefined;

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
    if (ancestors) {
      setExpandedItems(prev => [...new Set([...prev, 'project-documents', ...ancestors])]);
    }
  }, [getSplitTree]);

  // Sync the URL to the given item (the active one), or to the project root if none.
  const navigateToItem = useCallback((rawId: string | undefined) => {
    const oid = rawId ? idFromDocId(rawId) : undefined;
    const newPath = projectPath(project?.name ?? '', oid);
    if (location.pathname !== newPath) {
      navigate(newPath, { replace: true });
    }
  }, [project?.name, navigate]);

  // Opens the item in a tab. Only a double click (or Enter) does this.
  const openItem = useCallback((rawId: string | undefined) => {
    if (!rawId) return;
    navigateToItem(rawId);
    onSelectItem(rawId);
  }, [navigateToItem, onSelectItem]);

  // Only one of the two trees shows a selection at a time.
  const handleSelectionChange = useCallback((ids: string[]) => {
    setSelectedIds(ids);
    setRepoSelectedIds([]);   // documents and repos share a single active highlight
  }, []);

  // Same as handleSelectionChange, but for the separate repositories tree.
  const handleRepoSelectionChange = useCallback((event: React.SyntheticEvent | null, ids: string[]) => {
    const target = (event?.target as HTMLElement | null);
    if (target?.closest('[class*="iconContainer"]')) return;

    setRepoSelectedIds(ids);
    setSelectedIds([]);   // clear the documents selection so only one item is active
  }, []);

  // The central repo folder only toggles from its chevron, not by clicking the row.
  const handleRepoExpandedChange = useCallback((e: React.SyntheticEvent | null, itemIds: string[]) => {
    const wasExpanded = repoExpandedItems.includes(CENTRAL_REPO_FOLDER_ID);
    const willBeExpanded = itemIds.includes(CENTRAL_REPO_FOLDER_ID);
    if (wasExpanded !== willBeExpanded) {
      const target = (e?.target as HTMLElement | null);
      const fromChevron = !!target?.closest('[class*="iconContainer"]');
      if (!fromChevron) {
        const corrected = wasExpanded
          ? [...new Set([...itemIds, CENTRAL_REPO_FOLDER_ID])]
          : itemIds.filter(id => id !== CENTRAL_REPO_FOLDER_ID);
        setRepoExpandedItems(corrected);
        return;
      }
    }
    setRepoExpandedItems(itemIds);
  }, [repoExpandedItems]);

  // Forgets the selection and points the URL at the project root.
  const clearSelection = useCallback(() => {
    setSelectedIds([]);
    navigateToItem(undefined);
  }, [navigateToItem]);

  // Make the given document the selection (highlight, URL, open tab), or clear it when none.
  const selectDocument = useCallback((docOid?: string) => {
    if (!docOid) {
      clearSelection();
      return;
    }
    expandToDocument(docOid);
    const id = idDocId(docOid);
    setSelectedIds([id]);
    navigateToItem(id);
    onSelectItem(id);
  }, [clearSelection, expandToDocument, navigateToItem, onSelectItem]);

  // A click on the empty space of the panel clears both selections.
  const handleBackgroundClick = useCallback((event: React.MouseEvent) => {
    const target = event.target as HTMLElement | null;
    if (target?.closest('.MuiTreeItem-content, input, button')) return;
    setSelectedIds([]);
    setRepoSelectedIds([]);
  }, []);

  // Sync the selection from the URL on mount and when the project changes.
  useEffect(() => {
    const rawId = (docId && project?.documentIds.has(docId)) ? idDocId(docId) : undefined;
    setSelectedIds(rawId ? [rawId] : []);
    onSelectItem(rawId);
  }, [project?.name]);

  // The URL also changes when a tab is closed, so the highlight follows it.
  useEffect(() => {
    setSelectedIds(docId ? [idDocId(docId)] : []);
  }, [docId]);

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

  return (
    <Box sx={{ minHeight: '100%' }} onClick={handleBackgroundClick}>
    <TreeSearchBar value={search} onChange={setSearch} warning={searchWarning} />
    <Typography
      variant="overline"
      sx={{ display: 'block', mt: 1, mb: 1, color: 'text.secondary', fontWeight: 600, letterSpacing: 1 }}
    >
      Workspace Explorer
    </Typography>
    <DocumentTreeWhole
      project={project}
      docs={filteredDocs}
      tree={displayTree}
      depth={viewSettings.maxDepth}
      selectedIds={selectedIds}
      onSelectedIdsChange={handleSelectionChange}
      expandedItems={effectiveExpandedItems}
      onExpandedItemsChange={setExpandedItems}
      onOpenItem={openItem}
      onSelectDocument={selectDocument}
      onCleared={clearSelection}
    />
    <RepoTreeWhole
      selectedIds={repoSelectedIds}
      onSelectedItemsChange={handleRepoSelectionChange}
      expandedItems={repoExpandedItems}
      onExpandedItemsChange={handleRepoExpandedChange}
      onOpenItem={openItem}
    />
    </Box>
  );
};
