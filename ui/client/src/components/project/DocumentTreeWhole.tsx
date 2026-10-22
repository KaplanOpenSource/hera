import { ContentCopy, Folder } from '@mui/icons-material';
import { Stack, Tooltip, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { useRef, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId, isSplitId } from '../../shared/idDocId';
import { clickedRow } from '../../utils/clickedRow';
import { collectSubtreeKeys, findNodeByKey } from '../../utils/collectSubtreeKeys';
import { SplitTree } from '../../utils/splitTree';
import { DocumentSplitGroup } from './DocumentSplitGroup';
import { ProjectActionsButton } from './ProjectActionsButton';
import { MenuDocument, MenuPosition, TreeContextMenu } from './TreeContextMenu';
import { treeSelectionSx } from './treeSelectionSx';

// The documents section: the project tree and its right-click menu.
// Selection and expansion state are owned by the parent.
export const DocumentTreeWhole = ({
  project,
  docs,
  tree,
  depth,
  selectedIds,
  onSelectedIdsChange,
  expandedItems,
  onExpandedItemsChange,
  onOpenItem,
  onSelectDocument,
  onCleared,
}: {
  project: ProjectObj,
  docs: DocumentObj[],
  // The tree as shown, used to find what a branch row holds.
  tree: SplitTree,
  depth: number,
  selectedIds: string[],
  onSelectedIdsChange: (ids: string[]) => void,
  expandedItems: string[],
  onExpandedItemsChange: (ids: string[]) => void,
  onOpenItem: (rawId: string | undefined) => void,
  onSelectDocument: (docOid?: string) => void,
  // Called when the selected documents are gone, e.g. after a delete.
  onCleared: () => void,
}) => {
  const [menuPosition, setMenuPosition] = useState<MenuPosition | null>(null);
  const [menuDocs, setMenuDocs] = useState<MenuDocument[]>([]);
  // The row clicked last. A double click opens it.
  const lastClickedRef = useRef<string | undefined>(undefined);

  // Of the given tree ids, the documents, without the project config document.
  const documentsOf = (ids: string[]): MenuDocument[] => {
    const configDocId = project.configDocument?.docid;
    const oids = ids
      .map(idFromDocId)
      .filter((oid): oid is string => !!oid && oid !== configDocId);
    return oids.map(oid => {
      const doc = project.documents.find(d => d.docid === oid);
      return { oid, name: doc?.name ?? oid };
    });
  };

  // All the rows a branch row stands for, or just the row itself.
  const rowWithSubtree = (itemId: string) => {
    if (isSplitId(itemId)) {
      const node = findNodeByKey(tree.nodes, itemId);
      if (node) {
        return collectSubtreeKeys(node);
      }
    }
    return [itemId];
  };

  // A click only highlights rows. Clicking a branch takes everything under it.
  const handleSelectedItemsChange = (event: React.SyntheticEvent | null, ids: string[]) => {
    // Clicking the expand/collapse chevron shouldn't change the selection.
    const target = (event?.target as HTMLElement | null);
    if (target?.closest('[class*="iconContainer"]')) return;

    const clicked = clickedRow(ids, selectedIds);
    lastClickedRef.current = clicked;
    let nextIds = ids;
    if (clicked) {
      nextIds = [...new Set([...ids, ...rowWithSubtree(clicked)])];
    }
    onSelectedIdsChange(nextIds);
  };

  // Right click acts on the selection, after making sure the clicked row is in it.
  const handleRowContextMenu = (event: React.MouseEvent, itemId: string) => {
    event.preventDefault();
    event.stopPropagation();
    let actedOn = selectedIds;
    if (!selectedIds.includes(itemId)) {
      actedOn = rowWithSubtree(itemId);
      onSelectedIdsChange(actedOn);
      lastClickedRef.current = itemId;
    }
    setMenuDocs(documentsOf(actedOn));
    setMenuPosition({ x: event.clientX, y: event.clientY });
  };

  const filesDirectory = project.configDocument?.data.desc.filesDirectory;

  return (
    <>
      <SimpleTreeView
        expandedItems={expandedItems}
        onExpandedItemsChange={(_e, itemIds) => onExpandedItemsChange(itemIds)}
        selectedItems={selectedIds}
        onSelectedItemsChange={handleSelectedItemsChange}
        onDoubleClick={() => onOpenItem(lastClickedRef.current)}
        onKeyDown={e => e.key === 'Enter' && onOpenItem(lastClickedRef.current)}
        expansionTrigger={'content'}
        multiSelect
        sx={treeSelectionSx}
      >
        <TreeItem key={`project-documents`} itemId={`project-documents`}
          label={(
            <Stack direction='row' justifyContent="start" alignItems='center'>
              <Typography marginRight={1}>
                Project {project.name}
              </Typography>
              <ProjectActionsButton
                selectedIds={selectedIds}
                onSelectDocument={onSelectDocument}
              />
            </Stack>
          )}
        >
          <Stack direction={'row'} spacing={0} alignItems={'center'} justifyContent={'start'} sx={{ marginLeft: 5, width: 'fit-content' }}>
            <Folder sx={{ mr: 1 }} />
            <Tooltip title='Files directory where the project is located'>
              <Typography color={filesDirectory ? 'text.primary' : 'text.secondary'}>
                {filesDirectory || 'No directory'}
              </Typography>
            </Tooltip>
            {filesDirectory && (
              <ButtonTooltip
                title='Copy path'
                onClick={() => navigator.clipboard.writeText(filesDirectory ?? '')}
                sx={{ ml: 0.5 }}
              >
                <ContentCopy sx={{ fontSize: 16 }} />
              </ButtonTooltip>
            )}
          </Stack>
          <DocumentSplitGroup
            docs={docs}
            project={project}
            depth={depth}
            onDocumentDeleted={onCleared}
            onRowContextMenu={handleRowContextMenu}
          />
        </TreeItem>
      </SimpleTreeView>
      <TreeContextMenu
        position={menuPosition}
        docs={menuDocs}
        onClose={() => setMenuPosition(null)}
        onDeleted={onCleared}
      />
    </>
  );
};
