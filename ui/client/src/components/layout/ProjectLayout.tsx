import { Action, Actions, ITabRenderValues, Layout, TabNode } from 'flexlayout-react';
import { useCallback, useEffect, useState } from 'react';
import { ProjectObj } from '../../objects/ProjectObj';
import { classifyItemId, idFromDocId, ItemKind, normalizeSplitId } from '../../shared/idDocId';
import { classifyTab } from '../../shared/tabKind';
import { TAB_KIND_STYLES } from '../../shared/tabKindConfig';
import { useFlexlayoutTheme } from '../../theme';
import { hasPreview } from '../details/PreviewPanel';
import { DETAILS_TAB_PREFIX, LayoutModel } from './LayoutModel';
import { LayoutPanel } from './LayoutPanel';

// Fallback "item" shown when the tree selection isn't a document/repo/split.
const CONFIG_ITEM_ID = 'config';

export const ProjectLayout = ({
  project,
  treeCollapsed,
  resetSignal,
}: {
  project: ProjectObj,
  treeCollapsed: boolean,
  resetSignal: number,
}) => {
  useFlexlayoutTheme();
  const [activeShowItemId, setActiveShowItemId] = useState<string | undefined>(undefined);

  const [layout, setLayout] = useState(() => LayoutModel.create(!treeCollapsed));

  const activeDocId = activeShowItemId ? idFromDocId(activeShowItemId) : undefined;
  const activeDoc = activeDocId
    ? project.allDocuments.find(d => d.docid === activeDocId)
    : undefined;
  const previewAvailable = activeDoc ? hasPreview(activeDoc.data) : false;

  // When the reset button in the header is pressed, rebuild the layout in its
  // original left/right arrangement while keeping the currently open details tabs.
  useEffect(() => {
    if (resetSignal === 0) return;
    setLayout(layout.resetKeepingDetails(!treeCollapsed, activeShowItemId));
  }, [resetSignal]);

  useEffect(() => {
    layout.setTreeVisible(!treeCollapsed);
  }, [treeCollapsed, layout]);

  // The tree owns the selection; when it changes, open or focus that item's details tab.
  const handleSelectItem = useCallback((rawShowItemId: string | undefined) => {
    if (!rawShowItemId) return;

    const kind = classifyItemId(rawShowItemId);
    let showItemId = rawShowItemId;
    if (kind === ItemKind.Config) showItemId = CONFIG_ITEM_ID;
    else if (kind === ItemKind.Split) showItemId = normalizeSplitId(rawShowItemId, project.documents);

    layout.openOrFocusDetailsTab(showItemId, project);
    setActiveShowItemId(showItemId);
  }, [layout, project]);

  useEffect(() => {
    if (previewAvailable && activeDocId) {
      layout.setPreview(activeDocId, activeDoc!.name);
    } else {
      layout.setPreview();
    }
  }, [activeShowItemId, previewAvailable, activeDocId, layout]);

  // Keep open tabs in sync with the project: close tabs whose document was
  // deleted, then rename the survivors (e.g. a new notebook's tab once its
  // document loads).
  useEffect(() => {
    layout.closeMissingDocuments(project);
    layout.syncTabNames(project);
  }, [project, layout]);

  const handleAction = useCallback((action: Action) => {
    if (action.type === Actions.SELECT_TAB) {
      const tabId = action.data.tabNode as string;
      if (tabId?.startsWith(DETAILS_TAB_PREFIX)) {
        const node = layout.getTab(tabId);
        if (node) {
          setActiveShowItemId(node.getConfig()?.showItemId);
        }
      }
    }
    return action;
  }, [layout]);

  const onRenderTab = useCallback((node: TabNode, renderValues: ITabRenderValues) => {
    const showItemId = node.getConfig()?.showItemId as string | undefined;
    if (!showItemId) return;
    const kind = classifyTab(showItemId, project);
    if (!kind) return;
    const { icon: Icon, color } = TAB_KIND_STYLES[kind];
    renderValues.leading = <Icon sx={{ fontSize: 16, color }} />;
  }, [project.allDocuments]);

  const factory = (node: TabNode) => (
    <LayoutPanel
      component={node.getComponent()}
      config={node.getConfig()}
      project={project}
      onSelectItem={handleSelectItem}
    />
  );

  return (
    <Layout
      model={layout.model}
      factory={factory}
      onRenderTab={onRenderTab}
      onAction={handleAction}
    />
  );
};
