import { canvasTab } from './CanvasTab';
import { detailsTab } from './DetailsTab';
import { LayoutTab } from './LayoutTab';
import { outputTab } from './OutputTab';
import { previewTab } from './PreviewTab';
import { treeTab } from './TreeTab';

// Every kind of tab, for looking one up by a tab node's id or component.
export const ALL_TABS: LayoutTab[] = [treeTab, detailsTab, previewTab, canvasTab, outputTab];

// The tab kind that owns a tab id, or undefined for an id of no known kind.
export const tabForId = (tabId: string): LayoutTab | undefined => {
  return ALL_TABS.find(tab => tab.owns(tabId));
};

// The tab kind a flexlayout component id names.
export const tabForComponent = (component: string | undefined): LayoutTab | undefined => {
  return ALL_TABS.find(tab => tab.kind === component);
};
