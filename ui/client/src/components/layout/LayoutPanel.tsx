import { ReactNode } from 'react';
import { ProjectObj } from '../../objects/ProjectObj';
import { tabForComponent } from './tabs/allTabs';

// Renders the content for a single dock node, by asking the tab kind its
// flexlayout component id names.
export const LayoutPanel = ({
  component,
  config,
  project,
  onSelectItem,
}: {
  component: string | undefined,
  config: any,
  project: ProjectObj,
  onSelectItem: (rawId: string | undefined) => void,
}): ReactNode => {
  const tab = tabForComponent(component);
  if (!tab) {
    return null;
  }
  return tab.render({ config, project, onSelectItem });
};
