import { Paper } from '@mui/material';
import { ReactNode } from 'react';
import { ProjectObj } from '../../objects/ProjectObj';
import { DetailsViewPanel } from '../details/DetailsViewPanel';
import { PreviewPanel } from '../details/PreviewPanel';
import { ProjectTreeView } from '../project/ProjectTreeView';

// The component identifiers flexlayout stores on each tab and passes to the factory.
export enum LayoutComponent {
  Tree = 'tree',
  Details = 'details',
  Preview = 'preview',
}

// Renders the content for a single dock node, chosen by its flexlayout component id.
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
}) => {
  let content: ReactNode = null;
  switch (component) {
    case LayoutComponent.Tree:
      content = (
        <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
          <ProjectTreeView project={project} onSelectItem={onSelectItem} />
        </Paper>
      );
      break;
    case LayoutComponent.Details:
      content = (
        <Paper sx={{ height: '100%', overflow: 'hidden' }}>
          <DetailsViewPanel project={project} showItemId={config?.showItemId} />
        </Paper>
      );
      break;
    case LayoutComponent.Preview:
      content = <PreviewPanel docid={config?.docid} />;
      break;
  }
  return content;
};
