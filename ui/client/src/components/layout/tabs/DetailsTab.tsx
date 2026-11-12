import { Paper } from '@mui/material';
import { IJsonTabNode } from 'flexlayout-react';
import { ReactNode } from 'react';
import { ProjectObj } from '../../../objects/ProjectObj';
import { idFromDocId } from '../../../shared/idDocId';
import { tabKindClassName } from '../../../shared/tabKind';
import { detailsTabName, DetailsViewPanel } from '../../details/DetailsViewPanel';
import { LayoutTab, TabRenderContext } from './LayoutTab';
import { TabKind } from './TabKind';

// A tree item's details view. Keyed by the item id, which is not always a document.
export class DetailsTab extends LayoutTab {
  readonly kind = TabKind.Details;

  docIdOf(config: any): string | undefined {
    return idFromDocId(config?.showItemId ?? '');
  }

  // The tab's title. An item may be named before its document has loaded, so this
  // is recomputed as the project changes.
  tabName(showItemId: string, project: ProjectObj): string {
    return detailsTabName(showItemId, project);
  }

  node(showItemId: string, project: ProjectObj): IJsonTabNode {
    return {
      ...this.baseNode(showItemId),
      name: this.tabName(showItemId, project),
      className: tabKindClassName(showItemId, project),
      config: { showItemId },
    };
  }

  render({ config, project }: TabRenderContext): ReactNode {
    return (
      <Paper sx={{ height: '100%', overflow: 'hidden' }}>
        <DetailsViewPanel project={project} showItemId={config?.showItemId} />
      </Paper>
    );
  }
}

export const detailsTab = new DetailsTab();
