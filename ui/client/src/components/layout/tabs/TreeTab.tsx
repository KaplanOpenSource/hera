import { Paper } from '@mui/material';
import { IJsonTabNode } from 'flexlayout-react';
import { ReactNode } from 'react';
import { ProjectTreeView } from '../../project/ProjectTreeView';
import { LayoutTab, TabRenderContext } from './LayoutTab';
import { TabKind } from './TabKind';

// The single Workspace Explorer tab. It is static: it cannot be closed, renamed,
// or dragged out of place, and it has no key, so its id is just the kind.
export class TreeTab extends LayoutTab {
  readonly kind = TabKind.Tree;

  // The only tab with no key, so its id is the kind on its own.
  readonly tabId = TabKind.Tree;

  get prefix(): string {
    return this.tabId;
  }

  node(): IJsonTabNode {
    return {
      type: 'tab',
      id: this.tabId,
      component: this.kind,
      name: 'Workspace Explorer',
      enableClose: false,
      enableDrag: false,
      enableRename: false,
    };
  }

  render({ project, onSelectItem }: TabRenderContext): ReactNode {
    return (
      <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
        <ProjectTreeView project={project} onSelectItem={onSelectItem} />
      </Paper>
    );
  }
}

export const treeTab = new TreeTab();
