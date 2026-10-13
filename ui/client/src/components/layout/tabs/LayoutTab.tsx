import { IJsonTabNode } from 'flexlayout-react';
import { ReactNode } from 'react';
import { ProjectObj } from '../../../objects/ProjectObj';
import { TabKind } from './TabKind';

// What a tab needs to draw its panel. The config is the one stored on the tab node.
export interface TabRenderContext {
  config: any;
  project: ProjectObj;
  onSelectItem: (rawId: string | undefined) => void;
}

// One kind of dock tab: how its id is built, what its node looks like, which
// document it shows, and how its panel is drawn. Each kind is a single instance;
// it holds no state of its own, so it must not use hooks.
export abstract class LayoutTab {
  abstract readonly kind: TabKind;

  // Tab ids are "<kind>:<key>", where the key is a docid or a workflow name.
  get prefix(): string {
    return `${this.kind}:`;
  }

  id(key: string): string {
    return `${this.prefix}${key}`;
  }

  owns(tabId: string): boolean {
    return tabId.startsWith(this.prefix);
  }

  // The document this tab shows, or undefined when it shows none.
  docIdOf(config: any): string | undefined {
    return config?.docid;
  }

  // The parts every tab node shares; subclasses add the name and config.
  protected baseNode(key: string): IJsonTabNode {
    return { type: 'tab', id: this.id(key), component: this.kind };
  }

  abstract render(context: TabRenderContext): ReactNode;
}
