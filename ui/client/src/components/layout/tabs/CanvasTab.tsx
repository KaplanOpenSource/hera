import { IJsonTabNode } from 'flexlayout-react';
import { ReactNode } from 'react';
import { WorkflowCanvasPanel } from '../../workflow/WorkflowCanvasPanel';
import { LayoutTab, TabRenderContext } from './LayoutTab';
import { TabKind } from './TabKind';

// A workflow document's canvas. It opens and closes with the document's details tab.
export class CanvasTab extends LayoutTab {
  readonly kind = TabKind.Canvas;

  node(docid: string, docName: string): IJsonTabNode {
    return {
      ...this.baseNode(docid),
      name: `Canvas: ${docName}`,
      config: { docid },
    };
  }

  render({ config, project }: TabRenderContext): ReactNode {
    return <WorkflowCanvasPanel project={project} docid={config?.docid} />;
  }
}

export const canvasTab = new CanvasTab();
