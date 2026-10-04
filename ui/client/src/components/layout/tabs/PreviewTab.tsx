import { IJsonTabNode } from 'flexlayout-react';
import { ReactNode } from 'react';
import { PreviewPanel } from '../../details/PreviewPanel';
import { LayoutTab, TabRenderContext } from './LayoutTab';
import { TabKind } from './TabKind';

// A document's preview pane. Only one is open at a time.
export class PreviewTab extends LayoutTab {
  readonly kind = TabKind.Preview;

  node(docid: string, docName: string): IJsonTabNode {
    return {
      ...this.baseNode(docid),
      name: `Preview: ${docName}`,
      config: { docid },
    };
  }

  render({ config }: TabRenderContext): ReactNode {
    return <PreviewPanel docid={config?.docid} />;
  }
}

export const previewTab = new PreviewTab();
