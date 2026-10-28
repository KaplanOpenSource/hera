import { IJsonTabNode } from 'flexlayout-react';
import { ReactNode } from 'react';
import { WorkflowOutputPanel } from '../../workflow/log/WorkflowOutputPanel';
import { LayoutTab, TabRenderContext } from './LayoutTab';
import { TabKind } from './TabKind';

// A workflow run's output. Keyed by workflow name, not by document, so it shows no
// document of its own.
export class OutputTab extends LayoutTab {
  readonly kind = TabKind.Output;

  node(workflowName: string): IJsonTabNode {
    return {
      ...this.baseNode(workflowName),
      name: `Output: ${workflowName}`,
      config: { workflowName },
    };
  }

  render({ config }: TabRenderContext): ReactNode {
    return <WorkflowOutputPanel workflowName={config?.workflowName} />;
  }
}

export const outputTab = new OutputTab();
