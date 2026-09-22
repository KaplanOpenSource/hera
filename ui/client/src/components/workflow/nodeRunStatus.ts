import { WorkflowChunk } from '../../io/runWorkflow';
import { EVENT_PREFIX } from './log/classifyLog';
import { nodeNameFromTask } from './taskNodeName';

// How a node of the workflow did in the current run.
export enum NodeRunStatus {
  Pending = 'pending',
  Running = 'running',
  Success = 'success',
  Failure = 'failure',
}

export type NodeRunStatusMap = { [nodeName: string]: NodeRunStatus };

// The Luigi events that change a task's state, as printed by the server runner.
const EVENT_TO_STATUS: { [event: string]: NodeRunStatus } = {
  START: NodeRunStatus.Running,
  SUCCESS: NodeRunStatus.Success,
  FAILURE: NodeRunStatus.Failure,
  BROKEN: NodeRunStatus.Failure,
};

// "[luigi-event] SUCCESS taskname" or "[luigi-event] FAILURE taskname: message".
const EVENT_LINE = new RegExp(`^${EVENT_PREFIX.replace(/[[\]]/g, '\\$&')} (\\w+) ([^\\s:]+)`);

// Reads the run's output and says how each node is doing: the last event for a
// task wins, so a task that started and then failed shows as failed. Nodes with
// no event yet stay pending. Once the run is over, a task that started but never
// reported an end is counted as failed - the run stopped on it.
export const nodeRunStatuses = ({
  chunks,
  nodeNames,
  runFinished,
}: {
  chunks: WorkflowChunk[],
  nodeNames: string[],
  runFinished: boolean,
}): NodeRunStatusMap => {
  const statuses: NodeRunStatusMap = {};
  for (const chunk of chunks) {
    for (const line of chunk.text.split('\n')) {
      const match = EVENT_LINE.exec(line);
      if (!match) {
        continue;
      }
      const status = EVENT_TO_STATUS[match[1]];
      const nodeName = nodeNameFromTask(match[2], nodeNames);
      if (status && nodeName) {
        statuses[nodeName] = status;
      }
    }
  }
  if (runFinished) {
    for (const nodeName of Object.keys(statuses)) {
      if (statuses[nodeName] === NodeRunStatus.Running) {
        statuses[nodeName] = NodeRunStatus.Failure;
      }
    }
  }
  return statuses;
};
