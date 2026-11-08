// A run status as the server sends it, from start (busy) or poll (the rest).
export enum WorkflowServerStatus {
  Running = 'running',
  Done = 'done',
  Error = 'error',
  Stopped = 'stopped',
  Busy = 'busy',
  NotFound = 'not_found',
}
