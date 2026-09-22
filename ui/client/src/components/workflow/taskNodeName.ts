// Luigi names each generated task <node>_<index>, so a task name from the run
// output has to be mapped back to the node it came from. The longest matching
// node wins, so "Run" does not swallow a task belonging to "Run_extra".
export const taskBelongsToNode = (taskName: string, nodeName: string): boolean => {
  return taskName === nodeName || (taskName.startsWith(nodeName) && /^_\d+$/.test(taskName.slice(nodeName.length)));
};

export const nodeNameFromTask = (taskName: string, nodeNames: string[]): string | undefined => {
  const matches = nodeNames.filter((name) => { return taskBelongsToNode(taskName, name); });
  return matches.sort((a, b) => { return b.length - a.length; })[0];
};
