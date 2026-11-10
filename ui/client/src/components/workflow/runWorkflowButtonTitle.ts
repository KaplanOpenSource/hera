// The run button's tooltip text; also its aria-label, so it is a string, not a component.
export const runWorkflowButtonTitle = ({
  canStop,
  isStarting,
  disabled,
  disabledReason,
}: {
  canStop: boolean,
  isStarting: boolean,
  disabled?: boolean,
  disabledReason?: string,
}): string => {
  let title = 'Run the workflow as shown (right click for options)';
  if (canStop) {
    title = 'Stop the workflow';
  } else if (isStarting) {
    title = 'Workflow is starting…';
  } else if (disabled && disabledReason) {
    title = disabledReason;
  }
  return title;
};
