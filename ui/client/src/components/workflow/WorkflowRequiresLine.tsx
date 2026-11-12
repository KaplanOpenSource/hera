import { Typography } from '@mui/material';

// The quiet line naming the nodes a node waits for. Nothing when it waits for none.
export const WorkflowRequiresLine = ({
  requires,
}: {
  requires: string[],
}) => {
  return (
    <>
      {requires.length > 0 && (
        <Typography noWrap variant="caption" sx={{ color: 'text.disabled', fontStyle: 'italic' }}>
          {`depends on: ${requires.join(', ')}`}
        </Typography>
      )}
    </>
  );
};
