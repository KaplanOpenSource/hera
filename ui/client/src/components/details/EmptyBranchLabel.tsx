import { Typography } from '@mui/material';

// Placeholder shown under a branch (object-valued) row that has no children.
// `level` is the row's depth, used to line the text up with where children
// would appear.
export const EmptyBranchLabel = ({
  level,
}: {
  level: number,
}) => {
  return (
    <Typography
      variant="body2"
      sx={{
        fontStyle: 'italic',
        color: 'text.secondary',
        ml: `${30 + (12 * level)}px`,
      }}
    >
      (empty)
    </Typography>
  );
};
