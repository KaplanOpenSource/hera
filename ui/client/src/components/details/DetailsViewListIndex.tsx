import { Typography } from '@mui/material';
import { DraggableAttributes } from '@dnd-kit/core';
import { SyntheticListenerMap } from '@dnd-kit/core/dist/hooks/utilities';

// A list element's index. It is the drag handle for reordering the list.
export const DetailsViewListIndex = ({
  index,
  attributes,
  listeners,
}: {
  index: number,
  attributes: DraggableAttributes,
  listeners?: SyntheticListenerMap,
}) => {
  return (
    <Typography
      data-testid={`list-index-${index}`}
      onClick={e => e.stopPropagation()}
      sx={{
        fontFamily: 'monospace',
        color: 'text.secondary',
        whiteSpace: 'nowrap',
        flexShrink: 0,
        cursor: 'grab',
        touchAction: 'none',
      }}
      {...attributes}
      {...listeners}
    >
      {`[${index}]`}
    </Typography>
  );
};
