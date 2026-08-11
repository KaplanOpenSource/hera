import { Delete } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';

// Removes a field from its parent object by clearing its key. Renders nothing
// when the row can't be deleted (no `setItemKey`). Styled as a white trash icon
// on a red circle, meant to overlay the right side of the field name, so it
// never takes layout space.
export const DeleteFieldButton = ({
  itemKey,
  setItemKey,
}: {
  itemKey: string,
  setItemKey?: (newKey: string | undefined) => void,
}) => {
  return (
    setItemKey && (
      <ButtonTooltip
        title={'Delete ' + itemKey}
        onClick={() => setItemKey(undefined)}
        sx={{
          color: 'common.white',
          bgcolor: 'error.main',
          '&:hover': { bgcolor: 'error.dark' },
        }}
      >
        <Delete fontSize='small' />
      </ButtonTooltip>
    )
  );
};
