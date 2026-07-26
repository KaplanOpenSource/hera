import { Delete } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';

// Removes a field from its parent object by clearing its key. Renders nothing
// when the row can't be deleted (no `setItemKey`).
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
      >
        <Delete fontSize='small'/>
      </ButtonTooltip>
    )
  );
};
