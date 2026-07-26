import { Add } from '@mui/icons-material';
import { Box } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { EditAsJsonButton } from './EditAsJsonButton';

// Actions shown on a branch (object-valued) row: add a child, or edit the whole
// subtree as JSON. A child becomes a nested object via its own type chip.
export const DetailsViewItemBranchActions = ({
  itemValue,
  setItemValue,
}: {
  itemValue: any,
  setItemValue: (newVal: any) => void,
}) => {
  const addSubItem = (initialValue: any) => {
    let name = '';
    for (let i = 1; i < 1e5; i++) {
      const key = 'newItem_' + i;
      if (!(key in itemValue)) {
        name = key;
        break;
      }
    }
    setItemValue({ ...itemValue, [name]: initialValue });
  };

  return (
    <>
      <ButtonTooltip
        title={'Add item'}
        onClick={() => addSubItem('')}
      >
        <Add />
      </ButtonTooltip>
      <Box className="field-json">
        <EditAsJsonButton
          data={itemValue}
          setData={setItemValue}
        />
      </Box>
    </>
  );
};
