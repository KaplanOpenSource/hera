import { Add, CreateNewFolder } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { EditAsJsonButton } from './EditAsJsonButton';

// Actions shown on a branch (object-valued) row: add a child, add a sub
// structure, or edit the whole subtree as JSON.
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
      <ButtonTooltip
        title={'Add sub structure'}
        onClick={() => addSubItem({})}
      >
        <CreateNewFolder />
      </ButtonTooltip>
      <EditAsJsonButton
        data={itemValue}
        setData={setItemValue}
      />
    </>
  );
};
