import { Add, CreateNewFolderOutlined, Delete } from '@mui/icons-material';
import { Stack } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DetailsViewItemSingle } from './DetailsViewItemSingle';

export const keyForDetailsViewItem = (itemKey: string, level: number = 0, index: number = 0) => {
  return `___lvl${level}_idx${index}_${itemKey}`;
}

export const DetailsViewItem = ({
  itemKey,
  itemValue,
  setItemValue,
  setItemKey = undefined,
  level = 0,
  index,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  setItemKey?: (newKey: string | undefined) => void | undefined,
  level: number,
  index: number,
}) => {
  const key = keyForDetailsViewItem(itemKey, level, index);

  const addSubItem = () => {
    let name = '';
    for (let i = 1; i < 1e5; i++) {
      const key = 'newItem_' + i;
      if (!(key in itemValue)) {
        name = key;
        break;
      }
    }
    setItemValue({ ...itemValue, [name]: '' })
  }

  return typeof itemValue === 'object'
    ? (
      <TreeItem
        key={key}
        itemId={key}
        label={<Stack direction={'row'} alignItems={'center'}>
          {itemKey}
          <ButtonTooltip
            title={'Add item'}
            onClick={addSubItem}
          >
            <Add />
          </ButtonTooltip>
          {/* <ButtonTooltip
            title={'Add substructure'}
            onClick={() => { }}
          >
            <CreateNewFolderOutlined />
          </ButtonTooltip> */}
        </Stack>}
      >
        {Object.entries(itemValue).map(([k, v], i) => {
          const isDir = level === 1 && itemKey === 'desc' && k === 'filesDirectory';
          const changeKey = (newKey: string | undefined) => {
            const item = { ...itemValue };
            delete item[k];
            if (newKey !== undefined) {
              item[newKey] = v;
            }
            setItemValue(item);
          };
          return (
            <DetailsViewItem
              key={k}
              itemKey={k}
              itemValue={v}
              level={level + 1}
              index={i}
              setItemValue={newVal => setItemValue({ ...itemValue, [k]: newVal })}
              setItemKey={isDir ? undefined : changeKey}
            />
          )
        })}
      </TreeItem>
    )
    : (
      <TreeItem
        key={key}
        itemId={key}
        label={(
          <Stack direction={'row'} alignItems={'center'}>
            <DetailsViewItemSingle
              itemKey={itemKey}
              itemValue={itemValue}
              setItemValue={newVal => setItemValue(newVal)}
            />
            {setItemKey &&
              <ButtonTooltip
                title={'Delete ' + itemKey}
                onClick={() => setItemKey(undefined)}
              >
                <Delete />
              </ButtonTooltip>
            }
          </Stack>
        )}
      />
    )
}
