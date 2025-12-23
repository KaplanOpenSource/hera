import { Stack } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { DetailsViewItemSingle } from './DetailsViewItemSingle';

export const keyForDetailsViewItem = (itemKey: string, level: number = 0, index: number = 0) => {
  return `___lvl${level}_idx${index}_${itemKey}`;
}

export const DetailsViewItem = ({
  itemKey,
  itemValue,
  setItemValue,
  level = 0,
  index,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  level: number,
  index: number,
}) => {
  const key = keyForDetailsViewItem(itemKey, level, index);

  return typeof itemValue === 'object'
    ? (
      <TreeItem
        key={key}
        itemId={key}
        label={<Stack direction={'row'} alignItems={'center'}>
          {itemKey}
        </Stack>}
      >
        {Object.entries(itemValue).map(([k, v], i) => (
          <DetailsViewItem
            key={k}
            itemKey={k}
            itemValue={v}
            level={level + 1}
            index={i}
            setItemValue={newVal => setItemValue({ ...itemValue, [k]: newVal })}
          />
        ))}
      </TreeItem>
    )
    : (
      <TreeItem
        key={key}
        itemId={key}
        label={(
          <Stack direction={'row'}>
            <DetailsViewItemSingle
              itemKey={itemKey}
              itemValue={itemValue}
              setItemValue={newVal => setItemValue(newVal)}
            />
          </Stack>
        )}
      />
    )
}
