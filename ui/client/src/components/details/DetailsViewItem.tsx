import { Stack } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { ReactNode } from 'react';
import { DetailsViewItemSingle } from './DetailsViewItemSingle';

export const DetailsViewItem = ({
  itemKey,
  itemValue,
  setItemValue,
  level = 0,
  index,
  labelAdditions = null,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  level: number,
  index: number,
  labelAdditions?: ReactNode,
}) => {
  const key = `___lvl${level}_idx${index}_${itemKey}`
  return typeof itemValue === 'object'
    ? (
      <TreeItem
        key={key}
        itemId={key}
        label={<Stack direction={'row'} alignItems={'center'}>
          {itemKey}
          {labelAdditions}
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
            {labelAdditions}
          </Stack>
        )}
      />
    )
}
