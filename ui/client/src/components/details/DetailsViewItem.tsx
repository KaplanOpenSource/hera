import { Add, CreateNewFolder, Delete, Edit } from '@mui/icons-material';
import { Box, Stack, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { RenameField } from '../../elements/RenameField';
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
  const [editing, setEditing] = useState(false);

  const key = keyForDetailsViewItem(itemKey, level, index);
  const isTree = typeof itemValue === 'object';

  const addSubItem = (initialValue: any) => {
    let name = '';
    for (let i = 1; i < 1e5; i++) {
      const key = 'newItem_' + i;
      if (!(key in itemValue)) {
        name = key;
        break;
      }
    }
    setItemValue({ ...itemValue, [name]: initialValue })
  }

  return (
    <TreeItem
      key={key}
      itemId={key}
      label={(
        <Stack
          direction='row'
          spacing={1}
          justifyItems={'stretch'}
          alignItems={'center'}
          style={{ marginTop: 7 }}
        >

          <Box sx={{ marginRight: 1 }}>
            {(editing && setItemKey)
              ? (
                <RenameField
                  value={itemKey}
                  setValue={val => {
                    setItemKey(val);
                    setEditing(false);
                  }}
                />
              )
              : (
                <Typography
                  sx={{
                    whiteSpace: 'nowrap',
                    minWidth: '100px'
                  }}
                >
                  {itemKey}
                </Typography>
              )}
          </Box>

          {isTree && (<>
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
          </>)}

          {!isTree && (<>
            <DetailsViewItemSingle
              itemKey={itemKey}
              itemValue={itemValue}
              setItemValue={newVal => setItemValue(newVal)}
            />
          </>)}

          {setItemKey && (<>
            <ButtonTooltip
              title={'Rename ' + itemKey}
              onClick={() => setEditing(!editing)}
            >
              <Edit />
            </ButtonTooltip>
            <ButtonTooltip
              title={'Delete ' + itemKey}
              onClick={() => setItemKey(undefined)}
            >
              <Delete />
            </ButtonTooltip>
          </>)}
        </Stack>
      )}
    >
      {isTree && (<>
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
      </>)}
    </TreeItem>
  )
}
