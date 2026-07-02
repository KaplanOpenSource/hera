import { Add, CreateNewFolder, Delete } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DetailsViewItemName } from './DetailsViewItemName';
import { DetailsViewItemSingle } from './DetailsViewItemSingle';
import { FieldDef } from './fieldDef';
import { FieldSourceDot } from './FieldSourceDot';
import { EditAsJsonButton } from './EditAsJsonButton';
import { SelectDataFormat } from './SelectDataFormat';

export const keyForDetailsViewItem = (itemKey: string, parentKey?: string) => {
  return parentKey ? `${parentKey}/${itemKey}` : itemKey;
};

export const DetailsViewItem = ({
  itemKey,
  itemValue,
  setItemValue,
  setItemKey = undefined,
  parentKey,
  def = undefined,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  setItemKey?: (newKey: string | undefined) => void | undefined,
  parentKey?: string,
  // Definition of this field: `required` for the value editor, `children` for
  // the sub-fields below it.
  def?: FieldDef,
}) => {
  const key = keyForDetailsViewItem(itemKey, parentKey);
  const isTree = typeof itemValue === 'object' && itemValue !== null;
  const level = parentKey?.split('/').length || 0;

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
          // Bottom space on every row reserves room for a field's "required"
          // helper text, so it shows without moving anything.
          style={{ marginTop: 7, marginBottom: 14 }}
        >

          <FieldSourceDot source={def?.source} />

          <DetailsViewItemName
            itemKey={itemKey}
            setItemKey={setItemKey}
          />

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
            <EditAsJsonButton
              data={itemValue}
              setData={setItemValue}
            />
          </>)}

          {isTree
            ? null
            : (itemKey === 'dataFormat'
              ? (
                <SelectDataFormat
                  value={itemValue}
                  setValue={v => setItemValue(v)}
                />
              )
              : (
                <DetailsViewItemSingle
                  itemValue={itemValue}
                  setItemValue={newVal => setItemValue(newVal)}
                  def={def}
                />
              )
            )
          }

          {setItemKey && (<>
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
      {isTree && Object.keys(itemValue).length === 0 && (
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
      )}
      {isTree && (<>
        {Object.entries(itemValue).sort().map(([k, v]) => {
          const isDir = parentKey === undefined && itemKey === 'desc' && k === 'filesDirectory';

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
              parentKey={key}
              setItemValue={newVal => setItemValue({ ...itemValue, [k]: newVal })}
              setItemKey={isDir ? undefined : changeKey}
              def={def?.children?.[k]}
            />
          )
        })}
      </>)}
    </TreeItem>
  )
}
