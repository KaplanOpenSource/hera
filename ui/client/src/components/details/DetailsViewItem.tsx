import { Add, CreateNewFolder, Delete, Keyboard } from '@mui/icons-material';
import { Stack, TextField } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useDialog } from '../../elements/useDialog';
import { DetailsViewItemName } from './DetailsViewItemName';
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
  const isTree = typeof itemValue === 'object';
  const { DialogComponent, openDialog } = useDialog();

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
            <ButtonTooltip
              title={'Edit as Json Text'}
              onClick={async () => {
                const { confirmed, values } = await openDialog({
                  title: 'Edit as Json Text',
                  initialValues: { jsonText: JSON.stringify(itemValue, null, 2), ok: true },
                  render: ({ values, setValues }) => (
                    <TextField
                      label="Json"
                      fullWidth
                      multiline
                      minRows={4}
                      maxRows={20}
                      value={values.jsonText}
                      onChange={e => {
                        e.stopPropagation();
                        let ok = true;
                        try { JSON.parse(e.target.value) } catch { ok = false }
                        setValues({ ...values, jsonText: e.target.value, ok });
                      }}
                      error={!values.ok}
                      helperText={values.ok ? '' : 'Invalid Json'}
                    />
                  )
                });
                if (confirmed && values && values.ok) {
                  try {
                    const json = JSON.parse(values.jsonText);
                    setItemValue(json);
                  } catch { }
                }
              }}
            >
              <Keyboard />
              {DialogComponent}
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
        {Object.entries(itemValue).sort().map(([k, v], i) => {
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
