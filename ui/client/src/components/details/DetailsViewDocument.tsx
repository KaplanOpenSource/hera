import { FormControl, InputLabel, MenuItem, Select, Stack, TextField, Typography } from '@mui/material';
import { SimpleTreeView, TreeItem } from '@mui/x-tree-view';
import { Case, SwitchCase } from '../../elements/SwitchCase';
import { useServerConstants } from '../../stores/useServerConstants';

const DetailsViewItemSingle = ({
  itemKey,
  itemValue,
}: {
  itemKey: string,
  itemValue: any,
}) => {
  const { dataTypes } = useServerConstants();
  return (
    <Stack direction='row' spacing={1} justifyItems={'center'} alignItems={'center'}>
      <Typography>
        {itemKey}
      </Typography>
      <SwitchCase test={itemKey}>
        <Case isDefault>
          <TextField
            size='small'
            value={
              JSON.stringify(itemValue)
            }
          />
        </Case>
        <Case value={'dataFormat'}>
          <FormControl style={{ marginTop: 10, minWidth: '100px' }}>
            <InputLabel id="demo-simple-select-label">dataFormat</InputLabel>
            <Select
              labelId="demo-simple-select-label"
              id="demo-simple-select"
              value={itemValue}
              label="dataFormat"
              size='small'
            // onChange={(event: SelectChangeEvent) => setAge(event.target.value as string)}
            >
              {Object.entries(dataTypes).map(([_upcasename, name]) => (
                <MenuItem key={name} value={name}>{name}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Case>
      </SwitchCase>
    </Stack>
  )
}

const DetailsViewItem = ({
  itemKey,
  itemValue,
  level = 0,
  index,
}: {
  itemKey: string,
  itemValue: any,
  level: number,
  index: number,
}) => {
  const key = `___lvl${level}_idx${index}_${itemKey}`
  return typeof itemValue === 'object'
    ? (
      <TreeItem
        key={key}
        itemId={key}
        label={itemKey}
      >
        {Object.entries(itemValue).map(([k, v], i) => (
          <DetailsViewItem
            key={k}
            itemKey={k}
            itemValue={v}
            level={level + 1}
            index={i}
          />
        ))}
      </TreeItem>
    )
    : (
      <TreeItem
        key={key}
        itemId={key}
        label={(
          <DetailsViewItemSingle
            itemKey={itemKey}
            itemValue={itemValue}
          />
        )}
      />
    )
}

export const DetailsViewDocument = ({
  doc,
}: {
  doc: any,
}) => {
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;
  const items: [string, any][] = Object.entries(doc);

  return (
    <SimpleTreeView>
      <DetailsViewItem
        itemKey={name}
        itemValue={doc}
        level={0}
        index={0}
      />
    </SimpleTreeView>
  )
}