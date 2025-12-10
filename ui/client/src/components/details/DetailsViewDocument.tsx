import { SimpleTreeView, TreeItem } from '@mui/x-tree-view';
import { DetailsViewItemSingle } from './DetailsViewItemSingle';
import { useState } from 'react';

const DetailsViewItem = ({
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
          <DetailsViewItemSingle
            itemKey={itemKey}
            itemValue={itemValue}
            setItemValue={newVal => setItemValue(newVal)}
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
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc)));
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;

  console.log(shownDoc)
  return (
    <SimpleTreeView>
      <DetailsViewItem
        itemKey={name}
        itemValue={shownDoc}
        level={0}
        index={0}
        setItemValue={newVal => setShownDoc(newVal)}
      />
    </SimpleTreeView>
  )
}