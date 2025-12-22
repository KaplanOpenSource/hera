import { Close, Done } from '@mui/icons-material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: any,
  setDoc: (newVal: any) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc)));
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;

  return (
    <SimpleTreeView
      defaultExpandedItems={[keyForDetailsViewItem(name)]}
    >
      <DetailsViewItem
        itemKey={name}
        itemValue={shownDoc}
        level={0}
        index={0}
        setItemValue={newVal => setShownDoc(newVal)}
        labelAdditions={JSON.stringify(doc) === JSON.stringify(shownDoc) ? null : (
          <>
            <ButtonTooltip
              title='Update Document'
              onClick={() => setDoc(shownDoc)}
            >
              <Done />
            </ButtonTooltip>
            <ButtonTooltip
              title='Revert Document'
              onClick={() => setShownDoc(JSON.parse(JSON.stringify(doc)))}
            >
              <Close />
            </ButtonTooltip>
          </>
        )}
      />
    </SimpleTreeView>
  )
}
