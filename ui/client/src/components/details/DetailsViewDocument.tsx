import { Close, Done } from '@mui/icons-material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';
import { DocumentObj } from '../../objects/ProjectObj';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc.data)));

  const isChanged = JSON.stringify(doc.data) === JSON.stringify(shownDoc);
  return (
    <SimpleTreeView
      defaultExpandedItems={[keyForDetailsViewItem(doc.name)]}
    >
      <DetailsViewItem
        itemKey={doc.name}
        itemValue={shownDoc}
        level={0}
        index={0}
        setItemValue={newVal => setShownDoc(newVal)}
        labelAdditions={isChanged
          ? null
          : (
            <>
              <ButtonTooltip
                title='Update Document'
                onClick={() => setDoc(new DocumentObj(shownDoc, doc.project))}
              >
                <Done />
              </ButtonTooltip>
              <ButtonTooltip
                title='Revert Document'
                onClick={() => setShownDoc(JSON.parse(JSON.stringify(doc.data)))}
              >
                <Close />
              </ButtonTooltip>
            </>
          )
        }
      />
    </SimpleTreeView>
  )
}
