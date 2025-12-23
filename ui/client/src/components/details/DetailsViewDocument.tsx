import { Close, Done } from '@mui/icons-material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';
import { DocumentObj } from '../../objects/ProjectObj';
import { Stack, Typography } from '@mui/material';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc.data)));

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  const isChanged = JSON.stringify(doc.data) !== JSON.stringify(shownDoc);
  return (
    <>
      <Stack direction={'row'} alignItems={'center'}>
        <Typography variant='h6'>
          {doc.isConfig ? doc.project.name + ' config' : doc.name}
        </Typography>
        {isChanged
          ? (<>
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
          </>)
          : null}
      </Stack>
      <SimpleTreeView
        defaultExpandedItems={[keyForDetailsViewItem('desc', 1, 3)]}
      >
        {Object.entries(shownDoc).map(([k, v], i) => {
          if (doc.isConfig && ['projectName', '_id'].includes(k)) {
            return null;
          }
          return (
            <DetailsViewItem
              key={k}
              itemKey={k}
              itemValue={v}
              level={1}
              index={i}
              setItemValue={newVal => setShownDoc({ ...shownDoc, [k]: newVal })}
            />
          )
        })}
      </SimpleTreeView>
    </>
  )
}
