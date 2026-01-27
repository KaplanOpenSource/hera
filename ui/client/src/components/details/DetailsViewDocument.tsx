import { Close, Done, DynamicForm } from '@mui/icons-material';
import { Grid, Stack, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj } from '../../objects/ProjectObj';
import { FORBIDDEN_FIELDS } from '../../shared/constants';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';
import { copyWithout, reorderEntries } from '../../utils/utils';

const HIDE_ON_DESC = ['datasourceName', 'toolkit', 'version'];

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc.data)));
  const [showFormulated, setShowFormulated] = useState(false);

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  const isChanged = JSON.stringify(doc.data) !== JSON.stringify(shownDoc);
  return (
    <>
      <Stack direction={'row'} alignItems={'center'} justifyItems={'center'}>
        <Typography variant='h6' sx={{ marginRight: 1 }}>
          {doc.isConfig ? doc.project.name + ' config' : doc.name}
        </Typography>
        <ButtonTooltip
          title={'Show Formulated'}
          onClick={() => setShowFormulated(!showFormulated)}
        >
          <DynamicForm color={showFormulated ? 'primary' : 'inherit'} />
        </ButtonTooltip>
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
      <Grid container spacing={1}>
        <Grid size={1}>
          <Typography sx={{ fontSize: 12 }}>
            Id:
          </Typography>
        </Grid>
        <Grid size={11}>
          <Typography sx={{ fontSize: 12 }}>
            {doc.docid}
          </Typography>
        </Grid>
        <Grid size={1}>
          <Typography sx={{ fontSize: 12 }}>
            Cls:
          </Typography>
        </Grid>
        <Grid size={11}>
          <Typography sx={{ fontSize: 12 }}>
            {doc.data._cls}
          </Typography>
        </Grid>
        {!showFormulated ? null : (<>
          <Grid size={1}>
            <Typography sx={{ fontSize: 12 }}>
              Toolkit:
            </Typography>
          </Grid>
          <Grid size={11}>
            <Typography sx={{ fontSize: 12 }}>
              {doc.data.desc.toolkit || 'None'}
            </Typography>
          </Grid>
          <Grid size={1}>
            <Typography sx={{ fontSize: 12 }}>
              Version:
            </Typography>
          </Grid>
          <Grid size={11}>
            <Typography sx={{ fontSize: 12 }}>
              {(doc.data.desc.version || []).join('.')}
            </Typography>
          </Grid>
        </>)}
      </Grid>
      <SimpleTreeView
        defaultExpandedItems={[keyForDetailsViewItem('desc', 1, 3)]}
      >
        {reorderEntries(Object.entries(shownDoc), ['desc', 'resource']).map(([k, v], i) => {
          if (FORBIDDEN_FIELDS.includes(k)) {
            return null;
          }
          return (
            <DetailsViewItem
              key={k}
              itemKey={k}
              itemValue={!showFormulated || k !== 'desc' ? v : copyWithout(v, HIDE_ON_DESC)}
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
