import { Grid, Typography } from '@mui/material';
import { DocumentObj } from '../../objects/ProjectObj';
import { VersionFields } from './VersionFields';

export const DetailsViewDocumentHeader = ({
  doc,
  showFormulated,
  shownDoc,
  setShownDoc
}: {
  doc: DocumentObj,
  showFormulated: boolean,
  shownDoc: any,
  setShownDoc: (v: any) => void,
}) => {
  return (
    <Grid container spacing={1} alignItems={'center'}>
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
        {!doc.data.desc.version ? null : (<>
          <Grid size={1}>
            <Typography sx={{ fontSize: 12 }}>
              Version:
            </Typography>
          </Grid>
          <Grid size={11}>
            <VersionFields
              projectDoc={shownDoc}
              setProjectDoc={setShownDoc} />
          </Grid>
        </>)}
      </>)}
    </Grid>
  )
}

