import { Grid, Typography } from '@mui/material';
import { ProjectDocument } from '../../shared/types';
import { VersionFields } from './VersionFields';

export const DetailsViewDocumentHeader = ({
  docid,
  shownDoc,
  setShownDoc,
  showFormulated,
}: {
  docid: string,
  shownDoc: ProjectDocument,
  setShownDoc: (v: ProjectDocument) => void,
  showFormulated: boolean,
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
          {docid}
        </Typography>
      </Grid>
      <Grid size={1}>
        <Typography sx={{ fontSize: 12 }}>
          Cls:
        </Typography>
      </Grid>
      <Grid size={11}>
        <Typography sx={{ fontSize: 12 }}>
          {shownDoc._cls}
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
            {shownDoc.desc.toolkit || 'None'}
          </Typography>
        </Grid>
        {!shownDoc.desc.version ? null : (<>
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

