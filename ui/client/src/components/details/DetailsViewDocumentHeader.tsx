import { Grid, Typography } from '@mui/material';
import { Fragment } from 'react';
import { ProjectDocument } from '../../shared/types';
import { VersionFields } from './VersionFields';
export const DetailsViewDocumentHeader = ({
  docid,
  shownDoc,
  setShownDoc,
  showFormulated,
  extraFields = [],
}: {
  docid: string,
  shownDoc: ProjectDocument,
  setShownDoc: (v: ProjectDocument) => void,
  showFormulated: boolean,
  extraFields?: { name: string, value: string }[],
}) => {
  return (
    <Grid container spacing={1} alignItems={'center'}>
      <Grid key="id-label" size={2}>
        <Typography sx={{ fontSize: 12 }}>
          Id:
        </Typography>
      </Grid>
      <Grid key="id-value" size={10}>
        <Typography sx={{ fontSize: 12 }}>
          {docid}
        </Typography>
      </Grid>
      <Grid key="cls-label" size={2}>
        <Typography sx={{ fontSize: 12 }}>
          Cls:
        </Typography>
      </Grid>
      <Grid key="cls-value" size={10}>
        <Typography sx={{ fontSize: 12 }}>
          {shownDoc._cls}
        </Typography>
      </Grid>
      {!showFormulated ? null : (<>
        <Grid key="toolkit-label" size={2}>
          <Typography sx={{ fontSize: 12 }}>
            Toolkit:
          </Typography>
        </Grid>
        <Grid key="toolkit-value" size={10}>
          <Typography sx={{ fontSize: 12 }}>
            {shownDoc.desc.toolkit || 'None'}
          </Typography>
        </Grid>
        {!shownDoc.desc.version ? null : (<>
          <Grid key="version-label" size={2}>
            <Typography sx={{ fontSize: 12 }}>
              Version:
            </Typography>
          </Grid>
          <Grid key="version-value" size={10}>
            <VersionFields
              projectDoc={shownDoc}
              setProjectDoc={setShownDoc} />
          </Grid>
        </>)}
      </>)}
      {extraFields.map(({ name, value }) => (
        <Fragment key={`extra-${name}`}>
          <Grid size={2}>
            <Typography sx={{ fontSize: 12 }}>
              {name}
            </Typography>
          </Grid>
          <Grid size={10}>
            <Typography sx={{ fontSize: 12 }}>
              {value}
            </Typography>
          </Grid>
        </Fragment>
      ))}
    </Grid>
  )
}
