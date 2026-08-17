import { ExpandMore } from '@mui/icons-material';
import { Accordion, AccordionDetails, AccordionSummary, Grid, Typography } from '@mui/material';
import { Fragment, useState } from 'react';
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
  const [open, setOpen] = useState(true);
  return (
    <Accordion
      expanded={open}
      onChange={() => setOpen(!open)}
      disableGutters
      elevation={0}
      sx={{
        mt: 1,
        width: 'fit-content',
        maxWidth: '100%',
        bgcolor: 'background.paper',
        border: 1,
        borderColor: 'divider',
        borderRadius: 1,
        '&:before': { display: 'none' },
      }}
    >
      <AccordionSummary
        expandIcon={<ExpandMore />}
        sx={{
          minHeight: 0,
          px: 1,
          gap: 0.5,
          flexDirection: 'row-reverse',
          '& .MuiAccordionSummary-content': { my: 0.75 },
        }}
      >
        <Typography variant="overline" sx={{ fontSize: 10, color: 'text.secondary', fontWeight: 600, letterSpacing: 1 }}>
          Node Metadata Attributes
        </Typography>
      </AccordionSummary>
      <AccordionDetails sx={{ px: 1, pt: 0, pb: 1 }}>
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
          {showFormulated && (<>
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
            {shownDoc.desc.version && (<>
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
      </AccordionDetails>
    </Accordion>
  )
}
