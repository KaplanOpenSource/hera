import { Article, ReceiptLong } from '@mui/icons-material';
import { Paper, Stack, Typography } from '@mui/material';
import type { ProjectEntire } from '../../shared/types';
import { useEffect, useState } from 'react';
import { idFromDocId } from '../../shared/idDocId';
import { execPython } from '../../io/execPython';

export const DetailsViewPanel = ({
  project,
  selectedItemsIds,
}: {
  project: ProjectEntire,
  selectedItemsIds: string[],
}) => {
  const docid = idFromDocId(selectedItemsIds[0]);
  const [doc, setDoc] = useState<any>(undefined);

  useEffect(() => {
    (async () => {
      console.log(docid);
      if (docid) {
        const { data } = await execPython(`
import json
from hera.datalayer import All
docs = All.getDocumentByID('${docid}')
result = docs.asDict(with_id=True)
`);
        if (data) (
          setDoc(data)
        )
      } else {
        setDoc(undefined);
      }
    })()
  }, [docid])

  return (
    <Paper sx={{ p: 2, height: '100%' }}>
      {doc
        ? <pre>{JSON.stringify(doc, null, 2)}</pre>
        : (<>
          <Typography variant="h6">{project.name}</Typography>
          {/* <Typography>ID: {project.name}</Typography> */}
          <Typography>Documents: {project.documents?.length ?? 'N/A'}</Typography>
          {project.documents?.map(d => (
            <Stack key={d._id.$oid} direction='row'>
              {d.desc.datasourceName
                ? (<>
                  <Article />
                  <Typography>
                    {d.desc.datasourceName}
                  </Typography>
                </>)
                : (<>
                  <ReceiptLong />
                  <Typography>
                    {d.type ?? 'N/A'}
                  </Typography>
                </>)}
            </Stack>
          ))}
          {/* <Typography>Cache Documents: {project.documents?.cache ?? 'N/A'}</Typography>
      <Typography>Toolkits: {project.toolkitCount ?? 'N/A'}</Typography> */}
        </>)}
    </Paper>
  );
};
