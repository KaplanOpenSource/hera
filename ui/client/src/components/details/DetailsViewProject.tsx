import { Article, ReceiptLong } from "@mui/icons-material"
import { Typography, Stack } from "@mui/material"
import { ProjectEntire } from "@shared/types"

export const DetailsViewProject = ({
  project,
}: {
  project: ProjectEntire,
}) => {
  return (<>
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
  </>)
}