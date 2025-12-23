import { Article, ReceiptLong } from "@mui/icons-material"
import { Stack, Typography } from "@mui/material"
import { ProjectObj } from "../../objects/ProjectObj"

export const DetailsViewProject = ({
  project,
}: {
  project: ProjectObj,
}) => {
  return (<>
    <Typography variant="h6">{project.name}</Typography>
    {/* <Typography>ID: {project.name}</Typography> */}
    <Typography>Documents: {project.documents?.length ?? 'N/A'}</Typography>
    {project.allDocuments?.map(d => (
      <Stack key={d.docid} direction='row'>
        {!d.isConfig
          ? (<>
            <Article />
            <Typography>
              {d.name}
            </Typography>
          </>)
          : (<>
            <ReceiptLong />
            <Typography>
              {d.name ?? 'N/A'}
            </Typography>
          </>)}
      </Stack>
    ))}
    {/* <Typography>Cache Documents: {project.documents?.cache ?? 'N/A'}</Typography>
      <Typography>Toolkits: {project.toolkitCount ?? 'N/A'}</Typography> */}
  </>)
}