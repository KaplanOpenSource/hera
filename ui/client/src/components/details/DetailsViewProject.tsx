import { Article, ReceiptLong } from "@mui/icons-material"
import { Stack, Typography } from "@mui/material"
import { ProjectObj } from "../../objects/ProjectObj"
import { DetailsViewDocument } from "./DetailsViewDocument"

export const DetailsViewProject = ({
  project,
}: {
  project: ProjectObj,
}) => {
  const doc = project.configDocument;
  return (<>
    <Typography variant="h6">{project.name}</Typography>
    {/* <Typography>ID: {project.name}</Typography> */}
    {doc ? <>
      <DetailsViewDocument
        doc={doc.data}
        setDoc={(newVal) => {
          // changeDocument(newVal)
        }}
      />
    </>
      : <Typography>
        No config document
      </Typography>}
    <Typography>Documents: {project.documents?.length ?? 'N/A'}</Typography>
    {project.documents.map(d => (
      <Stack key={d.docid} direction='row'>
        <Article />
        <Typography>
          {d.name}
        </Typography>
      </Stack>
    ))}
    {/* <Typography>Cache Documents: {project.documents?.cache ?? 'N/A'}</Typography>
      <Typography>Toolkits: {project.toolkitCount ?? 'N/A'}</Typography> */}
  </>)
}