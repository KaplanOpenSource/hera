import { ProjectEntire, ProjectName } from "@shared/types";
import { useProjectStore } from "../stores/useProjectStore";
import { FetchPython } from "./FetchPython";

export const FetchProjects = ({ }) => {
  const {
    projectNames,
    setProjectNames,
    setCurrentProject,
    currProjectName,
    currProject,
  } = useProjectStore();

  return (
    <>
      <FetchPython
        code={`
from hera.datalayer.project import getProjectList;
result = [{"name": proj} for proj in getProjectList()]
        `}
        onSuccess={(data: ProjectName[] | undefined) => {
          setProjectNames(data || [])
        }}
      />
      <FetchPython
        code={`
import json
from hera.datalayer import All
docList = []
for doc in All.getDocuments('${currProjectName}'):
    docDict = doc.asDict()
    if ('docid' not in docDict['desc']):
        docDict['desc']['docid'] = str(doc.id)
    docList.append(docDict)
project = {"name": '${currProjectName}', "documents": docList}
result = json.dumps(project,indent=4)
        `}
        onSuccess={(data: any) => {
          const project = JSON.parse(data) as ProjectEntire;
          setCurrentProject(project);
        }}
      />
    </>
  )
}

