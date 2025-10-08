import { ProjectEntire, ProjectName } from "@shared/types";
import { execPython } from "./execPython";
import { useEffect } from "react";
import { NO_PROJECT, useProjectStore } from "../stores/useProjectStore";

export const fetchProjectNames = async (): Promise<{ data: ProjectName[] | undefined; problem: string | undefined; }> => {
  const { data, problem } = await execPython(`
from hera.datalayer.project import getProjectList;
result = [{"name": proj} for i, proj in enumerate(getProjectList())]
          `)
  return { data: data as ProjectName[] | undefined, problem };
}


export const FetcherProjectNames = ({ }) => {
  const { projectNames, setProjectNames, setCurrentProject, currProjectName, currProject } = useProjectStore();

  useEffect(() => {
    (async () => {
      console.log('loading project names started');
      const { data, problem } = await fetchProjectNames();
      if (data) {
        setProjectNames(data || []);
      } else {
        console.log('problem loading:', problem);
        setProjectNames([]);
      }
      console.log('loading project names done');
    })();
  }, []);

  useEffect(() => {
    (async () => {
      if (currProjectName !== NO_PROJECT && currProjectName !== currProject?.name) {
        console.log(`loading project ${currProjectName} started`);
        if (currProjectName) {
          const { data, problem } = await execPython(`
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
          `)
          if (data) {
            console.log(data);
            const project = data as ProjectEntire;
            setCurrentProject(project);
          } else {
            console.log('problem loading:', problem);
            setProjectNames([]);
          }
        }
        console.log(`loading project ${currProjectName} done`);
      }
    })();
  }, [currProjectName, currProject])

  return null;
}

