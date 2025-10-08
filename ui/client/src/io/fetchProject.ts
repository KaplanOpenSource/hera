import { ProjectName } from "@shared/types";
import { execPython } from "./execPython";
import { useEffect } from "react";
import { useProjectStore } from "../stores/useProjectStore";

export const fetchProjectNames = async (): Promise<{ data: ProjectName[] | undefined; problem: string | undefined; }> => {
  const { data, problem } = await execPython(`
from hera.datalayer.project import getProjectList;
result = [{"id": "p-" + str(i), "name": proj} for i, proj in enumerate(getProjectList())]
          `)
  return { data: data as ProjectName[] | undefined, problem };
}


export const FetcherProjectNames = ({ }) => {
  const { setProjectNames } = useProjectStore();

  useEffect(() => {
    console.log('loading started');
    (async () => {
      const { data, problem } = await fetchProjectNames();
      if (data) {
        setProjectNames(data || []);
      } else {
        console.log('problem loading:', problem);
        setProjectNames([]);
      }
      console.log('loading done');
    })();
  }, []);

  return null;
}

