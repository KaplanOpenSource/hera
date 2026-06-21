import { useEffect, useState } from "react";
import { fetchPython } from "../../io/fetchPython";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { RepoTreeDisplay } from "./RepoTreeDisplay";

export const DetailsViewRepo = ({
  repoPath,
}: {
  repoPath: string,
}) => {
  const isTempRepo = repoPath.includes(TEMP_REPO_NAME);
  const [tree, setTree] = useState<any>(undefined);

  useEffect(() => {
    (async () => {
      if (!isTempRepo) {
        const isFilePath = repoPath.includes('/') || repoPath.endsWith('.json');
        const { data } = await fetchPython({
          results: ['jsonData'],
          label: 'repo JSON',
          code: isFilePath
            ? `
import json
with open('${repoPath}', 'r') as fjson:
  jsonData = json.load(fjson)
`
            : `
from hera.utils.data.toolkit import dataToolkit
jsonData = dataToolkit().getRepository('${repoPath}')
`,
        });
        setTree(data?.jsonData ?? undefined);
      }
    })()
  }, [repoPath])

  return (
    <RepoTreeDisplay
      tree={tree}
      setTree={setTree}
      label={repoPath.split('/').pop() ?? repoPath}
      itemId={idRepoId(repoPath)}
      showStrDefault={isTempRepo}
    />
  );
};
