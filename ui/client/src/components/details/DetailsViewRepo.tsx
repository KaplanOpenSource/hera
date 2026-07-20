import { Link, LinkOff } from "@mui/icons-material";
import { useEffect, useState } from "react";
import { fetchPython } from "../../io/fetchPython";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { RepoTreeDisplay } from "./RepoTreeDisplay";

export const DetailsViewRepo = ({
  repoPath,
}: {
  repoPath: string,
}) => {
  const isTempRepo = repoPath.includes(TEMP_REPO_NAME);
  const [tree, setTree] = useState<any>(undefined);
  const [origTree, setOrigTree] = useState<any>(undefined);
  const [peekOriginal, setPeekOriginal] = useState(false);

  useEffect(() => {
    (async () => {
      if (!isTempRepo) {
        const isFilePath = repoPath.includes('/') || repoPath.endsWith('.json');
        const { data } = await fetchPython({
          results: ['jsonData', 'jsonDataResolved'],
          label: 'repo JSON',
          code: isFilePath
            ? `
import json, os
from hera.utils.data.toolkit import dataToolkit
with open('${repoPath}', 'r') as fjson:
  jsonData = json.load(fjson)
jsonDataResolved = dataToolkit.resolveDataSourcePaths(jsonData, os.path.dirname('${repoPath}'))
`
            : `
from hera.utils.data.toolkit import dataToolkit
jsonData = dataToolkit().getRepository('${repoPath}')
jsonDataResolved = jsonData
`,
        });
        setOrigTree(data?.jsonData ?? undefined);
        setTree(data?.jsonDataResolved ?? data?.jsonData ?? undefined);
      }
    })()
  }, [repoPath])

  const pathsResolved = JSON.stringify(origTree) !== JSON.stringify(tree);
  const displayedTree = (peekOriginal && origTree) ? origTree : tree;

  return (
    <RepoTreeDisplay
      tree={displayedTree}
      setTree={setTree}
      label={repoPath.split('/').pop() ?? repoPath}
      itemId={idRepoId(repoPath)}
      showStrDefault={isTempRepo}
      labelIcons={pathsResolved && (
        <ButtonTooltip
          title={peekOriginal ? 'Original (relative) paths' : 'Resolved (absolute) paths — hover to see original'}
          onClick={() => { }}
          onMouseEnter={() => setPeekOriginal(true)}
          onMouseLeave={() => setPeekOriginal(false)}
        >
          {peekOriginal ? <LinkOff fontSize="small" /> : <Link fontSize="small" />}
        </ButtonTooltip>
      )}
    />
  );
};
