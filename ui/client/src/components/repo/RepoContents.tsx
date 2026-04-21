import { Stack } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ReactNode, useState } from "react";
import { JsonTreeNode } from "../../elements/JsonTreeNode";
import { fetchPython } from "../../io/fetchPython";
import { Repository } from "../../shared/types";
import { RepoLabel } from "./RepoLabel";

export const RepoContents = ({
  repo,
  actions,
}: {
  repo: Repository,
  actions?: ReactNode,
}) => {
  const [contents, setContents] = useState<Record<string, any> | null>(null);
  const [loaded, setLoaded] = useState(false);

  const handleExpand = async () => {
    if (loaded) return;
    const { data } = await fetchPython({
      results: ['repoData'],
      label: `repo ${repo.datasourceName}`,
      code: `
from hera.utils.data.toolkit import dataToolkit
repoData = dataToolkit().getRepository('${repo.datasourceName}')
`,
    });
    if (data?.repoData) {
      setContents(data.repoData);
    }
    setLoaded(true);
  };

  return (
    <TreeItem
      itemId={repo.datasourceName}
      label={actions
        ? (
          <Stack direction="row" alignItems="center">
            <RepoLabel repo={repo} />
            {actions}
          </Stack>
        )
        : (
          <RepoLabel repo={repo} />
        )
      }
      onClick={handleExpand}
    >
      {!loaded
        ? <TreeItem itemId={`${repo.datasourceName}/__loading`} label="Loading..." />
        : contents && Object.entries(contents).map(([key, val]) => (
          <JsonTreeNode
            key={key}
            label={key}
            parentKey={`${repo.datasourceName}.${key}`}
            value={val}
          />
        ))
      }
    </TreeItem>
  )
}
