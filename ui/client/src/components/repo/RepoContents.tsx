import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { JsonTreeNode } from "../../elements/JsonTreeNode";
import { fetchPython } from "../../io/fetchPython";
import { Repository } from "../../shared/types";
import { RepoLabel } from "./RepoLabel";

export const RepoContents = ({ repo }: { repo: Repository }) => {
  const [contents, setContents] = useState<Record<string, any> | null>(null);
  const [loaded, setLoaded] = useState(false);

  const handleExpand = async () => {
    if (loaded) return;
    const { data, problem } = await fetchPython({
      results: ['repoData'],
      code: `
from hera.utils.data.toolkit import dataToolkit
repoData = dataToolkit().getRepository('${repo.datasourceName}')
`,
    });
    if (!problem && data?.repoData) {
      setContents(data.repoData);
    }
    setLoaded(true);
  };

  return (
    <TreeItem
      itemId={repo.datasourceName}
      label={<RepoLabel repo={repo} />}
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
