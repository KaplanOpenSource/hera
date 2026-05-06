import { Edit } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { ReactNode, useCallback, useEffect, useMemo, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { JsonTreeNode } from "../../elements/JsonTreeNode";
import { SECTION_DATASOURCE, SECTION_EXPERIMENT } from "./RepoJsonMerger";
import { RepoTreeAddButton } from "./RepoTreeAddButton";
import { TextFieldWithApply } from "./TextFieldWithApply";
import { VisibilityToggleButton } from "./VisibilityToggleButton";

export const RepoTreeDisplay = ({
  tree,
  setTree,
  label,
  itemId,
  tooltipContent,
  showStrDefault = false,
  overrides,
  repoJsons,
}: {
  tree: any;
  setTree: (next: any) => void;
  label: string;
  itemId: string;
  tooltipContent?: ReactNode;
  showStrDefault?: boolean;
  overrides?: { [path: string]: string[] };
  repoJsons?: { [path: string]: { [key: string]: any } };
}) => {
  const [repoStr, setRepoStr] = useState<string>('');
  const [showStr, setShowStr] = useState(showStrDefault);
  const [hiddenPaths, setHiddenPaths] = useState<Set<string>>(new Set());
  const [showHidden, setShowHidden] = useState(false);

  const canToggleHidden = useCallback((treePath: string[]) => {
    if (treePath.length === 1) return true;
    if (treePath[0].toLowerCase() === SECTION_EXPERIMENT) {
      return treePath.length === 2;
    }
    if (treePath.length === 3 && treePath[1].toLowerCase() === SECTION_DATASOURCE) {
      return true;
    }
    return false;
  }, []);

  const onToggleHidden = (path: string) => {
    setHiddenPaths((prev) => {
      const next = new Set(prev);
      if (next.has(path)) {
        next.delete(path);
      } else {
        next.add(path);
      }
      return next;
    });
  };

  useEffect(() => {
    const newText = JSON.stringify(tree, undefined, 2);
    if (newText !== repoStr) {
      setRepoStr(JSON.stringify(tree, undefined, 2));
    }
  }, [tree]);

  const expandedByDefault = useMemo(() => {
    const ids = [itemId];
    if (typeof tree === 'object' && tree !== null) {
      for (const key of Object.keys(tree)) {
        ids.push(`${itemId}.${key}.${key}`);
      }
    }
    return ids;
  }, [tree, itemId]);

  return (
    <>
      <SimpleTreeView
        defaultExpandedItems={expandedByDefault}
        key={expandedByDefault.join(',')}
      >
        <TreeItem
          itemId={itemId}
          label={(
            <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
              <Typography marginRight={1}>
                {label}
              </Typography>
              <RepoTreeAddButton tree={tree} hiddenPaths={hiddenPaths} />
              <ButtonTooltip
                title={'Edit Repository as Json text'}
                onClick={() => setShowStr(!showStr)}
              >
                <Edit />
              </ButtonTooltip>
              <VisibilityToggleButton
                hasHidden={hiddenPaths.size > 0}
                showHidden={showHidden}
                setShowHidden={setShowHidden}
                onRestoreAll={() => setHiddenPaths(new Set())}
              />
            </Stack>
          )}
        >
          {typeof tree === "object" && tree !== null
            ? Object.entries(tree).map(([key, val]) => (
              <JsonTreeNode
                key={key}
                label={key}
                parentKey={`${itemId}.${key}`}
                value={val as any}
                overrides={overrides}
                treePath={[key]}
                repoJsons={repoJsons}
                hiddenPaths={hiddenPaths}
                showHidden={showHidden}
                onToggleHidden={onToggleHidden}
                canToggleHidden={canToggleHidden}
              />
            ))
            : null}
        </TreeItem>
      </SimpleTreeView>
      {showStr
        ? (
          <TextFieldWithApply
            value={repoStr}
            setValue={v => setRepoStr(v)}
            onApply={() => {
              const t = JSON.parse(repoStr);
              setTree(t);
            }}
            textTitle="Repository Json (as text)"
            buttonTitle={'Apply Json text to tree (locally, not affecting the file)'}
          />
        )
        : null}
    </>
  );
};
