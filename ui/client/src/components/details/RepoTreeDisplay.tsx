import { Edit, Visibility, VisibilityOff } from "@mui/icons-material";
import { Stack, Typography } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { ReactNode, useEffect, useMemo, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { JsonTreeNode } from "../../elements/JsonTreeNode";
import { DatasourceAddButton } from "./DatasourceAddButton";
import { ExperimentAddButton } from "./ExperimentAddButton";
import { SECTION_DATASOURCE, SECTION_EXPERIMENT } from "./RepoJsonMerger";
import { ToolkitAddButton } from "./ToolkitAddButton";
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

  const canToggleHidden = (treePath: string[]) => {
    if (treePath.length === 1) return true;
    if (treePath[0].toLowerCase() === SECTION_EXPERIMENT) {
      return treePath.length === 2;
    }
    if (treePath.length === 3 && treePath[1].toLowerCase() === SECTION_DATASOURCE) {
      return true;
    }
    return false;
  };

  const onToggleHidden = (treePathKey: string) => {
    setHiddenPaths((prev) => {
      const next = new Set(prev);
      if (next.has(treePathKey)) {
        next.delete(treePathKey);
      } else {
        next.add(treePathKey);
      }
      return next;
    });
  };

  const isDatasource = (treePath: string[]) =>
    treePath.length === 3 && treePath[1].toLowerCase() === SECTION_DATASOURCE;

  const isExperiment = (treePath: string[]) =>
    treePath.length === 2 && treePath[0].toLowerCase() === SECTION_EXPERIMENT;

  const isToolkit = (treePath: string[]) =>
    treePath.length === 1 && treePath[0].toLowerCase() !== SECTION_EXPERIMENT;

  const nodeActions = (treePath: string[]) => {
    const actions = [];
    if (canToggleHidden(treePath)) {
      const treePathKey = treePath.join('/');
      const isHidden = hiddenPaths.has(treePathKey);
      actions.push(
        <ButtonTooltip
          key="visibility"
          title={isHidden ? 'Show' : 'Hide'}
          onClick={() => onToggleHidden(treePathKey)}
        >
          {isHidden ? <VisibilityOff fontSize="inherit" /> : <Visibility fontSize="inherit" />}
        </ButtonTooltip>
      );
    }
    if (isDatasource(treePath)) {
      actions.push(
        <DatasourceAddButton key="add-ds" tree={tree} treePath={treePath} />
      );
    }
    if (isExperiment(treePath)) {
      actions.push(
        <ExperimentAddButton key="add-exp" tree={tree} treePath={treePath} />
      );
    }
    if (isToolkit(treePath)) {
      actions.push(
        <ToolkitAddButton key="add-toolkit" tree={tree} toolkitName={treePath[0]} />
      );
    }
    return actions.length > 0 ? <>{actions}</> : null;
  };

  useEffect(() => {
    const newText = JSON.stringify(tree, undefined, 2);
    if (newText !== repoStr) {
      setRepoStr(JSON.stringify(tree, undefined, 2));
    }
  }, [tree]);

  const allToggleablePaths = useMemo(() => {
    const paths: string[] = [];
    const collect = (obj: any, prefix: string[]) => {
      if (typeof obj !== 'object' || obj === null) return;
      for (const key of Object.keys(obj)) {
        const tp = [...prefix, key];
        if (canToggleHidden(tp)) {
          paths.push(tp.join('/'));
        }
        collect(obj[key], tp);
      }
    };
    collect(tree, []);
    return paths;
  }, [tree, canToggleHidden]);

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
                setHiddenPaths={setHiddenPaths}
                allToggleablePaths={allToggleablePaths}
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
                nodeActions={nodeActions}
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
