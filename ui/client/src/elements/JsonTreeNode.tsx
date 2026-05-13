import DeleteIcon from "@mui/icons-material/Delete";
import { IconButton, Stack, Tooltip, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ReactNode } from "react";
import { OverrideLabel } from "./OverrideLabel";

export type JsonValue =
  | string
  | number
  | boolean
  | null
  | JsonValue[]
  | { [key: string]: JsonValue };

export const JsonTreeNode = ({
  label,
  value,
  setData,
  parentKey,
  rootLabelComponents = null,
  overrides,
  treePath = [],
  repoJsons,
  hiddenPaths,
  showHidden,
  nodeActions,
}: {
  label: string;
  value: JsonValue;
  setData?: (next: JsonValue | undefined) => void;
  parentKey: string,
  rootLabelComponents?: ReactNode | null,
  overrides?: { [path: string]: string[] };
  treePath?: string[];
  repoJsons?: { [path: string]: { [key: string]: any } };
  hiddenPaths?: Set<string>;
  showHidden?: boolean;
  nodeActions?: (treePath: string[]) => ReactNode;
}) => {
  const treePathKey = treePath.join('/');
  const isHidden = hiddenPaths?.has(treePathKey) ?? false;

  if (hiddenPaths && !showHidden && isHidden) {
    return null;
  }

  const isOverridden = treePathKey in (overrides ?? {});
  const hasOverriddenDescendant = !isOverridden
    && overrides !== undefined
    && Object.keys(overrides).some(k => k.startsWith(treePathKey + '/'));

  const onDelete = (e: React.MouseEvent) => {
    e.stopPropagation();
    if (setData) {
      setData(undefined);
    }
  };

  const labelNode = (
    <Stack direction="row" spacing={1} alignItems="center">
      <OverrideLabel
        isOverridden={isOverridden}
        treePathKey={treePathKey}
        overrides={overrides}
        repoJsons={repoJsons}
      >
        {label}
      </OverrideLabel>
      {hasOverriddenDescendant
        ? (
          <Tooltip title="Some children have overrides">
            <Typography variant="body2" color="error">●</Typography>
          </Tooltip>
        )
        : null}
      {nodeActions?.(treePath)}
      {setData === undefined ? null :
        <IconButton size="small" onClick={onDelete}>
          <DeleteIcon fontSize="inherit" />
        </IconButton>
      }
    </Stack>
  );

  if (Array.isArray(value)) {
    return (
      <TreeItem itemId={label} label={labelNode}>
        {value.map((item, index) => (
          <JsonTreeNode
            key={index}
            label={`[${index}]`}
            parentKey={`${parentKey}[${index}]`}
            value={item}
            overrides={overrides}
            treePath={[...treePath, String(index)]}
            repoJsons={repoJsons}
            hiddenPaths={hiddenPaths}
            showHidden={showHidden}
            nodeActions={nodeActions}
            setData={setData === undefined ? undefined :
              (next) =>
                setData(
                  next === undefined
                    ? value.filter((_, i) => i !== index)
                    : value.map((v, i) => (i === index ? next : v))
                )
            }
          />
        ))}
      </TreeItem>
    );
  }

  if (typeof value === "object" && value !== null) {
    return (
      <TreeItem
        itemId={`${parentKey}.${label}`}
        label={(<Stack direction={'row'}>
          {labelNode}
          {rootLabelComponents}
        </Stack>)}
      >
        {Object.entries(value).map(([key, val]) => (
          <JsonTreeNode
            key={key}
            label={key}
            parentKey={`${parentKey}.${label}[${key}]`}
            value={val}
            setData={setData === undefined ? undefined :
              (next) =>
                setData(
                  next === undefined
                    ? Object.fromEntries(
                      Object.entries(value).filter(([k]) => k !== key)
                    )
                    : { ...value, [key]: next }
                )
            }
            overrides={overrides}
            treePath={[...treePath, key]}
            repoJsons={repoJsons}
            hiddenPaths={hiddenPaths}
            showHidden={showHidden}
            nodeActions={nodeActions}
          />
        ))}
      </TreeItem>
    );
  }

  return (
    <TreeItem
      itemId={`${parentKey}.${label}`}
      label={
        <Stack direction="row" spacing={1} alignItems="center">
          <OverrideLabel
            isOverridden={isOverridden}
            treePathKey={treePathKey}
            overrides={overrides}
            repoJsons={repoJsons}
          >
            {label}: {String(value)}
          </OverrideLabel>
          {nodeActions?.(treePath)}
          {setData === undefined ? null :
            <IconButton size="small" onClick={onDelete}>
              <DeleteIcon fontSize="inherit" />
            </IconButton>
          }
        </Stack>
      }
    />
  );
};