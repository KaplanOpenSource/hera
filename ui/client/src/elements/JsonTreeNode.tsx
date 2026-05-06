import DeleteIcon from "@mui/icons-material/Delete";
import { IconButton, Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ReactNode } from "react";

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
  overridePath,
}: {
  label: string;
  value: JsonValue;
  setData?: (next: JsonValue | undefined) => void;
  parentKey: string,
  rootLabelComponents?: ReactNode | null,
  overrides?: Map<string, string[]>;
  overridePath?: string;
}) => {
  const isOverridden = overridePath !== undefined && overrides?.has(overridePath);
  const hasOverriddenDescendant = !isOverridden
    && overridePath !== undefined
    && overrides !== undefined
    && Array.from(overrides.keys()).some(k => k.startsWith(overridePath + '/'));

  const onDelete = (e: React.MouseEvent) => {
    e.stopPropagation();
    if (setData) {
      setData(undefined);
    }
  };

  const labelNode = (
    <Stack direction="row" spacing={1} alignItems="center">
      <Typography variant="body2" color={isOverridden ? 'error' : undefined}>{label}</Typography>
      {hasOverriddenDescendant
        ? <Typography variant="body2" color="error">●</Typography>
        : null}
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
            overridePath={overridePath !== undefined ? `${overridePath}/${key}` : undefined}
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
          <Typography variant="body2" color={isOverridden ? 'error' : undefined}>
            {label}: {String(value)}
          </Typography>
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