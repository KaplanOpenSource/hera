import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { IconButton, Stack, Typography, Box } from "@mui/material";
import DeleteIcon from "@mui/icons-material/Delete";
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
}: {
  label: string;
  value: JsonValue;
  setData?: (next: JsonValue | undefined) => void;
  parentKey: string,
  rootLabelComponents?: ReactNode | null,
}) => {
  const onDelete = (e: React.MouseEvent) => {
    e.stopPropagation();
    if (setData) {
      setData(undefined);
    }
  };

  const labelNode = (
    <Stack direction="row" spacing={1} alignItems="center">
      <Typography variant="body2">{label}</Typography>
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
          <Typography variant="body2">
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

export const JsonTreeView = ({
  data,
  setData,
}: {
  data: JsonValue;
  setData: (next: JsonValue) => void;
}) => {
  return (
    <Box sx={{ width: "100%", height: "100%" }}>
      <SimpleTreeView
        onClick={e => e.stopPropagation()}
      >
        <JsonTreeNode
          label="root"
          parentKey="root"
          value={data}
          setData={(next) => {
            if (next !== undefined) setData(next);
          }}
        />
      </SimpleTreeView>
    </Box>
  );
};
