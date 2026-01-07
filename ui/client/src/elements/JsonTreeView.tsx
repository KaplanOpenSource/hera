import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { IconButton, Stack, Typography, Box } from "@mui/material";
import DeleteIcon from "@mui/icons-material/Delete";

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
}: {
  label: string;
  value: JsonValue;
  setData: (next: JsonValue | undefined) => void;
}) => {
  const onDelete = (e: React.MouseEvent) => {
    e.stopPropagation();
    setData(undefined);
  };

  const labelNode = (
    <Stack direction="row" spacing={1} alignItems="center">
      <Typography variant="body2">{label}</Typography>
      <IconButton size="small" onClick={onDelete}>
        <DeleteIcon fontSize="inherit" />
      </IconButton>
    </Stack>
  );

  if (Array.isArray(value)) {
    return (
      <TreeItem itemId={label} label={labelNode}>
        {value.map((item, index) => (
          <JsonTreeNode
            key={index}
            label={`[${index}]`}
            value={item}
            setData={(next) =>
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
      <TreeItem itemId={label} label={labelNode}>
        {Object.entries(value).map(([key, val]) => (
          <JsonTreeNode
            key={key}
            label={key}
            value={val}
            setData={(next) =>
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
      itemId={label}
      label={
        <Stack direction="row" spacing={1} alignItems="center">
          <Typography variant="body2">
            {label}: {String(value)}
          </Typography>
          <IconButton size="small" onClick={onDelete}>
            <DeleteIcon fontSize="inherit" />
          </IconButton>
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
          value={data}
          setData={(next) => {
            if (next !== undefined) setData(next);
          }}
        />
      </SimpleTreeView>
    </Box>
  );
};
