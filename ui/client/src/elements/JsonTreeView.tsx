import { useMemo } from "react";
import { RichTreeView, TreeViewBaseItem } from "@mui/x-tree-view";
import { Box, Typography } from "@mui/material";

export type JsonValue =
  | string
  | number
  | boolean
  | null
  | JsonValue[]
  | { [key: string]: JsonValue };

function buildTree(
  value: JsonValue,
  path: string = "root"
): TreeViewBaseItem[] {
  if (Array.isArray(value)) {
    return value.map((item, index) => {
      const id = `${path}[${index}]`;

      if (typeof item === "object" && item !== null) {
        return {
          id,
          label: `[${index}]`,
          children: buildTree(item, id),
        };
      }

      return {
        id,
        label: `[${index}]: ${String(item)}`,
      };
    });
  }

  if (typeof value === "object" && value !== null) {
    return Object.entries(value).map(([key, val]) => {
      const id = `${path}.${key}`;

      if (typeof val === "object" && val !== null) {
        return {
          id,
          label: key,
          children: buildTree(val, id),
        };
      }

      return {
        id,
        label: `${key}: ${String(val)}`,
      };
    });
  }

  return [
    {
      id: path,
      label: String(value),
    },
  ];
}

export const JsonTreeView = ({ data }: { data: JsonValue }) => {
  const items = useMemo<TreeViewBaseItem[]>(
    () => [
      {
        id: "root",
        label: "root",
        children: buildTree(data),
      },
    ],
    [data]
  );

  return (
    <Box sx={{ width: "100%", height: "100%" }}>
      {!data || Object.keys(data).length === 0
        ? (
          <Typography>
            Json tree is empty
          </Typography>
        ) : (
          <RichTreeView
            items={items}
            onClick={e => e.stopPropagation()}
          />
        )}
    </Box>
  );
};
