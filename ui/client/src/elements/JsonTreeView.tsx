import { Box } from "@mui/material";
import { SimpleTreeView } from "@mui/x-tree-view";
import { JsonTreeNode, JsonValue } from "./JsonTreeNode";

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
