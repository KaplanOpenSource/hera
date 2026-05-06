import CheckIcon from "@mui/icons-material/Check";
import { Paper, Table, TableBody, TableCell, TableHead, TableRow, Typography } from "@mui/material";
import { SECTION_EXPERIMENT } from "../components/details/RepoJsonMerger";

const lookupValue = (
  doc: { [key: string]: any },
  overridePath: string,
) => {
  const parts = overridePath.split('/');
  const isExperiment = parts[0].toLowerCase() === SECTION_EXPERIMENT;
  const lookupParts = isExperiment ? [parts[0], ...parts.slice(2)] : parts;
  let current: any = doc;
  for (const part of lookupParts) {
    if (typeof current !== 'object' || current === null) return undefined;
    current = current[part];
  }
  return current;
};

export const OverrideTooltipTable = ({
  overridePath,
  filePaths,
  repoJsons,
}: {
  overridePath: string;
  filePaths: string[];
  repoJsons: { [path: string]: { [key: string]: any } };
}) => {
  const rows = filePaths.map(filePath => ({
    path: filePath,
    value: lookupValue(repoJsons[filePath], overridePath),
  }));

  const chosenPath = filePaths[filePaths.length - 1];

  return (
    <Paper elevation={4} sx={{ p: 1 }}>
      <Typography variant="h6" sx={{ mb: 1, px: 1 }}>
        {overridePath}
      </Typography>
      <Table size="small">
        <TableHead>
          <TableRow>
            <TableCell />
            <TableCell>File</TableCell>
            <TableCell>Value</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {rows.map(row => (
            <TableRow key={row.path}>
              <TableCell sx={{ width: 24, p: 0 }}>
                {row.path === chosenPath ? <CheckIcon fontSize="small" color="success" /> : null}
              </TableCell>
              <TableCell>{row.path}</TableCell>
              <TableCell>
                {typeof row.value === 'object'
                  ? JSON.stringify(row.value)
                  : String(row.value)}
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Paper>
  );
};
