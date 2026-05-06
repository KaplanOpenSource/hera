import { Table, TableBody, TableCell, TableHead, TableRow } from "@mui/material";

const lookupValue = (
  doc: { [key: string]: any },
  overridePath: string,
) => {
  const parts = overridePath.split('/');
  const isExperiment = parts[0].toLowerCase() === 'experiment';
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

  return (
    <Table size="small">
      <TableHead>
        <TableRow>
          <TableCell sx={{ color: 'inherit' }}>File</TableCell>
          <TableCell sx={{ color: 'inherit' }}>Value</TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {rows.map(row => (
          <TableRow key={row.path}>
            <TableCell sx={{ color: 'inherit' }}>{row.path}</TableCell>
            <TableCell sx={{ color: 'inherit' }}>
              {typeof row.value === 'object'
                ? JSON.stringify(row.value)
                : String(row.value)}
            </TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
};
