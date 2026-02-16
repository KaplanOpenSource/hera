import { Stack, TextField, Typography } from "@mui/material";
import { ProjectDocument } from "../../shared/types";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { Edit } from "@mui/icons-material";

export const VersionFields = ({
  projectDoc,
  setProjectDoc,
}: {
  projectDoc: ProjectDocument,
  setProjectDoc: (changedDoc: ProjectDocument) => void,
}) => {
  const [editable, setEditable] = useState(false);
  const rawVersion = projectDoc.desc.version;
  const isString = typeof rawVersion === 'string';
  const version = Array.isArray(rawVersion)
    ? rawVersion
    : isString
      ? rawVersion.split('.').map(Number)
      : [];

  return (
    <Stack direction={'row'} alignItems={'center'}>
      {editable
        ? version.map((vers, ivers) => (
            <TextField key={ivers}
              size='small'
              sx={{
                width: 75, fontSize: 6, height: '20px',
                '& .MuiInputBase-input': {
                  paddingRight: '4px',
                },
              }}
              slotProps={{
                input: {
                  sx: { fontSize: 10, padding: '2px -50px', height: 24, margin: 0 },
                },
              }}
              type='number'
              value={vers}
              onChange={e => {
                e.stopPropagation();
                const newVersion = [...version];
                const parsed = parseFloat(e.target.value);
                newVersion[ivers] = isNaN(parsed) ? 0 : parsed;
                const versionValue = isString ? newVersion.join('.') : newVersion;
                setProjectDoc({ ...projectDoc, desc: { ...projectDoc.desc, version: versionValue } });
              }}
            />
          ))
        : (
          <Typography sx={{ fontSize: 12 }}>
            {version.join('.')}
          </Typography>
        )
      }
      <ButtonTooltip onClick={() => setEditable(!editable)}>
        <Edit fontSize="small" />
      </ButtonTooltip>
    </Stack>
  )
}