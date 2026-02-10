import { Stack, Typography, Paper, Chip, Tooltip, IconButton, TextField, Button } from "@mui/material";
import { useState } from "react";
import { InjuryLevelType, makeDefaultLevelParams, LevelParamsEditor } from "./LevelParamsEditor";
import { Add, Delete } from "@mui/icons-material";

export const InjuryLevelsEditor = ({
  levelType,
  levels,
  parameters,
  onChange,
}: {
  levelType: InjuryLevelType;
  levels: string[];
  parameters: { [key: string]: any };
  onChange: (levels: string[], parameters: { [key: string]: any }) => void;
}) => {
  const [newLevelName, setNewLevelName] = useState("");

  const addLevel = () => {
    const name = newLevelName.trim();
    if (!name || levels.includes(name)) return;
    onChange([...levels, name], { ...parameters, [name]: makeDefaultLevelParams(levelType) });
    setNewLevelName("");
  };

  const removeLevel = (name: string) => {
    const { [name]: _, ...rest } = parameters;
    onChange(
      levels.filter((l) => l !== name),
      rest
    );
  };

  const reorderLevel = (index: number, direction: -1 | 1) => {
    const target = index + direction;
    if (target < 0 || target >= levels.length) return;
    const next = [...levels];
    [next[index], next[target]] = [next[target], next[index]];
    onChange(next, parameters);
  };

  const updateLevelParams = (name: string, params: any) => {
    onChange(levels, { ...parameters, [name]: params });
  };

  return (
    <Stack spacing={1}>
      <Typography variant="caption" color="text.secondary">
        Levels are ordered highest severity first
      </Typography>

      {levels.map((name, i) => (
        <Paper key={name} variant="outlined" sx={{ p: 1.5 }}>
          <Stack spacing={1.5}>
            <Stack direction="row" alignItems="center" justifyContent="space-between">
              <Stack direction="row" alignItems="center" spacing={1}>
                <Chip label={name} size="small" color="primary" variant="outlined" />
                <Typography variant="caption" color="text.secondary">
                  #{i + 1}
                </Typography>
              </Stack>
              <Stack direction="row" spacing={0.5}>
                <Tooltip title="Move up (higher severity)">
                  <span>
                    <IconButton size="small" disabled={i === 0} onClick={() => reorderLevel(i, -1)}>
                      ▲
                    </IconButton>
                  </span>
                </Tooltip>
                <Tooltip title="Move down (lower severity)">
                  <span>
                    <IconButton
                      size="small"
                      disabled={i === levels.length - 1}
                      onClick={() => reorderLevel(i, 1)}
                    >
                      ▼
                    </IconButton>
                  </span>
                </Tooltip>
                <Tooltip title="Remove level">
                  <IconButton
                    size="small"
                    onClick={() => removeLevel(name)}
                  >
                    <Delete fontSize="small" />
                  </IconButton>
                </Tooltip>
              </Stack>
            </Stack>
            <LevelParamsEditor
              levelType={levelType}
              params={parameters[name] ?? makeDefaultLevelParams(levelType)}
              onChange={(p) => updateLevelParams(name, p)}
            />
          </Stack>
        </Paper>
      ))}

      <Stack direction="row" spacing={1} alignItems="center">
        <TextField
          label="New level name"
          size="small"
          value={newLevelName}
          onChange={(e) => setNewLevelName(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && addLevel()}
          sx={{ flex: 1 }}
        />
        <Button
          variant="outlined"
          size="small"
          startIcon={<Add />}
          onClick={addLevel}
          disabled={!newLevelName.trim() || levels.includes(newLevelName.trim())}
        >
          Add
        </Button>
      </Stack>
    </Stack>
  );
};
