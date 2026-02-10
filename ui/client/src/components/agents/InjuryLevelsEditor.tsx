import { Stack, Typography, Paper, Tooltip, IconButton, TextField, Button } from "@mui/material";
import { useState } from "react";
import { InjuryLevelType, makeDefaultLevelParams, LevelParamsEditor } from "./LevelParamsEditor";
import { Add, Delete, DragIndicator } from "@mui/icons-material";
import {
  DndContext,
  closestCenter,
  KeyboardSensor,
  PointerSensor,
  useSensor,
  useSensors,
  DragEndEvent,
} from "@dnd-kit/core";
import {
  SortableContext,
  sortableKeyboardCoordinates,
  verticalListSortingStrategy,
  useSortable,
  arrayMove,
} from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import { restrictToVerticalAxis, restrictToParentElement } from "@dnd-kit/modifiers";
import { DetailsViewItemName } from "../details/DetailsViewItemName";

const SortableLevelItem = ({
  name,
  index,
  levelType,
  parameters,
  removeLevel,
  updateLevelParams,
  renameLevel,
}: {
  name: string;
  index: number;
  levelType: InjuryLevelType;
  parameters: { [key: string]: any };
  removeLevel: (name: string) => void;
  updateLevelParams: (name: string, params: any) => void;
  renameLevel: (oldName: string, newName: string) => void;
}) => {
  const { attributes, listeners, setNodeRef, transform, transition, isDragging } = useSortable({
    id: name,
  });

  const style = {
    transform: CSS.Transform.toString(transform),
    transition,
    opacity: isDragging ? 0.5 : 1,
  };

  return (
    <Paper ref={setNodeRef} style={style} variant="outlined" sx={{ p: 1.5 }}>
      <Stack spacing={1.5}>
        <Stack direction="row" alignItems="center" justifyContent="space-between">
          <Stack direction="row" alignItems="center" spacing={1} sx={{ flex: 1, minWidth: 0 }}>
            <IconButton
              size="small"
              sx={{ cursor: "grab", touchAction: "none" }}
              {...attributes}
              {...listeners}
            >
              <DragIndicator fontSize="small" />
            </IconButton>
            <DetailsViewItemName
              itemKey={name}
              setItemKey={(newName) => {
                if (newName && newName !== name) {
                  renameLevel(name, newName);
                }
              }}
            />
            <Typography variant="caption" color="text.secondary">
              #{index + 1}
            </Typography>
          </Stack>
          <Tooltip title="Remove level">
            <IconButton size="small" onClick={() => removeLevel(name)}>
              <Delete fontSize="small" />
            </IconButton>
          </Tooltip>
        </Stack>
        <LevelParamsEditor
          levelType={levelType}
          params={parameters[name] ?? makeDefaultLevelParams(levelType)}
          onChange={(p) => updateLevelParams(name, p)}
        />
      </Stack>
    </Paper>
  );
};

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

  const sensors = useSensors(
    useSensor(PointerSensor),
    useSensor(KeyboardSensor, { coordinateGetter: sortableKeyboardCoordinates })
  );

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

  const renameLevel = (oldName: string, newName: string) => {
    const trimmed = newName.trim();
    if (!trimmed || (trimmed !== oldName && levels.includes(trimmed))) return;
    const newLevels = levels.map((l) => (l === oldName ? trimmed : l));
    const { [oldName]: params, ...rest } = parameters;
    onChange(newLevels, { ...rest, [trimmed]: params });
  };

  const handleDragEnd = (event: DragEndEvent) => {
    const { active, over } = event;
    if (over && active.id !== over.id) {
      const oldIndex = levels.indexOf(active.id as string);
      const newIndex = levels.indexOf(over.id as string);
      onChange(arrayMove(levels, oldIndex, newIndex), parameters);
    }
  };

  const updateLevelParams = (name: string, params: any) => {
    onChange(levels, { ...parameters, [name]: params });
  };

  return (
    <Stack spacing={1}>
      <Typography variant="caption" color="text.secondary">
        Levels are ordered highest severity first — drag to reorder
      </Typography>
      <DndContext
        sensors={sensors}
        collisionDetection={closestCenter}
        onDragEnd={handleDragEnd}
        modifiers={[restrictToVerticalAxis, restrictToParentElement]}
      >
        <SortableContext items={levels} strategy={verticalListSortingStrategy}>
          {levels.map((name, i) => (
            <SortableLevelItem
              key={name}
              name={name}
              index={i}
              levelType={levelType}
              parameters={parameters}
              removeLevel={removeLevel}
              updateLevelParams={updateLevelParams}
              renameLevel={renameLevel}
            />
          ))}
        </SortableContext>
      </DndContext>
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
