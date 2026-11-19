import { Upload } from "@mui/icons-material";
import { ChangeEvent, useRef } from "react";
import { useNavigate } from "react-router-dom";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { projectApi } from "../../api/api";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { ProjectName } from "../../shared/types";

export const LoadProjectButton = () => {
  const inputRef = useRef<HTMLInputElement>(null);
  const navigate = useNavigate();

  const loadProject = async (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    // Cleared so picking the same file again still fires a change.
    e.target.value = '';
    if (!file) {
      return;
    }
    const data = await projectApi.load(await file.arrayBuffer(), file.name);
    if (!data) {
      return;
    }
    // Same as after create: the project travels with the navigation. See AddProjectButton.
    useProjectListStore.getState().setProjectNames(data.projects as ProjectName[]);
    navigate('/' + encodeURIComponent(data.project.name), { state: { project: data.project } });
  }

  return (<>
    <ButtonTooltip
      title='Load project'
      onClick={() => inputRef.current?.click()}
    >
      <Upload />
    </ButtonTooltip>
    <input
      ref={inputRef}
      type="file"
      accept=".zip"
      hidden
      data-testid="load-project-input"
      onChange={loadProject}
    />
  </>);
}
