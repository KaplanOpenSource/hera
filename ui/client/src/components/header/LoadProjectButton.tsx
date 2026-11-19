import { Upload } from "@mui/icons-material";
import { ChangeEvent, useRef, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { LoadProjectDialog } from "./LoadProjectDialog";

// Opens the file picker first; the dialog shows once a zip is picked.
export const LoadProjectButton = ({
  onLoaded,
}: {
  onLoaded?: () => void,
}) => {
  const [file, setFile] = useState<File | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  const pickFile = (e: ChangeEvent<HTMLInputElement>) => {
    const picked = e.target.files?.[0];
    // Cleared so picking the same file again still fires a change.
    e.target.value = '';
    if (!picked) {
      return;
    }
    setFile(picked);
  }

  return (<>
    <ButtonTooltip
      button
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
      onChange={pickFile}
    />
    {file && (
      <LoadProjectDialog
        file={file}
        onPickAnother={() => inputRef.current?.click()}
        onClose={() => setFile(null)}
        onLoaded={onLoaded}
      />
    )}
  </>);
}
