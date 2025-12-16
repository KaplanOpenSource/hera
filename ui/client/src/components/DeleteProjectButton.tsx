import { Add, Delete } from "@mui/icons-material";
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { useConfirm } from "../elements/useConfirm";
import { useProjectStore } from "../stores/useProjectStore";

export const DeleteProjectButton = ({ }) => {
  const { confirmOpen, ConfirmDialog } = useConfirm()
  const { currProjectName } = useProjectStore();

  return (<>
    <ButtonTooltip
      title='Delete project'
      onClick={async () => {
        const { confirmed, text } = await confirmOpen({
          title: `Are you sure you want to delete ${currProjectName}?`,
          requireText: true,
          textLabel: 'Type project name to confirm',
          textPlaceholder: currProjectName,
          textValidate: (text) => text === currProjectName,
        });
        if (confirmed && text === currProjectName) {
          // console.l
        }
      }}
    >
      <Delete />
      {ConfirmDialog}
    </ButtonTooltip>
  </>)
}