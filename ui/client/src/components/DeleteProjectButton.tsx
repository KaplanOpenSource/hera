import { Add, Delete } from "@mui/icons-material";
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { useConfirm } from "../elements/useConfirm";
import { useProjectStore } from "../stores/useProjectStore";
import { execPython } from "../io/execPython";

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
          await execPython(`
import json
from hera.datalayer import All
docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)['documents']
docs = sorted(docs, key=lambda d: d['type'] == '${currProjectName}__config__')
for d in docs:
    All.deleteDocumentByID(d['_id']['$oid'])
            `);
        }
      }}
    >
      <Delete />
      {ConfirmDialog}
    </ButtonTooltip>
  </>)
}