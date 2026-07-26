import { FindInPage } from '@mui/icons-material';
import { Stack } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchPython } from '../../io/fetchPython';
import { pushInfo } from '../../io/snackbar';
import { ProjectEntire } from '../../shared/types';
import { useProjectStore } from '../../stores/useProjectStore';
import { buildDetectNotebooksCode } from './buildDetectNotebooksCode';

// Scans the project's files directory for .ipynb files and registers any that
// are not already tracked as notebook documents.
export const DetectNotebooksButton = ({
  onDetected,
}: {
  onDetected?: () => void,
}) => {
  const filesDir = useProjectStore.getState().getProject()?.configDocument?.data.desc.filesDirectory ?? '';

  const detectNotebooks = async () => {
    const projectName = useProjectStore.getState().currProjectName;
    const { data } = await fetchPython({
      results: ['project', 'addedCount'],
      label: 'detect notebooks',
      code: buildDetectNotebooksCode({ projectName, filesDir }),
    });
    if (!data) {
      return;
    }
    useProjectStore.getState().setCurrentProject(data.project as ProjectEntire);
    const added = data.addedCount as number;
    pushInfo(added > 0 ? `Detected ${added} new notebook${added === 1 ? '' : 's'}` : 'No new notebooks found');
    onDetected?.();
  };

  return (
    <ButtonTooltip
      button
      title="Detect notebooks"
      disabled={!filesDir}
      onClick={detectNotebooks}
    >
      <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
        <FindInPage fontSize="small" />
      </Stack>
    </ButtonTooltip>
  );
};
