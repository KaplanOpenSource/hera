import { installReporter } from '../api/api';
import { Reporter } from '../api/Reporter';
import { dismiss, pushError, pushRunning } from './snackbar';

// The api package knows nothing about notistack; this joins the two.
const snackbarReporter: Reporter = {
  started: (label: string) => {
    return pushRunning(label);
  },
  finished: (ticket) => {
    dismiss(ticket);
  },
  failed: (message: string) => {
    pushError(message);
  },
};

export const installSnackbarReporter = () => {
  installReporter(snackbarReporter);
};
