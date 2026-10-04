import { create } from 'zustand';
import { fetchPythonClean } from '../io/fetchPython';

type ServerUserStore = {
  username: string | null;
  inDocker: boolean;
  loadServerUser: () => Promise<void>;
};

let loadStarted = false;

// The user that started the server, read once and shared by every view that
// shows it (the project chooser label, the user indicator).
export const useServerUserStore = create<ServerUserStore>()((set) => ({
  username: null,
  inDocker: false,
  loadServerUser: async () => {
    if (loadStarted) return;
    loadStarted = true;
    const response = await fetchPythonClean({
      results: ['username', 'inDocker'],
      code: "import getpass, os; username = getpass.getuser(); inDocker = os.path.exists('/.dockerenv')",
    });
    set({
      username: (response.data?.username as string) ?? null,
      inDocker: Boolean(response.data?.inDocker),
    });
  },
}));
