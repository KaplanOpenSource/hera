import { create } from 'zustand';

type ServerInfo = {
  heraVersion?: string;
  username?: string;
  inDocker?: boolean;
};

type ServerUserStore = {
  username: string | null;
  inDocker: boolean;
  heraVersion: string | null;
  setServerInfo: (info: ServerInfo) => void;
};

// About the server itself: the user that started it, whether it runs in docker,
// and the hera version. Read once at startup (see fetchInit).
export const useServerUserStore = create<ServerUserStore>()((set) => ({
  username: null,
  inDocker: false,
  heraVersion: null,
  setServerInfo: (info) => {
    set({
      username: info.username ?? null,
      inDocker: Boolean(info.inDocker),
      heraVersion: info.heraVersion ?? null,
    })
  },
}));
