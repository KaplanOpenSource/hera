// A ticket for one running call, handed back so the call can be dismissed.
export type Ticket = string | number;

// How the api package tells the user a call is running or has failed. The
// client installs one at startup; without it the calls simply stay quiet.
export type Reporter = {
  started: (label: string) => Ticket,
  finished: (ticket: Ticket) => void,
  failed: (message: string) => void,
};
