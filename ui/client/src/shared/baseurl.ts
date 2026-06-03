const apiHost = import.meta.env.VITE_API_HOST || 'localhost';
const apiPort = import.meta.env.VITE_API_PORT || '8000';
// Production: the server serves the bundle, so the API is same-origin (works on any host).
// Dev: Vite serves on a different port, so target the API host:port explicitly.
export const BASEURL = import.meta.env.DEV ? `http://${apiHost}:${apiPort}` : '';
