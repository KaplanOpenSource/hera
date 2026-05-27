const localhost = ['localhost', '127.0.0.1'];
const apiPort = import.meta.env.VITE_API_PORT || '8000';
const isLocal = localhost.includes(window.location.hostname) && !import.meta.env.DEV;
export const BASEURL = isLocal ? '' : `http://localhost:${apiPort}`;
