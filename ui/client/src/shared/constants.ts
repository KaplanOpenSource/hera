const localhost = ['localhost', '127.0.0.1'];
const isLocal = localhost.includes(window.location.hostname) && window.location.port === '8000';
export const API_BASE = isLocal ? '' : 'http://localhost:8000';
