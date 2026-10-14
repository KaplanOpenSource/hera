import { SnackbarKey } from 'notistack';
import { BASEURL } from '../shared/baseurl';
import { WARMING_UP } from './fetchPython';
import { dismiss, pushError, pushRunning } from './snackbar';

const SHORT_ERROR_MAX = 120;

export type EntrypointOptions = {
  // Shown in the running toast, and in front of the error if the call fails.
  label: string,
  // No toasts at all, for background calls like the periodic reload.
  silent?: boolean,
};

export type QueryParams = { [key: string]: string };

const shorten = (text: string) => {
  if (text.length <= SHORT_ERROR_MAX) {
    return text;
  }
  return text.slice(0, SHORT_ERROR_MAX - 1) + '…';
};

const report = (options: EntrypointOptions, problem: string) => {
  console.error(`${options.label}: ${problem}`);
  // Stay quiet while hera warms up; the startup gate handles that wait.
  if (!options.silent && problem !== WARMING_UP) {
    pushError(`${options.label}: ${shorten(problem)}`);
  }
};

// One call to a named server entrypoint. Returns its JSON answer, or null when
// the call failed.
const call = async (path: string, init: RequestInit, options: EntrypointOptions): Promise<any | null> => {
  let key: SnackbarKey | null = null;
  if (!options.silent) {
    key = pushRunning(options.label);
  }
  try {
    // no-store: the answer must never come from the browser cache.
    const response = await fetch(`${BASEURL}${path}`, { cache: 'no-store', ...init });
    const contentType = response.headers.get('content-type');
    // A server without this route answers with the SPA's HTML, not an error.
    if (!contentType?.includes('application/json')) {
      report(options, `${response.status} ${contentType}`);
      return null;
    }
    const body = await response.json();
    if (!response.ok) {
      report(options, String(body?.detail ?? body?.error ?? response.status));
      return null;
    }
    return body;
  } catch (e: any) {
    // The server restarting must not break the UI, so a dead connection is
    // reported like any other failure instead of throwing at the caller.
    report(options, `network error: ${e?.message ?? e}`);
    return null;
  } finally {
    if (key !== null) {
      dismiss(key);
    }
  }
};

export const getJson = async (path: string, params: QueryParams, options: EntrypointOptions) => {
  const query = new URLSearchParams(params).toString();
  let suffix = '';
  if (query) {
    suffix = `?${query}`;
  }
  return await call(`${path}${suffix}`, {}, options);
};

export const postJson = async (path: string, body: object, options: EntrypointOptions) => {
  return await call(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  }, options);
};
