import { BASEURL } from '../shared/baseurl';
import { ExecRequest, ExecResponse } from '../shared/types';
import { pushRunning, pushError, dismiss } from './snackbar';

export type PythonCommand = {
  results: string[];
  code: string;
};

const SHORT_ERROR_MAX = 120;

const shortenError = (error: string) => {
  const firstLine = error.split('\n')[0];
  return firstLine.length > SHORT_ERROR_MAX
    ? firstLine.slice(0, SHORT_ERROR_MAX - 1) + '…'
    : firstLine;
};

const getCallerLabel = (): string => {
  const stack = new Error().stack ?? '';
  const lines = stack.split('\n');
  const caller = lines.find((line, i) => i > 0 && !line.includes('/fetchPython') && !line.includes('/snackbar'));
  if (!caller) return 'Python';
  const match = caller.match(/at\s+(?:async\s+)?([^\s(]+)/);
  const name = match?.[1];
  if (!name || name.startsWith('http') || name === '<anonymous>') return 'Python';
  return name;
};

const fetchPythonDirect = async (code: string): Promise<ExecResponse> => {
  console.log('executing', code);
  const payload: ExecRequest = { code };
  let text: string | undefined;
  try {
    await new Promise(r => setTimeout(r, 1500));
    const r = await fetch(`${BASEURL}/exec`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    text = await r.text();
    const parsed = JSON.parse(text) as ExecResponse;
    if (parsed.problem) {
      console.error('python error:', parsed.problem.error, parsed.problem.traceback);
    } else {
      console.log('result =', parsed.data);
    }
    return parsed;
  } catch (e: any) {
    const error = text !== undefined
      ? `unexpected response: ${text.slice(0, 500)}`
      : `network error: ${e?.message ?? e}`;
    console.error(error);
    return {
      data: null,
      problem: { error, traceback: '' },
    };
  }
};

export const fetchPython = async (...commands: PythonCommand[]): Promise<{ data: any }> => {
  const caller = getCallerLabel();
  const lines: string[] = [];
  const resultVars: string[] = [];

  for (const cmd of commands) {
    lines.push(cmd.code);
    resultVars.push(...cmd.results);
  }

  lines.push('result = {}');
  for (const name of resultVars) {
    lines.push(`result["${name}"] = ${name}`);
  }

  const key = pushRunning(caller);
  try {
    const response = await fetchPythonDirect(lines.join('\n'));
    if (response.problem) {
      pushError(shortenError(response.problem.error));
      return { data: null };
    }
    return { data: response.data };
  } finally {
    dismiss(key);
  }
};
