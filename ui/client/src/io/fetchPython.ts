import { BASEURL } from '../shared/baseurl';
import { ExecRequest, ExecResponse } from '../shared/types';
import { pushRunning, pushError, dismiss } from './snackbar';

export type PythonCommand = {
  results: string[];
  code: string;
  label?: string;
};

const SHORT_ERROR_MAX = 120;

const shortenError = (error: string) => {
  const firstLine = error.split('\n')[0];
  return firstLine.length > SHORT_ERROR_MAX
    ? firstLine.slice(0, SHORT_ERROR_MAX - 1) + '…'
    : firstLine;
};

const fetchPythonDirect = async (code: string): Promise<ExecResponse> => {
  console.log('executing', code);
  const payload: ExecRequest = { code };
  let text: string | undefined;
  try {
    // await new Promise(r => setTimeout(r, 1500));
    const r = await fetch(`${BASEURL}/exec`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    text = await r.text();
    const parsed = JSON.parse(text) as ExecResponse;
    if (parsed.problem) {
      console.error(
        'python error:', parsed.problem.error,
        '\ntrace:\n', parsed.problem.traceback,
        '\npayload:\n', code);
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

const assembleCode = (commands: PythonCommand[]) => {
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

  return lines.join('\n');
};

export const fetchPythonClean = async (...commands: PythonCommand[]): Promise<ExecResponse> => {
  return fetchPythonDirect(assembleCode(commands));
};

export const fetchPython = async (...commands: PythonCommand[]): Promise<{ data: any }> => {
  const labels = commands.map(c => c.label).filter(Boolean);
  const label = labels.length > 0 ? labels.join(', ') : 'Python';
  const key = pushRunning(label);
  try {
    const response = await fetchPythonClean(...commands);
    if (response.problem) {
      pushError(`${label}: ${shortenError(response.problem.error)}`);
      return { data: null };
    }
    return { data: response.data };
  } finally {
    dismiss(key);
  }
};
