import { BASEURL } from '../shared/baseurl';
import { ExecRequest, ExecResponse } from '../shared/types';

export type PythonCommand = {
  results: string[];
  code: string;
};

const fetchPythonDirect = async (code: string): Promise<ExecResponse> => {
  console.log('executing', code);
  const payload: ExecRequest = { code };
  let text: string | undefined;
  try {
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

  const response = await fetchPythonDirect(lines.join('\n'));
  return { data: response.problem ? null : response.data };
};
