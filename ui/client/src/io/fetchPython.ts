import { BASEURL } from '../shared/baseurl';
import { ExecRequest } from '../shared/types';

export type PythonCommand = {
  results: string[];
  code: string;
};

const execRaw = async (code: string): Promise<{ data: any; problem: undefined | string }> => {
  try {
    console.log('executing', code);
    const payload: ExecRequest = {
      code,
    };
    const r = await fetch(`${BASEURL}/exec`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    const data = await r.json();
    console.log('result =', data);
    return { data, problem: undefined };
  } catch (e: any) {
    const problem = e?.message ?? 'Failed to run';
    console.trace('problem:', problem);
    return { data: undefined, problem };
  }
};

export const fetchPython = async (...commands: PythonCommand[]): Promise<{ data: any; problem: string | undefined }> => {
  const lines: string[] = [];
  const resultVars: string[] = [];

  for (const cmd of commands) {
    lines.push(cmd.code);
    resultVars.push(...cmd.results);
  }

  if (resultVars.length > 0) {
    lines.push('result = {}');
    for (const name of resultVars) {
      lines.push(`result["${name}"] = ${name}`);
    }
  }

  const { data, problem } = await execRaw(lines.join('\n'));
  return { data: problem ? undefined : data, problem };
};
