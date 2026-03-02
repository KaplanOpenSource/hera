import { execPython } from './execPython';

export type PythonCommand = {
  results: string[];
  code: string;
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

  const { data, problem } = await execPython(lines.join('\n'));
  return { data: problem ? undefined : data, problem };
};
