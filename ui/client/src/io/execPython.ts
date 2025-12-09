import { BASEURL } from "../shared/baseurl";
import { ExecRequest } from "../shared/types";

export const execPython = async (code: string): Promise<{ data: any; problem: undefined | string; }> => {
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
    console.log('got', data);
    return { data, problem: undefined };
  } catch (e: any) {
    const problem = e?.message ?? 'Failed to run';
    console.log('problem:', problem);
    return { data: undefined, problem };
  }
}
