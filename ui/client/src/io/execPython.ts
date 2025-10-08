import { API_BASE } from "../shared/constants";
import { ExecRequest } from "../shared/types";

export const execPython = async (code: string): Promise<{ data: any; problem: undefined | string; }> => {
  try {
    const payload: ExecRequest = {
      code,
    };
    const r = await fetch(`${API_BASE}/exec`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    const data = await r.json();
    return { data, problem: undefined };
  } catch (e: any) {
    const problem = e?.message ?? 'Failed to run';
    return { data: undefined, problem };
  }
}
