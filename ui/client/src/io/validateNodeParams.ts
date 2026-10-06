import { BASEURL } from '../shared/baseurl';

export interface NodeParamsValidation {
  ok: boolean;
  // Hermes's own message for the first bad parameter. Empty when there is
  // nothing to say (all good, unknown type, or a node with no checks).
  message: string;
}

// Tests one node's parameter values with the node's own Hermes checks. The
// server reaches hera and MongoDB for this, so call it sparingly (debounced).
export const validateNodeParams = async ({
  type,
  params,
}: {
  type: string,
  params: { [param: string]: unknown },
}): Promise<NodeParamsValidation> => {
  const response = await fetch(`${BASEURL}/node-params/validate`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ type, params }),
  });
  const text = await response.text();
  if (!response.ok) {
    const problem = JSON.parse(text);
    throw new Error(problem.error ?? text);
  }
  return JSON.parse(text);
};
