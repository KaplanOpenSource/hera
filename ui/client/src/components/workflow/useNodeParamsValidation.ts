import { useEffect, useState } from 'react';
import { validateNodeParams } from '../../io/validateNodeParams';
import { WorkflowNode } from '../../shared/types';

// Settle time before asking the server, so typing a value doesn't send a request
// per keystroke. Each check hits hera and MongoDB.
const VALIDATE_SETTLE_MS = 500;

// Hermes's complaint about this node's parameter values, or '' when it has none.
// One message at a time: the node's testParamValues reports only its first bad
// parameter, and names it in the text rather than in a field, so the caller shows
// this on the node.
export const useNodeParamsValidation = (node: WorkflowNode): string => {
  const [message, setMessage] = useState('');
  const type = node.type ?? '';
  const params = node.Execution?.input_parameters ?? {};
  // The values themselves are the trigger: re-check when any of them changes.
  const paramsKey = JSON.stringify(params);

  useEffect(() => {
    if (!type) {
      setMessage('');
      return;
    }
    let current = true;
    const timer = setTimeout(() => {
      validateNodeParams({ type, params })
        .then(result => {
          if (current) {
            setMessage(result.message);
          }
        })
        .catch(error => console.error('node params validation', error));
    }, VALIDATE_SETTLE_MS);
    // A newer edit landed (or the node went away): drop this answer.
    return () => {
      current = false;
      clearTimeout(timer);
    };
  }, [type, paramsKey]);

  return message;
};
