import { useEffect, useState } from 'react';
import { validateNodeParams } from '../../io/validateNodeParams';
import { WorkflowNode } from '../../shared/types';

const VALIDATE_SETTLE_MS = 500;

// Hermes's complaint about this node's parameter values, or '' when it has none.
export const useNodeParamsValidation = (node: WorkflowNode): string => {
  const [checked, setChecked] = useState({ request: '', message: '' });
  const type = node.type ?? '';
  const params = node.Execution?.input_parameters ?? {};
  // A string, so it compares by value between renders.
  const request = JSON.stringify({ type, params });

  useEffect(() => {
    const timer = setTimeout(async () => {
      if (!type) {
        setChecked({ request, message: '' });
        return;
      }
      try {
        const result = await validateNodeParams({ type, params });
        setChecked({ request, message: result.message });
      } catch (error) {
        console.error('node params validation', error);
      }
    }, VALIDATE_SETTLE_MS);
    return () => clearTimeout(timer);
  }, [request]);

  if (checked.request !== request) {
    return '';
  }
  return checked.message;
};
