import { useEffect, useState } from 'react';
import { validateNodeParams } from '../../io/validateNodeParams';
import { WorkflowNode } from '../../shared/types';

const VALIDATE_SETTLE_MS = 500;

// Hermes's complaint about this node's parameter values, or '' when it has none.
export const useNodeParamsValidation = (node: WorkflowNode): string => {
  const [checked, setChecked] = useState({ request: '', message: '' });
  const request = JSON.stringify({
    type: node.type ?? '',
    params: node.Execution?.input_parameters ?? {},
  });

  useEffect(() => {
    const timer = setTimeout(async () => {
      const { type, params } = JSON.parse(request);
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
