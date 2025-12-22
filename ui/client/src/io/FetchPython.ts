import { useEffect } from "react";
import { execPython } from "./execPython";

export const FetchPython = ({
  code,
  onSuccess,
  onFail,
}: {
  code: string,
  onSuccess: (data: any) => void,
  onFail?: (problem: string) => void,
}) => {
  useEffect(() => {
    (async () => {
      const { data, problem } = await execPython(code);
      if (!problem) {
        onSuccess(data);
      } else {
        if (onFail) {
          onFail(problem);
        }
      }
    })();
  }, [code]);

  return null;
}

