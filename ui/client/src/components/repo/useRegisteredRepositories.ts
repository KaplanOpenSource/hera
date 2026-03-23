import { useEffect, useState } from "react";
import { fetchPython } from "../../io/fetchPython";
import { Repository } from "../../shared/types";

export const useRegisteredRepositories = () => {
  const [repositories, setRepositories] = useState<Repository[]>([]);

  const fetchRepositories = async () => {
    const { data, problem } = await fetchPython({
      results: ['repos'],
      code: `
from hera.utils.data.toolkit import dataToolkit
repos = dataToolkit().getRepositoryTable().to_dict(orient='records')
`,
    });
    if (!problem && data?.repos) {
      setRepositories(data.repos);
    }
  };

  useEffect(() => {
    fetchRepositories();
  }, []);

  return { repositories, fetchRepositories };
}
