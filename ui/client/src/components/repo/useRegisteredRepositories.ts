import { useEffect, useState } from "react";
import { fetchRegisteredRepositories } from "../../io/fetchRegisteredRepositories";
import { Repository } from "../../shared/types";

export const useRegisteredRepositories = () => {
  const [repositories, setRepositories] = useState<Repository[]>([]);

  const fetchRepositories = async () => {
    setRepositories(await fetchRegisteredRepositories());
  };

  useEffect(() => {
    fetchRepositories();
  }, []);

  return { repositories, fetchRepositories };
}
