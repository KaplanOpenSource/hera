import { useEffect, useState } from "react";
import { repositoryApi } from "../../api/api";
import { Repository } from "../../shared/types";

export const useRegisteredRepositories = () => {
  const [repositories, setRepositories] = useState<Repository[]>([]);

  const fetchRepositories = async () => {
    setRepositories(await repositoryApi.registered());
  };

  useEffect(() => {
    fetchRepositories();
  }, []);

  return { repositories, fetchRepositories };
}
