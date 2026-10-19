import { Repository } from "../shared/types";
import { getJson } from "./entrypoint";

// The repositories registered in the default project, for showing their names.
export const fetchRegisteredRepositories = async (): Promise<Repository[]> => {
  const data = await getJson('/repository/registered', {}, { label: 'registered repositories' });
  if (!data) {
    return [];
  }
  return (data.repositories ?? []) as Repository[];
};
