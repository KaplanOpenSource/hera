import { Repository } from "../../shared/types";

export const RepositoriesInProject = ({
  repositories,
}: {
  repositories: Repository[],
}) => {
  if (repositories.length === 0) {
    return null;
  }
  return (
    <ul style={{ margin: '8px 0', paddingLeft: 24 }}>
      {repositories.map(r => (
        <li key={r.datasourceName}>
          <strong>{r.datasourceName}</strong> — {r.resource}
        </li>
      ))}
    </ul>
  )
}
