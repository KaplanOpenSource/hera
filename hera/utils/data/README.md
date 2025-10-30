# Hera Repository Attacher CLI

A tiny command-line tool that **attaches a repository/toolkit to every user** in your Hera MongoDB.
It writes/updates a single record in a `repositories` collection and appends the repo id to each user document
(via `$addToSet` on the field `repos`).

> **Scope**: This tool targets one or more MongoDB databases that you control. It does **not** require (or attempt) to connect to any external/remote multi-tenant MongoDB.

---

## Features

- Idempotent writes (`$addToSet` for users, upsert for repository record).
- Works across **multiple DBs** (e.g., `hera_*`) or an **explicit DB list**.
- Supports a **dry-run** mode.
- Simple wrapper **CLI** (`hera_repo_cli.py`) built on top of the lower-level script `add_repo_to_all_users.py` (which you keep in your repo).

---

## Prerequisites

- Python 3.8+
- `pymongo` and `click` installed in your environment:
  ```bash
  pip install pymongo click
  ```
- MongoDB reachable from your machine.
- A MongoDB user with sufficient privileges to read/write to your target DB(s).

> **Example roles**: for DB `dbhera`, a user with `readWrite` there; if you also write into DB `hera` you can grant `readWrite` on `hera` as well.

---

## Files in this folder

- `hera_repo_cli.py` — the user-facing CLI wrapper (invokes `add_repo_to_all_users.py` under the hood).
- `hera_repo_guide.ipynb` — a runnable step-by-step notebook (copy into your project or open locally).
- `README.md` — this file (copy/adapt into your repository).

> The script `add_repo_to_all_users.py` is expected to live in your repo (e.g., `hera/hera/utils/data/add_repo_to_all_users.py`).
> The CLI calls it as a subprocess.

---

## Quick Start

1. **Place the files** in your project (any folder is fine). Make sure `add_repo_to_all_users.py` is available in the same folder or adjust the CLI path accordingly.
2. **Install deps**:
   ```bash
   pip install pymongo click
   ```
3. **Add a repository to all users** (scan all DBs with prefix `hera`):
   ```bash
   python hera_repo_cli.py add      --repo-name IMS      --repo-path /opt/hera/toolkits/IMS      --username <mongo-user>      --password <mongo-pass>      --auth-db dbhera      --db-prefix hera      --verbose
   ```
   Replace `--username/--password/--auth-db` according to your MongoDB.

4. **Verify**:
   ```bash
   python hera_repo_cli.py check      --username <mongo-user>      --password <mongo-pass>      --auth-db dbhera
   ```

---

## CLI Usage

### `add`
Attach (or upsert) a repo and add its id to all users across matching DBs.

```bash
python hera_repo_cli.py add   --repo-name <NAME>   --repo-path </path/to/repo>   --username <mongo-user>   --password <mongo-pass>   --auth-db <auth-db>   --db-prefix hera   [--host 127.0.0.1] [--port 27017]   [--dry-run] [--verbose]
```

**Common flags:**

- `--repo-name` – human-friendly name (e.g., `IMS`).
- `--repo-path` – filesystem/git path (e.g., `/opt/hera/toolkits/IMS`).
- `--username`, `--password`, `--auth-db` – MongoDB auth.
- `--db-prefix` – when scanning all DBs, restrict to DBs starting with this prefix.
- `--dry-run` – print what will happen without modifying data.
- `--verbose` – extra logs.

### `check`
Prints all users (from DB `hera`) and existing repository records.

```bash
python hera_repo_cli.py check   --username <mongo-user>   --password <mongo-pass>   --auth-db <auth-db>   [--host 127.0.0.1] [--port 27017]
```

> You can easily extend `check` to read from a different DB or to print additional diagnostics.

---

## Examples

**Dry run:**
```bash
python hera_repo_cli.py add   --repo-name IMS   --repo-path /opt/hera/toolkits/IMS   --username ilay   --password ilay-strong-pass   --auth-db dbhera   --db-prefix hera   --dry-run --verbose
```

**Apply changes:**
```bash
python hera_repo_cli.py add   --repo-name IMS   --repo-path /opt/hera/toolkits/IMS   --username ilay   --password ilay-strong-pass   --auth-db dbhera   --db-prefix hera   --verbose
```

**Check state:**
```bash
python hera_repo_cli.py check   --username ilay   --password ilay-strong-pass   --auth-db dbhera
```

---

## Troubleshooting

- **Authentication failed**: verify `--username/--password/--auth-db` and that the user has the right roles on the target DBs.
- **No users found**: make sure your `users` collection is populated in the DB(s) you target.
- **Nothing updated**: if you already ran once, a second run may update `0` users (since `$addToSet` is idempotent). This is expected.

---

## License

MIT (adapt to your project).