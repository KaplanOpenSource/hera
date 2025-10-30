#!/usr/bin/env python3
"""
add_repo_to_all_users_multidb.py  (refactored)
----------------------------------------------
Admin script for environments where users are spread across multiple MongoDB databases
OR a single shared DB. Attaches a repository to all users (idempotent via $addToSet).

Examples (NO URI, like your local flow):
  python add_repo_to_all_users_multidb.py \
    --host 127.0.0.1 --port 27017 \
    --username ilay --password ilay2899 --auth-db dbhera \
    --db-list hera \
    --users-collection users \
    --repos-collection repositories \
    --user-repos-field repos \
    --repo-path /opt/hera/toolkits/IMS \
    --repo-name IMS \
    --verbose

Using URI (optional alternative):
  python add_repo_to_all_users_multidb.py \
    --mongo-uri "mongodb://Admin:Admin@127.0.0.1:27017/?authSource=admin" \
    --db-list hera \
    --repo-path /opt/hera/toolkits/IMS \
    --repo-name IMS
"""

import argparse
import hashlib
import re
import sys
from datetime import datetime, timezone
from typing import List, Set, Optional, Dict, Any

from pymongo import MongoClient, ASCENDING
from pymongo.errors import PyMongoError, OperationFailure

SYSTEM_DBS: Set[str] = {"admin", "local", "config"}

def stable_id_from_path(repo_path: str) -> str:
    norm = repo_path.strip().rstrip('/')
    h = hashlib.sha1(norm.encode('utf-8')).hexdigest()[:12]
    return f"repo-{h}"

def filter_dbs(all_dbs: List[str], prefix: Optional[str] = None, regex: Optional[str] = None) -> List[str]:
    candidates = [d for d in all_dbs if d not in SYSTEM_DBS]
    if prefix:
        candidates = [d for d in candidates if d.startswith(prefix)]
    if regex:
        r = re.compile(regex)
        candidates = [d for d in candidates if r.search(d)]
    return sorted(candidates)

def build_client_from_args(args: argparse.Namespace) -> MongoClient:
    if args.mongo_uri:
        return MongoClient(args.mongo_uri)
    # “כמו שעשית מקודם”: host/port/username/password/authSource
    kwargs: Dict[str, Any] = dict(host=args.host, port=args.port)
    if args.username:
        kwargs["username"] = args.username
    if args.password:
        kwargs["password"] = args.password
    if args.auth_db:
        kwargs["authSource"] = args.auth_db
    return MongoClient(**kwargs)

def main():
    p = argparse.ArgumentParser(description="Attach a repository toolkit to users across MongoDB DBs (idempotent).")

    # Connection (either --mongo-uri OR the discrete params below)
    p.add_argument("--mongo-uri", help="MongoDB URI with sufficient privileges (e.g., admin).")

    p.add_argument("--host", default="127.0.0.1", help="MongoDB host (default: 127.0.0.1)")
    p.add_argument("--port", type=int, default=27017, help="MongoDB port (default: 27017)")
    p.add_argument("--username", help="MongoDB username (optional if no auth)")
    p.add_argument("--password", help="MongoDB password (optional if no auth)")
    p.add_argument("--auth-db", help="Authentication database (e.g., admin / dbhera)")

    # Target DBs
    p.add_argument("--scan-all", action="store_true", help="Scan all non-system DBs (filter with --db-prefix / --db-regex).")
    p.add_argument("--db-prefix", default=None, help="Filter DBs by prefix (with --scan-all).")
    p.add_argument("--db-regex", default=None, help="Filter DBs by regex (with --scan-all).")
    p.add_argument("--db-list", default=None, help="Comma-separated DB names (e.g., hera_user1,hera_user2).")

    # Collections & fields
    p.add_argument("--users-collection", default="users", help='Users collection name (default: "users")')
    p.add_argument("--repos-collection", default="repositories", help='Repositories registry collection (default: "repositories")')
    p.add_argument("--user-repos-field", default="repos", help='Name of array field on user docs to store repo ids (default: "repos")')

    # Repo inputs
    p.add_argument("--repo-path", required=True, help="Filesystem path or git clone location of the repository")
    p.add_argument("--repo-name", default=None, help="Optional friendly name")
    p.add_argument("--repo-id", default=None, help="Optional stable repo id; if omitted derived from --repo-path")

    # Execution mode
    p.add_argument("--dry-run", action="store_true", help="Plan only, do not write")
    p.add_argument("--verbose", action="store_true", help="Verbose logs")

    args = p.parse_args()

    # Sanity: must choose *either* URI or (host/port) style; both are okay too, URI wins.
    if not args.mongo_uri and not args.host:
        print("[error] Provide --mongo-uri OR at least --host/--port (and auth if needed).", file=sys.stderr)
        sys.exit(2)

    # Figure repo_id
    repo_id = args.repo_id or stable_id_from_path(args.repo_path)
    user_repos_field = args.user_repos_field

    try:
        client = build_client_from_args(args)

        # Determine target DBs
        if args.db_list:
            target_dbs = [d.strip() for d in args.db_list.split(",") if d.strip()]
        elif args.scan_all:
            all_dbs = client.list_database_names()
            target_dbs = filter_dbs(all_dbs, prefix=args.db_prefix, regex=args.db_regex)
        else:
            print("[error] You must specify either --scan-all (optionally with filters) or --db-list.", file=sys.stderr)
            sys.exit(2)

        if not target_dbs:
            print("[info] No target databases after filtering. Nothing to do.")
            sys.exit(0)

        grand_total_users = 0
        grand_updated_users = 0
        touched_dbs = 0

        for dbname in target_dbs:
            try:
                db = client[dbname]
                users_col = db[args.users_collection]
                repos_col = db[args.repos_collection]

                # quick existence check
                users_count_est = users_col.estimated_document_count()
                if users_count_est == 0:
                    if args.verbose:
                        print(f"[skip] DB '{dbname}': no users found in '{args.users_collection}'.")
                    continue

                # ensure unique index on repositories.repo_id
                try:
                    repos_col.create_index([("repo_id", ASCENDING)], unique=True, name="uniq_repo_id")
                except Exception as e:
                    if args.verbose:
                        print(f"[warn] DB '{dbname}': could not create index uniq_repo_id: {e}", file=sys.stderr)

                now = datetime.now(timezone.utc)
                repo_doc = {
                    "repo_id": repo_id,
                    "path": args.repo_path.strip().rstrip('/'),
                    **({"name": args.repo_name} if args.repo_name else {}),
                    "updatedAt": now,
                }

                if args.verbose:
                    print(f"[info] DB '{dbname}': Upserting repo {repo_doc}")

                if not args.dry_run:
                    repos_col.update_one(
                        {"repo_id": repo_id},
                        {"$set": repo_doc, "$setOnInsert": {"createdAt": now}},
                        upsert=True
                    )

                # AddToSet repo_id into users.<user_repos_field>
                update_op = {"$addToSet": {user_repos_field: repo_id}}

                # For dry-run estimation (who is missing it)
                query_missing = {
                    "$or": [
                        {user_repos_field: {"$exists": False}},
                        {user_repos_field: {"$nin": [repo_id]}},
                    ]
                }
                missing_count = users_col.count_documents(query_missing)

                if args.dry_run:
                    if args.verbose:
                        print(f"[dry-run] DB '{dbname}': users needing update: {missing_count} / ~{users_count_est}")
                else:
                    result = users_col.update_many({}, update_op)
                    if args.verbose:
                        print(f"[info] DB '{dbname}': matched={result.matched_count}, modified={result.modified_count}")
                    grand_updated_users += result.modified_count

                grand_total_users += users_count_est
                touched_dbs += 1

            except OperationFailure as e:
                print(f"[warn] Skipping DB '{dbname}' due to OperationFailure: {e}", file=sys.stderr)
            except PyMongoError as e:
                print(f"[warn] Skipping DB '{dbname}' due to PyMongoError: {e}", file=sys.stderr)

        # summary
        mode = "DRY-RUN" if args.dry_run else "APPLIED"
        print("\n=== SUMMARY ({}) ===".format(mode))
        print(f"Target DBs considered: {len(target_dbs)} | DBs with users processed: {touched_dbs}")
        print(f"Repo: repo_id='{repo_id}', path='{args.repo_path}'" + (f", name='{args.repo_name}'" if args.repo_name else ""))
        print(f"Users total (estimated across DBs): ~{grand_total_users}")
        if args.dry_run:
            print("No writes performed.")
        else:
            print(f"Users updated across DBs: {grand_updated_users}")
        print("Done.")

    except PyMongoError as e:
        print(f"[error] MongoDB operation failed at top-level: {e}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()
