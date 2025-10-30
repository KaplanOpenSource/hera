import click
import pymongo
import subprocess
from datetime import datetime


@click.group()
def cli():
    """CLI לניהול ריפוזיטוריז ומשתמשים במערכת Hera"""
    pass


@cli.command()
@click.option("--repo-name", required=True, help="שם הריפו (למשל IMS)")
@click.option("--repo-path", required=True, help="נתיב הריפו (למשל /opt/hera/toolkits/IMS)")
@click.option("--host", default="127.0.0.1")
@click.option("--port", default=27017)
@click.option("--username", required=True)
@click.option("--password", required=True)
@click.option("--auth-db", default="dbhera")
@click.option("--db-prefix", default="hera")
@click.option("--dry-run", is_flag=True, help="להריץ כ־dry run בלבד")
@click.option("--verbose", is_flag=True)
def add(repo_name, repo_path, host, port, username, password, auth_db, db_prefix, dry_run, verbose):
    """הוספת ריפו לכל המשתמשים"""
    cmd = [
        "python", "add_repo_to_all_users.py",
        "--host", host,
        "--port", str(port),
        "--username", username,
        "--password", password,
        "--auth-db", auth_db,
        "--scan-all", "--db-prefix", db_prefix,
        "--users-collection", "users",
        "--repos-collection", "repositories",
        "--user-repos-field", "repos",
        "--repo-path", repo_path,
        "--repo-name", repo_name
    ]
    if dry_run:
        cmd.append("--dry-run")
    if verbose:
        cmd.append("--verbose")

    subprocess.run(cmd)


@cli.command()
@click.option("--host", default="127.0.0.1")
@click.option("--port", default=27017)
@click.option("--username", required=True)
@click.option("--password", required=True)
@click.option("--auth-db", default="dbhera")
def check(host, port, username, password, auth_db):
    """בודקת שה־repo נוסף לכל המשתמשים"""
    client = pymongo.MongoClient(
        host=host,
        port=port,
        username=username,
        password=password,
        authSource=auth_db
    )
    db = client["hera"]

    print("\n📋 רשימת משתמשים והריפוזיטוריז שלהם:\n")
    for u in db.users.find({}, {"username": 1, "repos": 1, "_id": 0}):
        print(f"👤 {u['username']}: {u.get('repos', [])}")

    print("\n📦 רשימת ריפוזיטוריז קיימים:\n")
    for r in db.repositories.find({}, {"_id": 0}):
        r["createdAt"] = r["createdAt"].strftime("%Y-%m-%d %H:%M:%S") if "createdAt" in r else ""
        print(f"• {r['name']} | {r['path']} | {r.get('repo_id','')}")

    print("\n✅ בדיקה הסתיימה בהצלחה!\n")


if __name__ == "__main__":
    cli()
