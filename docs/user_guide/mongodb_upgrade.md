# Upgrading MongoDB to 9.0

Hera now runs on MongoDB 9.0. If you already have a Hera database made by
**MongoDB 5.0** (what `init_with_mongo.sh`, the `Makefile` and `docker-compose.yml`
started before #1125) or **MongoDB 6.0**, the data has to be moved once with
`scripts/mongo/migrate_to_mongo9.sh`. Both source versions are tested end to end
into `mongo:9.0.4`.

---

## When you need this

You have an existing data dir (default `~/mongo-db-datadir`) and you switch to the
9.0 image. The 9.0 server will not open it and exits right away with:

```
UPGRADE PROBLEM: Found an invalid featureCompatibilityVersion document ... version: "5.0"
```

A new install with an empty data dir does not need this page.

MongoDB 9.0 only opens data at feature compatibility version (FCV) 8.0. Upgrading in place
would take one binary-and-FCV step per release (5.0 → 6.0 → 7.0 → 8.0 → 9.0). We can afford
downtime, so we do a single dump and restore into a new, empty 9.0 data dir instead. That
also gives fresh data files and freshly built indexes.

Only MongoDB is migrated, and MongoDB holds Hera's metadata documents. The data files that
documents point to through `resource` (parquet, netCDF, GeoTIFF, ...) stay where they are,
and their paths do not change.

## Downtime

Downtime starts when the script stops the old container and ends when the new one is
up. It reads every document three times (dump, source check, target check) and writes
it once, so it grows linearly with the number of documents. In the test run, 4,210 documents
over two databases took **9 seconds** from start to finish. To time a large database,
rehearse on a copy of the data dir first (see [Rehearsal](#rehearsal-optional)).

## Before you start

- **Stop everything that writes to Hera**: notebooks, `hera-server`, workflows, cron jobs.
  The script stops the MongoDB container anyway, but a client that is in the middle of
  writing will get errors.
- **Disk space**: about 2× the size of the data dir, on the disk that holds the new data
  dir: one for the new data files and one for the dump archive. The script checks this and
  refuses to start if there is not enough.
- **The 9.0 image**: use `mongo:9.0.4` or newer. Earlier 9.0.x crash on host kernels
  ≥ 6.19. If you already have an older `mongo:9.0` locally, run `docker pull mongo:9.0` first.
- Run the script from a checkout of this branch. It uses `mongo-init.d/` from the repo.

## Run it

```bash
cd ~/hera                                    # your Hera checkout
./scripts/mongo/migrate_to_mongo9.sh
```

All settings are environment variables. Their defaults match `init_with_mongo.sh` and the
`Makefile`:

| Variable | Default | Meaning |
|---|---|---|
| `MONGO_CONTAINER` | `hera-mongo` | container that serves Hera today |
| `MONGO_DATA` | `~/mongo-db-datadir` | old data dir, only read |
| `NEW_MONGO_DATA` | `$MONGO_DATA-9.0` | new data dir; must not exist yet |
| `MONGO_SOURCE_IMAGE` | image of `$MONGO_CONTAINER` | server that wrote `MONGO_DATA`, e.g. `mongo:5.0` |
| `MONGO_IMAGE` | `mongo:9.0` | new server; its tools do the dump, restore and checks |
| `MONGO_PORT` | `27017` | host port of the new container |
| `MONGO_BACKUP_DIR` | `~/mongo-backups` | where the dump archive is kept |

If `hera-mongo` no longer exists, or was already re-created with the 9.0 image (for example
by running the new `init_with_mongo.sh`, which removes the old container), set the source
version yourself:

```bash
MONGO_SOURCE_IMAGE=mongo:5.0 ./scripts/mongo/migrate_to_mongo9.sh
```

What the script does:

1. **Checks before it touches anything.** It refuses to run if `NEW_MONGO_DATA` exists, if
   `hera-mongo-old` exists (a migration already ran), or if there is not enough disk space.
2. **Stops writes.** It stops `hera-mongo` and starts the source image on the old data dir
   on a private Docker network with no published port, so nothing else can write to it.
3. **Dumps** everything (all databases and the users) to
   `~/mongo-backups/hera-mongo-<timestamp>.archive.gz` (`mongodump --archive --gzip`).
4. **Restores** into a new 9.0 server on `NEW_MONGO_DATA` with `--noIndexRestore`.
5. **Rebuilds the indexes** and **validates** (see below).
6. **Switches over**, only if the validation passed. It renames the old container to
   `hera-mongo-old` (stopped, kept for rollback) and starts a new `hera-mongo` from
   `MONGO_IMAGE` on `NEW_MONGO_DATA`, with the same port and `mongo-init.d` mount as
   `init_with_mongo.sh`.

If any step fails, the script removes its temporary containers and network, leaves the old
data dir and the old container as they were, and prints how to resume
(`docker start hera-mongo`). Before you run it again, remove the partial `NEW_MONGO_DATA`
with the command it prints.

`~/.pyhera/config.json` does not change: same host, same port, same users and passwords.

## Validate

The script validates on its own and exits non-zero on any difference. For every collection
in every database (except `local`, `config` and the FCV document in
`admin.system.version`) it compares the source and the target:

- the **document count**;
- a **sha256 of the raw BSON of every document in `_id` order**. This is a byte-for-byte
  check, so it also catches type-only changes such as an int that became a double;
- the **index specs** (keys and options).

It also prints the indexes on the target and its FCV:

```
rebuilt 2 index(es) on olymp.metadata: _cls_1_projectName_1, type_1
rebuilt 1 index(es) on olymp_shared.metadata: _cls_1_projectName_1
OK       admin.system.users  count=3  sha256=7f06545f4927d7c8  indexes=_id_,user_1_db_1
OK       olymp.metadata  count=3608  sha256=13f9ceb9706c87e8  indexes=_cls_1_projectName_1,_id_,type_1
OK       olymp_shared.metadata  count=602  sha256=e6d47ac72eeb427e  indexes=_cls_1_projectName_1,_id_
target FCV: 9.0
validation passed: every collection has the same documents and indexes
```

Then check through Hera. For a stronger check, run this before the migration as well and
compare the two outputs:

```python
from hera import Project
from hera.datalayer import getProjectList

for name in sorted(getProjectList()):
    project = Project(projectName=name)
    print(name, len(project.getMeasurementsDocuments()),
          len(project.getSimulationsDocuments()), len(project.getCacheDocuments()))
```

Also call `getData()` on a document or two that you know.

## Roll back

The old data dir was never written by 9.0, so it still opens with the old image. The script
prints the exact command:

```bash
docker rm -f hera-mongo && docker rename hera-mongo-old hera-mongo && docker start hera-mongo
```

Anything written to the 9.0 server after the migration is not in the old data dir. The
dump archive is a second copy of the state before the migration. To restore it, run
`mongorestore --archive --gzip` into an empty server of the old version.

## Finish

Until you are sure you will not roll back, keep both data dirs and use
`docker start/stop hera-mongo` or `make mongo-up` / `make mongo-down`. These reuse the
container and its mounts. Do **not** run these two in that period:

- `init_with_mongo.sh` re-creates `hera-mongo` on `~/mongo-db-datadir`, which is still the
  old data;
- `make mongo-clean` deletes `~/mongo-db-datadir`, which is your rollback.

When you no longer need the rollback, put the new data at the default path:

```bash
docker rm hera-mongo-old
docker stop hera-mongo && docker rm hera-mongo
mv ~/mongo-db-datadir ~/mongo-db-datadir-retired     # delete it later
mv ~/mongo-db-datadir-9.0 ~/mongo-db-datadir
make mongo-up
```

## Gotchas

**Feature compatibility version.** A new 9.0 data dir starts at FCV 9.0, and
`mongorestore` does not restore the FCV document from the dump. The validation prints
`target FCV: 9.0`. This is also why the rollback works: only the old data dir has the old FCV.

**Users and the `admin` database.** The dump includes `admin.system.users`. `mongorestore`
merges those users into the new server, so every user keeps its password, including users
added by hand. For this reason the restore target is started *without* `mongo-init.d`: the
users come from the dump, not from the init script. The validation compares
`admin.system.users` byte for byte. The script assumes Hera's default server, which runs
without `--auth` (that is how `init_with_mongo.sh`, the `Makefile` and `docker-compose.yml`
start it). On a server with `--auth`, the tool calls would need credentials, and the script
does not handle that.

**Index rebuild.** Hera's code declares one index: `indexes: ['projectName']` on the
`Metadata` document in `hera/datalayer/document/__init__.py`. Because the class uses
`allow_inheritance`, mongoengine stores it as `{_cls: 1, projectName: 1}`
(`_cls_1_projectName_1`). mongoengine also creates it by itself on first access
(`auto_create_indexes`). The script rebuilds the indexes **from the source's index
definitions, not from the code**:

- Data is restored with `--noIndexRestore`. Then every index is built new on 9.0 with
  `createIndexes`, once all the data is in. No index files or index versions come over
  from the old server.
- Options that only describe how the old server stored the index are dropped: `v`, `ns`,
  and `background` (mongoengine writes `background: false`, and 9.0 ignores it). The rebuilt
  indexes are therefore clean 9.0 indexes.
- Taking the definitions from the source and not the code keeps indexes that an admin
  added by hand for speed (for example on `type`). Rebuilding only from the code would drop
  them without telling anyone, and Hera would recreate its own index anyway. The validation
  checks that the target has exactly the source's indexes.

**"Cross-version dump & restore is unsupported".** `mongorestore` prints this warning
because the archive comes from 5.0/6.0. The byte-for-byte validation is what guards against
it. The warning is expected, and the script does not switch over unless every document
matches.

**Wrong source image.** If `MONGO_SOURCE_IMAGE` is newer than the data, the source server
will not start. The script stops at once and shows the `UPGRADE PROBLEM` log line, and
nothing is changed. Set `MONGO_SOURCE_IMAGE` to the version that wrote the data.

**docker compose.** The compose file runs MongoDB as the `mongo` service on
`./mongo-db-datadir`. To migrate, stop that service and run the script under a container
name that does not exist, so the script leaves the compose container alone. Then remove the
plain container the script started, and swap the dirs:

```bash
docker compose stop mongo
MONGO_CONTAINER=hera-mongo-migration MONGO_SOURCE_IMAGE=mongo:5.0 \
  MONGO_DATA="$PWD/mongo-db-datadir" MONGO_PORT=27099 ./scripts/mongo/migrate_to_mongo9.sh
docker rm -f hera-mongo-migration
mv mongo-db-datadir mongo-db-datadir-5.0 && mv mongo-db-datadir-9.0 mongo-db-datadir
docker compose up -d
```

To roll back, swap the two dirs again and start the stack with `MONGO_IMAGE=mongo:5.0`.

## Rehearsal (optional)

To measure the downtime or practise without touching production, copy the stopped data dir
and migrate the copy under another container name and port:

```bash
docker stop hera-mongo
sudo cp -a ~/mongo-db-datadir /tmp/rehearsal-data
docker start hera-mongo
time MONGO_CONTAINER=hera-mongo-rehearsal MONGO_SOURCE_IMAGE=mongo:5.0 \
  MONGO_DATA=/tmp/rehearsal-data MONGO_PORT=27099 ./scripts/mongo/migrate_to_mongo9.sh
docker rm -f hera-mongo-rehearsal
```
