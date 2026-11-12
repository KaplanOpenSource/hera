#!/usr/bin/env bash
# Migrate Hera's MongoDB from a 5.0 or 6.0 data dir to a fresh 9.0 data dir (issue #1126).
#
# MongoDB 9.0 cannot open a data dir written by 5.0/6.0, so the data goes through
# mongodump/mongorestore. The old data dir is only read, never changed: it is the rollback.
# See docs/user_guide/mongodb_upgrade.md for the full procedure.
#
# Usage:
#   ./scripts/mongo/migrate_to_mongo9.sh
#
# Settings (env vars; defaults match init_with_mongo.sh and the Makefile):
#   MONGO_CONTAINER     hera-mongo              container that serves Hera now
#   MONGO_DATA          $HOME/mongo-db-datadir  old data dir (read only)
#   NEW_MONGO_DATA      $MONGO_DATA-9.0         new data dir, must not exist yet
#   MONGO_SOURCE_IMAGE  image of MONGO_CONTAINER  server that wrote MONGO_DATA, e.g. mongo:5.0
#   MONGO_IMAGE         mongo:9.0               new server; its tools do all the work
#   MONGO_PORT          27017                   host port of the new container
#   MONGO_BACKUP_DIR    $HOME/mongo-backups     where the dump archive is kept
#
# Steps: stop the old container, dump it, restore into a new 9.0 container without
# indexes, rebuild the indexes, compare every collection (count + sha256 of the raw
# documents in _id order + index specs), and only then start the new 9.0 container
# as $MONGO_CONTAINER. The old container is kept, stopped, as $MONGO_CONTAINER-old.
# Any failure leaves the old container stopped but untouched: `docker start $MONGO_CONTAINER`.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INIT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)/mongo-init.d"

MONGO_CONTAINER="${MONGO_CONTAINER:-hera-mongo}"
MONGO_DATA="${MONGO_DATA:-${HOME}/mongo-db-datadir}"
NEW_MONGO_DATA="${NEW_MONGO_DATA:-${MONGO_DATA}-9.0}"
MONGO_IMAGE="${MONGO_IMAGE:-mongo:9.0}"
MONGO_PORT="${MONGO_PORT:-27017}"
MONGO_BACKUP_DIR="${MONGO_BACKUP_DIR:-${HOME}/mongo-backups}"

OLD_CONTAINER="${MONGO_CONTAINER}-old"
NET="${MONGO_CONTAINER}-migrate"
SRC="${MONGO_CONTAINER}-migrate-src"
DST="${MONGO_CONTAINER}-migrate-dst"
ARCHIVE="${MONGO_BACKUP_DIR}/${MONGO_CONTAINER}-$(date +%Y%m%d-%H%M%S).archive.gz"

# Runs in mongosh connected to the target; reads the source at SRC_URI.
# Rebuilds the source's indexes on the target, then compares both servers.
CHECK_JS=$(cat <<'EOF'
const crypto = require('crypto');
const source = new Mongo(SRC_URI), target = db.getMongo();

function snapshot(mongo) {
  const out = {};
  mongo.getDB('admin').adminCommand({listDatabases: 1, nameOnly: true}).databases
    .map(d => d.name).filter(name => name !== 'local' && name !== 'config').sort()
    .forEach(dbName => {
      const database = mongo.getDB(dbName);
      database.getCollectionInfos({type: 'collection'}).map(c => c.name).sort().forEach(name => {
        const ns = `${dbName}.${name}`;
        if (ns === 'admin.system.version') return;  // holds the FCV, differs on purpose
        const coll = database.getCollection(name), hash = crypto.createHash('sha256');
        let count = 0;
        // raw BSON bytes: catches type-only changes (int32 vs double) that parsed documents hide
        coll.find({}, {}, {raw: true}).sort({_id: 1}).forEach(bson => { hash.update(bson); count++; });
        // v/ns/background are storage details of the old server, not part of the index.
        const indexes = coll.getIndexes()
          .map(({v, ns, background, ...spec}) => Object.fromEntries(Object.entries(spec).sort()))
          .sort((a, b) => a.name.localeCompare(b.name));
        out[ns] = {count, sha256: hash.digest('hex'), indexes};
      });
    });
  return out;
}

const before = snapshot(source);
for (const [ns, {indexes}] of Object.entries(before)) {
  const [dbName, ...rest] = ns.split('.'), collName = rest.join('.');
  const specs = indexes.filter(spec => spec.name !== '_id_');
  if (collName.startsWith('system.') || specs.length === 0) continue;  // system indexes are the server's
  target.getDB(dbName).runCommand({createIndexes: collName, indexes: specs});
  print(`rebuilt ${specs.length} index(es) on ${ns}: ${specs.map(s => s.name).join(', ')}`);
}

const after = snapshot(target);
let bad = 0;
for (const ns of [...new Set([...Object.keys(before), ...Object.keys(after)])].sort()) {
  const a = EJSON.stringify(before[ns] ?? null), b = EJSON.stringify(after[ns] ?? null);
  const row = after[ns] || before[ns];
  print(`${a === b ? 'OK      ' : 'MISMATCH'} ${ns}  count=${row.count}  sha256=${row.sha256.slice(0, 16)}  indexes=${row.indexes.map(i => i.name).join(',')}`);
  if (a !== b) { bad++; print(`  source: ${a}\n  target: ${b}`); }
}
print(`target FCV: ${target.getDB('admin').adminCommand({getParameter: 1, featureCompatibilityVersion: 1}).featureCompatibilityVersion.version}`);
if (bad) { print(`${bad} namespace(s) differ`); quit(1); }
print('validation passed: every collection has the same documents and indexes');
EOF
)

log() { echo "=== $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
exists() { docker ps -a --format '{{.Names}}' | grep -qx "$1"; }

wait_ready() {  # the entrypoint execs mongod as PID 1 only after any init scripts ran
    for _ in $(seq 1 60); do
        if [ "$(docker exec "$1" cat /proc/1/comm 2>/dev/null)" = mongod ] &&
           docker exec "$1" mongosh --quiet --norc --eval 'db.runCommand({ping: 1}).ok' >/dev/null 2>&1; then
            return 0
        fi
        [ "$(docker inspect -f '{{.State.Running}}' "$1")" = true ] || break
        sleep 1
    done
    docker logs "$1" 2>&1 | grep -E '"s":"[EF]"' | tail -5 >&2 || docker logs --tail 20 "$1" >&2
    die "$1 did not become ready (wrong MONGO_SOURCE_IMAGE? see the log lines above)"
}

cleanup() {
    local status=$?
    docker rm -f "$SRC" "$DST" >/dev/null 2>&1 || true
    docker network rm "$NET" >/dev/null 2>&1 || true
    if [ "$status" -ne 0 ] && [ "${switched:-no}" = yes ]; then
        echo "The new server did not come up. Roll back with:"
        echo "  docker rm -f ${MONGO_CONTAINER} && docker rename ${OLD_CONTAINER} ${MONGO_CONTAINER} && docker start ${MONGO_CONTAINER}"
    elif [ "$status" -ne 0 ]; then
        echo "Migration aborted. ${MONGO_DATA} was not changed."
        echo "  Resume the old server:  docker start ${MONGO_CONTAINER}"
        if [ -e "$NEW_MONGO_DATA" ]; then
            echo "  Before re-running, remove the partial ${NEW_MONGO_DATA}:"
            echo "    docker run --rm -v ${NEW_MONGO_DATA}:/data/db ${MONGO_IMAGE} bash -c 'rm -rf /data/db/*' && rmdir ${NEW_MONGO_DATA}"
        fi
    fi
    exit "$status"
}

# --- preflight: refuse anything that could overwrite data ---
[ -d "$MONGO_DATA" ] || die "old data dir ${MONGO_DATA} not found"
[ -e "$NEW_MONGO_DATA" ] && die "${NEW_MONGO_DATA} already exists; refusing to overwrite it"
exists "$OLD_CONTAINER" && die "container ${OLD_CONTAINER} exists: a migration already ran"
for c in "$SRC" "$DST"; do exists "$c" && die "leftover container $c; remove it first"; done
[ -d "$INIT_DIR" ] || die "${INIT_DIR} not found"
if [ -z "${MONGO_SOURCE_IMAGE:-}" ] && exists "$MONGO_CONTAINER"; then
    MONGO_SOURCE_IMAGE="$(docker inspect -f '{{.Config.Image}}' "$MONGO_CONTAINER")"
fi
[ -n "${MONGO_SOURCE_IMAGE:-}" ] || die "set MONGO_SOURCE_IMAGE to the image that wrote ${MONGO_DATA}, e.g. mongo:5.0"

data_bytes=$(docker run --rm -v "${MONGO_DATA}:/data/db:ro" --entrypoint du "$MONGO_SOURCE_IMAGE" -sb /data/db | cut -f1)
mkdir -p "$MONGO_BACKUP_DIR" "$(dirname "$NEW_MONGO_DATA")"
free_bytes=$(df -B1 --output=avail "$(dirname "$NEW_MONGO_DATA")" | tail -1)
[ "$free_bytes" -gt $((2 * data_bytes)) ] ||
    die "need about 2x ${data_bytes} bytes free next to ${NEW_MONGO_DATA}, have ${free_bytes}"

trap cleanup EXIT

# --- 1. stop writes: stop the old server and serve its data dir on a private network only ---
if exists "$MONGO_CONTAINER"; then
    log "stopping ${MONGO_CONTAINER} (downtime starts)"
    docker stop -t 60 "$MONGO_CONTAINER" >/dev/null
fi
docker network create "$NET" >/dev/null
log "starting ${MONGO_SOURCE_IMAGE} on ${MONGO_DATA} (no published port)"
docker run -d --name "$SRC" --network "$NET" -v "${MONGO_DATA}:/data/db" "$MONGO_SOURCE_IMAGE" >/dev/null
wait_ready "$SRC"
echo "source server $(docker exec "$SRC" mongosh --quiet --norc --eval \
    'db.version() + ", FCV " + db.adminCommand({getParameter: 1, featureCompatibilityVersion: 1}).featureCompatibilityVersion.version')"

# --- 2. dump everything (users included) ---
log "dumping to ${ARCHIVE}"
docker run --rm --network "$NET" "$MONGO_IMAGE" mongodump --uri "mongodb://${SRC}/" --archive --gzip > "$ARCHIVE"

# --- 3. restore into a fresh 9.0 data dir, data and users only ---
log "starting ${MONGO_IMAGE} on new ${NEW_MONGO_DATA}"
mkdir "$NEW_MONGO_DATA"
docker run -d --name "$DST" --network "$NET" -v "${NEW_MONGO_DATA}:/data/db" "$MONGO_IMAGE" >/dev/null
wait_ready "$DST"
log "restoring (indexes are rebuilt in the next step)"
docker run --rm -i --network "$NET" "$MONGO_IMAGE" \
    mongorestore --uri "mongodb://${DST}/" --archive --gzip --noIndexRestore < "$ARCHIVE"

# --- 4. rebuild indexes and validate ---
log "rebuilding indexes and validating"
docker run --rm --network "$NET" "$MONGO_IMAGE" \
    mongosh --quiet --norc "mongodb://${DST}/" --eval "const SRC_URI = 'mongodb://${SRC}/'; ${CHECK_JS}"

# --- 5. switch: keep the old container as the rollback, start 9.0 as Hera's server ---
docker stop -t 60 "$SRC" "$DST" >/dev/null
exists "$MONGO_CONTAINER" && docker rename "$MONGO_CONTAINER" "$OLD_CONTAINER"
switched=yes
log "starting ${MONGO_CONTAINER} (${MONGO_IMAGE}) on 127.0.0.1:${MONGO_PORT}"
docker run --name "$MONGO_CONTAINER" \
    -v "${NEW_MONGO_DATA}:/data/db" \
    -v "${INIT_DIR}:/docker-entrypoint-initdb.d" \
    -p "127.0.0.1:${MONGO_PORT}:27017" \
    -d "$MONGO_IMAGE" >/dev/null
wait_ready "$MONGO_CONTAINER"

log "done (downtime ends)"
echo "  new data dir: ${NEW_MONGO_DATA}"
echo "  dump archive: ${ARCHIVE}"
echo "  rollback:     docker rm -f ${MONGO_CONTAINER} && docker rename ${OLD_CONTAINER} ${MONGO_CONTAINER} && docker start ${MONGO_CONTAINER}"
