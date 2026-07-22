
// use admin;
db = db.getSiblingDB("admin");

var adminUser = process.env.MONGO_ADMIN_USER || "MathAdmin";
var adminPwd  = process.env.MONGO_ADMIN_PWD  || "MathAdmin";
var heraUser  = process.env.MONGO_HERA_USER  || "hera";
var heraPwd   = process.env.MONGO_HERA_PWD   || "heracles";

db.getUser(adminUser) || db.createUser(
  {
    user: adminUser,
    pwd: adminPwd,
    roles: [ { role: "userAdminAnyDatabase", db: "admin" } , "readWriteAnyDatabase"]
  }
);

// use admin;
db.getUser(heraUser) || db.createUser(
  {
    user: heraUser,
    pwd:  heraPwd,
    roles: [ { role: "readWrite", db: "olymp" } ]
  }
);

printjson(db.getUsers())
