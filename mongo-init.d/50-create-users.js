
// use admin;
db = db.getSiblingDB("admin");

db.getUser("MathAdmin") || db.createUser(
  {
    user: "MathAdmin",
    pwd: "MathAdmin",
    roles: [ { role: "userAdminAnyDatabase", db: "admin" } , "readWriteAnyDatabase"]
  }
);

// use admin;
db.getUser("hera") || db.createUser(
  {
    user: "hera",
    pwd:  "heracles",   
    roles: [ { role: "readWrite", db: "olymp" } ]

  }
);

printjson(db.getUsers())
