const express = require("express");
const https = require("https");
const path = require("path");
const fs = require("fs");

const app = express();

app.get("/", (req, res) => {
  res.send("Hello, HTTPS World!");
});

const sslOptions = {
  key: fs.readFileSync(path.resolve(__dirname, "ssl", "server.key")),
  cert: fs.readFileSync(path.resolve(__dirname, "ssl", "server.crt")),
  passphrase: "test123", // Ук
};

https.createServer(sslOptions, app).listen(443, () => {
  console.log("Express HTTPS server listening on port 443");
});
