const express = require("express");
const fs = require("fs");
const path = require("path");
const https = require("https");

const app = express();

const sslOptions = {
  key: fs.readFileSync(path.resolve(__dirname, "ssl", "server.key")),
  cert: fs.readFileSync(path.resolve(__dirname, "ssl", "server.crt")),
  passphrase: "test123", // Ук
};

app.use((req, res, next) => {
  res.setHeader("Cross-Origin-Opener-Policy", "same-origin");
  res.setHeader("Cross-Origin-Embedder-Policy", "require-corp");
  res.setHeader("Cross-Origin-Resource-Policy", "same-site");
  next();
});

app.use(express.static(path.join(__dirname, "dist")));

// app.get("/", (req, res, next) => {
//   const headerValue = req.headers["custom-header"];

//   //   console.log(req.headers);

//   if (headerValue === "value2") {
//     res.sendFile(path.join(__dirname, "index2.html"));
//   } else {
//     res.sendFile(path.join(__dirname, "index1.html"));
//     res.setHeader("Cross-Origin-Opener-Policy", "same-origin");
//     res.setHeader("Cross-Origin-Embedder-Policy", "require-corp");
//   }
// });

app.get("/your-proxy-route/", async (req, res) => {
  res.setHeader("Cross-Origin-Opener-Policy", "unsafe-none");
  res.setHeader("Cross-Origin-Embedder-Policy", "unsafe-none");
  res.sendFile(path.join(__dirname, "index1.html"));
});

// app.listen(80, () => {
//   console.log(`Сервер запущен на http://127.0.0.1:${80}`);
// });

https.createServer(sslOptions, app).listen(443, () => {
  console.log(`Сервер запущен на https://127.0.0.1:${443}`);
});
