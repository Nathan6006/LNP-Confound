#!/usr/bin/env python3
"""Localhost receiver: the browser POSTs file bytes here and we write them.

Why this exists: the Playwright MCP server intercepts and cancels downloads
inside browser_run_code_unsafe, and the publisher session cookies that would
let curl fetch directly are httpOnly. So the browser -- which has the cookies
-- fetches the asset and POSTs the bytes to this process, which writes them.
No session secret and no multi-megabyte payload passes through the model's
context.

    python3 -I receiver.py <root-dir> <port>

POST /put?path=<relative/path/under/root>   body = raw file bytes
"""
import json, os, sys, hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

ROOT = os.path.abspath(sys.argv[1])
PORT = int(sys.argv[2]) if len(sys.argv) > 2 else 8765
MANIFEST = os.path.join(ROOT, "_received.jsonl")


class H(BaseHTTPRequestHandler):
    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")

    def do_OPTIONS(self):
        self.send_response(200); self._cors(); self.end_headers()

    def do_POST(self):
        q = parse_qs(urlparse(self.path).query)
        rel = (q.get("path") or [""])[0]
        src = (q.get("src") or [""])[0]
        # Refuse anything that escapes the root. Downloaded files are
        # untrusted and so is a path that came back through a web page.
        dest = os.path.abspath(os.path.join(ROOT, rel))
        if not dest.startswith(ROOT + os.sep) or not rel:
            self.send_response(400); self._cors(); self.end_headers()
            self.wfile.write(b"bad path"); return
        n = int(self.headers.get("Content-Length") or 0)
        data = self.rfile.read(n)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with open(dest, "wb") as f:
            f.write(data)
        rec = {"path": rel, "bytes": len(data), "src": src,
               "sha256": hashlib.sha256(data).hexdigest(),
               "magic": data[:8].decode("latin1", "replace")}
        with open(MANIFEST, "a") as f:
            f.write(json.dumps(rec) + "\n")
        self.send_response(200); self._cors()
        self.send_header("Content-Type", "application/json"); self.end_headers()
        self.wfile.write(json.dumps(rec).encode())
        print(f"  <- {rel}  {len(data)} bytes", flush=True)

    def log_message(self, *a):
        pass


if __name__ == "__main__":
    os.makedirs(ROOT, exist_ok=True)
    print(f"receiver on 127.0.0.1:{PORT} -> {ROOT}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", PORT), H).serve_forever()
