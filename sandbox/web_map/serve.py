"""Serve index.html and the generated out/ packages on 127.0.0.1 only."""

import argparse
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent


class Handler(SimpleHTTPRequestHandler):
    """Serves "/" as index.html and files under out/; everything else is 404."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(HERE), **kwargs)

    def translate_path(self, path):
        if path.split("?")[0] == "/":
            path = "/index.html"
        resolved = Path(super().translate_path(path)).resolve()
        out = HERE / "out"
        if resolved == HERE / "index.html" or resolved.is_relative_to(out):
            return str(resolved)
        return str(HERE / "does-not-exist")

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    with ThreadingHTTPServer(("127.0.0.1", args.port), Handler) as server:
        print(f"Open http://127.0.0.1:{server.server_port}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
