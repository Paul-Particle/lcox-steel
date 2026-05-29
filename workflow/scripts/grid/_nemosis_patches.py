"""Monkey-patches for NEMOSIS to work around AEMO API changes.

Applied on import via `import _nemosis_patches` in download_nem.py.
Delete this file and that import once the upstream package ships fixes.

Patch 1: AEMO now blocks NEMOSIS's stale Chrome 80 User-Agent with HTTP 403.

Patch 2: MMSDM filenames from Aug-2024 contain literal '#' characters.
  requests percent-encodes '#' → '%23', which AEMO's Azure endpoint rejects (HTTP 400).
  Sending a literal '#' also fails (proxies strip it as a fragment). The only working
  encoding is double-encoding '%2523'; Azure decodes one layer to '%23' as expected.
  requests would re-normalise '%2523' → '%23', so we bypass it entirely with http.client.
"""

import http.client
import io
import ssl
import zipfile

import nemosis.downloader as _dl

_dl.USR_AGENT_HEADER.clear()
_dl.USR_AGENT_HEADER["User-Agent"] = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
)


def _download_unzip_csv_patched(url: str, down_load_to: str) -> None:
    from urllib.parse import urlsplit
    url_fixed = url.replace("#", "%2523").replace("%23", "%2523")
    u = urlsplit(url_fixed)
    path = u.path + (f"?{u.query}" if u.query else "")
    conn = http.client.HTTPSConnection(
        u.hostname, u.port or 443, context=ssl.create_default_context(), timeout=120  # pyright: ignore
    )
    try:
        conn.request("GET", path,
                     headers={**_dl.USR_AGENT_HEADER, "Host": u.hostname, "Accept": "*/*"})  # pyright: ignore
        resp = conn.getresponse()
        if resp.status != 200:
            raise IOError(f"GET {url} -> {resp.status} {resp.reason}")
        body = resp.read()
    finally:
        conn.close()
    zipfile.ZipFile(io.BytesIO(body)).extractall(down_load_to)


_dl.download_unzip_csv = _download_unzip_csv_patched
