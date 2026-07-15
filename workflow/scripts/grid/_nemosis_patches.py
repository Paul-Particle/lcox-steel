"""Monkey-patches for NEMOSIS to work around AEMO API changes.

Applied on import via `import _nemosis_patches` in download_nem.py.
Delete this file and that import once the upstream package ships fixes.

Patch 1: AEMO now blocks NEMOSIS's stale Chrome 80 User-Agent with HTTP 403.

Patch 2: MMSDM filenames from Aug-2024 contain literal '#' characters
  (e.g. PUBLIC_ARCHIVE#DISPATCHPRICE#FILE01#202408010000.zip). The nemweb.com.au
  MMSDM archive serves these when the '#' is percent-encoded once as '%23'
  (a literal '#' is treated as a URL fragment and stripped, giving a 404).
  requests re-normalises the path and strips the fragment, so we bypass it with
  http.client, which sends the pre-encoded path verbatim.
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
    """Download and extract a NEMOSIS MMSDM zip, single-encoding '#' as '%23' (Patch 2)."""
    from urllib.parse import urlsplit
    url_fixed = url.replace("#", "%23")
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
