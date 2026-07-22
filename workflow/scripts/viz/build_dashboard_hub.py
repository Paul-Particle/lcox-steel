#!/usr/bin/env python3
"""Combine the v2 scenario-comparison dashboard, the network schematic and the
workflow DAG into ONE body-only artifact: a work-in-progress banner and a tab bar
over three isolated iframes.

Inputs:
  * results/dashboard_v2.html            — built by build_dashboard_v2.py
  * hub_pages/network_schematic.html     — static source page (no generator)
  * hub_pages/workflow_dag.html          — static source page (no generator)

Each page is a body-only, theme-aware Plotly document carrying its own
<style>/<script>, so we isolate them in `iframe.srcdoc` documents to avoid
CSS/ID/Plotly collisions. The three docs are embedded as JS string literals
(JSON-encoded, with `</` neutralised to `<\\/` so a nested `</script>` can't
terminate the host <script> element) and assigned to the iframe lazily on first
tab activation. The outer wrapper propagates the viewer's light/dark theme into
the active iframe.

Output is written to results/dashboard_hub.html (git-ignored, like the other
built dashboards). Publish that file to the artifact.
"""
import json
import re
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]        # workflow/scripts/viz/ -> repo root
RESULTS = REPO / "results"
HUB_PAGES = Path(__file__).with_name("hub_pages")  # tracked static source pages
OUT_PATH = RESULTS / "dashboard_hub.html"

# (key, path, tab label). "compare" is the dashboard build artifact; the other two
# are static source pages tracked under hub_pages/.
FRAGMENTS = {
    "compare":   (RESULTS / "dashboard_v2.html",         "Scenario comparison"),
    "schematic": (HUB_PAGES / "network_schematic.html",  "Network schematic"),
    "dag":       (HUB_PAGES / "workflow_dag.html",        "Workflow DAG"),
}
TAB_ORDER = ["compare", "schematic", "dag"]
DEFAULT = "compare"


def wrap(fragment_html: str, title: str) -> str:
    """Wrap a body-only fragment in a minimal full HTML document for the iframe."""
    return (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
        f"<title>{title}</title>"
        "<style>html,body{margin:0;padding:0}</style></head>"
        f"<body>{fragment_html}</body></html>"
    )


def js_string(s: str) -> str:
    """A JS string literal whose nested `</script>` can't close the host script."""
    return json.dumps(s).replace("</", "<\\/")


def build_hub() -> str:
    build_date = date.today().strftime("%-d %b %Y")   # stamped into the WIP banner

    docs, raw = {}, {}
    for key, (path, title) in FRAGMENTS.items():
        raw[key] = path.read_text(encoding="utf-8")
        docs[key] = wrap(raw[key], title)

    # The outer wrapper (banner + tab bar) styles itself with 'Titillium Web', but
    # the font data only lives inside the iframes. Lift the @font-face declarations
    # out of a fragment so the wrapper chrome renders in Titillium too.
    font_faces = "\n".join(re.findall(r"@font-face\{[^}]*\}", raw["compare"]))

    docs_js = "{" + ",".join(f'"{k}":{js_string(v)}' for k, v in docs.items()) + "}"

    tabs_html = "\n".join(
        f'      <button class="tab" role="tab" data-key="{k}" '
        f'aria-selected="{"true" if k == DEFAULT else "false"}">{FRAGMENTS[k][1]}</button>'
        for k in TAB_ORDER
    )

    return f"""<title>Green steel model — scenarios, network & workflow</title>
<style>
{font_faces}
  :root {{
    --bg:#ffffff; --nav-bg:#eef2f5; --fg:#1b2a33; --muted:#5b6b75;
    --accent:#0A5680; --on-accent:#ffffff; --border:rgba(10,86,128,.14);
    --hover:rgba(10,86,128,.08);
  }}
  @media (prefers-color-scheme:dark) {{
    :root {{
      --bg:#0b1418; --nav-bg:#0e1a20; --fg:#e6edf1; --muted:#8ea6b0;
      --accent:#91C096; --on-accent:#0E1A20; --border:rgba(145,192,150,.22);
      --hover:rgba(145,192,150,.12);
    }}
  }}
  :root[data-theme="light"] {{
    --bg:#ffffff; --nav-bg:#eef2f5; --fg:#1b2a33; --muted:#5b6b75;
    --accent:#0A5680; --on-accent:#ffffff; --border:rgba(10,86,128,.14);
    --hover:rgba(10,86,128,.08);
  }}
  :root[data-theme="dark"] {{
    --bg:#0b1418; --nav-bg:#0e1a20; --fg:#e6edf1; --muted:#8ea6b0;
    --accent:#91C096; --on-accent:#0E1A20; --border:rgba(145,192,150,.22);
    --hover:rgba(145,192,150,.12);
  }}
  html,body {{ height:100%; margin:0; }}
  .app {{
    display:flex; flex-direction:column; height:100vh; height:100dvh;
    background:var(--bg); color:var(--fg);
    font-family:'Titillium Web',-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  }}
  .bar {{
    display:flex; align-items:center; gap:18px; flex-wrap:wrap;
    padding:10px 16px; background:var(--nav-bg);
    border-bottom:1px solid var(--border); flex:0 0 auto;
  }}
  .brand {{
    font-weight:700; font-size:15px; letter-spacing:.02em; color:var(--fg);
    white-space:nowrap; display:flex; align-items:center; gap:8px;
  }}
  .brand .dot {{ width:10px; height:10px; border-radius:50%; background:var(--accent); }}
  .tabs {{ display:flex; gap:6px; flex-wrap:wrap; }}
  .tab {{
    appearance:none; border:1px solid transparent; border-radius:999px;
    padding:7px 15px; font:inherit; font-size:13.5px; font-weight:600;
    color:var(--muted); background:transparent; cursor:pointer;
    transition:background .12s ease,color .12s ease;
  }}
  .tab:hover {{ background:var(--hover); color:var(--fg); }}
  .tab[aria-selected="true"] {{ background:var(--accent); color:var(--on-accent); }}
  .tab:focus-visible {{ outline:2px solid var(--accent); outline-offset:2px; }}
  .stage {{ position:relative; flex:1 1 auto; min-height:0; background:var(--bg); }}
  iframe {{
    position:absolute; inset:0; width:100%; height:100%; border:0;
    background:var(--bg);
  }}
  /* work-in-progress banner — FCA house style: magenta-red (#D75674) accent with
     the same 4px left rule the dashboard's scenario cards use. */
  .wip {{ --wip:#D75674; --wip-bg:rgba(215,86,116,.10); --wip-line:rgba(215,86,116,.32); }}
  @media (prefers-color-scheme:dark) {{
    .wip {{ --wip:#E68199; --wip-bg:rgba(215,86,116,.16); --wip-line:rgba(215,86,116,.42); }}
  }}
  :root[data-theme="light"] .wip {{ --wip:#D75674; --wip-bg:rgba(215,86,116,.10); --wip-line:rgba(215,86,116,.32); }}
  :root[data-theme="dark"] .wip {{ --wip:#E68199; --wip-bg:rgba(215,86,116,.16); --wip-line:rgba(215,86,116,.42); }}
  .wip {{ position:relative; flex:0 0 auto; padding:10px 16px 10px 20px;
    background:var(--wip-bg); border-bottom:1px solid var(--wip-line); }}
  .wip::before {{ content:""; position:absolute; left:0; top:0; bottom:0; width:4px; background:var(--wip); }}
  .wip .wl {{ display:flex; align-items:center; gap:8px; font-size:16px; font-weight:800;
    letter-spacing:.1em; text-transform:uppercase; color:var(--wip); line-height:1.1; }}
  .wip .wi {{ margin-top:3px; font-size:11.5px; color:var(--muted); line-height:1.4; }}
</style>

<div class="app">
  <div class="wip" role="alert">
    <div class="wl"><span aria-hidden="true">&#9888;</span> Work in progress</div>
    <div class="wi">Based on unchecked / placeholder assumptions — the pipeline is not fully validated as of {build_date}.</div>
  </div>
  <nav class="bar" role="tablist" aria-label="Views">
    <span class="brand"><span class="dot"></span>Green steel model</span>
    <div class="tabs">
{tabs_html}
    </div>
  </nav>
  <div class="stage">
    <iframe id="view" title="View" allowfullscreen></iframe>
  </div>
</div>

<script>
  const DOCS = {docs_js};
  const view = document.getElementById('view');
  const tabs = Array.from(document.querySelectorAll('.tab'));
  let current = null;

  function curTheme() {{
    const dt = document.documentElement.getAttribute('data-theme');
    if (dt === 'dark' || dt === 'light') return dt;
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches
      ? 'dark' : 'light';
  }}

  function applyThemeToFrame() {{
    try {{
      const doc = view.contentDocument;
      if (doc && doc.documentElement) doc.documentElement.setAttribute('data-theme', curTheme());
    }} catch (e) {{ /* cross-origin should not happen for srcdoc */ }}
  }}
  view.addEventListener('load', applyThemeToFrame);

  function activate(key) {{
    if (key === current) return;
    current = key;
    tabs.forEach(t => t.setAttribute('aria-selected', String(t.dataset.key === key)));
    view.srcdoc = DOCS[key];        // (re)load; theme applied on load
  }}

  tabs.forEach(t => t.addEventListener('click', () => activate(t.dataset.key)));

  // Keep iframe in sync when the viewer toggles the artifact theme.
  new MutationObserver(applyThemeToFrame).observe(
    document.documentElement, {{ attributes: true, attributeFilter: ['data-theme'] }});
  if (window.matchMedia) {{
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    (mq.addEventListener ? mq.addEventListener.bind(mq, 'change')
                         : mq.addListener.bind(mq))(applyThemeToFrame);
  }}

  activate('{DEFAULT}');
</script>
"""


def main():
    html = build_hub()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(html, encoding="utf-8")
    print(f"wrote {OUT_PATH} ({OUT_PATH.stat().st_size/1e6:.2f} MB, {html.count(chr(0))} NULs)")


if __name__ == "__main__":
    main()
