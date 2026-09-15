"""Build the v2 scenario-comparison dashboard.

Reuses build_dashboard's payload assembly verbatim (same cases/synth/gas +
label/colour maps) and only swaps in the v2 template, which is scenario-centric:
a full button matrix picks Scenario A, a compact override strip derives Scenario
B from it, and the cost-breakdown / capacity charts show A next to B. The
levelised-cost table sits below the charts and shows each scenario's route group.
"""
from pathlib import Path

from build_dashboard import HTML_DIR, build_html

TEMPLATE_V2 = Path(__file__).with_name("dashboard_v2_template.html")
OUT_V2 = HTML_DIR / "dashboard_v2.html"


def main():
    html, cases, geos = build_html(TEMPLATE_V2)
    OUT_V2.parent.mkdir(parents=True, exist_ok=True)
    OUT_V2.write_text(html)
    print(f"wrote {OUT_V2} ({OUT_V2.stat().st_size/1e6:.2f} MB) — {len(cases)} projects, {len(geos)} geos")


if __name__ == "__main__":
    main()
