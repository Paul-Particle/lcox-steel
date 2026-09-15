"""Build the Cost breakdown page — the hub's first tab.

The same page as build_cost_taxonomies, showing the by-purpose split alone
rather than beside the reported cost groups. One template serves both: the
`PURPOSE_ONLY` flag hides the second column and skips its plots, so the two
cannot drift apart.

Output: results/html/cost_breakdown.html — body-only, for the hub.
"""
from build_cost_taxonomies import TEMPLATE_HTML, attach
from build_dashboard import HTML_DIR, build_html

OUT_PATH = HTML_DIR / "cost_breakdown.html"


def main() -> None:
    html, cases, geos = build_html(TEMPLATE_HTML, augment=attach)
    html = (html
            .replace("/*PURPOSE_ONLY*/false", "/*PURPOSE_ONLY*/true")
            .replace("/*PURPOSE_ONLY_CLASS*/", "purpose-only"))
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(html, encoding="utf-8")
    print(f"wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1e6:.2f} MB) — "
          f"{len(cases)} projects, {len(geos)} geos, body-only hub fragment")


if __name__ == "__main__":
    main()
