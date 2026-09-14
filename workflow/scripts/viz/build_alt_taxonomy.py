"""Build the single-taxonomy cost-breakdown page.

The same page as build_cost_taxonomies, showing the alternative taxonomy alone
rather than beside the one the dashboard plots. One template serves both: the
`ALT_ONLY` flag hides the second column and skips its plots, so the two cannot
drift apart.

Output: results/html/alt_taxonomy.html — body-only, for the hub.
"""
from build_cost_taxonomies import TEMPLATE_HTML, attach
from build_dashboard import HTML_DIR, build_html

OUT_PATH = HTML_DIR / "alt_taxonomy.html"


def main() -> None:
    html, cases, geos = build_html(TEMPLATE_HTML, augment=attach)
    html = (html
            .replace("/*ALT_ONLY*/false", "/*ALT_ONLY*/true")
            .replace("/*ALT_ONLY_CLASS*/", "alt-only"))
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(html, encoding="utf-8")
    print(f"wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1e6:.2f} MB) — "
          f"{len(cases)} projects, {len(geos)} geos, body-only hub fragment")


if __name__ == "__main__":
    main()
