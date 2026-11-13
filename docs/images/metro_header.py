#!/usr/bin/env python3
"""Replace nf-metro's one-line title with the CECILIA header, then move the map down to fit it.

  nf-metro render docs/images/cecilia_metro.mmd --embed-font -o docs/images/cecilia_metro.svg
  python3 docs/images/metro_header.py docs/images/cecilia_metro.svg nextflow.config

nf-metro draws the title as a single line at one size, so the header is added afterwards:

  CECILIA                                        large, bold, in the logo's colours
  usearCh basEd ampliCon pIpeLine for Illumina dAta   descriptor
  v<manifest.version> · Pipeline overview        version from nextflow.config

The logo has gold lettering on black and green microbes. Gold is unreadable on the map's light
background, so the brand text uses the logo's gold in dark mode and its green in light mode,
switched by the same light-dark() CSS the map uses. Run it once on a fresh render: it stops if
the SVG has no nf-metro title, which also means it was already applied.
"""
import json
import re
import sys

GOLD = "#dabf7a"   # logo lettering
GREEN = "#215d37"  # logo microbes
FONT = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"
SHIFT = 110        # extra height for the header; the map moves down by this much
X = 22             # left edge of the section boxes

CSS = (
    f".cecilia-brand {{ fill: light-dark({GREEN}, {GOLD}); }}\n"
    ".cecilia-descriptor { fill: var(--nfm-map-label-color, light-dark(#242424, #e0d9f7)); }\n"
    ".cecilia-version { fill: var(--nfm-map-section-label-color, light-dark(#6c757d, #9d94b8)); }\n"
)


def version(config):
    with open(config) as fh:
        text = fh.read()
    block = re.search(r"manifest\s*{(.*?)}", text, re.S)
    m = block and re.search(r"version\s*=\s*['\"]([^'\"]+)['\"]", block.group(1))
    if not m:
        raise SystemExit(f"ERROR: no manifest version in {config}")
    return m.group(1)


def header(ver):
    # Fill attributes match the dark palette, as nf-metro's own elements do; the CSS classes
    # choose the colour for the viewer's mode.
    return (
        '<g class="cecilia-header">\n'
        f'<text x="{X}" y="64" font-size="56" font-weight="bold" letter-spacing="5" fill="{GOLD}" '
        f'font-family="{FONT}" class="cecilia-brand">CECILIA</text>\n'
        f'<text x="{X + 2}" y="96" font-size="19" fill="#e0d9f7" font-family="{FONT}" '
        'class="cecilia-descriptor">Amplicon pipeline for Illumina data, built on USEARCH</text>\n'
        f'<text x="{X + 2}" y="123" font-size="14" letter-spacing="0.5" font-family="{FONT}">'
        f'<tspan fill="#9d94b8" class="cecilia-version">v{ver}</tspan>'
        '<tspan fill="#9d94b8" class="cecilia-version"> · </tspan>'
        f'<tspan font-weight="bold" fill="{GOLD}" class="cecilia-brand">Pipeline overview</tspan></text>\n'
        "</g>\n"
    )


def one(pattern, text, what, flags=0):
    found = list(re.finditer(pattern, text, flags))
    if len(found) != 1:
        raise SystemExit(f"ERROR: expected one {what}, found {len(found)}")
    return found[0]


def main(svg_path, config):
    with open(svg_path, encoding="utf-8") as fh:
        svg = fh.read()
    ver = version(config)

    title = one(r'<text [^>]*class="nf-metro-title"[^>]*>.*?</text>\n?', svg, "nf-metro title")
    svg = svg[:title.start()] + svg[title.end():]

    root = one(r'<svg [^>]*width="([\d.]+)" height="([\d.]+)" viewBox="0 0 ([\d.]+) ([\d.]+)"', svg, "root size")
    w, h = root.group(1), float(root.group(2))
    new_h = f"{h + SHIFT:g}"
    svg = (svg[:root.start()]
           + root.group(0).replace(f'height="{root.group(2)}"', f'height="{new_h}"')
                          .replace(f'0 0 {root.group(3)} {root.group(4)}', f"0 0 {root.group(3)} {new_h}")
           + svg[root.end():])

    bg = one(r'<rect x="0" y="0" width="[\d.]+" height="[\d.]+"[^>]*class="nf-metro-bg" />\n?', svg, "background")
    bg_tag = re.sub(r'height="[\d.]+"', f'height="{new_h}"', bg.group(0))
    end = svg.rindex("</svg>")
    svg = (svg[:bg.start()] + bg_tag + header(ver)
           + f'<g transform="translate(0,{SHIFT})">\n' + svg[bg.end():end] + "</g>\n" + svg[end:])

    style_end = one(r"</style>", svg, "style block").start()
    svg = svg[:style_end] + CSS + svg[style_end:]

    # Keep nf-metro's embedded node coordinates in step with the drawing.
    meta = one(r'(<metadata id="diagram-manifest"[^>]*><!\[CDATA\[)(.*?)(\]\]></metadata>)', svg, "manifest", re.S)
    data = json.loads(meta.group(2))
    data["height"] = h + SHIFT
    data["title"] = "CECILIA"
    for node in data.get("nodes", []):
        node["y"] += SHIFT
    payload = json.dumps(data, separators=(",", ":"), sort_keys=True)
    svg = svg[:meta.start()] + meta.group(1) + payload + meta.group(3) + svg[meta.end():]

    with open(svg_path, "w", encoding="utf-8") as fh:
        fh.write(svg)
    print(f"{svg_path}: header added (v{ver}), {w} x {new_h}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
