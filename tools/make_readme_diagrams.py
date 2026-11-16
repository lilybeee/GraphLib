#!/usr/bin/env python3
"""
Draws the two diagrams in the README that the library cannot draw itself:

  docs/images/architecture.svg   how a call travels through the library
  docs/images/scale_1M.svg       time and peak memory of each algorithm at 1,000,000 vertices

The scale chart is built from the output of the scale suite:

    ./build/scale_test 1000000 0 > scale.txt        (list backend, 1 000 000 vertices)
    python3 tools/make_readme_diagrams.py scale.txt

Without an argument only the architecture diagram is written.
"""
import re
import sys

BLUE, GREEN, ORANGE, GREY, RED = "#dbeafe", "#dcfce7", "#ffedd5", "#f3f4f6", "#fee2e2"


def box(x, y, w, h, title, lines=(), fill=GREY, title_size=15):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{fill}" stroke="#374151" stroke-width="1.4"/>',
           f'<text x="{x + w / 2}" y="{y + 23}" text-anchor="middle" font-size="{title_size}" font-weight="bold" fill="#111827">{title}</text>']
    for i, line in enumerate(lines):
        out.append(f'<text x="{x + w / 2}" y="{y + 43 + 16 * i}" text-anchor="middle" font-size="12" fill="#374151">{line}</text>')
    return "\n".join(out)


def arrow(x1, y1, x2, y2, label=""):
    s = (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#374151" stroke-width="1.6" marker-end="url(#a)"/>')
    if label:
        s += f'<text x="{(x1 + x2) / 2 + 6}" y="{(y1 + y2) / 2}" font-size="11" fill="#6b7280">{label}</text>'
    return s


def architecture():
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="980" height="640" viewBox="0 0 980 640" font-family="sans-serif">',
             '<defs><marker id="a" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto">'
             '<path d="M0,0 L10,4 L0,8 z" fill="#374151"/></marker></defs>',
             '<rect width="100%" height="100%" fill="#ffffff"/>',
             '<text x="490" y="30" text-anchor="middle" font-size="20" font-weight="bold" fill="#111827">How a call travels through the library</text>']
    # user code
    parts.append(box(340, 50, 300, 62, "Your program", ["graph_add_edge(g, u, v, w)   graph_pagerank(g, ...)"], BLUE))
    parts.append(arrow(490, 112, 490, 150))
    # facade
    parts.append(box(190, 150, 600, 96, "Graph handle  (graph.h)",
                     ["one API; vertices are integers; every call returns a status or records graph_last_error()",
                      "forwards to whichever representation is active, then tells the observer what happened"], GREEN))
    # side boxes
    parts.append(box(20, 150, 150, 96, "Cache", ["versioned LRU of", "baseline results", "(BFS, Dijkstra, ...)"], GREY))
    parts.append(box(810, 150, 150, 96, "Memory budget", ["4 GiB default;", "too big = clean", "GRAPH_ERR_TOO_LARGE"], RED))
    parts.append(f'<line x1="170" y1="198" x2="190" y2="198" stroke="#374151" stroke-width="1.4"/>')
    parts.append(f'<line x1="790" y1="198" x2="810" y2="198" stroke="#374151" stroke-width="1.4"/>')
    # observer
    parts.append(arrow(340, 246, 250, 300, "counts operations"))
    parts.append(box(60, 300, 380, 110, "Observer  (policy.h, stats.h)",
                     ["workload counters (decayed) + graph shape",
                      "+ each backend's cost model decide it:",
                      "migrate when the gain is clear, the saving repays",
                      "the rebuild, and the target fits budget and size limit"], ORANGE))
    parts.append(arrow(440, 355, 560, 355, "migrate (migrate.h)"))
    # backends
    parts.append(arrow(640, 246, 640, 300, "GraphBackend vtable"))
    parts.append(box(560, 300, 190, 110, "List backend", ["sorted neighbour arrays", "~47 B / edge (23 unweighted)", "any size; grows with edges"], BLUE))
    parts.append(box(770, 300, 190, 110, "Matrix backend", ["bit rows + sparse weight tiles", "fast has_edge, dense graphs", "size-limited: ~131k / 176k"], BLUE))
    # algorithm layers
    parts.append(f'<text x="490" y="450" text-anchor="middle" font-size="15" font-weight="bold" fill="#111827">What you can call on a Graph</text>')
    cols = [("graph_algorithms.h", ["~90 algorithms:", "paths, flow, matching,", "centrality, communities,", "colouring, cliques, I/O"]),
            ("graph_advanced.h", ["graph_adv_*:", "isomorphism, Johnson,", "IDA*, JPS, widest path,", "dynamic connectivity, GNN"]),
            ("graph_social.h", ["per-user queries:", "within k hops,", "people you may know,", "bidirectional path"]),
            ("graph_stream.h", ["edge-list FILE analyses:", "scan, degrees,", "components in O(V)", "memory (beyond RAM)"]),
            ("graph_viz.h", ["pictures of any size:", "SVG / PNG,", "communities, degree,", "matrix heat-map"])]
    for i, (t, ls) in enumerate(cols):
        x = 20 + i * 190
        parts.append(box(x, 464, 178, 108, t, ls, GREEN, 14))
    parts.append(f'<text x="490" y="600" text-anchor="middle" font-size="12" fill="#6b7280">'
                 f'Each of these is built on the Graph handle, so it works the same on either representation.</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def scale_chart(path):
    rows = []
    for line in open(path):
        m = re.match(r"^(.*?)\s+list\s+1000000\s+(\d+) ms\s+(\d+) MB\s+ok", line)
        if m:
            rows.append((m.group(1).strip(), int(m.group(2)), int(m.group(3))))
    rows.sort(key=lambda r: -r[1])
    rows = rows[:24]
    w, row_h, left = 980, 24, 250
    h = 90 + row_h * len(rows)
    biggest = max(r[1] for r in rows)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="sans-serif">',
           '<rect width="100%" height="100%" fill="#ffffff"/>',
           f'<text x="{w / 2}" y="28" text-anchor="middle" font-size="19" font-weight="bold" fill="#111827">'
           'One million users: time and peak memory per algorithm (list backend)</text>',
           f'<text x="{w / 2}" y="48" text-anchor="middle" font-size="12" fill="#6b7280">'
           'about 4 million friendships, 4 GiB cap; each case checked against an independent answer; the whole run took 3.5 minutes</text>']
    for i, (name, ms, mb) in enumerate(rows):
        y = 70 + i * row_h
        bar = max(2, (w - left - 190) * ms / biggest)
        out.append(f'<text x="{left - 10}" y="{y + 15}" text-anchor="end" font-size="12" fill="#111827">{name}</text>')
        out.append(f'<rect x="{left}" y="{y + 2}" width="{bar:.1f}" height="{row_h - 8}" rx="3" fill="#3b82f6"/>')
        out.append(f'<text x="{left + bar + 8:.1f}" y="{y + 15}" font-size="12" fill="#374151">{ms / 1000:.1f} s   ({mb} MB)</text>')
    out.append('</svg>')
    return "\n".join(out)


if __name__ == "__main__":
    open("docs/images/architecture.svg", "w").write(architecture())
    print("wrote docs/images/architecture.svg")
    if len(sys.argv) > 1:
        open("docs/images/scale_1M.svg", "w").write(scale_chart(sys.argv[1]))
        print("wrote docs/images/scale_1M.svg")
