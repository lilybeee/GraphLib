/*
 * graph_viz.h — draw a Graph as an image.
 *
 *   graph_visualize(g, "out.png", NULL);     // format chosen from the extension
 *
 * Works on any size and either representation (it only uses the public Graph
 * API).  Three things scale with the graph:
 *
 *   Layout   force-directed (Fruchterman-Reingold) where each vertex is only
 *            repelled by vertices in nearby grid cells, so one iteration costs
 *            about O(V + E) instead of O(V^2).  Handles 100k+ vertices.
 *   Format   SVG (crisp, zoomable, optional labels) up to a few thousand
 *            vertices; a raster PNG (written by our own encoder, no libraries)
 *            for anything bigger.  VIZ_FORMAT_AUTO picks by size.
 *   Style    edges get fainter as there are more of them and overlap into
 *            darker "density", vertex size follows degree, colour follows
 *            component, community or degree.  VIZ_LAYOUT_MATRIX instead draws
 *            the adjacency matrix as a density heat-map, which stays readable
 *            when a node-link drawing would be a solid blob.
 *
 * Reading the picture
 *   Colour by group (VIZ_COLOR_COMMUNITY / _COMPONENT): every group gets a colour, a name tag
 *     ("Community 3") at the centre of its members, and a legend entry with its size and its
 *     most-connected member.  Only the 11 largest groups are coloured; the rest are grey.
 *   Colour by degree (VIZ_COLOR_DEGREE): a temperature scale from blue (few connections) to dark
 *     red (very popular) on a log scale, with a colour bar and a list of the most popular people.
 *     Popular vertices are drawn last, larger and with a glow.
 *   The most-connected people are named on the picture (label_top).  Text in PNG output uses a
 *   built-in bitmap font, so no libraries are needed.
 *
 * Arrowheads: SVG draws them for directed graphs; the PNG only for directed graphs of <= 300 vertices.
 *
 * Failures return a GraphStatus (GRAPH_ERR_INVALID_ARGUMENT, _IO,
 * _OUT_OF_MEMORY, _TOO_LARGE).  Drawing does not change the graph, but like
 * any other call on a non-frozen Graph it counts as an operation for the
 * adaptive observer.
 */

#ifndef GRAPH_VIZ_H
#define GRAPH_VIZ_H

#include <stdbool.h>
#include <stdio.h>

#include "graph.h"

typedef enum {
    VIZ_FORMAT_AUTO,     /* from the file extension (.svg / .png), else by graph size */
    VIZ_FORMAT_SVG,
    VIZ_FORMAT_PNG
} VizFormat;

typedef enum {
    VIZ_LAYOUT_FORCE,    /* node-link drawing with a force-directed layout (default) */
    VIZ_LAYOUT_MATRIX    /* adjacency-matrix heat-map, vertices in ID order          */
} VizLayout;

typedef enum {
    VIZ_COLOR_AUTO,      /* component colours if the graph has several, else community */
    VIZ_COLOR_COMPONENT,
    VIZ_COLOR_COMMUNITY, /* Louvain (falls back to degree if it cannot run, or above 50000 vertices) */
    VIZ_COLOR_DEGREE,    /* blue (low) to yellow (high)                                */
    VIZ_COLOR_NONE
} VizColor;

typedef struct {
    VizFormat format;
    VizLayout layout;
    VizColor  color;
    int       width, height;     /* pixels / SVG units, 64..8192; default 1200 x 1200      */
    int       iterations;        /* layout iterations; 0 = choose from the graph size      */
    unsigned long long seed;     /* layout start; same seed gives the same picture         */
    int       show_labels;       /* -1 = auto (SVG, <= 150 vertices), 0 = never, 1 = always */
    bool      show_weights;      /* SVG only, auto-disabled above 150 edges                  */
    const char *const *labels;   /* optional: labels[id] for each ID; default is the ID      */
    const char *title;           /* optional caption                                         */
    bool      legend;            /* legend panel and community name tags (group colours), or the
                                    temperature bar; default true                                  */
    int       label_top;         /* name this many most-connected vertices: -1 = auto (12; none when
                                    every vertex already has a label or show_labels is 0), 0 = none                     */
    const char *group_name;      /* what a colour group is called in the legend ("Community"); default "Group" */
} VizOptions;

/** 1200x1200, force layout, auto format and colour, auto labels, seed 1. */
VizOptions graph_viz_default_options(void);

/**
 * graph_visualize - draw @g into the file @path (created or overwritten).
 * @opt may be NULL for the defaults.  With VIZ_FORMAT_AUTO the extension wins:
 * ".svg" or ".png"; anything else is chosen by size (SVG up to 1500 vertices and
 * 8000 edges, else PNG).
 */
GraphStatus graph_visualize(Graph *g, const char *path, const VizOptions *opt);

/** Same, to an open stream (VIZ_FORMAT_AUTO then decides by size only). */
GraphStatus graph_visualize_to(Graph *g, FILE *out, const VizOptions *opt);

/**
 * graph_layout - just the force-directed layout.  @x and @y need
 * graph_id_bound() entries; vertices are placed in [0, 1] x [0, 1], dead IDs get
 * NAN.  Same @opt fields as above (iterations, seed).  Use it to draw the graph
 * yourself.  Complexity: about O(iterations * (V + E)).
 */
GraphStatus graph_layout(Graph *g, double *x, double *y, const VizOptions *opt);

#endif /* GRAPH_VIZ_H */
