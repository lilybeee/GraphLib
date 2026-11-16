/* graph_stream.c — see graph_stream.h. */

#include "graph_stream.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define LINE_MAX_BYTES 512

/* Reads the next line into buf; false at end of file.  A line longer than the buffer
 * is reported through *too_long (edge-list lines are short, so that means a bad file). */
static bool next_line(FILE *in, char *buf, bool *too_long)
{
    *too_long = false;
    if (fgets(buf, LINE_MAX_BYTES, in) == NULL) {
        return false;
    }
    const size_t len = strlen(buf);
    if (len == LINE_MAX_BYTES - 1 && buf[len - 1] != '\n') {
        *too_long = true;
    }
    return true;
}

static const char *skip_space(const char *p)
{
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') {
        p++;
    }
    return p;
}

static bool parse_id(const char *p, const char **end, int *id)
{
    char *e;
    errno = 0;
    const long long v = strtoll(p, &e, 10);
    if (e == p || errno != 0 || v < 0 || v >= INT_MAX) {
        return false;
    }
    *end = e;
    *id  = (int)v;
    return true;
}

/*
 * Walks the file.  For each header it calls on_header(kind, value), for each edge
 * on_edge(u, v).  Returns GRAPH_OK or the first error.
 */
typedef enum { H_DIRECTED, H_UNDIRECTED, H_VERTICES } HeaderKind;
typedef GraphStatus (*HeaderFn)(void *ctx, HeaderKind kind, size_t value);
typedef GraphStatus (*EdgeFn)(void *ctx, int u, int v);

static GraphStatus walk(FILE *in, HeaderFn on_header, EdgeFn on_edge, void *ctx)
{
    if (in == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    char buf[LINE_MAX_BYTES];
    bool too_long;
    while (next_line(in, buf, &too_long)) {
        if (too_long) {
            return GRAPH_ERR_IO;
        }
        const char *p = skip_space(buf);
        if (*p == '\0') {
            continue;
        }
        if (*p == '#') {
            p = skip_space(p + 1);
            GraphStatus st = GRAPH_OK;
            if (strncmp(p, "directed", 8) == 0) {
                st = on_header(ctx, H_DIRECTED, 0);
            } else if (strncmp(p, "undirected", 10) == 0) {
                st = on_header(ctx, H_UNDIRECTED, 0);
            } else if (strncmp(p, "vertices", 8) == 0) {
                char *end;
                errno = 0;
                const unsigned long long n = strtoull(p + 8, &end, 10);
                if (end != p + 8 && errno == 0 && n < INT_MAX) {
                    st = on_header(ctx, H_VERTICES, (size_t)n);
                }
            }
            if (st != GRAPH_OK) {
                return st;
            }
            continue;
        }
        int u, v;
        const char *end;
        if (!parse_id(p, &end, &u) || !parse_id(skip_space(end), &end, &v)) {
            return GRAPH_ERR_IO;
        }
        const GraphStatus st = on_edge(ctx, u, v);
        if (st != GRAPH_OK) {
            return st;
        }
    }
    return ferror(in) ? GRAPH_ERR_IO : GRAPH_OK;
}

/* ---- generic edge walk ---- */

static GraphStatus ignore_header(void *ctx, HeaderKind kind, size_t value)
{
    (void)ctx; (void)kind; (void)value;
    return GRAPH_OK;
}

typedef struct {
    GraphStreamEdgeFn fn;
    void             *ctx;
} EdgeCall;

static GraphStatus forward_edge(void *ctx, int u, int v)
{
    const EdgeCall *c = ctx;
    return c->fn(c->ctx, u, v);
}

GraphStatus graph_stream_edges(FILE *in, GraphStreamEdgeFn fn, void *ctx)
{
    if (in == NULL || fn == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    EdgeCall call = {fn, ctx};
    const GraphStatus st = walk(in, ignore_header, forward_edge, &call);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

/* ---- scan ---- */

static GraphStatus scan_header(void *ctx, HeaderKind kind, size_t value)
{
    GraphStreamInfo *info = ctx;
    if (kind == H_DIRECTED)   info->directed = true;
    if (kind == H_UNDIRECTED) info->directed = false;
    if (kind == H_VERTICES && value > info->vertices) info->vertices = value;
    return GRAPH_OK;
}

static GraphStatus scan_edge(void *ctx, int u, int v)
{
    GraphStreamInfo *info = ctx;
    const size_t need = (size_t)(u > v ? u : v) + 1;
    if (need > info->vertices) info->vertices = need;
    info->edges++;
    return GRAPH_OK;
}

GraphStatus graph_stream_scan(FILE *in, GraphStreamInfo *info)
{
    if (info == NULL) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    *info = (GraphStreamInfo){0};
    const GraphStatus st = walk(in, scan_header, scan_edge, info);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

/* ---- degrees ---- */

typedef struct {
    size_t  n;
    bool    directed;
    size_t *out;
    size_t *in;
} DegCtx;

static GraphStatus deg_header(void *ctx, HeaderKind kind, size_t value)
{
    (void)ctx; (void)kind; (void)value;
    return GRAPH_OK;
}

static GraphStatus deg_edge(void *ctx, int u, int v)
{
    DegCtx *c = ctx;
    if ((size_t)u >= c->n || (size_t)v >= c->n) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    c->out[u]++;
    if (c->directed) {
        if (c->in != NULL) c->in[v]++;
    } else if (u != v) {
        c->out[v]++;                           /* an undirected edge touches both ends */
    }
    return GRAPH_OK;
}

GraphStatus graph_stream_degrees(FILE *in, size_t n, bool directed, size_t *out_degree,
                                 size_t *in_degree)
{
    if (in == NULL || (n > 0 && out_degree == NULL)) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    memset(out_degree, 0, n * sizeof *out_degree);
    if (in_degree != NULL) {
        memset(in_degree, 0, n * sizeof *in_degree);
    }
    DegCtx c = {n, directed, out_degree, in_degree};
    const GraphStatus st = walk(in, deg_header, deg_edge, &c);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}

/* ---- components ---- */

typedef struct {
    size_t n;
    int   *parent;
    int   *size;
} UfCtx;

static int uf_find(UfCtx *c, int x)
{
    while (c->parent[x] != x) {
        c->parent[x] = c->parent[c->parent[x]];     /* path halving */
        x = c->parent[x];
    }
    return x;
}

static GraphStatus uf_header(void *ctx, HeaderKind kind, size_t value)
{
    (void)ctx; (void)kind; (void)value;
    return GRAPH_OK;
}

static GraphStatus uf_edge(void *ctx, int u, int v)
{
    UfCtx *c = ctx;
    if ((size_t)u >= c->n || (size_t)v >= c->n) {
        return GRAPH_ERR_INVALID_VERTEX;
    }
    int a = uf_find(c, u), b = uf_find(c, v);
    if (a != b) {
        if (c->size[a] < c->size[b]) { int t = a; a = b; b = t; }
        c->parent[b] = a;                            /* union by size keeps the trees shallow */
        c->size[a] += c->size[b];
    }
    return GRAPH_OK;
}

GraphStatus graph_stream_components(FILE *in, size_t n, int *component, size_t *count)
{
    if (in == NULL || (n > 0 && component == NULL) || n >= INT_MAX) {
        return graph_set_error(GRAPH_ERR_INVALID_ARGUMENT);
    }
    UfCtx c = {n, malloc((n > 0 ? n : 1) * sizeof(int)), malloc((n > 0 ? n : 1) * sizeof(int))};
    if (c.parent == NULL || c.size == NULL) {
        free(c.parent);
        free(c.size);
        return graph_set_error(GRAPH_ERR_OUT_OF_MEMORY);
    }
    for (size_t i = 0; i < n; i++) {
        c.parent[i] = (int)i;
        c.size[i]   = 1;
    }
    GraphStatus st = walk(in, uf_header, uf_edge, &c);
    if (st == GRAPH_OK) {
        /* Number the components 0.. by first member; c.size is reused as the label table. */
        for (size_t i = 0; i < n; i++) c.size[i] = -1;
        size_t labels = 0;
        for (size_t i = 0; i < n; i++) {
            const int root = uf_find(&c, (int)i);
            if (c.size[root] < 0) c.size[root] = (int)labels++;
            component[i] = c.size[root];
        }
        if (count != NULL) *count = labels;
    }
    free(c.parent);
    free(c.size);
    return st == GRAPH_OK ? GRAPH_OK : graph_set_error(st);
}
