/*
 * matrix_io.c — reading and writing graphs.  See matrix_algorithms.h.
 */

#include <ctype.h>
#include <errno.h>
#include <limits.h>

#include "matrix_util.h"

/* ------------------------------------------------------------------ */
/*  Writers                                                             */
/* ------------------------------------------------------------------ */

/** Print s with `"` and `\` escaped, for DOT and JSON string literals. */
static void put_escaped(FILE *out, const char *s)
{
    for (; *s != '\0'; ++s) {
        if (*s == '"' || *s == '\\') {
            fputc('\\', out);
        }
        if (*s == '\n') {
            fputs("\\n", out);
        } else {
            fputc(*s, out);
        }
    }
}

MatrixStatus matrix_write_dot(const MatrixGraph *g, FILE *out, const char *name,
                              const char *const *labels, bool show_weights)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    fprintf(out, "%s \"", directed ? "digraph" : "graph");
    put_escaped(out, name != NULL ? name : "G");
    fputs("\" {\n", out);
    for (size_t v = 0; v < n; ++v) {
        if (!matrix_graph_is_alive(g, (int)v)) {
            continue;
        }
        if (labels != NULL && labels[v] != NULL) {
            fprintf(out, "  %zu [label=\"", v);
            put_escaped(out, labels[v]);
            fputs("\"];\n", out);
        } else {
            fprintf(out, "  %zu;\n", v);
        }
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (!directed && v < (int)u) {
                continue;
            }
            fprintf(out, "  %zu %s %d", u, directed ? "->" : "--", v);
            if (show_weights) {
                fprintf(out, " [label=\"%g\"]", w);
            }
            fputs(";\n", out);
        }
    }
    fputs("}\n", out);
    return ferror(out) ? GRAPH_ERR_IO : GRAPH_OK;
}

MatrixStatus matrix_write_edge_list(const MatrixGraph *g, FILE *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    fprintf(out, "# %s\n# vertices %zu\n", directed ? "directed" : "undirected", n);
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (directed || v >= (int)u) {
                fprintf(out, "%zu %d %.17g\n", u, v, w);
            }
        }
    }
    return ferror(out) ? GRAPH_ERR_IO : GRAPH_OK;
}

MatrixStatus matrix_write_graphml(const MatrixGraph *g, FILE *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
          "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">\n"
          "  <key id=\"w\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>\n",
          out);
    fprintf(out, "  <graph id=\"G\" edgedefault=\"%s\">\n", directed ? "directed" : "undirected");
    for (size_t v = 0; v < n; ++v) {
        if (matrix_graph_is_alive(g, (int)v)) {
            fprintf(out, "    <node id=\"n%zu\"/>\n", v);
        }
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (directed || v >= (int)u) {
                fprintf(out, "    <edge source=\"n%zu\" target=\"n%d\"><data key=\"w\">%.17g</data></edge>\n",
                        u, v, w);
            }
        }
    }
    fputs("  </graph>\n</graphml>\n", out);
    return ferror(out) ? GRAPH_ERR_IO : GRAPH_OK;
}

MatrixStatus matrix_write_gml(const MatrixGraph *g, FILE *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    fprintf(out, "graph [\n  directed %d\n", directed ? 1 : 0);
    for (size_t v = 0; v < n; ++v) {
        if (matrix_graph_is_alive(g, (int)v)) {
            fprintf(out, "  node [ id %zu ]\n", v);
        }
    }
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (directed || v >= (int)u) {
                fprintf(out, "  edge [ source %zu target %d weight %.17g ]\n", u, v, w);
            }
        }
    }
    fputs("]\n", out);
    return ferror(out) ? GRAPH_ERR_IO : GRAPH_OK;
}

MatrixStatus matrix_write_json(const MatrixGraph *g, FILE *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    fprintf(out, "{\n  \"directed\": %s,\n  \"vertices\": [", directed ? "true" : "false");
    bool first = true;
    for (size_t v = 0; v < n; ++v) {
        if (matrix_graph_is_alive(g, (int)v)) {
            fprintf(out, "%s%zu", first ? "" : ", ", v);
            first = false;
        }
    }
    fputs("],\n  \"edges\": [", out);
    first = true;
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (directed || v >= (int)u) {
                fprintf(out, "%s\n    {\"source\": %zu, \"target\": %d, \"weight\": %.17g}",
                        first ? "" : ",", u, v, w);
                first = false;
            }
        }
    }
    fputs(first ? "]\n}\n" : "\n  ]\n}\n", out);
    return ferror(out) ? GRAPH_ERR_IO : GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Edge-list reading                                                   */
/* ------------------------------------------------------------------ */

/** Reads one line of any length into *buf (grown as needed); false at EOF. */
static bool read_line(FILE *in, char **buf, size_t *cap)
{
    size_t len = 0;
    int c;
    while ((c = fgetc(in)) != EOF) {
        if (len + 2 > *cap) {
            const size_t ncap = *cap > 0 ? *cap * 2 : 128;
            char *nb = realloc(*buf, ncap);
            if (nb == NULL) {
                return false;
            }
            *buf = nb;
            *cap = ncap;
        }
        (*buf)[len++] = (char)c;
        if (c == '\n') {
            break;
        }
    }
    if (len == 0) {
        return false;
    }
    (*buf)[len] = '\0';
    return true;
}

typedef struct {
    int    *u;
    int    *v;
    double *w;
    size_t  count;
    size_t  cap;
} EdgeBuf;

static bool edgebuf_push(EdgeBuf *b, int u, int v, double w)
{
    if (b->count == b->cap) {
        const size_t cap = b->cap > 0 ? b->cap * 2 : 64;
        int    *nu = realloc(b->u, cap * sizeof(int));
        if (nu == NULL) {
            return false;
        }
        b->u = nu;
        int    *nv = realloc(b->v, cap * sizeof(int));
        if (nv == NULL) {
            return false;
        }
        b->v = nv;
        double *nw = realloc(b->w, cap * sizeof(double));
        if (nw == NULL) {
            return false;
        }
        b->w = nw;
        b->cap = cap;
    }
    b->u[b->count] = u;
    b->v[b->count] = v;
    b->w[b->count] = w;
    ++b->count;
    return true;
}

static void edgebuf_free(EdgeBuf *b)
{
    free(b->u);
    free(b->v);
    free(b->w);
    memset(b, 0, sizeof(*b));
}

/** Build the graph from collected edges; duplicate edges keep their first weight. */
static MatrixStatus build_graph(const EdgeBuf *e, size_t vertices, bool directed, MatrixGraph **out)
{
    MatrixGraph *g = matrix_graph_create(vertices, directed);
    if (g == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < e->count; ++i) {
        if (!matrix_graph_has_edge(g, e->u[i], e->v[i], NULL) &&
            !matrix_graph_add_edge(g, e->u[i], e->v[i], e->w[i])) {
            matrix_graph_destroy(g);
            return GRAPH_ERR_IO;      /* non-finite weight */
        }
    }
    *out = g;
    return GRAPH_OK;
}

static char *skip_space(char *s)
{
    while (*s != '\0' && isspace((unsigned char)*s)) {
        ++s;
    }
    return s;
}

/** Parse a non-negative vertex ID that fits in an int. */
static bool parse_id(const char *s, char **end, int *id)
{
    errno = 0;
    const long v = strtol(s, end, 10);
    if (*end == s || errno != 0 || v < 0 || v >= INT_MAX) {
        return false;
    }
    *id = (int)v;
    return true;
}

MatrixStatus matrix_read_edge_list(FILE *in, bool directed, MatrixGraph **out)
{
    if (in == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    char *line = NULL;
    size_t cap = 0;
    EdgeBuf edges;
    memset(&edges, 0, sizeof(edges));
    size_t vertices = 0;
    MatrixStatus st = GRAPH_OK;

    while (st == GRAPH_OK && read_line(in, &line, &cap)) {
        char *p = skip_space(line);
        if (*p == '\0') {
            continue;
        }
        if (*p == '#') {
            p = skip_space(p + 1);
            if (strncmp(p, "directed", 8) == 0) {
                directed = true;
            } else if (strncmp(p, "undirected", 10) == 0) {
                directed = false;
            } else if (strncmp(p, "vertices", 8) == 0) {
                char *end;
                errno = 0;
                const unsigned long long n = strtoull(p + 8, &end, 10);
                if (end != p + 8 && errno == 0 && n < INT_MAX && n > vertices) {
                    vertices = (size_t)n;
                }
            }
            continue;
        }
        int u, v;
        char *end;
        double w = 1.0;
        if (!parse_id(p, &end, &u)) {
            st = GRAPH_ERR_IO;
            break;
        }
        p = skip_space(end);
        if (!parse_id(p, &end, &v)) {
            st = GRAPH_ERR_IO;
            break;
        }
        p = skip_space(end);
        if (*p != '\0') {
            w = strtod(p, &end);
            if (end == p || *skip_space(end) != '\0' || !isfinite(w)) {
                st = GRAPH_ERR_IO;
                break;
            }
        }
        if (!edgebuf_push(&edges, u, v, w)) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        const size_t need = (size_t)(u > v ? u : v) + 1;
        if (need > vertices) {
            vertices = need;
        }
    }
    free(line);
    if (st == GRAPH_OK && ferror(in)) {
        st = GRAPH_ERR_IO;
    }
    if (st == GRAPH_OK) {
        st = build_graph(&edges, vertices, directed, out);
    }
    edgebuf_free(&edges);
    return st;
}

/* ------------------------------------------------------------------ */
/*  Named edge lists                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    char  **names;
    size_t  count;
    size_t  cap;
    int    *slots;       /* open-addressing table of name indices, -1 = empty */
    size_t  nslots;
} NameTable;

static size_t hash_name(const char *s)
{
    size_t h = 5381;
    for (; *s != '\0'; ++s) {
        h = h * 33 + (unsigned char)*s;
    }
    return h;
}

static void names_free(NameTable *t)
{
    for (size_t i = 0; i < t->count; ++i) {
        free(t->names[i]);
    }
    free(t->names);
    free(t->slots);
    memset(t, 0, sizeof(*t));
}

static bool names_rehash(NameTable *t, size_t nslots)
{
    int *slots = malloc(nslots * sizeof(int));
    if (slots == NULL) {
        return false;
    }
    mu_fill_int(slots, nslots, -1);
    for (size_t i = 0; i < t->count; ++i) {
        size_t h = hash_name(t->names[i]) & (nslots - 1);
        while (slots[h] >= 0) {
            h = (h + 1) & (nslots - 1);
        }
        slots[h] = (int)i;
    }
    free(t->slots);
    t->slots  = slots;
    t->nslots = nslots;
    return true;
}

/** ID of `name`, adding it if new; -1 on allocation failure. */
static int names_intern(NameTable *t, const char *name)
{
    if (t->nslots == 0 && !names_rehash(t, 64)) {
        return -1;
    }
    size_t h = hash_name(name) & (t->nslots - 1);
    while (t->slots[h] >= 0) {
        if (strcmp(t->names[t->slots[h]], name) == 0) {
            return t->slots[h];
        }
        h = (h + 1) & (t->nslots - 1);
    }
    if (t->count == t->cap) {
        const size_t cap = t->cap > 0 ? t->cap * 2 : 32;
        char **nn = realloc(t->names, cap * sizeof(char *));
        if (nn == NULL) {
            return -1;
        }
        t->names = nn;
        t->cap   = cap;
    }
    const size_t len = strlen(name) + 1;
    char *copy = malloc(len);
    if (copy == NULL) {
        return -1;
    }
    memcpy(copy, name, len);
    const int id = (int)t->count;
    t->names[t->count++] = copy;
    t->slots[h] = id;
    if (t->count * 2 > t->nslots && !names_rehash(t, t->nslots * 2)) {
        return -1;
    }
    return id;
}

void matrix_names_free(char **names, size_t count)
{
    if (names == NULL) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        free(names[i]);
    }
    free(names);
}

MatrixStatus matrix_read_named_edge_list(FILE *in, bool directed, MatrixGraph **out,
                                         char ***names_out, size_t *count)
{
    if (in == NULL || out == NULL || names_out == NULL || count == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    char *line = NULL;
    size_t cap = 0;
    EdgeBuf edges;
    NameTable names;
    memset(&edges, 0, sizeof(edges));
    memset(&names, 0, sizeof(names));
    MatrixStatus st = GRAPH_OK;

    while (st == GRAPH_OK && read_line(in, &line, &cap)) {
        char *p = skip_space(line);
        if (*p == '\0') {
            continue;
        }
        if (*p == '#') {
            p = skip_space(p + 1);
            if (strncmp(p, "directed", 8) == 0) {
                directed = true;
            } else if (strncmp(p, "undirected", 10) == 0) {
                directed = false;
            }
            continue;
        }
        char *tok[3];
        size_t nt = 0;
        while (*p != '\0' && nt < 4) {
            if (nt == 3) {
                nt = 4;         /* more than three fields */
                break;
            }
            tok[nt++] = p;
            while (*p != '\0' && !isspace((unsigned char)*p)) {
                ++p;
            }
            if (*p != '\0') {
                *p++ = '\0';
            }
            p = skip_space(p);
        }
        if (nt < 2 || nt > 3) {
            st = GRAPH_ERR_IO;
            break;
        }
        double w = 1.0;
        if (nt == 3) {
            char *end;
            w = strtod(tok[2], &end);
            if (end == tok[2] || *end != '\0' || !isfinite(w)) {
                st = GRAPH_ERR_IO;
                break;
            }
        }
        const int u = names_intern(&names, tok[0]);
        const int v = names_intern(&names, tok[1]);
        if (u < 0 || v < 0 || !edgebuf_push(&edges, u, v, w)) {
            st = GRAPH_ERR_OUT_OF_MEMORY;
        }
    }
    free(line);
    if (st == GRAPH_OK && ferror(in)) {
        st = GRAPH_ERR_IO;
    }
    if (st == GRAPH_OK) {
        st = build_graph(&edges, names.count, directed, out);
    }
    edgebuf_free(&edges);
    if (st != GRAPH_OK) {
        names_free(&names);
        return st;
    }
    free(names.slots);
    *names_out = names.names;      /* ownership passes to the caller */
    *count     = names.count;
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Binary snapshot                                                     */
/* ------------------------------------------------------------------ */

static const unsigned char BIN_MAGIC[4] = { 'M', 'X', 'G', '1' };

static void put_u64(unsigned char *b, uint64_t x)
{
    for (int i = 0; i < 8; ++i) {
        b[i] = (unsigned char)(x >> (8 * i));
    }
}

static uint64_t get_u64(const unsigned char *b)
{
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) {
        x |= (uint64_t)b[i] << (8 * i);
    }
    return x;
}

static bool write_u64(FILE *out, uint64_t x)
{
    unsigned char b[8];
    put_u64(b, x);
    return fwrite(b, 1, 8, out) == 8;
}

static bool read_u64(FILE *in, uint64_t *x)
{
    unsigned char b[8];
    if (fread(b, 1, 8, in) != 8) {
        return false;
    }
    *x = get_u64(b);
    return true;
}

MatrixStatus matrix_save_binary(const MatrixGraph *g, FILE *out)
{
    if (g == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = matrix_graph_id_bound(g);
    const bool directed = matrix_graph_is_directed(g);
    bool ok = fwrite(BIN_MAGIC, 1, 4, out) == 4 && fputc(directed ? 1 : 0, out) != EOF &&
              write_u64(out, n);
    for (size_t v = 0; ok && v < n; ++v) {
        ok = fputc(matrix_graph_is_alive(g, (int)v) ? 1 : 0, out) != EOF;
    }
    uint64_t stored = 0;
    for (size_t u = 0; u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        while (matrix_graph_next_neighbor(g, (int)u, &cur, &v, NULL)) {
            stored += (directed || v >= (int)u) ? 1u : 0u;
        }
    }
    ok = ok && write_u64(out, stored);
    for (size_t u = 0; ok && u < n; ++u) {
        if (!matrix_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (ok && matrix_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (directed || v >= (int)u) {
                uint64_t bits;
                memcpy(&bits, &w, sizeof(bits));
                ok = write_u64(out, (uint64_t)u) && write_u64(out, (uint64_t)v) &&
                     write_u64(out, bits);
            }
        }
    }
    return ok && !ferror(out) ? GRAPH_OK : GRAPH_ERR_IO;
}

MatrixStatus matrix_load_binary(FILE *in, MatrixGraph **out)
{
    if (in == NULL || out == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    unsigned char magic[4];
    int dir;
    uint64_t n, stored;
    if (fread(magic, 1, 4, in) != 4 || memcmp(magic, BIN_MAGIC, 4) != 0 ||
        (dir = fgetc(in)) == EOF || (dir != 0 && dir != 1) || !read_u64(in, &n) ||
        n >= INT_MAX) {
        return GRAPH_ERR_IO;
    }
    unsigned char *alive = malloc(n > 0 ? (size_t)n : 1);
    if (alive == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    if (n > 0 && fread(alive, 1, (size_t)n, in) != (size_t)n) {
        free(alive);
        return GRAPH_ERR_IO;
    }
    MatrixGraph *g = matrix_graph_create((size_t)n, dir == 1);
    if (g == NULL || !read_u64(in, &stored)) {
        free(alive);
        matrix_graph_destroy(g);
        return g == NULL ? (graph_last_error() == GRAPH_ERR_TOO_LARGE ? GRAPH_ERR_TOO_LARGE
                                                                      : GRAPH_ERR_OUT_OF_MEMORY)
                         : GRAPH_ERR_IO;
    }
    MatrixStatus st = GRAPH_OK;
    for (size_t v = 0; v < n && st == GRAPH_OK; ++v) {
        if (alive[v] > 1 || (alive[v] == 0 && !matrix_graph_remove_vertex(g, (int)v))) {
            st = GRAPH_ERR_IO;
        }
    }
    for (uint64_t i = 0; i < stored && st == GRAPH_OK; ++i) {
        uint64_t u, v, bits;
        double w;
        if (!read_u64(in, &u) || !read_u64(in, &v) || !read_u64(in, &bits) || u >= n || v >= n) {
            st = GRAPH_ERR_IO;
            break;
        }
        memcpy(&w, &bits, sizeof(w));
        if (!matrix_graph_add_edge(g, (int)u, (int)v, w)) {
            /* Over the memory budget or out of memory are reported as such; anything else is a
             * bad file (dead endpoint, duplicate or non-finite weight). */
            const GraphStatus why = graph_last_error();
            st = (why == GRAPH_ERR_TOO_LARGE || why == GRAPH_ERR_OUT_OF_MEMORY) ? why : GRAPH_ERR_IO;
        }
    }
    free(alive);
    if (st != GRAPH_OK) {
        matrix_graph_destroy(g);
        return st;
    }
    *out = g;
    return GRAPH_OK;
}
