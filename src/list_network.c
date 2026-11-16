/*
 * list_network.c — link analysis and community structure.
 * See list_algorithms.h.
 */

#include "list_util.h"

/* ------------------------------------------------------------------ */
/*  PageRank                                                            */
/* ------------------------------------------------------------------ */

ListStatus list_pagerank(const ListGraph *g, double damping, double tolerance,
                             size_t max_iterations, double *out, size_t *iterations)
{
    if (g == NULL || out == NULL || !(damping > 0.0 && damping < 1.0) || !(tolerance > 0.0)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = list_graph_id_bound(g);
    const size_t vc = list_graph_vertex_count(g);
    mu_fill_double(out, n, 0.0);
    if (iterations != NULL) {
        *iterations = 0;
    }
    if (vc == 0) {
        return GRAPH_OK;
    }
    double *outw = calloc(n > 0 ? n : 1, sizeof(double));
    double *next = malloc((n > 0 ? n : 1) * sizeof(double));
    if (outw == NULL || next == NULL) {
        free(outw);
        free(next);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t u = 0; u < n; ++u) {
        if (!list_graph_is_alive(g, (int)u)) {
            continue;
        }
        size_t cur = 0;
        int    v;
        double w;
        while (list_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
            if (w < 0.0) {
                free(outw);
                free(next);
                return GRAPH_ERR_INVALID_ARGUMENT;
            }
            outw[u] += w;
        }
        out[u] = 1.0 / (double)vc;
    }
    ListStatus st = GRAPH_ERR_NO_CONVERGENCE;
    size_t it = 0;
    while (it < max_iterations) {
        double dangling = 0.0;
        for (size_t u = 0; u < n; ++u) {
            if (list_graph_is_alive(g, (int)u) && outw[u] <= 0.0) {
                dangling += out[u];
            }
        }
        const double base = (1.0 - damping) / (double)vc + damping * dangling / (double)vc;
        for (size_t u = 0; u < n; ++u) {
            next[u] = list_graph_is_alive(g, (int)u) ? base : 0.0;
        }
        for (size_t u = 0; u < n; ++u) {
            if (!list_graph_is_alive(g, (int)u) || outw[u] <= 0.0) {
                continue;
            }
            size_t cur = 0;
            int    v;
            double w;
            while (list_graph_next_neighbor(g, (int)u, &cur, &v, &w)) {
                next[v] += damping * out[u] * w / outw[u];
            }
        }
        double diff = 0.0;
        for (size_t u = 0; u < n; ++u) {
            diff += fabs(next[u] - out[u]);
            out[u] = next[u];
        }
        ++it;
        if (diff < tolerance) {
            st = GRAPH_OK;
            break;
        }
    }
    if (iterations != NULL) {
        *iterations = it;
    }
    free(outw);
    free(next);
    return st;
}

/* ------------------------------------------------------------------ */
/*  HITS and eigenvector centrality                                     */
/* ------------------------------------------------------------------ */

static void l2_normalize(double *x, size_t n)
{
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        s += x[i] * x[i];
    }
    s = sqrt(s);
    if (s > 0.0) {
        for (size_t i = 0; i < n; ++i) {
            x[i] /= s;
        }
    }
}

ListStatus list_hits(const ListGraph *g, double tolerance, size_t max_iterations,
                         double *hub, double *authority)
{
    if (g == NULL || hub == NULL || authority == NULL || !(tolerance > 0.0)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    double *h2 = malloc((n > 0 ? n : 1) * sizeof(double));
    double *a2 = malloc((n > 0 ? n : 1) * sizeof(double));
    if (h2 == NULL || a2 == NULL) {
        free(h2);
        free(a2);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        hub[i] = authority[i] = list_graph_is_alive(g, (int)i) ? 1.0 : 0.0;
    }
    l2_normalize(hub, n);
    l2_normalize(authority, n);
    ListStatus st = GRAPH_ERR_NO_CONVERGENCE;
    for (size_t it = 0; it < max_iterations; ++it) {
        mu_fill_double(a2, n, 0.0);
        mu_fill_double(h2, n, 0.0);
        for (size_t u = 0; u < n; ++u) {      /* authority: sum of hubs pointing in */
            size_t cur = 0;
            int    v;
            while (list_graph_is_alive(g, (int)u) &&
                   list_graph_next_neighbor(g, (int)u, &cur, &v, NULL)) {
                a2[v] += hub[u];
            }
        }
        l2_normalize(a2, n);
        for (size_t u = 0; u < n; ++u) {      /* hub: sum of authorities pointed at */
            size_t cur = 0;
            int    v;
            while (list_graph_is_alive(g, (int)u) &&
                   list_graph_next_neighbor(g, (int)u, &cur, &v, NULL)) {
                h2[u] += a2[v];
            }
        }
        l2_normalize(h2, n);
        double diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            diff += fabs(h2[i] - hub[i]) + fabs(a2[i] - authority[i]);
            hub[i]       = h2[i];
            authority[i] = a2[i];
        }
        if (diff < tolerance) {
            st = GRAPH_OK;
            break;
        }
    }
    free(h2);
    free(a2);
    return st;
}

ListStatus list_eigenvector_centrality(const ListGraph *g, double tolerance,
                                           size_t max_iterations, double *out)
{
    if (g == NULL || out == NULL || !(tolerance > 0.0)) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    double *next = malloc((n > 0 ? n : 1) * sizeof(double));
    if (next == NULL) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        out[i] = list_graph_is_alive(g, (int)i) ? 1.0 : 0.0;
    }
    l2_normalize(out, n);
    ListStatus st = GRAPH_ERR_NO_CONVERGENCE;
    for (size_t it = 0; it < max_iterations; ++it) {
        /* (A^T + I) x : the identity shift stops bipartite graphs from oscillating. */
        memcpy(next, out, n * sizeof(double));
        for (size_t u = 0; u < n; ++u) {
            size_t cur = 0;
            int    v;
            while (list_graph_is_alive(g, (int)u) &&
                   list_graph_next_neighbor(g, (int)u, &cur, &v, NULL)) {
                next[v] += out[u];
            }
        }
        l2_normalize(next, n);
        double diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            diff += fabs(next[i] - out[i]);
            out[i] = next[i];
        }
        if (diff < tolerance) {
            st = GRAPH_OK;
            break;
        }
    }
    free(next);
    return st;
}

/* ------------------------------------------------------------------ */
/*  k-core decomposition (Batagelj-Zaversnik)                           */
/* ------------------------------------------------------------------ */

ListStatus list_core_numbers(const ListGraph *g, int *core)
{
    if (g == NULL || core == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n  = list_graph_id_bound(g);
    const size_t vc = list_graph_vertex_count(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    size_t md = 0;
    int *deg  = malloc((n > 0 ? n : 1) * sizeof(int));
    int *pos  = malloc((n > 0 ? n : 1) * sizeof(int));
    int *vert = malloc((vc > 0 ? vc : 1) * sizeof(int));
    if (deg == NULL || pos == NULL || vert == NULL) {
        free(deg); free(pos); free(vert);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t v = 0; v < n; ++v) {
        deg[v] = (int)mu_sym_degree(&sym, (int)v);
        core[v] = -1;
        if (list_graph_is_alive(g, (int)v) && (size_t)deg[v] > md) {
            md = (size_t)deg[v];
        }
    }
    int *bin = calloc(md + 2, sizeof(int));
    if (bin == NULL) {
        free(deg); free(pos); free(vert);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t v = 0; v < n; ++v) {
        if (list_graph_is_alive(g, (int)v)) {
            ++bin[deg[v]];
        }
    }
    int start = 0;
    for (size_t d = 0; d <= md; ++d) {
        const int num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (size_t v = 0; v < n; ++v) {
        if (list_graph_is_alive(g, (int)v)) {
            pos[v] = bin[deg[v]];
            vert[pos[v]] = (int)v;
            ++bin[deg[v]];
        }
    }
    for (size_t d = md; d >= 1; --d) {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;
    for (size_t i = 0; i < vc; ++i) {
        const int v = vert[i];
        core[v] = deg[v];
        for (size_t k = sym.off[v]; k < sym.off[v + 1]; ++k) {
            const int u = sym.to[k];
            if (deg[u] > deg[v]) {
                const int du = deg[u];
                const int pu = pos[u];
                const int pw = bin[du];
                const int w  = vert[pw];
                if (u != w) {
                    pos[u]   = pw;
                    vert[pu] = w;
                    pos[w]   = pu;
                    vert[pw] = u;
                }
                ++bin[du];
                --deg[u];
            }
        }
    }
    free(deg); free(pos); free(vert); free(bin);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Label propagation                                                   */
/* ------------------------------------------------------------------ */

/** Renumber labels 0.. in order of first appearance among alive vertices. */
static size_t renumber(const ListGraph *g, int *label, size_t n)
{
    int *map = malloc((n > 0 ? n : 1) * sizeof(int));
    if (map == NULL) {
        return 0;
    }
    mu_fill_int(map, n, -1);
    size_t next = 0;
    for (size_t v = 0; v < n; ++v) {
        if (!list_graph_is_alive(g, (int)v)) {
            label[v] = -1;
            continue;
        }
        if (map[label[v]] < 0) {
            map[label[v]] = (int)next++;
        }
        label[v] = map[label[v]];
    }
    free(map);
    return next;
}

ListStatus list_label_propagation(const ListGraph *g, int *community, size_t *count)
{
    if (g == NULL || community == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build(g, &sym)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t k = 0; k < sym.off[n]; ++k) {
        if (sym.w[k] < 0.0) {
            mu_sym_free(&sym);
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
    }
    double *lw      = calloc(n > 0 ? n : 1, sizeof(double));
    bool   *flag    = calloc(n > 0 ? n : 1, sizeof(bool));
    int    *touched = malloc((n > 0 ? n : 1) * sizeof(int));
    if (lw == NULL || flag == NULL || touched == NULL) {
        free(lw);
        free(flag);
        free(touched);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t v = 0; v < n; ++v) {
        community[v] = (int)v;
    }
    for (int round = 0; round < 100; ++round) {
        bool changed = false;
        for (size_t v = 0; v < n; ++v) {
            if (!list_graph_is_alive(g, (int)v) || mu_sym_degree(&sym, (int)v) == 0) {
                continue;
            }
            size_t nt = 0;
            for (size_t k = sym.off[v]; k < sym.off[v + 1]; ++k) {
                const int c = community[sym.to[k]];
                if (!flag[c]) {
                    flag[c] = true;
                    touched[nt++] = c;
                }
                lw[c] += sym.w[k];
            }
            /* Heaviest neighbouring label, ties to the smaller label. */
            int    best = -1;
            double bw   = -1.0;
            for (size_t i = 0; i < nt; ++i) {
                const int c = touched[i];
                if (lw[c] > bw || (lw[c] == bw && c < best)) {
                    best = c;
                    bw   = lw[c];
                }
            }
            const int cur = community[v];
            /* Move only on a strict gain: the total weight inside labels then
             * strictly increases, so the loop must terminate. */
            if (best != cur && bw > (flag[cur] ? lw[cur] : 0.0)) {
                community[v] = best;
                changed = true;
            }
            for (size_t i = 0; i < nt; ++i) {
                lw[touched[i]]   = 0.0;
                flag[touched[i]] = false;
            }
        }
        if (!changed) {
            break;
        }
    }
    const size_t k = renumber(g, community, n);
    if (count != NULL) {
        *count = k;
    }
    free(lw);
    free(flag);
    free(touched);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Modularity and Louvain                                              */
/* ------------------------------------------------------------------ */

ListStatus list_modularity(const ListGraph *g, const int *community, double *q)
{
    if (g == NULL || community == NULL || q == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    int maxc = -1;
    for (size_t v = 0; v < n; ++v) {
        if (!list_graph_is_alive(g, (int)v)) {
            continue;
        }
        if (community[v] < 0) {
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
        if (community[v] > maxc) {
            maxc = community[v];
        }
    }
    MuSym sym;
    if (!mu_sym_build(g, &sym)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    const size_t nc = (size_t)(maxc + 1);
    double *in  = calloc(nc + 1, sizeof(double));
    double *tot = calloc(nc + 1, sizeof(double));
    if (in == NULL || tot == NULL) {
        free(in);
        free(tot);
        mu_sym_free(&sym);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    double m2 = 0.0;
    for (size_t u = 0; u < n; ++u) {
        if (!list_graph_is_alive(g, (int)u)) {
            continue;
        }
        for (size_t k = sym.off[u]; k < sym.off[u + 1]; ++k) {
            const double w = sym.w[k];
            m2 += w;
            tot[community[u]] += w;
            if (community[sym.to[k]] == community[u]) {
                in[community[u]] += w;
            }
        }
    }
    double result = 0.0;
    if (m2 > 0.0) {
        for (size_t c = 0; c < nc; ++c) {
            result += in[c] / m2 - (tot[c] / m2) * (tot[c] / m2);
        }
    }
    *q = result;
    free(in);
    free(tot);
    mu_sym_free(&sym);
    return GRAPH_OK;
}

/** A weighted undirected graph level: CSR without self-loops plus loop weights. */
typedef struct {
    size_t  n;
    size_t *off;
    int    *to;
    double *w;
    double *loop;     /* internal weight of each node, each edge counted once */
    double *k;        /* weighted degree, loops counted twice                 */
} LvGraph;

static void lv_free(LvGraph *l)
{
    free(l->off);
    free(l->to);
    free(l->w);
    free(l->loop);
    free(l->k);
    memset(l, 0, sizeof(*l));
}

ListStatus list_louvain(const ListGraph *g, int *community, double *modularity, size_t *count)
{
    if (g == NULL || community == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build(g, &sym)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t k = 0; k < sym.off[n]; ++k) {
        if (sym.w[k] < 0.0) {
            mu_sym_free(&sym);
            return GRAPH_ERR_INVALID_ARGUMENT;
        }
    }

    /* Level 0 works on the original IDs directly (dead IDs are isolated nodes). */
    LvGraph cur;
    memset(&cur, 0, sizeof(cur));
    cur.n    = n;
    cur.off  = sym.off;  sym.off = NULL;
    cur.to   = sym.to;   sym.to  = NULL;
    cur.w    = sym.w;    sym.w   = NULL;
    cur.loop = calloc(n > 0 ? n : 1, sizeof(double));
    cur.k    = calloc(n > 0 ? n : 1, sizeof(double));
    int *node_of = malloc((n > 0 ? n : 1) * sizeof(int));     /* original ID -> current node */
    if (cur.loop == NULL || cur.k == NULL || node_of == NULL) {
        lv_free(&cur);
        free(node_of);
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < n; ++i) {
        node_of[i] = (int)i;
    }
    ListStatus st = GRAPH_OK;

    for (int level = 0; level < 64 && st == GRAPH_OK; ++level) {
        const size_t cn = cur.n;
        double m2 = 0.0;
        for (size_t i = 0; i < cn; ++i) {
            double s = 2.0 * cur.loop[i];
            for (size_t e = cur.off[i]; e < cur.off[i + 1]; ++e) {
                s += cur.w[e];
            }
            cur.k[i] = s;
            m2 += s;
        }
        if (m2 <= 0.0) {
            break;
        }
        int    *comm = malloc((cn > 0 ? cn : 1) * sizeof(int));
        double *tot  = malloc((cn > 0 ? cn : 1) * sizeof(double));
        double *nw   = calloc(cn > 0 ? cn : 1, sizeof(double));
        int    *touched = malloc((cn > 0 ? cn : 1) * sizeof(int));
        bool   *tflag   = calloc(cn > 0 ? cn : 1, sizeof(bool));
        if (comm == NULL || tot == NULL || nw == NULL || touched == NULL || tflag == NULL) {
            free(comm); free(tot); free(nw); free(touched); free(tflag);
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        for (size_t i = 0; i < cn; ++i) {
            comm[i] = (int)i;
            tot[i]  = cur.k[i];
        }
        bool improved = false;
        for (int pass = 0; pass < 100; ++pass) {
            bool moved = false;
            for (size_t i = 0; i < cn; ++i) {
                const int ci = comm[i];
                size_t nt = 0;
                for (size_t e = cur.off[i]; e < cur.off[i + 1]; ++e) {
                    const int c = comm[cur.to[e]];
                    if (!tflag[c]) {
                        tflag[c] = true;
                        touched[nt++] = c;
                    }
                    nw[c] += cur.w[e];
                }
                tot[ci] -= cur.k[i];
                int best = ci;
                double bgain = nw[ci] - tot[ci] * cur.k[i] / m2;
                for (size_t t = 0; t < nt; ++t) {
                    const int c = touched[t];
                    const double gain = nw[c] - tot[c] * cur.k[i] / m2;
                    if (gain > bgain + 1e-12) {
                        bgain = gain;
                        best  = c;
                    }
                }
                tot[best] += cur.k[i];
                comm[i] = best;
                if (best != ci) {
                    moved = improved = true;
                }
                for (size_t t = 0; t < nt; ++t) {
                    nw[touched[t]]    = 0.0;
                    tflag[touched[t]] = false;
                }
            }
            if (!moved) {
                break;
            }
        }
        free(tot);
        free(nw);
        free(touched);
        free(tflag);
        if (!improved) {
            free(comm);
            break;
        }

        /* Renumber communities and fold them into the next level's nodes. */
        int *label = malloc(cn * sizeof(int));
        if (label == NULL) {
            free(comm);
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        mu_fill_int(label, cn, -1);
        size_t nn = 0;
        for (size_t i = 0; i < cn; ++i) {
            if (label[comm[i]] < 0) {
                label[comm[i]] = (int)nn++;
            }
        }
        for (size_t i = 0; i < n; ++i) {
            node_of[i] = label[comm[node_of[i]]];
        }
        /* Members grouped by new node. */
        size_t *moff = calloc(nn + 2, sizeof(size_t));
        int    *mem  = malloc((cn > 0 ? cn : 1) * sizeof(int));
        LvGraph nxt;
        memset(&nxt, 0, sizeof(nxt));
        nxt.n    = nn;
        nxt.off  = calloc(nn + 1, sizeof(size_t));
        nxt.to   = malloc((cur.off[cn] > 0 ? cur.off[cn] : 1) * sizeof(int));
        nxt.w    = malloc((cur.off[cn] > 0 ? cur.off[cn] : 1) * sizeof(double));
        nxt.loop = calloc(nn > 0 ? nn : 1, sizeof(double));
        nxt.k    = calloc(nn > 0 ? nn : 1, sizeof(double));
        double *acc = calloc(nn > 0 ? nn : 1, sizeof(double));
        int    *tch = malloc((nn > 0 ? nn : 1) * sizeof(int));
        bool   *aflag = calloc(nn > 0 ? nn : 1, sizeof(bool));
        if (moff == NULL || mem == NULL || nxt.off == NULL || nxt.to == NULL || nxt.w == NULL ||
            nxt.loop == NULL || nxt.k == NULL || acc == NULL || tch == NULL || aflag == NULL) {
            free(moff); free(mem); free(acc); free(tch); free(aflag); free(label); free(comm);
            lv_free(&nxt);
            st = GRAPH_ERR_OUT_OF_MEMORY;
            break;
        }
        for (size_t i = 0; i < cn; ++i) {
            ++moff[label[comm[i]] + 1];
        }
        for (size_t c = 0; c < nn; ++c) {
            moff[c + 1] += moff[c];
        }
        {
            size_t *fill = malloc((nn + 1) * sizeof(size_t));
            if (fill == NULL) {
                free(moff); free(mem); free(acc); free(tch); free(aflag); free(label); free(comm);
                lv_free(&nxt);
                st = GRAPH_ERR_OUT_OF_MEMORY;
                break;
            }
            memcpy(fill, moff, nn * sizeof(size_t));
            for (size_t i = 0; i < cn; ++i) {
                mem[fill[label[comm[i]]]++] = (int)i;
            }
            free(fill);
        }
        size_t used = 0;
        for (size_t c = 0; c < nn; ++c) {
            size_t nt = 0;
            double internal = 0.0;
            for (size_t mi = moff[c]; mi < moff[c + 1]; ++mi) {
                const int i = mem[mi];
                nxt.loop[c] += cur.loop[i];
                for (size_t e = cur.off[i]; e < cur.off[i + 1]; ++e) {
                    const int d = label[comm[cur.to[e]]];
                    if ((size_t)d == c) {
                        internal += cur.w[e];     /* seen from both ends: halve below */
                    } else {
                        if (!aflag[d]) {
                            aflag[d] = true;
                            tch[nt++] = d;
                        }
                        acc[d] += cur.w[e];
                    }
                }
            }
            nxt.loop[c] += internal / 2.0;
            nxt.off[c] = used;
            for (size_t t = 0; t < nt; ++t) {
                nxt.to[used] = tch[t];
                nxt.w[used]  = acc[tch[t]];
                acc[tch[t]]  = 0.0;
                aflag[tch[t]] = false;
                ++used;
            }
        }
        nxt.off[nn] = used;
        free(moff); free(mem); free(acc); free(tch); free(aflag); free(label); free(comm);
        lv_free(&cur);
        cur = nxt;
    }
    lv_free(&cur);
    mu_sym_free(&sym);
    if (st != GRAPH_OK) {
        free(node_of);
        return st;
    }
    for (size_t i = 0; i < n; ++i) {
        community[i] = list_graph_is_alive(g, (int)i) ? node_of[i] : -1;
    }
    free(node_of);
    const size_t k = renumber(g, community, n);
    if (count != NULL) {
        *count = k;
    }
    if (modularity != NULL) {
        return list_modularity(g, community, modularity);
    }
    return GRAPH_OK;
}

/* ------------------------------------------------------------------ */
/*  Assortativity                                                       */
/* ------------------------------------------------------------------ */

ListStatus list_degree_assortativity(const ListGraph *g, double *r)
{
    if (g == NULL || r == NULL) {
        return GRAPH_ERR_INVALID_ARGUMENT;
    }
    const size_t n = list_graph_id_bound(g);
    MuSym sym;
    if (!mu_sym_build_ex(g, &sym, false)) {
        return GRAPH_ERR_OUT_OF_MEMORY;
    }
    double m = 0.0, sjk = 0.0, sj = 0.0, sj2 = 0.0;
    for (size_t u = 0; u < n; ++u) {
        for (size_t e = sym.off[u]; e < sym.off[u + 1]; ++e) {
            const size_t v = (size_t)sym.to[e];
            if (v <= u) {
                continue;       /* each undirected pair once */
            }
            const double j = (double)mu_sym_degree(&sym, (int)u);
            const double k = (double)mu_sym_degree(&sym, (int)v);
            m   += 1.0;
            sjk += j * k;
            sj  += j + k;
            sj2 += j * j + k * k;
        }
    }
    mu_sym_free(&sym);
    if (m == 0.0) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    const double mean = sj / (2.0 * m);
    const double den  = sj2 / (2.0 * m) - mean * mean;
    if (den <= 1e-15) {
        return GRAPH_ERR_UNSUPPORTED;
    }
    *r = (sjk / m - mean * mean) / den;
    return GRAPH_OK;
}
