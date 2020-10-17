/*
 * minorfinder.c: Game about finding a specific minor in a bigger
 * graph. To extract the minor from the bigger graph the player
 * can either contract edges, delete edges or delete vertices with
 * degree zero, i.e. vertices that aren't adjacent to any other
 * vertices of the graph.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include <time.h>

#include "puzzles.h"
#include "tree234.h"

/* debug mode */
#define DEBUG true

#define BENCHMARKS

/* enable or disable console log */
#if DEBUG
#define LOG(x) (printf x)
#else
#define LOG(x)
#endif

/* type aliases */
#define uint8 unsigned char
#define uint unsigned int
#define ulong unsigned long

enum {
    COL_SYSBACKGROUND,
    COL_BACKGROUND,
    COL_OUTLINE,
    COL_GRIDBORDER,
    COL_BASEPOINT,
    COL_MINORPOINT,
    COL_DRAGPOINT,
    COL_DELPOINT,
    COL_HIDEPOINT,
    COL_POINTOUTLINE,
    COL_EDGE,
    COL_DELEDGE,
    COL_HIDEEDGE,
    COL_FLASH,
    COL_FLASH2,
    COL_TEXT,
    COL_TEXTBACKGROUND,
#if DEBUG
    COL_SUBPOINT,
#endif
    NCOLOURS
};

/*
 * A grid, defines a region of a window
 */
typedef struct grid {

    /* grid offset */
    int x_off, y_off;

    /* grid size relativ to coordinate limit */
    float relsize;

} grid;

/*
 * A point in a grid, two rational coordinates and a denominator determine
 * the position of a point in a grid.
 */
typedef struct point {
    
    /* rational coordinates */
    long x, y;

    /* denominator - should always be 1 here */
    long d;
    
} point;

/*
 * A vertex that corresponds to a point that belongs to a graph
 */
typedef struct vertex {

     /* index in points array */
    int idx;

    /* number of edges that are incident to the vertex */
    int deg;

} vertex;

/*
 * An edge that connects two vertices of a graph. Despite the fact that we
 * use src and tgt as identifiers for the vertices that are incident to an
 * edge, edges do not have a direction.
 */
typedef struct edge {

    /* vertices that are incident to the edge */
    int src, tgt;

} edge;

/*
 * An undirected graph that consists of a set of points (vertices) and edges
 * that connect these vertices.
 */
typedef struct graph {

    /* number of references to the graph - for deallocation */
    int refcount;

    /*
     * the grid in which the graph is drawn - determines its coordinate
     * offset
     */
    grid grid;

    /* array of points - static size, current point coordinates */
    point* points;
    /* array of vertices - static size, current vertex degrees */
    vertex* vtcs;
    /* 234-tree of vertices - current visible vertices */
    tree234* vertices;

    /* 234-tree of edges - current edges */
    tree234* edges;

} graph;

enum game_mode {
    NORMAL,
    WAGNER
};

struct game_params {

    /* either NORMAL or WAGNER mode */
    enum game_mode mode;

    /* number of base graph points */
    int n_base;

    /* number of minor graph points */
    int n_min;

};

const struct game_params normal_presets[] = {
    { NORMAL, 15, 4 },
    { NORMAL, 19, 5 },
    { NORMAL, 23, 6 }
};

const struct game_params wagner_presets[] = {
    { WAGNER, 19, 5 },
    { WAGNER, 23, 6 }
};

#define DEFAULT_PRESET normal_presets[0]

struct game_state {

    game_params params;

    graph* base;
    graph* minor;

    /* player solved game */
    bool solved;
    /* player used solve function */
    bool cheated;

};

static game_params *default_params(void)
{
    game_params *ret = snew(game_params);

    *ret = DEFAULT_PRESET;

    return ret;
}

struct preset_menu* preset_menu(void)
{
    char buf[80];
    int i;
    struct preset_menu* ret = preset_menu_new();

    sprintf(buf, "Normal");
    struct preset_menu* normal = preset_menu_add_submenu(ret, dupstr(buf));
    for (i = 0; i < lenof(normal_presets); i++)
    {
        game_params* params = default_params();
        *params = normal_presets[i];
        sprintf(buf, "%d base, %d minor points", params->n_base, params->n_min);
        preset_menu_add_preset(normal, dupstr(buf), params);
    }

    sprintf(buf, "Wagner");
    struct preset_menu* wagner = preset_menu_add_submenu(ret, dupstr(buf));
    for (i = 0; i < lenof(wagner_presets); i++)
    {
        game_params* params = default_params();
        *params = wagner_presets[i];
        switch (params->n_min)
        {
            case 5:
                sprintf(buf, "K_5");
                break;
            case 6:
                sprintf(buf, "K_3,3");
                break;
            default:;
        }
        preset_menu_add_preset(wagner, dupstr(buf), params);
    }

    return ret;
}

static void free_params(game_params *params)
{
    sfree(params);
}

static game_params *dup_params(const game_params *params)
{
    game_params *ret = snew(game_params);
    *ret = *params;		       /* structure copy */
    return ret;
}

static void decode_params(game_params *params, char const *string)
{
    int mode;
    if (sscanf(string, "%d-%dx%d", &mode, &params->n_base, &params->n_min) != 3)
    {
        /* params encoding was incorrect */
        *params = DEFAULT_PRESET;
    }
    else
    {
        params->mode = mode;
    }
}

static char *encode_params(const game_params *params, bool full)
{
    char buf[80];

    sprintf(buf, "%d-%dx%d", params->mode, params->n_base, params->n_min);
    
    return dupstr(buf);
}

static config_item *game_configure(const game_params *params)
{
    return NULL;
}

static game_params *custom_params(const config_item *cfg)
{
    return NULL;
}

static const char *validate_params(const game_params *params, bool full)
{
    switch (params->mode)
    {
        case NORMAL:
            if (params->n_base < normal_presets[0].n_base
                || params->n_base > normal_presets[lenof(normal_presets)-1].n_base)
                return "Number of base graph points is invalid";
            else if (params->n_min < normal_presets[0].n_min
                || params->n_min > normal_presets[lenof(normal_presets)-1].n_min)
                return "Number of minor points is invalid";
            break;
        case WAGNER:
            if (params->n_base < wagner_presets[0].n_base
                || params->n_base > wagner_presets[lenof(wagner_presets)-1].n_base)
                return "Number of base graph points is invalid";
            else if (params->n_min < wagner_presets[0].n_min
                || params->n_min > wagner_presets[lenof(wagner_presets)-1].n_min)
                return "Number of minor points is invalid";
            break;
        default:;
    }
    
    return NULL;
}

/*
 * 
 * The functions and structures below are copied from the untangle
 * backend.
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

static ulong squarert(ulong n) {
    ulong d, a, b, di;

    d = n;
    a = 0;
    b = 1L << 30;		       /* largest available power of 4 */
    do {
        a >>= 1;
        di = 2*a + b;
        if (di <= d) {
            d -= di;
            a += b;
        }
        b >>= 2;
    } while (b);

    return a;
}

static void addedge(tree234 *edges, int s, int t)
{
    edge *e = snew(edge);

#if DEBUG
    assert(s != t);
#endif

    e->src = min(s, t);
    e->tgt = max(s, t);

    add234(edges, e);
}

static bool isedge(tree234 *edges, int s, int t)
{
    edge e;

#if DEBUG
    assert(s != t);
#endif

    e.src = min(s, t);
    e.tgt = max(s, t);

    return find234(edges, &e, NULL) != NULL;
}

static int edgecmpC(const void *av, const void *bv)
{
    const edge *a = (const edge *)av;
    const edge *b = (const edge *)bv;

    if (a->src < b->src)
	return -1;
    else if (a->src > b->src)
	return +1;
    else if (a->tgt < b->tgt)
	return -1;
    else if (a->tgt > b->tgt)
	return +1;
    return 0;
}

static int edgecmp(void *av, void *bv)
{
    return edgecmpC(av, bv);
}

static int vertcmpC(const void *av, const void *bv)
{
    const vertex *a = (vertex *)av;
    const vertex *b = (vertex *)bv;
    
    if (a->deg < b->deg)
	return -1;
    else if (a->deg > b->deg)
	return +1;
    else if (a->idx < b->idx)
	return -1;
    else if (a->idx > b->idx)
	return +1;
    return 0;
}

static int vertcmp(void *av, void *bv)
{
    return vertcmpC(av, bv);
}

/*
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * The functions and structures above are copied from the untangle
 * backend.
 * 
 */

#define POINTRADIUS 5
#define CROSSPOINT_THRESHOLD (POINTRADIUS / 2.0)

#define square(x) ((x) * (x))

/*
 * Check whether the edge between s and t crosses the point p
 */
static bool crosspoint(point s, point t, point p)
{
    double dist_st = sqrt(square(t.x - s.x) + square(t.y - s.y));
    double dist_sp = sqrt(square(p.x - s.x) + square(p.y - s.y));
    double dist_pt = sqrt(square(t.x - p.x) + square(t.y - p.y));

    return dist_sp + dist_pt - dist_st < CROSSPOINT_THRESHOLD;
}

/*
 * Add edges between the vertices in the range [offa, offa + cnta] and vertices in the range
 * [offb, offb + cntb]. Make sure that edges don't cross other edges or points and the degree
 * of the involved vertices doesn't increase beyond max_deg.
 */
static void addedges(tree234* edges, vertex* vertices, point* points, int offa, int cnta,
                    int maxdega, int offb, int cntb, int maxdegb, const int n, random_state* rs)
{
    bool* contains;
    int i;
#if DEBUG
    int j;
#endif
    vertex* vxa;
    vertex* vxb;
    tree234* vtcs;
    tree234* vtcsa;
    tree234** vtcsb;

#if DEBUG
    assert(offa >= 0 && offa < n);
    assert(offb >= 0 && offb < n);
    assert(cnta >= 0 || cntb >= 0);
    assert(cnta <= n - offa && cntb <= n - offb);
    assert(maxdega > 0 || maxdegb > 0);
    assert(maxdega < n && maxdegb < n);
#endif
    if (cnta < 0) cnta = offb + cntb - offa;
    else if (cntb < 0) cntb = offa + cnta - offb;
    if (cnta < cntb) {
        int tmp = cnta;
        cnta = cntb;
        cntb = tmp;
        tmp = offa;
        offa = offb;
        offb = tmp;
    }
    if (maxdega < 0) maxdega = maxdegb;
    else if (maxdegb < 0) maxdegb = maxdega;

    /* add all vertices in range a to a 234-tree */
    vtcsa = newtree234(vertcmp);
    LOG(("Initially added vertices "));
    for (i = offa; i < offa + cnta; i++)
    {
        add234(vtcsa, vertices + i);
        LOG(("%d, ", vertices[i].idx));
    }
    LOG(("to range a vertices\n"));

    /* for every vertex in range a add all vertices in range b to a 234-tree */
    vtcsb = snewn(cnta, tree234*);
    *vtcsb = newtree234(vertcmp);
    LOG(("Initially added vertices "));
    for (i = offb; i < offb + cntb; i++)
    {
        add234(*vtcsb, vertices + i);
        LOG(("%d, ", vertices[i].idx));
    }
    LOG(("to range b vertices\n"));
    for (i = 1; i < cnta; i++)
    {
        vtcsb[i] = copytree234(*vtcsb, NULL, NULL);
#if DEBUG
        assert(count234(vtcsb[i]) == count234(*vtcsb));
        for (j = 0; j < count234(vtcsb[i]); j++)
        {
            vxa = index234(vtcsb[i], j);
            vxb = index234(*vtcsb, j);
            assert(vxa->idx == vxb->idx);
        }
#endif
    }
    
    /*
     * Start adding edges. Pick the lowest degree vertex from the range a 234-tree and a
     * random vertex from the range b 234-tree and try to add an edge between them.
     * If the edge can't be added delete the range b vertex from the corresponding 234-tree.
     * Otherwise add the edge and update the involved vertices in the range a 234-tree and
     * all range b 234-trees. Repeat until there are no more edges to add.
     */
    while (count234(vtcsa))
    {
        vxa = index234(vtcsa, 0);
        vtcs = vtcsb[vxa->idx - offa];
        if (!count234(vtcs))
        {
            del234(vtcsa, vxa);
            LOG(("Removed vertex %d from range a vertices, no more edges to add\n",
                vxa->idx));
            continue;
        }

        next_vxb:
        vxb = index234(vtcs, random_upto(rs, count234(vtcs)));
        LOG(("Trying to add edge between vertices %d and %d\n", vxa->idx, vxb->idx));
        if (vxa->idx == vxb->idx || isedge(edges, vxa->idx, vxb->idx))
        {
            del234(vtcs, vxb);
            LOG(("Removed vertex %d from range b vertices of vertex %d,"\
                " the edge can't be added\n", vxb->idx, vxa->idx));
            if (count234(vtcs))
            {
                goto next_vxb;
            }
            else
            {
                del234(vtcsa, vxa);
                LOG(("Removed vertex %d from range a vertices, no more edges to add\n",
                    vxa->idx));
                continue;
            }
        }

        /* check for crossing points */
        for (i = 0; i < n; i++)
        {
            if (i == vxa->idx || i == vxb->idx)
            {
                continue;
            }
            else if (crosspoint(points[vxa->idx], points[vxb->idx], points[i]))
            {
                del234(vtcs, vxb);
                LOG(("Removed vertex %d from range b vertices of vertex %d,"\
                    " the edge crosses the point %d\n", vxb->idx, vxa->idx, i));
                if (count234(vtcs))
                {
                    goto next_vxb;
                }
                else
                {
                    del234(vtcsa, vxa);
                    LOG(("Removed vertex %d from range a vertices, no more edges to add\n",
                        vxa->idx));
                    goto next_vxa;
                }
            }
        }

        addedge(edges, vxa->idx, vxb->idx);
        LOG(("Added edge between vertices %d and %d\n", vxa->idx, vxb->idx));
        contains = snewn(cnta, bool);
        del234(vtcsa, vxa);
        for (i = 0; i < cnta; i++)
        {
            contains[i] = del234(vtcsb[i], vxa);
        }
        vxa->deg++;
        if (vxa->deg < maxdega)
        {
            add234(vtcsa, vxa);
            for (i = 0; i < cnta; i++)
            {
                if (contains[i]) add234(vtcsb[i], vxa);
            }
        }
        contains = sresize(contains, cnta + 1, bool);
        *contains = del234(vtcsa, vxb);
        for (i = 0; i < cnta; i++)
        {
            contains[i+1] = del234(vtcsb[i], vxb);
        }
        vxb->deg++;
        if (vxb->deg < maxdegb)
        {
            if (*contains) add234(vtcsa, vxb);
            for (i = 0; i < cnta; i++)
            {
                if (contains[i+1] && i != vxa->idx - offa)
                    add234(vtcsb[i], vxb);
            }
        }
        sfree(contains);
        LOG(("Updated degrees of vertices %d to %d and %d to %d\n", vxa->idx, vxa->deg,
            vxb->idx, vxb->deg));

        next_vxa:;
    }

    freetree234(vtcsa);
    for (i = 0; i < cnta; i++) freetree234(vtcsb[i]);
    sfree(vtcsb);
}

/*
 * Add edges to a graph with 5 vertices such that it becomes the K_5
 */
static void make_K_5_edges(tree234* edges, vertex* vertices, int n)
{
    int i, j;

#if DEBUG
    assert(n == 5);
#endif

    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            addedge(edges, i, j);
            vertices[i].deg++;
            vertices[j].deg++;
        }
    }
}

/*
 * Add edges to a graph with 6 vertices such that it becomes the K_33
 */
static void make_K_33_edges(tree234* edges, vertex* vertices, int n)
{
    int i;

#if DEBUG
    assert(n == 6);
#endif

    addedge(edges, 5, 2);
    addedge(edges, 5, 3);
    addedge(edges, 5, 4);

    addedge(edges, 0, 2);
    addedge(edges, 0, 3);
    addedge(edges, 0, 4);

    addedge(edges, 1, 2);
    addedge(edges, 1, 3);
    addedge(edges, 1, 4);

    for (i = 0; i < n; i++)
    {
        vertices[i].deg = 3;
    }
}

/*
 * These parameters are highly sensitive, changing them may cause problems when
 * generating new game descriptions.
 */
#define COORDMARGIN 1
#define COORDLIMIT(n) ((n) + (2 * COORDMARGIN))
#define COORDUNIT 16

#define MINORRADIUS(n) ((n) / 3)
#define SUBGRAPH_DISTANCE (2 * POINTRADIUS) + 1
#define SUBGRAPH_POINTENTROPY 2
#define OVERLAYPOINT_TRESHOLD_GAMEGEN square(4 * POINTRADIUS)

/*
 * Arrange the points of a K_33 in two rows of three points each, such that it
 * becomes the default presentation of the K_33.
 */
static void make_K_33_points(point* points, long clim, int n)
{
    int i;

#if DEBUG
    assert(n == 6);
#endif

    for (i = 0; i < n; i++)
    {
        switch (i)
        {
            case 0:
            case 3:
                points[i].x = clim * COORDUNIT / 2;
                break;
            case 1:
            case 2:
                points[i].x = clim * COORDUNIT / 5;
                break;
            case 4:
            case 5:
                points[i].x = clim * COORDUNIT * 4 / 5;
                break;
            default:;
        }
        switch (i)
        {
            case 5:
            case 0:
            case 1:
                points[i].y = clim * COORDUNIT / 4;
                break;
            case 2:
            case 3:
            case 4:
                points[i].y = clim * COORDUNIT * 3 / 4;
                break;
            default:;
        }
    }
}

/*
 * Extended edge. Stores a value proportional and the degree sum of its incident
 * vertices.
 */
typedef struct edge_ext {
    edge e;
    int degsum;
} edge_ext;

/*
 *  Compares two constant extended edges by their degsum and edge
 */
static int eextcmpC(const void *av, const void  *bv)
{
    const edge_ext *a = (edge_ext *)av;
    const edge_ext *b = (edge_ext *)bv;

    if (a->degsum < b->degsum) return -1;
    else if (a->degsum > b->degsum) return 1;
    else return edgecmpC(&a->e, &b->e);
}

/*
 * Compares two extended edges with eextcmpC
 */
static int eextcmp(void *av, void *bv)
{
    return eextcmpC(av, bv);
}

#ifdef BENCHMARKS
/*
 * Costs for visible changes to the graph to transform one into the other. Although edges are
 * moved together with their incident vertices we count that as two visible changes. One to the
 * vertex itself and one to the set of its incident edges.
 */
#define COST_SWITCHIDCS 1
#define COST_MOVEPOINT 2
#define COST_ADDEDGE 1
#define COST_DELEDGE 1

static tree234* LastBaseEdges = NULL;
static point* LastBasePts = NULL;
static int LastNBase = 0;
static int LastNMin = 0;

static int GraphDissim = 0;

/*
 * Calculate the dissimilarity of the base graphs between two game generations by comparing their
 * edges and vertices and summing up the predefined costs for the three operations add edge, delete
 * edge and move point that are required to transform the previous graph into the current.
 */
static void calc_basegraph_dissim(tree234* curr_base_edges, point* curr_base_pts, int curr_n_base,
                            int curr_n_min)
{
    int i, j, k;
    point* last_pt;
    point* curr_pt;
    edge* e;

    GraphDissim = 0;

    if (LastBaseEdges && LastBasePts && LastNBase && LastNMin
        && LastNBase == curr_n_base && LastNMin == curr_n_min)
    {
        int n_sub;
        for (i = 0; (e = index234(curr_base_edges, i)) != NULL; i++)
        {
            if (find234(LastBaseEdges, e, NULL) == NULL)
                GraphDissim += COST_ADDEDGE;
        }
        for (i = 0; (e = index234(LastBaseEdges, i)) != NULL; i++)
        {
            if (find234(curr_base_edges, e, NULL) == NULL)
                GraphDissim += COST_DELEDGE;
        }
        n_sub = LastNBase / LastNMin;
        for (i = 0; i < LastNMin; i++)
        {
            for (j = i * n_sub; j < (i + 1) * n_sub; j++)
            {
                bool moved = true;
                last_pt = LastBasePts + j;
                for (k = i * n_sub; k < (i + 1) * n_sub; k++)
                {
                    curr_pt = curr_base_pts + k;
                    if (curr_pt->x == last_pt->x
                        && curr_pt->y == last_pt->y)
                    {
                        moved = false;
                        break;
                    }
                }
                if (moved) GraphDissim += COST_MOVEPOINT;
                else GraphDissim += COST_SWITCHIDCS;
            }
        }
        for (i = LastNMin * n_sub; i < LastNBase; i++)
        {
            bool moved = true;
            last_pt = LastBasePts + i;
            for (j = LastNMin * n_sub; j < LastNBase; j++)
            {
                curr_pt = curr_base_pts + j;
                if (curr_pt->x == last_pt->x
                    && curr_pt->y == last_pt->y)
                {
                    moved = false;
                    break;
                }
            }
            if (moved) GraphDissim += COST_MOVEPOINT;
            else GraphDissim += COST_SWITCHIDCS;
        }
        while ((e = delpos234(LastBaseEdges, 0)) != NULL) sfree(e);
        freetree234(LastBaseEdges);
        sfree(LastBasePts);
    }

    LastBaseEdges = curr_base_edges;
    LastBasePts = curr_base_pts;
    LastNBase = curr_n_base;
    LastNMin = curr_n_min;
}
#endif

static char *new_game_desc(const game_params *params, random_state *rs,
			   char **aux, bool interactive)
{
    char* ret;
    
    const int n_min = params->n_min;
    const int n_base = params->n_base;
    const int n_sub = n_base / n_min;

    long i, j, k, l;
    long tmp, tmp2, tmp3;
    long coord_lim;
    long circle_rad;
    long* coords_x;
    long* coords_y;
    long* radii;
    
    double* angles;

    point* pt;
    point* pts_min;
    point* pts_base;

    vertex* vx;
    vertex* vtcs_min;
    vertex* vtcs_base;

    edge* e;
    tree234* edges_min_234;
    tree234* edges_base_234;

    coord_lim = COORDLIMIT(n_base);
    circle_rad = MINORRADIUS(n_base);

    tmp = coord_lim - (2 * COORDMARGIN);
    pts_min = snewn(n_min, point);
    /* Arrange the minor points in a circle with radius circle_rad */
    if (params->mode != WAGNER || n_min < 6)
    {
        for (i = 0; i < n_min; i++)
        {
            double angle = ((double) i * 2.0 * PI) / (double) n_min;
            pt = pts_min + i;
            pt->x = (((double) tmp / 2.0) + ((double) circle_rad * sin(angle)) + COORDMARGIN)
                    * COORDUNIT;
            pt->y = (((double) tmp / 2.0) + ((double) circle_rad * cos(angle)) + COORDMARGIN)
                    * COORDUNIT;
            pt->d = 1;
            LOG(("Assigned coordinates x:%ld, y:%ld and denominator %ld to minor point %ld\n",
                pt->x, pt->y, pt->d, i));
        }
    }
    /* Arrange the minor points in two rows of three points each */
    else
    {
        make_K_33_points(pts_min, coord_lim, n_min);
    }

    /* Add edges to the minor */
    vtcs_min = snewn(n_min, vertex);
    for (i = 0; i < n_min; i++)
    {
        vx = vtcs_min + i;
        vx->idx = i;
        vx->deg = 0;
    }
    edges_min_234 = newtree234(edgecmp);
    switch (params->mode)
    {
        case NORMAL:
            addedges(edges_min_234, vtcs_min, pts_min, 0, n_min, n_min - 2, 0, -1, -1, n_min,
                    rs);
            break;
        case WAGNER:
            switch (n_min) {
                case 5:
                    make_K_5_edges(edges_min_234, vtcs_min, n_min);
                    break;
                case 6:
                    make_K_33_edges(edges_min_234, vtcs_min, n_min);
                    break;
                default:;
            }
            break;
        default:;
    }

    /*
     * To create the base graph we need to replace all minor points by subgraphs.
     * Therefore we need to determine the areas in which we can place the subgraphs
     * which depend on the shortest distance between the minor points and their
     * neirest neighbours and the grid border respectively. For every minor point
     * the minimum of these two values is used to calculate the radius of a circular
     * subgraph area around the minor point.
     */
    tmp = COORDMARGIN * COORDUNIT;
    tmp2 = (coord_lim - COORDMARGIN) * COORDUNIT;
    radii = snewn(n_min, long);
    if (params->mode != WAGNER || n_min < 6)
    {
        for (i = 0; i < n_min; i++)
        {
            pt = pts_min + i;
            radii[i] = min(min(min(pt->x - tmp, pt->y - tmp), min(tmp2 - pt->x, tmp2 - pt->y)),
                            (squarert(square(pts_min[1].x - pts_min[0].x) + square(pts_min[1].y
                            - pts_min[0].y)) / 2)) - (SUBGRAPH_DISTANCE / 2);
            LOG(("Assigned subgraph radius %ld to subgraph %ld\n", radii[i], i));
        }
    }
    else
    {
        for (i = 0; i < n_min; i++)
        {
            pt = pts_min + i;
            radii[i] = min(min(pt->x - tmp, pt->y - tmp), min(tmp2 - pt->x, tmp2 - pt->y))
                        - (SUBGRAPH_DISTANCE / 2);
            LOG(("Assigned subgraph radius %ld to subgraph %ld\n", radii[i], i));
        }
        for (i = 0; i < n_min - 1; i++)
        {
            for (j = i + 1; j < n_min; j++)
            {
                long radius = (squarert(square(pts_min[i].x - pts_min[j].x) + square(pts_min[i].y
                                - pts_min[j].y)) / 2) - (SUBGRAPH_DISTANCE / 2);
                if (radius < radii[i])
                {
                    radii[i] = radius;
                    LOG(("Reassigned subgraph radius %ld to subgraph %ld\n", radii[i], i));
                }
                if (radius < radii[j])
                {
                    radii[j] = radius;
                    LOG(("Reassigned subgraph radius %ld to subgraph %ld\n", radii[j], j));
                }
            }
        }
    }

    /*
     * Assign coordinates to the subgraphs. The coordinates must lie in the previously
     * calculated subgraph areas.
     */
    tmp = n_sub * SUBGRAPH_POINTENTROPY;
    angles = snewn(tmp, double);
    for (i = 0; i < tmp; i++)
    {
        angles[i] = ((double) i * 2.0 * PI) / (double) tmp;
    }
    pts_base = snewn(n_base, point);
    for (i = 0; i < n_min; i++)
    {
        double rotation = (double) random_upto(rs, (100.0 * 2.0 * PI)
                            / (double) n_sub) / 100.0;
        shuffle(angles, tmp, sizeof(double), rs);
        for (j = 0; j < n_sub; j++)
        {
            pt = pts_base + (i * n_sub) + j;
            pt->x = (double) pts_min[i].x + ((double) radii[i] * sin(angles[j] + rotation));
            pt->y = (double) pts_min[i].y + ((double) radii[i] * cos(angles[j] + rotation));
            pt->d = 1;
            LOG(("Assigned coordinates x:%ld, y:%ld and denominator %ld to subgraph"\
                " point %ld of subgraph %ld (base graph point %ld)\n", pt->x, pt->y, pt->d,
                j, i, (i * n_sub) + j));
        }
    }
    sfree(radii);
    sfree(angles);

    /* Add edges to the subgraphs */
    vtcs_base = snewn(n_base, vertex);
    for (i = 0; i < n_base; i++)
    {
        vx = vtcs_base + i;
        vx->idx = i;
        vx->deg = 0;
    }
    edges_base_234 = newtree234(edgecmp);
    for (i = 0; i < n_min; i++)
    {
        addedges(edges_base_234, vtcs_base, pts_base,  i * n_sub, n_sub, 2 * n_sub / 3,
                i * n_sub, -1, -1, n_min * n_sub, rs);
    }

    /*
     * For every minor edge add an edge between random points of the subgraphs that
     * correspond to the adjacent minor points.
     */
    tmp = count234(edges_min_234);
    for (i = 0; (e = index234(edges_min_234, i)) != NULL; i++)
    {
        bool added = false;
        edge beste;
        edge_ext* eext;
        tree234* eexts = newtree234(eextcmp);
        beste.src = beste.tgt = -1;
        for (j = e->src * n_sub; j < (e->src + 1) * n_sub; j++)
        {
            for (k = e->tgt * n_sub; k < (e->tgt + 1) * n_sub; k++)
            {
                eext = snew(edge_ext);
                eext->e = (edge) { j, k };
                eext->degsum = vtcs_base[j].deg + vtcs_base[k].deg;
                add234(eexts, eext);
            }
        }
#if DEBUG
        for (j = 0; (eext = index234(eexts, j)) != NULL; j++)
        {
            LOG(("Added extended edge with source %d, target %d and degree"\
                "sum %d at position %ld to the 234-tree of extended edges\n",
                eext->e.src, eext->e.tgt, eext->degsum, j));
        }
#endif
        for (j = 0; (eext = index234(eexts, j)) != NULL; j++)
        {
            beste = eext->e;
            /* check for crossing points */
            for (l = 0; l < n_min * n_sub; l++)
            {
                if (l == beste.src || l == beste.tgt)
                    continue;
                else if (crosspoint(pts_base[beste.src], pts_base[beste.tgt], pts_base[l]))
                    goto next_subedge; /* this edge crosses a point => next target */
            }
            addedge(edges_base_234, beste.src, beste.tgt);
            added = true;
            goto next_minedge;
            next_subedge:;
        }
        next_minedge:
        if (!added)
        {
            beste = ((edge_ext*) index234(eexts, 0))->e;
            addedge(edges_base_234, beste.src, beste.tgt);
        }
        LOG(("Added edge between subgraphs %d and %d (base graph vertices %d and %d)\n",
            e->src, e->tgt, beste.src, beste.tgt));
        vtcs_base[beste.src].deg++;
        vtcs_base[beste.tgt].deg++;
        LOG(("Updated degrees of vertices %d to %d and %d to %d\n", vtcs_base[beste.src].idx,
            vtcs_base[beste.src].deg, vtcs_base[beste.tgt].idx, vtcs_base[beste.tgt].deg));
        while ((eext = delpos234(eexts, 0)) != NULL) sfree(eext);
        freetree234(eexts);
    }

    /*
     * Assign coordinates to the remaining points. The coordinates may not overlay
     * other points or edges.
     */
    tmp = coord_lim - (2 * COORDMARGIN) + 1;
    coords_x = snewn(tmp, long);
    coords_y = snewn(tmp, long);
    for (i = 0; i < tmp; i++)
    {
        coords_x[i] = (i + COORDMARGIN) * COORDUNIT;
        coords_y[i] = (i + COORDMARGIN) * COORDUNIT;
    }
    shuffle(coords_x, tmp, sizeof(long), rs);
    shuffle(coords_y, tmp, sizeof(long), rs);
    tmp2 = n_min * n_sub;
    tmp3 = 0;
    for (i = 0; i < tmp; i++)
    {
        point p;
        p.d = 1;
        p.x = coords_x[i];
        p.y = coords_y[i];
        /* check for a crossing with an edge */
        for (k = 0; (e = index234(edges_base_234, k)) != NULL; k++)
        {
            if (crosspoint(pts_base[e->src], pts_base[e->tgt], p))
                goto next_coords; /* an edge crosses the point => next coords */
        }
        /* check for an overlaying with a point */
        for (k = 0; k < tmp2 + tmp3; k++)
        {
            pt = pts_base + k;
            if (square(pt->x - p.x) + square(pt->y - p.y) < OVERLAYPOINT_TRESHOLD_GAMEGEN)
                goto next_coords; /* the point overlays another point => next coords */
        }
        pt = pts_base + tmp2 + tmp3;
        *pt = p;
        LOG(("Assigned coordinates x:%ld, y:%ld and denominator %ld to base graph point"\
            " %ld\n", pt->x, pt->y, pt->d, tmp2 + tmp3));
        if (tmp2 + ++tmp3 >= n_base) break;
        next_coords:;
    }
    sfree(coords_x);
    sfree(coords_y);

    /* Add edges to the remaining points */
    addedges(edges_base_234, vtcs_base, pts_base, 0, n_base, n_base - 1, n_min * n_sub, -1, 3,
            n_base, rs);
    /* Add more edges to confuse the player */
    addedges(edges_base_234, vtcs_base, pts_base, 0, n_base, 5, 0, -1, -1, n_base, rs);

#if DEBUG
    LOG(("Base graph vertices have  index - degree "));
    for (i = 0; i < n_base; i++) {
        vx = vtcs_base + i;
        LOG(("%d-%d, ", vx->idx, vx->deg));
        assert(vx->deg > 0);
    }
    LOG(("\n"));
#endif

    /*
     * The generation of a new game description is finished. Now we need to encode
     * the description in a dynamically allocated string and connect this string to
     * the return value.
     */
    ret = NULL;
    {
    const char* sep = ",";
    char buf[80];
    int len = 0;
    int off = 0;
    long count_min = count234(edges_min_234);
    long count_base = count234(edges_base_234);
    edge* edges_min = snewn(count_min, edge);
    edge* edges_base = snewn(count_base, edge);

    /*
     * Calculate the length of the game description. It contains information about
     * the points and edges of the minor and base graph. Points and edges are encoded
     * in the following way:
     * 
     * (1) point: <index> - <degree> - <x coordinate> - <y coordinate>
     * (2) edge: <source index> - <target index>
     * 
     * Single points or edges are separated by a comma while sets of points or edges
     * are separated by a semicolon:
     * 
     * (1),(1),...;(2),(2),...;(1),(1),...;(2),(2);
     */
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        len += (sprintf(buf, "%d-%d-%ld-%ld", vtcs_min[i].idx, vtcs_min[i].deg,
                        pt->x, pt->y) + 1);
    }
    for (i = 0; (e = index234(edges_min_234, i)) != NULL; i++)
    {
        edges_min[i] = *e;
        len += (sprintf(buf, "%d-%d", e->src, e->tgt) + 1);
    }
    for (i = 0; i < n_base; i++)
    {
        pt = pts_base + i;
        len += (sprintf(buf, "%d-%d-%ld-%ld", vtcs_base[i].idx, vtcs_base[i].deg,
                        pt->x, pt->y) + 1);
    }
    for (i = 0; (e = index234(edges_base_234, i)) != NULL; i++)
    {
        edges_base[i] = *e;
        len += (sprintf(buf, "%d-%d", e->src, e->tgt) + 1);
    }

    /*
     * Allocate memory for len+1 chars, that is exactly the length of the game
     * description including a trailing '\0'.
     */
    ret = snewn(++len, char);
    
    /*
     * Now encode the game description and write it into the allocated string
     * that will be connected to the return value.
     */
    for (i = 0; i < n_min - 1; i++)
    {
        pt = pts_min + vtcs_min[i].idx;
        off += sprintf(ret + off, "%d-%d-%ld-%ld%s", vtcs_min[i].idx, vtcs_min[i].deg,
                        pt->x, pt->y, sep);
    }
    sep = ";";
    pt = pts_min + vtcs_min[n_min-1].idx;
    off += sprintf(ret + off, "%d-%d-%ld-%ld%s", vtcs_min[n_min-1].idx, vtcs_min[n_min-1].deg,
                    pt->x, pt->y, sep);
    sep = ",";
    for (i = 0; i < count_min - 1; i++)
    {
        e = edges_min + i;
        off += sprintf(ret + off, "%d-%d%s",e->src, e->tgt, sep);
    }
    sep = ";";
    e = edges_min + count_min - 1;
    off += sprintf(ret + off, "%d-%d%s", e->src, e->tgt, sep);
    sep = ",";
    for (i = 0; i < n_base - 1; i++)
    {
        pt = pts_base + vtcs_base[i].idx;
        off += sprintf(ret + off, "%d-%d-%ld-%ld%s", vtcs_base[i].idx, vtcs_base[i].deg,
                        pt->x, pt->y, sep);
    }
    sep = ";";
    pt = pts_base + vtcs_base[n_base-1].idx;
    off += sprintf(ret + off, "%d-%d-%ld-%ld%s", vtcs_base[n_base-1].idx, vtcs_base[n_base-1].deg,
                    pt->x, pt->y, sep);
    sep = ",";
    for (i = 0; i < count_base - 1; i++)
    {
        e = edges_base + i;
        off += sprintf(ret + off, "%d-%d%s", e->src, e->tgt, sep);
    }
    sep = ";";
    e = edges_base + count_base - 1;
    sprintf(ret + off, "%d-%d%s", e->src, e->tgt, sep);

    sfree(edges_min);
    sfree(edges_base);
    }

    /* The aux string is not required and therefore it is set to NULL */
    *aux = NULL;

    sfree(vtcs_min);
    sfree(pts_min);
    while ((e = delpos234(edges_min_234, 0)) != NULL) sfree(e);
    freetree234(edges_min_234);

    sfree(vtcs_base);
#ifndef BENCHMARKS
    sfree(pts_base);
    while ((e = delpos234(edges_base_234, 0)) != NULL) sfree(e);
    freetree234(edges_base_234);
#else
    calc_basegraph_dissim(edges_base_234, pts_base, n_base, n_min);
    if (GraphDissim)
        printf("Base graph dissimiliarity between the last two game generations is %d\n",
                GraphDissim);
#endif

    return ret;
}

/*
 * Validate a graph description, i.e. a string that specifies a set of points by
 * their index position and degree and a set of edges between those points by the
 * point indices of their incident vertices.
 */
static const char* validate_graph(const char** desc, int n, long lim, long mar)
{
    int idx;
    int src, tgt;
    long x, y;
    while (**desc)
    {
        idx = atoi(*desc);
        if(idx < 0 || idx >= n)
            return "Point index out of range in game description";
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after point index in game description";
        (*desc)++;

        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after point degree in game description";
        (*desc)++;

        x = atol(*desc);
        if (x < mar || x > lim)
            return "X-coordinate out of range in game description";
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after point index in game description";
        (*desc)++;

        y = atol(*desc);
        if (y < mar || y > lim)
            return "Y-coordinate out of range in game description";
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        if (**desc != ',' && **desc != ';')
            return "Expected ',' or ';' after y-coordinate in game description";
        if (*((*desc)++) == ';') break;
    }
    while (**desc)
    {
        src = atoi(*desc);
        if (src < 0 || src >= n)
            return "Edge source index out of range in game description";
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after edge source index in game description";
        (*desc)++;

        tgt = atoi(*desc);
        if (tgt < 0 || tgt >= n)
            return "Edge target index out of range in game description";
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        if (**desc != ',' && **desc != ';')
            return "Expected ',' or ';' after edge target in game description";
        if (*((*desc)++) == ';') break;
    }    

    return NULL;
}

static const char *validate_desc(const game_params *params, const char *desc)
{
    const char* _desc = desc; /* pointer copy */
    const char* err;
    long coord_lim = COORDLIMIT(params->n_base) * COORDUNIT;
    long coord_mar = COORDMARGIN * COORDUNIT;
    if ((err = validate_graph(&_desc, params->n_min, coord_lim - coord_mar, coord_mar)) != NULL)
        return err;
    else if ((err = validate_graph(&_desc, params->n_base, coord_lim - coord_mar, coord_mar)) != NULL)
        return err;
    else
        return NULL;
}

/*
 * Parse a graph description, i.e. a string that specifies a set of points by their
 * index, position and degree and a set of edges between those points by the point
 * indices of their incident vertices.
 */
static graph* parse_graph(const char** desc, int n, long lim, long mar)
{
    int idx, deg;
    int src, tgt;
    long x, y;
    point* pt;
    vertex* vx;
    graph* ret = snew(graph);
    ret->refcount = 1;
    ret->points = snewn(n, point);
    ret->vtcs = snewn(n, vertex);
    ret->vertices = newtree234(vertcmp);
    ret->edges = newtree234(edgecmp);
    do
    {
        idx = atoi(*desc);
#if DEBUG
        assert(idx >= 0 && idx < n);
#endif
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

#if DEBUG
        assert(**desc == '-');
#endif
        (*desc)++;

        deg = atoi(*desc);
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

#if DEBUG
        assert(**desc == '-');
#endif
        (*desc)++;

        x = atol(*desc);
#if DEBUG
        assert(x >= mar && x <= lim);
#endif
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

#if DEBUG
        assert(**desc == '-');
#endif
        (*desc)++;

        y = atol(*desc);
#if DEBUG
        assert(y >= mar && y <= lim);
#endif
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        pt = ret->points + idx;
        pt->x = x;
        pt->y = y;
        pt->d = 1;

        vx = ret->vtcs + idx;
        vx->idx = idx;
        vx->deg = deg;
        add234(ret->vertices, vx);

#if DEBUG
        assert(**desc == ',' || **desc == ';');
#endif
    }
    while (*((*desc)++) != ';');
    do
    {
        src = atoi(*desc);
#if DEBUG
        assert(src >= 0 && src < n);
#endif
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

#if DEBUG
        assert(**desc == '-');
#endif
        (*desc)++;

        tgt = atoi(*desc);
#if DEBUG
        assert(tgt >= 0 && tgt < n);
#endif
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        addedge(ret->edges, src, tgt);

#if DEBUG
        assert(**desc == ',' || **desc == ';');
#endif
    }
    while (*((*desc)++) != ';');

    return ret;
}

static game_state *new_game(midend *me, const game_params *params,
                            const char *desc)
{
    const char* _desc = desc; /* pointer copy */
    game_state *state = snew(game_state);
    state->params = *params;
    long coord_lim = COORDLIMIT(params->n_base) * COORDUNIT;
    long coord_mar = COORDMARGIN * COORDUNIT;
    state->minor = parse_graph(&_desc, params->n_min, coord_lim - coord_mar, coord_mar);
    state->base = parse_graph(&_desc, params->n_base, coord_lim - coord_mar, coord_mar);

    grid minor_grid = {
        0,
        0,
        (float) COORDLIMIT(params->n_min) / (float) COORDLIMIT(params->n_base)
    };
    state->minor->grid = minor_grid;

    grid base_grid = {
        minor_grid.x_off + minor_grid.relsize * coord_lim,
        0,
        1.0F
    };
    state->base->grid = base_grid;

    state->solved = false;
    state->cheated = false;

    return state;
}

/*
 * Copy a single edge
 */
static void* edgecpy(void* state, void* elem)
{
    edge* e = (edge*) elem;
    edge* ecpy = snew(edge);
    *ecpy = *e;
    
    return (void*) ecpy;
}

/*
 * Duplicate a graph structure. The duplicates refcount will be 1.
 */
static graph* dup_graph(const graph* gr, int n)
{
    int i;
    vertex* vx;
    graph* ret = snew(graph);

    ret->refcount = 1;
    ret->grid = gr->grid;
    ret->points = snewn(n, point);
    memcpy(ret->points, gr->points, n * sizeof(point));
    ret->vtcs = snewn(n, vertex);
    memcpy(ret->vtcs, gr->vtcs, n * sizeof(vertex));

    ret->vertices = newtree234(vertcmp);
    for (i = 0; (vx = index234(gr->vertices, i)) != NULL; i++)
    {
        add234(ret->vertices, ret->vtcs + vx->idx);
#if DEBUG
    assert(ret->vtcs[vx->idx].idx == vx->idx);
#endif
    }

    ret->edges = copytree234(gr->edges, edgecpy, NULL);

    return ret;
}

static game_state *dup_game(const game_state *state)
{
    game_state *ret = snew(game_state);

    ret->params = state->params;
    ret->minor = state->minor;
    ret->minor->refcount++;
    ret->base = dup_graph(state->base, state->params.n_base);
    ret->solved = state->solved;
    ret->cheated = state->cheated;

    return ret;
}

/*
 * Free the memory that points to a graph structure if its refcount is equal to or
 * smaller than 0.
 */
static void free_graph(graph* gr, int n)
{
    edge* e;
    (gr->refcount)--;
    if (gr->refcount <= 0)
    {
        sfree(gr->points);
        sfree(gr->vtcs);
        freetree234(gr->vertices);
        while((e = delpos234(gr->edges, 0)) != NULL) sfree(e);
        freetree234(gr->edges);
        sfree(gr);
    }
}

static void free_game(game_state *state)
{
    free_graph(state->base, state->params.n_base);
    free_graph(state->minor, state->params.n_min);
    sfree(state);
}

/*
 * Replace the given edge in the given edge 234-tree. The new edge will have new_src
 * and new_tgt as source and target respetively.
 */
static void replace_edge(tree234* edges, edge* e, tree234* vertices, vertex* vtcs,
                        int new_src, int new_tgt)
{
    del234(edges, e);
    del234(vertices, vtcs + e->src);
    del234(vertices, vtcs + e->tgt);
    vtcs[e->src].deg--;
    vtcs[e->tgt].deg--;
    add234(vertices, vtcs + e->src);
    add234(vertices, vtcs + e->tgt);

    if (new_src != new_tgt && !isedge(edges, new_src, new_tgt))
    {
        e->src = min(new_src, new_tgt);
        e->tgt = max(new_src, new_tgt);
        add234(edges, e);
        del234(vertices, vtcs + e->src);
        del234(vertices, vtcs + e->tgt);
        vtcs[e->src].deg++;
        vtcs[e->tgt].deg++;
        add234(vertices, vtcs + e->src);
        add234(vertices, vtcs + e->tgt);
    }
    else
    {
        sfree(e);
    }
}

/* 
 * Contract an edge from a graph, i.e. merge its incident vertices such that no edges
 * are lost except for the contracted edge itself.
 */
static void contract_edge(graph* graph, int dom, int rec)
{
    int i;
    edge* e;
    edge* ecpy;
    tree234* edgescpy = copytree234(graph->edges, edgecpy, NULL);
     
    for (i = 0; (e = index234(graph->edges, i)) != NULL; i++)
    {
        ecpy = find234(edgescpy, e, edgecmp);
        if (rec == ecpy->src)
            replace_edge(edgescpy, ecpy, graph->vertices, graph->vtcs,
                        dom, ecpy->tgt);
        else if (rec == ecpy->tgt)
            replace_edge(edgescpy, ecpy, graph->vertices, graph->vtcs,
                        ecpy->src, dom);
    }

    del234(graph->vertices, graph->vtcs + rec);
#if DEBUG
    assert(graph->vtcs[rec].idx == rec);
#endif
    while((e = delpos234(graph->edges, 0)) != NULL) sfree(e);
    freetree234(graph->edges);
    graph->edges = edgescpy;
}

/*
 * Delete an edge from a graph
 */
static void delete_edge(graph* graph, edge e)
{
    edge* _e;

    /* delete the edge */
    _e = del234(graph->edges, &e);
    sfree(_e);

    /* update the involved vertices */
    del234(graph->vertices, graph->vtcs + e.src);
    del234(graph->vertices, graph->vtcs + e.tgt);
    graph->vtcs[e.src].deg--;
    graph->vtcs[e.tgt].deg--;
    add234(graph->vertices, graph->vtcs + e.src);
    add234(graph->vertices, graph->vtcs + e.tgt);
}

/*
 * Delete a vertex from a graph
 */
static void delete_vertex(graph* graph, vertex v)
{
    int i;
    edge* e;
    tree234* edgescpy = copytree234(graph->edges, edgecpy, NULL);

    for (i = 0; (e = index234(edgescpy, i)) != NULL; i++)
    {
        /* delete all edges that are incident to the vertex */
        if (e->src == v.idx || e->tgt == v.idx)
            delete_edge(graph, *e);
    }

    /* delete the vertex */
    del234(graph->vertices, &graph->vtcs[v.idx]);

    while((e = delpos234(edgescpy, 0)) != NULL) sfree(e);
    freetree234(edgescpy);
}

/*
 * A node of a graph isomorphism search tree
 */
typedef struct node node;

struct node {
    node* children;
    vertex** cells;
    int* cellsizes;
    int ncells;
    int nchildren;
    bool isleaf;
};

/*
 * Expand a node in an isomorphism search tree or if it is a leaf parse it into a
 * permutation that is candidate to be an isomorphism between two graphs.
 */
static vertex* expand_node(node* n)
{
    int i, j;
    
    for (i = 0; i < n->ncells; i++)
    {
        if (n->cellsizes[i] > 1)
        {
            n->isleaf = false;
            n->nchildren = n->cellsizes[i];
            break;
        }
    }
    if (n->isleaf)
    {
        vertex* permu = snewn(n->ncells, vertex);
        LOG(("Reached leaf and created corresponding permutation with size %d\n"\
            "The permutation is (", n->ncells));
        for (i = 0; i < n->ncells; i++)
        {
            permu[i] = *n->cells[i];
            LOG(("%d%s", permu[i].idx, (i < n->ncells - 1) ? ", " : ")\n"));
        }
        return permu;
    }
    else
    {
        node* child;
        n->children = snewn(n->nchildren, node);
        LOG(("Expanded node and created %d child nodes\n",
            n->nchildren));
        for (i = 0; i < n->nchildren; i++)
        {
            child = n->children + i;
            child->ncells = n->ncells + 1;
            child->cells = snewn(child->ncells, vertex*);
            child->cellsizes = snewn(child->ncells, int);
            child->isleaf = true;
            LOG(("Initialized child node %d with %d cells\nCreating new refinement"\
                " of parent cells\n", i, child->ncells));
            j = 0;
            while (n->cellsizes[j] == 1)
            {
                child->cells[j] = snew(vertex);
                *child->cells[j] = *n->cells[j];
                child->cellsizes[j++] = 1;
                LOG(("Initialized child cell %d with size 1\n", j - 1));
            }
            child->cells[j] = snew(vertex);
            *child->cells[j] = n->cells[j][i];
            child->cellsizes[j++] = 1;
            LOG(("Initialized child cell %d with size 1\n", j - 1));
            child->cellsizes[j] = n->cellsizes[j-1] - 1;
            child->cells[j] = snewn(child->cellsizes[j], vertex);
            if (i > 0)
                memcpy(child->cells[j], n->cells[j-1], i * sizeof(vertex));
            if (i < n->cellsizes[j-1] - 1)
                memcpy(child->cells[j] + i, n->cells[j-1] + i + 1,
                        (n->cellsizes[j-1] - 1 - i) * sizeof(vertex));
            LOG(("Initialized child cell %d with size %d\n", j, child->cellsizes[j]));
            for (j++; j < child->ncells; j++)
            {
                child->cells[j] = snewn(n->cellsizes[j-1], vertex);
                memcpy(child->cells[j], n->cells[j-1],
                        n->cellsizes[j-1] * sizeof(vertex));
                child->cellsizes[j] = n->cellsizes[j-1];
                LOG(("Initialized child cell %d with size %d\n", j, child->cellsizes[j]));
            }
            LOG(("Created new refinement of parent cells for child %d\n", i));
        }
        return NULL;
    }
}

/*
 * Free a node structure and all its child nodes recursively. Do only free the node
 * itself if it is the root node of an isomorphism search tree.
 */
static void free_node(node* n, bool isroot)
{
    int i;

    if (!n->isleaf)
    {
        for (i = 0; i < n->nchildren; i++)
        {
            free_node(n->children + i, false);
        }
        sfree(n->children);
    }
    for (i = 0; i < n->ncells; i++) sfree(n->cells[i]);
    sfree(n->cells);
    sfree(n->cellsizes);
    if (isroot) sfree(n);
}

/*
 * A key value mapping in a map
 */
typedef struct mapping {
    int key;
    int value;
} mapping;

/*
 * Compare a key value mapping by its key
 */
static int mappingcmpC(const void* av, const void* bv)
{
    const mapping* a = (mapping*) av;
    const mapping* b = (mapping*) bv;

    if (a->key < b->key) return -1;
    else if (a->key > b->key) return 1;
    else return 0;
}

static int mappingcmp(void* av, void* bv)
{
    return mappingcmpC(av, bv);
}

/*
 * Copy a single mapping
 */
static void* mappingcpy(void* state, void* elem)
{
    mapping* m = (mapping*) elem;
    mapping* mcpy = snew(mapping);
    *mcpy = *m;
    
    return (void*) mcpy;
}

/*
 * Check whether there exists an isomorphism between a test graph and a comparison
 * graph. The algorithm builds on the idea that vertices with different degree can't
 * be a vertex pair of a permutation that is an isomorphism between two graphs.
 */
static bool isomorphism_degheuristic(const graph* the_graph, const graph* cmp_graph,
                                    tree234** solution)
{
    bool found;
    int i, j;
    int tmp;
    int nvtcs = count234(the_graph->vertices);
    vertex* vx;
    vertex* vtcs_the;
    vertex* vtcs_cmp;
    vertex* permu;
    edge* e;
    node* n;
    node* root;
    tree234* lifo;

#ifdef BENCHMARKS
    clock_t begin = clock();
#endif
    
    /* Check whether the graphs have the same number of vertices */
    if (nvtcs != (tmp = count234(cmp_graph->vertices)))
    {
        LOG(("The graphs have different amounts of vertices (%d and %d)\n",
            nvtcs, tmp));
        return false;
    }

    /* Check whether the graphs have identical vertex degree distributions */
    vtcs_the = snewn(nvtcs, vertex);
    for (i = 0; (vx = index234(the_graph->vertices, i)) != NULL; i++)
    {
        vtcs_the[i] = *vx;
#if DEBUG
        if (i > 0) assert(vtcs_the[i-1].deg <= vtcs_the[i].deg);
#endif
    }
    vtcs_cmp = snewn(nvtcs, vertex);
    for (i = 0; (vx = index234(cmp_graph->vertices, i)) != NULL; i++)
    {
        vtcs_cmp[i] = *vx;
#if DEBUG
        if (i > 0) assert(vtcs_cmp[i-1].deg <= vtcs_cmp[i].deg);
#endif
    }

    tmp = 1;
    for (i = 0; i < nvtcs; i++)
    {
        LOG(("Vertices at position %d have indices %d and %d and degrees %d and %d\n",
            i, vtcs_the[i].idx, vtcs_cmp[i].idx, vtcs_the[i].deg, vtcs_cmp[i].deg));
        if (vtcs_the[i].deg != vtcs_cmp[i].deg)
        {
            LOG(("The graphs have different vertex degree distributions\n"));
            sfree(vtcs_the);
            sfree(vtcs_cmp);
            return false;
        }
        else if (i < nvtcs - 1 && vtcs_the[i+1].deg != vtcs_the[i].deg)
            tmp++;
    }

    /*
     * Initiate a search tree by adding vertices with the same number of incident
     * edges to the same cell of a node. Each vertex is unique in the union of all
     * cells of a node and the union of all cells of a node must contain all vertices
     * of the graph for which we create the search tree - the comparison graph.
     * Then we refine for each node the cells that contain more than one vertex.
     * Each refinement gives us another set of child nodes. We reach a leaf when
     * we refine a node such that all cells of the arising child do only contain
     * a single vertex. This means we just found a permutation for the comparison
     * graph which could suit the test graph. If the permutation suits the test graph,
     * i.e. if we apply it to the test graph and it turns out that the edge sets match,
     * then we've found an isomorphism between the graphs.
     */
    root = snew(node);
    root->ncells = tmp;
    root->cells = snewn(root->ncells, vertex*);
    root->cellsizes = snewn(root->ncells, int);
    root->isleaf = true;
    LOG(("Initialized root node with %d cells\n", root->ncells));
    for (i = 0; i < root->ncells; i++)
    {
        root->cellsizes[i] = 1;
    }
    tmp = 0;
    for (i = 1; i < nvtcs; i++)
    {
        if (vtcs_the[i].deg == vtcs_the[i-1].deg)
            root->cellsizes[tmp]++;
        else
            tmp++;
    }
#if DEBUG
    assert(tmp == root->ncells - 1);
#endif
    tmp = 0;
    for (i = 0; i < root->ncells; i++)
    {
        root->cells[i] = snewn(root->cellsizes[i], vertex);
        for (j = 0; j < root->cellsizes[i]; j++)
        {
            vx = root->cells[i] + j;
            *vx = vtcs_cmp[tmp+j];
            LOG(("Added vertex %d with index %d and degree %d to root cell %d at"\
                " position %d\n", tmp + j, vtcs_cmp[tmp+j].idx, vtcs_cmp[tmp+j].deg,
                i, j));
        }
        tmp += root->cellsizes[i];
    }
    sfree(vtcs_cmp);
    tmp = 0;
    lifo = newtree234(NULL);
    addpos234(lifo, root, tmp++);
    LOG(("Initialized lifo queue with root node\nStarting depth first search for"\
        " an isomorphism between the graphs\n"));
    while (tmp)
    {
        n = delpos234(lifo, --tmp);
        LOG(("Fetched node at position %d from the lifo queue\n", tmp));
        if((permu = expand_node(n)))
        {
            mapping* mappings;
            tree234* map = newtree234(mappingcmp);
            found = true;
            LOG(("Found vertex permuation that could be an ismomorphism"\
                " between the graphs\n"));
            /*
             * We don't need to check whether the test graph has additional edges that
             * are missing in the comparison graph since we already compared the vertex
             * degree distibutions of both graphs. Thus the sum of all vertex degrees
             * must be equal for both graphs and the number of edges respectively.
             */
            mappings = snewn(nvtcs, mapping);
            for (j = 0; j < nvtcs; j++)
            {
                mappings[j].key = permu[j].idx;
                mappings[j].value = vtcs_the[j].idx;
                add234(map, mappings + j);
            }
            sfree(permu);
            for (j = 0; (e = index234(cmp_graph->edges, j)) != NULL; j++)
            {
                mapping msrc, mtgt;
                msrc.key = e->src;
                mtgt.key = e->tgt;
                msrc.value = mtgt.value = -1;
                msrc.value = ((mapping*) find234(map, &msrc, NULL))->value;
                mtgt.value = ((mapping*) find234(map, &mtgt, NULL))->value;
                if (!isedge(the_graph->edges, msrc.value, mtgt.value))
                {
                    found = false;
                    LOG(("Permutation is no isomorphism between the graphs,\n"\
                        " missing edge %d-%d in the test graph\n", msrc.value,
                        mtgt.value));
                    sfree(mappings);
                    freetree234(map);
                    break;
                }
            }
            if (found)
            {
                if (solution) *solution = copytree234(map, mappingcpy, NULL);
                LOG(("Permutation is an isomorphism between the graphs\n"));
                sfree(mappings);
                freetree234(map);
                break;
            }
        }
        else
        {
            for (i = n->nchildren - 1; i >= 0; i--)
            {
                addpos234(lifo, n->children + i, tmp++);
                LOG(("Added child node %d at position %d to lifo queue\n", i,
                    tmp - 1));
            }
        }
    }

#ifdef BENCHMARKS
    clock_t end = clock();
    double duration = ((double) (end - begin) * 1000.0) / CLOCKS_PER_SEC;
    printf("Finished isomorphism test with degree heuristic, duration: %lf\n",
        duration);
#endif

    sfree(vtcs_the);
    free_node(root, true);
    freetree234(lifo);

    return found;
}

#ifdef BENCHMARKS
/*
 * Check whether there exists an isomorphism between a test graph and a comparison
 * graph. The algorithm is a simple bruteforce algorithm that checks for all possible
 * permutations whether it is an isomorphism between the graphs.
 */
static bool isomorphism_bruteforce(const graph* the_graph, const graph* cmp_graph)
{
    bool found;
    int i, j;
    int tmp;
    int nvtcs = count234(the_graph->vertices);
    vertex* vx;
    vertex* vtcs_the;
    vertex* vtcs_cmp;
    vertex* permu;
    edge* e;
    node* n;
    node* root;
    tree234* lifo;

    clock_t begin = clock();
    
    /* Check whether the graphs have the same number of vertices */
    if (nvtcs != (tmp = count234(cmp_graph->vertices)))
    {
        LOG(("The graphs have different amounts of vertices (%d and %d)\n",
            nvtcs, tmp));
        return false;
    }

    /* Check whether the graphs have identical vertex degree distributions */
    vtcs_the = snewn(nvtcs, vertex);
    for (i = 0; (vx = index234(the_graph->vertices, i)) != NULL; i++)
    {
        vtcs_the[i] = *vx;
#if DEBUG
        if (i > 0) assert(vtcs_the[i-1].deg <= vtcs_the[i].deg);
#endif
    }
    vtcs_cmp = snewn(nvtcs, vertex);
    for (i = 0; (vx = index234(cmp_graph->vertices, i)) != NULL; i++)
    {
        vtcs_cmp[i] = *vx;
#if DEBUG
        if (i > 0) assert(vtcs_cmp[i-1].deg <= vtcs_cmp[i].deg);
#endif
    }

    for (i = 0; i < nvtcs; i++)
    {
        LOG(("Vertices at position %d have indices %d and %d and degrees %d and %d\n",
            i, vtcs_the[i].idx, vtcs_cmp[i].idx, vtcs_the[i].deg, vtcs_cmp[i].deg));
        if (vtcs_the[i].deg != vtcs_cmp[i].deg)
        {
            LOG(("The graphs have different vertex degree distributions\n"));
            sfree(vtcs_the);
            sfree(vtcs_cmp);
            return false;
        }
    }

    root = snew(node);
    root->ncells = 1;
    root->cells = snewn(root->ncells, vertex*);
    root->cellsizes = snewn(root->ncells, int);
    root->isleaf = true;
    LOG(("Initialized root node with %d cells\n", root->ncells));
    *root->cellsizes = nvtcs;
    *root->cells = snewn(*root->cellsizes, vertex);
    for (j = 0; j < *root->cellsizes; j++)
    {
        vx = (*root->cells) + j;
        *vx = vtcs_cmp[j];
        LOG(("Added vertex %d with index %d and degree %d to root cell 0 at"\
            " position %d\n", j, vtcs_cmp[j].idx, vtcs_cmp[j].deg, j));
    }
    sfree(vtcs_cmp);
    tmp = 0;
    lifo = newtree234(NULL);
    addpos234(lifo, root, tmp++);
    LOG(("Initialized lifo queue with root node\nStarting depth first search for"\
        " an isomorphism between the graphs\n"));
    while (tmp)
    {
        n = delpos234(lifo, --tmp);
        LOG(("Fetched node at position %d from the lifo queue\n", tmp));
        if((permu = expand_node(n)))
        {
            mapping* mappings;
            tree234* map = newtree234(mappingcmp);
            found = true;
            LOG(("Found vertex permuation that could be an ismomorphism"\
                " between the graphs\n"));
            /*
             * We don't need to check whether the test graph has additional edges that
             * are missing in the comparison graph since we already compared the vertex
             * degree distibutions of both graphs. Thus the sum of all vertex degrees
             * must be equal for both graphs and the number of edges respectively.
             */
            mappings = snewn(nvtcs, mapping);
            for (j = 0; j < nvtcs; j++)
            {
                mappings[j].key = permu[j].idx;
                mappings[j].value = vtcs_the[j].idx;
                add234(map, mappings + j);
            }
            sfree(permu);
            for (j = 0; (e = index234(cmp_graph->edges, j)) != NULL; j++)
            {
                mapping msrc, mtgt;
                msrc.key = e->src;
                mtgt.key = e->tgt;
                msrc.value = mtgt.value = -1;
                msrc.value = ((mapping*) find234(map, &msrc, NULL))->value;
                mtgt.value = ((mapping*) find234(map, &mtgt, NULL))->value;
                if (!isedge(the_graph->edges, msrc.value, mtgt.value))
                {
                    found = false;
                    LOG(("Permutation is no isomorphism between the graphs,\n"\
                        " missing edge %d-%d in the test graph\n", msrc.value,
                        mtgt.value));
                    break;
                }
            }
            sfree(mappings);
            freetree234(map);
            if (found)
            {
                LOG(("Permutation is an isomorphism between the graphs\n"));
                break;
            }
        }
        else
        {
            for (i = n->nchildren - 1; i >= 0; i--)
            {
                addpos234(lifo, n->children + i, tmp++);
                LOG(("Added child node %d at position %d to lifo queue\n", i,
                    tmp - 1));
            }
        }
    }

    clock_t end = clock();
    double duration = ((double) (end - begin) * 1000.0) / CLOCKS_PER_SEC;
    printf("Finished bruteforce isomorphism test, duration: %lf\n", duration);

    sfree(vtcs_the);
    free_node(root, true);
    freetree234(lifo);

    return found;
}
#endif

/*
 * A move that can be performed by a player
 */
enum move {
    MOVE_IDLE = 0x0,
    MOVE_DRAGPOINT = 0x1,
    MOVE_CONTREDGE = 0x2,
    MOVE_DELPOINT = 0x4,
    MOVE_DELEDGE = 0x8
};

/*
 * Start from currstate and recursively check for all possible contraction sequences whether
 * they lead to a solved state or not. If yes the method should return a non-NULL pointer to
 * a dynamic string, otherwise it should return NULL. If the algorithm runs longer than time-
 * out milliseconds without finding a solution it will return NULL either.
 */
static char* solve_bruteforce(const game_state* currstate, game_state** solvedstate, tree234** solution,
                                int* movessize, int* moveslen, clock_t begin, int timeout)
{
    int i;
    edge* e;

    if (((double) (clock() - begin) * 1000.0) / CLOCKS_PER_SEC > timeout)
        return NULL;

    if (count234(currstate->base->vertices) > currstate->params.n_min)
    {
        LOG(("Bruteforce solver - Base graph has more vertices than minor graph left\n"));
        for (i = 0; (e = index234(currstate->base->edges, i)) != NULL; i++)
        {
            char* moves;
            game_state* nextstate = dup_game(currstate);
            contract_edge(nextstate->base, e->src, e->tgt);
            LOG(("Bruteforce solver - Contracted edge %d-%d\n", e->src, e->tgt));
            if ((moves = solve_bruteforce(nextstate, solvedstate, solution, movessize, moveslen, begin, timeout)))
            {
                char buf[80];
                char* oldmoves = NULL;
                int movesoff;
                free_game(nextstate);
                LOG(("Bruteforce solver - Current moveslen: %d, current movessize: %d\n", *moveslen, *movessize));
                LOG(("Bruteforce solver - First character pointed to by moves, should be 'S': %c\n", *moves));
                LOG(("Bruteforce solver - Encode last contraction %d-%d\n", e->src, e->tgt));
                /*
                 * Since we start encoding the very last move first, we always have to put
                 * subsequent move encodings in front of the previous encodings. This should
                 * ensure that our moves happen to be in the correct order when the most outter
                 * recursive call returns.
                 */
                if (*moveslen > 1)
                {
                    oldmoves = dupstr(moves+1);
                    LOG(("Bruteforce solver - oldmoves: %s\n", oldmoves));
                }
                movesoff = sprintf(buf, "%d:%d-%d;", MOVE_CONTREDGE, e->src, e->tgt);
                if ((*moveslen) + movesoff >= *movessize)
                {
                    *movessize = (*moveslen) + movesoff + 256;
                    moves = sresize(moves, *movessize, char);
                }
                strcpy(moves+1, buf);
                if (oldmoves)
                    strcpy(moves + movesoff + 1, oldmoves);
                sfree(oldmoves);
                LOG(("Bruteforce solver - newmoves: %s\n", moves));
                (*moveslen) += movesoff;
                return moves;
            }

            free_game(nextstate);
        }  
    }
    else if (isomorphism_degheuristic(currstate->minor, currstate->base, solution))
    {
        char* moves = snewn(*movessize, char);
        LOG(("Bruteforce solver - Found solution\n"));
#ifdef BENCHMARKS
        assert(isomorphism_bruteforce(currstate->minor, currstate->base));
#endif
        *moves = 'S';
        *moveslen = 1;
        *solvedstate = dup_game(currstate);
        return moves;
    }
    
    return NULL;
}

static char *solve_game(const game_state *state, const game_state *currstate,
                        const char *aux, const char **error)
{
    char buf[80];
    char* ret;

    int retsize = 256;
    int retlen = 0;
    int retoff;

    int n_base = currstate->params.n_base;
    int n_min = currstate->params.n_min;
    int n_1sub = n_base / n_min;
    int n_nsub = n_min * n_1sub;

    int i;

    point* pt;
    vertex* vx;
    mapping* m;

    game_state* solved;
    tree234* solution;

    /*
     * Delete all edges that are incident to the remaining points and then delete
     * the remaining points themselves.
     */
    ret = snewn(retsize, char);

    if (currstate->solved)
    {
        retoff = sprintf(buf, "%d:;", MOVE_IDLE);
        strcpy(ret + retlen, buf);
        retlen += retoff;

        LOG(("Idled\n"));

        return ret;
    }

    ret[retlen++] = 'S';

    solved = dup_game(currstate);
    for (i = n_nsub; i < n_base; i++)
    {
        vx = solved->base->vtcs + i;
        if (find234(solved->base->vertices, vx, NULL) != NULL)
        {
            retoff = sprintf(buf, "%d:%d;", MOVE_DELPOINT, vx->idx);
            if (retlen + retoff >= retsize)
            {
                retsize = retlen + retoff + 256;
                ret = sresize(ret, retsize, char);
            }
            strcpy(ret + retlen, buf);
            retlen += retoff;

            LOG(("Number of visible vertices is %d\n", count234(solved->base->vertices)));
            delete_vertex(solved->base, *vx);
            LOG(("Deleted point %d\n", vx->idx));
            LOG(("Number of visible vertices is %d\n", count234(solved->base->vertices)));
        }
    }

    /*
     * Contract the edges in all subgraphs such that every subgraph only consists of
     * a single point.
     */
    for (i = 0; i < n_min; i++)
    {
        bool contr = false;
        edge contre;
        contre.src = i * n_1sub;
        contre.tgt = contre.src + 1;
        while (contre.tgt < (i + 1) * n_1sub)
        {
            if (isedge(solved->base->edges, contre.src, contre.tgt))
            { 
                retoff = sprintf(buf, "%d:%d-%d;", MOVE_CONTREDGE, contre.src, contre.tgt);
                if (retlen + retoff >= retsize)
                {
                    retsize = retlen + retoff + 256;
                    ret = sresize(ret, retsize, char);
                }
                strcpy(ret + retlen, buf);
                retlen += retoff;

                LOG(("Number of visible vertices is %d\n", count234(solved->base->vertices)));
                contract_edge(solved->base, contre.src, contre.tgt);
                LOG(("Contracted edge %d-%d\n", contre.src, contre.tgt));
                LOG(("Number of visible vertices is %d\n", count234(solved->base->vertices)));

                contr = true;
                contre.tgt = contre.src + 1;
            }
            else if (contre.tgt == ((i + 1) * n_1sub) - 1 && !contr)
            {
                contre.src++;
                contre.tgt = contre.src + 1;
            }
            else
            {
                contre.tgt++;
            }
        }
    }

    /*
     * Check whether the found solution is valid. Doing it at this point may not be optimal
     * performance wise but we can ensure that if we are able to find a solution we will
     * always find it.
     */
    solution = NULL;
    if (!(solved->solved = isomorphism_degheuristic(solved->minor, solved->base, &solution)))
    {
#ifdef BENCHMARKS
        assert(!isomorphism_bruteforce(solved->minor, solved->base));
#endif
        sfree(ret);
        ret = NULL;
        free_game(solved);
        solved = NULL;
        retsize = 256;
        retlen = 0;
        clock_t begin = clock();
        if (!(ret = solve_bruteforce(currstate, &solved, &solution, &retsize, &retlen, begin, 100)))
        {
            sfree(solved);
            *error = "Solution not known for the current puzzle state";
            return NULL;
        }
        else
        {
            LOG(("Used bruteforce solver to solve puzzle\n"));
        }
    }

    for (i = 0; (vx = index234(solved->base->vertices, i)) != NULL; i++)
    {
        mapping mvx;
        mvx.key = vx->idx;
        mvx.value = -1;
        mvx.value = ((mapping*) find234(solution, &mvx, NULL))->value;
        LOG(("Vertex %d in the base graph corresponds to vertex %d in the minor\n", vx->idx,
            mvx.value));
        pt = solved->base->points + vx->idx;
        pt->x = solved->minor->points[mvx.value].x;
        pt->y = solved->minor->points[mvx.value].y;
        pt->d = 1;
        retoff = sprintf(buf, "%d:%d-%ld-%ld;", MOVE_DRAGPOINT, vx->idx, pt->x, pt->y);
        if (retlen + retoff >= retsize)
        {
            retsize = retlen + retoff + 256;
            ret = sresize(ret, retsize, char);
        }
        strcpy(ret + retlen, buf);
        retlen += retoff;
    }

    while ((m = delpos234(solution, 0)) != NULL) sfree(m);
    freetree234(solution);
    free_game(solved);
    
    return ret;
}

static bool game_can_format_as_text_now(const game_params *params)
{
    return true;
}

static char *game_text_format(const game_state *state)
{
    return NULL;
}

struct game_ui {

    /* new position of the drag point */
    point newpt;
    /* index of the drag point */
    int dragpt;

    /* index of vertex to delete */
    int delvx;

    /* index of the resulting point of the merge - the dominant point */
    int mergept_dom;
    /* index of the recessive point */
    int mergept_rec;

    /* edge that should be deleted */
    edge deledge;

    /* the currently ongoing move */
    enum move current_move;

};

static game_ui *new_ui(const game_state *state)
{
    game_ui* ret = snew(game_ui);
    ret->newpt.d = 1;
    ret->dragpt = -1;
    ret->delvx = -1;
    ret->mergept_dom = -1;
    ret->mergept_rec = -1;
    ret->deledge.src = -1;
    ret->deledge.tgt = -1;
    ret->current_move = MOVE_IDLE;

    return ret;
}

static void free_ui(game_ui *ui)
{
    sfree(ui);
}

static char *encode_ui(const game_ui *ui)
{
    return NULL;
}

static void decode_ui(game_ui *ui, const char *encoding)
{
}

static void game_changed_state(game_ui *ui, const game_state *oldstate,
                               const game_state *newstate)
{
    switch (ui->current_move)
    {
        case MOVE_CONTREDGE:
            ui->mergept_dom = -1;
            ui->mergept_rec = -1;
        case MOVE_DRAGPOINT:
            ui->dragpt = -1;
            break;
        case MOVE_DELPOINT:
            ui->delvx = -1;
            break;
        case MOVE_DELEDGE:
            ui->deledge.src = -1;
            ui->deledge.tgt = -1;
            break;
        case MOVE_IDLE:
        default:;
    }
    ui->current_move = MOVE_IDLE;
}

struct game_drawstate {

    int tilesize;

    long coord_lim;
    long coord_mar;

    /* game window size */
    long width, height;

};

#define HOVER_POINTRADIUS (POINTRADIUS + 4)

#define POINT_TRESHOLD square(HOVER_POINTRADIUS)
#define OVERLAYPOINT_TRESHOLD_MOVEINT square(2 * POINTRADIUS)
#define EDGE_TRESHOLD 2.0

/* heuristic to determine whether a point has been clicked */
#define point_heuristic(px, py, cx, cy) (square((px) - (cx)) + square((py) - (cy)))

/*
 * Heuristic to determine whether an edge has been clicked
 */
static double edge_heuristic(long esx, long esy, long etx, long ety,
                            long cx, long cy)
{
    double dist_st = sqrt(square(etx - esx) + square(ety - esy));
    double dist_sc = sqrt(square(cx - esx) + square(cy - esy));
    double dist_ct = sqrt(square(etx - cx) + square(ety - cy));
    return dist_sc + dist_ct - dist_st;
}

static char *interpret_move(const game_state *state, game_ui *ui,
                            const game_drawstate *ds,
                            int x, int y, int button)
{
    int i;

    /*
     * Since one can only perform moves on the base graph which is placed in GRID_RIGHT
     * and the point coordinates of a graph are given in relation to its grid one needs
     * to subtract the grid size from the x-coordinate of the mouse event. Make sure the
     * coordinates are correct, even if the player resizes the window.
     */
    long realx = (x * COORDUNIT / ds->tilesize) - state->base->grid.x_off;
    long realy = y * COORDUNIT / ds->tilesize;

    point* pt;
    point* pts = state->base->points;
    vertex* vx;
    tree234* vertices = state->base->vertices;
    tree234* edges = state->base->edges;

    if (IS_MOUSE_DOWN(button))
    {
        int i;
        int bestpt_idx = -1;
        long ptheur;
        long bestpt_heur = POINT_TRESHOLD;
        double eheur;
        double beste_heur = EDGE_TRESHOLD;
        edge beste;
        edge* e;
        beste.src = beste.tgt = -1;

        /*
         * Discard the move if the button that has been pressed is neither LEFT_BUTTON nor
         * RIGHT_BUTTON.
         */
        if (!(button == LEFT_BUTTON || button == RIGHT_BUTTON))
            return NULL;

        /*
         * Get the index of the point with the shortest point heuristic. Do only take points
         * into account with a heuristic smaller than POINT_TRESHOLD.
         */
        for (i = 0; (vx = index234(vertices, i)) != NULL; i++)
        {
            pt = pts + vx->idx;
            ptheur = point_heuristic(pt->x, pt->y, realx, realy);
            if (ptheur < bestpt_heur)
            {
                bestpt_heur = ptheur;
                bestpt_idx = vx->idx;
            }
        }
        /*
         * Check whether there is any point with a point heuristic smaller than the
         * treshold. If yes => game_ui requires update.
         */
        if (bestpt_heur < POINT_TRESHOLD)
        {
            if (button == LEFT_BUTTON)
            {
                ui->current_move = MOVE_DRAGPOINT;
                ui->dragpt = bestpt_idx;
                ui->newpt.x = realx;
                ui->newpt.y = realy;
                LOG(("Updated position of point %d to x:%ld, x:%ld\n", ui->dragpt, ui->newpt.x,
                    ui->newpt.y));
                return UI_UPDATE;
            }
            else
            {
                ui->current_move = MOVE_DELPOINT;
                ui->delvx = bestpt_idx;
                LOG(("Selected point %d to delete\n", ui->delvx));
                return UI_UPDATE;
            }
        }
        else if (button == RIGHT_BUTTON)
        {
            /*
            * Get the index of the edge with the shortest edge heuristic. Do only take edges
            * into account with a heuristic smaller than EDGE_TRESHOLD.
            */
            for (i = 0; (e = index234(edges, i)) != NULL; i++)
            {
                eheur = edge_heuristic(pts[e->src].x, pts[e->src].y, pts[e->tgt].x,
                                        pts[e->tgt].y, realx, realy);
                if (eheur < beste_heur)
                {
                    LOG(("New nearest edge %d-%d with distance %lf to the mouse position\n",
                        e->src, e->tgt, eheur));
                    beste_heur = eheur;
                    beste = *e;
                }
            }
            /*
             * Check whether there is at least one edge with an edge heuristic smaller
             * than the treshold. If yes => game_ui requires update.
             */
            if (beste_heur < EDGE_TRESHOLD && button == RIGHT_BUTTON)
            {
                ui->current_move = MOVE_DELEDGE;
                ui->deledge.src = beste.src;
                ui->deledge.tgt = beste.tgt;
                LOG(("Selected edge %d-%d to delete\n", ui->deledge.src, ui->deledge.tgt));
                return UI_UPDATE;
            }
        }
    }
    else if (IS_MOUSE_DRAG(button))
    {
        if (ui->current_move > MOVE_IDLE && ui->current_move <= MOVE_CONTREDGE)
        {
            ui->newpt.x = realx;
            ui->newpt.y = realy;
            LOG(("Updated position of drag point to x:%ld, y:%ld\n", ui->newpt.x, ui->newpt.y));
            for (i = 0; (vx = index234(vertices, i)) != NULL; i++)
            {
                if (vx->idx == ui->dragpt) continue;
                pt = pts + vx->idx;
                if (square(ui->newpt.x - pt->x) + square(ui->newpt.y - pt->y) < OVERLAYPOINT_TRESHOLD_MOVEINT)
                {
                    if (ui->current_move == MOVE_DRAGPOINT
                        && isedge(edges, ui->dragpt, vx->idx))
                    {
                        ui->current_move = MOVE_CONTREDGE;
                        ui->mergept_rec = ui->dragpt;
                        ui->mergept_dom = vx->idx;
                        LOG(("Converted drag of %d into contraction of %d-%d\n", ui->dragpt,
                            ui->mergept_dom, ui->mergept_rec));
                    }
                    else if (ui->current_move == MOVE_CONTREDGE)
                    {
                        ui->mergept_dom = vx->idx;
                    }
                    return UI_UPDATE;
                }
            }
            if (ui->current_move == MOVE_CONTREDGE)
            {
                ui->current_move = MOVE_DRAGPOINT;
                ui->dragpt = ui->mergept_rec;
                LOG(("Converted contraction of %d-%d to drag of %d to position x:%ld, y:%ld\n",
                    ui->mergept_dom, ui->mergept_rec, ui->dragpt, ui->newpt.x, ui->newpt.y));
                ui->mergept_rec = -1;
                ui->mergept_dom = -1;
            }
            return UI_UPDATE;
        }
    }
    else if (IS_MOUSE_RELEASE(button))
    {
        switch (ui->current_move)
        {
            /*
             * Check for an ongoing drag. If yes, check whether the player wants to discard
             * the drag. If no => make move, else => discard and update game_ui.
             */
            case MOVE_DRAGPOINT:
                if (ui->newpt.x < 2 * POINTRADIUS
                    || ui->newpt.x + (2 * POINTRADIUS) > ds->coord_lim * COORDUNIT / ds->tilesize
                    || ui->newpt.y < 2 * POINTRADIUS
                    || ui->newpt.y + (2 * POINTRADIUS) > ds->coord_lim * COORDUNIT / ds->tilesize)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->dragpt = -1;
                    LOG(("Unselected drag point\n"));
                    return UI_UPDATE;
                }
                else
                {
                    char buf[80];
                    sprintf(buf, "%d:%d-%ld-%ld;", ui->current_move, ui->dragpt, ui->newpt.x,
                            ui->newpt.y);
                    LOG(("Dragging point %d to position x:%ld, y:%ld\n", ui->dragpt, ui->newpt.x,
                        ui->newpt.y));
                    return dupstr(buf);
                }
            /*
             * Check for an ongoing contraction. If yes, check whether the player wants
             * to discard the contraction. If no => make move, else => discard and update
             * game_ui.
             */
            case MOVE_CONTREDGE:
                if (realx < 0 || realx > ds->coord_lim * COORDUNIT / ds->tilesize
                    || y < 0 || y > ds->coord_lim)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->mergept_dom = -1;
                    ui->mergept_rec = -1;
                    LOG(("position x:%ld, y:%ld; x-offset:%d; coordinate limit:%ld\n", realx, realy,
                        state->base->grid.x_off, ds->coord_lim));
                    LOG(("Unselected edge to contract\n"));
                    return UI_UPDATE;
                }
                else
                {
                    char buf[80];
                    sprintf(buf, "%d:%d-%d;", ui->current_move, ui->mergept_dom, ui->mergept_rec);
                    LOG(("Contracting edge between vertices %d and %d\n", ui->mergept_dom,
                        ui->mergept_rec));
                    return dupstr(buf);
                }
            /*
             * Check for an ongoing point deletion. If yes, check whether the player wants
             * to discard the deletion. If no => make move, else => discard and update
             * game_ui.
             */
            case MOVE_DELPOINT:
                if (realx < 0 || realx > ds->coord_lim * COORDUNIT / ds->tilesize
                    || y < 0 || y > ds->coord_lim)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->delvx = -1;
                    LOG(("Unselected point to delete\n"));
                    return UI_UPDATE;
                }
                else
                {
                    char buf[80];
                    sprintf(buf, "%d:%d;", ui->current_move, ui->delvx);
                    LOG(("Deleting point %d\n", ui->delvx));
                    return dupstr(buf);
                }
            /*
             * Check for an ongoing edge deletion. If yes, check whether the player wants
             * to discard the deletion. If no => make move, else => discard and update
             * game_ui.
             */
            case MOVE_DELEDGE:
                if (realx < 0 || realx > ds->coord_lim * COORDUNIT / ds->tilesize
                    || y < 0 || y > ds->coord_lim)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->deledge.src = -1;
                    ui->deledge.tgt = -1;
                    LOG(("Unselected edge to delete\n"));
                    return UI_UPDATE;
                }
                else
                {
                    char buf[80];
                    sprintf(buf, "%d:%d-%d;", ui->current_move, ui->deledge.src, ui->deledge.tgt);
                    LOG(("Deleting edge between vertices %d and %d\n", ui->deledge.src,
                        ui->deledge.tgt));
                    return dupstr(buf);
                }
            case MOVE_IDLE:
            default:;
        }
    }

    return NULL;
}

static game_state *execute_move(const game_state *state, const char *move)
{
    int off;
    int idx;
    int current_move;
    long coord_lim = COORDLIMIT(state->params.n_base) * COORDUNIT;
    point p;
    vertex v;
    edge e;
    game_state* ret;
    v.deg = 0;

    ret = dup_game(state);
    if (*move == 'S')
    {
        if (!ret->solved) ret->cheated = true;
        ret->solved = true;
        move++;
    }
    if (!(*move && isdigit((uint8) *move)))
    {
        LOG(("Failed to scan move, empty move string\n"));
        free_game(ret);
        return NULL;
    }

    /*
     * Parse the move description. Return NULL if either the description is incorrect
     * or the move is invalid.
     */
    do
    {
        if (sscanf(move, "%d:%n", &current_move, &off) != 1)
        {
            LOG(("Failed to scan new current move\n"));
            free_game(ret);
            return NULL;
        }

        move += off;

        switch (current_move)
        {
            case MOVE_DRAGPOINT:
                idx = -1;
                p.x = p.y = -1;
                if (sscanf(move, "%d-%ld-%ld;%n", &idx, &p.x, &p.y, &off) != 3
                    || idx < 0 || idx >= state->params.n_base
                    || p.x < 0 || p.x > coord_lim || p.y < 0 || p.y > coord_lim)
                {
                    LOG(("Failed to scan drag move,"));
                    if (idx >= 0) LOG((" point %d\n", idx));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }

                /* Assign the new coordinates to the dragged point */
                ret->base->points[idx].x = p.x;
                ret->base->points[idx].y = p.y;
                LOG(("Dragged point %d to position x:%ld, y:%ld\n", idx, p.x, p.y));

                move += off;
                break;
            case MOVE_CONTREDGE:
                e.src = e.tgt = -1;
                if (sscanf(move, "%d-%d;%n", &e.src, &e.tgt, &off) != 2
                    || e.src == e.tgt || !isedge(ret->base->edges, e.src, e.tgt))
                {
                    LOG(("Failed to scan contraction move,"));
                    if (e.src >= 0 && e.tgt >= 0) LOG((" edge %d-%d\n", e.src, e.tgt));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }
                
                /* Contract the edge between src and tgt */
                contract_edge(ret->base, e.src, e.tgt);
                LOG(("Contracted edge %d-%d\n", e.src, e.tgt));
                if (ret->solved) ret->solved = false;

                move += off;
                break;
            case MOVE_DELPOINT:
                if (sscanf(move, "%d;%n", &v.idx, &off) != 1
                    || (v.idx < 0 && v.idx >= state->params.n_base))
                {
                    LOG(("Failed to scan point deletion move,"));
                    if (v.idx >= 0) LOG((" point %d\n", v.idx));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }

                delete_vertex(ret->base, v);
                LOG(("Deleted point %d\n", v.idx));
                if (ret->solved) ret->solved = false;

                move += off;
                break;
            case MOVE_DELEDGE:
                e.src = e.tgt = -1;
                if (sscanf(move, "%d-%d;%n", &e.src, &e.tgt, &off) != 2
                    || e.src == e.tgt || !isedge(ret->base->edges, e.src, e.tgt))
                {
                    LOG(("Failed to scan edge deletion move,"));
                    if (e.src >= 0 && e.tgt >= 0) LOG((" edge %d-%d\n", e.src, e.tgt));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }

                delete_edge(ret->base, e);
                LOG(("Deleted edge %d-%d\n", e.src, e.tgt));
                if (ret->solved) ret->solved = false;

                move += off;
                break;
            case MOVE_IDLE:
                sscanf(move, ";%n", &off);
                if (off != 1)
                {
                    LOG(("Failed to scan idle move\n"));
                    free_game(ret);
                    return NULL;
                }

                LOG(("Idled\n"));

                move += off;
                break;
            default:
                free_game(ret);
                return NULL;
        }

        if (!(current_move == MOVE_DRAGPOINT || state->solved || ret->solved))
        {
            ret->solved = isomorphism_degheuristic(ret->minor, ret->base, NULL);
#ifdef BENCHMARKS
            assert(ret->solved == isomorphism_bruteforce(ret->minor, ret->base));
#endif
        }
    }
    while (*move && isdigit((uint8) *move));

    return ret;
}

/* ----------------------------------------------------------------------
 * Drawing routines.
 */

static void game_compute_size(const game_params *params, int tilesize,
                              int *x, int *y)
{
    int coord_lim = COORDLIMIT(params->n_base) * tilesize;
    *x = coord_lim + (((float) COORDLIMIT(params->n_min)
        / (float) COORDLIMIT(params->n_base)) * (float) coord_lim);
    *y = coord_lim;
}

static void game_set_size(drawing *dr, game_drawstate *ds,
                          const game_params *params, int tilesize)
{
    ds->tilesize = tilesize;
    ds->coord_lim = COORDLIMIT(params->n_base) * tilesize;
    ds->coord_mar = COORDMARGIN * tilesize;
    ds->width = ds->coord_lim + (((float) COORDLIMIT(params->n_min)
                / (float) COORDLIMIT(params->n_base)) * (float) ds->coord_lim);
    ds->height = ds->coord_lim;
}

static float *game_colours(frontend *fe, int *ncolours)
{
    float *ret = snewn(3 * NCOLOURS, float);

    game_mkhighlight(fe, ret, COL_BACKGROUND, -1, COL_SYSBACKGROUND);

    /* dark grey */
    ret[(COL_OUTLINE * 3) + 0] = 0.3F;
    ret[(COL_OUTLINE * 3) + 1] = 0.3F;
    ret[(COL_OUTLINE * 3) + 2] = 0.3F;

    /* dark grey */
    ret[(COL_GRIDBORDER * 3) + 0] = 0.3F;
    ret[(COL_GRIDBORDER * 3) + 1] = 0.3F;
    ret[(COL_GRIDBORDER * 3) + 2] = 0.3F;

    /* blue */
    ret[(COL_BASEPOINT * 3) + 0] = 0.0F;
    ret[(COL_BASEPOINT * 3) + 1] = 0.0F;
    ret[(COL_BASEPOINT * 3) + 2] = 1.0F;

    /* red */
    ret[(COL_MINORPOINT * 3) + 0] = 1.0F;
    ret[(COL_MINORPOINT * 3) + 1] = 0.0F;
    ret[(COL_MINORPOINT * 3) + 2] = 0.0F;

    /* green */
    ret[(COL_DRAGPOINT * 3) + 0] = 0.0F;
    ret[(COL_DRAGPOINT * 3) + 1] = 1.0F;
    ret[(COL_DRAGPOINT * 3) + 2] = 0.0F;

    /* orange */
    ret[(COL_DELPOINT * 3) + 0] = 1.0F;
    ret[(COL_DELPOINT * 3) + 1] = 0.5F;
    ret[(COL_DELPOINT * 3) + 2] = 0.0F;

    /* light blue */
    ret[(COL_HIDEPOINT * 3) + 0] = 0.3F;
    ret[(COL_HIDEPOINT * 3) + 1] = 0.7F;
    ret[(COL_HIDEPOINT * 3) + 2] = 1.0F;

    /* black */
    ret[(COL_POINTOUTLINE * 3) + 0] = 0.0F;
    ret[(COL_POINTOUTLINE * 3) + 1] = 0.0F;
    ret[(COL_POINTOUTLINE * 3) + 2] = 0.0F;

    /* black */
    ret[(COL_EDGE * 3) + 0] = 0.0F;
    ret[(COL_EDGE * 3) + 1] = 0.0F;
    ret[(COL_EDGE * 3) + 2] = 0.0F;

    /* orange */
    ret[(COL_DELEDGE * 3) + 0] = 1.0F;
    ret[(COL_DELEDGE * 3) + 1] = 0.5F;
    ret[(COL_DELEDGE * 3) + 2] = 0.0F;

    /* grey */
    ret[(COL_HIDEEDGE * 3) + 0] = 0.4F;
    ret[(COL_HIDEEDGE * 3) + 1] = 0.4F;
    ret[(COL_HIDEEDGE * 3) + 2] = 0.4F;

    /* white */
    ret[(COL_FLASH * 3) + 0] = 1.0F;
    ret[(COL_FLASH * 3) + 1] = 1.0F;
    ret[(COL_FLASH * 3) + 2] = 1.0F;

    /* grey */
    ret[(COL_FLASH2 * 3) + 0] = 0.5F;
    ret[(COL_FLASH2 * 3) + 1] = 0.5F;
    ret[(COL_FLASH2 * 3) + 2] = 0.5F;

#if DEBUG
    /* cyan */
    ret[(COL_SUBPOINT * 3) + 0] = 0.0F;
    ret[(COL_SUBPOINT * 3) + 1] = 1.0F;
    ret[(COL_SUBPOINT * 3) + 2] = 1.0F;
#endif

    *ncolours = NCOLOURS;
    return ret;
}

static game_drawstate *game_new_drawstate(drawing *dr, const game_state *state)
{
    struct game_drawstate *ds = snew(struct game_drawstate);

    ds->tilesize = COORDUNIT;
    ds->coord_lim = COORDLIMIT(state->params.n_base) * COORDUNIT;
    ds->coord_mar = COORDMARGIN * COORDUNIT;
    ds->width = ds->coord_lim + (((float) COORDLIMIT(state->params.n_min)
                / (float) COORDLIMIT(state->params.n_base)) * (float) ds->coord_lim);
    ds->height = ds->coord_lim;

    return ds;
}

static void game_free_drawstate(drawing *dr, game_drawstate *ds)
{
    sfree(ds);
}

#define ANIM_LENGTH 0.7F
#define FLASH_LENGTH 0.3F

static void game_redraw(drawing *dr, game_drawstate *ds,
                        const game_state *oldstate, const game_state *state,
                        int dir, const game_ui *ui,
                        float animtime, float flashtime)
{
    int i;
    int bg_color;
    float r_anim;

    int mxoff = state->minor->grid.x_off;
    int myoff = state->minor->grid.y_off;
    float mrelsize = state->minor->grid.relsize;

    int bxoff = state->base->grid.x_off;
    int byoff = state->base->grid.y_off;
    float brelsize = state->base->grid.relsize;

    edge* e;
    vertex* vx;
    point* esrc;
    point* etgt;
    point* pts;

    /*
     * Check whether game has been recently solved and solve function
     * hasn't been used.
     */
    if (!flashtime)
        bg_color = COL_BACKGROUND;
    else if (flashtime < (FLASH_LENGTH / 3.0F) ||
            flashtime > 2.0F * (FLASH_LENGTH / 3.0F))
        bg_color = COL_FLASH;
    else
        bg_color = COL_FLASH2;

    if (oldstate) r_anim = animtime / ANIM_LENGTH;
    else r_anim = 1.0F;
    
    /*
     * The initial contents of the window are not guaranteed and
     * can vary with front ends. To be on the safe side, all games
     * should start by drawing a big background-colour rectangle
     * covering the whole window.
     */
    draw_rect(dr, 0, 0, ds->width, ds->height, bg_color);
    draw_rect_outline(dr, 0, 0, ds->width, ds->height, COL_OUTLINE);

    /* Draw the minor grid outline */
    draw_rect_outline(dr, mxoff, myoff, ds->coord_lim * mrelsize,
                    ds->coord_lim * mrelsize, COL_OUTLINE);
    
    /* Draw the minor edges in the intended grid */
    pts = state->minor->points;
    for (i = 0; (e = index234(state->minor->edges, i)) != NULL; i++)
    {
        point psrc, ptgt;
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        psrc.x = (esrc->x + mxoff) * mrelsize * ds->tilesize / COORDUNIT;
        psrc.y = (esrc->y + myoff) * mrelsize * ds->tilesize / COORDUNIT;
        ptgt.x = (etgt->x + mxoff) * mrelsize * ds->tilesize / COORDUNIT;
        ptgt.y = (etgt->y + myoff) * mrelsize * ds->tilesize / COORDUNIT;
        draw_line(dr, psrc.x, psrc.y, ptgt.x, ptgt.y, COL_EDGE);
    }
    /* Draw the minor points in the intended grid */
    for (i = 0; (vx = index234(state->minor->vertices, i)) != NULL; i++)
    {
        long r;
        point p;
        p.x = (pts[vx->idx].x + mxoff) * mrelsize * ds->tilesize / COORDUNIT;
        p.y = (pts[vx->idx].y + myoff) * mrelsize * ds->tilesize / COORDUNIT;
        r = POINTRADIUS * ds->tilesize / COORDUNIT;
        draw_circle(dr, p.x, p.y, r, COL_MINORPOINT, COL_POINTOUTLINE);
    }

    /* Draw the base graph edges in the intended grid */
    pts = state->base->points;
    for (i = 0; (e = index234((dir > 0 && oldstate) ? oldstate->base->edges :
        state->base->edges, i)) != NULL; i++)
    {
        /*
         * Check whether there is a solve animation ongoing after which
         * the current edge will disappear.
         */
        bool hide = find234((dir < 0 && oldstate) ? oldstate->base->edges :
                            state->base->edges, e, NULL) == NULL;
        point psrc, ptgt;
        point* oesrc;
        point* oetgt;
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        /*
         * Check whether the edge source vertex is being dragged. Use the
         * vertex's new coordinates if so.
         */
        if (ui->dragpt == e->src)
        {
            psrc.x = (ui->newpt.x + bxoff) * brelsize * ds->tilesize / COORDUNIT;
            psrc.y = (ui->newpt.y + byoff) * brelsize * ds->tilesize / COORDUNIT;
        }
        /*
         * Check whether there is a solve animation ongoing. Calculate the
         * vertex's current animated position if so.
         */
        else if (oldstate)
        {
            oesrc = oldstate->base->points + e->src;
            psrc.x = (oesrc->x + ((float) (esrc->x - oesrc->x) * r_anim) + bxoff)
                    * brelsize * ds->tilesize / COORDUNIT;
            psrc.y = (oesrc->y + ((float) (esrc->y - oesrc->y) * r_anim) + byoff)
                    * brelsize * ds->tilesize / COORDUNIT;
        }
        else
        {
            psrc.x = (esrc->x + bxoff) * brelsize * ds->tilesize / COORDUNIT;
            psrc.y = (esrc->y + byoff) * brelsize * ds->tilesize / COORDUNIT;
        }
        if (ui->dragpt == e->tgt)
        {
            ptgt.x = (ui->newpt.x + bxoff) * brelsize * ds->tilesize / COORDUNIT;
            ptgt.y = (ui->newpt.y + byoff) * brelsize * ds->tilesize / COORDUNIT;
        }
        else if (oldstate)
        {
            oetgt = oldstate->base->points + e->tgt;
            ptgt.x = (oetgt->x + ((float) (etgt->x - oetgt->x) * r_anim) + bxoff)
                    * brelsize * ds->tilesize / COORDUNIT;
            ptgt.y = (oetgt->y + ((float) (etgt->y - oetgt->y) * r_anim) + byoff)
                    * brelsize * ds->tilesize / COORDUNIT;
        }
        else
        {
            ptgt.x = (etgt->x + bxoff) * brelsize * ds->tilesize / COORDUNIT;
            ptgt.y = (etgt->y + byoff) * brelsize * ds->tilesize / COORDUNIT;
        }
        draw_line(dr, psrc.x, psrc.y, ptgt.x, ptgt.y,
                    ((e->src == ui->deledge.src && e->tgt == ui->deledge.tgt)
                    || e->src == ui->delvx || e->tgt == ui->delvx) ?
                    COL_DELEDGE : (hide ? COL_HIDEEDGE : COL_EDGE));
    }
    /* Draw the base graph points in the intended grid */
    for (i = 0; (vx = index234((dir > 0 && oldstate) ? oldstate->base->vertices :
        state->base->vertices, i)) != NULL; i++)
    {
        bool hide = find234((dir < 0 && oldstate) ? oldstate->base->vertices :
                            state->base->vertices, ((dir < 0 && oldstate) ?
                            oldstate->base->vtcs : state->base->vtcs) + vx->idx, NULL)
                    == NULL;
        long r;
        point p;
        point* op;
        if (vx->idx == ui->dragpt)
        {
            continue;
        }
        else if (oldstate)
        {
            op = oldstate->base->points + vx->idx;
            p.x = (op->x + ((float) (pts[vx->idx].x - op->x) * r_anim) + bxoff)
                    * brelsize * ds->tilesize / COORDUNIT;
            p.y = (op->y + ((float) (pts[vx->idx].y - op->y) * r_anim) + byoff)
                    * brelsize * ds->tilesize / COORDUNIT;
        }
        else
        {
            p.x = (pts[vx->idx].x + bxoff) * brelsize * ds->tilesize / COORDUNIT;
            p.y = (pts[vx->idx].y + byoff) * brelsize * ds->tilesize / COORDUNIT;
        }
        if (vx->idx == ui->mergept_dom)
            r = HOVER_POINTRADIUS * ds->tilesize / COORDUNIT;
        else
            r = POINTRADIUS * ds->tilesize / COORDUNIT;
        draw_circle(dr, p.x, p.y, r, (vx->idx == ui->delvx) ? COL_DELPOINT :
                    ( hide ? COL_HIDEPOINT :
#if DEBUG
                    ((vx->idx < state->params.n_min * (state->params.n_base
                    / state->params.n_min)) ? COL_SUBPOINT : COL_BASEPOINT)),
#else
                    COL_BASEPOINT),
#endif
                    COL_POINTOUTLINE);
    }
    /*
     * If a drag is ongoing draw the drag point. The drag point is drawn last
     * such that when we hover over another point it doesn't hide behind that
     * point.
     */
    if (ui->dragpt > -1)
    {
        long r;
        point p;
        p.x = (ui->newpt.x + bxoff) * brelsize * ds->tilesize / COORDUNIT;
        p.y = (ui->newpt.y + byoff) * brelsize * ds->tilesize / COORDUNIT;
        r = POINTRADIUS * ds->tilesize / COORDUNIT;
        draw_circle(dr, p.x, p.y, r, COL_DRAGPOINT, COL_POINTOUTLINE);
    }

    draw_update(dr, 0, 0, ds->width, ds->height);
}

static float game_anim_length(const game_state *oldstate,
                              const game_state *newstate, int dir, game_ui *ui)
{
    if (((dir < 0) ? oldstate : newstate)->solved
        && !((dir > 0) ? oldstate : newstate)->solved
        && ((dir < 0) ? oldstate : newstate)->cheated)
        return ANIM_LENGTH;
    else
        return 0.0F;
}

static float game_flash_length(const game_state *oldstate,
                               const game_state *newstate, int dir, game_ui *ui)
{
    /*
     * Check whether game has been recently solved and solve function
     * hasn't been used.
     */
    if (newstate->solved && !oldstate->solved && !newstate->cheated)
        return FLASH_LENGTH;
    else
        return 0.0F;
}

static int game_status(const game_state *state)
{
    return state->solved ? 1 : 0;
}

static bool game_timing_state(const game_state *state, game_ui *ui)
{
    return true;
}

static void game_print_size(const game_params *params, float *x, float *y)
{
}

static void game_print(drawing *dr, const game_state *state, int tilesize)
{
}

#ifdef COMBINED
#define thegame minorfinder
#endif

const struct game thegame = {
    "Minor Finder", "games.minorfinder", "minorfinder",
    default_params,                                                     /* done */
    NULL, preset_menu,                                                  /* done */
    decode_params,                                                      /* done */
    encode_params,                                                      /* done */
    free_params,                                                        /* done */
    dup_params,                                                         /* done */
    false, game_configure, custom_params,
    validate_params,                                                    /* done */
    new_game_desc,                                                      /* done */
    validate_desc,                                                      /* done */
    new_game,                                                           /* done */
    dup_game,                                                           /* done */
    free_game,                                                          /* done */
    true, solve_game,                                                   /* done */
    false, game_can_format_as_text_now, game_text_format,
    new_ui,                                                             /* done */
    free_ui,                                                            /* done */
    encode_ui,
    decode_ui,
    NULL, /* game_request_keys */
    game_changed_state,                                                 /* done */
    interpret_move,                                                     /* done */
    execute_move,                                                       /* done */
    COORDUNIT, game_compute_size, game_set_size,                        /* done */
    game_colours,                                                       /* done */
    game_new_drawstate,                                                 /* done */
    game_free_drawstate,                                                /* done */
    game_redraw,                                                        /* done */
    game_anim_length,
    game_flash_length,                                                  /* done */
    game_status,                                                        /* done */
    false, false, game_print_size, game_print,
    false,			       /* wants_statusbar */
    false, game_timing_state,
    SOLVE_ANIMATES,				       /* flags */
};
