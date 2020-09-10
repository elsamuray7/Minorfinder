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
    COL_CONTREDGE,
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

enum default_params {
    DEFAULT_N_BASE = 15,
    DEFAULT_N_MIN = 4
};

/*
 * A grid, defines a region of a window
 */
enum grid {
    GRID_LEFT,
    GRID_RIGHT,
    NGRIDS
};

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
    enum grid grid;

    /* array of points - remains the same throughout the game */
    point* points;
    /* array of vertices - remains the same throughout the game */
    vertex* vtcs;
    /* 234-tree of vertices - maps the current game state */
    tree234* vertices;
    /* array of 234-trees of eaten vertices during contractions */
    tree234** eaten_vertices;

    /* 234-tree of edges - maps the current game state */
    tree234* edges;

} graph;

struct game_params {

    /* number of base graph points */
    int n_base;

    /* number of minor graph points */
    int n_min;

};

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

    ret->n_base = DEFAULT_N_BASE;
    ret->n_min = DEFAULT_N_MIN;

    return ret;
}

static bool game_fetch_preset(int i, char **name, game_params **params)
{
    game_params* ret;
    int n_base;
    int n_min;
    char buf[80];

    switch (i)
    {
        case 0:
            n_base = DEFAULT_N_BASE;
            n_min = DEFAULT_N_MIN;
            break;
        case 1:
            n_base = 19;
            n_min = 5;
            break;
        case 2:
            n_base = 23;
            n_min = 6;
            break;
        default:
            return false;
    }

    sprintf(buf, "%d base, %d minor points", n_base, n_min);
    *name = dupstr(buf);

    ret = snew(game_params);
    ret->n_base = n_base;
    ret->n_min = n_min;
    *params = ret;
    
    return true;
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
    if (*string == 'b')
    {
        string++;
        if (*string && isdigit((uint8) *string))
        {
            params->n_base = atoi(string);
            while (*string && isdigit((uint8) *string)) string++;
            if (*string == '-')
            {
                string++;
                if (*string == 'm')
                {
                    string++;
                    if (*string && isdigit((uint8) *string))
                    {
                        params->n_min = atoi(string);
                        return; /* params were correctly encoded */
                    }
                }
            }
        }
    }

    /* params encoding was incorrect */
    params->n_base = DEFAULT_N_BASE;
    params->n_min = DEFAULT_N_MIN;
}

static char *encode_params(const game_params *params, bool full)
{
    char buf[80];

    sprintf(buf, "b%d-m%d", params->n_base, params->n_min);
    
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
    if (params->n_base < DEFAULT_N_BASE || params->n_base > 30)
        return "Number of base graph points is invalid";
    if (params->n_min < DEFAULT_N_MIN || params->n_min > 15)
        return "Number of minor points is invalid";
    
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

/* ----------------------------------------------------------------------
 * Small number of 64-bit integer arithmetic operations, to prevent
 * integer overflow at the very core of cross().
 */

typedef struct {
    long hi;
    ulong lo;
} int64;

#define greater64(i,j) ( (i).hi>(j).hi || ((i).hi==(j).hi && (i).lo>(j).lo))
#define sign64(i) ((i).hi < 0 ? -1 : (i).hi==0 && (i).lo==0 ? 0 : +1)

static int64 mulu32to64(ulong x, ulong y)
{
    ulong a, b, c, d, t;
    int64 ret;

    a = (x & 0xFFFF) * (y & 0xFFFF);
    b = (x & 0xFFFF) * (y >> 16);
    c = (x >> 16) * (y & 0xFFFF);
    d = (x >> 16) * (y >> 16);

    ret.lo = a;
    ret.hi = d + (b >> 16) + (c >> 16);
    t = (b & 0xFFFF) << 16;
    ret.lo += t;
    if (ret.lo < t)
	ret.hi++;
    t = (c & 0xFFFF) << 16;
    ret.lo += t;
    if (ret.lo < t)
	ret.hi++;

#ifdef DIAGNOSTIC_VIA_LONGLONG
    assert(((ulong long)ret.hi << 32) + ret.lo ==
	   (ulong long)x * y);
#endif

    return ret;
}

static int64 mul32to64(long x, long y)
{
    int sign = +1;
    int64 ret;
#ifdef DIAGNOSTIC_VIA_LONGLONG
    long long realret = (long long)x * y;
#endif

    if (x < 0)
	x = -x, sign = -sign;
    if (y < 0)
	y = -y, sign = -sign;

    ret = mulu32to64(x, y);

    if (sign < 0) {
	ret.hi = -ret.hi;
	ret.lo = -ret.lo;
	if (ret.lo)
	    ret.hi--;
    }

#ifdef DIAGNOSTIC_VIA_LONGLONG
    assert(((ulong long)ret.hi << 32) + ret.lo == realret);
#endif

    return ret;
}

static int64 dotprod64(long a, long b, long p, long q)
{
    int64 ab, pq;

    ab = mul32to64(a, b);
    pq = mul32to64(p, q);
    ab.hi += pq.hi;
    ab.lo += pq.lo;
    if (ab.lo < pq.lo)
	ab.hi++;
    return ab;
}

/*
 * Determine whether the line segments between a1 and a2, and
 * between b1 and b2, intersect. We count it as an intersection if
 * any of the endpoints lies _on_ the other line.
 */
static bool cross(point a1, point a2, point b1, point b2)
{
    long b1x, b1y, b2x, b2y, px, py;
    int64 d1, d2, d3;

    /*
     * The condition for crossing is that b1 and b2 are on opposite
     * sides of the line a1-a2, and vice versa. We determine this
     * by taking the dot product of b1-a1 with a vector
     * perpendicular to a2-a1, and similarly with b2-a1, and seeing
     * if they have different signs.
     */

    /*
     * Construct the vector b1-a1. We don't have to worry too much
     * about the denominator, because we're only going to check the
     * sign of this vector; we just need to get the numerator
     * right.
     */
    b1x = b1.x * a1.d - a1.x * b1.d;
    b1y = b1.y * a1.d - a1.y * b1.d;
    /* Now construct b2-a1, and a vector perpendicular to a2-a1,
     * in the same way. */
    b2x = b2.x * a1.d - a1.x * b2.d;
    b2y = b2.y * a1.d - a1.y * b2.d;
    px = a1.y * a2.d - a2.y * a1.d;
    py = a2.x * a1.d - a1.x * a2.d;
    /* Take the dot products. Here we resort to 64-bit arithmetic. */
    d1 = dotprod64(b1x, px, b1y, py);
    d2 = dotprod64(b2x, px, b2y, py);
    /* If they have the same non-zero sign, the lines do not cross. */
    if ((sign64(d1) > 0 && sign64(d2) > 0) ||
	(sign64(d1) < 0 && sign64(d2) < 0))
	return false;

    /*
     * If the dot products are both exactly zero, then the two line
     * segments are collinear. At this point the intersection
     * condition becomes whether or not they overlap within their
     * line.
     */
    if (sign64(d1) == 0 && sign64(d2) == 0) {
	/* Construct the vector a2-a1. */
	px = a2.x * a1.d - a1.x * a2.d;
	py = a2.y * a1.d - a1.y * a2.d;
	/* Determine the dot products of b1-a1 and b2-a1 with this. */
	d1 = dotprod64(b1x, px, b1y, py);
	d2 = dotprod64(b2x, px, b2y, py);
	/* If they're both strictly negative, the lines do not cross. */
	if (sign64(d1) < 0 && sign64(d2) < 0)
	    return false;
	/* Otherwise, take the dot product of a2-a1 with itself. If
	 * the other two dot products both exceed this, the lines do
	 * not cross. */
	d3 = dotprod64(px, px, py, py);
	if (greater64(d1, d3) && greater64(d2, d3))
	    return false;
    }

    /*
     * We've eliminated the only important special case, and we
     * have determined that b1 and b2 are on opposite sides of the
     * line a1-a2. Now do the same thing the other way round and
     * we're done.
     */
    b1x = a1.x * b1.d - b1.x * a1.d;
    b1y = a1.y * b1.d - b1.y * a1.d;
    b2x = a2.x * b1.d - b1.x * a2.d;
    b2y = a2.y * b1.d - b1.y * a2.d;
    px = b1.y * b2.d - b2.y * b1.d;
    py = b2.x * b1.d - b1.x * b2.d;
    d1 = dotprod64(b1x, px, b1y, py);
    d2 = dotprod64(b2x, px, b2y, py);
    if ((sign64(d1) > 0 && sign64(d2) > 0) ||
	(sign64(d1) < 0 && sign64(d2) < 0))
	return false;

    /*
     * The lines must cross.
     */
    return true;
}

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

#define POINTRADIUS 6
#define CROSSPOINT_THRESHOLD (POINTRADIUS / 2)

#define square(x) ((x) * (x))

/*
 * Check whether the edge between s and t crosses the point p
 */
static bool crosspoint(point s, point t, point p)
{
    long dist_st = squarert(square(t.x - s.x) + square(t.y - s.y));
    long dist_sp = squarert(square(p.x - s.x) + square(p.y - s.y));
    long dist_pt = squarert(square(t.x - p.x) + square(t.y - p.y));

    return dist_sp + dist_pt - dist_st < CROSSPOINT_THRESHOLD;
}

/*
 * Add edges between the vertices in the source range and vertices in the target range.
 * Make sure that edges don't cross other edges or points and the degree of the involved
 * vertices doesn't increase beyond __max_deg.
 */
static void addedges(tree234* edges, vertex* vertices, point* points, const int __off_src_vtcs,
                    const int __n_src_vtcs, const int __off_tgt_vtcs, const int __n_tgt_vtcs,
                    const int n_pts, const int __max_deg, random_state* rs)
{
    bool* contains;
    int i;
#if DEBUG
    int j;
#endif
    int off_src_vtcs;
    int n_src_vtcs;
    int off_tgt_vtcs;
    int n_tgt_vtcs;
    int max_deg;
    vertex* vxa;
    vertex* vxb;
    tree234* vtcs;
    tree234* src_vtcs;
    tree234** tgt_vtcs;
    edge* e;

    if (__off_src_vtcs <= __off_tgt_vtcs)
    {
        off_src_vtcs = __off_src_vtcs;
        n_src_vtcs = __n_src_vtcs;
        off_tgt_vtcs = __off_tgt_vtcs;
        n_tgt_vtcs = __n_tgt_vtcs;
    }
    else
    {
        off_src_vtcs = __off_tgt_vtcs;
        n_src_vtcs = __n_tgt_vtcs;
        off_tgt_vtcs = __off_src_vtcs;
        n_tgt_vtcs = __n_src_vtcs;
    }
    if (n_src_vtcs < 0)
        n_src_vtcs = off_tgt_vtcs + n_tgt_vtcs - off_src_vtcs;
    else if (n_tgt_vtcs < 0)
        n_tgt_vtcs = off_src_vtcs + n_src_vtcs - off_tgt_vtcs;
    if (__max_deg < 0) max_deg = max(n_src_vtcs, n_tgt_vtcs) - 1;
    else max_deg = __max_deg;

    /* add all source vertices to a 234-tree */
    src_vtcs = newtree234(vertcmp);
    LOG(("Initially added vertices "));
    for (i = off_src_vtcs; i < off_src_vtcs + n_src_vtcs; i++)
    {
        add234(src_vtcs, vertices + i);
        LOG(("%d, ", vertices[i].idx));
    }
    LOG(("to source vertices\n"));

    /* for every source vertex add all target vertices to a 234-tree */
    tgt_vtcs = snewn(n_src_vtcs, tree234*);
    *tgt_vtcs = newtree234(vertcmp);
    LOG(("Initially added vertices "));
    for (i = off_tgt_vtcs; i < off_tgt_vtcs + n_tgt_vtcs; i++)
    {
        add234(*tgt_vtcs, vertices + i);
        LOG(("%d, ", vertices[i].idx));
    }
    LOG(("to target vertices\n"));
    for (i = 1; i < n_src_vtcs; i++)
    {
        tgt_vtcs[i] = copytree234(*tgt_vtcs, NULL, NULL);
#if DEBUG
        assert(count234(tgt_vtcs[i]) == count234(*tgt_vtcs));
        for (j = 0; j < count234(tgt_vtcs[i]); j++)
        {
            vxa = index234(tgt_vtcs[i], j);
            vxb = index234(*tgt_vtcs, j);
            assert(vxa->idx == vxb->idx);
        }
#endif
    }
    
    /*
     * Start adding edges. Randomly pick two vertices and try to add an edge between them.
     * If the edge can't be added delete the target vertex from the corresponding target
     * 234-tree. Otherwise add the edge and update the involved vertices in the source
     * 234-tree and all target 234-trees. Repeat until there are no more edges to add.
     */
    while (count234(src_vtcs))
    {
        vxa = index234(src_vtcs, random_upto(rs, count234(src_vtcs)));
        vtcs = tgt_vtcs[vxa->idx - off_src_vtcs];
        if (!count234(vtcs))
        {
            del234(src_vtcs, vxa);
            LOG(("Removed vertex %d from source vertices, no more edges to add\n",
                vxa->idx));
            continue;
        }
        vxb = index234(vtcs, random_upto(rs, count234(vtcs)));
        LOG(("Trying to add edge between vertices %d and %d\n", vxa->idx, vxb->idx));

        if (vxb->idx <= vxa->idx)
        {
            del234(vtcs, vxb);
            LOG(("Removed vertex %d from target vertices of vertex %d,"\
                " the edge can't be added\n", vxb->idx, vxa->idx));
            continue;
        }
        /* check for crossing edges */
        for (i = 0; (e = index234(edges, i)) != NULL; i++)
        {
            if (vxa->idx == e->src || vxa->idx == e->tgt || vxb->idx == e->src ||
                vxb->idx == e->tgt)
            {
                continue;
            }
            else if (cross(points[vxa->idx], points[vxb->idx], points[e->src],
                            points[e->tgt]))
            {
                del234(vtcs, vxb);
                LOG(("Removed vertex %d from target vertices of vertex %d,"\
                    " the edge crosses the edge %d-%d\n", vxb->idx, vxa->idx,
                    e->src, e->tgt));
                goto next_vertices; /* this edge crosses another edge => next vertex pair */
            }
        }
        /* check for crossing points */
        for (i = 0; i < n_pts; i++)
        {
            if (i == vxa->idx || i == vxb->idx)
            {
                continue;
            }
            else if (crosspoint(points[vxa->idx], points[vxb->idx], points[i]))
            {
                del234(vtcs, vxb);
                LOG(("Removed vertex %d from target vertices of vertex %d,"\
                    " the edge crosses the point %d\n", vxb->idx, vxa->idx, i));
                goto next_vertices; /* this edge crosses a point => next vertex pair */
            }
        }

        addedge(edges, vxa->idx, vxb->idx);
        LOG(("Added edge between vertices %d and %d\n", vxa->idx, vxb->idx));
        contains = snewn(n_src_vtcs, bool);
        del234(src_vtcs, vxa);
        for (i = 0; i < n_src_vtcs; i++)
        {
            contains[i] = del234(tgt_vtcs[i], vxa);
        }
        vxa->deg++;
        if (vxa->deg < max_deg)
        {
            add234(src_vtcs, vxa);
            for (i = 0; i < n_src_vtcs; i++)
            {
                if (contains[i]) add234(tgt_vtcs[i], vxa);
            }
        }
        contains = sresize(contains, n_src_vtcs + 1, bool);
        *contains = del234(src_vtcs, vxb);
        for (i = 0; i < n_src_vtcs; i++)
        {
            contains[i+1] = del234(tgt_vtcs[i], vxb);
        }
        vxb->deg++;
        if (vxb->deg < max_deg)
        {
            if (*contains) add234(src_vtcs, vxb);
            for (i = 0; i < n_src_vtcs; i++)
            {
                if (contains[i+1] && i != vxa->idx - off_src_vtcs)
                    add234(tgt_vtcs[i], vxb);
            }
        }
        sfree(contains);
        LOG(("Updated degrees of vertices %d to %d and %d to %d\n", vxa->idx, vxa->deg,
            vxb->idx, vxb->deg));

        next_vertices:;
    }

    freetree234(src_vtcs);
    for (i = 0; i < n_src_vtcs; i++) freetree234(tgt_vtcs[i]);
    sfree(tgt_vtcs);
}

/*
 * These parameters are highly sensitive, changing them may cause problems when
 * generating new game descriptions.
 */
#define COORDMARGIN 0.5
#define COORDLIMIT(n) ((n) + (2 * COORDMARGIN))
#define COORDUNIT 24

#define MINORRADIUS(n) ((4 * (n)) / 11)
#define SUBGRAPH_DISTANCE (4 * POINTRADIUS)
#define SUBGRAPH_POINTENTROPY 3
#define OVERLAYPOINT_TRESHOLD square(4 * POINTRADIUS)

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

    /* Arrange the minor points in a circle with radius circle_rad */
    tmp = coord_lim - (2 * COORDMARGIN);
    pts_min = snewn(n_min, point);
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

    /* Add edges to the minor */
    vtcs_min = snewn(n_min, vertex);
    for (i = 0; i < n_min; i++)
    {
        vx = vtcs_min + i;
        vx->idx = i;
        vx->deg = 0;
    }
    edges_min_234 = newtree234(edgecmp);
    addedges(edges_min_234, vtcs_min, pts_min, 0, n_min, 0, -1, n_min, -1, rs);

    /*
     * To create the orginal graph we need to replace all minor points by subgraphs.
     * To determine the areas in which we can place the subgraphs we need to know the
     * shortest distance between any two adjacent minor points and for every point the
     * shortest distance to the grid border. For every minor point the minimum of these
     * two values is used to calculate the radius of a circular subgraph area around
     * the minor point.
     */
    tmp = COORDMARGIN * COORDUNIT;
    tmp2 = (coord_lim - COORDMARGIN) * COORDUNIT;
    radii = snewn(n_min, long);
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        radii[i] = min(min(min(pt->x - tmp, pt->y - tmp), min(tmp2 - pt->x, tmp2 - pt->y)),
                        (squarert(square(pts_min[1].x - pts_min[0].x) + square(pts_min[1].y
                        - pts_min[0].y)) / 2)) - (SUBGRAPH_DISTANCE / 2);
        LOG(("Assigned subgraph radius %ld to subgraph %ld\n", radii[i], i));
    }

    /*
     * Assign coordinates to the subgraphs. The coordinates must lie in the previously
     * calculated subgraph areas.
     */
    tmp = SUBGRAPH_POINTENTROPY * n_sub;
    angles = snewn(tmp, double);
    for (i = 0; i < tmp; i++)
    {
        angles[i] = ((double) i * 2.0 * PI) / (double) tmp;
    }
    pts_base = snewn(n_base, point);
    for (i = 0; i < n_min; i++)
    {
        shuffle(angles, tmp, sizeof(double), rs);
        for (j = 0; j < n_sub; j++)
        {
            pt = pts_base + (i * n_sub) + j;
            pt->x = (double) pts_min[i].x + ((double) radii[i] * sin(angles[j]));
            pt->y = (double) pts_min[i].y + ((double) radii[i] * cos(angles[j]));
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
        addedges(edges_base_234, vtcs_base, pts_base,  i * n_sub, n_sub, i * n_sub, -1,
                n_min * n_sub, -1, rs);
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
        edge* _e;
        tree234* vtcs_src = newtree234(vertcmp);
        tree234* vtcs_tgt = newtree234(vertcmp);
        beste.src = beste.tgt = -1;
        for (j = e->src * n_sub; j < (e->src + 1) * n_sub; j++)
        {
            add234(vtcs_src, vtcs_base + j);
        }
        for (j = e->tgt * n_sub; j < (e->tgt + 1) * n_sub; j++)
        {
            add234(vtcs_tgt, vtcs_base + j);
        }
        for (j = 0; j < n_sub; j++)
        {
            beste.src = ((vertex*) index234(vtcs_src, j))->idx;
            for (k = 0; k < n_sub; k++)
            {
                beste.tgt = ((vertex*) index234(vtcs_tgt, k))->idx;
                /* check for crossing edges */
                for (l = 0; (_e = index234(edges_base_234, l)) != NULL; l++)
                {
                    if (_e->src == beste.src || _e->tgt == beste.src
                        || _e->src == beste.tgt || _e->tgt == beste.tgt)
                        continue;
                    else if (cross(pts_base[beste.src], pts_base[beste.tgt],
                                    pts_base[_e->src], pts_base[_e->tgt]))
                        goto next_tgt; /* this edge crosses another edge => next target */
                }
                /* check for crossing points */
                for (l = 0; l < n_min * n_sub; l++)
                {
                    if (l == beste.src || l == beste.tgt)
                        continue;
                    else if (crosspoint(pts_base[beste.src], pts_base[beste.tgt], pts_base[l]))
                        goto next_tgt; /* this edge crosses a point => next target */
                }
                addedge(edges_base_234, beste.src, beste.tgt);
                added = true;
                goto next_edge;
                next_tgt:;
            }
        }
        next_edge:
        if (!added)
        {
            int sqdist;
            int best_sqdist = square(coord_lim * COORDUNIT);
            for (j = e->src * n_sub; j < (e->src + 1) * n_sub; j++)
            {
                for (k = e->tgt * n_sub; k < (e->tgt + 1) * n_sub; k++)
                {
                    sqdist = square(pts_base[k].x - pts_base[j].x)
                            + square(pts_base[k].y - pts_base[j].y);
                    if (sqdist < best_sqdist)
                    {
                        best_sqdist = sqdist;
                        beste.src = j;
                        beste.tgt = k;
                    }
                }
            }
            addedge(edges_base_234, beste.src, beste.tgt);
        }
        LOG(("Added edge between subgraphs %d and %d (base graph vertices %d and %d)\n",
            e->src, e->tgt, beste.src, beste.tgt));
        vtcs_base[beste.src].deg++;
        vtcs_base[beste.tgt].deg++;
        LOG(("Updated degrees of vertices %d to %d and %d to %d\n", vtcs_base[beste.src].idx,
            vtcs_base[beste.src].deg, vtcs_base[beste.tgt].idx, vtcs_base[beste.tgt].deg));
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
            if (square(pt->x - p.x) + square(pt->y - p.y) < OVERLAYPOINT_TRESHOLD)
                goto next_coords; /* the point overlays another point => next coords */
        }
        pt = pts_base + tmp2 + tmp3++;
        *pt = p;
        LOG(("Assigned coordinates x:%ld, y:%ld and denominator %ld to base graph point %ld\n",
            pt->x, pt->y, pt->d, tmp2 + tmp3 - 1));
        if (tmp2 + tmp3 >= n_base) break;
        next_coords:;
    }
    sfree(coords_x);
    sfree(coords_y);

    /* Add edges to the remaining points */
    addedges(edges_base_234, vtcs_base, pts_base, 0, n_base, n_min * n_sub, -1, n_base, -1, rs);

#if DEBUG
    for (i = n_min * n_sub; i < n_base; i++)
    {
        vx = vtcs_base + i;
        assert(vx->deg > 0);
    }
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

    sfree(pts_min);
    sfree(pts_base);
    sfree(vtcs_min);
    sfree(vtcs_base);
    while ((e = delpos234(edges_min_234, 0)) != NULL) sfree(e);
    freetree234(edges_min_234);
    while ((e = delpos234(edges_base_234, 0)) != NULL) sfree(e);
    freetree234(edges_base_234);

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
    long grid_size = COORDLIMIT(params->n_base) * COORDUNIT;
    long grid_margin = COORDMARGIN * COORDUNIT;
    if ((err = validate_graph(&_desc, params->n_min, grid_size - grid_margin, grid_margin)) != NULL)
        return err;
    else if ((err = validate_graph(&_desc, params->n_base, grid_size - grid_margin, grid_margin)) != NULL)
        return err;
    else
        return NULL;
}

/*
 * Parse a graph description, i.e. a string that specifies a set of points by their
 * index, position and degree and a set of edges between those points by the point
 * indices of their incident vertices.
 */
static graph* parse_graph(const char** desc, enum grid grid, int n, long lim, long mar)
{
    int i;
    int idx, deg;
    int src, tgt;
    long x, y;
    point* pt;
    vertex* vx;
    graph* ret = snew(graph);
    ret->refcount = 1;
    ret->grid = grid;
    ret->points = snewn(n, point);
    ret->vtcs = snewn(n, vertex);
    ret->vertices = newtree234(vertcmp);
    ret->eaten_vertices = snewn(n, tree234*);
    for (i = 0; i < n; i++)
    {
        ret->eaten_vertices[i] = newtree234(vertcmp);
    }
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
    long grid_size = COORDLIMIT(params->n_base) * COORDUNIT;
    long grid_margin = COORDMARGIN * COORDUNIT;
    state->minor = parse_graph(&_desc, GRID_LEFT, params->n_min, grid_size - grid_margin, grid_margin);
    state->base = parse_graph(&_desc, GRID_RIGHT, params->n_base, grid_size - grid_margin, grid_margin);
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
    int i, j;
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

    ret->eaten_vertices = snewn(n, tree234*);
    for (i = 0; i < n; i++)
    {
        ret->eaten_vertices[i] = newtree234(vertcmp);
        for (j = 0; (vx = index234(gr->eaten_vertices[i], j)) != NULL; j++)
        {
            add234(ret->eaten_vertices[i], ret->vtcs + vx->idx);
#if DEBUG
    assert(ret->vtcs[vx->idx].idx == vx->idx);
#endif
        }
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
    int i;
    edge* e;
    (gr->refcount)--;
    if (gr->refcount <= 0)
    {
        sfree(gr->points);
        sfree(gr->vtcs);
        freetree234(gr->vertices);
        for (i = 0; i < n; i++) freetree234(gr->eaten_vertices[i]);
        sfree(gr->eaten_vertices);
        
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

    graph->points[dom].x += (graph->points[rec].x - graph->points[dom].x) / 2;
    graph->points[dom].y += (graph->points[rec].y - graph->points[dom].y) / 2;
    del234(graph->vertices, graph->vtcs + rec);
    add234(graph->eaten_vertices[dom], graph->vtcs + rec);
#if DEBUG
    assert(graph->vtcs[rec].idx == rec);
#endif
    for (i = 0; (e = delpos234(graph->edges, 0)) != NULL; i++) sfree(e);
    freetree234(graph->edges);
    graph->edges = edgescpy;
}

/*
 * Delete an edge from a graph
 */
static void delete_edge(graph* graph, edge e)
{
    /* delete the edge */
    del234(graph->edges, &e);

    /* update the involved vertices */
    del234(graph->vertices, graph->vtcs + e.src);
    del234(graph->vertices, graph->vtcs + e.tgt);
    graph->vtcs[e.src].deg--;
    graph->vtcs[e.tgt].deg--;
    add234(graph->vertices, graph->vtcs + e.src);
    add234(graph->vertices, graph->vtcs + e.tgt);
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
 * Check whether there exists an isomorphism between a test graph and a comparison
 * graph. The algorithm builds on the idea that vertices with different degree can't
 * be a vertex pair of a permutation that is an isomorphism between two graphs.
 */
static bool isomorphism_degheuristic(const graph* the_graph, const graph* cmp_graph)
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
                msrc.value = ((mapping*) find234(map, &msrc, mappingcmp))->value;
                mtgt.value = ((mapping*) find234(map, &mtgt, mappingcmp))->value;
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
                msrc.value = ((mapping*) find234(map, &msrc, mappingcmp))->value;
                mtgt.value = ((mapping*) find234(map, &mtgt, mappingcmp))->value;
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
static char* solve_bruteforce(const game_state* currstate, game_state** solvedstate,
                                int* movessize, int* moveslen, clock_t begin, int timeout)
{
#ifdef BENCHMARKS
    bool solved;
#endif
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
            if ((moves = solve_bruteforce(nextstate, solvedstate, movessize, moveslen, begin, timeout)))
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
                LOG(("Bruteforce solver - newmoves: %s\n", moves));
                (*moveslen) += movesoff;
                return moves;
            }

            free_game(nextstate);
        }  
    }
#ifdef BENCHMARKS
    else if ((solved = isomorphism_degheuristic(currstate->minor, currstate->base)))
    {
#else
    else if (isomorphism_degheuristic(currstate->minor, currstate->base))
    {
#endif
        LOG(("Bruteforce solver - Found solution\n"));
        char* moves = snewn(*movessize, char);
#ifdef BENCHMARKS
        assert(solved == isomorphism_bruteforce(currstate->minor, currstate->base));
#endif
        *moves = 'S';
        *moveslen = 1;
        *solvedstate = dup_game(currstate);
        return moves;
    }
    
    return NULL;
}

#ifdef SHARE_ADJACENT_VERTEX_HEURISTIC
/*
 * Check whether two vertices from different subgraphs are adjacent to the same
 * vertex.
 */
static bool share_adjacent_vertex(const game_state* state, int vxa, int vxb, int n1)
{
    int i;
    edge* e;
    if (vxa / n1 == vxb / n1) return false;
    for (i = 0; (e = index234(state->base->edges, i)) != NULL; i++)
    {
        if ((e->src == vxa && e->tgt != vxb && isedge(state->base->edges, vxb, e->tgt))
            || (e->tgt == vxa && e->tgt != vxb && isedge(state->base->edges, e->src, vxb)))
            return true;
    }
    return false;
}
#endif

#ifdef REPRESENT_MINOR_EDGE_HEURISTIC
/*
 * Check whether two adjacent vertices from different subgraphs represent an
 * edge of the minor.
 */
static bool represent_minor_edge(const game_state* state, int vxa, int vxb, int n1)
{
    if (vxa / n1 == vxb / n1) return false;
    return isedge(state->minor->edges, vxa / n1, vxb / n1);
}
#endif

#define EDGE_COUNT_DIFFERENCE_HEURISTIC

#ifdef EDGE_COUNT_DIFFERENCE_HEURISTIC
static int edge_count_difference(const game_state* state, int vxa, int vxb)
{
    int diff;
    game_state* tmp = dup_game(state);
    contract_edge(tmp->base, vxa, vxb);
    diff = count234(tmp->base->edges) - count234(tmp->minor->edges);
    free_game(tmp);
    return diff;
}
#endif

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

    int i, j;
    int tmp;

    point* pt;
    vertex* vx;
    edge* e;
    tree234* eaten_vtcs;

    game_state* solved;

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
    for (i = 0; (e = index234(currstate->base->edges, i)) != NULL; i++)
    {
        if (e->tgt >= n_nsub)
        {
            retoff = sprintf(buf, "%d:%d-%d;", MOVE_DELEDGE, e->src, e->tgt);
            if (retlen + retoff >= retsize)
            {
                retsize = retlen + retoff + 256;
                ret = sresize(ret, retsize, char);
            }
            strcpy(ret + retlen, buf);
            retlen += retoff;

            delete_edge(solved->base, *e);
            LOG(("Deleted edge %d-%d\n", e->src, e->tgt));
        }
    }
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
            del234(solved->base->vertices, vx);
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
        edge contre;
        contre.src = i * n_1sub;
        for (contre.src = i * n_1sub; contre.src < ((i + 1) * n_1sub) - 1; contre.src++)
        {
            for (contre.tgt = contre.src + 1; contre.tgt < (i + 1) * n_1sub; contre.tgt++)
            {
                if (isedge(solved->base->edges, contre.src, contre.tgt))
                {
                    eaten_vtcs = solved->base->eaten_vertices[contre.src];
                    if ((tmp = count234(eaten_vtcs)))
                    {
                        for (j = 0; (vx = index234(eaten_vtcs, j)) != NULL; j++)
                        {
                            if (vx->idx < n_nsub && vx->idx >= (i + 1) * n_1sub)
                            {
                                LOG(("Vertex %d has eaten vertex %d from different subgraph\n",
                                    contre.src, vx->idx));
                                goto next_src;
                            }
                        }
                    }
                    eaten_vtcs = solved->base->eaten_vertices[contre.tgt];
                    if ((tmp = count234(eaten_vtcs)))
                    {
                        for (j = 0; (vx = index234(eaten_vtcs, j)) != NULL; j++)
                        {
                            if (vx->idx < n_nsub && vx->idx >= (i + 1) * n_1sub)
                            {
                                LOG(("Vertex %d has eaten vertex %d from different subgraph\n",
                                    contre.tgt, vx->idx));
                                goto next_tgt;
                            }
                        }
                    }
                    
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
                }
                next_tgt:;
            }
            next_src:;
        }
    }
    for (i = 0; (e = index234(solved->base->edges, i)) != NULL; i++)
    {
        eaten_vtcs = solved->base->eaten_vertices[e->src];
        if ((tmp = count234(eaten_vtcs)))
        {
            for (j = 0; (vx = index234(eaten_vtcs, j)) != NULL; j++)
            {
                if (vx->idx < n_nsub
                    && ((e->src / n_1sub != e->tgt / n_1sub && vx->idx / n_1sub == e->tgt / n_1sub
                            && (edge_count_difference(solved, e->src, e->tgt) >= 0))
                        || (solved->base->vtcs[e->src].deg == 1 || solved->base->vtcs[e->tgt].deg == 1)))
                {
                    LOG(("Vertex %d has eaten vertex %d from same subgraph as vertex %d\n",
                        e->src, vx->idx, e->tgt));
                    retoff = sprintf(buf, "%d:%d-%d;", MOVE_CONTREDGE, e->src, e->tgt);
                    if (retlen + retoff >= retsize)
                    {
                        retsize = retlen + retoff + 256;
                        ret = sresize(ret, retsize, char);
                    }
                    strcpy(ret + retlen, buf);
                    retlen += retoff;

                    LOG(("Number of visible vertices is %d\n", count234(solved->base->vertices)));
                    LOG(("Contracting edge %d-%d\n", e->src, e->tgt));
                    contract_edge(solved->base, e->src, e->tgt);
                    LOG(("Number of visible vertices is %d\n", count234(solved->base->vertices)));
                    break;
                }
            }
        }
    }

    /*
     * Check whether the found solution is valid. Doing it at this point may not be optimal
     * performance wise but we can ensure that if we are able to find a solution we will
     * always find it.
     * 
     * TODO:
     * Implement another algorithm that can detect minors without any knowledge about the
     * subgraphs.
     */
    if (!(solved->solved = isomorphism_degheuristic(solved->minor, solved->base)))
    {
#ifdef BENCHMARKS
        assert(solved->solved == isomorphism_bruteforce(solved->minor, solved->base));
#endif
        sfree(ret);
        ret = NULL;
        free_game(solved);
        solved = NULL;
        retsize = 256;
        retlen = 0;
        clock_t begin = clock();
        if ((ret = solve_bruteforce(currstate, &solved, &retsize, &retlen, begin, 100)))
        {
            long coord_lim = COORDLIMIT(currstate->params.n_base);
            long circle_rad = MINORRADIUS(currstate->params.n_base);
            tmp = coord_lim - (2 * COORDMARGIN);
            LOG(("Used bruteforce solver to solve puzzle\n"));
            for (i = 0; (vx = index234(solved->base->vertices, i)) != NULL; i++)
            {
                double angle = ((double) i * 2.0 * PI) / (double) n_min;
                pt = solved->base->points + vx->idx;
                pt->x = (((double) tmp / 2.0) + ((double) circle_rad * sin(angle)) + COORDMARGIN)
                        * COORDUNIT;
                pt->y = (((double) tmp / 2.0) + ((double) circle_rad * cos(angle)) + COORDMARGIN)
                        * COORDUNIT;
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
            return ret;
        }
        *error = "Solution not known for the current puzzle state";
        return NULL;
    }

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

    /* index of the point to delete */
    int delpt;

    /* index of the resulting point of the merge - the dominant point */
    int mergept_dom;
    /* index of the recessive point */
    int mergept_rec;

    /* indices of the vertices that are incident to the edge that should be deleted */
    int deledge_src, deledge_tgt;

    /* the currently ongoing move */
    enum move current_move;

};

static game_ui *new_ui(const game_state *state)
{
    game_ui* ret = snew(game_ui);
    ret->newpt.d = 1;
    ret->dragpt = -1;
    ret->delpt = -1;
    ret->mergept_dom = -1;
    ret->mergept_rec = -1;
    ret->deledge_src = -1;
    ret->deledge_tgt = -1;
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
        case MOVE_DRAGPOINT:
            ui->dragpt = -1;
            break;
        case MOVE_CONTREDGE:
            ui->mergept_dom = -1;
            ui->mergept_rec = -1;
            break;
        case MOVE_DELPOINT:
            ui->delpt = -1;
            break;
        case MOVE_DELEDGE:
            ui->deledge_src = -1;
            ui->deledge_tgt = -1;
            break;
        case MOVE_IDLE:
        default:;
    }
    ui->current_move = MOVE_IDLE;
}

struct game_drawstate {

    int tilesize;

    long grid_size;
    long grid_margin;

    long headline_height;

    /* game window size */
    long width, height;

};

#define POINT_TRESHOLD square(POINTRADIUS + 2)
#define EDGE_TRESHOLD 2

/* heuristic to determine whether a point has been clicked */
#define point_heuristic(px, py, cx, cy) (square((px) - (cx)) + square((py) - (cy)))

/*
 * Heuristic to determine whether an edge has been clicked
 */
static double edge_heuristic(long esx, long esy, long etx, long ety,
                            long cx, long cy)
{
    long dist_st = squarert(square(etx - esx) + square(ety - esy));
    long dist_sc = squarert(square(cx - esx) + square(cy - esy));
    long dist_ct = squarert(square(etx - cx) + square(ety - cy));
    return dist_sc + dist_ct - dist_st;
}

static char *interpret_move(const game_state *state, game_ui *ui,
                            const game_drawstate *ds,
                            int x, int y, int button)
{
    /*
     * Since one can only perform moves on the base graph which is placed in GRID_RIGHT
     * and the point coordinates of a graph are given in relation to its grid one needs
     * to subtract the grid size from the x-coordinate of the mouse event. Make sure the
     * coordinates are correct, even if the player resizes the window.
     */
    long realx = (x - (state->base->grid * ds->grid_size)) * COORDUNIT / ds->tilesize;
    long realy = (y - ds->headline_height) * COORDUNIT / ds->tilesize;

    if (IS_MOUSE_DOWN(button))
    {
        int i;
        int bestpt_idx = -1;
        long ptheur;
        long bestpt_heur = POINT_TRESHOLD;
        double eheur;
        double beste_heur = EDGE_TRESHOLD;
        point* pt;
        point* pts;
        vertex* vx;
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
        pts = state->base->points;
        for (i = 0; (vx = index234(state->base->vertices, i)) != NULL; i++)
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
            else if (state->base->vtcs[bestpt_idx].deg == 0)
            {
                ui->current_move = MOVE_DELPOINT;
                ui->delpt = bestpt_idx;
                LOG(("Selected point %d to delete\n", ui->delpt));
                return UI_UPDATE;
            }
        }
        else
        {
            /*
            * Get the index of the edge with the shortest edge heuristic. Do only take edges
            * into account with a heuristic smaller than EDGE_TRESHOLD.
            */
            for (i = 0; (e = index234(state->base->edges, i)) != NULL; i++)
            {
                eheur = edge_heuristic(pts[e->src].x, pts[e->src].y, pts[e->tgt].x,
                                        pts[e->tgt].y, realx, realy);
                if (eheur < beste_heur)
                {
                    beste_heur = eheur;
                    beste = *e;
                }
            }
            /*
             * Check whether there is at least one edge with an edge heuristic smaller
             * than the treshold. If yes => game_ui requires update.
             */
            if (beste_heur < EDGE_TRESHOLD)
            {
                if (button == LEFT_BUTTON)
                {
                    ui->current_move = MOVE_CONTREDGE;
                    ui->mergept_dom = beste.src;
                    ui->mergept_rec = beste.tgt;
                    LOG(("Selected edge %d-%d to contract\n", ui->mergept_dom, ui->mergept_rec));
                    return UI_UPDATE;
                }
                else
                {
                    ui->current_move = MOVE_DELEDGE;
                    ui->deledge_src = beste.src;
                    ui->deledge_tgt = beste.tgt;
                    LOG(("Selected edge %d-%d to delete\n", ui->deledge_src, ui->deledge_tgt));
                    return UI_UPDATE;
                }
            }
        }
    }
    else if (IS_MOUSE_DRAG(button))
    {
        /* Check for an ongoing drag. If yes => game_ui requires update. */
        if (ui->current_move == MOVE_DRAGPOINT)
        {
            ui->newpt.x = realx;
            ui->newpt.y = realy;
            LOG(("Updated position of drag point to x:%ld, y:%ld\n", ui->newpt.x, ui->newpt.y));
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
                if (ui->newpt.x < 0 || ui->newpt.x > ds->grid_size * COORDUNIT / ds->tilesize
                    || ui->newpt.y < 0 || ui->newpt.y > ds->grid_size * COORDUNIT / ds->tilesize)
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
                if (x < state->base->grid * ds->grid_size
                    || x > (state->base->grid * ds->grid_size) + ds->grid_size
                    || y < ds->headline_height || y > ds->height)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->mergept_dom = -1;
                    ui->mergept_rec = -1;
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
                if (x < state->base->grid * ds->grid_size
                    || x > (state->base->grid * ds->grid_size) + ds->grid_size
                    || y < ds->headline_height || y > ds->height)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->delpt = -1;
                    LOG(("Unselected point to delete\n"));
                    return UI_UPDATE;
                }
                else
                {
                    char buf[80];
                    sprintf(buf, "%d:%d;", ui->current_move, ui->delpt);
                    LOG(("Deleting point %d\n", ui->delpt));
                    return dupstr(buf);
                }
            /*
             * Check for an ongoing edge deletion. If yes, check whether the player wants
             * to discard the deletion. If no => make move, else => discard and update
             * game_ui.
             */
            case MOVE_DELEDGE:
                if (x < state->base->grid * ds->grid_size
                    || x > (state->base->grid * ds->grid_size) + ds->grid_size
                    || y < ds->headline_height || y > ds->height)
                {
                    ui->current_move = MOVE_IDLE;
                    ui->deledge_src = -1;
                    ui->deledge_tgt = -1;
                    LOG(("Unselected edge to delete\n"));
                    return UI_UPDATE;
                }
                else
                {
                    char buf[80];
                    sprintf(buf, "%d:%d-%d;", ui->current_move, ui->deledge_src, ui->deledge_tgt);
                    LOG(("Deleting edge between vertices %d and %d\n", ui->deledge_src,
                        ui->deledge_tgt));
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
    int src, tgt;
    int current_move;
    long grid_size = COORDLIMIT(state->params.n_base) * COORDUNIT;
    long x, y;
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
                x = y = -1;
                if (sscanf(move, "%d-%ld-%ld;%n", &idx, &x, &y, &off) != 3
                    || idx < 0 || idx >= state->params.n_base
                    || x < 0 || x > grid_size || y < 0 || y > grid_size)
                {
                    LOG(("Failed to scan drag move,"));
                    if (idx >= 0) LOG((" point %d\n", idx));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }

                /* Assign the new coordinates to the dragged point */
                ret->base->points[idx].x = x;
                ret->base->points[idx].y = y;
                LOG(("Dragged point %d to position x:%ld, y:%ld\n", idx, x, y));

                move += off;
                break;
            case MOVE_CONTREDGE:
                src = tgt = -1;
                if (sscanf(move, "%d-%d;%n", &src, &tgt, &off) != 2
                    || src == tgt || !isedge(ret->base->edges, src, tgt))
                {
                    LOG(("Failed to scan contraction move,"));
                    if (src >= 0 && tgt >= 0) LOG((" edge %d-%d\n", src, tgt));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }
                
                /* Contract the edge between src and tgt */
                contract_edge(ret->base, src, tgt);
                LOG(("Contracted edge %d-%d\n", src, tgt));
                if (ret->solved) ret->solved = false;

                move += off;
                break;
            case MOVE_DELPOINT:
                idx = -1;
                if (sscanf(move, "%d;%n", &idx, &off) != 1
                    || (idx < 0 && idx >= state->params.n_base))
                {
                    LOG(("Failed to scan point deletion move,"));
                    if (idx >= 0) LOG((" point %d\n", idx));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }

                v.idx = idx;
                del234(ret->base->vertices, &v);
                LOG(("Deleted point %d\n", v.idx));
                if (ret->solved) ret->solved = false;

                move += off;
                break;
            case MOVE_DELEDGE:
                src = tgt = -1;
                if (sscanf(move, "%d-%d;%n", &src, &tgt, &off) != 2
                    || src == tgt || !isedge(ret->base->edges, src, tgt))
                {
                    LOG(("Failed to scan edge deletion move,"));
                    if (src >= 0 && tgt >= 0) LOG((" edge %d-%d\n", src, tgt));
                    else LOG(("\n"));
                    free_game(ret);
                    return NULL;
                }

                e.src = src;
                e.tgt = tgt;
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
            ret->solved = isomorphism_degheuristic(ret->minor, ret->base);
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

#define HEADLINE_HEIGHT(tilesize) (1.5F * (float) (tilesize))

static void game_compute_size(const game_params *params, int tilesize,
                              int *x, int *y)
{
    int grid_size = COORDLIMIT(params->n_base) * tilesize;
    *x = NGRIDS * grid_size;
    *y = grid_size + HEADLINE_HEIGHT(tilesize);
}

static void game_set_size(drawing *dr, game_drawstate *ds,
                          const game_params *params, int tilesize)
{
    ds->tilesize = tilesize;
    ds->grid_size = COORDLIMIT(params->n_base) * tilesize;
    ds->grid_margin = COORDMARGIN * tilesize;
    ds->width = NGRIDS * ds->grid_size;
    ds->headline_height = HEADLINE_HEIGHT(tilesize);
    ds->height = ds->grid_size + ds->headline_height;
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

    /* violet */
    ret[(COL_DELPOINT * 3) + 0] = 0.7F;
    ret[(COL_DELPOINT * 3) + 1] = 0.0F;
    ret[(COL_DELPOINT * 3) + 2] = 1.0F;

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
    ret[(COL_CONTREDGE * 3) + 0] = 1.0F;
    ret[(COL_CONTREDGE * 3) + 1] = 0.5F;
    ret[(COL_CONTREDGE * 3) + 2] = 0.0F;

    /* violet */
    ret[(COL_DELEDGE * 3) + 0] = 0.7F;
    ret[(COL_DELEDGE * 3) + 1] = 0.0F;
    ret[(COL_DELEDGE * 3) + 2] = 1.0F;

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

    /* black */
    ret[(COL_TEXT * 3) + 0] = 0.1F;
    ret[(COL_TEXT * 3) + 1] = 0.1F;
    ret[(COL_TEXT * 3) + 2] = 0.1F;

    /* light grey */
    ret[(COL_TEXTBACKGROUND * 3) + 0] = 0.7F;
    ret[(COL_TEXTBACKGROUND * 3) + 1] = 0.7F;
    ret[(COL_TEXTBACKGROUND * 3) + 2] = 0.7F;

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
    ds->grid_size = COORDLIMIT(state->params.n_base) * COORDUNIT;
    ds->grid_margin = COORDMARGIN * COORDUNIT;
    ds->width = NGRIDS * ds->grid_size;
    ds->headline_height = HEADLINE_HEIGHT(COORDUNIT);
    ds->height = ds->grid_size + ds->headline_height;

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
    long x_off;
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
    draw_rect(dr, 0, ds->headline_height, ds->width, ds->height - ds->headline_height,
                bg_color);

    /*
     * Draw the headline area, separate the two grids for the base
     * graph and minor from each other and draw an outline around
     * the area where the game is drawn.
     */
    draw_rect(dr, 0, 0, ds->width, ds->headline_height, COL_TEXTBACKGROUND);
    draw_thick_line(dr, 2.0F, (float) ds->grid_size, 0.0F, (float) ds->grid_size,
                    (float) (ds->height - 1), COL_GRIDBORDER);
    draw_rect_outline(dr, 0, 0, ds->width, ds->height, COL_OUTLINE);
    
    /*
     * Draw texts that make clear which of the two grids belongs to
     * the base graph and which to the minor.
     */
    draw_text(dr, ds->grid_size / 2, ds->headline_height / 2, FONT_FIXED,
                ds->headline_height / 2, ALIGN_VCENTRE | ALIGN_HCENTRE,
                COL_TEXT, "MINOR");
    draw_text(dr, ds->grid_size + (ds->grid_size / 2), ds->headline_height / 2,
                FONT_FIXED, ds->headline_height / 2, ALIGN_VCENTRE | ALIGN_HCENTRE,
                COL_TEXT, "ORIGINAL");
    
    /* Draw the minor edges in the intended grid */
    pts = state->minor->points;
    x_off = state->minor->grid * ds->grid_size;
    for (i = 0; (e = index234(state->minor->edges, i)) != NULL; i++)
    {
        point psrc, ptgt;
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        psrc.x = (esrc->x * ds->tilesize / COORDUNIT) + x_off;
        psrc.y = (esrc->y * ds->tilesize / COORDUNIT) + ds->headline_height;
        ptgt.x = (etgt->x * ds->tilesize / COORDUNIT) + x_off;
        ptgt.y = (etgt->y * ds->tilesize / COORDUNIT) + ds->headline_height;
        draw_line(dr, psrc.x, psrc.y, ptgt.x, ptgt.y, COL_EDGE);
    }
    /* Draw the minor points in the intended grid */
    for (i = 0; (vx = index234(state->minor->vertices, i)) != NULL; i++)
    {
        long r;
        point p;
        p.x = (pts[vx->idx].x * ds->tilesize / COORDUNIT) + x_off;
        p.y = (pts[vx->idx].y * ds->tilesize / COORDUNIT) + ds->headline_height;
        r = POINTRADIUS * ds->tilesize / COORDUNIT;
        draw_circle(dr, p.x, p.y, r, COL_MINORPOINT, COL_POINTOUTLINE);
    }

    /* Draw the base graph edges in the intended grid */
    pts = state->base->points;
    x_off = state->base->grid * ds->grid_size;
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
            psrc.x = (ui->newpt.x * ds->tilesize / COORDUNIT) + x_off;
            psrc.y = (ui->newpt.y * ds->tilesize / COORDUNIT) + ds->headline_height;
        }
        /*
         * Check whether there is a solve animation ongoing. Calculate the
         * vertex's current animated position if so.
         */
        else if (oldstate)
        {
            oesrc = oldstate->base->points + e->src;
            psrc.x = ((oesrc->x + ((float) (esrc->x - oesrc->x) * r_anim)) * ds->tilesize
                    / COORDUNIT) + x_off;
            psrc.y = ((oesrc->y + ((float) (esrc->y - oesrc->y) * r_anim)) * ds->tilesize
                    / COORDUNIT) + ds->headline_height;
        }
        else
        {
            psrc.x = (esrc->x * ds->tilesize / COORDUNIT) + x_off;
            psrc.y = (esrc->y * ds->tilesize / COORDUNIT) + ds->headline_height;
        }
        if (ui->dragpt == e->tgt)
        {
            ptgt.x = (ui->newpt.x * ds->tilesize / COORDUNIT) + + x_off;
            ptgt.y = (ui->newpt.y * ds->tilesize / COORDUNIT) + ds->headline_height;
        }
        else if (oldstate)
        {
            oetgt = oldstate->base->points + e->tgt;
            ptgt.x = ((oetgt->x + ((float) (etgt->x - oetgt->x) * r_anim)) * ds->tilesize
                    / COORDUNIT) + x_off;
            ptgt.y = ((oetgt->y + ((float) (etgt->y - oetgt->y) * r_anim)) * ds->tilesize
                    / COORDUNIT) + ds->headline_height;
        }
        else
        {
            ptgt.x = (etgt->x * ds->tilesize / COORDUNIT) + x_off;
            ptgt.y = (etgt->y * ds->tilesize / COORDUNIT) + ds->headline_height;
        }
        draw_line(dr, psrc.x, psrc.y, ptgt.x, ptgt.y,
                    (e->src == ui->mergept_dom && e->tgt == ui->mergept_rec) ?
                    COL_CONTREDGE :
                    ((e->src == ui->deledge_src && e->tgt == ui->deledge_tgt) ?
                    COL_DELEDGE :
                    ( hide ? 
                    COL_HIDEEDGE : COL_EDGE)));
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
            p.x = (ui->newpt.x * ds->tilesize / COORDUNIT) + x_off;
            p.y = (ui->newpt.y * ds->tilesize / COORDUNIT) + ds->headline_height;
        }
        else if (oldstate)
        {
            op = oldstate->base->points + vx->idx;
            p.x = ((op->x + ((float) (pts[vx->idx].x - op->x) * r_anim)) * ds->tilesize
                    / COORDUNIT) + x_off;
            p.y = ((op->y + ((float) (pts[vx->idx].y - op->y) * r_anim)) * ds->tilesize
                    / COORDUNIT) + ds->headline_height;
        }
        else
        {
            p.x = (pts[vx->idx].x * ds->tilesize / COORDUNIT) + x_off;
            p.y = (pts[vx->idx].y * ds->tilesize / COORDUNIT) + ds->headline_height;
        }
        r = POINTRADIUS * ds->tilesize / COORDUNIT;
        draw_circle(dr, p.x, p.y, r,
                    (vx->idx == ui->dragpt) ? COL_DRAGPOINT :
                    ((vx->idx == ui->delpt) ? COL_DELPOINT :
                    ( hide ? COL_HIDEPOINT :
#if DEBUG
                    (vx->idx < state->params.n_min * (state->params.n_base
                    / state->params.n_min)) ? COL_SUBPOINT : COL_BASEPOINT)),
#else
                    COL_BASEPOINT)),
#endif
                    COL_POINTOUTLINE);
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
    game_fetch_preset, NULL,                                            /* done */
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
