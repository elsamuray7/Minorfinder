/*
 * nullgame.c [FIXME]: Template defining the null game (in which no
 * moves are permitted and nothing is ever drawn). This file exists
 * solely as a basis for constructing new game definitions - it
 * helps to have something which will compile from the word go and
 * merely doesn't _do_ very much yet.
 * 
 * Parts labelled FIXME actually want _removing_ (e.g. the dummy
 * field in each of the required data structures, and this entire
 * comment itself) when converting this source file into one
 * describing a real game.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include "puzzles.h"
#include "tree234.h"

/* type alias for an unsigned integer */
#define uint unsigned int

/* type alias for an unsigned long */
#define ulong unsigned long

#define PREFERRED_TILESIZE 64

enum {
    COL_BACKGROUND,
    COL_BASEPOINT,
    COL_MINORPOINT,
    COL_MERGEPOINT,
    COL_EDGE,
    COL_CONTREDGE,
    COL_FLASH,
    NCOLOURS
};

/*
 * A grid, defines a region of a window
 */
enum grid {
    GRID_LEFT,
    GRID_RIGHT,
    NGRIDS
};

enum {
    DEFAULT_N_BASE = 12,
    DEFAULT_N_MIN = 5
};

/*
 * A point in a grid, two rational coordinates and a denominator determine the position of
 * a point in a grid.
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
 * An edge that connects two vertices of a graph. Despite the fact that we use src and tgt as
 * identifiers for the vertices that are incident to an edge, edges do not have a direction.
 */
typedef struct edge {

    /* vertices that are incident to the edge */
    int src, tgt;

} edge;

/*
 * An undirected graph that consists of a set of points (vertices) and edges that connect
 * these vertices.
 */
typedef struct graph {

    /* number of references to the graph - for deallocation */
    int refcount;

    /* the grid in which the graph is drawn - determines its coordinate offset */
    enum grid grid;

    /* array of points - remains the same throughout the game */
    point* points;
    /* array of point indices - remains the same throughout the game */
    int* indices;
    /* 234-tree of point indices - maps the current game state */
    tree234* indices_234;

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

    bool solved;

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
            n_base = 16;
            n_min = 5;
            break;
        case 2:
            n_base = 20;
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
    if (*(string++) == 'b') {
        params->n_base = atoi(string);
        while (*string && isdigit((unsigned char) *string)) string++;
    }
    else
    {
        params->n_base = DEFAULT_N_BASE;
    }
    if (*(string++) == 'm')
    {
        params->n_min = atoi(string);
    }
    else
    {
        params->n_min = DEFAULT_N_MIN;
    }
}

static char *encode_params(const game_params *params, bool full)
{
    char buf[80];

    sprintf(buf, "b%dm%d", params->n_base, params->n_min);
    
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
    if (params->n_base < 12)
        return "Number of base graph points must be at least 12";
    if (params->n_min < 5)
        return "Number of minor points must be at least 5";
    /*
     * TODO:
     * If we experience performance issues we could also think of an upper limit
     * for our params.
     */
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

/*
 * Our solutions are arranged on a square grid big enough that n
 * points occupy about 1/POINTDENSITY of the grid.
 */
#define POINTDENSITY 5
#define MAXDEGREE 4
#define COORDLIMIT(n) squarert((n) * POINTDENSITY)

static void addedge(tree234 *edges, int s, int t)
{
    edge *e = snew(edge);

    assert(s != t);

    e->src = min(s, t);
    e->tgt = max(s, t);

    add234(edges, e);
}

static bool isedge(tree234 *edges, int s, int t)
{
    edge e;

    assert(s != t);

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

    if (a->deg > b->deg)
	return -1;
    else if (a->deg < b->deg)
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

static int intcmpC(const void* av, const void* bv)
{
    const int* a = (int*) av;
    const int* b = (int*) bv;

    if (*a < *b) return -1;
    else if (*a > *b) return 1;
    else return 0;
}

static int intcmp(void* av, void* bv)
{
    return intcmpC(av, bv);
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

static void addedges(tree234* edges, tree234* vertices, point* points, long* cnt)
{
    int i;
    vertex* vxa = delpos234(vertices, 0);
    (*cnt)--;

    for (i = (*cnt) - 1; i >= 0;)
    {
        bool crossing = false;
        int j;
        vertex* vxb = index234(vertices, i);
        edge* e;

        if (isedge(edges, vxa->idx, vxb->idx)) 
        {
            i--;
            continue;
        }

        for (j = 0; (e = index234(edges, j)) != NULL; j++)
        {
            if (vxa->idx == e->src || vxa->idx == e->tgt ||
                vxb->idx == e->src || vxb->idx == e->tgt)
            {
                continue;
            }
            else if (cross(points[vxa->idx], points[vxb->idx],
                            points[e->src], points[e->tgt]))
            {
                crossing = true;
                break;
            }
        }

        if (!crossing)
        {
            addedge(edges, vxa->idx, vxb->idx);
            del234(vertices, vxb);
            vxb->deg++;
            if (vxb->deg < MAXDEGREE)
            {
                add234(vertices, vxb);
            }
            else
            {
                (*cnt)--;
                i--;
            }
            vxa->deg++;
            if (vxa->deg >= MAXDEGREE) break;
        }
        else
        {
            i--;
        }
    }
}

/* denominator must divide PREFFERED_TILESIZE, i.e. 64 */
#define COORDMARGIN(l) ((l) / 16)
/* must divide the denominator of COORDMARGIN, i.e. 16 and must be a multiple of 2 */
#define COORDUNIT 8
#define square(x) ((x) * (x))

static char *new_game_desc(const game_params *params, random_state *rs,
			   char **aux, bool interactive)
{
    char* ret;
    
    const int n_min = params->n_min;
    const int n_base = params->n_base;
    const int n_sub = n_base / n_min;
    int* sub_sizes;
    int* sub_offsets;

    long i, j;
    long margin_min = 4 * COORDUNIT;
    long size;
    long count;
    long g_size;
    long g_margin;
    long* coords_x;
    long* coords_y;
    long* radii;

    point* pt;
    point* pts_min;
    point* pts_base;
    point* sub;
    point** subs;

    vertex* vx;
    vertex* vtcs_min;
    vertex* vtcs_base;
    tree234* vtcs_234;

    edge* e;
    tree234* edges_min_234;
    tree234* edges_base_234;

    /*
     * Set the grid size and margin. The grid size depends on the number of
     * points that our base graph has. It determines the actual size of our
     * grids. Together with the grid margin we get the exact area in which
     * we can place our points.
     */
    g_size = COORDLIMIT(n_base) * PREFERRED_TILESIZE;
    g_margin = COORDMARGIN(g_size);

    /*
     * Generate random coordinates for the points of the minor. The coordinates will
     * be in the range (g_margin + margin, g_size - g_margin - margin).
     */
    size = g_size - (NGRIDS * (g_margin + margin_min));
    count = (size / COORDUNIT) + 1;
    coords_x = snewn(count, long);
    coords_y = snewn(count, long);
    for (i = 0; i <= size; i += COORDUNIT)
    {
        int idx = i / COORDUNIT;
        coords_x[idx] = i + g_margin + margin_min;
        coords_y[idx] = i + g_margin + margin_min;
    }
    shuffle(coords_x, count, sizeof(*coords_x), rs);
    shuffle(coords_y, count, sizeof(*coords_y), rs);

    /* Allocate memory for the points of the minor */
    pts_min = snewn(n_min, point);

    /* Assign random coordinates to the points of the minor */
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        pt->x = coords_x[i];
        pt->y = coords_y[i];
        pt->d = 1;
    }
    sfree(coords_x);
    sfree(coords_y);

    /* Create the 234-tree that stores the edges of the minor */
    edges_min_234 = newtree234(edgecmp);
    
    /*
     * Add edges to the minor. Make sure that new edges don't cross existing ones
     * and that the degree of the vertices doesn't increase beyond MAXDEGREE. Keep
     * adding new edges untill there are less than two vertices left in the vertex
     * 234-tree.
     */
    vtcs_min = snewn(n_min, vertex);
    vtcs_234 = newtree234(vertcmp);
    for (i = 0; i < n_min; i++)
    {
        vx = vtcs_min + i;
        vx->idx = i;
        vx->deg = 0;
        add234(vtcs_234, vx);
    }
    count = n_min;
    while (count >= 2)
    {
        addedges(edges_min_234, vtcs_234, pts_min, &count);
    }
    sfree(vtcs_min);
    freetree234(vtcs_234);

    /*
     * To get the orginal graph out of the minor we need to replace all points
     * of the minor by subgraphs. To determine the areas in which we can place
     * the subgraphs we need to find for every point of the minor the distance
     * to its nearest neighbour. The distances are used to calculate the radii
     * of circular subgraph areas around the points of our minor.
     */
    size = g_size - g_margin;
    radii = snewn(n_min, long);
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        radii[i] = min(min(pt->x - g_margin, size - pt->x),
                        min(pt->y - g_margin, size - pt->y));
    }
    for (i = 0; i < n_min - 1; i++)
    {
        pt = pts_min + i;
        for (j = i + 1; j < n_min; j++)
        {
            point* ptb = pts_min + j;
            long dist = squarert(square(pt->x - ptb->x) + square(pt->y - ptb->y)) / 2;
            if (dist < radii[i]) radii[i] = dist;
            if (dist < radii[j]) radii[j] = dist;
        }
    }

    /* Allocare memory for the subgraphs that replace the points of the minor */
    subs = snewn(n_min, point*);
    sub_sizes = snewn(n_min, int);
    sub_offsets = snewn(n_min + 1, int);

    /*
     * Generate random coordinates for the points of our subgraphs. The points
     * of a subgraph must lay in the previously calculated circle with a minor
     * point as center and the distance to the neirest neighbour of the center
     * as radius.
     */
    count = 4 * n_sub;
    *sub_offsets = 0;
    for (i = 0; i < n_min; i++)
    {
        if (radii[i] < margin_min)
        {
            subs[i] = snew(point);
            sub_sizes[i] = 1;
            sub_offsets[i+1] = sub_offsets[i] + 1;
            sub = subs[i];
            *sub = pts_min[i];
        }
        else
        {
            double* angles = snewn(count, double);
            subs[i] = snewn(n_sub, point);
            sub_sizes[i] = n_sub;
            sub_offsets[i+1] = sub_offsets[i] + n_sub;
            sub = subs[i];
            for (j = 0; j < count; j++)
            {
                angles[j] = (double) j * (2.0 * PI) / (double) count;
            }
            shuffle(angles, count, sizeof(double), rs);
            for (j = 0; j < n_sub; j++)
            {
                long r = random_upto(rs, (radii[i] - ((3 * COORDUNIT) / 2)) + 1) + COORDUNIT;
                pt = sub + j;
                pt->x = pts_min[i].x + (r * sin(angles[j]));
                pt->y = pts_min[i].y + (r * cos(angles[j]));
                pt->d = 1;
            }
            sfree(angles);
        }
    }
    sfree(radii);

    /* Allocate memory for the base graph points */
    pts_base = snewn(n_base, point);

    /* Copy all subgraph points into the array of base points */
    for (i = 0; i < n_min; i++)
    {
        sub = subs[i];
        for (j = sub_offsets[i]; j < sub_offsets[i+1]; j++)
        {
            pts_base[j] = sub[j-sub_offsets[i]];
        }
    }
    for (i = 0; i < n_min; i++) sfree(subs[i]);
    sfree(subs);

    /* Create the 234-tree that stores the edges of the base graph */
    edges_base_234 = newtree234(edgecmp);

    /*
     * Add edges between the subgraphs in such a way that if the subgraphs would
     * be replaced by the minor points again the outcome would be our minor graph.
     */
    vtcs_base = snewn(n_base, vertex);
    for (i = 0; i < n_base; i++)
    {
        vx = vtcs_base + i;
        vx->idx = i;
        vx->deg = 0;
    }
    for (i = 0; (e = index234(edges_min_234, i)) != NULL; i++)
    {
        int idxa = random_upto(rs, sub_sizes[e->src]) + sub_offsets[e->src];
        int idxb = random_upto(rs, sub_sizes[e->tgt]) + sub_offsets[e->tgt];
        addedge(edges_base_234, idxa, idxb);
        vtcs_base[idxa].deg++;
        vtcs_base[idxb].deg++;
    }

    /*
     * Add edges to the subgraphs of our base graph. Make sure that edges don't
     * cross and the degree of the vertices does not increase beyond MAXDEGREE.
     */
    for (i = 0; i < n_min; i++)
    {
        vtcs_234 = newtree234(vertcmp);
        for (j = sub_offsets[i]; j < sub_offsets[i+1]; j++)
        {
            add234(vtcs_234, vtcs_base + j);
        }
        count = sub_sizes[i];
        while (count >= 2)
        {
            addedges(edges_base_234, vtcs_234, pts_base, &count);
        }
        freetree234(vtcs_234);
    }

    /*
     * TODO:
     * Currently the remaining points and the subgraph points can overlap.
     * To avoid this we would have to check the remaining points for their
     * distance to any of the subgraph points. This might be complicated
     * and expensive, hence I want to see whether it's really neccessary
     * first.
     */

    /*
     * Generate random coordinates for the remaining points. The coordinates
     * will be in the range (g_margin, g_size - g_margin).
     */
    size = g_size - (NGRIDS * g_margin);
    count = (size / COORDUNIT) + 1;
    coords_x = snewn(count, long);
    coords_y = snewn(count, long);
    for (i = 0; i <= size; i += COORDUNIT)
    {
        int idx = i / COORDUNIT;
        coords_x[idx] = i + g_margin;
        coords_y[idx] = i + g_margin;
    }
    shuffle(coords_x, count, sizeof(*coords_x), rs);
    shuffle(coords_y, count, sizeof(*coords_y), rs);

    /* Assign random coordinates to the remaining points */
    count = n_base - sub_offsets[n_min];
    for (i = 0; i < count; i++)
    {
       pt = pts_base + i + sub_offsets[n_min];
       pt->x = coords_x[i];
       pt->y = coords_y[i];
       pt->d = 1;
    }
    sfree(coords_x);
    sfree(coords_y);

    /*
     * Add edges to the base graph including the remaining points. Make sure that
     * edges don't cross and the degree of the vertices does not increase beyond
     * MAXDEGREE.
     */
    vtcs_234 = newtree234(vertcmp);
    for (i = 0; i < n_base; i++)
    {
        add234(vtcs_234, vtcs_base + i);
    }
    count = n_base;
    while (count >= 2)
    {
        addedges(edges_base_234, vtcs_234, pts_base, &count);
    }
    sfree(vtcs_base);
    freetree234(vtcs_234);

    /*
     * The generation of a new game description is finished. Now we need to encode
     * the description in a dynamically allocated string and connect this string to
     * the return value of our method.
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
     * Calculate the length of our game description. The game description contains
     * information about the points and edges of our minor and base graph. The points
     * and edges are encoded like this:
     * 
     * (1) point: <index> - <x coordinate> - <y coordinate>
     * (2) edge: <source index> - <target index>
     * 
     * There we have our game description, a concatenation of encoded point and edge
     * data. We have four collections of either points or edges that go into our
     * description. The first two sets belong to the minor, the last two belong to
     * the base graph. Single points and edges are separated by commas whereas sets
     * are separated by semicolons:
     * 
     * (1),(1),...;(2),(2),...;(1),(1),...;(2),(2),...
     */
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        len += (sprintf(buf, "%ld-%ld-%ld", i, pt->x, pt->y) + 1);
    }
    for (i = 0; (e = index234(edges_min_234, i)) != NULL; i++)
    {
        edges_min[i] = *e;
        len += (sprintf(buf, "%d-%d", e->src, e->tgt) + 1);
    }
    for (i = 0; i < n_base; i++)
    {
        pt = pts_base + i;
        len += (sprintf(buf, "%ld-%ld-%ld", i, pt->x, pt->y) + 1);
    }
    for (i = 0; (e = index234(edges_base_234, i)) != NULL; i++)
    {
        edges_base[i] = *e;
        len += (sprintf(buf, "%d-%d", e->src, e->tgt) + 1);
    }

    /*
     * Allocate memory for len+1 chars, that is exactly the length of our game
     * description including a trailing '\0'.
     */
    ret = snewn(++len, char);
    
    /*
     * Now encode the game description and write it into the allocated string
     * that will be connected to our return value.
     */
    for (i = 0; i < n_min - 1; i++)
    {
        pt = pts_min + i;
        off += sprintf(ret + off, "%ld-%ld-%ld%s", i, pt->x, pt->y, sep);
    }
    sep = ";";
    pt = pts_min + n_min - 1;
    off += sprintf(ret + off, "%d-%ld-%ld%s", n_min - 1, pt->x, pt->y, sep);
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
        pt = pts_base + i;
        off += sprintf(ret + off, "%ld-%ld-%ld%s", i, pt->x, pt->y, sep);
    }
    sep = ";";
    pt = pts_base + n_base - 1;
    off += sprintf(ret + off, "%d-%ld-%ld%s", n_base - 1, pt->x, pt->y, sep);
    sep = ",";
    for (i = 0; i < count_base - 1; i++)
    {
        e = edges_base + i;
        off += sprintf(ret + off, "%d-%d%s", e->src, e->tgt, sep);
    }
    sep = ";";
    e = edges_base + count_base - 1;
    off += sprintf(ret + off, "%d-%d%s", e->src, e->tgt, sep);

    sfree(edges_min);
    sfree(edges_base);
    }

    /*
     * Last but not least we have to fill our aux string with some useful
     * information that will help us to solve the game, i.e. information
     * about the subgraphs that replaced the minor points in our base graph.
     */
    *aux = NULL;
    {
    const char* sep;
    char* str;
    char buf[80];
    int len = 0;
    int off = 0;

    /*
     * Calculate the length of our aux string. It will contain a sequence
     * of encoded subgraph points. The subgraph points are separated by
     * commas and the subgraphs themselves are separated by semicolons.
     * The encoded string looks like this:
     * 
     * <subgraph index> : <index in subgraph> - <index in base points> , ... ; ...
     */
    for (i = 0; i < n_min; i++)
    {
        len += (sprintf(buf, "%ld", i) + 1);
        for (j = sub_offsets[i]; j < sub_offsets[i+1]; j++)
        {
            len += (sprintf(buf, "%ld-%ld", j - sub_offsets[i], j) + 1);
        }
    }

    /*
     * Allocate memory for len+1 chars, that is exactly the length of our aux
     * string including a trailing '\0'.
     */
    str = snewn(++len, char);

    /*
     * Now encode the subgraph data and write it into the allocated aux string.
     */
    for (i = 0; i < n_min; i++)
    {
        sep = ":";
        off += sprintf(str + off, "%ld%s", i, sep);
        sep = ",";
        for (j = sub_offsets[i]; j < sub_offsets[i+1] - 1; j++)
        {
            off += sprintf(str + off, "%ld-%ld%s", j - sub_offsets[i], j, sep);
        }
        sep = ";";
        off += sprintf(str + off, "%d-%d%s", (sub_offsets[i+1] - 1) - sub_offsets[i], sub_offsets[i+1] - 1, sep);
    }

    *aux = str;
    }

    sfree(pts_min);
    sfree(pts_base);
    sfree(sub_sizes);
    sfree(sub_offsets);
    while ((e = delpos234(edges_min_234, 0)) != NULL) sfree(e);
    freetree234(edges_min_234);
    while ((e = delpos234(edges_base_234, 0)) != NULL) sfree(e);
    freetree234(edges_base_234);

    return ret;
}

/*
 * Validate a graph description, i.e. a string that specifies a set of points by their index and
 * position and a set of edges between those points by the point indices of their incident vertices.
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
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after point index in game description";
        (*desc)++;

        x = atol(*desc);
        if (x < mar || x > lim)
            return "X-coordinate out of range in game description";
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after point index in game description";
        (*desc)++;

        y = atol(*desc);
        if (y < mar || y > lim)
            return "Y-coordinate out of range in game description";
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        if (**desc != ',' && **desc != ';')
            return "Expected ',' or ';' after y-coordinate in game description";
        if (*((*desc)++) == ';') break;
    }
    while (**desc)
    {
        src = atoi(*desc);
        if (src < 0 || src >= n)
            return "Edge source index out of range in game description";
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        if (**desc != '-')
            return "Expected '-' after edge source index in game description";
        (*desc)++;

        tgt = atoi(*desc);
        if (tgt < 0 || tgt >= n)
            return "Edge target index out of range in game description";
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

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
    long g_size = COORDLIMIT(params->n_base) * PREFERRED_TILESIZE;
    long g_margin = COORDMARGIN(g_size);
    if ((err = validate_graph(&_desc, params->n_min, g_size - g_margin, g_margin)) != NULL)
        return err;
    else if ((err = validate_graph(&_desc, params->n_base, g_size - g_margin, g_margin)) != NULL)
        return err;
    else
        return NULL;
    
}

/*
 * Parse a graph description that has been validated by validate_graph. The graph description consists
 * of a description of a set of points (vertices) and a description of a set of edges.
 */
static graph* parse_graph(const char** desc, int n, long lim, long mar)
{
    graph* ret = snew(graph);
    point* pt;
    int idx;
    int src, tgt;
    int* ix;
    long x, y;
    ret->refcount = 1;
    ret->points = snewn(n, point);
    ret->indices = snewn(n, int);
    ret->indices_234 = newtree234(intcmp);
    ret->edges = newtree234(edgecmp);
    do
    {
        idx = atoi(*desc);
        /**/printf("%d-", idx);/**/
        assert(idx >= 0 && idx < n);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        x = atol(*desc);
        /**/printf("%ld-", x);/**/
        assert(x >= mar && x <= lim);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        y = atol(*desc);
        /**/printf("%ld\n", y);/**/
        assert(y >= mar && y <= lim);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        pt = ret->points + idx;
        pt->x = x;
        pt->y = y;
        pt->d = 1;

        ix = ret->indices + idx;
        *ix = idx;
        add234(ret->indices_234, ix);

        assert(**desc == ',' || **desc == ';');
    }
    while (*((*desc)++) != ';');
    do
    {
        src = atoi(*desc);
        /**/printf("%d-", src);/**/
        assert(src >= 0 && src < n);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        tgt = atoi(*desc);
        /**/printf("%d\n", tgt);/**/
        assert(tgt >= 0 && tgt < n);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        addedge(ret->edges, src, tgt);

        assert(**desc == ',' || **desc == ';');
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
    long g_size = COORDLIMIT(params->n_base) * PREFERRED_TILESIZE;
    long g_margin = COORDMARGIN(g_size);
    state->minor = parse_graph(&_desc, params->n_min, g_size - g_margin, g_margin);
    state->base = parse_graph(&_desc, params->n_base, g_size - g_margin, g_margin);
    state->solved = false;

    return state;
}

static game_state *dup_game(const game_state *state)
{
    game_state *ret = snew(game_state);

    ret->params = state->params;
    ret->minor = state->minor;
    ret->minor->refcount++;
    ret->base = state->base;
    ret->base->refcount++;
    ret->solved = state->solved;

    return ret;
}

/*
 * Free the memory that points to a graph structure if its refcount is equal to or smaller than 0,
 * inlcudes freeing the point array and the vertex and edge 234-trees.
 */
static void free_graph(graph* graph)
{
    edge* e;
    graph->refcount--;
    if (graph->refcount <= 0)
    {
        sfree(graph->points);
        sfree(graph->indices);
        freetree234(graph->indices_234);
        while((e = delpos234(graph->edges, 0)) != NULL) sfree(e);
        freetree234(graph->edges);
        sfree(graph);
    }
}

static void free_game(game_state *state)
{
    free_graph(state->minor);
    free_graph(state->base);
    sfree(state);
}

/*
 * Replace the given edge in the given edge 234-tree. The new edge will have new_src and new_tgt 
 * as source and target respetively.
 */
/**/static void replace_edge(tree234* edges, edge* e, int new_src, int new_tgt)
{
    del234(edges, e);

    if (new_src == new_tgt)
    {
        sfree(e);
    }
    else
    {
        int s = min(new_src, new_tgt);
        int t = max(new_src, new_tgt);
        e->src = s;
        e->tgt = t;
        if (find234(edges, e, NULL) == NULL)
            add234(edges, e);
        else
            sfree(e);
    }
}/**/


/* 
 * Contract an edge from our graph, i.e. merge its incident vertices in such a way that no edges
 * are lost except for the contracted edge itself.
 */
/**/static void contract_edge(graph* graph, int src, int tgt)
{
    int i;
    int s, t;
    edge* e;

    assert(src != tgt);
    assert(isedge(graph->edges, src, tgt));

    s = min(src, tgt);
    t = max(src, tgt);
     
    for (i = 0; (e = index234(graph->edges, i)) != NULL; i++)
    {
        if (t == e->src)
        {
            replace_edge(graph->edges, e, s, e->tgt);
        }
        else if (t == e->tgt)
        {
            replace_edge(graph->edges, e, e->src, s);
        }
    }

    del234(graph->indices_234, graph->indices + t);
}/**/

static char *solve_game(const game_state *state, const game_state *currstate,
                        const char *aux, const char **error)
{
    return NULL;
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

    /* game window size */
    long width, height;

    /* indices of the points that should be merged */
    int mergept1, mergept2;
    /* the new coordinates of the resulting point */
    long new_x, new_y;

    bool just_merged;

};

static game_ui *new_ui(const game_state *state)
{
    game_ui* ret = snew(game_ui);
    long g_size = COORDLIMIT(state->params.n_base) * PREFERRED_TILESIZE;
    ret->width = NGRIDS * g_size;
    ret->height = g_size;
    ret->mergept1 = -1;
    ret->mergept2 = -1;
    ret->new_x = -1;
    ret->new_y = -1;
    ret->just_merged = false;

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
    ui->mergept1 = -1;
    ui->mergept2 = -1;
    ui->new_x = -1;
    ui->new_y = -1;
}

struct game_drawstate {

    int tilesize;

    /* indices of the points that should be merged */
    long mergept1, mergept2;

    graph* base;
    graph* minor;

};

static char *interpret_move(const game_state *state, game_ui *ui,
                            const game_drawstate *ds,
                            int x, int y, int button)
{
    return NULL;
}

static game_state *execute_move(const game_state *state, const char *move)
{
    return NULL;
}

/* ----------------------------------------------------------------------
 * Drawing routines.
 */

static void game_compute_size(const game_params *params, int tilesize,
                              int *x, int *y)
{
    int g_size = COORDLIMIT(params->n_base) * tilesize;
    *x = NGRIDS * g_size;
    *y = g_size;
}

static void game_set_size(drawing *dr, game_drawstate *ds,
                          const game_params *params, int tilesize)
{
    ds->tilesize = tilesize;
}

static float *game_colours(frontend *fe, int *ncolours)
{
    float *ret = snewn(3 * NCOLOURS, float);

    frontend_default_colour(fe, &ret[COL_BACKGROUND * 3]);

    /* blue */
    ret[(COL_BASEPOINT * 3) + 0] = 0.0F;
    ret[(COL_BASEPOINT * 3) + 1] = 0.0F;
    ret[(COL_BASEPOINT * 3) + 2] = 1.0F;

    /* red */
    ret[(COL_MINORPOINT * 3) + 0] = 0.0F;
    ret[(COL_MINORPOINT * 3) + 1] = 0.0F;
    ret[(COL_MINORPOINT * 3) + 2] = 1.0F;

    /* yellow */
    ret[(COL_MERGEPOINT * 3) + 0] = 1.0F;
    ret[(COL_MERGEPOINT * 3) + 1] = 0.8F;
    ret[(COL_MERGEPOINT * 3) + 2] = 0.0F;

    /* black */
    ret[(COL_EDGE * 3) + 0] = 0.0F;
    ret[(COL_EDGE * 3) + 1] = 0.0F;
    ret[(COL_EDGE * 3) + 2] = 0.0F;

    /* grey */
    ret[(COL_CONTREDGE * 3) + 0] = 0.5F;
    ret[(COL_CONTREDGE * 3) + 1] = 0.5F;
    ret[(COL_CONTREDGE * 3) + 2] = 0.5F;

    /* white */
    ret[(COL_FLASH * 3) + 0] = 1.0F;
    ret[(COL_FLASH * 3) + 1] = 1.0F;
    ret[(COL_FLASH * 3) + 2] = 1.0F;

    *ncolours = NCOLOURS;
    return ret;
}

static game_drawstate *game_new_drawstate(drawing *dr, const game_state *state)
{
    struct game_drawstate *ds = snew(struct game_drawstate);

    ds->tilesize = PREFERRED_TILESIZE;
    ds->mergept1 = -1;
    ds->mergept2 = -1;
    ds->base = state->base;
    state->base->refcount++;
    ds->minor = state->minor;
    state->minor->refcount++;

    return ds;
}

static void game_free_drawstate(drawing *dr, game_drawstate *ds)
{
    free_graph(ds->minor);
    free_graph(ds->base);
    sfree(ds);
}

static void game_redraw(drawing *dr, game_drawstate *ds,
                        const game_state *oldstate, const game_state *state,
                        int dir, const game_ui *ui,
                        float animtime, float flashtime)
{
    /*
     * The initial contents of the window are not guaranteed and
     * can vary with front ends. To be on the safe side, all games
     * should start by drawing a big background-colour rectangle
     * covering the whole window.
     */
    draw_rect(dr, 0, 0, 10*ds->tilesize, 10*ds->tilesize, COL_BACKGROUND);
    draw_update(dr, 0, 0, 10*ds->tilesize, 10*ds->tilesize);
}

static float game_anim_length(const game_state *oldstate,
                              const game_state *newstate, int dir, game_ui *ui)
{
    return 0.0F;
}

static float game_flash_length(const game_state *oldstate,
                               const game_state *newstate, int dir, game_ui *ui)
{
    return 0.0F;
}

static int game_status(const game_state *state)
{
    return 0;
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
    default_params, /* done */
    game_fetch_preset, NULL, /* done */
    decode_params, /* done */
    encode_params, /* done */
    free_params, /* done */
    dup_params, /* done */
    false, game_configure, custom_params,
    validate_params, /* done */
    new_game_desc, /* done */
    validate_desc, /* done */
    new_game, /* done */
    dup_game, /* done */
    free_game, /* done */
    false, solve_game,
    false, game_can_format_as_text_now, game_text_format,
    new_ui,
    free_ui,
    encode_ui,
    decode_ui,
    NULL, /* game_request_keys */
    game_changed_state,
    interpret_move,
    execute_move,
    PREFERRED_TILESIZE, game_compute_size, game_set_size,
    game_colours,
    game_new_drawstate,
    game_free_drawstate,
    game_redraw,
    game_anim_length,
    game_flash_length,
    game_status,
    false, false, game_print_size, game_print,
    false,			       /* wants_statusbar */
    false, game_timing_state,
    0,				       /* flags */
};
