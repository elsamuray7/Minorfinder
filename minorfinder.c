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

enum {
    COL_BACKGROUND,
    COL_OUTLINE,
    COL_GRIDBORDER,
    COL_BASEPOINT,
    COL_MINORPOINT,
    COL_MERGEPOINT,
    COL_POINTOUTLINE,
    COL_EDGE,
    COL_CONTREDGE,
    COL_FLASH,
    COL_TEXT,
    NCOLOURS
};

enum default_params {
    DEFAULT_N_BASE = 8,
    DEFAULT_N_MIN = 3
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

    /* number of copies of the graph - for deallocation */
    int* cpycount;

    /* the grid in which the graph is drawn - determines its coordinate offset */
    enum grid grid;

    /* array of points - remains the same throughout the game */
    point* points;
    /* array of point indices - remains the same throughout the game */
    int* idcs;
    /* 234-tree of point indices - maps the current game state */
    tree234* indices;

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
            n_base = 11;
            n_min = 4;
            break;
        case 2:
            n_base = 15;
            n_min = 4;
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

#define MAXDEGREE 4

static void addedges(tree234* edges, tree234* vertices, point* points, int n_pts, long* cnt)
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
        /* check for crossing edges */
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
        /* check for crossing points */
        for (j = 0; j < n_pts; j++)
        {
            if (vxa->idx == j || vxb->idx == j)
            {
                continue;
            }
            else if (cross(points[vxa->idx], points[vxb->idx], points[j], points[j]))
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

#define COORDMARGIN_BASE 1
#define COORDMARGIN_MIN 2
#define COORDDENSITY_MIN 2
#define COORDLIMIT(n) ((n) - ((n) % COORDDENSITY_MIN) + (2 * COORDMARGIN_BASE) + 1)
#define COORDUNIT 32

#define POINTRADIUS 5

#define square(x) ((x) * (x))

static char *new_game_desc(const game_params *params, random_state *rs,
			   char **aux, bool interactive)
{
    char* ret;
    
    const int n_min = params->n_min;
    const int n_base = params->n_base;
    const int n_sub = n_base / n_min;

    long i, j, k, l;
    long count;
    long upper_lim, lower_lim;
    long coord_lim;
    long* coords_x;
    long* coords_y;
    long* radii;

    double* angles;

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

    coord_lim = COORDLIMIT(n_base);

    /*
     * Generate random coordinates for the points of the minor. The coordinates will
     * be in the range (g_margin + margin, g_size - g_margin - margin).
     */
    upper_lim = coord_lim - (2 * COORDMARGIN_MIN);
    count = (upper_lim / COORDDENSITY_MIN) + 1;
    coords_x = snewn(count, long);
    coords_y = snewn(count, long);
    for (i = 0; i <= upper_lim; i += COORDDENSITY_MIN)
    {
        int idx = i / COORDDENSITY_MIN;
        coords_x[idx] = i + COORDMARGIN_MIN;
        coords_y[idx] = i + COORDMARGIN_MIN;
    }
    shuffle(coords_x, count, sizeof(*coords_x), rs);
    shuffle(coords_y, count, sizeof(*coords_y), rs);

    /* Allocate memory for the points of the minor */
    pts_min = snewn(n_min, point);

    /* Assign random coordinates to the points of the minor */
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        pt->x = coords_x[i] * COORDUNIT;
        pt->y = coords_y[i] * COORDUNIT;
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
        addedges(edges_min_234, vtcs_234, pts_min, n_min, &count);
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
    upper_lim = (coord_lim - COORDMARGIN_BASE) * COORDUNIT;
    lower_lim = COORDMARGIN_BASE * COORDUNIT;
    radii = snewn(n_min, long);
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        radii[i] = min(min(pt->x - lower_lim, upper_lim - pt->x),
                        min(pt->y - lower_lim, upper_lim - pt->y));
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

    /*
     * Generate random coordinates for the points of our subgraphs. The points
     * of a subgraph must lie in the previously calculated circle with a minor
     * point as center and the distance to the neirest neighbour of the center
     * as radius.
     */
    count = 2 * n_sub;
    angles = snewn(count, double);
    for (i = 0; i < count; i++)
    {
        angles[i] = (double) i * (2.0 * PI) / (double) count;
    }
    for (i = 0; i < n_min; i++)
    {
        subs[i] = snewn(n_sub, point);
        sub = subs[i];
        shuffle(angles, count, sizeof(double), rs);
        for (j = 0; j < n_sub; j++)
        {
            long r = random_upto(rs, radii[i] - (3 * POINTRADIUS) + 1) + (2 * POINTRADIUS);
            pt = sub + j;
            pt->x = pts_min[i].x + (r * sin(angles[j]));
            pt->y = pts_min[i].y + (r * cos(angles[j]));
            pt->d = 1;
        }
    }
    sfree(angles);
    sfree(radii);

    /* Allocate memory for the base graph points */
    pts_base = snewn(n_base, point);

    /* Copy all subgraph points into the array of base points */
    for (i = 0; i < n_min; i++)
    {
        sub = subs[i];
        for (j = i * n_sub; j < (i + 1) * n_sub; j++)
        {
            pts_base[j] = sub[j-(i*n_sub)];
        }
    }
    for (i = 0; i < n_min; i++) sfree(subs[i]);
    sfree(subs);

    /* Create the 234-tree that stores the edges of the base graph */
    edges_base_234 = newtree234(edgecmp);

    /*
     * Add edges between the subgraphs such that if the subgraphs would be replaced
     * by the minor points again the outcome would be our minor graph. Also make sure
     * that the edges don't cross existing ones. Therefore we just iterate over all
     * vertices of the source and target subgraph respectively until we find a non-
     * corssing edge. This works because the circular areas in which the subgraphs
     * lie don't intersect. Hence there must exist a non-corssing edge between two
     * subgraphs if the corresponding minor points share an edge.
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
        for (j = e->src * n_sub; j < (e->src + 1) * n_sub; j++)
        {
            for (k = e->tgt * n_sub; k < (e->tgt + 1) * n_sub; k++)
            {
                bool crossing = false;
                edge* eb;
                /* check for crossing edges */
                for (l = 0; (eb = index234(edges_base_234, l)) != NULL; l++)
                {
                    if (j == e->src || j == eb->tgt || k == eb->src || k == eb->tgt)
                    {
                        continue;
                    }
                    else if (cross(pts_base[j], pts_base[k], pts_base[eb->src], pts_base[eb->tgt]))
                    {
                        crossing = true;
                        break;
                    }
                }
                /* check for crossing points */
                for (l = 0; l < n_min * n_sub; l++)
                {
                    if (j == l || k == l)
                    {
                        continue;
                    }
                    else if (cross(pts_base[j], pts_base[k], pts_base[l], pts_base[l]))
                    {
                        crossing = true;
                        break;
                    }
                }
                if (!crossing)
                {
                    addedge(edges_base_234, j, k);
                    vtcs_base[j].deg++;
                    vtcs_base[k].deg++;
                    goto next_iter;
                }
            }
        }
        next_iter:;
    }

    /*
     * Add edges to the subgraphs of our base graph. Make sure that edges don't
     * cross and the degree of the vertices does not increase beyond MAXDEGREE.
     */
    for (i = 0; i < n_min; i++)
    {
        vtcs_234 = newtree234(vertcmp);
        for (j = i * n_sub; j < (i + 1) * n_sub; j++)
        {
            add234(vtcs_234, vtcs_base + j);
        }
        count = n_sub;
        while (count >= 2)
        {
            addedges(edges_base_234, vtcs_234, pts_base, n_min * n_sub, &count);
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
    count = coord_lim - (2 * COORDMARGIN_BASE) + 1;
    coords_x = snewn(count, long);
    coords_y = snewn(count, long);
    for (i = 0; i < count; i++)
    {
        coords_x[i] = i + COORDMARGIN_BASE;
        coords_y[i] = i + COORDMARGIN_BASE;
    }
    shuffle(coords_x, count, sizeof(*coords_x), rs);
    shuffle(coords_y, count, sizeof(*coords_y), rs);

    /* Assign random coordinates to the remaining points */
    count = n_base - (n_min * n_sub);
    for (i = 0; i < count; i++)
    {
       pt = pts_base + (n_base - count + i);
       pt->x = coords_x[i] * COORDUNIT;
       pt->y = coords_y[i] * COORDUNIT;
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
        addedges(edges_base_234, vtcs_234, pts_base, n_base, &count);
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
     * Furthermore the game description contain information about the subgraph offsets.
     * This information is required to determine whether the player has solved the game
     * or not. The subgraph offsets are encoded like this:
     * 
     * (3) <subgraph index> - <subgraph offset>
     * 
     * There we have our game description, a concatenation of encoded point and edge
     * data followed by an encoding of the subgraph offsets. We have four collections
     * of either points or edges that go into our description. The first two sets
     * belong to the minor, the latter two belong to the base graph. These four sets
     * are followed by a collection of subgraph offsets. Single points, edges and
     * subgraph offsets are separated by commas whereas whole sets are separated by
     * semicolons:
     * 
     * (1),(1),...;(2),(2),...;(1),(1),...;(2),(2),...;(3),(3),...;
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
    sprintf(ret + off, "%d-%d%s", e->src, e->tgt, sep);

    sfree(edges_min);
    sfree(edges_base);
    }

    /* The aux string is not required anymore. Therefore it is set to NULL. */
    *aux = NULL;

    sfree(pts_min);
    sfree(pts_base);
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
    long g_size = COORDLIMIT(params->n_base) * COORDUNIT;
    long g_margin = COORDMARGIN_BASE * COORDUNIT;
    if ((err = validate_graph(&_desc, params->n_min, g_size - g_margin, g_margin)) != NULL)
        return err;
    else if ((err = validate_graph(&_desc, params->n_base, g_size - g_margin, g_margin)) != NULL)
        return err;
    else
        return NULL;
}

/*
 * Parse a graph description that has been validated by validate_graph
 */
static graph* parse_graph(const char** desc, enum grid grid, int n, long lim, long mar)
{
    int idx;
    int src, tgt;
    int* ix;
    long x, y;
    point* pt;
    graph* ret = snew(graph);
    ret->refcount = 1;
    ret->cpycount = snew(int);
    *ret->cpycount = 1;
    ret->grid = grid;
    ret->points = snewn(n, point);
    ret->idcs = snewn(n, int);
    ret->indices = newtree234(intcmp);
    ret->edges = newtree234(edgecmp);
    do
    {
        idx = atoi(*desc);
        assert(idx >= 0 && idx < n);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        x = atol(*desc);
        assert(x >= mar && x <= lim);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        y = atol(*desc);
        assert(y >= mar && y <= lim);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        pt = ret->points + idx;
        pt->x = x;
        pt->y = y;
        pt->d = 1;

        ix = ret->idcs + idx;
        *ix = idx;
        add234(ret->indices, ix);

        assert(**desc == ',' || **desc == ';');
    }
    while (*((*desc)++) != ';');
    do
    {
        src = atoi(*desc);
        assert(src >= 0 && src < n);
        while (**desc && isdigit((unsigned char) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        tgt = atoi(*desc);
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
    long g_size = COORDLIMIT(params->n_base) * COORDUNIT;
    long g_margin = COORDMARGIN_BASE * COORDUNIT;
    state->minor = parse_graph(&_desc, GRID_LEFT, params->n_min, g_size - g_margin, g_margin);
    state->base = parse_graph(&_desc, GRID_RIGHT, params->n_base, g_size - g_margin, g_margin);
    state->solved = false;
    state->cheated = false;

    return state;
}

/*
 * Duplicate a graph structure. Set its refcount to 1 and increment its copycount.
 */
static graph* dup_graph(const graph* gr)
{
    int i;
    int* ix;
    edge* e;
    edge* ecpy;
    graph* ret = snew(graph);

    ret->refcount = 1;
    ret->cpycount = gr->cpycount;
    (*ret->cpycount)++;
    ret->grid = gr->grid;
    ret->points = gr->points;
    ret->idcs = gr->idcs;
    ret->indices = newtree234(intcmp);
    for (i = 0; (ix = index234(gr->indices, i)) != NULL; i++) add234(ret->indices, ix);
    ret->edges = newtree234(edgecmp);
    for (i = 0; (e = index234(gr->edges, i)) != NULL; i++)
    {
        ecpy = snew(edge);
        *ecpy = *e;
        add234(ret->edges, ecpy);
    }

    return ret;
}

static game_state *dup_game(const game_state *state)
{
    game_state *ret = snew(game_state);

    ret->params = state->params;
    ret->minor = dup_graph(state->minor);
    ret->base = dup_graph(state->base);
    ret->solved = state->solved;
    ret->cheated = state->cheated;

    return ret;
}

/*
 * Free the memory that points to a graph structure if its refcount is equal to or smaller than 0,
 * inlcudes freeing the point array and the vertex and edge 234-trees.
 */
static void free_graph(graph* gr, bool iscpy)
{
    edge* e;
    if (iscpy) (*gr->cpycount)--;
    (gr->refcount)--;
    if (gr->refcount <= 0)
    {
        if (*gr->cpycount <= 0)
        {
            sfree(gr->cpycount);
            sfree(gr->points);
            sfree(gr->idcs);
        }
        freetree234(gr->indices);
        while((e = delpos234(gr->edges, 0)) != NULL) sfree(e);
        freetree234(gr->edges);
        sfree(gr);
    }
}

static void free_game(game_state *state)
{
    free_graph(state->base, true);
    free_graph(state->minor, true);
    sfree(state);
}

/*
 * Replace the given edge in the given edge 234-tree. The new edge will have new_src and new_tgt 
 * as source and target respetively.
 */
/*static void replace_edge(tree234* edges, edge* e, int new_src, int new_tgt)
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
}*/


/* 
 * Contract an edge from our graph, i.e. merge its incident vertices in such a way that no edges
 * are lost except for the contracted edge itself.
 */
/*static void contract_edge(graph* graph, int src, int tgt)
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
            replace_edge(graph->edges, e, s, e->tgt);
        else if (t == e->tgt)
            replace_edge(graph->edges, e, e->src, s);
    }

    del234(graph->indices, graph->idcs + t);
}*/

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

    /* grid size */
    long g_size;
    /* grid margin*/
    long g_margin;

    /* game window size */
    long width, height;

    /* index of the resulting point of the points merge - the dominant point */
    int mergept_dom;
    /* index of the recessive point */
    int mergept_rec;

    bool just_merged;

};

static game_ui *new_ui(const game_state *state)
{
    game_ui* ret = snew(game_ui);
    ret->g_size = COORDLIMIT(state->params.n_base) * COORDUNIT;
    ret->g_margin = COORDMARGIN_BASE * COORDUNIT;
    ret->width = NGRIDS * ret->g_size;
    ret->height = ret->g_size;
    ret->mergept_dom = -1;
    ret->mergept_rec = -1;
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
    ui->mergept_dom = -1;
    ui->mergept_rec = -1;
}

struct game_drawstate {

    int tilesize;

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

    /* black */
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
    ret[(COL_MERGEPOINT * 3) + 0] = 0.0F;
    ret[(COL_MERGEPOINT * 3) + 1] = 1.0F;
    ret[(COL_MERGEPOINT * 3) + 2] = 0.0F;

    /* black */
    ret[(COL_POINTOUTLINE * 3) + 0] = 0.0F;
    ret[(COL_POINTOUTLINE * 3) + 1] = 0.0F;
    ret[(COL_POINTOUTLINE * 3) + 2] = 0.0F;

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

    /* black */
    ret[(COL_TEXT * 3) + 0] = 0.1F;
    ret[(COL_TEXT * 3) + 1] = 0.1F;
    ret[(COL_TEXT * 3) + 2] = 0.1F;

    *ncolours = NCOLOURS;
    return ret;
}

static game_drawstate *game_new_drawstate(drawing *dr, const game_state *state)
{
    struct game_drawstate *ds = snew(struct game_drawstate);

    ds->tilesize = COORDUNIT;
    ds->base = state->base;
    (ds->base->refcount)++;
    ds->minor = state->minor;
    (ds->minor->refcount)++;

    return ds;
}

static void game_free_drawstate(drawing *dr, game_drawstate *ds)
{
    free_graph(ds->minor, false);
    free_graph(ds->base, false);
    sfree(ds);
}

static void game_redraw(drawing *dr, game_drawstate *ds,
                        const game_state *oldstate, const game_state *state,
                        int dir, const game_ui *ui,
                        float animtime, float flashtime)
{
    int i;
    int idx;
    long x_off;
    edge* e;
    point* esrc;
    point* etgt;
    point* pts;
    
    /*
     * The initial contents of the window are not guaranteed and
     * can vary with front ends. To be on the safe side, all games
     * should start by drawing a big background-colour rectangle
     * covering the whole window.
     */
    draw_rect(dr, 0, 0, ui->width, ui->height, COL_BACKGROUND);

    draw_rect_outline(dr, 0, 0, ui->width, ui->height, COL_OUTLINE);
    draw_thick_line(dr, 2.0F, (float) ui->g_size, 0.0F, (float) ui->g_size,
                    (float) ui->g_size, COL_GRIDBORDER);
    draw_text(dr, ui->g_size / 2, ui->g_margin / 2, FONT_FIXED, ui->g_margin / 2,
                ALIGN_VCENTRE | ALIGN_HCENTRE, COL_TEXT, "MINOR");
    draw_text(dr, ui->g_size + (ui->g_size / 2), ui->g_margin / 2, FONT_FIXED, ui->g_margin / 2,
                ALIGN_VCENTRE | ALIGN_HCENTRE, COL_TEXT, "ORIGINAL");

    pts = ds->base->points;
    x_off = ds->base->grid * ui->g_size;
    for (i = 0; (e = index234(ds->base->edges, i)) != NULL; i++)
    {
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        draw_line(dr, esrc->x + x_off, esrc->y, etgt->x + x_off, etgt->y, COL_EDGE);
    }
    for (i = 0; index234(ds->base->indices, i) != NULL; i++)
    {
        idx = *((int*) index234(ds->base->indices, i));
        draw_circle(dr, pts[idx].x + x_off, pts[idx].y, POINTRADIUS,
                    COL_BASEPOINT, COL_POINTOUTLINE);
    }

    pts = ds->minor->points;
    x_off = ds->minor->grid * ui->g_size;
    for (i = 0; (e = index234(ds->minor->edges, i)) != NULL; i++)
    {
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        draw_line(dr, esrc->x + x_off, esrc->y, etgt->x + x_off, etgt->y, COL_EDGE);
    }
    for (i = 0; index234(ds->minor->indices, i) != NULL; i++)
    {
        idx = *((int*) index234(ds->minor->indices, i));
        draw_circle(dr, pts[idx].x + x_off, pts[idx].y, POINTRADIUS,
                    COL_MINORPOINT, COL_POINTOUTLINE);
    }

    draw_update(dr, 0, 0, ui->width, ui->height);
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
    false, solve_game,
    false, game_can_format_as_text_now, game_text_format,
    new_ui,                                                             /* done */
    free_ui,                                                            /* done */
    encode_ui,
    decode_ui,
    NULL, /* game_request_keys */
    game_changed_state,
    interpret_move,
    execute_move,
    COORDUNIT, game_compute_size, game_set_size,                        /* done */
    game_colours,                                                       /* done */
    game_new_drawstate,                                                 /* done */
    game_free_drawstate,                                                /* done */
    game_redraw,
    game_anim_length,
    game_flash_length,
    game_status,
    false, false, game_print_size, game_print,
    false,			       /* wants_statusbar */
    false, game_timing_state,
    0,				       /* flags */
};
