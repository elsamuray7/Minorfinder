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

/* type alias for an unsigned char */
#define uint8 unsigned char

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
    COL_SOLVEFLASH,
    COL_LOSEFLASH,
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

    /* indicates whether this graph has been created by duplicating another graph */
    bool iscpy;

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
    /* player lost game */
    bool lost;

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
    if (params->n_base < DEFAULT_N_BASE)
        return "Number of base graph points is too low";
    if (params->n_min < DEFAULT_N_MIN)
        return "Number of minor points is too low";
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

#define MAXDEGREE 5
#define POINTRADIUS 5
#define POINT_CROSSCHECK_ACCURACY 50

static void addedges(tree234* edges, tree234* vertices, point* points, int n_pts, long* cnt)
{
    int i;
    vertex* vxa = delpos234(vertices, 0);
    (*cnt)--;

    for (i = (*cnt) - 1; i >= 0;)
    {
        int j, k;
        vertex* vxb = index234(vertices, i);
        edge* e;

        if (isedge(edges, vxa->idx, vxb->idx))
            goto next_vertex; /* this edge already exists => next vertex */
        /* check for crossing edges */
        for (j = 0; (e = index234(edges, j)) != NULL; j++)
        {
            if (vxa->idx == e->src || vxa->idx == e->tgt ||
                vxb->idx == e->src || vxb->idx == e->tgt)
                continue;
            else if (cross(points[vxa->idx], points[vxb->idx],
                            points[e->src], points[e->tgt]))
                goto next_vertex; /* this edge crosses another edge => next vertex */
        }
        /* check for crossing points */
        for (j = 0; j < n_pts; j++)
        {
            if (vxa->idx == j || vxb->idx == j)
            {
                continue;
            }
            else
            {
                for (k = 0; k < POINT_CROSSCHECK_ACCURACY; k++)
                {
                    double a = (double) k * (2.0 * PI) / (double) (2 * POINT_CROSSCHECK_ACCURACY);
                    point pta;
                    point ptb;
                    pta.x = points[j].x + ((POINTRADIUS + 2) * sin(a));
                    pta.y = points[j].y + ((POINTRADIUS + 2) * cos(a));
                    pta.d = 1;
                    ptb.x = points[j].x - ((POINTRADIUS + 2) * sin(a));
                    ptb.y = points[j].y - ((POINTRADIUS + 2) * cos(a));
                    ptb.d = 1;
                    if (cross(points[vxa->idx], points[vxb->idx], pta, ptb))
                        goto next_vertex; /* this edge crosses a point => next vertex */
                }
            }
        }

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
        continue;

        /*
         * I'm very sorry, I had to decide between multiple breaks, continues,
         * further effort for checking the conditions for every break and con-
         * tinue respectively and goto. The latter seemed to be the better evil.
         * Here we decrement the loop variable and continue the loop to get the
         * next vertex.
         */
        next_vertex:
        i--;
    }
}

/*
 * These parameters are highly sensitive, changing them may cause problems when
 * generating new game descriptions.
 */
#define COORDMARGIN_BASE 1
#define COORDMARGIN_MIN 3
#define COORDDENSITY_MIN 2
#define COORDLIMIT(n) ((n) - ((n) % COORDDENSITY_MIN) + (2 * COORDMARGIN_BASE) + 1)
#define COORDUNIT 32

#define square(x) ((x) * (x))

static char *new_game_desc(const game_params *params, random_state *rs,
			   char **aux, bool interactive)
{
    char* ret;
    
    const int n_min = params->n_min;
    const int n_base = params->n_base;
    const int n_sub = n_base / n_min;

    long i, j, k, l, m;
    long tmp, tmp2, tmp3;
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
     * Generate random coordinates for the minor points. The coordinates will
     * be in the following range:
     * 
     * [COORDMARGIN_MIN, coord_lim - COORDMARGIN_MIN]
     */
    tmp = coord_lim - (2 * COORDMARGIN_MIN);
    tmp2 = (tmp / COORDDENSITY_MIN) + 1;
    coords_x = snewn(tmp2, long);
    coords_y = snewn(tmp2, long);
    for (i = 0; i <= tmp; i += COORDDENSITY_MIN)
    {
        int idx = i / COORDDENSITY_MIN;
        coords_x[idx] = i + COORDMARGIN_MIN;
        coords_y[idx] = i + COORDMARGIN_MIN;
    }
    shuffle(coords_x, tmp2, sizeof(*coords_x), rs);
    shuffle(coords_y, tmp2, sizeof(*coords_y), rs);

    /* Allocate memory for the minor points */
    pts_min = snewn(n_min, point);

    /* Assign random coordinates to the minor points */
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        pt->x = coords_x[i] * COORDUNIT;
        pt->y = coords_y[i] * COORDUNIT;
        pt->d = 1;
    }
    sfree(coords_x);
    sfree(coords_y);

    /* Create the 234-tree that stores the minor edges */
    edges_min_234 = newtree234(edgecmp);
    
    /*
     * Add edges to the minor. Make sure that new edges don't cross existing ones
     * and that the degree of the vertices doesn't increase beyond MAXDEGREE. Keep
     * adding new edges until there are less than two vertices left in the vertex
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
    tmp = n_min;
    while (tmp >= 2)
    {
        addedges(edges_min_234, vtcs_234, pts_min, n_min, &tmp);
    }
    sfree(vtcs_min);
    freetree234(vtcs_234);

    /*
     * To get the orginal graph out of the minor we need to replace all points
     * of the minor by subgraphs. To determine the areas in which we can place
     * the subgraphs we need to find for every point of the minor the distance
     * to its nearest neighbour. The distances are used to calculate the radii
     * of circular subgraph areas around the minor points.
     */
    tmp = COORDMARGIN_BASE * COORDUNIT;
    tmp2 = (coord_lim - COORDMARGIN_BASE) * COORDUNIT;
    radii = snewn(n_min, long);
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        radii[i] = min(min(pt->x - tmp, tmp2 - pt->x),
                        min(pt->y - tmp, tmp2 - pt->y));
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

    /* Allocare memory for the subgraphs that replace the minor points */
    subs = snewn(n_min, point*);

    /*
     * Generate random coordinates for the subgraph points. The points of a
     * subgraph must lie in the previously calculated circle with the cor-
     * responding minor point as center and half the distance to the centers
     * neirest neighbour as radius.
     */
    tmp = 2 * n_sub;
    angles = snewn(tmp, double);
    for (i = 0; i < tmp; i++)
    {
        angles[i] = (double) i * (2.0 * PI) / (double) tmp;
    }
    for (i = 0; i < n_min; i++)
    {
        subs[i] = snewn(n_sub, point);
        sub = subs[i];
        shuffle(angles, tmp, sizeof(double), rs);
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

    /* Copy all subgraph points into the array of base graph points */
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

    /* Create the 234-tree that stores the base graph edges */
    edges_base_234 = newtree234(edgecmp);

    /*
     * Add edges between the subgraphs such that if the subgraphs would be replaced
     * by the minor points again the outcome would be the minor graph. Also make sure
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
                edge* eb;
                /* check for crossing edges */
                for (l = 0; (eb = index234(edges_base_234, l)) != NULL; l++)
                {
                    if (j == e->src || j == eb->tgt || k == eb->src || k == eb->tgt)
                        continue;
                    else if (cross(pts_base[j], pts_base[k], pts_base[eb->src], pts_base[eb->tgt]))
                        goto next_vertex; /* the edge crosses another edge => next vertex */
                }
                /* check for crossing points */
                for (l = 0; l < n_min * n_sub; l++)
                {
                    if (j == l || k == l)
                    {
                        continue;
                    }
                    else
                    {
                        for (m = 0; m < POINT_CROSSCHECK_ACCURACY; m++)
                        {
                            double a = (double) m * (2.0 * PI) / (double) (2 * POINT_CROSSCHECK_ACCURACY);
                            point pta;
                            point ptb;
                            pta.x = pts_base[l].x + ((POINTRADIUS + 2) * sin(a));
                            pta.y = pts_base[l].y + ((POINTRADIUS + 2) * cos(a));
                            pta.d = 1;
                            ptb.x = pts_base[l].x - ((POINTRADIUS + 2) * sin(a));
                            ptb.y = pts_base[l].y - ((POINTRADIUS + 2) * cos(a));
                            ptb.d = 1;
                            if (cross(pts_base[j], pts_base[k], pta, ptb))
                                goto next_vertex; /* the edge crosses a point => next vertex */
                        }
                    }
                }
                addedge(edges_base_234, j, k);
                vtcs_base[j].deg++;
                vtcs_base[k].deg++;
                goto next_edge;
                /*
                 * For consistency I decided to use goto whenever it was possible.
                 * This seemed to be the better of two evils. Here we just continue
                 * the loop to get the next vertex.
                 */
                next_vertex:;
            }
        }
        /*
         * For consistency I decided to use goto whenever it was possible.
         * This seemed to be the better of two evils. Here we just continue
         * the loop to get the next edge.
         */
        next_edge:;
    }

    /*
     * Add edges to the subgraphs. Make sure that edges don't cross and the degree
     * of the vertices does not increase beyond MAXDEGREE.
     */
    for (i = 0; i < n_min; i++)
    {
        vtcs_234 = newtree234(vertcmp);
        for (j = i * n_sub; j < (i + 1) * n_sub; j++)
        {
            add234(vtcs_234, vtcs_base + j);
        }
        tmp = n_sub;
        while (tmp >= 2)
        {
            addedges(edges_base_234, vtcs_234, pts_base, n_min * n_sub, &tmp);
        }
        freetree234(vtcs_234);
    }

    /*
     * Generate random coordinates for the remaining points. The coordinates will
     * be in the following range:
     * 
     * [COORDMARGIN_BASE, coord_lim - COORDMARGIN_BASE]
     */
    tmp = coord_lim - (2 * COORDMARGIN_BASE) + 1;
    coords_x = snewn(tmp, long);
    coords_y = snewn(tmp, long);
    for (i = 0; i < tmp; i++)
    {
        coords_x[i] = i + COORDMARGIN_BASE;
        coords_y[i] = i + COORDMARGIN_BASE;
    }
    shuffle(coords_x, tmp, sizeof(*coords_x), rs);
    shuffle(coords_y, tmp, sizeof(*coords_y), rs);

    /* Assign random coordinates to the remaining points */
    tmp2 = n_min * n_sub;
    tmp3 = 0;
    for (i = 0; i < tmp; i++)
    {
        long x = coords_x[i] * COORDUNIT;
        for (j = 0; j < tmp; j++)
        {
            long y = coords_y[j] * COORDUNIT;
            /* check for an overlaying with an edge */
            for (k = 0; (e = index234(edges_base_234, k)) != NULL; k++)
            {
                for (l = 0; l < POINT_CROSSCHECK_ACCURACY; l++)
                {
                    double a = (double) l * (2.0 * PI) / (double) (2 * POINT_CROSSCHECK_ACCURACY);
                    point pta;
                    point ptb;
                    pta.x = x + ((POINTRADIUS + 2) * sin(a));
                    pta.y = y + ((POINTRADIUS + 2) * cos(a));
                    pta.d = 1;
                    ptb.x = x - ((POINTRADIUS + 2) * sin(a));
                    ptb.y = y - ((POINTRADIUS + 2) * cos(a));
                    ptb.d = 1;
                    if (cross(pts_base[e->src], pts_base[e->tgt], pta, ptb))
                        goto next_coord_y; /* the point overlays an edge => next y-coordinate */
                }
            }
            /* check for an overlaying with a point */
            for (k = 0; k < tmp2 + tmp3; k++)
            {
                pt = pts_base + k;
                if (square(pt->x - x) + square(pt->y - y) < square((2 * POINTRADIUS) + 2))
                    goto next_coord_y; /* the point overlays another point => next y-coordinate */
            }
            pt = pts_base + (tmp2 + tmp3++);
            pt->x = x;
            pt->y = y;
            pt->d = 1;
            goto next_coord_x; /* succesfully assigned coordinates => find next coordinate pair */
            /*
             * For consistency I decided to use goto whenever it was possible.
             * This seemed to be the better of two evils. Here we just continue
             * the loop to get the next y-coordinate.
             */
            next_coord_y:;
        }
        /*
         * For consistency I decided to use goto whenever it was possible.
         * This seemed to be the better of two evils. Here we just continue
         * the loop to get the next x-coordinate.
         */
        next_coord_x:
        if (tmp2 + tmp3 >= n_base) break;
    }
    sfree(coords_x);
    sfree(coords_y);

    /*
     * Add edges to the base graph including the remaining points. Make sure that
     * edges don't cross and the degree of the vertices doesn't increase beyond
     * MAXDEGREE.
     */
    vtcs_234 = newtree234(vertcmp);
    for (i = 0; i < n_base; i++)
    {
        add234(vtcs_234, vtcs_base + i);
    }
    tmp = n_base;
    while (tmp >= 2)
    {
        addedges(edges_base_234, vtcs_234, pts_base, n_base, &tmp);
    }
    sfree(vtcs_base);
    freetree234(vtcs_234);

    /*
     * The generation of a new game description is finished. Now we need to encode
     * the description in a dynamically allocated string and connect this string to
     * the methods return value.
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
     * Calculate the length of the game description. The game description contains
     * information about the points and edges of the minor and base graph. The points
     * and edges are encoded like this:
     * 
     * (1) point: <index> - <x coordinate> - <y coordinate>
     * (2) edge: <source index> - <target index>
     * 
     * There we have our game description, a concatenation of encoded point and edge
     * data followed by an encoding of the subgraph offsets. We have four collections
     * of either points or edges that go into our description. The first two sets be-
     * long to the minor, the latter two belong to the base graph.
     * 
     * (1),(1),...;(2),(2),...;(1),(1),...;(2),(2);
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
     * that will be connected to the methods return value.
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

    /* The aux string is not required. Therefore it is set to NULL. */
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
 * Validate a graph description, i.e. a string that specifies a set of points by
 * their index and position and a set of edges between those points by the point
 * indices of their incident vertices.
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
 * Parse a graph description, i.e. a string that specifies a set of points by their
 * index and position and a set of edges between those points by the point indices
 * of their incident vertices.
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
    ret->grid = grid;
    ret->points = snewn(n, point);
    ret->idcs = snewn(n, int);
    ret->indices = newtree234(intcmp);
    ret->edges = newtree234(edgecmp);
    ret->iscpy = false;
    do
    {
        idx = atoi(*desc);
        assert(idx >= 0 && idx < n);
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        x = atol(*desc);
        assert(x >= mar && x <= lim);
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        y = atol(*desc);
        assert(y >= mar && y <= lim);
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

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
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

        assert(**desc == '-');
        (*desc)++;

        tgt = atoi(*desc);
        assert(tgt >= 0 && tgt < n);
        while (**desc && isdigit((uint8) (**desc))) (*desc)++;

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
    state->lost = false;

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
 * Duplicate a graph structure. The duplicates refcount will be 1 and it will be
 * marked as copy.
 */
static graph* dup_graph(const graph* gr)
{
    graph* ret = snew(graph);

    ret->refcount = 1;
    ret->grid = gr->grid;
    ret->points = gr->points;
    ret->idcs = gr->idcs;
    ret->indices = copytree234(gr->indices, NULL, NULL);
    ret->edges = copytree234(gr->edges, edgecpy, NULL);
    ret->iscpy = true;

    return ret;
}

static game_state *dup_game(const game_state *state)
{
    game_state *ret = snew(game_state);

    ret->params = state->params;
    ret->minor = state->minor;
    ret->minor->refcount++;
    ret->base = dup_graph(state->base);
    ret->solved = state->solved;
    ret->cheated = state->cheated;
    ret->lost = state->lost;

    return ret;
}

/*
 * Free the memory that points to a graph structure if its refcount is equal to or
 * smaller than 0. If it hasn't been created by duplicating another graph also free
 * the point array and the vertex and edge 234-trees.
 */
static void free_graph(graph* gr)
{
    edge* e;
    (gr->refcount)--;
    if (gr->refcount <= 0)
    {
        if (!gr->iscpy)
        {
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
    free_graph(state->base);
    free_graph(state->minor);
    sfree(state);
}

/*
 * Replace the given edge in the given edge 234-tree. The new edge will have
 * new_src and new_tgt as source and target respetively.
 */
static void replace_edge(tree234* edges, edge* e, int new_src, int new_tgt)
{
    del234(edges, e);

    if (new_src != new_tgt && !isedge(edges, new_src, new_tgt))
    {
        e->src = min(new_src, new_tgt);
        e->tgt = max(new_src, new_tgt);
        add234(edges, e);
    }
    else
    {
        sfree(e);
    }
}

/* 
 * Contract an edge from a graph, i.e. merge its incident vertices such that
 * no edges are lost except for the contracted edge itself.
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
            replace_edge(edgescpy, ecpy, dom, ecpy->tgt);
        else if (rec == ecpy->tgt)
            replace_edge(edgescpy, ecpy, ecpy->src, dom);
    }

    del234(graph->indices, graph->idcs + rec);
    graph->edges = edgescpy;
}

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
    if (!ui->just_merged) return;
    ui->mergept_dom = -1;
    ui->mergept_rec = -1;
    ui->just_merged = false;
}

struct game_drawstate {

    int tilesize;

};

#define MOUSEPOS_TRESHOLD (POINTRADIUS + 2)

static char *interpret_move(const game_state *state, game_ui *ui,
                            const game_drawstate *ds,
                            int x, int y, int button)
{
    if (IS_MOUSE_DOWN(button))
    {
        int i;
        int opt_idx;
        int* ix;
        long sq_dist;
        long sq_opt_dist = square(MOUSEPOS_TRESHOLD);
        point* pt;

        /*
         * Get the index of the point with the shortest distance to the coordinates
         * of the mouse click. Do only take points into account with a distance smaller
         * than MOUSEPOS_TRESHOLD.
         */
        x -= state->base->grid * ui->g_size;
        for (i = 0; (ix = (int*) index234(state->base->indices, i)) != NULL; i++)
        {
            pt = state->base->points + (*ix);
            sq_dist = square(pt->x - x) + square(pt->y - y);
            if (sq_dist < sq_opt_dist)
            {
                sq_opt_dist = sq_dist;
                opt_idx = *ix;
            }
        }

        /*
         * Check whether there is any point with a distance smaller than the treshold.
         * If yes, either update the game_ui or make a move, otherwise do nothing.
         */
        if (sq_opt_dist < square(MOUSEPOS_TRESHOLD))
        {
            if (ui->mergept_rec == -1)
            {
                ui->mergept_rec = opt_idx;
                return UI_UPDATE;
            }
            else if (ui->mergept_dom == -1)
            {
                if (ui->mergept_rec == opt_idx)
                {
                    ui->mergept_rec = -1;
                    return UI_UPDATE;
                }
                else if (isedge(state->base->edges, opt_idx, ui->mergept_rec))
                {
                    char buf[80];
                    ui->mergept_dom = opt_idx;
                    ui->just_merged = true;
                    sprintf(buf, "d%d-r%d", ui->mergept_dom, ui->mergept_rec);
                    return dupstr(buf);
                }
            }
        }
    }

    return NULL;
}

static game_state *execute_move(const game_state *state, const char *move)
{
    int dom, rec;

    /* Parse the move description. Return NULL if it is incorrect. */
    if (*move != 'd') return NULL;
    move++;
    if (!(*move || isdigit((uint8) *move))) return NULL;
    dom = atoi(move);
    while(*move && isdigit((uint8) *move)) move++;
    if (*move != '-') return NULL;
    move++;
    if (*move != 'r') return NULL;
    move++;
    if (!(*move || isdigit((uint8) *move))) return NULL;
    rec = atoi(move);

    /* Return NULL if the move is invalid */
    if (dom == rec || !isedge(state->base->edges, dom, rec)) return NULL;

    /* Copy the current game state and contract the edge between dom and rec */
    game_state* ret = dup_game(state);
    contract_edge(ret->base, dom, rec);

    return ret;
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
    ret[(COL_SOLVEFLASH * 3) + 0] = 1.0F;
    ret[(COL_SOLVEFLASH * 3) + 1] = 1.0F;
    ret[(COL_SOLVEFLASH * 3) + 2] = 1.0F;

    /* light red */
    ret[(COL_LOSEFLASH * 3) + 0] = 1.0F;
    ret[(COL_LOSEFLASH * 3) + 1] = 0.3F;
    ret[(COL_LOSEFLASH * 3) + 2] = 0.3F;

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

    return ds;
}

static void game_free_drawstate(drawing *dr, game_drawstate *ds)
{
    sfree(ds);
}

static void game_redraw(drawing *dr, game_drawstate *ds,
                        const game_state *oldstate, const game_state *state,
                        int dir, const game_ui *ui,
                        float animtime, float flashtime)
{
    int i;
    int bg_color;
    int* ix;
    long x_off;
    edge* e;
    point* esrc;
    point* etgt;
    point* pts;

    /*
     * Check whether game has been recently solved and solve function
     * hasn't been used.
     */
    if (state->solved && !oldstate->solved && !state->cheated)
        bg_color = COL_SOLVEFLASH;
    /* Check whether the game has been recently lost */
    else if (state->lost && !oldstate->lost)
        bg_color = COL_LOSEFLASH;
    else
        bg_color = COL_BACKGROUND;
    
    /*
     * The initial contents of the window are not guaranteed and
     * can vary with front ends. To be on the safe side, all games
     * should start by drawing a big background-colour rectangle
     * covering the whole window.
     */
    draw_rect(dr, 0, 0, ui->width, ui->height, bg_color);

    /*
     * Draw an outline around the area where the game is drawn and
     * separate the two grids for the base graph and minor by a thick
     * line.
     */
    draw_rect_outline(dr, 0, 0, ui->width, ui->height, COL_OUTLINE);
    draw_thick_line(dr, 2.0F, (float) ui->g_size, 0.0F, (float) ui->g_size,
                        (float) ui->g_size, COL_GRIDBORDER);
    /*
     * Draw texts that make clear which of the two grids belongs to
     * the base graph and which to the minor.
     */
    draw_text(dr, ui->g_size / 2, ui->g_margin / 2, FONT_FIXED, ui->g_margin / 2,
                ALIGN_VCENTRE | ALIGN_HCENTRE, COL_TEXT, "MINOR");
    draw_text(dr, ui->g_size + (ui->g_size / 2), ui->g_margin / 2, FONT_FIXED,
                ui->g_margin / 2, ALIGN_VCENTRE | ALIGN_HCENTRE, COL_TEXT, "ORIGINAL");

    /* Draw the minor edges in the intended grid */
    pts = state->minor->points;
    x_off = state->minor->grid * ui->g_size;
    for (i = 0; (e = index234(state->minor->edges, i)) != NULL; i++)
    {
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        draw_line(dr, esrc->x + x_off, esrc->y, etgt->x + x_off, etgt->y, COL_EDGE);
    }
    /* Draw the minor points in the intended grid */
    for (i = 0; (ix = (int*) index234(state->minor->indices, i)) != NULL; i++)
    {
        draw_circle(dr, pts[*ix].x + x_off, pts[*ix].y, POINTRADIUS, COL_MINORPOINT,
                    COL_POINTOUTLINE);
    }

    /* Draw the base graph edges in the intended grid */
    pts = state->base->points;
    x_off = state->base->grid * ui->g_size;
    for (i = 0; (e = index234(state->base->edges, i)) != NULL; i++)
    {
        esrc = pts + e->src;
        etgt = pts + e->tgt;
        draw_line(dr, esrc->x + x_off, esrc->y, etgt->x + x_off, etgt->y, COL_EDGE);
    }
    /* Draw the base graph points in the intended grid */
    for (i = 0; (ix = (int*) index234(state->base->indices, i)) != NULL; i++)
    {
        draw_circle(dr, pts[*ix].x + x_off, pts[*ix].y, POINTRADIUS, (*ix == ui->mergept_rec) ?
                    COL_MERGEPOINT : COL_BASEPOINT, COL_POINTOUTLINE);
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
    /*
     * Check whether game has been recently solved and solve function
     * hasn't been used.
     */
    if (newstate->solved && !oldstate->solved && !newstate->cheated)
        return 0.3F;
    /* Check whether the game has been recently lost */
    else if (newstate->lost && !oldstate->lost)
        return 0.3F;
    else
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
