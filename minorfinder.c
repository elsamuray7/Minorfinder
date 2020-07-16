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
    NCOLOURS
};

/*
 * A point in a grid, spefified by its coordinates and its
 * respective denominator.
 */
typedef struct point {

    long x, y;

    long d; /* denominator */

} point;

/* A vertex that corresponds to a point in a graph. */
typedef struct vertex {

    int idx; /* points array index */

    int deg;

} vertex;

/*
 * An edge that connects to points of a graph.
 * Edges are stored in such a way that src < tgt
 * holds true.
 */
typedef struct edge {

    int src;
    int tgt;

} edge;

/*
 * An undirected graph that consists of a a set of points/vertices
 * and a set of edges that interconnect these points.
 */
typedef struct graph {

    int refcount; /* for deallocation */

    point* points;
    tree234* edges;

} graph;

struct game_params {

    int n_base; /* number of base graph points */
    int n_min; /* number of minor graph points */

};

struct game_state {

    game_params params;

    int width;
    int height;

    graph* base;
    graph* minor;

    bool solved;

};

struct game_ui {
    /* not implemented yet */
};

static game_params *default_params(void)
{
    game_params *ret = snew(game_params);

    ret->n_base = 12;
    ret->n_min = 5;

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
            n_base = 12;
            n_min = 5;
            break;
        case 1:
            n_base = 18;
            n_min = 5;
            break;
        case 2:
            n_base = 24;
            n_min = 5;
            break;
        default:
            return false;
    }

    sprintf(buf, "%d+%d points", n_base, n_min);
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
}

static char *encode_params(const game_params *params, bool full)
{
    return dupstr("FIXME");
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
#define POINTDENSITY 3
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

static void _addedge(tree234* edges, tree234* vertices, point* points, long* cnt)
{
    int i;
    vertex* vxa = delpos234(vertices, 0);
    (*cnt)--;
    /**/printf("Creating edges from vertex %d\n", vxa->idx);/**/

    for (i = (*cnt) - 1; i >= 0;)
    {
        bool crossing = false;
        int j;
        vertex* vxb = index234(vertices, i);
        edge* e;
        /**/printf("Creating edges to vertex %d\n", vxb->idx);/**/

        if (isedge(edges, vxa->idx, vxb->idx)) 
        {
            /**/printf("Found existing edge from vertex %d to vertex %d\n",
                        vxa->idx, vxb->idx);/**/
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
                /**/printf("Found crossing between edges %d-%d and %d-%d\n",
                            vxa->idx, vxb->idx, e->src, e->tgt);/**/
                crossing = true;
                break;
            }
        }

        if (crossing)
        {
            i--;
            continue;
        }
        else
        {
            /**/printf("Add edge from vertex %d to vertex %d to the edges 234-tree\n",
                        vxa->idx, vxb->idx);/**/
            addedge(edges, vxa->idx, vxb->idx);
            del234(vertices, vxb);
            vxb->deg++;
            if (vxb->deg < MAXDEGREE)
            {
                /**/printf("Update vertex %d in the vertices 234-tree\n",
                            vxb->idx);/**/
                add234(vertices, vxb);
            }
            else
            {
                /**/printf("Delete vertex %d in the vertices 234-tree\n",
                            vxb->idx);/**/
                (*cnt)--;
                i--;
            }
            vxa->deg++;
            if (vxa->deg >= MAXDEGREE) break;
        }
    }

    /**/printf("Finished creating edges from vertex %d\n", vxa->idx);/**/
}

/* denominator must divide PREFFERED_TILESIZE, i.e. 64 */
#define COORDMARGIN(l) ((l) / 16)
/* must divide the denominator of COORDMARGIN, i.e. 16 and must be a multiple of 2 */
#define COORDUNIT 8
#define square(x) ((x) * (x))

static char *new_game_desc(const game_params *params, random_state *rs,
			   char **aux, bool interactive)
{
    /*char* ret;*/
    
    const int n_min = params->n_min;
    const int n_base = params->n_base;
    const int n_sub = ((5 * n_base) / 6) / n_min;
    int* sub_sizes;

    long i, j;
    long tmp, tmp2;
    long g_size;
    long g_margin;
    long* coords_x;
    long* coords_y;
    long* radii;
    long* sub_offsets;

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
    tree234* edges_min;
    tree234* edges_base;

    /*
     * Set the grid size and margin. The grid size depends on the number of
     * points that our base graph has. It determines the actual size of our
     * grids. Together with the grid margin we get the exact area in which
     * we can place our points.
     */
    g_size = COORDLIMIT(n_base) * PREFERRED_TILESIZE;
    g_margin = COORDMARGIN(g_size);
    /**/printf("Set grid size: %d and margin: %d\n", (int) g_size, (int) g_margin);/**/

    /*
     * Generate random coordinates for the points of the minor and the base graph.
     * The coordinates will be in the range (coord_mar, coord_lim - coord_mar).
     */
    tmp = g_size - (2 * g_margin); /* coordinate range */
    tmp2 = (tmp / COORDUNIT) + 1; /* number of possible coordinates */
    coords_x = snewn(tmp2, long);
    coords_y = snewn(tmp2, long);
    for (i = 0; i <= tmp; i += COORDUNIT)
    {
        int idx = i / COORDUNIT;
        coords_x[idx] = i + g_margin;
        coords_y[idx] = i + g_margin;
    }
    shuffle(coords_x, tmp2, sizeof(*coords_x), rs);
    shuffle(coords_y, tmp2, sizeof(*coords_y), rs);

    /* Allocate memory for the points of the minor */
    pts_min = snewn(n_min, point);

    /* Assign random coordinates to the points of the minor */
    for (i = 0; i < n_min; i++)
    {
        pt = pts_min + i;
        pt->x = coords_x[i];
        pt->y = coords_y[i];
        pt->d = 1;
        /**/printf("Assigned minor point %d with coordinates x: %d, y: %d\n",
                    (int) i, (int) pt->x, (int) pt->y);/**/
    }

    /* Create the 234-tree that stores the edges of the minor */
    edges_min = newtree234(edgecmp);
    
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
    tmp = n_min; /* number of minor vertices */
    while (tmp >= 2)
    {
        _addedge(edges_min, vtcs_234, pts_min, &tmp);
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
    tmp = g_size - (2 * g_margin); /* maximal radius */
    radii = snewn(n_min, long);
    for (i = 0; i < n_min; i++)
    {
        radii[i] = tmp;
    }
    for (i = 0; i < n_min - 1; i++)
    {
        point* pta = pts_min + i;
        for (j = i + 1; j < n_min; j++)
        {
            point* ptb = pts_min + j;
            long dist = squarert(square(pta->x - ptb->x) + square(pta->y - ptb->y)) / 2;
            if (dist < radii[i])
            {
                radii[i] = dist;
                /**/printf("Assigned new minimal distance %d from minor point %d to %d\n",
                        (int) radii[i], (int) i, (int) j);/**/
            }
            if (dist < radii[j])
            {
                radii[j] = dist;
                /**/printf("Assigned new minimal distance %d from minor point %d to %d\n",
                        (int) radii[j], (int) j, (int) i);/**/
            }
        }
    }

    /* Allocare memory for the subgraphs that replace the points of the minor */
    subs = snewn(n_min, point*);

    /*
     * Generate random coordinates for the points of our subgraphs. The points
     * of a subgraph must lay in the previously calculated circle with a minor
     * point as center and the distance to the neirest neighbour of the center
     * as radius.
     */
    tmp = 4 * n_sub; /* number of possible angles */
    tmp2 = 4 * COORDUNIT; /* minimal radius */
    sub_sizes = snewn(n_min, int);
    for (i = 0; i < n_min; i++)
    {
        double angles[tmp];
        /**/printf("Creating subgraph that replaces minor point %d\n", (int) i);/**/
        if (radii[i] < tmp2 ||
            (pts_min[i].x - g_margin) < tmp2 || (g_size - g_margin - pts_min[i].x) < tmp2 ||
            (pts_min[i].y - g_margin) < tmp2 || (g_size - g_margin - pts_min[i].y) < tmp2)
        {
            subs[i] = snew(point);
            sub_sizes[i] = 1;
            sub = subs[i];
            *sub = pts_min[i];
            /**/printf("Copied minor point into subgraph and finished creating subgraph %d\n",
                        (int) i);/**/
            continue;
        }
        subs[i] = snewn(n_sub, point);
        sub_sizes[i] = n_sub;
        sub = subs[i];
        for (j = 0; j < tmp; j++)
        {
            angles[j] = (double) j * (2.0 * PI) / (double) tmp;
        }
        shuffle(angles, tmp, sizeof(double), rs);
        for (j = 0; j < n_sub; j++)
        {
            long x, y;
            do
            {
                long r = random_upto(rs, (radii[i] - ((3 * COORDUNIT) / 2)) + 1) + COORDUNIT;
                x = pts_min[i].x + (r * sin(angles[j]));
                y = pts_min[i].y + (r * cos(angles[j]));
            }
            while (x < g_margin || x > (g_size - g_margin) ||
                    y < g_margin || y > (g_size - g_margin));
            pt = sub + j;
            pt->x = x;
            pt->y = y;
            pt->d = 1;
            /**/printf("Created new subgraph point with coordinates x: %d, y: %d\n",
                    (int) pt->x, (int) pt->y);/**/
        }
        /**/printf("Finished creating subgraph %d\n", (int) i);/**/
    }
    sfree(radii);

    /* Create the 234-tree that stores the edges of the base graph */
    edges_base = newtree234(edgecmp);

    vtcs_base = snewn(n_base, vertex);

    /*
     * Add edges to the subgraphs of our base graph. Make sure that edges don't
     * cross and the degree of the vertices does not increase beyond MAXDEGREE.
     */
    sub_offsets = snewn(n_min + 1, long);
    *sub_offsets = 0;
    for (i = 0; i < n_min; i++)
    {
        vtcs_234 = newtree234(vertcmp);
        for (j = sub_offsets[i]; j < sub_offsets[i] + sub_sizes[i]; j++)
        {
            vx = vtcs_base + j;
            vx->idx = j;
            vx->deg = 0;
            add234(vtcs_234, vx);
        }
        sub_offsets[i+1] = sub_offsets[i] + sub_sizes[i];
        tmp = sub_sizes[i]; /* number of subgraph vertices */
        while (tmp >= 2)
        {
            _addedge(edges_base, vtcs_234, subs[i], &tmp);
        }
        freetree234(vtcs_234);
    }

    /* Allocate memory for the base graph points */
    pts_base = snewn(n_base, point);

    /* Copy all subgraph points into the array of base points */
    for (i = 0; i < n_min; i++)
    {
        sub = subs[i];
        for (j = 0; j < sub_sizes[i]; j++)
        {
            pts_base[sub_offsets[i]+j] = sub[j];
            /**/printf("Copied subgraph point %d of subgraph %d into base graph points at %d\n",
                        (int) j, (int) i, (int) (sub_offsets[i] + j));/**/
        }
    }

    /*
     * TODO:
     * Currently the remaining points and the subgraph points can overlap.
     * To avoid this we would have to check the remaining points for their
     * distance to any of the subgraph points. This might be complicated
     * and expensive, hence I want to see whether it's really neccessary
     * first.
     */

    /* Assign random coordinates to the remaining points */
    for (i = sub_offsets[n_min]; i < n_base; i++)
    {
       pt = pts_base + i;
       pt->x = coords_x[i + n_min - sub_offsets[n_min]];
       pt->y = coords_y[i + n_min - sub_offsets[n_min]];
       pt->d = 1;
       /**/printf("Assigned remaining point %d with coordinates x: %d, y: %d\n",
                    (int) i, (int) pt->x, (int) pt->y);/**/
    }
    sfree(coords_x);
    sfree(coords_y);

    /*
     * Add edges between the subgraphs in such a way that if the subgraphs would
     * be replaced by the minor points again the outcome would be our minor graph.
     */
    for (i = 0; (e = index234(edges_min, i)) != NULL; i++)
    {
        long idxa = random_upto(rs, sub_sizes[e->src]) + sub_offsets[e->src];
        long idxb = random_upto(rs, sub_sizes[e->tgt]) + sub_offsets[e->tgt];
        addedge(edges_base, idxa, idxb);
        /**/printf("Added edge point %d (vertex %d), subgraph %d to p:%d (v:%d), s:%d to base graph edges\n",
                    (int) (idxa - sub_offsets[e->src]), (int) idxa, e->src,
                    (int) (idxb - sub_offsets[e->tgt]), (int) idxb, e->tgt);/**/
        vtcs_base[idxa].deg++;
        vtcs_base[idxb].deg++;
    }
    sfree(sub_sizes);
    sfree(sub_offsets);

    /*
     * Add edges to the base graph including the remaining points. Make sure that
     * edges don't cross and the degree of the vertices does not increase beyond
     * MAXDEGREE.
     */
    vtcs_234 = newtree234(vertcmp);
    for (i = 0; i < sub_offsets[n_min]; i++)
    {
        vx = vtcs_base + i;
        add234(vtcs_234, vx);
    }
    for (i = sub_offsets[n_min]; i < n_base; i++)
    {
        vx = vtcs_base + i;
        vx->idx = i;
        vx->deg = 0;
        add234(vtcs_234, vx);
    }
    tmp = n_base;
    while (tmp >= 2)
    {
        _addedge(edges_base, vtcs_234, pts_base, &tmp);
    }
    sfree(vtcs_base);
    freetree234(vtcs_234);

    /*ret = NULL;
    cnt = count234(edges_min);
    {
    int retlen = 0;
    char buf[80];
    int off;
    edge* edges = snewn(cnt, edge);

    for (i = 0; (e = index234(edges_min, i)) != NULL; i++)
    {
        edges[i] = *e;
        retlen += (sprintf(buf, "%d-%d", e->src, e->tgt) + 1);
    }
    ret = snewn(retlen, char);

    e = edges;
    off = sprintf(ret, "%d-%d", e->src, e->tgt);
    for (i = 1; i < cnt; i++)
    {
        e = edges + i;
        off += sprintf(ret + off, ",%d-%d", e->src, e->tgt);
    }
    sfree(edges);
    }*/

    sfree(pts_min);
    sfree(pts_base);
    for (i = 0; i < n_min; i++) sfree(subs[i]);
    sfree(subs);
    while ((e = delpos234(edges_min, 0)) != NULL) sfree(e);
    freetree234(edges_min);
    while ((e = delpos234(edges_base, 0)) != NULL) sfree(e);
    freetree234(edges_base);

    return dupstr("FIXME");
}

static const char *validate_desc(const game_params *params, const char *desc)
{
    return NULL;
}

static game_state *new_game(midend *me, const game_params *params,
                            const char *desc)
{
    game_state *state = snew(game_state);

    /* state->FIXME = 0; */

    return state;
}

static game_state *dup_game(const game_state *state)
{
    game_state *ret = snew(game_state);

    /* ret->FIXME = state->FIXME; */

    return ret;
}

static void free_game(game_state *state)
{
    sfree(state);
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

static game_ui *new_ui(const game_state *state)
{
    return NULL;
}

static void free_ui(game_ui *ui)
{
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
}

struct game_drawstate {
    int tilesize;
    int FIXME;
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
    *x = 2 * g_size;
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

    *ncolours = NCOLOURS;
    return ret;
}

static game_drawstate *game_new_drawstate(drawing *dr, const game_state *state)
{
    struct game_drawstate *ds = snew(struct game_drawstate);

    ds->tilesize = PREFERRED_TILESIZE;
    ds->FIXME = 0;

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
    default_params,
    game_fetch_preset, NULL,
    decode_params,
    encode_params,
    free_params,
    dup_params,
    false, game_configure, custom_params,
    validate_params,
    new_game_desc,
    validate_desc,
    new_game,
    dup_game,
    free_game,
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
