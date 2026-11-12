/* sr3d.c -- 3D scale-rich model with plates: iterated division of convex
 * polyhedra (single-threaded).
 *
 * Step t = 1..T: pick a free cell with probability proportional to its volume,
 * a uniform point x in it and a plane normal n (isotropic, or one of the three
 * axes); the plate is the section of the cell by the plane through x.  It ends
 * on every face it crosses (contacts along lines).  With thickness
 * lam_t = lam0 t^-alpha the two daughters are the parent clipped by the planes
 * offset by +-lam_t/2 (the free region shrinks; empty daughters are dropped).
 * Dual graph: plate = node, contact = edge.  Faces carry the label of the plate
 * they lie on; cube faces are labelled -1..-6.
 *
 * Faces are kept with outward (counter-clockwise from outside) orientation, so
 * each crossed face contributes exactly one "exit" point (inside -> outside)
 * and the m exit points are the m vertices of the section polygon.
 *
 * usage: sr3d T seed prefix [orient=iso|axes] [lam0=0] [alpha=0]
 * output (binary, little endian):
 *   prefix.steps : T records of 10 float64:
 *       t, depth, nfaces_parent, m (faces crossed), ncube (cube faces crossed),
 *       plate area, plate perimeter, parent volume, u (volume fraction of the
 *       daughter on the +n side; zero thickness), free volume after the step
 *     (t = -1 marks a rejected step)
 *   prefix.edges : records int32 t, int32 label, float32 face area,
 *       float32 face perimeter, float32 cell volume, int32 cell nfaces
 *   prefix.cells : final cells, records int32 nfaces, int32 nvert_total
 *       (sum of face sizes = 2E), float64 volume, int32 touches_cube, int32 depth
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/* ---------- RNG: xoshiro256** ---------- */
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t next(void) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}
static double urand(void) { return (next() >> 11) * 0x1.0p-53; }
static double nrand(void) {
    double u1 = urand(), u2 = urand();
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}
static void seed_rng(uint64_t seed) {
    uint64_t z = seed * 0x9E3779B97F4A7C15ULL + 12345;
    for (int i = 0; i < 4; i++) {
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t x = z;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        s[i] = x ^ (x >> 31);
    }
}

/* ---------- Fenwick tree ---------- */
static int fn;
static double *ft, *fw, ftot;
static void fen_init(int n) { fn = n; ft = calloc(n + 1, sizeof(double)); fw = calloc(n, sizeof(double)); ftot = 0; }
static void fen_set(int i, double w) {
    double d = w - fw[i]; fw[i] = w; ftot += d;
    for (i++; i <= fn; i += i & -i) ft[i] += d;
}
static int fen_find(double u) {
    int pos = 0, step = 1;
    while (step * 2 <= fn) step *= 2;
    for (; step; step >>= 1) {
        int nx = pos + step;
        if (nx <= fn && ft[nx] <= u) { pos = nx; u -= ft[nx]; }
    }
    return pos < fn ? pos : fn - 1;
}

/* ---------- cells ---------- */
typedef struct {
    int nf, nv;      /* faces, total vertex slots */
    int *lab, *fs;   /* labels, face sizes */
    double *x;       /* 3*nv coordinates, faces concatenated */
    double vol;
    int depth;
} Cell;

static Cell *cells;
static int ncells;

#define MAXF 256
#define MAXV 4096
/* scratch for building a cell */
static int b_nf, b_nv, b_lab[MAXF], b_fs[MAXF];
static double b_x[3 * MAXV];

static void cell_from_build(Cell *c) {
    c->nf = b_nf; c->nv = b_nv;
    c->lab = malloc(sizeof(int) * 2 * b_nf);
    c->fs = c->lab + b_nf;
    memcpy(c->lab, b_lab, sizeof(int) * b_nf);
    memcpy(c->fs, b_fs, sizeof(int) * b_nf);
    c->x = malloc(sizeof(double) * 3 * b_nv);
    memcpy(c->x, b_x, sizeof(double) * 3 * b_nv);
}
static void cell_free(Cell *c) { free(c->lab); free(c->x); c->lab = NULL; c->x = NULL; c->nf = c->nv = 0; c->vol = 0; }

static inline double dot(const double *a, const double *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline void cross(const double *a, const double *b, double *o) {
    o[0] = a[1]*b[2] - a[2]*b[1]; o[1] = a[2]*b[0] - a[0]*b[2]; o[2] = a[0]*b[1] - a[1]*b[0];
}

static void centroid(int nv, const double *x, double *c) {
    c[0] = c[1] = c[2] = 0;
    for (int i = 0; i < nv; i++) { c[0] += x[3*i]; c[1] += x[3*i+1]; c[2] += x[3*i+2]; }
    c[0] /= nv; c[1] /= nv; c[2] /= nv;
}

/* volume by fan from the vertex centroid */
static double volume(int nf, const int *fs, const double *x, int nv) {
    double c[3]; centroid(nv, x, c);
    double V = 0; int off = 0;
    for (int f = 0; f < nf; f++) {
        const double *a = x + 3 * off;
        for (int j = 1; j < fs[f] - 1; j++) {
            const double *b = x + 3 * (off + j), *d = x + 3 * (off + j + 1);
            double u[3] = {a[0]-c[0], a[1]-c[1], a[2]-c[2]};
            double v[3] = {b[0]-c[0], b[1]-c[1], b[2]-c[2]};
            double w[3] = {d[0]-c[0], d[1]-c[1], d[2]-c[2]};
            double cr[3]; cross(v, w, cr);
            V += fabs(dot(u, cr)) / 6.0;
        }
        off += fs[f];
    }
    return V;
}

static void sample_point(const Cell *K, double *p) {
    double c[3]; centroid(K->nv, K->x, c);
    double u = urand() * K->vol, acc = 0;
    int off = 0, found = 0;
    const double *A = NULL, *B = NULL, *D = NULL;
    for (int f = 0; f < K->nf && !found; f++) {
        const double *a = K->x + 3 * off;
        for (int j = 1; j < K->fs[f] - 1; j++) {
            const double *b = K->x + 3 * (off + j), *d = K->x + 3 * (off + j + 1);
            double uu[3] = {a[0]-c[0], a[1]-c[1], a[2]-c[2]};
            double v[3] = {b[0]-c[0], b[1]-c[1], b[2]-c[2]};
            double w[3] = {d[0]-c[0], d[1]-c[1], d[2]-c[2]};
            double cr[3]; cross(v, w, cr);
            acc += fabs(dot(uu, cr)) / 6.0;
            A = a; B = b; D = d;   /* last tet seen: fallback for rounding */
            if (u <= acc) { found = 1; break; }
        }
        off += K->fs[f];
    }
    double s1 = urand(), t1 = urand(), w1 = urand();
    if (s1 + t1 > 1) { s1 = 1 - s1; t1 = 1 - t1; }
    if (t1 + w1 > 1) { double tmp = w1; w1 = 1 - s1 - t1; t1 = 1 - tmp; }
    else if (s1 + t1 + w1 > 1) { double tmp = w1; w1 = s1 + t1 + w1 - 1; s1 = 1 - t1 - tmp; }
    for (int k = 0; k < 3; k++)
        p[k] = c[k] + s1 * (A[k] - c[k]) + t1 * (B[k] - c[k]) + w1 * (D[k] - c[k]);
}

/* Section of K by plane n.x = c0: exit points of crossed faces (one per face).
 * Fills pts (3*m), lab (m), returns m.  Exit = from d>=0 to d<0 in the face
 * orientation. */
static int section(const Cell *K, const double *n, double c0, double *pts, int *labs, int *fidx) {
    int off = 0, m = 0;
    for (int f = 0; f < K->nf; f++) {
        int sz = K->fs[f];
        const double *X = K->x + 3 * off;
        for (int i = 0; i < sz; i++) {
            int j = (i + 1) % sz;
            double di = dot(n, X + 3*i) - c0, dj = dot(n, X + 3*j) - c0;
            if (di >= 0 && dj < 0) {
                double t = di / (di - dj);
                for (int k = 0; k < 3; k++) pts[3*m + k] = X[3*i + k] + t * (X[3*j + k] - X[3*i + k]);
                labs[m] = K->lab[f]; fidx[m] = f; m++;
                break; /* convex face: one exit */
            }
        }
        off += sz;
    }
    return m;
}

/* sort section points CCW about n (or about -n if rev) */
static void order_poly(double *P, int m, const double *n, int rev, int *labs, int *fidx) {
    double c[3]; centroid(m, P, c);
    double e1[3] = {P[0]-c[0], P[1]-c[1], P[2]-c[2]};
    double l = sqrt(dot(e1, e1)) + 1e-300; e1[0] /= l; e1[1] /= l; e1[2] /= l;
    double e2[3]; cross(n, e1, e2);
    double ang[MAXF]; int idx[MAXF];
    for (int i = 0; i < m; i++) {
        double d[3] = {P[3*i]-c[0], P[3*i+1]-c[1], P[3*i+2]-c[2]};
        ang[i] = atan2(dot(d, e2), dot(d, e1)); idx[i] = i;
    }
    for (int i = 1; i < m; i++) { /* insertion sort */
        int k = idx[i]; double a = ang[k]; int j = i - 1;
        while (j >= 0 && ang[idx[j]] > a) { idx[j+1] = idx[j]; j--; }
        idx[j+1] = k;
    }
    double Q[3 * MAXF]; int L[MAXF], FI[MAXF];
    for (int i = 0; i < m; i++) {
        int k = rev ? idx[m - 1 - i] : idx[i];
        memcpy(Q + 3*i, P + 3*k, 3 * sizeof(double));
        if (labs) L[i] = labs[k];
        if (fidx) FI[i] = fidx[k];
    }
    memcpy(P, Q, 3 * m * sizeof(double));
    if (labs) memcpy(labs, L, m * sizeof(int));
    if (fidx) memcpy(fidx, FI, m * sizeof(int));
}

static void poly_area_perim(const double *P, int m, double *A, double *per) {
    double s[3] = {0, 0, 0}; *per = 0;
    for (int i = 1; i < m - 1; i++) {
        double u[3] = {P[3*i]-P[0], P[3*i+1]-P[1], P[3*i+2]-P[2]};
        double v[3] = {P[3*(i+1)]-P[0], P[3*(i+1)+1]-P[1], P[3*(i+1)+2]-P[2]};
        double c[3]; cross(u, v, c); s[0] += c[0]; s[1] += c[1]; s[2] += c[2];
    }
    *A = 0.5 * sqrt(dot(s, s));
    for (int i = 0; i < m; i++) {
        int j = (i + 1) % m;
        double d[3] = {P[3*j]-P[3*i], P[3*j+1]-P[3*i+1], P[3*j+2]-P[3*i+2]};
        *per += sqrt(dot(d, d));
    }
}

/* clip K to n.x >= c0; new face gets label newlab.  Result in build buffers.
 * returns 1 if a non-empty polyhedron results. */
static int clip(const Cell *K, const double *n, double c0, int newlab) {
    b_nf = 0; b_nv = 0;
    int off = 0, m = 0;
    double cut[3 * MAXF];
    for (int f = 0; f < K->nf; f++) {
        int sz = K->fs[f];
        const double *X = K->x + 3 * off;
        int start = b_nv;
        for (int i = 0; i < sz; i++) {
            int j = (i + 1) % sz;
            double di = dot(n, X + 3*i) - c0, dj = dot(n, X + 3*j) - c0;
            if (di >= 0) { memcpy(b_x + 3*b_nv, X + 3*i, 3*sizeof(double)); b_nv++; }
            if ((di >= 0) != (dj >= 0)) {
                double t = di / (di - dj);
                double *q = b_x + 3*b_nv;
                for (int k = 0; k < 3; k++) q[k] = X[3*i + k] + t * (X[3*j + k] - X[3*i + k]);
                b_nv++;
                if (di >= 0 && m < MAXF) { memcpy(cut + 3*m, q, 3*sizeof(double)); m++; }
            }
        }
        int cnt = b_nv - start;
        if (cnt >= 3) { b_lab[b_nf] = K->lab[f]; b_fs[b_nf] = cnt; b_nf++; }
        else b_nv = start;
        off += sz;
        if (b_nv > MAXV - 64 || b_nf > MAXF - 4) { fprintf(stderr, "overflow\n"); exit(1); }
    }
    if (m >= 3) {
        double nn[3] = {-n[0], -n[1], -n[2]};
        order_poly(cut, m, nn, 0, NULL, NULL); /* CCW about outward normal -n */
        memcpy(b_x + 3*b_nv, cut, 3*m*sizeof(double));
        b_lab[b_nf] = newlab; b_fs[b_nf] = m; b_nf++; b_nv += m;
    }
    return b_nf >= 4;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: sr3d T seed prefix [iso|axes] [lam0] [alpha]\n"); return 1; }
    long T = atol(argv[1]);
    uint64_t seed = strtoull(argv[2], NULL, 10);
    const char *pre = argv[3];
    int axes = argc > 4 && strcmp(argv[4], "axes") == 0;
    double lam0 = argc > 5 ? atof(argv[5]) : 0.0;
    double alpha = argc > 6 ? atof(argv[6]) : 0.0;
    seed_rng(seed);
    char fn_[1024];
    snprintf(fn_, sizeof fn_, "%s.steps", pre); FILE *fs_ = fopen(fn_, "wb");
    snprintf(fn_, sizeof fn_, "%s.edges", pre); FILE *fe = fopen(fn_, "wb");

    long cap = T + 2;
    cells = calloc(cap, sizeof(Cell));
    fen_init((int)cap);
    /* unit cube, outward CCW faces */
    double V[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    int F[6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},{2,3,7,6},{1,2,6,5},{0,4,7,3}};
    b_nf = 6; b_nv = 24;
    for (int f = 0; f < 6; f++) { b_lab[f] = -1 - f; b_fs[f] = 4;
        for (int j = 0; j < 4; j++) memcpy(b_x + 3*(4*f + j), V[F[f][j]], 3*sizeof(double)); }
    cell_from_build(&cells[0]); cells[0].vol = 1.0; cells[0].depth = 0; ncells = 1;
    fen_set(0, 1.0);

    double pts[3 * MAXF]; int labs[MAXF], fidx[MAXF];
    long nrej = 0;
    for (long t = 1; t <= T; t++) {
        double rec[10];
        if (ftot < 1e-14) { /* jammed */
            for (long tt = t; tt <= T; tt++) { double r[10] = {-2,0,0,0,0,0,0,0,0,ftot}; fwrite(r, sizeof r, 1, fs_); }
            break;
        }
        int ci = fen_find(urand() * ftot);
        Cell *K = &cells[ci];
        if (K->nf == 0) { t--; continue; }
        double p[3]; sample_point(K, p);
        double n[3];
        if (axes) { int a = (int)(urand() * 3); n[0] = n[1] = n[2] = 0; n[a] = 1; }
        else { do { n[0] = nrand(); n[1] = nrand(); n[2] = nrand(); } while (dot(n, n) < 1e-20);
               double l = sqrt(dot(n, n)); n[0] /= l; n[1] /= l; n[2] /= l; }
        double c0 = dot(n, p);
        int m = section(K, n, c0, pts, labs, fidx);
        if (m < 3) { nrej++; double r[10] = {-1,0,0,0,0,0,0,0,0,ftot}; fwrite(r, sizeof r, 1, fs_); continue; }
        order_poly(pts, m, n, 0, labs, fidx);
        double area, per; poly_area_perim(pts, m, &area, &per);
        double lam = lam0 > 0 ? lam0 * pow((double)t, -alpha) : 0.0;
        /* face area/perimeter of the hit faces (before the split) */
        int ncube = 0;
        {
            int offs[MAXF]; int o = 0;
            for (int f = 0; f < K->nf; f++) { offs[f] = o; o += K->fs[f]; }
            for (int i = 0; i < m; i++) {
                int f = fidx[i]; double fa, fp;
                poly_area_perim(K->x + 3*offs[f], K->fs[f], &fa, &fp);
                int32_t a32[2] = {(int32_t)t, labs[i]};
                float f32[3] = {(float)fa, (float)fp, (float)K->vol};
                int32_t nfK = K->nf;
                fwrite(a32, 4, 2, fe); fwrite(f32, 4, 3, fe); fwrite(&nfK, 4, 1, fe);
                if (labs[i] < 0) ncube++;
            }
        }
        double Vp = K->vol; int dep = K->depth, nfp = K->nf;
        Cell D1 = {0}, D2 = {0};
        int ok1 = clip(K, n, c0 + lam / 2, (int)t);
        if (ok1) { cell_from_build(&D1); D1.vol = volume(D1.nf, D1.fs, D1.x, D1.nv); D1.depth = dep + 1; }
        double nm[3] = {-n[0], -n[1], -n[2]};
        int ok2 = clip(K, nm, -c0 + lam / 2, (int)t);
        if (ok2) { cell_from_build(&D2); D2.vol = volume(D2.nf, D2.fs, D2.x, D2.nv); D2.depth = dep + 1; }
        double v1 = ok1 ? D1.vol : 0, v2 = ok2 ? D2.vol : 0;
        cell_free(K);
        if (ok1) { cells[ci] = D1; fen_set(ci, D1.vol); } else fen_set(ci, 0);
        if (ok2) { cells[ncells] = D2; fen_set(ncells, D2.vol); } else { cells[ncells].nf = 0; }
        ncells++;
        rec[0] = t; rec[1] = dep; rec[2] = nfp; rec[3] = m; rec[4] = ncube;
        rec[5] = area; rec[6] = per; rec[7] = Vp; rec[8] = (v1 + v2) > 0 ? v1 / (v1 + v2) : 0; rec[9] = ftot;
        fwrite(rec, sizeof rec, 1, fs_);
        if (t % 200000 == 0) { fprintf(stderr, "t=%ld free=%.4g rej=%ld\n", t, ftot, nrej); }
    }
    fclose(fs_); fclose(fe);
    snprintf(fn_, sizeof fn_, "%s.cells", pre); FILE *fc = fopen(fn_, "wb");
    for (int i = 0; i < ncells; i++) {
        Cell *K = &cells[i];
        if (K->nf == 0) continue;
        int tc = 0; for (int f = 0; f < K->nf; f++) if (K->lab[f] < 0) tc = 1;
        int32_t a[2] = {K->nf, K->nv}; double v = K->vol; int32_t b[2] = {tc, K->depth};
        fwrite(a, 4, 2, fc); fwrite(&v, 8, 1, fc); fwrite(b, 4, 2, fc);
    }
    fclose(fc);
    fprintf(stderr, "done T=%ld rejected=%ld cells=%d\n", T, nrej, ncells);
    return 0;
}
