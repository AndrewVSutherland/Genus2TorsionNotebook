/*
 * General 5-torsion sieve on the compact two-parameter M(12) chart.
 *
 * This is deliberately broader than the point-contact-5 cover.  For
 *
 *   L = b + (2b-1)x,
 *   A = x + w(1+bx),
 *   F = L*(L*A^2 + 4b(1+x)^2*(wL-x^2)),
 *
 * every good rational (b,w) carries the marked M(12) class.  Any rational
 * 5-torsion class (with arbitrary Mumford support) forces
 *
 *   5 | #Jac(y^2=F)(F_p)
 *
 * at every good p != 2,5.  Bad reductions and denominator residues are
 * passed conservatively.  The program computes all masks independently by
 * counting C(F_p) and C(F_{p^2}), then scans reduced rational pairs.
 *
 * Usage:
 *   clang -O3 -Wall -Wextra -o /tmp/m12_general5_bw_sieve \
 *       code/m12_general5_bw_height_sieve.c
 *   /tmp/m12_general5_bw_sieve 100
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXP 64
#define MAXPR 32
#define MAXDEG 12

typedef struct {
    int p, nonsquare;
    unsigned char *ok;       /* p^2 entries indexed by b*p+w */
    unsigned char *goodres;  /* distinguishes a genuine mask from bad reduction */
    int good, allowed, bad;
} Mask;

typedef struct {
    int n, d;
} Rat;

static const int prime_list[] = {
    3,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
    67,71,73,79,83,89,97,101,103,107,109,113
};
static const int nprimes = (int)(sizeof(prime_list)/sizeof(prime_list[0]));

static int imod(int64_t a, int p) {
    int r = (int)(a % p);
    return r < 0 ? r + p : r;
}

static int ipow(int a, int64_t e, int p) {
    int64_t r = 1, x = imod(a,p);
    while (e) {
        if (e & 1) r = r*x % p;
        x = x*x % p;
        e >>= 1;
    }
    return (int)r;
}

static int invmod(int a, int p) {
    return ipow(a,p-2,p);
}

static int gcd_int(int a, int b) {
    if (a < 0) a = -a;
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

static int legendre(int a, int p) {
    a = imod(a,p);
    if (!a) return 0;
    return ipow(a,(p-1)/2,p)==1 ? 1 : -1;
}

static int first_nonsquare(int p) {
    for (int a=2;a<p;a++) if (legendre(a,p)<0) return a;
    return -1;
}

/* F_{p^2}=F_p[u]/(u^2-ns), encoded as a+p*b. */
static int f2_add(int x, int y, int p) {
    int a=(x%p+y%p)%p, b=(x/p+y/p)%p;
    return a+p*b;
}

static int f2_mul(int x, int y, int p, int ns) {
    int a=x%p,b=x/p,c=y%p,d=y/p;
    int aa=imod((int64_t)a*c+(int64_t)b*d*ns,p);
    int bb=imod((int64_t)a*d+(int64_t)b*c,p);
    return aa+p*bb;
}

static int f2_pow(int x, int64_t e, int p, int ns) {
    int r=1, y=x;
    while (e) {
        if (e&1) r=f2_mul(r,y,p,ns);
        y=f2_mul(y,y,p,ns);
        e>>=1;
    }
    return r;
}

static int f2_char(int x, int p, int ns) {
    if (!x) return 0;
    return f2_pow(x,((int64_t)p*p-1)/2,p,ns)==1 ? 1 : -1;
}

static void poly_zero(int *a) {
    for (int i=0;i<MAXDEG;i++) a[i]=0;
}

static void poly_add_to(int *a, const int *b, int p) {
    for (int i=0;i<MAXDEG;i++) a[i]=imod((int64_t)a[i]+b[i],p);
}

static void poly_mul(const int *a, int da, const int *b, int db, int *c, int p) {
    poly_zero(c);
    for (int i=0;i<=da;i++) for (int j=0;j<=db;j++)
        c[i+j]=imod((int64_t)c[i+j]+(int64_t)a[i]*b[j],p);
}

static void poly_scale(int *a, int d, int c, int p) {
    for (int i=0;i<=d;i++) a[i]=imod((int64_t)a[i]*c,p);
}

/* Build the compact M(12) quintic F. */
static void family_poly(int bv, int wv, int p, int *f) {
    int L[MAXDEG], A[MAXDEG], A2[MAXDEG], LA2[MAXDEG];
    int one2[MAXDEG], wLx2[MAXDEG], prod[MAXDEG], inner[MAXDEG];
    poly_zero(L); poly_zero(A); poly_zero(one2); poly_zero(wLx2);
    L[0]=bv; L[1]=imod(2*bv-1,p);
    A[0]=wv; A[1]=imod(1+(int64_t)bv*wv,p);
    one2[0]=1; one2[1]=2%p; one2[2]=1;
    wLx2[0]=imod((int64_t)wv*L[0],p);
    wLx2[1]=imod((int64_t)wv*L[1],p);
    wLx2[2]=imod(-1,p);
    poly_mul(A,1,A,1,A2,p);
    poly_mul(L,1,A2,2,LA2,p);
    poly_mul(one2,2,wLx2,2,prod,p);
    poly_scale(prod,4,imod(4*bv,p),p);
    memcpy(inner,LA2,sizeof(inner));
    poly_add_to(inner,prod,p);
    poly_mul(L,1,inner,4,f,p);
}

static int poly_degree(const int *a, int maxd, int p) {
    for (int d=maxd;d>=0;d--) if (imod(a[d],p)) return d;
    return -1;
}

static void poly_rem(int *a, int *da, const int *b, int db, int p) {
    int ib=invmod(b[db],p);
    while (*da>=db && *da>=0) {
        int c=imod((int64_t)a[*da]*ib,p), sh=*da-db;
        for (int j=0;j<=db;j++)
            a[j+sh]=imod((int64_t)a[j+sh]-(int64_t)c*b[j],p);
        *da=poly_degree(a,*da,p);
    }
}

static int squarefree_quintic(const int *f, int p) {
    if (poly_degree(f,5,p)!=5) return 0;
    int a[MAXDEG], b[MAXDEG], r[MAXDEG];
    memcpy(a,f,sizeof(a));
    poly_zero(b);
    for (int i=1;i<=5;i++) b[i-1]=imod((int64_t)i*f[i],p);
    int da=5, db=poly_degree(b,4,p);
    while (db>0) {
        memcpy(r,a,sizeof(r));
        int dr=da;
        poly_rem(r,&dr,b,db,p);
        memcpy(a,b,sizeof(a)); da=db;
        memcpy(b,r,sizeof(b)); db=dr;
    }
    return db==0;
}

static int eval_fp(const int *f, int x, int p) {
    int r=0;
    for (int i=5;i>=0;i--) r=imod((int64_t)r*x+f[i],p);
    return r;
}

static int eval_f2(const int *f, int x, int p, int ns) {
    int r=0;
    for (int i=5;i>=0;i--) r=f2_add(f2_mul(r,x,p,ns),f[i],p);
    return r;
}

static int jac_order(const int *f, int p, int ns) {
    int n1=1;
    for (int x=0;x<p;x++) n1 += 1+legendre(eval_fp(f,x,p),p);
    int n2=1;
    for (int x=0;x<p*p;x++) n2 += 1+f2_char(eval_f2(f,x,p,ns),p,ns);
    int c1=n1-p-1;
    int c2=(n2-p*p-1+c1*c1)/2;
    return 1+c1+c2+p*c1+p*p;
}

static void build_mask(Mask *m, int p) {
    m->p=p; m->nonsquare=first_nonsquare(p);
    m->ok=(unsigned char*)calloc((size_t)p*p,1);
    m->goodres=(unsigned char*)calloc((size_t)p*p,1);
    m->good=m->allowed=m->bad=0;
    if (!m->ok || !m->goodres || m->nonsquare<0) { fprintf(stderr,"mask allocation failure\n"); exit(2); }
    for (int b=0;b<p;b++) for (int w=0;w<p;w++) {
        int idx=b*p+w, f[MAXDEG];
        /* Explicit chart boundaries and every singular reduction pass. */
        if (!b || b==1 || imod(2*b-1,p)==0 || !w) {
            m->ok[idx]=1; m->bad++; continue;
        }
        family_poly(b,w,p,f);
        if (!squarefree_quintic(f,p)) {
            m->ok[idx]=1; m->bad++; continue;
        }
        m->good++;
        m->goodres[idx]=1;
        int nj=jac_order(f,p,m->nonsquare);
        /* Positive control: the marked order-12 point injects at p != 2,3. */
        if (p!=3 && nj%12!=0) {
            fprintf(stderr,"ORDER12 CONTROL FAILED p=%d b=%d w=%d #J=%d\n",p,b,w,nj);
            exit(3);
        }
        if (nj%5==0) {
            m->ok[idx]=1; m->allowed++;
        }
    }
    printf("MASK p=%d good=%d allowed5=%d bad=%d killed=%d\n",
           p,m->good,m->allowed,m->bad,m->good-m->allowed);
    fflush(stdout);
}

static Rat *rational_parameters(int H, int *count) {
    size_t cap=(size_t)2*H*H+16, n=0;
    Rat *v=(Rat*)malloc(cap*sizeof(Rat));
    if (!v) { fprintf(stderr,"parameter allocation failure\n"); exit(2); }
    for (int d=1;d<=H;d++) for (int a=-H;a<=H;a++) {
        if (gcd_int(a,d)!=1) continue;
        if (n==cap) { cap*=2; v=(Rat*)realloc(v,cap*sizeof(Rat)); if(!v) exit(2); }
        v[n++]=(Rat){a,d};
    }
    *count=(int)n;
    return v;
}

static int exact_b_boundary(Rat q) {
    return q.n==0 || q.n==q.d || 2*q.n==q.d;
}

static int exact_linear_discriminant(Rat b, Rat w) {
    /* w*(b-1)+1=0, the visible multiplicity-six discriminant component. */
    return (int64_t)w.n*(b.n-b.d)+(int64_t)w.d*b.d==0;
}

static int passes(const Mask *masks, int nm, Rat b, Rat w, int *good_checked) {
    *good_checked=0;
    for (int i=0;i<nm;i++) {
        int p=masks[i].p;
        int bd=imod(b.d,p), wd=imod(w.d,p);
        if (!bd || !wd) continue;
        int bv=imod((int64_t)imod(b.n,p)*invmod(bd,p),p);
        int wv=imod((int64_t)imod(w.n,p)*invmod(wd,p),p);
        int idx=bv*p+wv;
        /* ok includes bad residues; only count genuinely good constraints. */
        if (!masks[i].ok[idx]) return 0;
        if (masks[i].goodres[idx]) (*good_checked)++;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc<2) { fprintf(stderr,"usage: %s HEIGHT [NPRIMES] [MIN_GOOD_TO_PRINT]\n",argv[0]); return 2; }
    int H=atoi(argv[1]);
    int nm=argc>=3 ? atoi(argv[2]) : nprimes;
    int min_good=argc>=4 ? atoi(argv[3]) : 0;
    if (H<1 || nm<1 || nm>nprimes) return 2;
    Mask masks[MAXPR];
    for (int i=0;i<nm;i++) build_mask(&masks[i],prime_list[i]);

    int nr=0; Rat *par=rational_parameters(H,&nr);
    printf("SCAN H=%d parameters=%d pairs=%llu primes=%d\n",
           H,nr,(unsigned long long)nr*(unsigned long long)nr,nm);
    uint64_t checked=0, survivors=0, printed=0;
    for (int ib=0;ib<nr;ib++) {
        Rat b=par[ib];
        if (exact_b_boundary(b)) continue;
        for (int iw=0;iw<nr;iw++) {
            Rat w=par[iw];
            if (w.n==0) continue;
            if (exact_linear_discriminant(b,w)) continue;
            checked++;
            int ngood=0;
            if (!passes(masks,nm,b,w,&ngood)) continue;
            survivors++;
            if (ngood>=min_good) {
                printed++;
                printf("SURVIVOR b=%d/%d w=%d/%d good_masks=%d\n",
                       b.n,b.d,w.n,w.d,ngood);
            }
        }
        if ((ib&1023)==0) {
            fprintf(stderr,"PROGRESS b_index=%d/%d checked=%llu survivors=%llu\n",
                    ib,nr,(unsigned long long)checked,(unsigned long long)survivors);
        }
    }
    printf("DONE H=%d checked=%llu survivors=%llu printed_min_good_%d=%llu\n",H,
           (unsigned long long)checked,(unsigned long long)survivors,min_good,
           (unsigned long long)printed);
    for (int i=0;i<nm;i++) { free(masks[i].ok); free(masks[i].goodres); }
    free(par);
    return 0;
}
