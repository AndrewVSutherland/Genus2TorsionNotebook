/*
 * Memory-bounded modular sieve for cyclic order 60 in the split-B Hermite
 * family with a marked typical 5-torsion class.
 *
 * Put
 *
 *   B=x(x-1),  q=x^2+(t^2-s^2-1)x+s^2,
 *
 * and let A be the monic quintic satisfying
 *
 *   A(0)=s^5, A'(0)=5*s^3*q'(0)/2,
 *   A(1)=t^5, A'(1)=5*t^3*q'(1)/2.
 *
 * Its x^4 coefficient m is free.  The quotient
 *
 *   F=(A^2-q^5)/(x^2(x-1)^2)
 *
 * is a quintic, and [q,A mod q] has order dividing 5 on y^2=F.  On the
 * squarefree/coprime open it is the marked nonzero 5-class.  A rational
 * cyclic order-60 subgroup therefore forces 60 | #J(F_p) at every prime of
 * good reduction p != 2,3,5.
 *
 * Bad reductions and reductions leaving the marked open are passed
 * conservatively.  Masks are computed from first principles using C(F_p)
 * and C(F_{p^2}).  The integer scanner combines the first SEED masks by CRT
 * in the m-coordinate, so its main cost is quadratic rather than cubic in
 * the height once the local conditions become selective.
 *
 * Build and use, for example:
 *
 *   clang -O3 -Wall -Wextra -o /tmp/hermite60 \
 *       code/hermite5_order60_modular_sieve.c
 *   /tmp/hermite60 integer 1000 8 3
 *   /tmp/hermite60 rational 20 8
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXDEG 12
#define MAXPR 16

typedef struct {
    int p, nonsquare;
    unsigned char *ok;       /* p^3 entries, s-major then t then m */
    unsigned char *goodres;  /* a genuine good/open residue */
    int good, allowed, bad;
} Mask;

typedef struct { int n, d; } Rat;

static const int prime_list[] = { 7,11,13,17,19,23,29,31,37,41,43 };
static const int nprimes=(int)(sizeof(prime_list)/sizeof(prime_list[0]));

static int imod(int64_t a,int p) { int r=(int)(a%p); return r<0?r+p:r; }
static int ipow(int a,int64_t e,int p) {
    int64_t r=1,x=imod(a,p);
    while(e) { if(e&1) r=r*x%p; x=x*x%p; e>>=1; }
    return (int)r;
}
static int invmod(int a,int p) { return ipow(a,p-2,p); }
static int gcd_int(int a,int b) {
    if(a<0)a=-a;
    if(b<0)b=-b;
    while(b){int z=a%b;a=b;b=z;} return a;
}
static int legendre(int a,int p) {
    a=imod(a,p); if(!a)return 0;
    return ipow(a,(p-1)/2,p)==1?1:-1;
}
static int first_nonsquare(int p) {
    for(int a=2;a<p;a++) {
        if(legendre(a,p)<0)return a;
    }
    return -1;
}

/* F_(p^2)=F_p[u]/(u^2-ns), encoded as a+p*b. */
static int f2_add(int x,int y,int p) {
    return (x%p+y%p)%p+p*((x/p+y/p)%p);
}
static int f2_mul(int x,int y,int p,int ns) {
    int a=x%p,b=x/p,c=y%p,d=y/p;
    return imod((int64_t)a*c+(int64_t)b*d*ns,p)
        +p*imod((int64_t)a*d+(int64_t)b*c,p);
}
static int f2_pow(int x,int64_t e,int p,int ns) {
    int r=1; while(e){if(e&1)r=f2_mul(r,x,p,ns);x=f2_mul(x,x,p,ns);e>>=1;}return r;
}
static int f2_char(int x,int p,int ns) {
    if(!x)return 0;
    return f2_pow(x,((int64_t)p*p-1)/2,p,ns)==1?1:-1;
}

static void pzero(int *a){for(int i=0;i<MAXDEG;i++)a[i]=0;}
static int pdegree(const int *a,int d,int p){for(;d>=0;d--)if(imod(a[d],p))return d;return -1;}
static void pmul(const int *a,int da,const int *b,int db,int *c,int p){
    pzero(c); for(int i=0;i<=da;i++)for(int j=0;j<=db;j++)
        c[i+j]=imod((int64_t)c[i+j]+(int64_t)a[i]*b[j],p);
}
static void prem(int *a,int *da,const int *b,int db,int p){
    int ib=invmod(b[db],p);
    while(*da>=db&&*da>=0){
        int c=imod((int64_t)a[*da]*ib,p),sh=*da-db;
        for(int j=0;j<=db;j++)a[j+sh]=imod((int64_t)a[j+sh]-(int64_t)c*b[j],p);
        *da=pdegree(a,*da,p);
    }
}
static int pgcd_degree(const int *aa,int da,const int *bb,int db,int p){
    int a[MAXDEG],b[MAXDEG],r[MAXDEG]; memcpy(a,aa,sizeof(a));memcpy(b,bb,sizeof(b));
    da=pdegree(a,da,p);db=pdegree(b,db,p);
    while(db>=0){memcpy(r,a,sizeof(r));int dr=da;prem(r,&dr,b,db,p);
        memcpy(a,b,sizeof(a));da=db;memcpy(b,r,sizeof(b));db=dr;}
    return da;
}
static int squarefree_quintic(const int *f,int p){
    if(pdegree(f,5,p)!=5)return 0;
    int df[MAXDEG];pzero(df);for(int i=1;i<=5;i++)df[i-1]=imod((int64_t)i*f[i],p);
    return pgcd_degree(f,5,df,4,p)==0;
}
static int eval_fp(const int *f,int x,int p){int r=0;for(int i=5;i>=0;i--)r=imod((int64_t)r*x+f[i],p);return r;}
static int eval_f2(const int *f,int x,int p,int ns){int r=0;for(int i=5;i>=0;i--)r=f2_add(f2_mul(r,x,p,ns),f[i],p);return r;}
static int jac_order(const int *f,int p,int ns,int *roots1,int *roots2){
    int n1=1,n2=1;*roots1=0;*roots2=0;
    for(int x=0;x<p;x++){
        int z=eval_fp(f,x,p);if(!z)(*roots1)++;
        n1+=1+legendre(z,p);
    }
    for(int x=0;x<p*p;x++){
        int z=eval_f2(f,x,p,ns);if(!z)(*roots2)++;
        n2+=1+f2_char(z,p,ns);
    }
    int c1=n1-p-1,c2=(n2-p*p-1+c1*c1)/2;
    return 1+c1+c2+p*c1+p*p;
}

/* Build F by exact polynomial arithmetic modulo p. */
static void family_poly(int s,int t,int m,int p,int *f,int *qout){
    int inv2=invmod(2,p),s2=imod((int64_t)s*s,p),t2=imod((int64_t)t*t,p);
    int s3=imod((int64_t)s2*s,p),t3=imod((int64_t)t2*t,p);
    int s5=imod((int64_t)s3*s2,p),t5=imod((int64_t)t3*t2,p);
    int e=imod((int64_t)t2-s2-1,p);
    int q[MAXDEG],q2[MAXDEG],q4[MAXDEG],q5[MAXDEG];pzero(q);
    q[0]=s2;q[1]=e;q[2]=1;
    int A[MAXDEG],A2[MAXDEG],N[MAXDEG];pzero(A);
    A[0]=s5;
    A[1]=imod((int64_t)5*inv2%p*s3%p*e,p);
    A[2]=imod((int64_t)m+2LL*s5-5LL*s3*t2+5LL*s3
        +(int64_t)5*inv2%p*s2%p*t3+(int64_t)inv2*t5
        -(int64_t)5*inv2%p*t3+2,p);
    A[3]=imod(-2LL*m-(int64_t)inv2*s5+(int64_t)5*inv2%p*s3%p*t2
        -(int64_t)5*inv2%p*s3-(int64_t)5*inv2%p*s2%p*t3
        +(int64_t)inv2*t5+(int64_t)5*inv2%p*t3-3,p);
    A[4]=imod(m,p);A[5]=1;
    pmul(A,5,A,5,A2,p);pmul(q,2,q,2,q2,p);pmul(q2,4,q2,4,q4,p);pmul(q4,8,q,2,q5,p);
    for(int i=0;i<MAXDEG;i++)N[i]=imod((int64_t)A2[i]-q5[i],p);
    /* Divide by x^2(x-1)^2=x^4-2x^3+x^2. */
    int den[MAXDEG];pzero(den);den[2]=1;den[3]=imod(-2,p);den[4]=1;
    pzero(f);int dn=pdegree(N,10,p);
    while(dn>=4){int c=N[dn],sh=dn-4;f[sh]=c;
        for(int j=2;j<=4;j++)N[j+sh]=imod((int64_t)N[j+sh]-(int64_t)c*den[j],p);
        dn=pdegree(N,dn,p);
    }
    if(dn>=0){fprintf(stderr,"HERMITE DIVISION FAILURE p=%d s=%d t=%d m=%d remdeg=%d\n",p,s,t,m,dn);exit(3);}
    if(qout)memcpy(qout,q,sizeof(q));
}

static size_t midx(int s,int t,int m,int p){return ((size_t)s*p+t)*p+m;}
static void build_mask(Mask *M,int p){
    M->p=p;M->nonsquare=first_nonsquare(p);size_t n=(size_t)p*p*p;
    M->ok=(unsigned char*)calloc(n,1);M->goodres=(unsigned char*)calloc(n,1);
    M->good=M->allowed=M->bad=0;if(!M->ok||!M->goodres||M->nonsquare<0){fprintf(stderr,"allocation failure\n");exit(2);}
    for(int s=0;s<p;s++)for(int t=0;t<p;t++)for(int m=0;m<p;m++){
        size_t z=midx(s,t,m,p);int f[MAXDEG],q[MAXDEG];
        /* s*t=0, repeated q, q meeting F, degree drop, and bad curves pass. */
        if(!s||!t){M->ok[z]=1;M->bad++;continue;}
        family_poly(s,t,m,p,f,q);
        int dq[MAXDEG];pzero(dq);dq[0]=q[1];dq[1]=imod(2,p);
        if(pgcd_degree(q,2,dq,1,p)>0||pgcd_degree(q,2,f,5,p)>0||!squarefree_quintic(f,p)){
            M->ok[z]=1;M->bad++;continue;
        }
        M->good++;M->goodres[z]=1;int roots1,roots2;
        int nj=jac_order(f,p,M->nonsquare,&roots1,&roots2);
        if(nj%5){fprintf(stderr,"ORDER5 CONTROL FAILED p=%d s=%d t=%d m=%d J=%d\n",p,s,t,m,nj);exit(3);}
        /* For a squarefree quintic, dim J(F_p)[2] is the number of
         * irreducible factors minus one.  The 2-Sylow contains an element
         * of order 4 exactly when v_2(#J) is larger than this dimension. */
        int nquad=(roots2-roots1)/2;
        int remaining=5-roots1-2*nquad;
        int factors=roots1+nquad+(remaining>0);
        int r2=factors-1,v2=0,z2=nj;
        while(!(z2&1)){v2++;z2>>=1;}
        if(nj%60==0&&v2>r2){M->ok[z]=1;M->allowed++;}
    }
    printf("MASK p=%d good=%d allowed60=%d bad=%d killed=%d\n",p,M->good,M->allowed,M->bad,M->good-M->allowed);fflush(stdout);
}

static int exact_open_integer(int s,int t,int m){
    (void)m;if(!s||!t)return 0;
    /* Disc(q)=(t^2-s^2-1)^2-4s^2. */
    int64_t e=(int64_t)t*t-(int64_t)s*s-1;
    return e*e-4LL*s*s!=0;
}
static int passes_remaining(const Mask *M,int nm,int start,int s,int t,int m,int *ngood){
    for(int i=start;i<nm;i++){
        int p=M[i].p,sr=imod(s,p),tr=imod(t,p),mr=imod(m,p);size_t z=midx(sr,tr,mr,p);
        if(!M[i].ok[z])return 0;
        if(M[i].goodres[z])(*ngood)++;
    }return 1;
}

static int cmp_int(const void *a,const void *b){int x=*(const int*)a,y=*(const int*)b;return (x>y)-(x<y);}

/* CRT-combine the allowed m residues for fixed integer s,t at seed masks. */
static int crt_m_residues(const Mask *M,int seed,int s,int t,int **out,int *modout){
    int cap=64,n=1,mod=1;int *v=(int*)malloc((size_t)cap*sizeof(int));if(!v)exit(2);v[0]=0;
    for(int i=0;i<seed;i++){
        int p=M[i].p,sr=imod(s,p),tr=imod(t,p),vals[64],nv=0;
        for(int m=0;m<p;m++)if(M[i].ok[midx(sr,tr,m,p)])vals[nv++]=m;
        if(!nv){free(v);*out=NULL;*modout=mod*p;return 0;}
        int nn=n*nv;if(nn>cap){while(cap<nn)cap*=2;v=(int*)realloc(v,(size_t)cap*sizeof(int));if(!v)exit(2);}
        int oldn=n,oldmod=mod,inv=invmod(oldmod%p,p);
        /* Fill backwards is unsafe with expansion; use a temporary array. */
        int *w=(int*)malloc((size_t)nn*sizeof(int));if(!w)exit(2);int k=0;
        for(int j=0;j<oldn;j++)for(int h=0;h<nv;h++){
            int z=imod((int64_t)(vals[h]-v[j])*inv,p);
            w[k++]=v[j]+oldmod*z;
        }
        free(v);v=w;n=nn;mod*=p;cap=nn;
    }
    qsort(v,(size_t)n,sizeof(int),cmp_int);*out=v;*modout=mod;return n;
}

static void scan_integer(const Mask *M,int nm,int H,int seed){
    uint64_t pairs=0,crtclasses=0,lifts=0,survivors=0;int maxclasses=0;
    printf("INTEGER_SCAN H=%d primes=%d seed=%d\n",H,nm,seed);fflush(stdout);
    for(int s=-H;s<=H;s++){
        if(!s)continue;
        for(int t=-H;t<=H;t++){
            if(!t)continue;
            pairs++;
            int *classes=NULL,mod=1,nc=crt_m_residues(M,seed,s,t,&classes,&mod);
            crtclasses+=(uint64_t)nc;if(nc>maxclasses)maxclasses=nc;
            for(int j=0;j<nc;j++){
                int r=classes[j];
                int64_t k0=(-((int64_t)H)-r+mod-1)/mod; /* ceil((-H-r)/mod), corrected below */
                int64_t num=-(int64_t)H-r;
                k0=num>=0?(num+mod-1)/mod:-((-num)/mod);
                for(int64_t k=k0;;k++){
                    int64_t mm=(int64_t)r+k*mod;if(mm>H)break;if(mm< -H)continue;
                    int m=(int)mm;if(!exact_open_integer(s,t,m))continue;lifts++;
                    int ngood=0;
                    for(int i=0;i<seed;i++)if(M[i].goodres[midx(imod(s,M[i].p),imod(t,M[i].p),imod(m,M[i].p),M[i].p)])ngood++;
                    if(!passes_remaining(M,nm,seed,s,t,m,&ngood))continue;
                    survivors++;printf("SURVIVOR s=%d t=%d m=%d good_masks=%d\n",s,t,m,ngood);
                }
            }
            free(classes);
        }
        if(((s+H)&127)==0){fprintf(stderr,"PROGRESS s=%d pairs=%llu lifts=%llu survivors=%llu\n",s,(unsigned long long)pairs,(unsigned long long)lifts,(unsigned long long)survivors);}
    }
    printf("INTEGER_DONE H=%d pairs=%llu crt_classes=%llu max_classes_per_pair=%d candidate_lifts=%llu survivors=%llu\n",
        H,(unsigned long long)pairs,(unsigned long long)crtclasses,maxclasses,(unsigned long long)lifts,(unsigned long long)survivors);
}

static Rat *rational_parameters(int H,int *count){
    size_t cap=(size_t)2*H*H+16,n=0;Rat *v=(Rat*)malloc(cap*sizeof(Rat));if(!v)exit(2);
    for(int d=1;d<=H;d++)for(int a=-H;a<=H;a++)if(gcd_int(a,d)==1){
        if(n==cap){cap*=2;v=(Rat*)realloc(v,cap*sizeof(Rat));if(!v)exit(2);}v[n++]=(Rat){a,d};
    }*count=(int)n;return v;
}
static int rat_res(Rat a,int p,int *ok){int d=imod(a.d,p);if(!d){*ok=0;return 0;}*ok=1;return imod((int64_t)imod(a.n,p)*invmod(d,p),p);}
static int exact_q_disc_rat(Rat s,Rat t){
    /* Clear denominators in (t^2-s^2-1)^2-4s^2. */
    int64_t sd=s.d,td=t.d;
    /* __int128 keeps H up to several thousand safe. */
    __int128 e=(__int128)t.n*t.n*sd*sd-(__int128)s.n*s.n*td*td-(__int128)sd*sd*td*td;
    __int128 rhs=4*(__int128)s.n*s.n*sd*sd*td*td*td*td;
    return e*e!=rhs;
}
static void scan_rational(const Mask *M,int nm,int H){
    int nr=0;Rat *v=rational_parameters(H,&nr);uint64_t checked=0,survivors=0;
    printf("RATIONAL_SCAN H=%d parameters=%d formal_triples=%llu primes=%d\n",H,nr,(unsigned long long)nr*nr*nr,nm);fflush(stdout);
    for(int is=0;is<nr;is++){
        Rat s=v[is];if(!s.n)continue;
        for(int it=0;it<nr;it++){
            Rat t=v[it];if(!t.n||!exact_q_disc_rat(s,t))continue;
            for(int im=0;im<nr;im++){
                Rat m=v[im];checked++;int ngood=0,pass=1;
                for(int i=0;i<nm;i++){
                    int os,ot,om,p=M[i].p,sr=rat_res(s,p,&os),tr=rat_res(t,p,&ot),mr=rat_res(m,p,&om);
                    if(!os||!ot||!om)continue;
                    size_t z=midx(sr,tr,mr,p);
                    if(!M[i].ok[z]){pass=0;break;}if(M[i].goodres[z])ngood++;
                }
                if(pass){survivors++;printf("SURVIVOR s=%d/%d t=%d/%d m=%d/%d good_masks=%d\n",s.n,s.d,t.n,t.d,m.n,m.d,ngood);}
            }
        }
        if((is&31)==0)fprintf(stderr,"PROGRESS s_index=%d/%d checked=%llu survivors=%llu\n",is,nr,(unsigned long long)checked,(unsigned long long)survivors);
    }
    printf("RATIONAL_DONE H=%d checked=%llu survivors=%llu\n",H,(unsigned long long)checked,(unsigned long long)survivors);free(v);
}

int main(int argc,char **argv){
    if(argc<3){fprintf(stderr,"usage: %s integer|rational HEIGHT [NPRIMES] [SEED]\n",argv[0]);return 2;}
    int H=atoi(argv[2]),nm=argc>=4?atoi(argv[3]):8,seed=argc>=5?atoi(argv[4]):3;
    if(H<1||nm<1||nm>nprimes||seed<1||seed>nm||seed>4)return 2;
    Mask M[MAXPR];for(int i=0;i<nm;i++)build_mask(&M[i],prime_list[i]);
    if(!strcmp(argv[1],"integer"))scan_integer(M,nm,H,seed);
    else if(!strcmp(argv[1],"rational"))scan_rational(M,nm,H);
    else{fprintf(stderr,"unknown mode\n");return 2;}
    for(int i=0;i<nm;i++){free(M[i].ok);free(M[i].goodres);}return 0;
}
