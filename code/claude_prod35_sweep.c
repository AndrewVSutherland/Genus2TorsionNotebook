/* PRODUCTION C35 sweep on the Elkies A_1(5) chart (target Z/35).
 * Model: L=x, L'=1, Q=q2 x^2+q1 x+q0, Q'=Q-x, H=Q'-x*Q, curve z^2 = f = H^2 + 4 Q^2 Q'.
 * Universal 5-torsion; sieve kills smooth chart reductions with 7 !| #J(F_p)
 * (prime-to-p part of 35 at p=5,7).  Boundary/nonintegral = conservative pass.
 *
 * Production upgrades over code/claude_c35_sweep.c:
 *   - kill tables through p<=43 (optionally p<=61 with nfp=17), checked big-prime-first
 *   - degenerate planes q0=0, q2=0 excluded in-engine (globally singular models)
 *   - quotient by the involution (q0,q1,q2)->(q2,1-q1,q0): a whole q1-slice is skipped
 *     when 1-q1 is in the height box with smaller index; on the fixed slice q1=1/2
 *     only i0<=i2 is enumerated
 *   - chunked by kept q1-slice round-robin: ./prod35 search H nchunks chunkid [nfp]
 *
 * Usage:
 *   ./prod35 search Hmax nchunks chunkid [nfp(13|17)]
 *   ./prod35 tabcheck p q0 q1 q2          -- print table verdict + #J for one residue
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef int32_t i32;

static i32 P;
static i32 nres;
static int8_t *qr = 0;
static int8_t *qr2 = 0;

static i32 addm(i32 a, i32 b){ i32 s=a+b; return s>=P? s-P : s; }
static i32 subm(i32 a, i32 b){ i32 s=a-b; return s<0? s+P : s; }
static i32 mulm(i32 a, i32 b){ return (i32)(((i64)a*b)%P); }
static i32 powm(i32 a, i64 e){ i32 r=1; a%=P; if(a<0)a+=P; while(e){ if(e&1) r=mulm(r,a); a=mulm(a,a); e>>=1;} return r; }

static void setup_prime(i32 p){
  P=p;
  free(qr); free(qr2);
  qr = malloc((size_t)p); qr2 = malloc((size_t)p*p);
  memset(qr, -1, (size_t)p); qr[0]=0;
  for(i32 v=1;v<p;v++) qr[mulm(v,v)] = 1;
  nres=2; while(qr[nres]!=-1) nres++;
  memset(qr2, -1, (size_t)p*p);
  qr2[0]=0;
  for(i32 a=0;a<p;a++) for(i32 b=0;b<p;b++){
    if(a==0&&b==0) continue;
    i32 re = addm(mulm(a,a), mulm(nres, mulm(b,b)));
    i32 im = mulm(2, mulm(a,b));
    qr2[(size_t)re*p+im] = 1;
  }
}

static int poly_gcd_deg(i32 *Ain, int da, i32 *Bin, int db){
  i32 A[9], B[9], T[9];
  int i;
  for(i=0;i<=da;i++) A[i]=Ain[i];
  for(i=0;i<=db;i++) B[i]=Bin[i];
  while(db>=0){
    i32 inv = powm(B[db], P-2);
    int k;
    for(k=da;k>=db;k--){
      i32 c = mulm(A[k], inv);
      if(c){ int j; for(j=0;j<db;j++) A[k-db+j] = subm(A[k-db+j], mulm(c,B[j])); }
      A[k]=0;
    }
    int nd = db-1; while(nd>=0 && A[nd]==0) nd--;
    for(i=0;i<=db;i++) T[i]=B[i];
    for(i=0;i<=nd && nd>=0;i++) B[i]=A[i];
    for(i=0;i<=db;i++) A[i]=T[i];
    da = db; db = nd;
  }
  return da;
}

static void build_f(i32 q0, i32 q1, i32 q2, i32 *f){
  i32 Q[3]  = {q0, q1, q2};
  i32 Qp[3] = {q0, subm(q1,1), q2};
  i32 H[4];
  H[0]=q0; H[1]=subm(subm(q1,1),q0); H[2]=subm(q2,q1); H[3]=subm(0,q2);
  i32 Q2[5]; int i,j;
  for(i=0;i<5;i++) Q2[i]=0;
  for(i=0;i<3;i++) for(j=0;j<3;j++) Q2[i+j]=addm(Q2[i+j], mulm(Q[i],Q[j]));
  i32 QQ[7]; for(i=0;i<7;i++) QQ[i]=0;
  for(i=0;i<5;i++) for(j=0;j<3;j++) QQ[i+j]=addm(QQ[i+j], mulm(Q2[i],Qp[j]));
  i32 H2[7]; for(i=0;i<7;i++) H2[i]=0;
  for(i=0;i<4;i++) for(j=0;j<4;j++) H2[i+j]=addm(H2[i+j], mulm(H[i],H[j]));
  for(i=0;i<7;i++) f[i]=addm(H2[i], mulm(4%P, QQ[i]));
}

static int count_curve(i32 *f, i64 *Jout){
  int d = 6;
  if(f[6]==0){ if(f[5]==0) return 0; d=5; }
  i32 fp[8]; int i;
  int dd = -1;
  for(i=0;i<d;i++){ fp[i] = mulm((i+1)%P, f[i+1]); if(fp[i]) dd=i; }
  if(dd<0) return 0;
  if(poly_gcd_deg(f, d, fp, dd) > 0) return 0;
  i64 N1 = (d==6) ? (qr[f[6]]==1 ? 2 : 0) : 1;
  for(i32 x=0;x<P;x++){
    i32 v=f[d]; int k;
    for(k=d-1;k>=0;k--) v = addm(mulm(v,x), f[k]);
    N1 += 1 + qr[v];
  }
  i64 N2 = 2;
  if(d==5) N2 = 1;
  for(i32 a=0;a<P;a++) for(i32 b=0;b<P;b++){
    i32 vr=f[d], vi=0; int k;
    for(k=d-1;k>=0;k--){
      i32 nr = addm(mulm(vr,a), mulm(nres, mulm(vi,b)));
      i32 ni = addm(mulm(vr,b), mulm(vi,a));
      vr = addm(nr, f[k]); vi = ni;
    }
    N2 += 1 + ((vr==0&&vi==0)?0:qr2[(size_t)vr*P+vi]);
  }
  i64 c1 = N1 - (P+1);
  i64 c2 = (N2 - ((i64)P*P+1) + c1*c1)/2;
  *Jout = 1 + c1 + c2 + (i64)P*c1 + (i64)P*P;
  return 1;
}

/* ---------------- production search ---------------- */
#define MAXFP 17
/* big-prime-first check order for kill efficiency; table build order irrelevant */
static i32 FPALL[MAXFP] = {41,43,37,31,29,23,19,17,13,11,7,5,3,47,53,59,61};
static int NFPuse = 13;
static int8_t *tab[MAXFP];  /* 0 boundary-pass, 1 kill, 2 good-pass */

static void build_tables(void){
  int t;
  for(t=0;t<NFPuse;t++){
    i32 p = FPALL[t];
    setup_prime(p);
    i64 req = 35; if(p==5) req=7; if(p==7) req=5;
    tab[t] = malloc((size_t)p*p*p);
    for(i32 q2=0;q2<p;q2++) for(i32 q1=0;q1<p;q1++) for(i32 q0=0;q0<p;q0++){
      i32 f[7]; i64 J;
      build_f(q0,q1,q2,f);
      size_t idx = ((size_t)q2*p+q1)*p+q0;
      if(!count_curve(f,&J)){ tab[t][idx]=0; continue; }
      tab[t][idx] = (J%req==0) ? 2 : 1;
    }
    i64 k0=0,k1=0,k2=0;
    for(i64 i=0;i<(i64)p*p*p;i++){ if(tab[t][i]==0)k0++; else if(tab[t][i]==1)k1++; else k2++; }
    fprintf(stderr,"table p=%d boundary=%lld kill=%lld pass=%lld\n",p,(long long)k0,(long long)k1,(long long)k2);
  }
}

#define MAXR 12000
static i32 nums[MAXR], dens[MAXR];
static int R=0, izero=-1, ihalf=-1;
static i32 invq[MAXR];              /* index of 1-q, or -1 if outside box */
static int8_t res8[MAXFP][MAXR];    /* residue of rational i mod FPALL[t], -1 if p|den */

static void build_rationals(i32 H){
  for(i32 d=1;d<=H;d++) for(i32 n=-H;n<=H;n++){
    i32 a = n<0? -n:n, b=d, g;
    while(b){ g=a%b; a=b; b=g; }
    if(n==0){ if(d!=1) continue; }
    else if(a!=1) continue;
    if(R>=MAXR){ fprintf(stderr,"MAXR overflow\n"); exit(1); }
    nums[R]=n; dens[R]=d; R++;
  }
  /* index lookup: key = (d-1)*(2H+1) + (n+H) */
  i32 *lk = malloc(sizeof(i32)*(size_t)H*(2*H+1));
  for(i64 i=0;i<(i64)H*(2*H+1);i++) lk[i]=-1;
  for(int i=0;i<R;i++) lk[(i64)(dens[i]-1)*(2*H+1) + (nums[i]+H)] = i;
  for(int i=0;i<R;i++){
    /* 1 - n/d = (d-n)/d, then reduce */
    i32 n2 = dens[i]-nums[i], d2 = dens[i];
    i32 a = n2<0? -n2:n2, b=d2, g;
    while(b){ g=a%b; a=b; b=g; }
    if(a){ n2/=a; d2/=a; }
    else { d2=1; }                   /* n2==0 */
    invq[i] = (n2>=-((i32)H) && n2<=(i32)H && d2<=(i32)H) ? lk[(i64)(d2-1)*(2*H+1)+(n2+H)] : -1;
    if(nums[i]==0 && dens[i]==1) izero=i;
    if(nums[i]==1 && dens[i]==2) ihalf=i;
  }
  free(lk);
  for(int t=0;t<NFPuse;t++){
    i32 p=FPALL[t]; P=p;
    for(int i=0;i<R;i++){
      if(dens[i]%p==0){ res8[t][i]=-1; continue; }
      i32 nn = nums[i]%p; if(nn<0)nn+=p;
      res8[t][i] = (int8_t)mulm(nn, powm(dens[i]%p, p-2));
    }
  }
  fprintf(stderr,"rationals height<=%d: %d (izero=%d ihalf=%d)\n",H,R,izero,ihalf);
}

static void search(i32 H, int nchunks, int chunkid){
  build_tables();
  build_rationals(H);
  int nslices=0, nskip=0;
  i64 tot=0, surv=0;
  int actT[MAXFP]; i32 actP[MAXFP];
  int8_t *rowptr[MAXFP];
  int kept=0;
  for(int i1=0;i1<R;i1++){
    if(invq[i1]>=0 && invq[i1]<i1){ nskip++; continue; }   /* involution partner slice earlier */
    int sid = kept++;
    if(sid % nchunks != chunkid) continue;
    nslices++;
    int diag = (i1==ihalf);          /* fixed slice q1=1/2: enumerate i0<=i2 only */
    /* active primes for this slice (q1 integral at p) */
    int nact1=0; int act1[MAXFP];
    for(int t=0;t<NFPuse;t++) if(res8[t][i1]>=0) act1[nact1++]=t;
    for(int i2=0;i2<R;i2++){
      if(i2==izero) continue;
      int nact=0;
      for(int u=0;u<nact1;u++){
        int t=act1[u];
        int r2 = res8[t][i2];
        if(r2<0) continue;
        i32 p=FPALL[t];
        rowptr[nact] = tab[t] + ((size_t)r2*p + (i32)res8[t][i1])*p;
        actT[nact] = t; actP[nact]=p; (void)actP;
        nact++;
      }
      int i0max = diag ? i2 : R-1;
      for(int i0=0;i0<=i0max;i0++){
        if(i0==izero) continue;
        tot++;
        int ok=1;
        for(int u=0;u<nact;u++){
          int r0 = res8[actT[u]][i0];
          if(r0<0) continue;
          if(rowptr[u][r0]==1){ ok=0; break; }
        }
        if(ok){
          surv++;
          printf("SURV %d/%d %d/%d %d/%d\n", nums[i0],dens[i0], nums[i1],dens[i1], nums[i2],dens[i2]);
        }
      }
    }
    fprintf(stderr,"slice i1=%d (%d/%s) done tot=%lld surv=%lld\n", i1, nslices, "chunk", (long long)tot,(long long)surv);
    fflush(stderr); fflush(stdout);
  }
  fprintf(stderr,"DONE chunk %d/%d H=%d slices=%d skipped_slices=%d triples=%lld survivors=%lld\n",
    chunkid,nchunks,H,nslices,nskip,(long long)tot,(long long)surv);
}

int main(int argc, char **argv){
  if(argc>=5 && !strcmp(argv[1],"search")){
    i32 H = atoi(argv[2]);
    int nc = atoi(argv[3]), ci = atoi(argv[4]);
    if(argc>=6) NFPuse = atoi(argv[5]);
    if(NFPuse<1||NFPuse>MAXFP) NFPuse=13;
    search(H,nc,ci);
    return 0;
  }
  if(argc>=6 && !strcmp(argv[1],"tabcheck")){
    i32 p=atoi(argv[2]);
    setup_prime(p);
    i32 f[7]; i64 J;
    i32 r3=(i32)(((long)atoi(argv[3])%p+p)%p), r4=(i32)(((long)atoi(argv[4])%p+p)%p), r5=(i32)(((long)atoi(argv[5])%p+p)%p);
    build_f(r3, r4, r5, f);
    int st = count_curve(f,&J);
    printf("p=%d status=%s J=%lld Jmod35=%lld\n", p, st?"smooth":"boundary", st?(long long)J:0, st?(long long)(J%35):-1);
    return 0;
  }
  fprintf(stderr,"usage: %s search Hmax nchunks chunkid [nfp] | tabcheck p q0 q1 q2\n",argv[0]);
  return 1;
}
