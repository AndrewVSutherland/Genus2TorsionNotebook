/* C35 target: Elkies A_1(5) chart local densities + CRT-guided rational search.
 * Model: L=x, L'=1, Q=q2 x^2+q1 x+q0, Q'=Q-x, H=Q'-x*Q, curve z^2 = f = H^2 + 4 Q^2 Q'.
 * Usage:
 *   ./c35_sweep density p [nsample]     -- chart stats over F_p (full or sampled)
 *   ./c35_sweep search Hmax             -- rational triples height<=Hmax, filtered at
 *                                          p in {3,5,7,11,13,17,19,23,29}
 * Single-threaded.
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

/* gcd degree of two polys (arrays low->high), returns degree of gcd */
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
    /* swap: new A = old B, new B = remainder */
    for(i=0;i<=db;i++) T[i]=B[i];
    for(i=0;i<=nd && nd>=0;i++) B[i]=A[i];
    for(i=0;i<=db;i++) A[i]=T[i];
    da = db; db = nd;
  }
  return da;
}

/* build f = H^2 + 4 Q^2 Q' from (q0,q1,q2) mod P; f has 7 slots */
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

/* status: 0 = boundary (deg<5 or not squarefree), 1 = good, *Jout set */
static int count_curve(i32 *f, i64 *Jout){
  int d = 6;
  if(f[6]==0){ if(f[5]==0) return 0; d=5; }
  i32 fp[8]; int i;
  int dd = -1;
  for(i=0;i<d;i++){ fp[i] = mulm((i+1)%P, f[i+1]); if(fp[i]) dd=i; }
  if(dd<0) return 0;             /* f' = 0 mod p => not squarefree */
  if(poly_gcd_deg(f, d, fp, dd) > 0) return 0;
  i64 N1 = (d==6) ? (qr[f[6]]==1 ? 2 : 0) : 1;
  for(i32 x=0;x<P;x++){
    i32 v=f[d]; int k;
    for(k=d-1;k>=0;k--) v = addm(mulm(v,x), f[k]);
    N1 += 1 + qr[v];
  }
  i64 N2 = 2;                    /* deg6: lc always square in F_p^2; deg5: 1 */
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

static void density(i32 p, i64 nsample){
  setup_prime(p);
  i64 tot=0, bad=0, sm=0, n5=0, n7=0, n35=0;
  i64 total_space = (i64)p*p*p;
  uint64_t rng = 88172645463325252ULL;
  i64 iter = nsample>0 ? nsample : total_space;
  for(i64 it=0; it<iter; it++){
    i32 q0,q1,q2;
    if(nsample>0){
      rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17;
      i64 v = (i64)(rng % (uint64_t)total_space);
      q0 = v%p; q1=(v/p)%p; q2=(v/((i64)p*p))%p;
    } else {
      q0 = it%p; q1=(it/p)%p; q2=(it/((i64)p*p))%p;
    }
    tot++;
    i32 f[7]; i64 J;
    build_f(q0,q1,q2,f);
    if(!count_curve(f,&J)){ bad++; continue; }
    sm++;
    if(J%5==0) n5++;
    if(J%7==0) n7++;
    if(J%35==0) n35++;
  }
  printf("density p=%d %s tot=%lld boundary=%lld smooth=%lld 5|J=%lld 7|J=%lld 35|J=%lld frac35=%.5f\n",
    p, nsample>0?"(sampled)":"(full)", (long long)tot,(long long)bad,(long long)sm,
    (long long)n5,(long long)n7,(long long)n35, sm? (double)n35/sm : 0.0);
}

/* ---------------- search mode ---------------- */
#define NFP 9
static i32 FP[NFP] = {3,5,7,11,13,17,19,23,29};
static int8_t *tab[NFP];  /* 0 boundary-pass, 1 kill, 2 good-pass */

static void build_tables(void){
  int t;
  for(t=0;t<NFP;t++){
    i32 p = FP[t];
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

static void search(i32 H){
  build_tables();
  /* enumerate rationals num/den, |num|<=H, 1<=den<=H, gcd=1 */
  i32 nums[400], dens[400]; int R=0;
  for(i32 d=1;d<=H;d++) for(i32 n=-H;n<=H;n++){
    i32 a = n<0? -n:n, b=d, g;
    while(b){ g=a%b; a=b; b=g; }
    if(a!=1 && !(n==0&&d==1)) continue;
    if(n==0&&d!=1) continue;
    nums[R]=n; dens[R]=d; R++;
  }
  fprintf(stderr,"rationals height<=%d: %d\n",H,R);
  /* residue of each rational at each filter prime; -1 if den divisible by p */
  static i32 res[NFP][400];
  for(int t=0;t<NFP;t++){
    i32 p=FP[t]; P=p;
    for(int i=0;i<R;i++){
      if(dens[i]%p==0){ res[t][i]=-1; continue; }
      i32 nn = nums[i]%p; if(nn<0)nn+=p;
      res[t][i] = mulm(nn, powm(dens[i]%p, p-2));
    }
  }
  i64 tot=0, surv=0;
  for(int i0=0;i0<R;i0++) for(int i1=0;i1<R;i1++) for(int i2=0;i2<R;i2++){
    tot++;
    int ok=1;
    for(int t=0;t<NFP;t++){
      i32 p=FP[t];
      i32 r0=res[t][i0], r1=res[t][i1], r2=res[t][i2];
      if(r0<0||r1<0||r2<0) continue;              /* nonintegral: conservative pass */
      size_t idx = ((size_t)r2*p+r1)*p+r0;
      if(tab[t][idx]==1){ ok=0; break; }
    }
    if(ok){
      surv++;
      printf("SURV %d/%d %d/%d %d/%d\n", nums[i0],dens[i0], nums[i1],dens[i1], nums[i2],dens[i2]);
    }
  }
  fprintf(stderr,"search H=%d triples=%lld survivors=%lld\n",H,(long long)tot,(long long)surv);
}

int main(int argc, char **argv){
  if(argc>=3 && !strcmp(argv[1],"density")){
    i32 p = atoi(argv[2]);
    i64 ns = argc>=4? atoll(argv[3]) : 0;
    density(p, ns);
    return 0;
  }
  if(argc>=3 && !strcmp(argv[1],"search")){
    search(atoi(argv[2]));
    return 0;
  }
  fprintf(stderr,"usage: %s density p [nsample] | search Hmax\n",argv[0]);
  return 1;
}
