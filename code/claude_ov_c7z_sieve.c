// claude_ov_c7z_sieve.c -- Lane 2 contact-7 harvester, z-coordinate version.
//
// MATH.  The contact-7 chart h = 1 - (7/2)x + a x^2 + b x^3, f = (h^2+(x-1)^7)/x^2 has
// rational Weierstrass points r = 1-v^2 with v a root of
//     Q(v) = v^5 + c4 v^4 + (2c4-1) v^3 + (c4+c0-1/2) v^2 + 2 c0 v + c0 ,  c4=b+2, c0=5/2-(a+b).
// Q satisfies IDENTICALLY  Q(-1) = -1/2,  Q'(-1) = 3,  Q'(0) = 2 Q(0);  equivalently, with
//     z_i = 1/(1+v_i):      sum z_i = 6,     prod z_i = 2,     sum 1/(1-z_i) = 3.
// Given two of the five z's, alpha and beta, the other three are the roots of
//     w^3 - A w^2 + B w - C ,   A = 6-alpha-beta,  C = 2/(alpha*beta),
//     T = 3 - 1/(1-alpha) - 1/(1-beta),   B = (3 - 2A - T + T(A+C))/(T-1).
// 3 rational z  <-> quintic type [1,1,1,2] -> 2-rank 3 -> torsion [2,2,14]  (order 56)
// 5 rational z  <-> quintic type [1,1,1,1,1] -> 2-rank 4 -> torsion [2,2,2,14] (order 112, RECORD)
//
// SIEVE.  Two chains sharing the same modular (A,B,C):
//   ROOT  chain: "cubic has a root mod p"      (survival 2/3 per prime)  -> [2,2,14] candidates
//   SPLIT chain: "cubic splits completely mod p" = has a root AND disc is a QR
//                                                  (survival 1/6 per prime) -> order-112 candidates
// The split chain is exactly the disc-is-a-square Legendre test the brief asks for, fused with
// the root test.  ~78 primes are used, so a genuinely-irreducible cubic survives the root chain
// with probability (2/3)^78 = 3.4e-14.
//
// Usage: claude_ov_c7z_sieve H NTHREADS [SHARD NSHARD]

// IMPLEMENTATION NOTE (v2).  For a fixed prime p, A, B, C -- and therefore both the
// "cubic has a root mod p" and "cubic splits mod p" predicates -- depend ONLY on the residue
// pair (z_i mod p, z_j mod p).  So each prime is compiled once into a (p+1)^2 byte table
// (row/col index p = "z is 0, 1 or infinity mod p", i.e. prime unusable for that index)
// with bit0 = has-a-root, bit1 = splits completely.  The whole sieve is then one table
// lookup per prime per pair, which is ~15x faster than doing modular arithmetic inline.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef int32_t i32; typedef int64_t i64; typedef uint32_t u32; typedef uint16_t u16;

#define MAXP 160
static int NP=0; static i32 PR[MAXP];
static uint8_t *TAB[MAXP];        // (p+1)^2 bytes: bit0 hasroot, bit1 split
static u16 *ZC[MAXP];             // per-index residue code, prime-major

static i32 mpw(i32 a, i32 e, i32 p){ i64 r=1,x=a%p; while(e){ if(e&1) r=r*x%p; x=x*x%p; e>>=1;} return (i32)r; }
static int isprime(int n){ if(n<2) return 0; for(int d=2;(i64)d*d<=n;d++) if(n%d==0) return 0; return 1; }

static void buildtab(int k){
  i32 p = PR[k];
  i32 *inv = (i32*)malloc(sizeof(i32)*p); inv[0]=0;
  for(i32 t=1;t<p;t++) inv[t]=mpw(t,p-2,p);
  uint8_t *qr = (uint8_t*)calloc(p,1); qr[0]=1;
  for(i32 t=1;t<p;t++) qr[(i32)(((i64)t*t)%p)]=1;
  size_t n1 = (size_t)p+1;
  uint8_t *T = (uint8_t*)malloc(n1*n1);
  memset(T, 3, n1*n1);                       // default PASS both chains
  for(i32 zi=2; zi<p; zi++) for(i32 zj=2; zj<p; zj++){
    // z=0 and z=1 are the degenerate residues; loop starts at 2 so they keep the PASS default
    i64 A = (6 - zi - zj) % p; if(A<0) A+=p;
    i64 C = (2LL*inv[zi]) % p; C = (C*inv[zj]) % p;
    i64 Tt;
    i32 m1i = (i32)(((1-(i64)zi)%p+p)%p);
    i32 m1j = (i32)(((1-(i64)zj)%p+p)%p);
    if(m1i==0 || m1j==0) continue;            // z=1 mod p: PASS
    Tt = (3 - inv[m1i] - inv[m1j]) % p; if(Tt<0) Tt+=p;
    i64 Bd = (Tt-1+p) % p;
    if(Bd==0) continue;                       // PASS
    i64 Bn = (3 - 2*A - Tt + (Tt*((A+C)%p))%p) % p; if(Bn<0) Bn+=p;
    i64 B = (Bn * inv[(i32)Bd]) % p;
    i32 a2 = (i32)((p-A)%p), a1 = (i32)B, a0 = (i32)((p-C)%p);
    int has = 0;
    for(i32 w=0; w<p; w++){ i64 t = (((i64)w + a2)*w + a1) % p; t = (t*w + a0) % p; if(t==0){ has=1; break; } }
    // disc of w^3 - A w^2 + B w - C
    i64 A2=(A*A)%p, A3=(A2*A)%p, B2=(B*B)%p, B3=(B2*B)%p, C2=(C*C)%p;
    i64 d = (18*((A*B)%p)%p*C - 4*A3%p*C + A2*B2 - 4*B3 - 27*C2) % p; d%=p; if(d<0) d+=p;
    int spl = has && qr[(i32)d];
    T[(size_t)zi*n1 + zj] = (uint8_t)((has?1:0) | (spl?2:0));
  }
  TAB[k]=T; free(inv); free(qr);
}

int main(int argc, char **argv){
  int H      = argc>1 ? atoi(argv[1]) : 40;
  int nthr   = argc>2 ? atoi(argv[2]) : 1;
  int NPRIME = argc>3 ? atoi(argv[3]) : 78;
  int JB     = argc>4 ? atoi(argv[4]) : 2048;
  omp_set_num_threads(nthr);
  { int q=5; while(NP<NPRIME && NP<MAXP){ if(isprime(q)) PR[NP++]=q; q+=2; } }
#pragma omp parallel for schedule(dynamic,1)
  for(int k=0;k<NP;k++) buildtab(k);
  fprintf(stderr,"primes %d: %d ... %d\n", NP, PR[0], PR[NP-1]);

  int cap = 3*H*H + 64;
  i32 *ZN=(i32*)malloc(sizeof(i32)*cap), *ZD=(i32*)malloc(sizeof(i32)*cap);
  i32 *VN=(i32*)malloc(sizeof(i32)*cap), *VD=(i32*)malloc(sizeof(i32)*cap);
  int nv=0;
  for(int d=1;d<=H;d++) for(int n=-H;n<=H;n++){
    int a=n<0?-n:n, b=d; while(b){int t=a%b;a=b;b=t;}
    if(a!=1) continue;
    if(n==0||n==d||n==-d) continue;
    i32 zn=d, zd=n+d;
    if(zd<0){ zn=-zn; zd=-zd; }
    VN[nv]=n; VD[nv]=d; ZN[nv]=zn; ZD[nv]=zd; nv++;
  }
  fprintf(stderr,"H=%d #v=%d pairs=%.6g\n",H,nv,(double)nv*(nv-1)/2);
#pragma omp parallel for schedule(dynamic,1)
  for(int k=0;k<NP;k++){
    i32 p=PR[k];
    u16 *c=(u16*)malloc(sizeof(u16)*nv);
    i32 *inv=(i32*)malloc(sizeof(i32)*p); inv[0]=0;
    for(i32 t=1;t<p;t++) inv[t]=mpw(t,p-2,p);
    for(int t=0;t<nv;t++){
      i32 zn=(i32)(ZN[t]%p), zd=(i32)(ZD[t]%p);
      if(zn<0) zn+=p; if(zd<0) zd+=p;
      if(zd==0 || zn==0) { c[t]=(u16)p; continue; }
      i32 z=(i32)(((i64)zn*inv[zd])%p);
      c[t] = (z==0||z==1) ? (u16)p : (u16)z;
    }
    ZC[k]=c; free(inv);
  }

  int nbl=(nv+JB-1)/JB;
  long long tot_pairs=0, tot_root=0, tot_split=0, blkdone=0;
#pragma omp parallel for schedule(dynamic,1) reduction(+:tot_pairs,tot_root,tot_split)
  for(int tt=0; tt<nbl; tt++){
    int jb = nbl-1-tt;                       // largest blocks first (LPT)
    int j0=jb*JB, j1=j0+JB; if(j1>nv) j1=nv;
    long long lp=0, lr=0, ls=0;
    u16 ci[MAXP];
    for(int i=0;i<j1;i++){
      i32 vni=VN[i], vdi=VD[i];
      for(int k=0;k<NP;k++) ci[k]=ZC[k][i];
      int js=(i+1>j0)?i+1:j0;
      for(int j=js;j<j1;j++){
        if(VD[j]==vdi && (VN[j]==vni || VN[j]==-vni)) continue;   // v_j = +- v_i -> singular f
        lp++;
        int alive=1, spl=1, np=0, ns=0;
        for(int k=0;k<NP;k++){
          i32 p=PR[k];
          uint8_t e = TAB[k][(size_t)ci[k]*(p+1) + ZC[k][j]];
          if(ci[k]==p || ZC[k][j]==p) continue;
          np++;
          if(!(e&1)){ alive=0; spl=0; break; }
          if(spl){ if(e&2) ns++; else spl=0; }
        }
        if(alive){
          lr++; if(spl) ls++;
#pragma omp critical
          { printf("%s v1=%d/%d v2=%d/%d nprime=%d nsplit=%d\n", spl?"SPLITALL_CAND":"R3",
                   vni,vdi,VN[j],VD[j],np,ns); fflush(stdout); }
        }
      }
    }
    tot_pairs+=lp; tot_root+=lr; tot_split+=ls;
#pragma omp atomic
    blkdone++;
    if((blkdone & 255)==0){
#pragma omp critical
      { fprintf(stderr,"PROGRESS blocks=%lld/%d pairs=%.4g\n",blkdone,nbl,(double)tot_pairs); fflush(stderr); }
    }
  }
  printf("SEARCH_DONE H=%d pairs=%lld root_survivors=%lld split_survivors=%lld nprimes=%d\n",
         H,tot_pairs,tot_root,tot_split,NP);
  fflush(stdout);
  return 0;
}
