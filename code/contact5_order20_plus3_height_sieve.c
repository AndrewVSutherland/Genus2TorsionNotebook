/*
 * Fast exact finite-mask scan for the full contact-5/order-20 +3 lane.
 *
 * Usage:
 *   clang -O3 -o /tmp/order20_plus3_sieve code/contact5_order20_plus3_height_sieve.c
 *   /tmp/order20_plus3_sieve 100000 t
 *   /tmp/order20_plus3_sieve 100000 s
 *
 * The masks are the exact conditions 3 | #J_t(F_p), with every bad
 * residue passed conservatively.  They were independently point-counted
 * over F_p and F_{p^2}; all good residues also satisfy the positive control
 * 20 | #J_t(F_p).  The s-chart means s=(t-1)/(t+1).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXP 127
#define NMASK 27
#define NBASE 4

typedef struct { int p; unsigned char ok[MAXP]; } Mask;
typedef struct { int mod, n; int *a; } Classes;

static const char *spec[NMASK] = {
 "7:4,5,6",
 "11:0,2,8,10",
 "13:0,2,6,10,12",
 "17:2,3,8,13,14,15,16",
 "19:2,3,5,10,12,13,14,15,16,18",
 "23:0,6,7,9,10,11,13,15,16,17,20,22",
 "29:0,1,2,7,11,14,15,23,24,25,26,28",
 "31:3,5,9,10,11,13,15,16,17,21,26,28,29,30",
 "37:0,1,4,7,8,14,17,19,20,22,23,26,28,30,34,36",
 "41:0,3,6,7,14,18,29,31,32,38,40",
 "43:1,2,6,11,12,14,18,29,30,34,37,38,40,41,42",
 "47:1,2,3,4,6,7,9,10,11,12,20,23,27,28,30,35,38,40,42,44,46",
 "53:1,2,3,6,12,14,16,18,21,29,30,32,35,42,43,45,50,52",
 "59:1,3,5,8,11,13,14,17,19,20,21,22,23,34,36,38,40,44,46,48,51,52,55,56,58",
 "61:2,3,5,7,8,9,12,15,16,23,26,27,28,29,37,41,42,46,50,53,57,58,60",
 "67:3,7,12,15,16,22,26,34,35,38,39,42,43,44,53,62,63,64,65,66",
 "71:2,10,11,13,14,15,18,22,24,26,30,33,34,38,42,44,45,54,55,58,59,60,61,63,65,67,68,70",
 "73:2,5,13,14,15,16,17,25,26,36,37,41,42,43,44,46,50,54,56,57,58,69,70,72",
 "79:0,2,7,8,10,18,21,22,24,27,30,34,36,41,44,46,47,50,54,55,58,59,60,62,65,66,67,70,73,74,76,78",
 "83:4,6,7,9,10,11,12,13,15,20,25,28,30,32,34,35,37,38,44,45,46,47,48,50,54,55,61,64,65,66,68,71,74,75,76,77,80,82",
 "89:2,4,5,7,10,12,16,17,20,23,24,27,30,32,34,39,48,52,55,57,58,64,69,71,73,78,82,84,86,88",
 "97:0,2,3,4,7,8,9,10,12,18,19,20,22,25,28,29,35,37,39,40,41,43,44,45,50,57,58,61,70,73,79,84,89,91,94,96",
 "101:1,5,6,8,10,13,20,22,23,25,28,29,32,33,35,36,37,46,47,48,51,55,58,59,60,64,65,66,67,76,82,83,88,91,92,93,98,100",
 "103:6,10,11,15,17,20,21,22,23,25,32,36,38,39,40,42,45,48,51,56,58,63,66,67,68,74,79,84,85,93,94,96,100,101,102",
 "107:0,1,5,6,7,8,10,11,12,17,19,20,21,24,25,26,28,34,35,36,37,40,41,42,48,49,51,59,60,61,64,65,66,79,80,81,83,85,89,92,93,94,95,96,100,103,104,105,106",
 "109:0,3,4,6,8,11,12,13,15,17,18,19,24,27,31,36,40,42,45,48,51,52,59,61,63,64,65,66,69,71,74,76,77,78,80,82,85,92,97,98,99,103,104,106,107,108",
 "113:3,6,14,16,17,20,21,22,23,25,27,28,33,34,35,36,37,39,43,46,49,50,51,55,56,58,64,73,75,76,78,81,83,86,87,92,93,97,99,103,104,107,108,110,112"
};

static int base_index[NBASE] = {0,1,9,15}; /* p=7,11,41,67 */

static long long gcdll(long long a, long long b) {
  if (a < 0) a = -a;
  while (b) { long long r = a % b; a = b; b = r; }
  return a;
}

static int mod(int64_t a, int p) {
  int r = (int)(a % p); return r < 0 ? r + p : r;
}

static int invmod(int a, int p) {
  int t=0, nt=1, r=p, nr=mod(a,p);
  while (nr) { int q=r/nr, x=t-q*nt; t=nt; nt=x; x=r-q*nr; r=nr; nr=x; }
  return mod(t,p);
}

static void parse_masks(Mask *m) {
  for (int i=0;i<NMASK;i++) {
    char buf[1024]; strncpy(buf,spec[i],sizeof(buf)-1); buf[sizeof(buf)-1]=0;
    char *colon=strchr(buf,':'); *colon=0; m[i].p=atoi(buf);
    memset(m[i].ok,0,sizeof(m[i].ok));
    for (char *tok=strtok(colon+1,","); tok; tok=strtok(NULL,","))
      m[i].ok[atoi(tok)]=1;
  }
}

static Classes make_classes(Mask *m, int bits, int schart) {
  Classes c={1,1,malloc(sizeof(int))}; c.a[0]=0;
  for (int j=0;j<NBASE;j++) if (bits&(1<<j)) {
    int p=m[base_index[j]].p, vals[MAXP], nv=0;
    for (int t=0;t<p;t++) if (m[base_index[j]].ok[t]) {
      if (!schart) vals[nv++]=t;
      else if ((t+1)%p) vals[nv++]=mod((int64_t)(t-1)*invmod(t+1,p),p);
    }
    int *na=malloc((size_t)c.n*nv*sizeof(int)); int k=0, inv=invmod(c.mod,p);
    for (int i=0;i<c.n;i++) for (int z=0;z<nv;z++) {
      int delta=mod((int64_t)(vals[z]-c.a[i]%p)*inv,p);
      na[k++]=c.a[i]+c.mod*delta;
    }
    free(c.a); c.a=na; c.n=k; c.mod*=p;
  }
  return c;
}

static int passes(Mask *m, int64_t n, int64_t d, int schart) {
  int64_t tn=schart?d+n:n, td=schart?d-n:d;
  for (int i=0;i<NMASK;i++) {
    int p=m[i].p, dd=mod(td,p);
    if (!dd) continue;
    int tv=mod((int64_t)mod(tn,p)*invmod(dd,p),p);
    if (!m[i].ok[tv]) return 0;
  }
  return 1;
}

int main(int argc,char **argv) {
  if (argc<3) { fprintf(stderr,"usage: %s HEIGHT t|s\n",argv[0]); return 2; }
  int H=atoi(argv[1]), schart=argv[2][0]=='s';
  Mask masks[NMASK]; parse_masks(masks);
  Classes cache[16]; for (int b=0;b<16;b++) cache[b]=make_classes(masks,b,schart);
  uint64_t candidates=0, primitive=0, survivors=0;
  for (int64_t d=1;d<=H;d++) {
    int bits=0;
    for (int j=0;j<NBASE;j++) if (d%masks[base_index[j]].p) bits|=1<<j;
    Classes *c=&cache[bits];
    for (int i=0;i<c->n;i++) {
      int64_t rn=(d*c->a[i])%c->mod;
      int64_t lo=-H, rem=mod(lo,c->mod), delta=mod(rn-rem,c->mod);
      for (int64_t n=lo+delta;n<=H;n+=c->mod) {
        candidates++;
        if (gcdll(n,d)!=1) continue;
        primitive++;
        if (passes(masks,n,d,schart)) {
          survivors++;
          if (!schart) printf("SURVIVOR t=%lld/%lld\n",(long long)n,(long long)d);
          else printf("SURVIVOR s=%lld/%lld t=(%lld)/(%lld)\n",
                      (long long)n,(long long)d,(long long)(d+n),(long long)(d-n));
        }
      }
    }
  }
  printf("DONE chart=%c H=%d CRT_candidates=%llu primitive=%llu survivors=%llu\n",
         schart?'s':'t',H,(unsigned long long)candidates,
         (unsigned long long)primitive,(unsigned long long)survivors);
  for (int b=0;b<16;b++) free(cache[b].a);
  return 0;
}
