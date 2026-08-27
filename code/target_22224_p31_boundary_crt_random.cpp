// Boundary-conditioned Monte Carlo/CRT search for the q-square
// [2,2,2,24] cover.  Unlike the old square box, every normalized trial is
// drawn from one of the COMPLETE p=31 boundary residue classes; toric trials
// force a p=31 zero or pole.  Complete target-or-boundary masks at seven
// other primes are applied before exact 128-bit square tests.
//
// Build:
//   c++ -O3 -std=c++17 -o /tmp/target_22224_p31_boundary_crt_random \
//     code/target_22224_p31_boundary_crt_random.cpp
// Run (trials height seed output log):
//   /tmp/target_22224_p31_boundary_crt_random 100000000 10000 1 \
//     results/target_22224_p31_boundary_crt_candidates.tsv \
//     results/target_22224_p31_boundary_crt_random.log

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using i128=__int128_t; using u128=__uint128_t;

struct Frac { int n=1,d=1; };
struct Mask { int p; std::vector<unsigned char> target; };

static int norm(long long x,int p){int r=int(x%p);return r<0?r+p:r;}
static int pw(int a,int e,int p){long long r=1,b=norm(a,p);while(e){if(e&1)r=r*b%p;b=b*b%p;e>>=1;}return int(r);}
static int inv(int a,int p){return pw(a,p-2,p);}
static bool sq0(int a,int p){a=norm(a,p);return a==0||pw(a,(p-1)/2,p)==1;}

static std::vector<std::string> tabs(const std::string&s){std::vector<std::string>v;std::stringstream z(s);std::string x;while(std::getline(z,x,'\t'))v.push_back(x);return v;}
static Mask load(int p){
  std::string stem=p<=41?"results/target_22224_qsquare_crt_sieve":"results/target_22224_p31_boundary_crt_mask";
  std::ifstream in(stem+"_p"+std::to_string(p)+"_masks.tsv");
  if(!in){std::cerr<<"missing mask p="<<p<<"\n";std::exit(2);} Mask m{p,std::vector<unsigned char>(p*p)};std::string l;std::getline(in,l);
  while(std::getline(in,l)){auto f=tabs(l);if(f.size()>=11&&std::stoi(f[10])>0)m.target[std::stoi(f[1])*p+std::stoi(f[2])]=1;}return m;
}
static bool raw(const Mask&m,int a,int b){a=norm(a,m.p);b=norm(b,m.p);return a&&b&&m.target[a*m.p+b];}
static bool symtarget(const Mask&m,int a,int b,int c){
  int p=m.p,x[3]={a,b,c}; const int pe[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}; const int sg[4][3]={{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
  for(auto&s:sg){int y[3]={norm(s[0]*x[0],p),norm(s[1]*x[1],p),norm(s[2]*x[2],p)};for(auto&q:pe)if(raw(m,y[q[0]],y[q[1]]))return true;}return false;
}
enum Stat{BOUNDARY,TARGET,REJECT};
static Stat status(const Mask&m,const Frac&A,const Frac&B){
 int p=m.p,an=A.n%p,ad=A.d%p,bn=B.n%p,bd=B.d%p;if(!an||!ad||!bn||!bd)return BOUNDARY;
 int a=(long long)an*inv(ad,p)%p,b=(long long)bn*inv(bd,p)%p,c=inv((long long)a*b%p,p);int x=(long long)a*a%p,y=(long long)b*b%p,z=(long long)c*c%p;
 int R=norm(x+y+z-3,p),S=norm(inv(x,p)+inv(y,p)+inv(z,p)-3,p);if(!sq0(R,p)||!sq0(S,p))return REJECT;if(!R||!S)return BOUNDARY;
 int w=(long long)S*inv(R,p)%p;int rr[4]={x,y,z,w};std::sort(rr,rr+4);if(std::adjacent_find(rr,rr+4)!=rr+4)return BOUNDARY;
 return symtarget(m,a,b,c)?TARGET:REJECT;
}

struct RNG{uint64_t s;uint64_t next(){s^=s<<7;s^=s>>9;s^=s<<8;return s;}int pick(int n){return int(next()%uint64_t(n));}};
static Frac random_frac(RNG&r,int H,bool unit31=true){for(;;){Frac f{1+r.pick(H),1+r.pick(H)};if(std::gcd(f.n,f.d)!=1)continue;if(unit31&&(f.n%31==0||f.d%31==0))continue;return f;}}
static bool with_residue(RNG&r,int H,int residue,Frac&f){
 for(int tries=0;tries<20;++tries){int d=1+r.pick(H);if(d%31==0)continue;int want=(long long)residue*(d%31)%31;int first=want?want:31; if(first>H)continue;int count=(H-first)/31+1;int n=first+31*r.pick(count);if(std::gcd(n,d)==1){f={n,d};return true;}}return false;
}

static unsigned long long rootfloor(i128 x){if(x<=0)return 0;unsigned long long r=(unsigned long long)std::sqrt((long double)x);while((u128)r*r>(u128)x)--r;while((u128)(r+1)*(r+1)<=(u128)x)++r;return r;}
static bool square(i128 x,unsigned long long&r){if(x<0)return false;r=rootfloor(x);return(u128)r*r==(u128)x;}
struct Exact{bool ok=false;unsigned long long rn,rd,sn,sd;};
static Exact exact(const Frac&A,const Frac&B){
 i128 an=A.n,ad=A.d,bn=B.n,bd=B.d,an2=an*an,ad2=ad*ad,bn2=bn*bn,bd2=bd*bd,D=an2*ad2*bn2*bd2;
 i128 R=bd2*an2*an2*bn2+ad2*an2*bn2*bn2+ad2*ad2*bd2*bd2-3*D;
 i128 S=ad2*ad2*bd2*bn2+ad2*bd2*bd2*an2+an2*an2*bn2*bn2-3*D;unsigned long long rr,ss;if(!square(R,rr)||!square(S,ss)||!rr||!ss)return{};
 unsigned long long den=(unsigned long long)A.n*A.d*B.n*B.d,gr=std::gcd(rr,den),gs=std::gcd(ss,den);return{true,rr/gr,den/gr,ss/gs,den/gs};
}
static bool collision(const Frac&A,const Frac&B){
 i128 an=A.n,ad=A.d,bn=B.n,bd=B.d;
 if(an==ad||bn==bd||an*bn==ad*bd)return true;                 // A,B,C=1
 if(an*bd==bn*ad)return true;                                 // A=B
 if(an*an*bn==ad*ad*bd)return true;                           // A=C
 if(an*bn*bn==ad*bd*bd)return true;                           // B=C
 return false;
}

int main(int ac,char**av){
 uint64_t N=ac>1?std::stoull(av[1]):100000000ULL;int H=ac>2?std::stoi(av[2]):10000;uint64_t seed=ac>3?std::stoull(av[3]):1;
 std::string outn=ac>4?av[4]:"results/target_22224_p31_boundary_crt_candidates.tsv";
 std::string logn=ac>5?av[5]:"results/target_22224_p31_boundary_crt_random.log";
 const std::array<int,10> ps{29,37,41,43,47,53,59,61,67,71};
 std::array<Mask,10> ms{load(29),load(37),load(41),load(43),load(47),
                        load(53),load(59),load(61),load(67),load(71)};
 // Complete normalized p=31 boundary lists, including the four R=S=0 bases.
 std::array<std::vector<int>,31> allow31;
 for(int a=1;a<31;++a)for(int b=1;b<31;++b){int c=inv(a*b%31,31),x=a*a%31,y=b*b%31,z=c*c%31,R=norm(x+y+z-3,31),S=norm(inv(x,31)+inv(y,31)+inv(z,31)-3,31);if(!sq0(R,31)||!sq0(S,31))continue;bool bd=!R||!S;if(!bd){int w=(long long)S*inv(R,31)%31;int q[4]={x,y,z,w};std::sort(q,q+4);bd=std::adjacent_find(q,q+4)!=q+4;}if(bd)allow31[a].push_back(b);}
 RNG rng{seed?seed:1};std::ofstream out(outn),log(logn);out<<"A_num\tA_den\tB_num\tB_den\trho_num\trho_den\tsigma_num\tsigma_den\tmode\n";
 uint64_t made=0,modpass=0,exactn=0,nontriv=0;
 std::array<uint64_t,10> rej{},bd{},tg{};
 uint64_t normalized=0,toric=0,gcdfail=0;
 for(uint64_t it=0;it<N;++it){Frac A,B;int mode=rng.pick(10)<8?0:1;if(mode==0){A=random_frac(rng,H,true);int ar=(long long)(A.n%31)*inv(A.d%31,31)%31;if(allow31[ar].empty()||!with_residue(rng,H,allow31[ar][rng.pick(allow31[ar].size())],B)){++gcdfail;continue;}++normalized;}
   else{A=random_frac(rng,H,false);B=random_frac(rng,H,false);int which=rng.pick(4);if(which==0)A.n=31*(1+rng.pick(H/31));if(which==1)A.d=31*(1+rng.pick(H/31));if(which==2)B.n=31*(1+rng.pick(H/31));if(which==3)B.d=31*(1+rng.pick(H/31));if(A.n>H||A.d>H||B.n>H||B.d>H||std::gcd(A.n,A.d)!=1||std::gcd(B.n,B.d)!=1){++gcdfail;continue;}++toric;}
   ++made;bool pass=true;for(size_t k=0;k<ps.size();++k){Stat s=status(ms[k],A,B);if(s==REJECT){++rej[k];pass=false;break;}if(s==BOUNDARY)++bd[k];else++tg[k];}if(!pass)continue;++modpass;Exact e=exact(A,B);if(!e.ok)continue;++exactn;if(collision(A,B))continue;++nontriv;out<<A.n<<'\t'<<A.d<<'\t'<<B.n<<'\t'<<B.d<<'\t'<<e.rn<<'\t'<<e.rd<<'\t'<<e.sn<<'\t'<<e.sd<<'\t'<<(mode?"toric":"normalized")<<'\n';}
 std::ostringstream h;h<<"P31_BOUNDARY_CRT_RANDOM trials "<<N<<" height "<<H<<" seed "<<seed<<" generated "<<made<<" normalized "<<normalized<<" toric "<<toric<<" generation_fail "<<gcdfail<<" ten_prime_survivors "<<modpass<<" exact_double_square "<<exactn<<" exact_nontrivial "<<nontriv;std::cout<<h.str()<<"\n";log<<h.str()<<"\n";
 for(size_t k=0;k<ps.size();++k){std::ostringstream z;z<<"MASK p "<<ps[k]<<" reject "<<rej[k]<<" boundary "<<bd[k]<<" target "<<tg[k];std::cout<<z.str()<<"\n";log<<z.str()<<"\n";}log<<"OUTPUT "<<outn<<"\n";
}
