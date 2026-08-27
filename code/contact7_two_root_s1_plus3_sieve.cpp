// One-parameter local sieve for the cancelled s=1 branch of the
// contact-7 two-root surface.
//
//   a = 35/8,
//   b = -(8t^4+24t^3+48t^2+45t+15)/(8(t+1)^3),
//
// and
//
//   f=x^5+(b^2-7)x^4+(2ab+21)x^3
//       +(a^2-7b-35)x^2+(2b-7a+35)x+(2a-35/4).
//
// The roots x=0 and x=1-t^2 give visible Z/2 x Z/14.  A rational
// 3-class forces 84 (prime-to-p) to divide #J(F_p) at every good prime.
// Only good-reduction failures are discarded; t=0,+-1, denominator
// disks, and bad reduction are kept as unknown.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using std::int64_t;
using std::uint64_t;
using std::vector;

static int md(int64_t a,int p){ int r=int(a%p); return r<0?r+p:r; }
static int add(int a,int b,int p){ return md(int64_t(a)+b,p); }
static int sub(int a,int b,int p){ return md(int64_t(a)-b,p); }
static int mul(int a,int b,int p){ return int(int64_t(a)*b%p); }
static int pw(int a,int64_t n,int p){ int r=1; for(;n;n>>=1,a=mul(a,a,p)) if(n&1) r=mul(r,a,p); return r; }
static int inv(int a,int p){ return pw(a,p-2,p); }
static int chi(int a,int p){ if(!a)return 0; return pw(a,(p-1)/2,p)==1?1:-1; }

static void trim(vector<int>&a){ while(a.size()>1 && a.back()==0)a.pop_back(); }
static vector<int> rem(vector<int>a,vector<int>b,int p){
    trim(a);trim(b); int db=int(b.size())-1,ib=inv(b.back(),p);
    while(int(a.size())-1>=db && !(a.size()==1&&a[0]==0)){
        int d=int(a.size())-1-db,q=mul(a.back(),ib,p);
        for(int j=0;j<=db;++j)a[d+j]=sub(a[d+j],mul(q,b[j],p),p);
        trim(a);
    }
    return a;
}
static bool squarefree(const std::array<int,6>&f,int p){
    vector<int>a(f.begin(),f.end()),b(5);
    for(int i=1;i<=5;++i)b[i-1]=mul(i%p,f[i],p);
    trim(a);trim(b);
    while(!(b.size()==1&&b[0]==0)){auto r=rem(a,b,p);a.swap(b);b.swap(r);}
    return a.size()==1;
}

struct F2{int a,b;};
static F2 fadd(F2 x,F2 y,int p){return{add(x.a,y.a,p),add(x.b,y.b,p)};}
static F2 fmul(F2 x,F2 y,int d,int p){return{add(mul(x.a,y.a,p),mul(d,mul(x.b,y.b,p),p),p),add(mul(x.a,y.b,p),mul(x.b,y.a,p),p)};}
static F2 eval2(const std::array<int,6>&f,F2 x,int d,int p){F2 y{f[5],0};for(int i=4;i>=0;--i)y=fadd(fmul(y,x,d,p),{f[i],0},p);return y;}
static int chi2(F2 x,int d,int p){if(!x.a&&!x.b)return 0;return chi(sub(mul(x.a,x.a,p),mul(d,mul(x.b,x.b,p),p),p),p);}

static bool model(int t,int p,std::array<int,6>&f){
    if(p==2 || t==0 || t==1 || t==p-1)return false;
    int den=mul(8%p,pw(add(t,1,p),3,p),p); if(!den)return false;
    int t2=mul(t,t,p),t3=mul(t2,t,p),t4=mul(t2,t2,p);
    int num=add(add(mul(8%p,t4,p),mul(24%p,t3,p),p),add(mul(48%p,t2,p),add(mul(45%p,t,p),15%p,p),p),p);
    int b=mul(sub(0,num,p),inv(den,p),p), a=mul(35%p,inv(8%p,p),p);
    int half=inv(2,p),quarter=mul(half,half,p);
    f[0]=sub(mul(2,a,p),mul(35%p,quarter,p),p);
    f[1]=add(sub(mul(2,b,p),mul(7%p,a,p),p),35%p,p);
    f[2]=sub(sub(mul(a,a,p),mul(7%p,b,p),p),35%p,p);
    f[3]=add(mul(2,mul(a,b,p),p),21%p,p);
    f[4]=sub(mul(b,b,p),7%p,p); f[5]=1;
    return true;
}
static int jac(const std::array<int,6>&f,int p){
    int S1=0;for(int x=0;x<p;++x){int y=f[5];for(int i=4;i>=0;--i)y=add(mul(y,x,p),f[i],p);S1+=chi(y,p);} int a1=-S1;
    int d=2;while(chi(d,p)!=-1)++d; int S2=0;
    for(int a=0;a<p;++a)for(int b=0;b<p;++b)S2+=chi2(eval2(f,{a,b},d,p),d,p);
    int a2=(S2+a1*a1)/2;return p*p+1+a2-(p+1)*a1;
}
static int req(int n,int p){while(n%p==0)n/=p;return n;}

struct Mask{int p;vector<unsigned char>s;uint64_t chart=0,good=0,pass=0,fail=0,anomaly=0;};
static Mask mask(int p){
    Mask M{p,vector<unsigned char>(p,0)};
    for(int t=0;t<p;++t){std::array<int,6>f{};if(!model(t,p,f))continue;++M.chart;if(!squarefree(f,p))continue;++M.good;int n=jac(f,p);if(n%req(28,p))++M.anomaly;M.s[t]=(n%req(84,p)==0)?1:2;if(M.s[t]==1)++M.pass;else++M.fail;}
    return M;
}
static int residue(int n,int d,int p){int dd=md(d,p);return dd?mul(md(n,p),inv(dd,p),p):-1;}

int main(int argc,char**argv){
    int H=argc>1?std::atoi(argv[1]):1000; std::ofstream log,sv;
    if(argc>2)log.open(argv[2]); if(argc>3)sv.open(argv[3]); std::ostream&out=log.is_open()?log:std::cout;
    const vector<int>ps={5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251};
    vector<Mask>ms;out<<"CONTACT7_TWO_ROOT_S1_PLUS3_SIEVE\nheight "<<H<<"\n";
    for(int p:ps){Mask M=mask(p);out<<"prime "<<p<<" chart "<<M.chart<<" good "<<M.good<<" pass "<<M.pass<<" fail "<<M.fail<<" unknown "<<p-M.good<<" base_anomaly "<<M.anomaly<<" required "<<req(84,p)<<"\n";ms.push_back(std::move(M));}
    // Most selective masks first.
    std::stable_sort(ms.begin(),ms.end(),[](const Mask&A,const Mask&B){return (__int128)A.fail*B.p>(__int128)B.fail*A.p;});
    out<<"filter_order";for(auto&M:ms)out<<" "<<M.p;out<<"\n";
    uint64_t raw=0,boundary=0,tested=0,survived=0;vector<uint64_t>kills(ms.size());
    if(sv.is_open())sv<<"[\n";bool first=true;
    for(int d=1;d<=H;++d)for(int n=-H;n<=H;++n){if(std::gcd(std::abs(n),d)!=1)continue;++raw;if(n==0||n==d||n==-d){++boundary;continue;}++tested;size_t k=0;for(;k<ms.size();++k){int r=residue(n,d,ms[k].p);if(r<0)continue;if(ms[k].s[r]==2)break;}if(k<ms.size()){++kills[k];continue;}++survived;if(sv.is_open()){if(!first)sv<<",\n";first=false;sv<<"<"<<n<<","<<d<<">";}if(survived<=100)out<<"SURVIVOR "<<n<<"/"<<d<<"\n";}
    if(sv.is_open())sv<<"\n]\n";
    out<<"SUMMARY rational_parameters "<<raw<<" exact_boundary "<<boundary<<" tested "<<tested<<" survivors "<<survived<<"\nKILL_COUNTS";for(size_t k=0;k<ms.size();++k)if(kills[k])out<<" "<<ms[k].p<<":"<<kills[k];out<<"\n";
}
