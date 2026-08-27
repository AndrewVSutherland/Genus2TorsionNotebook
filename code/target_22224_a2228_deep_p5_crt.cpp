// CRT/rational-reconstruction search on the certified p^5 direct-contact
// tubes for [2,2,2,24].
//
// The p=11 bank contains 484 projected tangent states modulo 11^5 and the
// p=13 bank contains 52 modulo 13^5.  We quotient the p=11 side by a global
// signed permutation/projective unit, enumerate every relative signed
// permutation on the p=13 side, combine the three affine ratios by CRT, and
// rationally reconstruct them at the sharp sqrt(M/2) bound.  Every resulting
// rational A2228 tuple is independently filtered and exact-square tested.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using i64=long long;
using i128=__int128_t;

struct Options {
    std::string p11="results/target_22224_direct_contact_deep13_padic_tangent_p11.tsv";
    std::string p13="results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv";
    std::string mask_prefix="results/target_22224_full_family_halving_corrected_complete";
    std::string output="results/target_22224_a2228_deep_p5_crt.tsv";
    size_t start=0,stop=(size_t)-1,progress=25;
    int prime_stop=97;
};

static std::vector<std::string> split(const std::string&s,char sep){
    std::vector<std::string>z;std::stringstream ss(s);std::string x;
    while(std::getline(ss,x,sep))z.push_back(x);return z;
}
static i64 modll(i64 a,i64 m){a%=m;return a<0?a+m:a;}
static i64 invunit(i64 a,i64 m){
    i64 r0=m,r1=modll(a,m),s0=0,s1=1;
    while(r1){i64 q=r0/r1; i64 nr=r0-q*r1;r0=r1;r1=nr;
        i64 ns=s0-q*s1;s0=s1;s1=ns;}
    if(r0!=1)throw std::runtime_error("inverse of nonunit");
    return modll(s0,m);
}
static i64 absres(i64 a,i64 m){a=modll(a,m);return std::min(a,a?m-a:0);}

static std::vector<std::array<i64,4>> read_bases(const std::string&path){
    std::ifstream in(path);if(!in)throw std::runtime_error("cannot open "+path);
    std::string line;std::getline(in,line);auto h=split(line,'\t');
    int col[4]={-1,-1,-1,-1};const char*nm[4]={"a","b","c","d"};
    for(int j=0;j<4;j++)for(size_t i=0;i<h.size();i++)if(h[i]==nm[j])col[j]=(int)i;
    std::vector<std::array<i64,4>>out;
    while(std::getline(in,line)){if(line.empty())continue;auto f=split(line,'\t');
        std::array<i64,4>v{};for(int j=0;j<4;j++)v[j]=std::stoll(f[col[j]]);out.push_back(v);}
    return out;
}

static std::array<i64,4> canonical_abs_projective(
    const std::array<i64,4>&v,int p,i64 m){
    std::array<i64,4>best={m,m,m,m};bool found=false;
    for(int i=0;i<4;i++)if(v[i]%p){i64 q=invunit(v[i],m);std::array<i64,4>z{};
        for(int j=0;j<4;j++)z[j]=absres((i64)((i128)modll(v[j],m)*q%m),m);
        std::sort(z.begin(),z.end());if(!found||z<best){best=z;found=true;}}
    if(!found)throw std::runtime_error("deep base has no unit coordinate");return best;
}

static std::set<std::array<i64,3>> oriented_ratios(
    const std::vector<std::array<i64,4>>&rows,int p,i64 m){
    std::set<std::array<i64,3>>out;
    const int perms[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for(auto v:rows)for(int pivot=0;pivot<4;pivot++)if(v[pivot]%p){
        i64 q=invunit(v[pivot],m);std::array<i64,3>r{};int k=0;
        for(int j=0;j<4;j++)if(j!=pivot)r[k++]=(i64)((i128)modll(v[j],m)*q%m);
        for(auto &pi:perms)for(int mask=0;mask<8;mask++){
            std::array<i64,3>z{};for(int j=0;j<3;j++){
                z[j]=r[pi[j]];if(mask>>j&1)z[j]=z[j]?m-z[j]:0;
            }out.insert(z);
        }
    }return out;
}

struct Rat{i64 n=0,d=1;bool ok=false;};
static Rat ratrec(i64 a,i64 m,i64 B){
    i64 r0=m,r1=modll(a,m),t0=0,t1=1;
    while(r1>B){i64 q=r0/r1; i64 nr=(i64)((i128)r0-(i128)q*r1);
        i64 nt=(i64)((i128)t0-(i128)q*t1);r0=r1;r1=nr;t0=t1;t1=nt;}
    if(!t1)return {};
    i64 n=r1,d=t1;if(d<0){n=-n;d=-d;}i64 g=std::gcd(std::llabs(n),d);n/=g;d/=g;
    if(std::llabs(n)>B||d>B||std::gcd(d,m)!=1)return {};
    if(modll((i64)(((i128)n-(i128)a*d)%m),m))return {};
    return {n,d,true};
}

static i64 crt(i64 a,i64 m,i64 b,i64 n){
    i64 q=modll((i64)((i128)(b-a)*invunit(m,n)%n),n);
    return (i64)((i128)a+(i128)m*q);
}

struct LocalMask{
    int p;std::vector<unsigned char>a;
    size_t ix(int x,int y,int z,int w)const{return (((size_t)x*p+y)*p+z)*p+w;}
    bool get(int x,int y,int z,int w)const{return a[ix(x,y,z,w)]!=0;}
};
static bool boundary4(int a,int b,int c,int d,int p){
    int s[5]={0,a*a%p,b*b%p,c*c%p,d*d%p};for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)if(s[i]==s[j])return true;return false;
}
static LocalMask make_mask(int p,const std::string&pre){
    LocalMask M{p,std::vector<unsigned char>((size_t)p*p*p*p,0)};
    for(int a=0;a<p;a++)for(int b=0;b<p;b++)for(int c=0;c<p;c++)for(int d=0;d<p;d++)
        if(boundary4(a,b,c,d,p))M.a[M.ix(a,b,c,d)]=1;
    std::ifstream in(pre+"_p"+std::to_string(p)+".tsv");if(!in)throw std::runtime_error("mask open failed");
    std::string line;std::getline(in,line);while(std::getline(in,line)){if(line.empty())continue;auto f=split(line,'\t');
        int v[4];for(int i=0;i<4;i++)v[i]=std::stoi(f[i])%p;
        for(int l=1;l<p;l++)M.a[M.ix(l*v[0]%p,l*v[1]%p,l*v[2]%p,l*v[3]%p)]=1;}
    return M;
}

static bool square_mod(int x,int p){x%=p;if(x<0)x+=p;for(int r=0;r<p;r++)if(r*r%p==x)return true;return false;}
static bool modular_pass(const std::array<Rat,3>&q,int signs,const std::vector<LocalMask>&masks,
                         const std::vector<int>&square_primes,std::map<int,unsigned long long>&kills){
    for(int p:square_primes){bool unit=true;int x[4]={1,0,0,0};
        for(int j=0;j<3;j++){if(q[j].d%p==0){unit=false;break;}i64 z=(i64)((i128)modll(q[j].n,p)*invunit(q[j].d,p)%p);
            if(signs>>j&1)z=z?p-z:0;x[j+1]=(int)z;}
        if(!unit)continue;
        int a=x[0],b=x[1],c=x[2],d=x[3];int R[4]={
            (int)((i64)a*b%p*c%p*d%p),
            (int)((i64)a*(a+b)%p*(a+c)%p*(a+d)%p),
            (int)((i64)b*(b+a)%p*(b+c)%p*(b+d)%p),
            (int)((i64)c*(c+a)%p*(c+b)%p*(c+d)%p)};
        bool ok=true;for(int r:R)if(!square_mod(r,p)){ok=false;break;}
        if(!ok){kills[p]++;return false;}
        for(auto&M:masks)if(M.p==p&&!M.get(a,b,c,d)){kills[-p]++;return false;}
    }return true;
}

static bool smooth_offrectangle(const std::array<Rat,3>&q,int signs){
    std::array<i64,4>n={1,q[0].n,q[1].n,q[2].n},d={1,q[0].d,q[1].d,q[2].d};
    for(int j=1;j<4;j++)if(signs>>(j-1)&1)n[j]=-n[j];
    for(int i=0;i<4;i++)if(!n[i])return false;
    for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if((i128)n[i]*d[j]==(i128)n[j]*d[i]||(i128)n[i]*d[j]==-(i128)n[j]*d[i])return false;
    auto eqprod=[&](int i,int j,int k,int l){return (i128)std::llabs(n[i])*std::llabs(n[j])*d[k]*d[l]==(i128)std::llabs(n[k])*std::llabs(n[l])*d[i]*d[j];};
    return !(eqprod(0,1,2,3)||eqprod(0,2,1,3)||eqprod(0,3,1,2));
}

static Options parse(int ac,char**av){Options o;for(int i=1;i<ac;i++){std::string a=av[i];auto need=[&](){if(++i>=ac)throw std::runtime_error("missing option");return std::string(av[i]);};
    if(a=="--p11")o.p11=need();else if(a=="--p13")o.p13=need();else if(a=="--mask-prefix")o.mask_prefix=need();
    else if(a=="--output")o.output=need();else if(a=="--start")o.start=std::stoull(need());else if(a=="--stop")o.stop=std::stoull(need());
    else if(a=="--progress")o.progress=std::stoull(need());else if(a=="--prime-stop")o.prime_stop=std::stoi(need());
    else throw std::runtime_error("unknown option "+a);}return o;}

int main(int ac,char**av){try{
    Options o=parse(ac,av);const int p11=11,p13=13;const i64 m11=161051,m13=371293,M=m11*m13;
    i64 B=(i64)std::sqrt((long double)(M/2));while((B+1)*(B+1)<=M/2)++B;while(B*B>M/2)--B;
    auto A=read_bases(o.p11),C=read_bases(o.p13);std::set<std::array<i64,3>>left;
    for(auto v:A){auto z=canonical_abs_projective(v,p11,m11);if(z[0]!=1)throw std::runtime_error("canonical pivot failure");left.insert({z[1],z[2],z[3]});}
    auto right=oriented_ratios(C,p13,m13);std::vector<std::array<i64,3>>R(right.begin(),right.end());
    size_t stop=std::min(o.stop,R.size());if(o.start>stop)o.start=stop;
    std::vector<LocalMask>masks;for(int p:{17,19,23,29,31})masks.push_back(make_mask(p,o.mask_prefix));
    std::vector<int>sqpr={17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
    sqpr.erase(std::remove_if(sqpr.begin(),sqpr.end(),[&](int p){return p>o.prime_stop;}),sqpr.end());
    std::ofstream out(o.output);if(!out)throw std::runtime_error("cannot write output");
    out<<"left_index\tright_index\tsigns\tq1_num\tq1_den\tq2_num\tq2_den\tq3_num\tq3_den\n";
    unsigned long long pairs=0,reconstructed=0,modpass=0;std::map<int,unsigned long long>kills;
    std::cout<<"A2228_DEEP_P5_CRT_START modulus="<<M<<" ratrec_bound="<<B<<" p11_rows="<<A.size()
             <<" p11_profiles="<<left.size()<<" p13_rows="<<C.size()<<" p13_oriented="<<R.size()
             <<" right_range=["<<o.start<<","<<stop<<") prime_stop="<<o.prime_stop
             <<" model=independent_unit_scaling_signed_permutations\n";
    size_t li=0;for(auto L:left){for(size_t ri=o.start;ri<stop;ri++){++pairs;std::array<Rat,3>q{};bool ok=true;
        for(int j=0;j<3;j++){q[j]=ratrec(crt(L[j],m11,R[ri][j],m13),M,B);if(!q[j].ok){ok=false;break;}}
        if(!ok)continue;++reconstructed;
        for(int signs=0;signs<8;signs++){if(!smooth_offrectangle(q,signs))continue;
            if(!modular_pass(q,signs,masks,sqpr,kills))continue;++modpass;
            out<<li<<'\t'<<ri<<'\t'<<signs;for(auto z:q)out<<'\t'<<z.n<<'\t'<<z.d;out<<'\n';
        }}++li;if(o.progress&&li%o.progress==0){out.flush();std::cout<<"PROGRESS left="<<li<<" pairs="<<pairs<<" reconstructed="<<reconstructed<<" modpass="<<modpass<<"\n";}}
    std::cout<<"A2228_DEEP_P5_CRT_DONE pairs="<<pairs<<" reconstructed="<<reconstructed<<" modpass="<<modpass<<" kills=";
    for(auto z:kills)std::cout<<z.first<<":"<<z.second<<",";std::cout<<" output="<<o.output<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<"\n";return 2;}}
