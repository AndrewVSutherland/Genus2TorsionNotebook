// Resumable transverse integer-box enumerator for the full A(2,2,2,8)
// four-radicand cover, streamed into the certified p^5 contact tubes.
//
// We projectively order positive magnitudes 0<a<b<c<d<=B and force the
// first cover equation abcd=square by writing
//
//      d = sf(abc) k^2.
//
// The remaining three radicands are tested exactly for each of the eight
// signed presentations modulo global sign.  Full-cover hits are primitive
// normalized, checked against the corrected p=17,...,31 contact masks, and
// matched against the 484 certified 11^5 and 52 certified 13^5 tangent tubes
// up to independent unit scaling, permutations, and sign choices.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using i64=long long;using i128=__int128_t;using u128=__uint128_t;

struct Options{
    int B=400,a_start=1,a_stop=-1;bool all_signs=true;size_t progress=25;
    std::string bank="data/tor2228_bank.txt";
    std::string deep11="results/target_22224_direct_contact_deep13_padic_tangent_p11.tsv";
    std::string deep13="results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv";
    std::string mask_prefix="results/target_22224_full_family_halving_corrected_complete";
    std::string output="results/target_22224_a2228_deep_p5_box.tsv";
    std::string candidates="results/target_22224_a2228_deep_p5_box_candidates.tsv";
};
static std::vector<std::string>split(const std::string&s,char c){std::vector<std::string>z;std::stringstream q(s);std::string x;while(std::getline(q,x,c))z.push_back(x);return z;}
static i64 modll(i64 a,i64 m){a%=m;return a<0?a+m:a;}
static i64 invunit(i64 a,i64 m){i64 r0=m,r1=modll(a,m),s0=0,s1=1;while(r1){i64 q=r0/r1,nr=r0-q*r1,ns=s0-q*s1;r0=r1;r1=nr;s0=s1;s1=ns;}if(r0!=1)throw std::runtime_error("nonunit");return modll(s0,m);}
static i64 absres(i64 a,i64 m){a=modll(a,m);return std::min(a,a?m-a:0);}
static i64 gcd4(i64 a,i64 b,i64 c,i64 d){return std::gcd(std::gcd(std::llabs(a),std::llabs(b)),std::gcd(std::llabs(c),std::llabs(d)));}
static std::array<i64,4>primitive(std::array<i64,4>v){i64 g=gcd4(v[0],v[1],v[2],v[3]);for(auto&x:v)x/=g;for(auto x:v)if(x){if(x<0)for(auto&y:v)y=-y;break;}return v;}
static std::string key(std::array<i64,4>v){v=primitive(v);for(auto&x:v)x=std::llabs(x);std::sort(v.begin(),v.end());return std::to_string(v[0])+","+std::to_string(v[1])+","+std::to_string(v[2])+","+std::to_string(v[3]);}
static u128 isqrt128(u128 n){if(!n)return 0;u128 r=(u128)std::sqrt((long double)n);while((r+1)<=n/(r+1))++r;while(r>n/r)--r;return r;}
static bool square128(i128 x){if(x<0)return false;u128 n=(u128)x,r=isqrt128(n);return r*r==n;}
static i64 isqrt64(i64 n){i64 r=(i64)std::sqrt((long double)n);while((r+1)<=n/(r+1))++r;while(r>n/r)--r;return r;}
static i64 ceil_sqrt64(i64 n){i64 r=isqrt64(n);return r*r==n?r:r+1;}
static std::array<i128,4>rad(const std::array<i64,4>&v){i128 a=v[0],b=v[1],c=v[2],d=v[3];return {a*b*c*d,a*(a+b)*(a+c)*(a+d),b*(b+a)*(b+c)*(b+d),c*(c+a)*(c+b)*(c+d)};}
static bool rectangle(const std::array<i64,4>&v){i128 a=std::llabs(v[0]),b=std::llabs(v[1]),c=std::llabs(v[2]),d=std::llabs(v[3]);return a*b==c*d||a*c==b*d||a*d==b*c;}

struct LocalMask{int p;std::vector<unsigned char>a;size_t ix(int x,int y,int z,int w)const{return (((size_t)x*p+y)*p+z)*p+w;}bool get(const std::array<i64,4>&v)const{int z[4];for(int i=0;i<4;i++)z[i]=(int)modll(v[i],p);return a[ix(z[0],z[1],z[2],z[3])];}};
static bool boundary4(int a,int b,int c,int d,int p){int s[5]={0,a*a%p,b*b%p,c*c%p,d*d%p};for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)if(s[i]==s[j])return true;return false;}
static LocalMask make_mask(int p,const std::string&pre){LocalMask M{p,std::vector<unsigned char>((size_t)p*p*p*p,0)};for(int a=0;a<p;a++)for(int b=0;b<p;b++)for(int c=0;c<p;c++)for(int d=0;d<p;d++)if(boundary4(a,b,c,d,p))M.a[M.ix(a,b,c,d)]=1;std::ifstream in(pre+"_p"+std::to_string(p)+".tsv");if(!in)throw std::runtime_error("mask open");std::string line;std::getline(in,line);while(std::getline(in,line)){if(line.empty())continue;auto f=split(line,'\t');int v[4];for(int i=0;i<4;i++)v[i]=std::stoi(f[i])%p;for(int l=1;l<p;l++)M.a[M.ix(l*v[0]%p,l*v[1]%p,l*v[2]%p,l*v[3]%p)]=1;}return M;}

static std::array<i64,4>deep_key(const std::array<i64,4>&v,int p,i64 m){std::array<i64,4>best={m,m,m,m};bool found=false;for(int i=0;i<4;i++){i64 x=modll(v[i],m);if(x%p==0)continue;i64 q=invunit(x,m);std::array<i64,4>z{};for(int j=0;j<4;j++)z[j]=absres((i64)((i128)modll(v[j],m)*q%m),m);std::sort(z.begin(),z.end());if(!found||z<best){best=z;found=true;}}return best;}
struct DeepBank{int p;i64 m;size_t rows=0;std::set<std::array<i64,4>>keys;bool get(const std::array<i64,4>&v)const{return keys.count(deep_key(v,p,m));}};
static DeepBank read_deep(const std::string&path,int p){DeepBank B{p,1};for(int i=0;i<5;i++)B.m*=p;std::ifstream in(path);if(!in)throw std::runtime_error("deep open");std::string line;std::getline(in,line);auto h=split(line,'\t');int col[4]={-1,-1,-1,-1};const char*nm[4]={"a","b","c","d"};for(int j=0;j<4;j++)for(size_t i=0;i<h.size();i++)if(h[i]==nm[j])col[j]=(int)i;while(std::getline(in,line)){if(line.empty())continue;auto f=split(line,'\t');std::array<i64,4>v{};for(int j=0;j<4;j++)v[j]=std::stoll(f[col[j]]);B.rows++;B.keys.insert(deep_key(v,p,B.m));}return B;}
static std::set<std::string>read_bank(const std::string&path){std::ifstream in(path);std::set<std::string>out;std::string line;while(std::getline(in,line)){if(line.size()<3||line[0]!='[')continue;auto f=split(line.substr(1,line.size()-2),',');if(f.size()!=4)continue;std::array<i64,4>v{};for(int i=0;i<4;i++)v[i]=std::stoll(f[i]);out.insert(key(v));}return out;}

static Options parse(int ac,char**av){Options o;for(int i=1;i<ac;i++){std::string a=av[i];auto need=[&](){if(++i>=ac)throw std::runtime_error("missing option");return std::string(av[i]);};if(a=="--B")o.B=std::stoi(need());else if(a=="--a-start")o.a_start=std::stoi(need());else if(a=="--a-stop")o.a_stop=std::stoi(need());else if(a=="--positive-only")o.all_signs=false;else if(a=="--progress")o.progress=std::stoull(need());else if(a=="--bank")o.bank=need();else if(a=="--deep11")o.deep11=need();else if(a=="--deep13")o.deep13=need();else if(a=="--mask-prefix")o.mask_prefix=need();else if(a=="--output")o.output=need();else if(a=="--candidates")o.candidates=need();else throw std::runtime_error("unknown option "+a);}return o;}

int main(int ac,char**av){
 try{
    Options o=parse(ac,av);
    if(o.a_stop<0)o.a_stop=o.B-2;
    o.a_stop=std::min(o.a_stop,o.B-2);
    auto known=read_bank(o.bank);
    auto D11=read_deep(o.deep11,11),D13=read_deep(o.deep13,13);
    std::vector<LocalMask>masks;
    for(int p:{17,19,23,29,31})masks.push_back(make_mask(p,o.mask_prefix));
    std::vector<int>spf(o.B+1),sf(o.B+1,1);
    for(int i=2;i<=o.B;i++)if(!spf[i])
        for(int j=i;j<=o.B;j+=i)if(!spf[j])spf[j]=i;
    for(int n=2;n<=o.B;n++){
        int p=spf[n],m=n/p;sf[n]=(m%p==0)?sf[m/p]:sf[m]*p;
    }
    auto combine=[](i64 x,i64 y){i64 g=std::gcd(x,y);return (x/g)*(y/g);};
    std::ofstream out(o.output),cand(o.candidates);
    if(!out||!cand)throw std::runtime_error("output open");
    std::string header="a\tb\tc\td\tsigns\tprimitive_a\tprimitive_b\tprimitive_c\tprimitive_d\trectangle\tknown_bank\tlocal_contact\tdeep11\tdeep13\n";
    out<<header;cand<<header;
    unsigned long long triples=0,kvals=0,presentations=0,full=0,offrectangle=0;
    unsigned long long positive_full=0,known_positive=0,newcurves=0,local=0;
    unsigned long long deep11=0,deep13=0,candidates=0;
    std::set<std::string>seennew;
    std::cout<<"A2228_DEEP_P5_BOX_START B="<<o.B<<" a_range=["<<o.a_start<<","<<o.a_stop
             <<"] signs="<<(o.all_signs?8:1)<<" p11_tubes="<<D11.keys.size()
             <<" p13_tubes="<<D13.keys.size()<<" forcing=abcd_squarefree_kernel\n";
    for(int a=o.a_start;a<=o.a_stop;a++){
      for(int b=a+1;b<=o.B-2;b++){
        i64 sab=combine(sf[a],sf[b]);
        for(int c=b+1;c<=o.B-1;c++){
          ++triples;
          i64 s=combine(sab,sf[c]);
          if(s>o.B)continue;
          i64 k0=ceil_sqrt64((c+1+s-1)/s),k1=isqrt64(o.B/s);
          for(i64 k=k0;k<=k1;k++){
            ++kvals;i64 d=s*k*k;
            if(d<=c||d>o.B)continue;
            std::array<i64,4>mag={(i64)a,(i64)b,(i64)c,d};
            auto pmag=primitive(mag);
            bool d11=D11.get(pmag),d13=D13.get(pmag);
            for(int sm=0;sm<(o.all_signs?8:1);sm++){
              ++presentations;
              std::array<i64,4>v=mag;
              for(int j=1;j<4;j++)if(sm>>(j-1)&1)v[j]=-v[j];
              auto rr=rad(v);
              if(!square128(rr[1])||!square128(rr[2])||!square128(rr[3]))continue;
              ++full;if(sm==0)++positive_full;
              bool rect=rectangle(v);if(!rect)++offrectangle;
              auto pv=primitive(v);std::string ck=key(pv);bool old=known.count(ck);
              if(sm==0&&old)++known_positive;
              if(!old&&seennew.insert(ck).second)++newcurves;
              bool lok=!rect;for(auto&M:masks)if(lok&&!M.get(pv)){lok=false;break;}
              local+=lok;deep11+=d11;deep13+=d13;
              out<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<'\t'<<sm;
              for(auto x:pv)out<<'\t'<<x;
              out<<'\t'<<rect<<'\t'<<old<<'\t'<<lok<<'\t'<<d11<<'\t'<<d13<<'\n';
              if(lok&&d11&&d13){
                ++candidates;
                cand<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<'\t'<<sm;
                for(auto x:pv)cand<<'\t'<<x;
                cand<<'\t'<<rect<<'\t'<<old<<'\t'<<lok<<'\t'<<d11<<'\t'<<d13<<'\n';cand.flush();
                std::cout<<"DEEP_P5_FULL_COVER_CANDIDATE tuple=["<<pv[0]<<","<<pv[1]<<","<<pv[2]<<","<<pv[3]<<"] signs="<<sm<<" old="<<old<<"\n";
              }
            }
          }
        }
      }
      if(o.progress&&(a-o.a_start+1)%o.progress==0){
        out.flush();
        std::cout<<"PROGRESS a="<<a<<" triples="<<triples<<" kvals="<<kvals
                 <<" full="<<full<<" positive="<<positive_full<<" local="<<local
                 <<" deep_both="<<candidates<<"\n";
      }
    }
    std::cout<<"A2228_DEEP_P5_BOX_DONE triples="<<triples<<" kvals="<<kvals
             <<" presentations="<<presentations<<" full="<<full
             <<" offrectangle="<<offrectangle
             <<" positive_full="<<positive_full<<" known_positive="<<known_positive
             <<" new_curve_keys="<<newcurves<<" local="<<local<<" deep11="<<deep11
             <<" deep13="<<deep13<<" candidates="<<candidates<<" output="<<o.output
             <<" candidates_output="<<o.candidates<<"\n";
    return 0;
 }catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<"\n";return 2;}
}
