// Exhaustive bounded CRT/reconstruction lane for the fresh direct-contact chart.
//
// We cross every projected p^2 boundary class at p=11,13 with every generic
// open target class at p=17,19,23,29,31.  There are exactly
//
//     600 * 650 * 1 * 2 * 4 * 3 * 8 = 74,880,000
//
// combinations.  Coordinatewise classical rational reconstruction is performed
// at the full modulus.  The p=37,41,43 masks are then used as sound lookup
// filters: if every denominator is a unit at p, the reduced triple must occur
// in the union of the fresh open and boundary masks; otherwise the candidate is
// retained because this affine finite-prime filter does not apply.
//
// Compile, for example, with
//   c++ -O3 -std=c++17 -o /tmp/direct_exhaustive <this-file>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using I = std::int64_t;
using Triple = std::array<I,3>;

static const std::string ROOT = "results/";

static std::vector<std::string> split_tab(const std::string &s) {
    std::vector<std::string> z; std::stringstream ss(s); std::string x;
    while (std::getline(ss,x,'\t')) z.push_back(x);
    return z;
}

static std::vector<Triple> read_triples(const std::string &path, int offset) {
    std::ifstream f(path); if (!f) { std::cerr << "cannot open " << path << "\n"; std::exit(2); }
    std::string line; std::getline(f,line); std::set<Triple> ans;
    while (std::getline(f,line)) {
        if (line.empty()) continue; auto z=split_tab(line);
        if ((int)z.size()<offset+3) continue;
        ans.insert({std::stoll(z[offset]),std::stoll(z[offset+1]),std::stoll(z[offset+2])});
    }
    return std::vector<Triple>(ans.begin(),ans.end());
}

static I modnorm(__int128 a, I m) {
    I z=(I)(a%m); return z<0 ? z+m : z;
}

static I egcd(I a,I b,I &x,I &y) {
    if (!b) {x=1;y=0;return a;}
    I u,v,g=egcd(b,a%b,u,v);x=v;y=u-(a/b)*v;return g;
}

static I invmod(I a,I m) {
    I x,y,g=egcd(modnorm(a,m),m,x,y); if (g!=1) {std::cerr<<"nonunit inverse\n";std::exit(2);}
    return modnorm(x,m);
}

static I crt(I a,I m,I b,I n) {
    I k=modnorm((__int128)(b-a)*invmod(m,n),n);
    return (I)((__int128)a+(__int128)m*k);
}

static Triple crt_triple(const Triple&a,I m,const Triple&b,I n) {
    return {crt(a[0],m,b[0],n),crt(a[1],m,b[1],n),crt(a[2],m,b[2],n)};
}

struct Rat { I n,d; };

static bool ratrec(I a,I m,I bound,Rat &out) {
    I r0=m,r1=modnorm(a,m),t0=0,t1=1;
    while (r1>bound) {
        I q=r0/r1;
        I nr=r0-q*r1; r0=r1;r1=nr;
        I nt=t0-q*t1; t0=t1;t1=nt;
    }
    if (!t1) return false;
    I n=r1,d=t1; if (d<0) {n=-n;d=-d;}
    I g=std::gcd(n<0?-n:n,d);n/=g;d/=g;
    if ((n<0?-n:n)>bound || d>bound || std::gcd(d,m)!=1) return false;
    if (modnorm((__int128)n-(__int128)a*d,m)!=0) return false;
    out={n,d};return true;
}

static I triple_key(const Triple &a,I p) { return (a[0]*p+a[1])*p+a[2]; }

static std::unordered_set<I> allowed_at(I p) {
    std::set<Triple> z;
    auto b=read_triples(ROOT+"target_22224_direct_contact_deep13_boundary_p"+std::to_string(p)+".tsv",1);
    auto t=read_triples(ROOT+"target_22224_direct_contact_deep13_p"+std::to_string(p)+".tsv",1);
    z.insert(b.begin(),b.end());z.insert(t.begin(),t.end());
    std::unordered_set<I> ans;ans.reserve(2*z.size());for (auto a:z) ans.insert(triple_key(a,p));
    std::cerr << "FILTER_MASK p="<<p<<" boundary="<<b.size()<<" open="<<t.size()<<" union="<<z.size()<<"\n";
    return ans;
}

// False means a mathematically sound rejection.  A denominator divisible by p
// makes this affine chart filter inapplicable, so such a row is conservatively kept.
static bool passes_filter(const std::array<Rat,3>&q,I p,const std::unordered_set<I>&allowed,bool &nonunit) {
    Triple r; nonunit=false;
    for (int j=0;j<3;j++) {
        if (q[j].d%p==0) {nonunit=true;return true;}
        r[j]=modnorm((__int128)q[j].n*invmod(q[j].d,p),p);
    }
    return allowed.count(triple_key(r,p));
}

int main(int argc,char **argv) {
    std::string output=ROOT+"target_22224_direct_contact_deep13_exhaustive_survivors.tsv";
    std::string log=ROOT+"target_22224_direct_contact_deep13_exhaustive.log";
    std::size_t start=0,count=(std::size_t)-1;
    bool padic_tangent=false;
    bool expanded_bank=false;
    for (int i=1;i<argc;i++) {
        std::string a=argv[i];
        if (a=="--output" && i+1<argc) output=argv[++i];
        else if (a=="--log" && i+1<argc) log=argv[++i];
        else if (a=="--start-deep" && i+1<argc) start=std::stoull(argv[++i]);
        else if (a=="--count-deep" && i+1<argc) count=std::stoull(argv[++i]);
        else if (a=="--padic-tangent") padic_tangent=true;
        else if (a=="--expanded-bank") {expanded_bank=true;padic_tangent=true;}
        else {std::cerr<<"unknown/missing argument "<<a<<"\n";return 2;}
    }
    I mod11=padic_tangent?161051:121,mod13=padic_tangent?371293:169;
    std::string stem11=expanded_bank?"target_22224_direct_contact_deep13_padic_expanded_p11.tsv":
                       (padic_tangent?"target_22224_direct_contact_deep13_padic_tangent_p11.tsv":"target_22224_direct_contact_deep13_lift_p11.tsv");
    std::string stem13=expanded_bank?"target_22224_direct_contact_deep13_padic_expanded_p13.tsv":
                       (padic_tangent?"target_22224_direct_contact_deep13_padic_tangent_p13.tsv":"target_22224_direct_contact_deep13_lift_p13.tsv");
    auto m11=read_triples(ROOT+stem11,0);auto m13=read_triples(ROOT+stem13,0);
    std::vector<Triple> deep;deep.reserve(m11.size()*m13.size());
    for (auto a:m11) for(auto b:m13) deep.push_back(crt_triple(a,mod11,b,mod13));

    const std::array<I,5> ps={17,19,23,29,31};
    std::vector<Triple> good(1,Triple{0,0,0});I mgood=1;
    for (I p:ps) {
        auto z=read_triples(ROOT+"target_22224_direct_contact_deep13_p"+std::to_string(p)+".tsv",1);
        std::vector<Triple> next;next.reserve(good.size()*z.size());
        for(auto a:good) for(auto b:z) next.push_back(crt_triple(a,mgood,b,p));
        good.swap(next);mgood*=p;
    }
    auto a37=allowed_at(37),a41=allowed_at(41),a43=allowed_at(43);
    I mdeep=mod11*mod13,modulus=mdeep*mgood,bound=(I)std::sqrt((long double)(modulus/2));
    while ((__int128)(bound+1)*(bound+1)<=modulus/2) ++bound;
    while ((__int128)bound*bound>modulus/2) --bound;
    std::size_t end=std::min(deep.size(),start+std::min(count,deep.size()-std::min(start,deep.size())));
    if(start>deep.size()) start=end=deep.size();
    std::uint64_t total=0,r1=0,r2=0,r3=0,allrec=0,tzero=0,pass37=0,pass41=0,pass43=0;
    std::uint64_t nonunit37=0,nonunit41=0,nonunit43=0;
    std::ofstream out(output);out<<"u_num\tu_den\tt_num\tt_den\tw_num\tw_den\n";
    auto t0=std::chrono::steady_clock::now();
    for(std::size_t di=start;di<end;di++) for(const auto &g:good) {
        ++total;std::array<Rat,3> q;bool ok=true;
        for(int j=0;j<3;j++) {
            I a=crt(deep[di][j],mdeep,g[j],mgood);
            if(!ratrec(a,modulus,bound,q[j])) {ok=false;break;}
            if(j==0) ++r1; else if(j==1) ++r2; else ++r3;
        }
        if(!ok) continue;++allrec;if(q[1].n==0){++tzero;continue;}
        bool nu=false;if(!passes_filter(q,37,a37,nu))continue;++pass37;if(nu)++nonunit37;
        if(!passes_filter(q,41,a41,nu))continue;++pass41;if(nu)++nonunit41;
        if(!passes_filter(q,43,a43,nu))continue;++pass43;if(nu)++nonunit43;
        for(int j=0;j<3;j++) out<<q[j].n<<'\t'<<q[j].d<<(j==2?'\n':'\t');
    }
    out.close();double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::ostringstream s;
    s<<"DIRECT_CONTACT_EXHAUSTIVE\n"
     <<"MODE "<<(expanded_bank?"expanded_padic_bank":(padic_tangent?"padic_tangent":"projected_p2"))<<"\n"
     <<"MASK_SIZES p11sq "<<m11.size()<<" p13sq "<<m13.size()<<" open_product "<<good.size()<<"\n"
     <<"DEEP_RANGE "<<start<<" "<<end<<" of "<<deep.size()<<"\n"
     <<"MODULUS "<<modulus<<" RECONSTRUCTION_BOUND "<<bound<<"\n"
     <<"COMBINATIONS "<<total<<" RATREC_PREFIX_U "<<r1<<" UT "<<r2<<" UTW "<<r3
     <<" ALL_RECONSTRUCTED "<<allrec<<" TZERO "<<tzero<<"\n"
     <<"PASS_P37 "<<pass37<<" NONUNIT_P37 "<<nonunit37
     <<" PASS_P41 "<<pass41<<" NONUNIT_P41 "<<nonunit41
     <<" PASS_P43 "<<pass43<<" NONUNIT_P43 "<<nonunit43<<"\n"
     <<"SECONDS "<<sec<<" OUTPUT "<<output<<"\n";
    std::cout<<s.str();std::ofstream lf(log);lf<<s.str();
    return 0;
}
