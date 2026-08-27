// Resumable integer/fiber generator for new points on the full
// four-radicand A(2,2,2,8) cover, targeted at [2,2,2,24].
//
// This is the genuinely new 3D/integer lane.  It does not enumerate the
// known Filip/Adam P1 families.  Starting with every primitive three-coordinate
// fiber occurring in tor2228.txt, it extends the fourth coordinate far beyond
// the old B=16384 bank.  For fixed signed (a,b,c), abcd being a square forces
//
//     d = sign(abc) * sf(|abc|) * k^2.
//
// We enumerate these d, apply the corrected necessary local masks first, and
// only then test the other three integer square conditions.  At p=11,13 the
// allowed set is exactly branch-boundary reduction; at p=17,19,23,29,31 it is
// boundary union the corrected complete smooth 3-contact mask.

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
#include <unordered_set>
#include <utility>
#include <vector>

using i64 = long long;
using i128 = __int128_t;
using u128 = __uint128_t;

struct Options {
    std::string bank = "data/tor2228_bank.txt";
    std::string mask_prefix =
        "results/target_22224_full_family_halving_corrected_complete";
    std::string output =
        "results/target_22224_a2228_deep13_generator_integer_candidates.tsv";
    i64 d_min = 16385;
    i64 d_max = 10000000000LL;
    size_t fiber_start = 0;
    size_t fiber_stop = (size_t)-1;
    bool all_signs = true;
    size_t progress = 100;
    bool surface = false;
    bool rank8_charts = false;
    bool deep_p5_charts = false;
    std::string deep11 = "results/target_22224_direct_contact_deep13_padic_tangent_p11.tsv";
    std::string deep13 = "results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv";
    i64 t_num = 1;
    i64 t_den = 2;
    i64 r_min = -10000000;
    i64 r_max = 10000000;
};

static std::vector<std::string> split(const std::string &s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string z;
    while (std::getline(ss, z, sep)) out.push_back(z);
    return out;
}

static i64 gcd4(i64 a, i64 b, i64 c, i64 d) {
    return std::gcd(std::gcd(std::llabs(a), std::llabs(b)),
                    std::gcd(std::llabs(c), std::llabs(d)));
}

static std::array<i64,4> primitive(std::array<i64,4> v) {
    i64 g = gcd4(v[0], v[1], v[2], v[3]);
    if (g > 1) for (auto &z : v) z /= g;
    for (i64 z : v) {
        if (!z) continue;
        if (z < 0) for (auto &w : v) w = -w;
        break;
    }
    return v;
}

static std::string curve_key(std::array<i64,4> v) {
    v = primitive(v);
    for (auto &z : v) z = std::llabs(z);
    std::sort(v.begin(), v.end());
    return std::to_string(v[0]) + "," + std::to_string(v[1]) + "," +
           std::to_string(v[2]) + "," + std::to_string(v[3]);
}

static std::vector<std::array<i64,4>> read_bank(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open bank " + path);
    std::set<std::array<i64,4>> seen;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 3 || line.front() != '[' || line.back() != ']') continue;
        auto fields = split(line.substr(1, line.size()-2), ',');
        if (fields.size() != 4) continue;
        std::array<i64,4> v{};
        for (int i=0;i<4;i++) v[i] = std::stoll(fields[i]);
        seen.insert(primitive(v));
    }
    return std::vector<std::array<i64,4>>(seen.begin(), seen.end());
}

static u128 isqrt128(u128 n) {
    if (!n) return 0;
    u128 r = (u128)std::sqrt((long double)n);
    while ((r+1) <= n/(r+1)) ++r;
    while (r > n/r) --r;
    return r;
}

static bool square128(i128 z) {
    if (z < 0) return false;
    u128 n = (u128)z;
    u128 r = isqrt128(n);
    return r*r == n;
}

static i64 floor_sqrt64(i64 n) {
    i64 r = (i64)std::sqrt((long double)n);
    while ((r+1) <= n/(r+1)) ++r;
    while (r > n/r) --r;
    return r;
}

static i64 ceil_sqrt64(i64 n) {
    i64 r = floor_sqrt64(n);
    return r*r == n ? r : r+1;
}

static std::array<i128,4> radicands(const std::array<i64,4> &v) {
    i128 a=v[0], b=v[1], c=v[2], d=v[3];
    return {a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)};
}

static bool full_cover(const std::array<i64,4> &v) {
    auto rr = radicands(v);
    return square128(rr[0]) && square128(rr[1]) &&
           square128(rr[2]) && square128(rr[3]);
}

static bool smooth_tuple(const std::array<i64,4> &v) {
    for (i64 z : v) if (z == 0) return false;
    for (int i=0;i<4;i++) for (int j=i+1;j<4;j++)
        if (std::llabs(v[i]) == std::llabs(v[j])) return false;
    return true;
}

static bool rectangle_tuple(const std::array<i64,4> &v) {
    i128 a=std::llabs(v[0]),b=std::llabs(v[1]);
    i128 c=std::llabs(v[2]),d=std::llabs(v[3]);
    return a*b==c*d || a*c==b*d || a*d==b*c;
}

struct LocalMask {
    int p;
    std::vector<unsigned char> allowed; // indexed by (((a*p+b)*p+c)*p+d)

    size_t index(int a,int b,int c,int d) const {
        return (((size_t)a*p+b)*p+c)*p+d;
    }
    bool get(i64 a,i64 b,i64 c,i64 d) const {
        auto mod = [this](i64 z) { int r=(int)(z%p); return r<0?r+p:r; };
        return allowed[index(mod(a),mod(b),mod(c),mod(d))] != 0;
    }
};

static bool boundary4(int a,int b,int c,int d,int p) {
    int s[5] = {0, a*a%p, b*b%p, c*c%p, d*d%p};
    for (int i=0;i<5;i++) for (int j=i+1;j<5;j++) if (s[i]==s[j]) return true;
    return false;
}

static int invmod(int a,int p) {
    int r=1,e=p-2;
    while(e){ if(e&1)r=(int)((long long)r*a%p); a=(int)((long long)a*a%p); e>>=1; }
    return r;
}

static LocalMask make_mask(int p,const std::string &prefix) {
    LocalMask M{p,std::vector<unsigned char>((size_t)p*p*p*p,0)};
    for(int a=0;a<p;a++) for(int b=0;b<p;b++)
      for(int c=0;c<p;c++) for(int d=0;d<p;d++)
        if(boundary4(a,b,c,d,p)) M.allowed[M.index(a,b,c,d)]=1;
    if (p==11 || p==13) return M;
    std::string path = prefix + "_p" + std::to_string(p) + ".tsv";
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open corrected mask " + path);
    std::string line; std::getline(in,line);
    while(std::getline(in,line)) {
        if(line.empty()) continue;
        auto f=split(line,'\t');
        if(f.size()<4) continue;
        int v[4]; for(int i=0;i<4;i++)v[i]=std::stoi(f[i])%p;
        for(int lam=1;lam<p;lam++) {
            int a=lam*v[0]%p,b=lam*v[1]%p,c=lam*v[2]%p,d=lam*v[3]%p;
            M.allowed[M.index(a,b,c,d)]=1;
        }
    }
    return M;
}

static LocalMask make_projective_permutation_mask(
    int p, const std::array<int,4> &target) {
    LocalMask M{p,std::vector<unsigned char>((size_t)p*p*p*p,0)};
    std::array<int,4> perm=target;
    std::sort(perm.begin(),perm.end());
    do {
        for(int signs=0;signs<16;signs++)
          for(int lam=1;lam<p;lam++) {
            int z[4];
            for(int j=0;j<4;j++) {
                z[j]=lam*perm[j]%p;
                if((signs>>j)&1)z[j]=z[j]?p-z[j]:0;
            }
            M.allowed[M.index(z[0],z[1],z[2],z[3])]=1;
          }
    } while(std::next_permutation(perm.begin(),perm.end()));
    return M;
}

static bool raw_chart_prefilter(
    const LocalMask &M, i64 a, i64 b, i64 c, i64 d) {
    // If all coordinates vanish modulo p, primitive normalization may divide
    // out a p before the exact chart check.  Retain that rare case here.
    auto mod=[&](i64 z){int r=(int)(z%M.p);return r<0?r+M.p:r;};
    if(mod(a)==0 && mod(b)==0 && mod(c)==0 && mod(d)==0)return true;
    return M.get(a,b,c,d);
}

static i64 modll(i64 a,i64 m) { a%=m; return a<0?a+m:a; }

static i64 unit_inverse(i64 a,i64 m) {
    i64 old_r=a,r=m,old_s=1,s=0;
    while(r) {
        i64 q=old_r/r;
        i64 nr=old_r-q*r;old_r=r;r=nr;
        i64 ns=old_s-q*s;old_s=s;s=ns;
    }
    if(std::llabs(old_r)!=1)throw std::runtime_error("nonunit projective pivot");
    if(old_r<0)old_s=-old_s;
    return modll(old_s,m);
}

static std::array<i64,4> deep_projective_key(
    const std::array<i64,4> &v,i64 modulus,int p) {
    std::array<i64,4> best={modulus,modulus,modulus,modulus};
    bool found=false;
    for(int pivot=0;pivot<4;pivot++) {
        i64 x=modll(v[pivot],modulus);
        if(x%p==0)continue;
        i64 q=unit_inverse(x,modulus);
        std::array<i64,4> z{};
        for(int j=0;j<4;j++) {
            i64 r=(i64)((i128)modll(v[j],modulus)*q%modulus);
            z[j]=std::min(r,r?modulus-r:0);
        }
        std::sort(z.begin(),z.end());
        if(!found || z<best){best=z;found=true;}
    }
    if(!found)return {modulus,modulus,modulus,modulus};
    return best;
}

struct DeepChartBank {
    int p=0,exponent=0;
    i64 modulus=0;
    size_t rows=0;
    std::set<std::array<i64,4>> keys;

    bool get(const std::array<i64,4> &v) const {
        return keys.count(deep_projective_key(v,modulus,p))!=0;
    }
};

static DeepChartBank read_deep_charts(
    const std::string &path,int p,int exponent) {
    DeepChartBank B;B.p=p;B.exponent=exponent;B.modulus=1;
    for(int i=0;i<exponent;i++)B.modulus*=p;
    std::ifstream in(path);
    if(!in)throw std::runtime_error("cannot open deep chart bank "+path);
    std::string line;
    if(!std::getline(in,line))throw std::runtime_error("empty deep chart bank "+path);
    auto header=split(line,'\t');
    int col[4]={-1,-1,-1,-1};
    const std::string names[4]={"a","b","c","d"};
    for(int j=0;j<4;j++)for(size_t i=0;i<header.size();i++)
        if(header[i]==names[j]){col[j]=(int)i;break;}
    for(int j=0;j<4;j++)if(col[j]<0)
        throw std::runtime_error("missing base column in "+path);
    while(std::getline(in,line)) {
        if(line.empty())continue;
        auto f=split(line,'\t');
        std::array<i64,4> v{};bool ok=true;
        for(int j=0;j<4;j++) {
            if((size_t)col[j]>=f.size()){ok=false;break;}
            v[j]=std::stoll(f[col[j]]);
        }
        if(!ok)continue;
        ++B.rows;B.keys.insert(deep_projective_key(v,B.modulus,p));
    }
    if(B.keys.empty())throw std::runtime_error("no deep chart keys in "+path);
    return B;
}

static Options parse_options(int argc,char **argv) {
    Options o;
    for(int i=1;i<argc;i++) {
        std::string a=argv[i];
        auto need=[&](){if(i+1>=argc)throw std::runtime_error("missing value for "+a);return std::string(argv[++i]);};
        if(a=="--bank")o.bank=need();
        else if(a=="--mask-prefix")o.mask_prefix=need();
        else if(a=="--output")o.output=need();
        else if(a=="--d-min")o.d_min=std::stoll(need());
        else if(a=="--d-max")o.d_max=std::stoll(need());
        else if(a=="--fiber-start")o.fiber_start=std::stoull(need());
        else if(a=="--fiber-stop")o.fiber_stop=std::stoull(need());
        else if(a=="--progress")o.progress=std::stoull(need());
        else if(a=="--positive-only")o.all_signs=false;
        else if(a=="--surface")o.surface=true;
        else if(a=="--rank8-charts")o.rank8_charts=true;
        else if(a=="--deep-p5-charts")o.deep_p5_charts=true;
        else if(a=="--deep11")o.deep11=need();
        else if(a=="--deep13")o.deep13=need();
        else if(a=="--t-num")o.t_num=std::stoll(need());
        else if(a=="--t-den")o.t_den=std::stoll(need());
        else if(a=="--r-min")o.r_min=std::stoll(need());
        else if(a=="--r-max")o.r_max=std::stoll(need());
        else throw std::runtime_error("unknown option "+a);
    }
    return o;
}

int main(int argc,char **argv) {
  try {
    Options opt=parse_options(argc,argv);
    if(opt.deep_p5_charts)opt.rank8_charts=true;
    auto bank=read_bank(opt.bank);
    std::unordered_set<std::string> known;
    int maxcoord=0;
    for(auto v:bank){known.insert(curve_key(v));for(auto z:v)maxcoord=std::max(maxcoord,(int)std::llabs(z));}

    std::set<std::array<i64,3>> triple_set;
    size_t off_rectangle_bank=0;
    for(auto v:bank) {
      if(rectangle_tuple(v))continue;
      ++off_rectangle_bank;
      for(int omit=0;omit<4;omit++) {
        std::array<i64,3> t{};int q=0;
        for(int i=0;i<4;i++)if(i!=omit)t[q++]=std::llabs(v[i]);
        std::sort(t.begin(),t.end()); triple_set.insert(t);
      }
    }
    std::vector<std::array<i64,3>> triples;
    // Corrected high-survival off-rectangle charts first.
    const std::vector<std::array<i64,4>> priority_seeds = {
        {18,686,4932,8631}, {833,2529,4496,4913},
        {17,4208,4471,5329}, {25,264,936,6864}
    };
    std::set<std::array<i64,3>> priority_triples;
    for(auto v:priority_seeds) for(int omit=0;omit<4;omit++) {
        std::array<i64,3> t{};int q=0;
        for(int i=0;i<4;i++)if(i!=omit)t[q++]=v[i];
        std::sort(t.begin(),t.end());
        if(triple_set.count(t))priority_triples.insert(t);
    }
    for(auto t:priority_triples){triples.push_back(t);triple_set.erase(t);}
    triples.insert(triples.end(),triple_set.begin(),triple_set.end());
    size_t stop=std::min(opt.fiber_stop,triples.size());
    if(opt.fiber_start>stop)opt.fiber_start=stop;

    std::vector<int> sf(maxcoord+1);
    for(int n=1;n<=maxcoord;n++) {
        int x=n,k=1;
        for(int p=2;(long long)p*p<=x;p+=(p==2?1:2)) {
            int parity=0; while(x%p==0){x/=p;parity^=1;} if(parity)k*=p;
        }
        if(x>1)k*=x; sf[n]=k;
    }
    auto combine=[](i64 x,i64 y){i64 g=std::gcd(x,y);return (x/g)*(y/g);};

    std::vector<LocalMask> masks;
    for(int p:{11,13,17,19,23,29,31})masks.push_back(make_mask(p,opt.mask_prefix));
    // Corrected deep-lift charts supplied by the p^5 incidence audit.
    // We take their full projective permutation orbits independently at the
    // two primes, since the four branch labels are intrinsically unordered.
    const LocalMask chart11=make_projective_permutation_mask(11,{0,5,2,8});
    const LocalMask chart13=make_projective_permutation_mask(13,{7,7,8,5});
    DeepChartBank deep11,deep13;
    if(opt.deep_p5_charts) {
        deep11=read_deep_charts(opt.deep11,11,5);
        deep13=read_deep_charts(opt.deep13,13,5);
        std::cout<<"DEEP_P5_MODEL p11_modulus="<<deep11.modulus
                 <<" rows="<<deep11.rows<<" projective_keys="<<deep11.keys.size()
                 <<" path="<<opt.deep11
                 <<" p13_modulus="<<deep13.modulus
                 <<" rows="<<deep13.rows<<" projective_keys="<<deep13.keys.size()
                 <<" path="<<opt.deep13
                 <<" orbit=independent_unit_scaling_signed_permutations\n";
    }

    // Genuine rational surface supplied by the corrected full-cover audit:
    //   (a,b,c,d)=(1,r,s,rs),
    //   q=2t/(1-t^2), s=(1-q^2 r)/(q^2-r).
    // For t=m/n put A=(n^2-m^2)^2, B=(2mn)^2.  At integral r a common
    // projective multiple is (B-Ar, r(B-Ar), A-Br, r(A-Br)).
    if(opt.surface) {
        i128 mm=opt.t_num,nn=opt.t_den;
        i128 A0=(nn*nn-mm*mm)*(nn*nn-mm*mm);
        i128 B0=(2*mm*nn)*(2*mm*nn);
        if(A0>INT64_MAX || B0>INT64_MAX)throw std::runtime_error("surface t too large");
        i64 A=(i64)A0,B=(i64)B0;
        std::vector<std::pair<int,std::vector<unsigned char>>> rres;
        for(auto &M:masks) {
            std::vector<unsigned char> ok(M.p,0);
            for(int r=0;r<M.p;r++) {
                i64 sd=B-A*r,sn=A-B*r;
                if(M.get(sd,(i64)r*sd,sn,(i64)r*sn))ok[r]=1;
            }
            int cnt=std::accumulate(ok.begin(),ok.end(),0);
            rres.push_back({cnt,std::move(ok)});
        }
        std::sort(rres.begin(),rres.end(),[](auto &x,auto &y){return x.first<y.first;});
        std::ofstream out(opt.output);
        if(!out)throw std::runtime_error("cannot write "+opt.output);
        out<<"t_num\tt_den\tr\ta\tb\tc\td\tcurve_key\tknown_bank\n";
        unsigned long long tested=0,local=0,cover=0,news=0;
        std::unordered_set<std::string> emitted;
        std::cout<<"A2228_SURFACE_GENERATOR_START t="<<opt.t_num<<"/"<<opt.t_den
                 <<" r=["<<opt.r_min<<","<<opt.r_max<<"] profiles=";
        for(auto &z:rres)std::cout<<z.second.size()<<":"<<z.first<<",";
        std::cout<<"\n";
        for(i64 r=opt.r_min;r<=opt.r_max;r++) {
            ++tested; bool ok=true;
            for(auto &rec:rres){int p=(int)rec.second.size();int rr=(int)(r%p);if(rr<0)rr+=p;if(!rec.second[rr]){ok=false;break;}}
            if(!ok)continue;
            ++local;
            i128 sd0=(i128)B-(i128)A*r,sn0=(i128)A-(i128)B*r;
            i128 vals0[4]={sd0,(i128)r*sd0,sn0,(i128)r*sn0};
            bool fits=true;for(auto z:vals0)if(z<INT64_MIN||z>INT64_MAX)fits=false;
            if(!fits)continue;
            std::array<i64,4> v={(i64)vals0[0],(i64)vals0[1],(i64)vals0[2],(i64)vals0[3]};
            v=primitive(v);
            if(!smooth_tuple(v) || rectangle_tuple(v))continue;
            bool final_local=true;for(auto &M:masks)if(!M.get(v[0],v[1],v[2],v[3])){final_local=false;break;}
            if(!final_local)continue;
            if(!full_cover(v))throw std::runtime_error("surface full-cover identity failed");
            ++cover;
            std::string key=curve_key(v);bool old=known.count(key);
            if(!emitted.insert(key).second)continue;
            if(!old)++news;
            out<<opt.t_num<<'\t'<<opt.t_den<<'\t'<<r<<'\t'<<v[0]<<'\t'<<v[1]
               <<'\t'<<v[2]<<'\t'<<v[3]<<'\t'<<key<<'\t'<<(old?1:0)<<'\n';
            out.flush();
            std::cout<<"SURFACE_HIT r="<<r<<" tuple=["<<v[0]<<","<<v[1]<<","<<v[2]<<","<<v[3]
                     <<"] old="<<old<<"\n";
        }
        std::cout<<"A2228_SURFACE_GENERATOR_DONE tested="<<tested<<" local="<<local
                 <<" cover="<<cover<<" new="<<news<<" output="<<opt.output<<"\n";
        return 0;
    }

    std::ofstream out(opt.output);
    if(!out)throw std::runtime_error("cannot write "+opt.output);
    out<<"fiber\tsigns\tk\ta\tb\tc\td\tprimitive_a\tprimitive_b\tprimitive_c\tprimitive_d\tcurve_key\tsquare_mask\tfull_cover\tknown_bank\n";

    unsigned long long k_total=0,k_local=0,deep_tests=0,deep_survivors=0;
    unsigned long long exact_tests=0,cover_hits=0,new_hits=0;
    std::unordered_set<std::string> emitted;
    std::vector<std::array<int,3>> signs={{{1,1,1}}};
    if(opt.all_signs)signs={{{1,1,1}},{{1,1,-1}},{{1,-1,1}},{{1,-1,-1}}};

    std::cout<<"A2228_INTEGER_GENERATOR_START bank="<<bank.size()
             <<" off_rectangle_bank="<<off_rectangle_bank
             <<" fibers="<<triples.size()<<" range=["<<opt.fiber_start<<","<<stop
             <<") d=["<<opt.d_min<<","<<opt.d_max<<"] signs="<<signs.size()<<"\n";
    for(size_t fi=opt.fiber_start;fi<stop;fi++) {
      auto mag=triples[fi];
      i64 s=combine(combine(sf[mag[0]],sf[mag[1]]),sf[mag[2]]);
      if(s<=0 || s>opt.d_max)continue;
      i64 k0=ceil_sqrt64((opt.d_min+s-1)/s);
      i64 k1=floor_sqrt64(opt.d_max/s);
      if(k0<1)k0=1;
      for(size_t si=0;si<signs.size();si++) {
        i64 a=mag[0]*signs[si][0],b=mag[1]*signs[si][1],c=mag[2]*signs[si][2];
        int dsign=((a<0)^(b<0)^(c<0))?-1:1;
        std::vector<std::pair<int,std::vector<unsigned char>>> kres;
        for(auto &M:masks) {
            std::vector<unsigned char> ok(M.p,0);
            for(int r=0;r<M.p;r++) {
                i64 d=(i64)dsign*s*r*r;
                bool good=M.get(a,b,c,d);
                if(opt.rank8_charts && M.p==11)
                    good=raw_chart_prefilter(chart11,a,b,c,d);
                if(opt.rank8_charts && M.p==13)
                    good=raw_chart_prefilter(chart13,a,b,c,d);
                if(good)ok[r]=1;
            }
            int cnt=std::accumulate(ok.begin(),ok.end(),0);
            kres.push_back({cnt,std::move(ok)});
        }
        std::sort(kres.begin(),kres.end(),[](auto &x,auto &y){return x.first<y.first;});
        auto test_k = [&](i64 k) {
            ++k_total;
            bool lok=true;
            // kres lost p association; recover it from vector size.
            for(auto &rec:kres){int p=(int)rec.second.size();if(!rec.second[k%p]){lok=false;break;}}
            if(!lok)return;
            ++k_local;
            i64 d=dsign*s*k*k;
            std::array<i64,4> raw={a,b,c,d};
            auto v=primitive(raw);
            if(!smooth_tuple(v) || rectangle_tuple(v))return;
            bool final_local=true;
            for(auto &M:masks) {
                bool good=M.get(v[0],v[1],v[2],v[3]);
                if(opt.rank8_charts && M.p==11)
                    good=chart11.get(v[0],v[1],v[2],v[3]);
                if(opt.rank8_charts && M.p==13)
                    good=chart13.get(v[0],v[1],v[2],v[3]);
                if(!good){final_local=false;break;}
            }
            if(!final_local)return;
            if(opt.deep_p5_charts) {
                ++deep_tests;
                if(!deep11.get(v) || !deep13.get(v))return;
                ++deep_survivors;
            }
            ++exact_tests;
            auto rr=radicands(v);
            int square_mask=0;
            for(int j=0;j<4;j++)if(square128(rr[j]))square_mask|=1<<j;
            bool is_cover=square_mask==15;
            std::string key=curve_key(v);
            bool old=known.count(key);
            if(opt.rank8_charts || is_cover) {
                out<<fi<<'\t'<<si<<'\t'<<k<<'\t'<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d
                   <<'\t'<<v[0]<<'\t'<<v[1]<<'\t'<<v[2]<<'\t'<<v[3]
                   <<'\t'<<key<<'\t'<<square_mask<<'\t'<<(is_cover?1:0)
                   <<'\t'<<(old?1:0)<<'\n';
                out.flush();
            }
            if(!is_cover)return;
            ++cover_hits;
            if(!emitted.insert(key).second)return;
            if(!old)++new_hits;
            std::cout<<"FULL_COVER_HIT fiber="<<fi<<" signs="<<si<<" k="<<k
                     <<" tuple=["<<v[0]<<","<<v[1]<<","<<v[2]<<","<<v[3]
                     <<"] old="<<old<<"\n";
        };
        if(opt.rank8_charts) {
            // The two deep charts reduce k to a tiny set modulo 11*13.
            // Iterate only those CRT classes instead of scanning every k.
            std::vector<int> crt;
            for(int r=0;r<143;r++) {
                i64 d=(i64)dsign*s*r*r;
                if(raw_chart_prefilter(chart11,a,b,c,d) &&
                   raw_chart_prefilter(chart13,a,b,c,d))crt.push_back(r);
            }
            for(int r:crt) {
                i64 rem=(r-(k0%143)+143)%143;
                for(i64 k=k0+rem;k<=k1;k+=143)test_k(k);
            }
        } else {
            for(i64 k=k0;k<=k1;k++)test_k(k);
        }
      }
      if(opt.progress && (fi+1)%opt.progress==0)
        std::cout<<"PROGRESS fiber="<<(fi+1)<<" k_total="<<k_total
                 <<" local="<<k_local<<" deep_tests="<<deep_tests
                 <<" deep_survivors="<<deep_survivors<<" exact="<<exact_tests
                 <<" cover="<<cover_hits<<" new="<<new_hits<<"\n";
    }
    std::cout<<"A2228_INTEGER_GENERATOR_DONE fibers="<<(stop-opt.fiber_start)
             <<" k_total="<<k_total<<" local="<<k_local<<" exact="<<exact_tests
             <<" deep_tests="<<deep_tests<<" deep_survivors="<<deep_survivors
             <<" cover_hits="<<cover_hits<<" new_hits="<<new_hits
             <<" output="<<opt.output<<"\n";
    return 0;
  } catch(const std::exception &e) {
    std::cerr<<"ERROR "<<e.what()<<"\n"; return 2;
  }
}
