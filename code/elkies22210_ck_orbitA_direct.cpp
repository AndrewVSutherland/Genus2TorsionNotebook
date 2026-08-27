// Direct rational search on Orbit A of the Clebsch--Klein source-halving
// cover.  This does not enumerate the base CK surface first.
//
// Normalize the marked coordinate r1=1 and write
//
//   r_j = (n_j^2-m_j^2)/(n_j^2+m_j^2),   j=2,3,
//
// so 1-r_j^2 is automatically a rational square.  For x=r2, y=r3,
// eliminate u=r4 and v=r5 from the two CK equations.  They must satisfy
//
//   u+v = S = -1-x-y,
//   uv  = P = (x+y)(1+x)(1+y)/(1+x+y).
//
// Thus D=S^2-4P must be a rational square, after which u,v are recovered
// exactly and 1-u^2, 1-v^2 are tested for rational squareness.  By symmetry
// this is complete for cover points for which any two of t2,...,t5 have
// primitive numerator and denominator at most B.
//
// Compile:
//   c++ -O3 -std=c++17 code/elkies22210_ck_orbitA_direct.cpp \
//       -o /tmp/elkies22210_ck_orbitA_direct
// Run:
//   /tmp/elkies22210_ck_orbitA_direct 100

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

static u128 uabs128(i128 x) { return x < 0 ? (u128)(-x) : (u128)x; }

static u128 gcd128(u128 a, u128 b) {
    while (b) { u128 r=a%b; a=b; b=r; }
    return a;
}

static std::string str128(i128 x) {
    if (x==0) return "0";
    bool neg=x<0; u128 u=uabs128(x); std::string s;
    while (u) { s.push_back(char('0'+u%10)); u/=10; }
    if (neg) s.push_back('-');
    std::reverse(s.begin(),s.end()); return s;
}

struct Q {
    i128 n, d;
    Q(i128 nn=0, i128 dd=1): n(nn), d(dd) {
        if (d==0) { std::cerr << "zero denominator\n"; std::abort(); }
        if (d<0) { n=-n; d=-d; }
        u128 g=gcd128(uabs128(n),(u128)d);
        n/=i128(g); d/=i128(g);
    }
};

static Q add(const Q&a,const Q&b) {
    u128 g=gcd128((u128)a.d,(u128)b.d);
    i128 ad=a.d/i128(g), bd=b.d/i128(g);
    return Q(a.n*bd+b.n*ad,ad*b.d);
}
static Q neg(const Q&a) { return Q(-a.n,a.d); }
static Q sub(const Q&a,const Q&b) { return add(a,neg(b)); }
static Q mul(const Q&a,const Q&b) {
    u128 g1=gcd128(uabs128(a.n),(u128)b.d);
    u128 g2=gcd128(uabs128(b.n),(u128)a.d);
    return Q((a.n/i128(g1))*(b.n/i128(g2)),
             (a.d/i128(g2))*(b.d/i128(g1)));
}
static Q divide(const Q&a,const Q&b) {
    if (b.n==0) { std::cerr << "division by zero\n"; std::abort(); }
    return mul(a,Q(b.d,b.n));
}
static bool equal(const Q&a,const Q&b) { return a.n==b.n && a.d==b.d; }
static Q sqr(const Q&a) { return mul(a,a); }
static std::string qstr(const Q&a) { return str128(a.n)+"/"+str128(a.d); }

static bool sq64[64],sq63[63],sq65[65];

static void init_square_residues() {
    for (unsigned z=0;z<320;++z) {
        sq64[(z*z)%64]=true;
        sq63[(z*z)%63]=true;
        sq65[(z*z)%65]=true;
    }
}

static bool square128(u128 n, u128 &root) {
    long double d=std::sqrt((long double)n);
    std::uint64_t r=(std::uint64_t)d;
    while ((u128)r*r>n) --r;
    while ((u128)(r+1)*(r+1)<=n) ++r;
    root=r; return (u128)r*r==n;
}

static bool square_q(const Q&a,Q &root) {
    if (a.n<=0) return false; // all tested factors are nonzero on the open
    u128 un=(u128)a.n, ud=(u128)a.d;
    if (!sq64[un%64] || !sq64[ud%64]
        || !sq63[un%63] || !sq63[ud%63]
        || !sq65[un%65] || !sq65[ud%65]) return false;
    u128 rn,rd;
    if (!square128(un,rn) || !square128(ud,rd)) return false;
    root=Q(i128(rn),i128(rd)); return true;
}

struct CircleValue { i64 m,n; Q x; };

int main(int argc,char **argv) {
    int B=argc>1 ? std::atoi(argv[1]) : 100;
    if (B<2 || B>500) {
        std::cerr << "This exact __int128 implementation requires 2 <= B <= 500\n";
        return 2;
    }
    init_square_residues();
    std::vector<CircleValue> vals;
    std::set<std::pair<i64,i64>> seen;
    for (i64 n=1;n<=B;++n) for (i64 m=1;m<=B;++m) {
        if (m==n || std::gcd(m,n)!=1) continue; // x=0 or duplicate t
        Q x((i128)n*n-(i128)m*m,(i128)n*n+(i128)m*m);
        // x=+-1 cannot occur here; keep one parameter for each x value.
        std::pair<i64,i64> key={(i64)x.n,(i64)x.d};
        if (seen.insert(key).second) vals.push_back({m,n,x});
    }

    std::uint64_t pairs=0, distinct_input_squares=0, nonzero_den=0,
                  disc_positive=0, disc_square=0, uv_nonzero=0,
                  u_circle=0, v_circle=0, smooth=0;
    std::vector<std::string> hits;
    std::vector<std::string> rational_ck_completions;
    const Q one(1), two(2), four(4);
    for (std::size_t ix=0;ix<vals.size();++ix) {
        const Q &x=vals[ix].x;
        for (std::size_t iy=ix+1;iy<vals.size();++iy) {
            ++pairs;
            const Q &y=vals[iy].x;
            if (equal(sqr(x),sqr(y))) continue;
            ++distinct_input_squares;
            Q z=add(x,y), den=add(one,z);
            if (den.n==0) continue;
            ++nonzero_den;
            Q S=neg(den);
            Q P=divide(mul(mul(z,add(one,x)),add(one,y)),den);
            Q D=sub(sqr(S),mul(four,P)), sqrtD;
            if (D.n<=0) continue;
            ++disc_positive;
            if (!square_q(D,sqrtD)) continue;
            ++disc_square;
            Q u=divide(add(S,sqrtD),two);
            Q v=divide(sub(S,sqrtD),two);
            if (u.n==0 || v.n==0) continue;
            ++uv_nonzero;
            Q zu,zv;
            bool cu=square_q(sub(one,sqr(u)),zu);
            if (cu) ++u_circle;
            bool cv=square_q(sub(one,sqr(v)),zv);
            if (cv) ++v_circle;
            if (rational_ck_completions.size()<20) {
                std::ostringstream o;
                o << "CK_COMPLETION t2=" << vals[ix].m << "/" << vals[ix].n
                  << " t3=" << vals[iy].m << "/" << vals[iy].n
                  << " r=[1," << qstr(x) << "," << qstr(y) << ","
                  << qstr(u) << "," << qstr(v) << "]"
                  << " remaining_circle=" << cu << "," << cv;
                rational_ck_completions.push_back(o.str());
            }
            if (!(cu&&cv)) continue;
            std::array<Q,5> rr={one,x,y,u,v};
            bool ok=true;
            for (int i=0;i<5;++i) for (int j=i+1;j<5;++j)
                if (equal(sqr(rr[i]),sqr(rr[j]))) ok=false;
            if (!ok) continue;
            ++smooth;

            // Independent exact checks of both CK equations.
            Q p1(0),p3(0);
            for (const Q&r:rr) { p1=add(p1,r); p3=add(p3,mul(sqr(r),r)); }
            if (p1.n!=0 || p3.n!=0) {
                std::cerr << "internal CK identity failure\n"; return 3;
            }
            std::ostringstream o;
            o << "HIT t2=" << vals[ix].m << "/" << vals[ix].n
              << " t3=" << vals[iy].m << "/" << vals[iy].n
              << " r=[1," << qstr(x) << "," << qstr(y) << ","
              << qstr(u) << "," << qstr(v) << "]";
            hits.push_back(o.str());
        }
    }

    std::cout << "ELKIES22210_CK_ORBITA_DIRECT\n";
    std::cout << "parameter_height " << B << "\n";
    std::cout << "circle_values " << vals.size() << "\n";
    std::cout << "pairs " << pairs << "\n";
    std::cout << "distinct_input_squares " << distinct_input_squares << "\n";
    std::cout << "nonzero_elimination_denominator " << nonzero_den << "\n";
    std::cout << "disc_positive " << disc_positive << "\n";
    std::cout << "disc_square " << disc_square << "\n";
    std::cout << "uv_nonzero " << uv_nonzero << "\n";
    std::cout << "u_circle " << u_circle << " v_circle " << v_circle << "\n";
    for (const auto&s:rational_ck_completions) std::cout << s << "\n";
    std::cout << "smooth_cover_points " << smooth << "\n";
    std::cout << "hits " << hits.size() << "\n";
    for (const auto&s:hits) std::cout << s << "\n";
    std::cout << "DONE\n";
    return 0;
}
