// Fast local sieve for the two-root contact-7 surface with an added
// rational 3-torsion class.
//
// Absorb the two signs in the root conditions into signed parameters u,v:
//
//   r = 1-u^2,  w = 1-v^2,  h(r)=u^7,  h(w)=v^7,
//   h = 1-(7/2)x+a*x^2+b*x^3.
//
// In cancelled coordinates put
//
//   A(z)=(2z^5+4z^4+6z^3+8z^2+10z+5)/(2(z+1)^2),
//   b=(A(v)-A(u))/(u^2-v^2),  a=A(u)-b(1-u^2).
//
// This extends across u=1 and v=1; those are removable roots-at-zero
// branches.  Its genuine parameter poles are u=-1, v=-1, and u=-v.
// Away from those poles and the singular u=0, v=0, u=v loci, the resulting
// genus-2 curve
//
//   y^2 = f = (h^2+(x-1)^7)/x^2
//
// has visible torsion Z/2 x Z/14, of order 28.  If it also has a rational
// 3-class, then the prime-to-p part of 84 divides #J(F_p) at every good
// prime p.  For every affine residue chart this program distinguishes:
//
//   unknown : a chart denominator vanishes, the contact point degenerates,
//             or the displayed integral model has bad reduction;
//   pass    : good reduction and the required divisibility holds;
//   fail    : good reduction and the required divisibility fails.
//
// Only fail residues are discarded.  Thus the search treats all bad and
// boundary residue disks conservatively.
//
// The finite-field Jacobian order is computed independently from #C(F_p)
// and #C(F_{p^2}).  No Magma tables are embedded.
//
// Build/run:
//   c++ -O3 -std=c++17 code/contact7_two_root_plus3_sieve.cpp \
//       -o /tmp/contact7_two_root_plus3_sieve
//   /tmp/contact7_two_root_plus3_sieve 100 \
//       data/contact7_two_root_plus3_h100.txt \
//       data/contact7_two_root_plus3_h100_survivors.m
//   /tmp/contact7_two_root_plus3_sieve 10000 \
//       data/contact7_two_root_plus3_u1_h10000.txt \
//       data/contact7_two_root_plus3_u1_h10000_survivors.m branch1

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

static int mod(int64_t a, int p) {
    int r = int(a % p);
    return r < 0 ? r + p : r;
}
static int add(int a, int b, int p) { return mod(int64_t(a) + b, p); }
static int sub(int a, int b, int p) { return mod(int64_t(a) - b, p); }
static int mul(int a, int b, int p) { return int((int64_t(a) * b) % p); }
static int power(int a, int64_t n, int p) {
    int ans = 1;
    while (n) {
        if (n & 1) ans = mul(ans, a, p);
        a = mul(a, a, p);
        n >>= 1;
    }
    return ans;
}
static int inv(int a, int p) { return power(a, p - 2, p); }
static int legendre(int a, int p) {
    if (a == 0) return 0;
    return power(a, (p - 1) / 2, p) == 1 ? 1 : -1;
}

// Polynomial helpers, coefficients low to high.
static void trim(vector<int>& a) {
    while (a.size() > 1 && a.back() == 0) a.pop_back();
}
static vector<int> poly_mod(vector<int> a, const vector<int>& b, int p) {
    trim(a);
    vector<int> bb = b;
    trim(bb);
    const int db = int(bb.size()) - 1;
    const int ib = inv(bb.back(), p);
    while (int(a.size()) - 1 >= db && !(a.size() == 1 && a[0] == 0)) {
        const int d = int(a.size()) - 1 - db;
        const int q = mul(a.back(), ib, p);
        for (int j = 0; j <= db; ++j)
            a[d + j] = sub(a[d + j], mul(q, bb[j], p), p);
        trim(a);
    }
    return a;
}
static bool squarefree(const std::array<int,6>& f, int p) {
    vector<int> a(f.begin(), f.end());
    vector<int> b(5);
    for (int i = 1; i <= 5; ++i) b[i - 1] = mul(i % p, f[i], p);
    trim(a); trim(b);
    while (!(b.size() == 1 && b[0] == 0)) {
        vector<int> r = poly_mod(a, b, p);
        a.swap(b); b.swap(r);
    }
    trim(a);
    return a.size() == 1;
}

struct Fp2 { int a, b; }; // a+b*z, z^2=d
static Fp2 fp2_add(Fp2 x, Fp2 y, int p) {
    return {add(x.a,y.a,p), add(x.b,y.b,p)};
}
static Fp2 fp2_mul(Fp2 x, Fp2 y, int d, int p) {
    return {add(mul(x.a,y.a,p),mul(d,mul(x.b,y.b,p),p),p),
            add(mul(x.a,y.b,p),mul(x.b,y.a,p),p)};
}
static Fp2 eval_fp2(const std::array<int,6>& f, Fp2 x, int d, int p) {
    Fp2 y{f[5],0};
    for (int i=4;i>=0;--i) y=fp2_add(fp2_mul(y,x,d,p),{f[i],0},p);
    return y;
}
static int chi_fp2(Fp2 x, int d, int p) {
    if (x.a == 0 && x.b == 0) return 0;
    // The quadratic character over F_{p^2} is the character of the norm.
    int norm = sub(mul(x.a,x.a,p),mul(d,mul(x.b,x.b,p),p),p);
    return legendre(norm,p);
}

static int required_prime_to_p(int n, int p) {
    while (n % p == 0) n /= p;
    return n;
}

// Return false on a denominator/boundary of the affine (u,v)-chart.
static bool finite_model(int u, int v, int p, std::array<int,6>& f,
                         int& a, int& b) {
    if (p == 2) return false;
    const int half = inv(2,p);
    const int r = sub(1,mul(u,u,p),p);
    if (u==0 || v==0 || u==p-1 || v==p-1 || u==v || add(u,v,p)==0)
        return false;
    auto Aval=[&](int z) {
        int N=2%p;
        N=add(mul(N,z,p),4%p,p);
        N=add(mul(N,z,p),6%p,p);
        N=add(mul(N,z,p),8%p,p);
        N=add(mul(N,z,p),10%p,p);
        N=add(mul(N,z,p),5%p,p);
        int zp1=add(z,1,p);
        return mul(N,inv(mul(2%p,mul(zp1,zp1,p),p),p),p);
    };
    const int Au=Aval(u),Av=Aval(v);
    const int den=sub(mul(u,u,p),mul(v,v,p),p);
    if (den==0) return false;
    b=mul(sub(Av,Au,p),inv(den,p),p);
    a=sub(Au,mul(b,r,p),p);
    // h(1)=a+b-5/2.  Its vanishing is the degenerate contact-7 boundary.
    if (sub(add(a,b,p),mul(5 % p,half,p),p)==0) return false;
    // f = x^5+(b^2-7)x^4+(2ab+21)x^3
    //       +(a^2-7b-35)x^2+(2b-7a+35)x+(2a-35/4).
    const int quarter=mul(half,half,p);
    f[0]=sub(mul(2,a,p),mul(35 % p,quarter,p),p);
    f[1]=add(sub(mul(2,b,p),mul(7 % p,a,p),p),35 % p,p);
    f[2]=sub(sub(mul(a,a,p),mul(7 % p,b,p),p),35 % p,p);
    f[3]=add(mul(2,mul(a,b,p),p),21 % p,p);
    f[4]=sub(mul(b,b,p),7 % p,p);
    f[5]=1;
    return true;
}

static int jacobian_order(const std::array<int,6>& f, int p) {
    int S1=0;
    for (int x=0;x<p;++x) {
        int y=f[5];
        for (int i=4;i>=0;--i) y=add(mul(y,x,p),f[i],p);
        S1 += legendre(y,p);
    }
    const int a1=-S1;
    int d=2;
    while (legendre(d,p)!=-1) ++d;
    int S2=0;
    for (int aa=0;aa<p;++aa) for (int bb=0;bb<p;++bb)
        S2 += chi_fp2(eval_fp2(f,{aa,bb},d,p),d,p);
    const int numer=S2+a1*a1;
    if (numer & 1) {
        std::cerr << "internal parity failure at p=" << p << "\n";
        std::exit(3);
    }
    const int a2=numer/2;
    return p*p+1+a2-(p+1)*a1;
}

struct Mask {
    int p;
    vector<unsigned char> status; // 0 unknown, 1 pass, 2 fail
    uint64_t chart=0,good=0,pass=0,fail=0,base_anomaly=0;
};

static Mask make_mask(int p) {
    Mask M; M.p=p; M.status.assign(p*p,0);
    const int req=required_prime_to_p(84,p);
    const int base=required_prime_to_p(28,p);
    for (int u=0;u<p;++u) for (int v=0;v<p;++v) {
        std::array<int,6> f{}; int a=0,b=0;
        if (!finite_model(u,v,p,f,a,b)) continue;
        ++M.chart;
        if (!squarefree(f,p)) continue;
        ++M.good;
        const int n=jacobian_order(f,p);
        if (n % base != 0) ++M.base_anomaly;
        unsigned char st = n % req == 0 ? 1 : 2;
        M.status[u*p+v]=M.status[v*p+u]=st;
        if (st==1) ++M.pass; else ++M.fail;
    }
    return M;
}

struct Rat { int n,d; };
static bool rat_less(const Rat& x,const Rat& y) {
    return (__int128)x.n*y.d < (__int128)y.n*x.d;
}
static bool boundary(const Rat& x) {
    // u=0 is a singular/contact degeneration; u=-1 is an incompatible
    // pole.  The formerly excluded u=+1 branch is valid and is retained.
    return x.n==0 || x.n==-x.d;
}
static bool equal_rat(const Rat& x,const Rat& y) {
    return int64_t(x.n)*y.d == int64_t(y.n)*x.d;
}
static bool opposite(const Rat& x,const Rat& y) {
    return int64_t(x.n)*y.d == -int64_t(y.n)*x.d;
}
static int residue(const Rat& x,int p) {
    const int dd=mod(x.d,p);
    if (!dd) return -1;
    return mul(mod(x.n,p),inv(dd,p),p);
}

struct Out {
    std::ostream& a; std::ostream* b;
    template<class T> Out& operator<<(const T& x) { a<<x; if(b) *b<<x; return *this; }
    Out& operator<<(std::ostream& (*pf)(std::ostream&)) { pf(a); if(b) pf(*b); return *this; }
};

int main(int argc,char** argv) {
    const int H=argc>=2 ? std::atoi(argv[1]) : 50;
    const bool branch_one=argc>=5 && std::string(argv[4])=="branch1";
    if (H<1) return 2;
    std::ofstream log_file;
    if (argc>=3) log_file.open(argv[2]);
    Out out{std::cout,log_file.is_open()?&log_file:nullptr};
    std::ofstream survivors;
    if (argc>=4) survivors.open(argv[3]);
    if (survivors.is_open()) survivors << "[\n";

    // p=3 is retained as a diagnostic; it cannot test the added 3-primary
    // class.  All other masks are genuine target-84 filters.
    const vector<int> primes={3,5,7,11,13,17,19,23,29,31,37,41,43,47,
                              53,59,61,67,71,73,79,83};
    vector<Mask> masks;
    out << "CONTACT7_TWO_ROOT_PLUS3_SIEVE\nheight "<<H
        <<" mode "<<(branch_one?"u=1 branch":"full surface")<<"\n";
    for (int p:primes) {
        Mask M=make_mask(p);
        out << "prime "<<p<<" total "<<p*p<<" chart "<<M.chart
            <<" good "<<M.good<<" pass "<<M.pass<<" fail "<<M.fail
            <<" unknown "<<(uint64_t(p)*p-M.good)
            <<" base_anomaly "<<M.base_anomaly
            <<" required "<<required_prime_to_p(84,p)<<"\n";
        if (p!=3) masks.push_back(std::move(M));
    }
    std::stable_sort(masks.begin(),masks.end(),[](const Mask& A,const Mask& B){
        return (__int128)A.fail*(B.p+1)*(B.p+1)
             > (__int128)B.fail*(A.p+1)*(A.p+1);
    });
    out << "filter_order";
    for (const auto& M:masks) out << " "<<M.p;
    out << "\n";

    uint64_t raw=0, exact_boundary=0, tested=0, survived=0;
    vector<uint64_t> killed(masks.size(),0);
    bool first_survivor=true;
    auto TestPair = [&](const Rat& u,const Rat& v) {
        ++raw;
        if (boundary(u)||boundary(v)||equal_rat(u,v)||opposite(u,v)) {
            ++exact_boundary; return;
        }
        ++tested;
        size_t k=0;
        for (;k<masks.size();++k) {
            const int ur=residue(u,masks[k].p),vr=residue(v,masks[k].p);
            if (ur<0||vr<0) continue; // projective-infinity disk: unknown
            const unsigned char st=masks[k].status[ur*masks[k].p+vr];
            if (st==2) break;
        }
        if (k<masks.size()) { ++killed[k]; return; }
        ++survived;
        out << "SURVIVOR "<<u.n<<"/"<<u.d<<" "<<v.n<<"/"<<v.d<<"\n";
        if (survivors.is_open()) {
            if (!first_survivor) survivors << ",\n";
            first_survivor=false;
            survivors << "<"<<u.n<<","<<u.d<<","<<v.n<<","<<v.d<<">";
        }
    };
    if (branch_one) {
        uint64_t nvals=0;
        const Rat one{1,1};
        for (int d=1;d<=H;++d) for (int n=-H;n<=H;++n) {
            if (std::gcd(std::abs(n),d)!=1) continue;
            ++nvals;
            TestPair(one,{n,d});
        }
        out << "rational_parameters "<<nvals<<" branch_raw_parameters "<<raw<<"\n";
    } else {
        vector<Rat> vals;
        for (int d=1;d<=H;++d) for (int n=-H;n<=H;++n)
            if (std::gcd(std::abs(n),d)==1) vals.push_back({n,d});
        std::sort(vals.begin(),vals.end(),rat_less);
        out << "rational_parameters "<<vals.size()<<" unordered_raw_pairs "
            <<uint64_t(vals.size())*(vals.size()-1)/2<<"\n";
        for (size_t i=0;i<vals.size();++i)
            for (size_t j=i+1;j<vals.size();++j) TestPair(vals[i],vals[j]);
    }
    if (survivors.is_open()) survivors << "\n]\n";
    out << "SUMMARY raw_pairs "<<raw<<" exact_boundary "<<exact_boundary
        <<" tested "<<tested<<" survivors "<<survived<<"\n";
    out << "KILL_COUNTS";
    for (size_t k=0;k<masks.size();++k) if(killed[k])
        out << " "<<masks[k].p<<":"<<killed[k];
    out << "\n";
    return 0;
}
