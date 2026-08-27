// Projective rational-point sieve for the genus-6 trigonal quotient D of
// the contact-30 C3-root cover.
//
// We use the integral bidegree-(10,3) equation F(Z,W;X,Y)=0 documented in
// notes/contact30_g6_quotient_points_2026_07_11.md.  The two searches are:
//   z=Z/W of bounded projective height, asking for some rho=[X:Y] in P1(Q);
//   rho=X/Y of bounded projective height, asking for some z=[Z:W] in P1(Q).
// Reduction of an integral bihomogeneous equation is an unconditional local
// test even at bad/singular fibers: a rational point always reduces to a
// projective F_p-point.  Thus no special bad-disk exceptions are needed.
//
// Build/run:
//   c++ -O3 -std=c++17 code/contact30_g6_quotient_projective_sieve.cpp \
//       -o /tmp/contact30_g6_quotient_projective_sieve
//   /tmp/contact30_g6_quotient_projective_sieve 25000 quiet


#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using std::int64_t;
using std::string;
using std::vector;

static int mod(int64_t a, int p) {
    int r = int(a % p);
    return r < 0 ? r + p : r;
}
static int add(int a, int b, int p) { return mod(int64_t(a) + b, p); }
static int sub(int a, int b, int p) { return mod(int64_t(a) - b, p); }
static int mul(int a, int b, int p) { return int((int64_t(a) * b) % p); }
static int power(int a, int64_t n, int p) {
    int r = 1;
    while (n) {
        if (n & 1) r = mul(r, a, p);
        a = mul(a, a, p);
        n >>= 1;
    }
    return r;
}
static int inv(int a, int p) { return power(a, p - 2, p); }

// Low-to-high coefficients.
static const vector<int64_t> PA = {
    2886221168640LL, -8735943819264LL, 11658061873152LL,
    -9038973468672LL, 4514219274240LL, -1519391001600LL,
    349528189824LL, -54340244352LL, 5471086032LL,
    -322506900LL, 8461353LL
};
static const vector<int64_t> PB = {
    -13343981568LL, 31885688832LL, -32393066496LL, 18329174016LL,
    -6337639680LL, 1375090368LL, -183283224LL, 13749864LL,
    -445299LL
};
static const vector<int64_t> PC = {
    24776704LL, -40903680LL, 27411072LL, -9583360LL,
    1849992LL, -187500LL, 7813LL
};

static int eval_affine(const vector<int64_t>& f, int x, int p) {
    int v = 0;
    for (auto it = f.rbegin(); it != f.rend(); ++it)
        v = add(mul(v, x, p), mod(*it, p), p);
    return v;
}

// Coefficients [c0,c1,c2,c3] of F at a projective z-coordinate.
static std::array<int,4> fiber_at_z(int zidx, int p) {
    int d, pa, pb, pc;
    if (zidx == p) {
        d = 5 % p;
        pa = mod(PA.back(), p);
        pb = mod(PB.back(), p);
        pc = mod(PC.back(), p);
    } else {
        int z = zidx;
        d = sub(sub(mul(5 % p, mul(z,z,p), p), mul(20 % p,z,p), p),
                16 % p, p);
        pa = eval_affine(PA,z,p);
        pb = eval_affine(PB,z,p);
        pc = eval_affine(PC,z,p);
    }
    int d2=mul(d,d,p), d3=mul(d2,d,p), d4=mul(d2,d2,p), d5=mul(d4,d,p);
    int c3=mul(2 % p,d5,p);
    int c2=sub(pa,mul(3 % p,d5,p),p);
    int c1=mul(d,add(pb,mul(3 % p,d4,p),p),p);
    int c0=mul(d2,sub(pc,d3,p),p);
    return {c0,c1,c2,c3};
}

static int eval_rho(const std::array<int,4>& c, int ridx, int p) {
    if (ridx == p) return c[3]; // [X:Y]=[1:0]
    int v=0;
    for (int i=3;i>=0;--i) v=add(mul(v,ridx,p),c[i],p);
    return v;
}

struct Masks {
    int p;
    vector<unsigned char> z_allowed;
    vector<unsigned char> rho_allowed;
};

static Masks make_masks(int p) {
    Masks M;
    M.p=p;
    M.z_allowed.assign(p+1,0);
    M.rho_allowed.assign(p+1,0);
    vector<std::array<int,4>> fibers(p+1);
    for (int z=0;z<=p;++z) fibers[z]=fiber_at_z(z,p);
    for (int z=0;z<=p;++z) {
        for (int rho=0;rho<=p;++rho) if (eval_rho(fibers[z],rho,p)==0) {
            M.z_allowed[z]=1;
            M.rho_allowed[rho]=1;
        }
    }
    return M;
}

static int count_allowed(const vector<unsigned char>& a) {
    int n=0; for (auto x:a) n += x != 0; return n;
}

struct Survivor { int a,b; };

int main(int argc,char** argv) {
    const int H = argc>=2 ? std::atoi(argv[1]) : 10000;
    const bool verbose = !(argc>=3 && string(argv[2])=="quiet");
    if (H<1) return 2;
    const vector<int> primes = {
        7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
        79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
        157,163,167,173,179,181,191,193,197,199,211,223,227,229,
        233,239,241,251,257,263,269,271,277,281,283,293,307,311,
        313,317,331,337,347,349,353,359,367,373,379,383,389,397,401
    };
    vector<Masks> masks;
    for (int p:primes) masks.push_back(make_masks(p));
    vector<int> zorder(masks.size()), rorder(masks.size());
    std::iota(zorder.begin(),zorder.end(),0);
    std::iota(rorder.begin(),rorder.end(),0);
    std::stable_sort(zorder.begin(),zorder.end(),[&](int i,int j){
        return int64_t(count_allowed(masks[i].z_allowed))*(masks[j].p+1)
             < int64_t(count_allowed(masks[j].z_allowed))*(masks[i].p+1);
    });
    std::stable_sort(rorder.begin(),rorder.end(),[&](int i,int j){
        return int64_t(count_allowed(masks[i].rho_allowed))*(masks[j].p+1)
             < int64_t(count_allowed(masks[j].rho_allowed))*(masks[i].p+1);
    });

    std::cout << "CONTACT30_GENUS6_QUOTIENT_PROJECTIVE_SIEVE\n";
    std::cout << "height "<<H<<" primes "<<masks.size()<<"\n";
    if (verbose) for (const auto& M:masks)
        std::cout << "prime "<<M.p
                  << " z_allowed "<<count_allowed(M.z_allowed)<<"/"<<M.p+1
                  << " rho_allowed "<<count_allowed(M.rho_allowed)<<"/"<<M.p+1
                  << "\n";

    uint64_t primitive=0;
    vector<Survivor> zsurv, rsurv;
    vector<uint64_t> zdepth(masks.size()+1,0), rdepth(masks.size()+1,0);
    const size_t cap=100000;
    for (int b=1;b<=H;++b) {
        vector<int> bmod(masks.size()),binv(masks.size());
        for (size_t j=0;j<masks.size();++j) {
            int p=masks[j].p;
            bmod[j]=b%p;
            if (bmod[j]) binv[j]=inv(bmod[j],p);
        }
        for (int a=-H;a<=H;++a) {
            if (std::gcd(std::abs(a),b)!=1) continue;
            ++primitive;
            size_t iz=0;
            for (;iz<zorder.size();++iz) {
                int j=zorder[iz],p=masks[j].p;
                int idx=bmod[j] ? mul(mod(a,p),binv[j],p) : p;
                if (!masks[j].z_allowed[idx]) break;
            }
            ++zdepth[iz];
            if (iz==zorder.size() && zsurv.size()<cap) zsurv.push_back({a,b});
            size_t ir=0;
            for (;ir<rorder.size();++ir) {
                int j=rorder[ir],p=masks[j].p;
                int idx=bmod[j] ? mul(mod(a,p),binv[j],p) : p;
                if (!masks[j].rho_allowed[idx]) break;
            }
            ++rdepth[ir];
            if (ir==rorder.size() && rsurv.size()<cap) rsurv.push_back({a,b});
        }
    }
    std::cout << "SUMMARY primitive_parameters "<<primitive
              << " z_survivors "<<zsurv.size()
              << " rho_survivors "<<rsurv.size()<<"\n";
    for (const auto&s:zsurv) std::cout<<"Z_SURVIVOR "<<s.a<<"/"<<s.b<<"\n";
    for (const auto&s:rsurv) std::cout<<"RHO_SURVIVOR "<<s.a<<"/"<<s.b<<"\n";
    std::cout << "Z_INFINITY exact fiber is an irreducible cubic over Q\n";
    std::cout << "RHO_INFINITY requires 5*z^2-20*z-16=0, no rational z\n";
    std::cout << "Z_FILTER_DEPTH";
    for (size_t i=0;i<zdepth.size();++i) if(zdepth[i]) std::cout<<" "<<i<<":"<<zdepth[i];
    std::cout << "\nRHO_FILTER_DEPTH";
    for (size_t i=0;i<rdepth.size();++i) if(rdepth[i]) std::cout<<" "<<i<<":"<<rdepth[i];
    std::cout<<"\n";
    return 0;
}
