// Exact bounded search for rational points on the two source-halving
// covers of Elkies' Clebsch--Klein [2,2,2,10] family.
//
// Compile:
//   c++ -O3 -std=c++17 code/elkies22210_source_halving_search.cpp \
//       -o code/elkies22210_source_halving_search
// Run:
//   code/elkies22210_source_halving_search 1000
//
// The search is complete for primitive integral projective tuples with
// max |r_i| <= H.  Permutation and global-sign symmetry are removed by
// imposing 0 < |r1| < ... < |r5| and r1 > 0.
//
// For every orbit-12 marked pair, the search also applies the compatible
// local square conditions modulo 11^3, 19^2, and 23^2.  These are necessary
// for a rational cover point and constitute a CRT sieve: the same marked
// pair must lie in an allowed boundary disk at all three forced primes.

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
#include <unordered_set>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

static i64 iabs64(i64 x) { return x < 0 ? -x : x; }

static bool square_u64(std::uint64_t n, std::uint64_t &root) {
    long double d = std::sqrt((long double)n);
    std::uint64_t r = (std::uint64_t)d;
    while ((u128)r*r > n) --r;
    while ((u128)(r+1)*(r+1) <= n) ++r;
    root = r;
    return (u128)r*r == n;
}

static bool square_pos_i128(i128 n) {
    if (n <= 0) return false;
    u128 u = (u128)n;
    long double d = std::sqrt((long double)u);
    u128 r = (u128)d;
    while (r*r > u) --r;
    while ((r+1)*(r+1) <= u) ++r;
    return r*r == u;
}

// This is the exact Stoll/Zarhin radicand -D*(c-B)/(c-A), after
// cancelling A-c from D=Product_c(A-c).  Cross-ratio quotients alone
// would omit the common squareclass and produce false positives.
static bool pair_halves(const std::array<i64,6> &alpha, int i, int j) {
    int rem[4], nr=0;
    for (int k=0;k<6;++k) if (k!=i && k!=j) rem[nr++]=k;
    for (int z=0;z<4;++z) {
        const int k=rem[z];
        i128 value=(i128)alpha[k]-alpha[j];
        for (int w=0;w<4;++w) if (w!=z)
            value *= (i128)alpha[i]-alpha[rem[w]];
        if (!square_pos_i128(value)) return false;
    }
    return true;
}

static std::string tuple_string(const std::array<i64,5> &r) {
    std::ostringstream o;
    o << "[";
    for (int i=0;i<5;++i) { if (i) o << ","; o << r[i]; }
    o << "]";
    return o.str();
}

static std::string square_key(const std::array<i64,5> &r) {
    std::array<i64,5> q;
    for (int i=0;i<5;++i) q[i]=r[i]*r[i];
    std::sort(q.begin(),q.end());
    std::ostringstream o;
    for (auto z:q) o << z << ',';
    return o.str();
}

static bool boundary_mod_p(const std::array<i64,5> &r, i64 p) {
    for (int i=0;i<5;++i) {
        if ((r[i]%p+p)%p==0) return true;
        i64 ai=(i64)(((i128)r[i]*r[i])%p);
        for (int j=0;j<i;++j) {
            i64 aj=(i64)(((i128)r[j]*r[j])%p);
            if (ai==aj) return true;
        }
    }
    return false;
}

static i64 mod_i128(i128 x, i64 q) {
    i64 z=(i64)(x%q);
    return z<0 ? z+q : z;
}

static std::vector<unsigned char> square_residues(i64 q) {
    std::vector<unsigned char> sq((std::size_t)q,0);
    for (i64 z=0;z<q;++z) sq[(std::size_t)((i128)z*z%q)]=1;
    return sq;
}

// Exact orbit-12 radicands for the marked pair {r_i^2,r_j^2}, reduced
// modulo q.  They are homogeneous of degree 6, so projective rescaling by
// a local unit multiplies every radicand by a square and does not affect
// this test.  The three unmarked indices may occur in any order.
static bool orbit12_mod_q(const std::array<i64,5> &r, int i, int j,
                          i64 q, const std::vector<unsigned char> &sq) {
    std::array<i64,5> a;
    for (int z=0;z<5;++z) a[z]=mod_i128((i128)r[z]*r[z],q);
    int rem[3], nr=0;
    for (int z=0;z<5;++z) if (z!=i && z!=j) rem[nr++]=z;
    i128 g0=-1;
    for (int z=0;z<3;++z) g0*=a[i]-a[rem[z]];
    if (!sq[(std::size_t)mod_i128(g0,q)]) return false;
    for (int z=0;z<3;++z) {
        i128 g=(i128)a[rem[z]]-a[j];
        for (int w=0;w<3;++w) if (w!=z) g*=a[i]-a[rem[w]];
        if (!sq[(std::size_t)mod_i128(g,q)]) return false;
    }
    return true;
}

int main(int argc, char **argv) {
    int H = argc > 1 ? std::atoi(argv[1]) : 300;
    std::size_t candidate_print_limit = argc > 2 ?
        (std::size_t)std::strtoull(argv[2],nullptr,10) : 100;
    if (H < 5) { std::cerr << "H must be at least 5\n"; return 2; }

    const i64 q11=11*11*11, q19=19*19, q23=23*23;
    const auto local_sq11=square_residues(q11);
    const auto local_sq19=square_residues(q19);
    const auto local_sq23=square_residues(q23);
    // Resolved open Hensel seeds from the normalized boundary classifier.
    // These startup checks guard the indexing and exact modular formulas
    // used by the global sieve.
    const std::array<i64,5> seed11={1,242,959,9,120};
    const std::array<i64,5> seed19={1,248,98,3,11};
    const std::array<i64,5> seed23={1,392,118,8,10};
    if (!orbit12_mod_q(seed11,0,1,q11,local_sq11)
        || !orbit12_mod_q(seed19,0,1,q19,local_sq19)
        || !orbit12_mod_q(seed23,0,1,q23,local_sq23)) {
        std::cerr << "local seed self-check failed\n";
        return 3;
    }

    // Fast necessary square-residue filters for the quadratic
    // discriminant.  They only reject nonsquares.
    bool sq64[64]={false}, sq63[63]={false}, sq65[65]={false};
    for (int i=0;i<320;++i) {
        sq64[(i*i)%64]=true; sq63[(i*i)%63]=true; sq65[(i*i)%65]=true;
    }

    std::uint64_t triples=0, disc_nonnegative=0, disc_modpass=0,
                  disc_square=0, integral_roots=0, open_ck=0,
                  primitive_ck=0, unique_curves=0;
    std::uint64_t orbit0_points=0, orbit1_points=0,
                  orbit0_pairs=0, orbit1_pairs=0;
    std::uint64_t forced_boundary_11_19=0,
                  forced_boundary_orbit0=0,
                  forced_boundary_orbit1=0;
    std::uint64_t orbit12_pairs_tested=0,
                  local11_pairs=0, local11_19_pairs=0,
                  local11_19_23_pairs=0,
                  local11_points=0, local11_19_points=0,
                  local11_19_23_points=0;
    std::unordered_set<std::string> seen;
    std::vector<std::string> hits, crt_candidates;

    for (i64 a=1;a<=H;++a) {
        for (i64 b=a+1;b<=H;++b) {
            for (int sb : {-1,1}) {
                i64 r2=sb*b;
                for (i64 c=b+1;c<=H;++c) {
                    for (int sc : {-1,1}) {
                        ++triples;
                        i64 r1=a, r3=sc*c;
                        i64 s=r1+r2+r3;
                        if (s==0) continue;
                        i128 A=(i128)r1*r1*r1+(i128)r2*r2*r2+(i128)r3*r3*r3;
                        i128 disc128=3*(i128)s*(4*A-(i128)s*s*s);
                        if (disc128 < 0 || disc128 > (i128)UINT64_MAX) continue;
                        ++disc_nonnegative;
                        std::uint64_t disc=(std::uint64_t)disc128;
                        if (!sq64[disc%64] || !sq63[disc%63] || !sq65[disc%65]) continue;
                        ++disc_modpass;
                        std::uint64_t rootu;
                        if (!square_u64(disc,rootu)) continue;
                        ++disc_square;
                        i64 root=(i64)rootu;
                        i64 den=6*s;
                        i128 num=-(i128)3*s*s+root;
                        if (den==0 || num%den!=0) continue;
                        i64 r4=(i64)(num/den);
                        i64 r5=-s-r4;
                        ++integral_roots;
                        if (r4==0 || r5==0 || iabs64(r4)>H || iabs64(r5)>H) continue;
                        if (iabs64(r4)>iabs64(r5)) std::swap(r4,r5);
                        if (!(c<iabs64(r4) && iabs64(r4)<iabs64(r5))) continue;
                        ++open_ck;
                        std::array<i64,5> r={r1,r2,r3,r4,r5};
                        i64 g=0; for (i64 z:r) g=std::gcd(g,iabs64(z));
                        if (g!=1) continue;
                        ++primitive_ck;
                        std::string key=square_key(r);
                        if (!seen.insert(key).second) continue;
                        ++unique_curves;

                        bool b11=boundary_mod_p(r,11), b19=boundary_mod_p(r,19),
                             b23=boundary_mod_p(r,23);
                        bool b29=boundary_mod_p(r,29), b31=boundary_mod_p(r,31);
                        if (b11 && b19) ++forced_boundary_11_19;
                        if (b11 && b19 && b23) ++forced_boundary_orbit1;
                        if (b11 && b19 && b29 && b31) ++forced_boundary_orbit0;

                        std::array<i64,6> alpha={0,r1*r1,r2*r2,r3*r3,r4*r4,r5*r5};
                        int c0=0,c1=0;
                        for (int i=1;i<6;++i) {
                            if (pair_halves(alpha,0,i)) {
                                ++c0;
                                std::ostringstream o;
                                o << "HIT orbit=0i pair=[0," << i << "] rs=" << tuple_string(r);
                                hits.push_back(o.str());
                            }
                        }
                        int c11=0,c1119=0,c111923=0;
                        for (int i=1;i<5;++i) for (int j=i+1;j<6;++j) {
                            ++orbit12_pairs_tested;
                            bool l11=orbit12_mod_q(r,i-1,j-1,q11,local_sq11);
                            bool l19=l11 && orbit12_mod_q(r,i-1,j-1,q19,local_sq19);
                            bool l23=l19 && orbit12_mod_q(r,i-1,j-1,q23,local_sq23);
                            if (l11) { ++local11_pairs; ++c11; }
                            if (l19) { ++local11_19_pairs; ++c1119; }
                            if (l23) {
                                ++local11_19_23_pairs; ++c111923;
                                if (crt_candidates.size()<candidate_print_limit) {
                                    std::ostringstream o;
                                    o << "CRT_CANDIDATE pair=[" << i << "," << j
                                      << "] rs=" << tuple_string(r);
                                    crt_candidates.push_back(o.str());
                                }
                            }
                            bool exact=pair_halves(alpha,i,j);
                            if (exact && !l23) {
                                std::cerr << "internal CRT necessity failure pair=["
                                          << i << "," << j << "] rs="
                                          << tuple_string(r) << "\n";
                                return 4;
                            }
                            if (exact) {
                                ++c1;
                                std::ostringstream o;
                                o << "HIT orbit=ij pair=[" << i << "," << j << "] rs=" << tuple_string(r);
                                hits.push_back(o.str());
                            }
                        }
                        if (c11) ++local11_points;
                        if (c1119) ++local11_19_points;
                        if (c111923) ++local11_19_23_points;
                        orbit0_pairs += c0; orbit1_pairs += c1;
                        if (c0) ++orbit0_points;
                        if (c1) ++orbit1_points;
                    }
                }
            }
        }
        if (a%100==0) {
            std::cerr << "PROGRESS a=" << a << " triples=" << triples
                      << " unique=" << unique_curves << " hits=" << hits.size() << "\n";
        }
    }

    std::cout << "ELKIES22210_SOURCE_HALVING_SEARCH\n";
    std::cout << "height " << H << "\n";
    std::cout << "triples " << triples << "\n";
    std::cout << "disc_nonnegative " << disc_nonnegative << "\n";
    std::cout << "disc_modpass " << disc_modpass << "\n";
    std::cout << "disc_square " << disc_square << "\n";
    std::cout << "integral_roots " << integral_roots << "\n";
    std::cout << "open_ck " << open_ck << "\n";
    std::cout << "primitive_ck " << primitive_ck << "\n";
    std::cout << "unique_curves " << unique_curves << "\n";
    std::cout << "orbit0_points " << orbit0_points << " orbit0_pairs " << orbit0_pairs << "\n";
    std::cout << "orbit1_points " << orbit1_points << " orbit1_pairs " << orbit1_pairs << "\n";
    std::cout << "forced_boundary_11_19 " << forced_boundary_11_19 << "\n";
    std::cout << "forced_boundary_orbit0_11_19_29_31 " << forced_boundary_orbit0 << "\n";
    std::cout << "forced_boundary_orbit1_11_19_23 " << forced_boundary_orbit1 << "\n";
    std::cout << "orbit12_crt_moduli " << q11 << " " << q19 << " " << q23 << "\n";
    std::cout << "orbit12_crt_product " << q11*q19*q23 << "\n";
    std::cout << "orbit12_pairs_tested " << orbit12_pairs_tested << "\n";
    std::cout << "orbit12_local11_pairs " << local11_pairs
              << " points " << local11_points << "\n";
    std::cout << "orbit12_local11_19_pairs " << local11_19_pairs
              << " points " << local11_19_points << "\n";
    std::cout << "orbit12_local11_19_23_pairs " << local11_19_23_pairs
              << " points " << local11_19_23_points << "\n";
    std::cout << "crt_candidates_printed " << crt_candidates.size()
              << " limit " << candidate_print_limit << "\n";
    for (const auto &h:crt_candidates) std::cout << h << "\n";
    std::cout << "hits " << hits.size() << "\n";
    for (const auto &h:hits) std::cout << h << "\n";
    std::cout << "DONE\n";
    return 0;
}
