// Projective finite-prime sieve for a rational root of the cubic factor in
// the simultaneous contact-5/contact-6 order-30 family.
//
// A parameter is R=a/b with |a|,b <= H, b>0 and gcd(a,b)=1.  At each prime
// we evaluate the original rational family on P^1(F_p).  A residue is kept
// if either
//   (1) the genus-2 fiber is good and C3 has an F_p-root, or
//   (2) the family/model is on a boundary residue.
// Thus the sieve never rejects a global point merely because its reduction
// lies in a bad disk.  The two signs in the order-30 parameterization are
// tracked independently.
//
// Build/run:
//   c++ -O3 -std=c++17 code/contact30_c3root_projective_sieve.cpp \
//       -o /tmp/contact30_c3root_projective_sieve
//   /tmp/contact30_c3root_projective_sieve 10000

// This is an exhaustive parameter-height sieve, not a proof about all
// rational parameters.  Any survivor must still be checked over Q.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

using std::array;
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
static int divide(int a, int b, int p) { return mul(a, inv(b, p), p); }

// Low-to-high coefficient vectors over F_p.
static void trim(vector<int>& f) {
    while (!f.empty() && f.back() == 0) f.pop_back();
}

static vector<int> derivative(const vector<int>& f, int p) {
    vector<int> d;
    for (int i = 1; i < int(f.size()); ++i) d.push_back(mul(i % p, f[i], p));
    trim(d);
    return d;
}

static vector<int> remainder(vector<int> a, const vector<int>& b, int p) {
    trim(a);
    if (b.empty()) return {};
    const int db = int(b.size()) - 1;
    const int ib = inv(b.back(), p);
    while (!a.empty() && int(a.size()) - 1 >= db) {
        int da = int(a.size()) - 1;
        int c = mul(a.back(), ib, p);
        int shift = da - db;
        for (int i = 0; i <= db; ++i)
            a[i + shift] = sub(a[i + shift], mul(c, b[i], p), p);
        trim(a);
    }
    return a;
}

static int gcd_degree(vector<int> a, vector<int> b, int p) {
    trim(a);
    trim(b);
    while (!b.empty()) {
        vector<int> r = remainder(a, b, p);
        a.swap(b);
        b.swap(r);
    }
    return a.empty() ? -1 : int(a.size()) - 1;
}

enum Boundary : unsigned {
    B_NONE = 0,
    B_PARAMETER_POLE = 1u << 0,
    B_U_ZERO = 1u << 1,
    B_C_ZERO = 1u << 2,
    B_Q_POLE = 1u << 3,
    B_WRONG_DEGREE = 1u << 4,
    B_SINGULAR = 1u << 5
};

struct ResidueResult {
    bool root = false;
    unsigned boundary = B_NONE;
};

static ResidueResult family_residue(int R, bool infinity, int eps, int p) {
    // Homogeneous R=[a:b].  R is affine when infinity=false; [1:0]
    // otherwise.  Using homogeneous input prevents a missing infinity disk.
    const int a = infinity ? 1 : R;
    const int b = infinity ? 0 : 1;
    const int two_inv = inv(2, p);
    const int aa = mul(a, a, p), ab = mul(a, b, p), bb = mul(b, b, p);
    const int denR = sub(aa, mul(5 % p, bb, p), p);
    if (denR == 0) return {false, B_PARAMETER_POLE};

    int tnum = add(sub(mul(5 % p, aa, p), mul(20 % p, ab, p), p),
                   mul(19 % p, bb, p), p);
    int ynum = mod(-2LL * add(sub(mul(5 % p, aa, p), mul(22 % p, ab, p), p),
                              mul(25 % p, bb, p), p), p);
    int t = divide(tnum, denR, p);
    int Y = divide(ynum, denR, p);
    int u = power(t, 3, p);
    if (u == 0) return {false, B_U_ZERO};

    int s = 0;
    s = add(s, power(t, 5, p), p);
    s = add(s, power(t, 4, p), p);
    s = add(s, mul(mul(5 % p, two_inv, p), power(t, 3, p), p), p);
    s = add(s, mul(two_inv, t, p), p);
    int contact = mul(mul(mul(t, sub(t, two_inv, p), p), add(t, 1, p), p), Y, p);
    s = add(s, eps == 1 ? contact : mod(-contact, p), p);

    int u2 = mul(u, u, p);
    int denomC = mul(2 % p, u, p);
    int C = divide(add(u2, 1, p), denomC, p);
    int c = divide(sub(u2, 1, p), denomC, p);
    if (c == 0) return {false, B_C_ZERO};

    int u3 = mul(u2, u, p), u4 = mul(u2, u2, p);
    int u5 = mul(u4, u, p), u6 = mul(u3, u3, p);
    int s2 = mul(s, s, p), s3 = mul(s2, s, p);
    int denq = u6;
    denq = add(denq, mul(6 % p, mul(u4, s, p), p), p);
    denq = sub(denq, mul(2 % p, u4, p), p);
    denq = add(denq, mul(15 % p, mul(u3, s, p), p), p);
    denq = sub(denq, mul(u, s3, p), p);
    denq = add(denq, u2, p);
    if (denq == 0) return {false, B_Q_POLE};

    int numq = 0;
    numq = add(numq, mul(15 % p, u5, p), p);
    numq = add(numq, mul(90 % p, u4, p), p);
    numq = add(numq, mul(20 % p, mul(u3, s, p), p), p);
    numq = sub(numq, mul(6 % p, mul(u2, s2, p), p), p);
    numq = add(numq, mul(231 % p, u3, p), p);
    numq = add(numq, mul(2 % p, mul(u2, s, p), p), p);
    numq = sub(numq, mul(15 % p, mul(u, s2, p), p), p);
    numq = add(numq, mul(90 % p, u2, p), p);
    numq = sub(numq, mul(20 % p, mul(u, s, p), p), p);
    numq = add(numq, mul(15 % p, u, p), p);
    numq = sub(numq, mul(2 % p, s, p), p);
    int q = divide(numq, denq, p);
    int A = mul(add(s, q, p), two_inv, p);
    int B = mul(sub(15 % p, mul(s, q, p), p), two_inv, p);

    // Q2=h-(x-1)^3 and C3=h+(x-1)^3.
    vector<int> Q2 = {add(C, 1, p), sub(B, 3 % p, p), add(A, 3 % p, p)};
    vector<int> C3 = {sub(C, 1, p), add(B, 3 % p, p), sub(A, 3 % p, p), 2 % p};
    vector<int> f(6, 0);
    for (int i = 0; i < int(Q2.size()); ++i)
        for (int j = 0; j < int(C3.size()); ++j)
            f[i + j] = add(f[i + j], mul(Q2[i], C3[j], p), p);
    trim(f);
    if (int(f.size()) != 6) return {false, B_WRONG_DEGREE};
    if (gcd_degree(f, derivative(f, p), p) > 0) return {false, B_SINGULAR};

    for (int rho = 0; rho < p; ++rho) {
        int value = 0;
        for (int i = 3; i >= 0; --i) value = add(mul(value, rho, p), C3[i], p);
        if (value == 0) return {true, B_NONE};
    }
    return {false, B_NONE};
}

struct PrimeMask {
    int p;
    // Bit 0 is eps=-1, bit 1 is eps=+1.  A bit is set for a root or boundary.
    vector<unsigned char> allowed;
    array<int, 2> roots{{0, 0}};
    array<int, 2> boundaries{{0, 0}};
    array<vector<string>, 2> boundary_labels;
};

static string boundary_name(unsigned b) {
    if (b & B_PARAMETER_POLE) return "parameter_pole";
    if (b & B_U_ZERO) return "u_zero";
    if (b & B_C_ZERO) return "c_zero";
    if (b & B_Q_POLE) return "q_pole";
    if (b & B_WRONG_DEGREE) return "wrong_degree";
    if (b & B_SINGULAR) return "singular";
    return "unknown";
}

static PrimeMask make_mask(int p) {
    PrimeMask M;
    M.p = p;
    M.allowed.assign(p + 1, 0);
    for (int idx = 0; idx <= p; ++idx) {
        bool infinity = idx == p;
        for (int branch = 0; branch < 2; ++branch) {
            int eps = branch == 0 ? -1 : 1;
            ResidueResult r = family_residue(idx, infinity, eps, p);
            if (r.root) ++M.roots[branch];
            if (r.boundary != B_NONE) {
                ++M.boundaries[branch];
                M.boundary_labels[branch].push_back(
                    (infinity ? string("inf") : std::to_string(idx)) + ":" +
                    boundary_name(r.boundary));
            }
            if (r.root || r.boundary != B_NONE) M.allowed[idx] |= (1u << branch);
        }
    }
    return M;
}

struct Candidate { int a, b, eps; };

static bool is_global_c_zero_boundary(int a, int b) {
    // Over Q, c=0 iff u^2=1, hence t=+/-1.  Substitution gives
    //   t= 1: (R-2)(R-3)=0,
    //   t=-1: (R-1)(3R-7)=0.
    return a == b || a == 2*b || a == 3*b || 3*a == 7*b;
}

int main(int argc, char** argv) {
    int H = argc >= 2 ? std::atoi(argv[1]) : 10000;
    bool verbose_masks = !(argc >= 3 && string(argv[2]) == "quiet");
    if (H < 1) return 2;
    const vector<int> primes = {
        7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
        67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
        131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
        263, 269, 271, 277, 281, 283, 293, 307, 311
    };
    vector<PrimeMask> masks;
    for (int p : primes) masks.push_back(make_mask(p));
    // Low allowed density first gives the fastest exact intersection.
    std::stable_sort(masks.begin(), masks.end(), [](const PrimeMask& a, const PrimeMask& b) {
        int ca = 0, cb = 0;
        for (auto x : a.allowed) ca += x != 0;
        for (auto x : b.allowed) cb += x != 0;
        return int64_t(ca) * (b.p + 1) < int64_t(cb) * (a.p + 1);
    });

    std::cout << "CONTACT30_C3ROOT_PROJECTIVE_SIEVE\n";
    std::cout << "height " << H << " primes " << masks.size() << "\n";
    if (verbose_masks) for (const auto& M : masks) {
        int union_allowed = 0;
        for (auto x : M.allowed) union_allowed += x != 0;
        std::cout << "prime " << M.p << " union_allowed " << union_allowed
                  << "/" << M.p + 1;
        for (int j = 0; j < 2; ++j) {
            std::cout << " eps" << (j == 0 ? -1 : 1)
                      << " roots " << M.roots[j]
                      << " boundary " << M.boundaries[j] << " [";
            for (size_t k = 0; k < M.boundary_labels[j].size(); ++k) {
                if (k) std::cout << ",";
                std::cout << M.boundary_labels[j][k];
            }
            std::cout << "]";
        }
        std::cout << "\n";
    }

    uint64_t primitive_parameters = 0;
    array<uint64_t, 2> raw_survivors{{0, 0}};
    array<uint64_t, 2> global_boundary_survivors{{0, 0}};
    array<uint64_t, 2> open_survivors{{0, 0}};
    vector<uint64_t> killed_after(masks.size() + 1, 0);
    vector<Candidate> candidates;
    const size_t candidate_print_cap = 10000;

    for (int b = 1; b <= H; ++b) {
        vector<int> binv(masks.size(), 0), bmod(masks.size(), 0);
        for (size_t j = 0; j < masks.size(); ++j) {
            int p = masks[j].p;
            bmod[j] = b % p;
            if (bmod[j]) binv[j] = inv(bmod[j], p);
        }
        for (int a = -H; a <= H; ++a) {
            if (std::gcd(std::abs(a), b) != 1) continue;
            ++primitive_parameters;
            unsigned char alive = 3;
            size_t used = 0;
            for (; used < masks.size() && alive; ++used) {
                const auto& M = masks[used];
                int idx;
                if (bmod[used] == 0) {
                    // gcd(a,b)=1 implies a != 0 mod p, hence [a:b]=infinity.
                    idx = M.p;
                } else {
                    idx = mul(mod(a, M.p), binv[used], M.p);
                }
                alive &= M.allowed[idx];
            }
            ++killed_after[used];
            for (int branch = 0; branch < 2; ++branch) if (alive & (1u << branch)) {
                ++raw_survivors[branch];
                if (is_global_c_zero_boundary(a, b)) {
                    ++global_boundary_survivors[branch];
                } else {
                    ++open_survivors[branch];
                    if (candidates.size() < candidate_print_cap)
                        candidates.push_back({a, b, branch == 0 ? -1 : 1});
                }
            }
        }
    }

    std::cout << "SUMMARY primitive_parameters " << primitive_parameters
              << " branch_tests " << 2 * primitive_parameters
              << " raw_survivors_eps_minus " << raw_survivors[0]
              << " raw_survivors_eps_plus " << raw_survivors[1]
              << " global_c_zero_boundary_eps_minus " << global_boundary_survivors[0]
              << " global_c_zero_boundary_eps_plus " << global_boundary_survivors[1]
              << " open_survivors_eps_minus " << open_survivors[0]
              << " open_survivors_eps_plus " << open_survivors[1] << "\n";
    std::cout << "GLOBAL_BOUNDARY R=1,2,3,7/3 both signs; exact identity c=0\n";
    std::cout << "CANDIDATES_STORED " << candidates.size() << "\n";
    for (const auto& c : candidates)
        std::cout << "CANDIDATE " << c.a << "/" << c.b << " eps " << c.eps << "\n";
    std::cout << "FILTER_DEPTH";
    for (size_t i = 0; i < killed_after.size(); ++i)
        if (killed_after[i]) std::cout << " " << i << ":" << killed_after[i];
    std::cout << "\n";
    return (open_survivors[0] || open_survivors[1]) ? 1 : 0;
}
