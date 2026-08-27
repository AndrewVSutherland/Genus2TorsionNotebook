// Exact search on a rational two-parameter chart of the Clebsch--Klein
// cubic, followed by the two source-halving cover tests.
//
// Put
//
//   r1=x+y, r2=x-y, r3=z+w, r4=z-w, r5=-2(x+z).
//
// After x=1 and z=t, the CK cubic is the conic
//
//   y^2 + t*w^2 = (1+t)(t^2+3t+1),
//
// with rational section (y,w)=(-1,-t-2).  Intersecting the line
// w=-(t+2)+m(y+1) with the conic a second time gives the following
// projective CK point (m=0 is the other horizontal intersection):
//
//   R1 = 1+t(t+2)m,
//   R2 = tm(m-t-2),
//   R3 = -1+m+t(t+1)m^2,
//   R4 = 1+t-m-tm^2,
//   R5 = -(1+t)(1+tm^2).
//
// For t=a/b and m=c/d, multiplication by b^2*d^2 gives the integral
// formulas used below.  Thus every parameter pair is a rational CK point;
// unlike the quadratic base enumerator, no discriminant-square event is
// required first.
//
// Conversely, every smooth labelled CK point is in this chart:
//
//   t = (r3+r4)/(r1+r2),
//   y = (r1-r2)/(r1+r2),  w = (r3-r4)/(r1+r2),
//   m = (w+t+2)/(y+1).
//
// Indeed r1+r2=0 forces one of r3,r4,r5 to vanish by the CK equations,
// and y+1=0 is r1=0.  Both exceptional loci are outside the smooth open.
//
// On the real locus only three of the 15 nonzero 2-classes can halve:
// the orbit-01 class using the largest r_i^2, and the orbit-12 pairs in
// sorted positions {1,2} and {3,4}.  Testing those three exact
// Stoll--Zarhin radicand lists is therefore complete on this chart.
//
// Compile:
// On the current Apple-Silicon workspace, GMP is the x86_64 Homebrew build:
//   c++ -arch x86_64 -O3 -std=c++17 -I/usr/local/include \
//       code/elkies22210_ck_rational_param_search.cpp \
//       -L/usr/local/lib -lgmpxx -lgmp \
//       -o /tmp/elkies22210_ck_rational_param_search
// Run:
//   /tmp/elkies22210_ck_rational_param_search 50

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <gmpxx.h>

using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct Rat {
    i64 n;
    i64 d;
};

static bool square_residue(const mpz_class &n) {
    static const std::array<int, 11> mods =
        {64, 63, 65, 7, 11, 13, 17, 19, 23, 29, 31};
    static const auto residues = [] {
        std::array<std::array<bool, 65>, 11> table{};
        for (std::size_t i = 0; i < mods.size(); ++i)
            for (int x = 0; x < mods[i]; ++x)
                table[i][(x * x) % mods[i]] = true;
        return table;
    }();
    for (std::size_t i = 0; i < mods.size(); ++i) {
        const unsigned rem =
            unsigned(mpz_fdiv_ui(n.get_mpz_t(), unsigned(mods[i])));
        if (!residues[i][rem]) return false;
    }
    return true;
}

static std::vector<unsigned char> square_residues(unsigned q) {
    std::vector<unsigned char> out(q, 0);
    for (unsigned x = 0; x < q; ++x) out[(x * x) % q] = 1;
    return out;
}

static bool positive_square(const mpz_class &n) {
    if (n <= 0 || !square_residue(n)) return false;
    return mpz_perfect_square_p(n.get_mpz_t()) != 0;
}

static bool square_u128(u128 n) {
    long double d = std::sqrt((long double)n);
    std::uint64_t r = (std::uint64_t)d;
    while ((u128)r * r > n) --r;
    while ((u128)(r + 1) * (r + 1) <= n) ++r;
    return (u128)r * r == n;
}

static bool orbit01(const std::array<i64, 5> &r,
                    const std::array<int, 5> &ord) {
    const int i = ord[4];
    const i128 ai = i128(r[i]) * r[i];
    for (int z = 0; z < 4; ++z) {
        const i128 aj = i128(r[ord[z]]) * r[ord[z]];
        if (!square_u128(u128(ai - aj))) return false;
    }
    return true;
}

static i64 mod_i128(i128 x, i64 q) {
    i64 z = i64(x % q);
    return z < 0 ? z + q : z;
}

static bool orbit12_mod_q(const std::array<i64, 5> &r, int i, int j,
                          i64 q, const std::vector<unsigned char> &sq) {
    std::array<i64, 5> a;
    for (int z = 0; z < 5; ++z) a[z] = mod_i128(i128(r[z]) * r[z], q);
    int rem[3], nr = 0;
    for (int z = 0; z < 5; ++z) if (z != i && z != j) rem[nr++] = z;
    i128 g0 = -1;
    for (int z = 0; z < 3; ++z) g0 *= a[i] - a[rem[z]];
    if (!sq[std::size_t(mod_i128(g0, q))]) return false;
    for (int z = 0; z < 3; ++z) {
        i128 g = i128(a[rem[z]]) - a[j];
        for (int w = 0; w < 3; ++w) if (w != z) g *= a[i] - a[rem[w]];
        if (!sq[std::size_t(mod_i128(g, q))]) return false;
    }
    return true;
}

static std::array<mpz_class, 4>
orbit12_radicands(const std::array<mpz_class, 5> &aa, int i, int j) {
    int rem[3], nr = 0;
    for (int k = 0; k < 5; ++k) if (k != i && k != j) rem[nr++] = k;
    const int k = rem[0], l = rem[1], m = rem[2];
    return {
        -(aa[i] - aa[k]) * (aa[i] - aa[l]) * (aa[i] - aa[m]),
        (aa[k] - aa[j]) * (aa[i] - aa[l]) * (aa[i] - aa[m]),
        (aa[l] - aa[j]) * (aa[i] - aa[k]) * (aa[i] - aa[m]),
        (aa[m] - aa[j]) * (aa[i] - aa[k]) * (aa[i] - aa[l])
    };
}

static std::string tuple_string(const std::array<i64, 5> &r) {
    std::ostringstream out;
    out << "[";
    for (int i = 0; i < 5; ++i) {
        if (i) out << ",";
        out << r[i];
    }
    out << "]";
    return out.str();
}

int main(int argc, char **argv) {
    const int B = argc > 1 ? std::atoi(argv[1]) : 40;
    if (B < 1 || B > 500) {
        std::cerr << "require 1 <= B <= 500\n";
        return 2;
    }

    std::vector<Rat> vals;
    for (i64 d = 1; d <= B; ++d) {
        for (i64 n = -B; n <= B; ++n) {
            if (std::gcd(std::llabs(n), d) != 1) continue;
            vals.push_back({n, d});
        }
    }

    const auto sq1331 = square_residues(1331);
    const auto sq361 = square_residues(361);
    const auto sq529 = square_residues(529);

    std::uint64_t parameter_pairs = 0, nonzero_tuples = 0,
                  primitive_tuples = 0, smooth_tuples = 0,
                  orbit01_hits = 0, orbit12_hits = 0,
                  orbit12_pairs = 0, orbit12_mod1331 = 0,
                  orbit12_mod361 = 0, orbit12_mod529 = 0;
    std::array<std::uint64_t, 16> orbit12_exact_masks{};
    std::vector<std::string> hits;
    std::vector<std::string> crt_samples, partial_samples;

    for (const Rat &tv : vals) {
        const i64 a = tv.n, b = tv.d;
        for (const Rat &mv : vals) {
            ++parameter_pairs;
            const i64 c = mv.n, d = mv.d;

            // Common denominator b^2*d^2 for the five projective formulas.
            const i128 bb = i128(b) * b;
            const i128 dd = i128(d) * d;
            std::array<i128, 5> wide = {
                bb * dd + i128(a) * (a + 2 * b) * c * d,
                i128(a) * c * (i128(b) * c - i128(d) * (a + 2 * b)),
                -bb * dd + bb * c * d + i128(a) * (a + b) * c * c,
                i128(b) * (a + b) * dd - bb * c * d - i128(a) * b * c * c,
                -i128(a + b) * (i128(b) * dd + i128(a) * c * c)
            };

            if (std::any_of(wide.begin(), wide.end(), [](i128 z) { return z == 0; }))
                continue;
            ++nonzero_tuples;

            i128 g = 0;
            for (i128 z : wide) {
                i128 u = z < 0 ? -z : z;
                while (u) {
                    i128 r = g % u;
                    g = u;
                    u = r;
                }
            }
            if (g == 0) continue;
            std::array<i64, 5> r;
            bool fits = true;
            for (int i = 0; i < 5; ++i) {
                wide[i] /= g;
                if (wide[i] < i128(INT64_MIN) || wide[i] > i128(INT64_MAX)) fits = false;
                r[i] = fits ? i64(wide[i]) : 0;
            }
            if (!fits) continue;
            ++primitive_tuples;

            i128 p1 = 0, p3 = 0;
            for (i64 z : r) {
                p1 += z;
                p3 += i128(z) * z * z;
            }
            if (p1 != 0 || p3 != 0) {
                std::cerr << "CK identity failure at t=" << a << "/" << b
                          << " m=" << c << "/" << d << "\n";
                return 3;
            }

            std::array<int, 5> ord = {0, 1, 2, 3, 4};
            std::sort(ord.begin(), ord.end(), [&](int i, int j) {
                const i64 ai = std::llabs(r[i]), aj = std::llabs(r[j]);
                return ai < aj;
            });
            bool smooth = true;
            for (int i = 1; i < 5; ++i)
                if (std::llabs(r[ord[i - 1]]) == std::llabs(r[ord[i]])) smooth = false;
            if (!smooth) continue;
            ++smooth_tuples;

            static_assert(sizeof(long) >= sizeof(i64), "GMP conversion needs 64-bit long");
            if (orbit01(r, ord)) {
                ++orbit01_hits;
                std::ostringstream out;
                out << "HIT orbit=01 t=" << a << "/" << b
                    << " m=" << c << "/" << d << " r=" << tuple_string(r);
                hits.push_back(out.str());
            }
            for (const std::array<int, 2> pos :
                {std::array<int, 2>{0, 1}, std::array<int, 2>{2, 3}}) {
                const int i = ord[pos[0]], j = ord[pos[1]];
                ++orbit12_pairs;
                if (!orbit12_mod_q(r, i, j, 1331, sq1331)) continue;
                ++orbit12_mod1331;
                if (!orbit12_mod_q(r, i, j, 361, sq361)) continue;
                ++orbit12_mod361;
                if (!orbit12_mod_q(r, i, j, 529, sq529)) continue;
                ++orbit12_mod529;

                std::array<mpz_class, 5> aa;
                for (int z = 0; z < 5; ++z)
                    aa[z] = mpz_class(long(r[z])) * long(r[z]);
                const std::array<mpz_class, 4> gs = orbit12_radicands(aa, i, j);
                unsigned mask = 0;
                for (unsigned z = 0; z < 4; ++z)
                    if (positive_square(gs[z])) mask |= 1U << z;
                ++orbit12_exact_masks[mask];
                if (crt_samples.size() < 3) {
                    std::ostringstream out;
                    out << "CRT_SURVIVOR mask=" << mask
                        << " ranks=" << pos[0] + 1 << "," << pos[1] + 1
                        << " t=" << a << "/" << b << " m=" << c << "/" << d
                        << " r=" << tuple_string(r);
                    crt_samples.push_back(out.str());
                }
                if (mask != 0 && mask != 15 && partial_samples.size() < 20) {
                    std::ostringstream out;
                    out << "PARTIAL mask=" << mask
                        << " ranks=" << pos[0] + 1 << "," << pos[1] + 1
                        << " t=" << a << "/" << b << " m=" << c << "/" << d
                        << " r=" << tuple_string(r);
                    partial_samples.push_back(out.str());
                }
                if (mask == 15) {
                    ++orbit12_hits;
                    std::ostringstream out;
                    out << "HIT orbit=12 ranks=" << pos[0] + 1 << "," << pos[1] + 1
                        << " t=" << a << "/" << b << " m=" << c << "/" << d
                        << " r=" << tuple_string(r);
                    hits.push_back(out.str());
                }
            }
        }
    }

    std::cout << "ELKIES22210_CK_RATIONAL_PARAM_SEARCH\n";
    std::cout << "parameter_height " << B << "\n";
    std::cout << "rational_values " << vals.size() << "\n";
    std::cout << "parameter_pairs " << parameter_pairs << "\n";
    std::cout << "nonzero_tuples " << nonzero_tuples << "\n";
    std::cout << "primitive_tuples " << primitive_tuples << "\n";
    std::cout << "smooth_tuples " << smooth_tuples << "\n";
    std::cout << "orbit01_hits " << orbit01_hits << "\n";
    std::cout << "orbit12_real_eligible_pairs " << orbit12_pairs << "\n";
    std::cout << "orbit12_mod1331 " << orbit12_mod1331 << "\n";
    std::cout << "orbit12_mod361 " << orbit12_mod361 << "\n";
    std::cout << "orbit12_mod529 " << orbit12_mod529 << "\n";
    std::cout << "orbit12_exact_mask_counts";
    for (unsigned mask = 0; mask < 16; ++mask)
        if (orbit12_exact_masks[mask])
            std::cout << " " << mask << ":" << orbit12_exact_masks[mask];
    std::cout << "\n";
    for (const std::string &s : crt_samples) std::cout << s << "\n";
    for (const std::string &s : partial_samples) std::cout << s << "\n";
    std::cout << "orbit12_hits " << orbit12_hits << "\n";
    std::cout << "hits " << hits.size() << "\n";
    for (const std::string &s : hits) std::cout << s << "\n";
    std::cout << "DONE\n";
    return 0;
}
