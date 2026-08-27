// Search the canonical (2,4) model of the [2,2,2,12] square cover.
//
// Starting with
//   H = A^2-B^2+C^2,
//   L = (A^2-B^2)(B^2-C^2),
// the two equations in the scratch calculation factor as
//   y^2 = H L,                 z^2 = H(H^2+L).
// On H != 0 put t=y/H and w=z/H and clear denominators.  The resulting
// primitive projective integral point (a,c,t;b,w) satisfies
//
//   b^2+w^2 = a^2+c^2+t^2,
//   b^2w^2  = a^2c^2+a^2t^2+c^2t^2.                 (1)
//
// If R=w^2-a^2-c^2 and S=b^2-a^2-c^2, then (1) is equivalent to
//
//   R S = (ac)^2,
//   w^2 = a^2+c^2+R,
//   b^2 = a^2+c^2+S,
//   t^2 = a^2+c^2+R+S.                              (2)
//
// This program enumerates the signed divisor pairs in (2).  It is complete
// for primitive rational points of projective height at most T on (1).
// Independent signs are quotiented out, as are S_3 on {a,c,t} and S_2 on
// {b,w}.  For every surviving orbit, all 6*2 original (A,B,C;t,w) charts are
// checked against the full nondegeneracy product from the scratch section.
//
// Build, for example:
//   g++ -O3 -std=c++17 -fopenmp code/search_22212_canonical.cpp -o search_22212_canonical
//
// Examples:
//   ./search_22212_canonical --self-test-only
//   ./search_22212_canonical --T 1200 --threads 8
//   ./search_22212_canonical --T 10000 --a-min 1 --a-max 500 --threads 8

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128;

namespace {

struct Options {
    u64 T = 1200;
    int threads = 1;
    u64 a_min = 1;
    u64 a_max = 0;  // 0 means T
    u64 c_min = 1;
    u64 c_max = 0;  // 0 means T
    bool self_test_only = false;
    bool show_all_charts = true;
};

struct Stats {
    u64 pairs = 0;
    u64 divisor_pairs = 0;
    u64 signed_candidates = 0;
    u64 bounds_pass = 0;
    u64 residue_pass = 0;
    u64 square_triples = 0;
    u64 primitive = 0;
    u64 equation_verified = 0;
    u64 nondegenerate_orbits = 0;
    u64 nondegenerate_charts = 0;

    Stats& operator+=(const Stats& o) {
        pairs += o.pairs;
        divisor_pairs += o.divisor_pairs;
        signed_candidates += o.signed_candidates;
        bounds_pass += o.bounds_pass;
        residue_pass += o.residue_pass;
        square_triples += o.square_triples;
        primitive += o.primitive;
        equation_verified += o.equation_verified;
        nondegenerate_orbits += o.nondegenerate_orbits;
        nondegenerate_charts += o.nondegenerate_charts;
        return *this;
    }
};

struct Chart {
    u64 A = 0, B = 0, C = 0, t = 0, w = 0;

    auto key() const { return std::tie(A, B, C, t, w); }
    bool operator<(const Chart& o) const { return key() < o.key(); }
    bool operator==(const Chart& o) const { return key() == o.key(); }
};

struct Hit {
    // left is sorted a<=c<=t; right is sorted b<=w.
    std::array<u64, 3> left{};
    std::array<u64, 2> right{};
    std::vector<Chart> charts;

    auto key() const {
        return std::tie(left[0], left[1], left[2], right[0], right[1]);
    }
    bool operator<(const Hit& o) const { return key() < o.key(); }
};

std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x != 0) {
        s.push_back(static_cast<char>('0' + x % 10));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::string to_string_i128(i128 x) {
    if (x >= 0) return to_string_u128(static_cast<u128>(x));
    // The values used here are far from the signed minimum.
    return "-" + to_string_u128(static_cast<u128>(-x));
}

u128 gcd_u128(u128 a, u128 b) {
    while (b != 0) {
        const u128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

std::string format_fraction(i128 numerator, i128 denominator) {
    if (denominator == 0) return "undefined";
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    const u128 an = numerator < 0 ? static_cast<u128>(-numerator)
                                  : static_cast<u128>(numerator);
    const u128 g = gcd_u128(an, static_cast<u128>(denominator));
    numerator /= static_cast<i128>(g);
    denominator /= static_cast<i128>(g);
    if (denominator == 1) return to_string_i128(numerator);
    return to_string_i128(numerator) + "/" + to_string_i128(denominator);
}

u64 parse_u64(const std::string& s, const char* name) {
    std::size_t pos = 0;
    unsigned long long v = 0;
    try {
        v = std::stoull(s, &pos);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("invalid ") + name + ": " + s);
    }
    if (pos != s.size())
        throw std::runtime_error(std::string("invalid ") + name + ": " + s);
    return static_cast<u64>(v);
}

int parse_int(const std::string& s, const char* name) {
    std::size_t pos = 0;
    long v = 0;
    try {
        v = std::stol(s, &pos);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("invalid ") + name + ": " + s);
    }
    if (pos != s.size() || v < 1 || v > std::numeric_limits<int>::max())
        throw std::runtime_error(std::string("invalid ") + name + ": " + s);
    return static_cast<int>(v);
}

void usage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " [T] [options]\n"
        << "  --T N                 canonical projective-height bound (default 1200)\n"
        << "  --threads N           OpenMP threads (default 1)\n"
        << "  --a-min N, --a-max N  inclusive range for smallest left coordinate\n"
        << "  --c-min N, --c-max N  inclusive range for middle left coordinate\n"
        << "  --first-chart-only    print one nondegenerate chart per orbit\n"
        << "  --self-test-only      validate the algebra and the two known points\n"
        << "  --help                show this message\n";
}

Options parse_options(int argc, char** argv) {
    Options o;
    bool positional_T_seen = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* name) -> std::string {
            if (++i >= argc) throw std::runtime_error(std::string("missing value for ") + name);
            return argv[i];
        };
        if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            std::exit(0);
        } else if (arg == "--T") {
            o.T = parse_u64(need_value("--T"), "T");
            positional_T_seen = true;
        } else if (arg == "--threads" || arg == "-j") {
            o.threads = parse_int(need_value("--threads"), "thread count");
        } else if (arg == "--a-min") {
            o.a_min = parse_u64(need_value("--a-min"), "a-min");
        } else if (arg == "--a-max") {
            o.a_max = parse_u64(need_value("--a-max"), "a-max");
        } else if (arg == "--c-min") {
            o.c_min = parse_u64(need_value("--c-min"), "c-min");
        } else if (arg == "--c-max") {
            o.c_max = parse_u64(need_value("--c-max"), "c-max");
        } else if (arg == "--self-test-only") {
            o.self_test_only = true;
        } else if (arg == "--first-chart-only") {
            o.show_all_charts = false;
        } else if (!arg.empty() && arg[0] != '-' && !positional_T_seen) {
            o.T = parse_u64(arg, "T");
            positional_T_seen = true;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    if (o.T < 1) throw std::runtime_error("T must be positive");
    // This also keeps all degree-four checks safely inside signed 128 bits.
    // An SPF array at this limit would already need about 4 GB.
    if (o.T > 1000000000ULL)
        throw std::runtime_error("T > 10^9 is unsupported");
    if (o.a_max == 0 || o.a_max > o.T) o.a_max = o.T;
    if (o.c_max == 0 || o.c_max > o.T) o.c_max = o.T;
    o.a_min = std::max<u64>(1, o.a_min);
    o.c_min = std::max<u64>(1, o.c_min);
    if (o.a_min > o.a_max) throw std::runtime_error("empty a range");
    if (o.c_min > o.c_max) throw std::runtime_error("empty c range");
    return o;
}

struct SquareSieve {
    static constexpr std::array<unsigned, 4> moduli{{64, 63, 65, 187}};
    std::array<std::vector<std::uint8_t>, 4> qr;

    SquareSieve() {
        for (std::size_t j = 0; j < moduli.size(); ++j) {
            const unsigned m = moduli[j];
            qr[j].assign(m, 0);
            for (unsigned x = 0; x < m; ++x) qr[j][(x * x) % m] = 1;
        }
    }

    bool passes(u128 n) const {
        for (std::size_t j = 0; j < moduli.size(); ++j) {
            if (!qr[j][static_cast<unsigned>(n % moduli[j])]) return false;
        }
        return true;
    }
};

bool square_root_bounded(u128 n, u64 bound, u64& root) {
    const long double approx = std::sqrt(static_cast<long double>(n));
    u64 r = static_cast<u64>(approx);
    if (r > bound) r = bound;
    while (r < bound && static_cast<u128>(r + 1) * (r + 1) <= n) ++r;
    while (static_cast<u128>(r) * r > n) --r;
    if (static_cast<u128>(r) * r != n) return false;
    root = r;
    return true;
}

std::vector<std::uint32_t> make_spf(u64 T) {
    // One uint32_t per integer; no divisor tables or per-pair factorization
    // caches are retained.  Memory use is therefore 4(T+1) bytes plus small
    // per-thread recursion and hit buffers.
    std::vector<std::uint32_t> spf(static_cast<std::size_t>(T) + 1, 0);
    if (T >= 1) spf[1] = 1;
    for (u64 i = 2; i <= T; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) continue;
        spf[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(i);
        if (i > T / i) continue;
        for (u64 j = i * i; j <= T; j += i) {
            auto& x = spf[static_cast<std::size_t>(j)];
            if (x == 0) x = static_cast<std::uint32_t>(i);
        }
    }
    return spf;
}

using Factorization = std::vector<std::pair<u64, unsigned>>;

void add_factorization(u64 n, const std::vector<std::uint32_t>& spf,
                       Factorization& f) {
    while (n > 1) {
        const u64 p = spf[static_cast<std::size_t>(n)];
        unsigned e = 0;
        do {
            n /= p;
            ++e;
        } while (n > 1 && n % p == 0);
        auto it = std::find_if(f.begin(), f.end(),
                               [&](const auto& pe) { return pe.first == p; });
        if (it == f.end()) f.emplace_back(p, e);
        else it->second += e;
    }
}

u64 gcd5(u64 a, u64 b, u64 c, u64 t, u64 w) {
    return std::gcd(std::gcd(std::gcd(std::gcd(a, b), c), t), w);
}

bool verify_surface(const std::array<u64, 3>& left,
                    const std::array<u64, 2>& right) {
    const i128 a = left[0], c = left[1], t = left[2];
    const i128 b = right[0], w = right[1];
    const i128 a2 = a * a, c2 = c * c, t2 = t * t;
    const i128 b2 = b * b, w2 = w * w;
    if (b2 + w2 != a2 + c2 + t2) return false;
    return b2 * w2 == a2 * c2 + a2 * t2 + c2 * t2;
}

bool verify_chart(const Chart& ch) {
    const i128 A = ch.A, B = ch.B, C = ch.C, t = ch.t, w = ch.w;
    const i128 A2 = A * A, B2 = B * B, C2 = C * C;
    const i128 H = A2 - B2 + C2;
    const i128 L = (A2 - B2) * (B2 - C2);
    return L == H * t * t && w * w == H + t * t;
}

bool discriminant_nonzero(u64 Au, u64 Bu, u64 Cu) {
    const i128 A = Au, B = Bu, C = Cu;
    if (A == 0 || B == 0 || C == 0) return false;
    // The paired linear factors are sign-independent after quotienting signs.
    if (A == B || A == C || B == C) return false;
    const i128 A2 = A * A, B2 = B * B, C2 = C * C;

    const i128 f10 = -A2 + B2 - 2 * C2;
    const i128 f11 = -A2 + B2 - C2;
    const i128 f12 = -A2 + B2 - A * C - C2;
    const i128 f13 = -A2 + B2 + A * C - C2;
    const i128 f14 = -A2 * B2 + B2 * B2 + 2 * A2 * C2
                    - 2 * B2 * C2 + C2 * C2;
    const i128 f15 = -A2 * A2 + A2 * B2 - A2 * C2
                    + B2 * C2 - C2 * C2;
    const i128 f16 = A2 * A2 - 2 * A2 * B2 + B2 * B2
                    + 2 * A2 * C2 - B2 * C2;
    return f10 != 0 && f11 != 0 && f12 != 0 && f13 != 0
        && f14 != 0 && f15 != 0 && f16 != 0;
}

std::vector<Chart> valid_charts(const std::array<u64, 3>& left,
                                const std::array<u64, 2>& right) {
    static constexpr int perms[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2},
        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };
    std::vector<Chart> out;
    out.reserve(12);
    for (const auto& p : perms) {
        for (int j = 0; j < 2; ++j) {
            Chart ch{left[p[0]], right[j], left[p[1]], left[p[2]], right[1-j]};
            if (!verify_chart(ch))
                throw std::runtime_error("internal chart-equation failure");
            if (discriminant_nonzero(ch.A, ch.B, ch.C)) out.push_back(ch);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

bool known_point_test(const std::array<u64, 3>& left,
                      const std::array<u64, 2>& right,
                      const Chart& expected) {
    if (!std::is_sorted(left.begin(), left.end())) return false;
    if (!std::is_sorted(right.begin(), right.end())) return false;
    if (!verify_surface(left, right)) return false;
    const i128 a = left[0], c = left[1];
    const i128 b = right[0], w = right[1];
    const i128 R = w * w - a * a - c * c;
    const i128 S = b * b - a * a - c * c;
    if (R * S != a * a * c * c) return false;
    const auto charts = valid_charts(left, right);
    return std::find(charts.begin(), charts.end(), expected) != charts.end();
}

void self_test() {
    // Canonical representatives of the two currently known nondegenerate
    // orbits.  The first has canonical height 1015; the second has height 266.
    const std::array<u64, 3> l1{{143, 408, 1015}};
    const std::array<u64, 2> r1{{437, 1013}};
    const Chart c1{408, 437, 143, 1015, 1013};
    const std::array<u64, 3> l2{{120, 143, 266}};
    const std::array<u64, 2> r2{{218, 241}};
    const Chart c2{120, 218, 143, 266, 241};
    if (!known_point_test(l1, r1, c1) || !known_point_test(l2, r2, c2))
        throw std::runtime_error("known-point self-test failed");

    // The alternate denominator-two chart (60,109,133) clears to
    // (120,266,143;218,241), hence the same second canonical orbit.
    std::array<u64, 3> alt_left{{120, 266, 143}};
    std::sort(alt_left.begin(), alt_left.end());
    if (alt_left != l2) throw std::runtime_error("denominator-chart self-test failed");
}

struct PairProcessor {
    const Options& opt;
    const SquareSieve& sieve;
    Stats& stats;
    std::vector<Hit>& hits;
    u64 a, c;
    u128 base, T2bound, product;

    void process_signed(u128 rmag, u128 smag, int sign) {
        ++stats.signed_candidates;
        u128 w2 = 0, b2 = 0, t2 = 0;
        if (sign > 0) {
            w2 = base + rmag;
            b2 = base + smag;
            t2 = base + rmag + smag;
            if (w2 > T2bound || b2 > T2bound || t2 > T2bound) return;
        } else {
            if (rmag >= base || smag >= base || rmag + smag >= base) return;
            w2 = base - rmag;
            b2 = base - smag;
            t2 = base - rmag - smag;
            if (w2 > T2bound || b2 > T2bound || t2 > T2bound) return;
        }
        // Since a,c must be the two smallest left coordinates in the chosen
        // canonical representative, t>=c.  Positivity also excludes every
        // chart that is automatically on the discriminant boundary.
        if (t2 < static_cast<u128>(c) * c || w2 == 0 || b2 == 0) return;
        ++stats.bounds_pass;
        if (!sieve.passes(w2) || !sieve.passes(b2) || !sieve.passes(t2)) return;
        ++stats.residue_pass;

        u64 wr = 0, br = 0, tr = 0;
        if (!square_root_bounded(w2, opt.T, wr)
            || !square_root_bounded(b2, opt.T, br)
            || !square_root_bounded(t2, opt.T, tr)) return;
        ++stats.square_triples;

        std::array<u64, 3> left{{a, c, tr}};
        std::sort(left.begin(), left.end());
        if (left[0] != a || left[1] != c) return;
        std::array<u64, 2> right{{br, wr}};
        std::sort(right.begin(), right.end());
        if (gcd5(left[0], right[0], left[1], left[2], right[1]) != 1) return;
        ++stats.primitive;
        if (!verify_surface(left, right))
            throw std::runtime_error("internal surface-equation failure");
        ++stats.equation_verified;

        auto charts = valid_charts(left, right);
        if (charts.empty()) return;
        ++stats.nondegenerate_orbits;
        stats.nondegenerate_charts += static_cast<u64>(charts.size());
        hits.push_back(Hit{left, right, std::move(charts)});
    }

    void divisor_leaf(u128 d) {
        // d<=ac is enforced during recursion, so this is one unordered pair.
        const u128 e = product / d;
        if (d > e) return;
        ++stats.divisor_pairs;
        // Any surviving radicand is <=T^2.  This early bound avoids carrying
        // enormous complementary divisors into the signed arithmetic.
        if (e > 2 * T2bound) return;
        process_signed(d, e, +1);
        process_signed(d, e, -1);
    }
};

void enumerate_divisors(const Factorization& f, std::size_t i, u128 current,
                        u128 limit, PairProcessor& pp) {
    if (i == f.size()) {
        pp.divisor_leaf(current);
        return;
    }
    const u64 p = f[i].first;
    const unsigned maxe = 2 * f[i].second;
    u128 value = current;
    for (unsigned e = 0; e <= maxe; ++e) {
        if (value > limit) break;
        enumerate_divisors(f, i + 1, value, limit, pp);
        if (e != maxe) value *= p;
    }
}

std::vector<Hit> run_search(const Options& opt, Stats& total) {
    const auto spf = make_spf(opt.T);
    const SquareSieve sieve;

    int actual_threads = 1;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(opt.threads);
    actual_threads = opt.threads;
#else
    (void)opt;
#endif
    std::vector<Stats> thread_stats(static_cast<std::size_t>(actual_threads));
    std::vector<std::vector<Hit>> thread_hits(static_cast<std::size_t>(actual_threads));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        Stats& stats = thread_stats[static_cast<std::size_t>(tid)];
        auto& hits = thread_hits[static_cast<std::size_t>(tid)];

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (long long ai = static_cast<long long>(opt.a_min);
             ai <= static_cast<long long>(opt.a_max); ++ai) {
            const u64 a = static_cast<u64>(ai);
            const u64 c0 = std::max<u64>(a, opt.c_min);
            for (u64 c = c0; c <= opt.c_max; ++c) {
                ++stats.pairs;
                Factorization f;
                f.reserve(16);
                add_factorization(a, spf, f);
                add_factorization(c, spf, f);
                std::sort(f.begin(), f.end());

                const u128 ac = static_cast<u128>(a) * c;
                PairProcessor pp{opt, sieve, stats, hits, a, c,
                                 static_cast<u128>(a) * a + static_cast<u128>(c) * c,
                                 static_cast<u128>(opt.T) * opt.T, ac * ac};
                enumerate_divisors(f, 0, 1, ac, pp);
            }
        }
    }

    std::vector<Hit> hits;
    for (std::size_t i = 0; i < thread_stats.size(); ++i) {
        total += thread_stats[i];
        hits.insert(hits.end(), std::make_move_iterator(thread_hits[i].begin()),
                    std::make_move_iterator(thread_hits[i].end()));
    }
    std::sort(hits.begin(), hits.end());
    hits.erase(std::unique(hits.begin(), hits.end(),
                           [](const Hit& x, const Hit& y) { return x.key() == y.key(); }),
               hits.end());
    // The counter is accumulated before this defensive deduplication.
    total.nondegenerate_orbits = static_cast<u64>(hits.size());
    total.nondegenerate_charts = 0;
    for (const auto& h : hits) total.nondegenerate_charts += h.charts.size();
    return hits;
}

void print_hit(const Hit& h, bool all_charts) {
    const i128 a = h.left[0], c = h.left[1];
    const i128 b = h.right[0], w = h.right[1];
    const i128 R = w * w - a * a - c * c;
    const i128 S = b * b - a * a - c * c;
    std::cout << "HIT left=[" << h.left[0] << ',' << h.left[1] << ',' << h.left[2]
              << "] right=[" << h.right[0] << ',' << h.right[1] << "]"
              << " R=" << to_string_i128(R) << " S=" << to_string_i128(S)
              << " valid_charts=" << h.charts.size() << '\n';
    const std::size_t n = all_charts ? h.charts.size() : std::min<std::size_t>(1, h.charts.size());
    for (std::size_t i = 0; i < n; ++i) {
        const auto& ch = h.charts[i];
        const i128 A = ch.A, B = ch.B, C = ch.C;
        const i128 A2 = A * A, B2 = B * B, C2 = C * C;
        const i128 H = A2 - B2 + C2;
        const i128 s_num = A * (A2 - B2 + 2 * C2);
        const i128 s_den = 2 * B * H;
        const i128 m_num = -2 * A2 * B2 + 2 * B2 * B2 + 2 * A2 * C2
                         - 6 * B2 * C2 + 4 * C2 * C2;
        const i128 m_den = 2 * A * B * H;
        std::cout << "  CHART A=" << ch.A << " B=" << ch.B << " C=" << ch.C
                  << " t=" << ch.t << " w=" << ch.w
                  << " s=" << format_fraction(s_num, s_den)
                  << " m=" << format_fraction(m_num, m_den)
                  << " n=" << format_fraction(A, B) << '\n';
    }
}

void print_stats(const Stats& s) {
    std::cout
        << "STATS pairs=" << s.pairs
        << " divisor_pairs=" << s.divisor_pairs
        << " signed_candidates=" << s.signed_candidates
        << " bounds_pass=" << s.bounds_pass
        << " residue_pass=" << s.residue_pass
        << " square_triples=" << s.square_triples
        << " primitive=" << s.primitive
        << " equation_verified=" << s.equation_verified
        << " nondegenerate_orbits=" << s.nondegenerate_orbits
        << " nondegenerate_charts=" << s.nondegenerate_charts << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parse_options(argc, argv);
        self_test();
        std::cout << "# self_test=ok\n";
        if (opt.self_test_only) return 0;

        std::cout << "# canonical_search T=" << opt.T << " threads=" << opt.threads
                  << " a_range=[" << opt.a_min << ',' << opt.a_max << ']'
                  << " c_range=[" << opt.c_min << ',' << opt.c_max << "]\n";
        std::cout << "# height=max(a,c,t,b,w); signs, S3(left), S2(right) quotiented\n";

        Stats stats;
        const auto hits = run_search(opt, stats);
        for (const auto& h : hits) print_hit(h, opt.show_all_charts);
        print_stats(stats);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
