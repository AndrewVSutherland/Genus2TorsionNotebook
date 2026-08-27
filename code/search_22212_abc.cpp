// Exhaustive bounded-height search on the [2,2,2,12] K3 fibre product.
//
// We seek primitive positive triples (A,B,C), up to a supplied height, for
// which both sextics in the accompanying notes are integer squares and every
// factor in the stated discriminant is nonzero.  The two sextics are evaluated
// through the identities
//
//   D = A^2 - B^2 + C^2,
//   F = (A^2-B^2)(B^2-C^2)D,
//   G = F + D^3.
//
// The equations are invariant under A <-> C.  We therefore enumerate A<C
// and test the discriminant in both orientations.  Positivity leaves exactly
// the two C-ranges described in scan_B().  Six precomputed 64-point block
// masks discard nonsquares modulo small primes without visiting individual C's.
// The remaining candidates pass four scalar prime filters and a 2-adic filter.
//
// Build:
//   g++ -O3 -march=native -std=c++17 -fopenmp code/search_22212_abc.cpp
//       -o /tmp/search_22212_abc
//
// Examples:
//   /tmp/search_22212_abc 500 --threads 8
//   /tmp/search_22212_abc 5000 --threads 8 --a-min 1 --a-max 1000
//       --output data/22212_abc_H5000_a1_1000.tsv

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using i128 = __int128_t;
using u128 = __uint128_t;

namespace {

constexpr std::size_t kMaxScalarFilters = 8;
constexpr std::uint64_t kMaximumHeight = 1000000;

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
    // All quantities in this program are far from the signed i128 minimum.
    return "-" + to_string_u128(static_cast<u128>(-x));
}

i128 abs_i128(i128 x) { return x < 0 ? -x : x; }

i128 gcd_i128(i128 a, i128 b) {
    a = abs_i128(a);
    b = abs_i128(b);
    while (b != 0) {
        const i128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

std::uint64_t floor_sqrt_u128(u128 n) {
    if (n == 0) return 0;
    const long double approximation =
        std::sqrt(static_cast<long double>(n));
    std::uint64_t r = static_cast<std::uint64_t>(approximation);
    while (static_cast<u128>(r) * r > n) --r;
    while (r != std::numeric_limits<std::uint64_t>::max()) {
        const std::uint64_t next = r + 1;
        if (static_cast<u128>(next) * next > n) break;
        r = next;
    }
    return r;
}

bool is_square_abs(i128 n, std::uint64_t& root) {
    const u128 a = static_cast<u128>(abs_i128(n));
    root = floor_sqrt_u128(a);
    return static_cast<u128>(root) * root == a;
}

std::uint64_t floor_sqrt_u64(std::uint64_t n) {
    return floor_sqrt_u128(static_cast<u128>(n));
}

int mod_norm(std::int64_t x, int m) {
    int r = static_cast<int>(x % m);
    return r < 0 ? r + m : r;
}

bool same_nonzero_sign(i128 a, i128 b, i128 c) {
    return (a > 0 && b > 0 && c > 0) ||
           (a < 0 && b < 0 && c < 0);
}

struct Fraction {
    i128 numerator = 0;
    i128 denominator = 1;
};

struct OrbitKey {
    std::array<std::uint64_t, 3> left{};
    std::array<std::uint64_t, 2> right{};
};

Fraction reduced_fraction(i128 numerator, i128 denominator) {
    if (denominator == 0) throw std::runtime_error("zero rational denominator");
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    const i128 g = gcd_i128(numerator, denominator);
    return {numerator / g, denominator / g};
}

struct Hit {
    std::uint64_t A = 0;
    std::uint64_t B = 0;
    std::uint64_t C = 0;
    std::uint64_t y = 0;
    std::uint64_t z = 0;
    std::uint64_t canonical_A = 0;
    Fraction s;
    Fraction m;
    Fraction n;
    OrbitKey orbit;
};

bool operator<(const Hit& lhs, const Hit& rhs) {
    return std::tie(lhs.A, lhs.B, lhs.C, lhs.y, lhs.z) <
           std::tie(rhs.A, rhs.B, rhs.C, rhs.y, rhs.z);
}

Hit make_hit(std::uint64_t A, std::uint64_t B, std::uint64_t C,
             std::uint64_t yroot, std::uint64_t zroot,
             std::uint64_t canonical_A, const OrbitKey& orbit) {
    const i128 a = A;
    const i128 b = B;
    const i128 c = C;
    const i128 a2 = a * a;
    const i128 b2 = b * b;
    const i128 c2 = c * c;
    const i128 d = a2 - b2 + c2;

    const i128 s_num = a * (a2 - b2 + 2 * c2);
    const i128 s_den = 2 * b * d;
    const i128 m_num =
        -2 * a2 * b2 + 2 * b2 * b2 + 2 * a2 * c2 -
        6 * b2 * c2 + 4 * c2 * c2;
    const i128 m_den = 2 * a * b * d;

    Hit h;
    h.A = A;
    h.B = B;
    h.C = C;
    h.y = yroot;
    h.z = zroot;
    h.canonical_A = canonical_A;
    h.s = reduced_fraction(s_num, s_den);
    h.m = reduced_fraction(m_num, m_den);
    h.n = reduced_fraction(a, b);
    h.orbit = orbit;
    return h;
}

// The full nonvanishing product from the question, tested factor by factor so
// that the (very large) product itself is never formed.
bool discriminant_nonzero(std::uint64_t A, std::uint64_t B,
                          std::uint64_t C) {
    const i128 a = A;
    const i128 b = B;
    const i128 c = C;
    if (a == 0 || b == 0 || c == 0) return false;
    if (b == c || a == b || a == c) return false;

    const i128 x = a * a;
    const i128 y = b * b;
    const i128 z = c * c;
    const i128 d = x - y + z;
    const i128 ac = a * c;

    if (-x + y - 2 * z == 0) return false;
    if (-d == 0) return false;
    if (-d - ac == 0) return false;
    if (-d + ac == 0) return false;

    const i128 q5 = -x * y + y * y + 2 * x * z - 2 * y * z + z * z;
    const i128 q6 = -x * x + x * y - x * z + y * z - z * z;
    const i128 q7 = x * x - 2 * x * y + y * y + 2 * x * z - y * z;
    return q5 != 0 && q6 != 0 && q7 != 0;
}

// Direct expanded formulas, used only by startup regression tests.
std::pair<i128, i128> direct_forms(std::uint64_t A, std::uint64_t B,
                                   std::uint64_t C) {
    const i128 x = static_cast<i128>(A) * A;
    const i128 y = static_cast<i128>(B) * B;
    const i128 z = static_cast<i128>(C) * C;
    const i128 F =
        x * x * y - 2 * x * y * y + y * y * y - x * x * z +
        3 * x * y * z - 2 * y * y * z - x * z * z + y * z * z;
    const i128 G =
        x * x * x - 2 * x * x * y + x * y * y + 2 * x * x * z -
        3 * x * y * z + y * y * z + 2 * x * z * z -
        2 * y * z * z + z * z * z;
    return {F, G};
}

std::pair<i128, i128> factored_forms(std::uint64_t A, std::uint64_t B,
                                     std::uint64_t C) {
    const i128 x = static_cast<i128>(A) * A;
    const i128 y = static_cast<i128>(B) * B;
    const i128 z = static_cast<i128>(C) * C;
    const i128 d = x - y + z;
    const i128 F = (x - y) * (y - z) * d;
    return {F, F + d * d * d};
}

bool startup_self_test() {
    const std::array<std::array<std::uint64_t, 3>, 8> samples{{
        {{1, 2, 3}}, {{2, 5, 4}}, {{7, 3, 11}}, {{8, 13, 9}},
        {{120, 218, 143}}, {{408, 437, 143}}, {{143, 437, 408}},
        {{499, 317, 211}}
    }};
    for (const auto& t : samples) {
        if (direct_forms(t[0], t[1], t[2]) !=
            factored_forms(t[0], t[1], t[2])) {
            std::cerr << "self-test failed: expanded/factored mismatch at "
                      << t[0] << ',' << t[1] << ',' << t[2] << '\n';
            return false;
        }
    }

    // Check that the positivity parametrization is exhaustive after A<->C.
    for (std::uint64_t aa = 1; aa <= 9; ++aa) {
        for (std::uint64_t bb = 1; bb <= 9; ++bb) {
            for (std::uint64_t cc = 1; cc <= 9; ++cc) {
                if (aa == cc) continue;
                const std::uint64_t a = std::min(aa, cc);
                const std::uint64_t c = std::max(aa, cc);
                const auto fg = factored_forms(a, bb, c);
                const bool positive = fg.first > 0 && fg.second > 0;
                const bool in_ranges =
                    bb > a && ((c > bb) ||
                               (c < bb && a * a + c * c < bb * bb));
                if (positive != in_ranges) {
                    std::cerr << "self-test failed: positivity ranges at "
                              << a << ',' << bb << ',' << c << '\n';
                    return false;
                }
            }
        }
    }

    struct Known {
        std::uint64_t A, B, C, y, z;
    };
    const std::array<Known, 2> known{{
        {120, 218, 143, 3371550, 3054675},
        {408, 437, 143, 4116840, 4108728}
    }};
    for (const Known& k : known) {
        const auto fg = factored_forms(k.A, k.B, k.C);
        if (fg.first != static_cast<i128>(k.y) * k.y ||
            fg.second != static_cast<i128>(k.z) * k.z ||
            !discriminant_nonzero(k.A, k.B, k.C)) {
            std::cerr << "self-test failed: known point " << k.A << ','
                      << k.B << ',' << k.C << '\n';
            return false;
        }
    }
    return true;
}

std::vector<unsigned char> square_residues(int modulus) {
    std::vector<unsigned char> result(static_cast<std::size_t>(modulus), 0);
    for (int r = 0; r < modulus; ++r) {
        result[static_cast<std::size_t>(1LL * r * r % modulus)] = 1;
    }
    return result;
}

bool double_square_residue(int a, int b, int c, int modulus,
                           const std::vector<unsigned char>& qr) {
    const int x = static_cast<int>(1LL * a * a % modulus);
    const int y = static_cast<int>(1LL * b * b % modulus);
    const int z = static_cast<int>(1LL * c * c % modulus);
    const int d = mod_norm(static_cast<std::int64_t>(x) - y + z, modulus);
    const int f = mod_norm(1LL * mod_norm(x - y, modulus) *
                               mod_norm(y - z, modulus) % modulus * d,
                           modulus);
    const int g = mod_norm(static_cast<std::int64_t>(f) +
                               1LL * d * d % modulus * d,
                           modulus);
    return qr[static_cast<std::size_t>(f)] &&
           qr[static_cast<std::size_t>(g)];
}

// For fixed (A mod p,B mod p,C0 mod p), table[] is a 64-bit mask whose j-th
// bit says whether C=C0+j satisfies both square conditions modulo p.
struct BlockPrime {
    int modulus = 0;
    std::vector<std::uint64_t> table;
    u128 allowed_bits = 0;

    BlockPrime(int p, int threads) : modulus(p) {
        const std::size_t p3 = static_cast<std::size_t>(p) * p * p;
        table.resize(p3);
        const auto qr = square_residues(p);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
        for (int ab = 0; ab < p * p; ++ab) {
            const int a = ab / p;
            const int b = ab % p;
            std::array<unsigned char, 128> allowed{};
            for (int c = 0; c < p; ++c) {
                allowed[static_cast<std::size_t>(c)] =
                    static_cast<unsigned char>(
                        double_square_residue(a, b, c, p, qr));
            }
            for (int c0 = 0; c0 < p; ++c0) {
                std::uint64_t mask = 0;
                int c = c0;
                for (unsigned bit = 0; bit < 64; ++bit) {
                    if (allowed[static_cast<std::size_t>(c)]) {
                        mask |= std::uint64_t{1} << bit;
                    }
                    if (++c == p) c = 0;
                }
                table[(static_cast<std::size_t>(ab) * p) + c0] = mask;
            }
        }

        for (std::uint64_t mask : table) {
            allowed_bits += static_cast<unsigned>(__builtin_popcountll(mask));
        }
    }

    std::uint64_t mask(std::uint64_t A, std::uint64_t B,
                       std::uint64_t C0) const {
        const std::size_t a = static_cast<std::size_t>(A % modulus);
        const std::size_t b = static_cast<std::size_t>(B % modulus);
        const std::size_t c = static_cast<std::size_t>(C0 % modulus);
        return table[((a * modulus + b) * modulus) + c];
    }
};

struct ScalarFilter {
    int modulus = 0;
    std::vector<unsigned char> qr;

    explicit ScalarFilter(int m) : modulus(m), qr(square_residues(m)) {}

    bool accepts(std::uint64_t A, std::uint64_t B, std::uint64_t C) const {
        return double_square_residue(static_cast<int>(A % modulus),
                                     static_cast<int>(B % modulus),
                                     static_cast<int>(C % modulus),
                                     modulus, qr);
    }
};

struct Stats {
    u128 positive_points = 0;
    u128 blocks = 0;
    std::array<u128, 6> block_zero{};
    u128 block_mask_survivors = 0;
    std::array<u128, kMaxScalarFilters> scalar_reject{};
    u128 scalar_survivors = 0;
    u128 imprimitive = 0;
    u128 exact_tested = 0;
    u128 h_square = 0;
    u128 hl_square = 0;
    u128 exact_double_square = 0;
    u128 discriminant_tests = 0;
    u128 discriminant_valid = 0;

    Stats& operator+=(const Stats& rhs) {
        positive_points += rhs.positive_points;
        blocks += rhs.blocks;
        for (std::size_t i = 0; i < block_zero.size(); ++i)
            block_zero[i] += rhs.block_zero[i];
        block_mask_survivors += rhs.block_mask_survivors;
        for (std::size_t i = 0; i < scalar_reject.size(); ++i)
            scalar_reject[i] += rhs.scalar_reject[i];
        scalar_survivors += rhs.scalar_survivors;
        imprimitive += rhs.imprimitive;
        exact_tested += rhs.exact_tested;
        h_square += rhs.h_square;
        hl_square += rhs.hl_square;
        exact_double_square += rhs.exact_double_square;
        discriminant_tests += rhs.discriminant_tests;
        discriminant_valid += rhs.discriminant_valid;
        return *this;
    }
};

// Padding keeps independently updated per-thread counters off the same cache
// line.  These are the only mutable objects shared by worker threads.
struct alignas(64) PaddedStats {
    Stats value;
};

struct SearchData {
    std::uint64_t H = 0;
    const std::vector<std::uint64_t>* squares = nullptr;
    const std::vector<BlockPrime>* block_primes = nullptr;
    const std::vector<ScalarFilter>* scalar_filters = nullptr;
};

void test_candidate(std::uint64_t A, std::uint64_t B, std::uint64_t C,
                    const SearchData& data, Stats& stats,
                    std::vector<Hit>& hits) {
    const auto& filters = *data.scalar_filters;
    for (std::size_t k = 0; k < filters.size(); ++k) {
        if (!filters[k].accepts(A, B, C)) {
            ++stats.scalar_reject[k];
            return;
        }
    }
    ++stats.scalar_survivors;

    const std::uint64_t gcd_ab = std::gcd(A, B);
    if (std::gcd(gcd_ab, C) != 1) {
        ++stats.imprimitive;
        return;
    }

    ++stats.exact_tested;
    const i128 x = static_cast<i128>(A) * A;
    const i128 y = static_cast<i128>(B) * B;
    const i128 z = static_cast<i128>(C) * C;
    const i128 Hform = x - y + z;
    const i128 Lform = (x - y) * (y - z);
    const i128 Qform = Lform + Hform * Hform;
    const i128 common = gcd_i128(Hform, Lform);
    const i128 h = Hform / common;
    const i128 l = Lform / common;
    const i128 q = Qform / common;
    if (!same_nonzero_sign(h, l, q)) return;

    std::uint64_t rh = 0, rl = 0, rq = 0;
    if (!is_square_abs(h, rh)) return;
    ++stats.h_square;
    if (!is_square_abs(l, rl)) return;
    ++stats.hl_square;
    if (!is_square_abs(q, rq)) return;
    ++stats.exact_double_square;

    const u128 yroot128 = static_cast<u128>(common) * rh * rl;
    const u128 zroot128 = static_cast<u128>(common) * rh * rq;
    if (yroot128 > std::numeric_limits<std::uint64_t>::max() ||
        zroot128 > std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("square root does not fit uint64_t");
    }
    const std::uint64_t yroot = static_cast<std::uint64_t>(yroot128);
    const std::uint64_t zroot = static_cast<std::uint64_t>(zroot128);
    const i128 F = Hform * Lform;
    const i128 G = Hform * Qform;
    if (F != static_cast<i128>(yroot) * yroot ||
        G != static_cast<i128>(zroot) * zroot) {
        throw std::logic_error("coprime-factor square test/root reconstruction");
    }

    const u128 rA = static_cast<u128>(rh) * A;
    const u128 rB = static_cast<u128>(rh) * B;
    const u128 rC = static_cast<u128>(rh) * C;
    if (rA > std::numeric_limits<std::uint64_t>::max() ||
        rB > std::numeric_limits<std::uint64_t>::max() ||
        rC > std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("surface-orbit key does not fit uint64_t");
    }
    OrbitKey orbit{{static_cast<std::uint64_t>(rA),
                    static_cast<std::uint64_t>(rC), rl},
                   {static_cast<std::uint64_t>(rB), rq}};
    std::sort(orbit.left.begin(), orbit.left.end());
    std::sort(orbit.right.begin(), orbit.right.end());

    // A<->C preserves F and G, but not every displayed discriminant factor.
    ++stats.discriminant_tests;
    if (discriminant_nonzero(A, B, C)) {
        ++stats.discriminant_valid;
        hits.push_back(make_hit(A, B, C, yroot, zroot, A, orbit));
    }
    ++stats.discriminant_tests;
    if (discriminant_nonzero(C, B, A)) {
        ++stats.discriminant_valid;
        hits.push_back(make_hit(C, B, A, yroot, zroot, A, orbit));
    }
}

void scan_interval(std::uint64_t A, std::uint64_t B, std::uint64_t lo,
                   std::uint64_t hi, const SearchData& data, Stats& stats,
                   std::vector<Hit>& hits) {
    if (lo > hi) return;
    stats.positive_points += static_cast<u128>(hi - lo + 1);
    const auto& primes = *data.block_primes;

    for (std::uint64_t start = lo; start <= hi;) {
        const std::uint64_t remaining = hi - start + 1;
        const unsigned length =
            static_cast<unsigned>(std::min<std::uint64_t>(64, remaining));
        std::uint64_t mask = length == 64
                                 ? ~std::uint64_t{0}
                                 : ((std::uint64_t{1} << length) - 1);
        ++stats.blocks;
        for (std::size_t k = 0; k < primes.size(); ++k) {
            mask &= primes[k].mask(A, B, start);
            if (mask == 0) {
                ++stats.block_zero[k];
                break;
            }
        }
        stats.block_mask_survivors +=
            static_cast<unsigned>(__builtin_popcountll(mask));
        while (mask != 0) {
            const unsigned bit = static_cast<unsigned>(__builtin_ctzll(mask));
            test_candidate(A, B, start + bit, data, stats, hits);
            mask &= mask - 1;
        }
        if (remaining <= 64) break;
        start += 64;
    }
}

void scan_B(std::uint64_t A, std::uint64_t B, const SearchData& data,
            Stats& stats, std::vector<Hit>& hits) {
    const auto& sq = *data.squares;
    // Once A<C is imposed, F>0 and G>0 force B>A.  There are two cases:
    //   A < B < C;
    //   A < C < B and A^2+C^2 < B^2 (equivalently D<0).
    if (B < data.H) scan_interval(A, B, B + 1, data.H, data, stats, hits);

    const std::uint64_t bound2 = sq[B] - sq[A] - 1;
    std::uint64_t upper = floor_sqrt_u64(bound2);
    if (upper >= B) upper = B - 1;
    if (A + 1 <= upper)
        scan_interval(A, B, A + 1, upper, data, stats, hits);
}

struct Options {
    std::uint64_t height = 0;
    std::uint64_t a_min = 1;
    std::uint64_t a_max = 0;
    int threads = 1;
    std::string output = "-";
    bool self_test_only = false;
};

void usage(const char* program) {
    std::cerr
        << "usage: " << program << " HEIGHT [options]\n"
        << "options:\n"
        << "  --threads N       OpenMP worker count (default min(system,8))\n"
        << "  --a-min N         first canonical min(A,C) (default 1)\n"
        << "  --a-max N         last canonical min(A,C) (default HEIGHT-1)\n"
        << "  --output FILE     candidate TSV (default stdout; '-' means stdout)\n"
        << "  --self-test-only  run algebra/regression tests and stop\n"
        << "  --help            show this message\n";
}

std::uint64_t parse_u64(const std::string& text, const char* what) {
    std::size_t used = 0;
    const unsigned long long result = std::stoull(text, &used);
    if (used != text.size()) throw std::invalid_argument(what);
    return static_cast<std::uint64_t>(result);
}

Options parse_options(int argc, char** argv) {
    if (argc < 2) throw std::invalid_argument("missing HEIGHT");
    if (std::string(argv[1]) == "--help") {
        usage(argv[0]);
        std::exit(0);
    }
    Options options;
    options.height = parse_u64(argv[1], "invalid HEIGHT");
#ifdef _OPENMP
    options.threads = std::min(omp_get_max_threads(), 8);
#else
    options.threads = 1;
#endif
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* name) -> std::string {
            if (++i >= argc) throw std::invalid_argument(std::string("missing ") + name);
            return argv[i];
        };
        if (arg == "--threads") {
            options.threads = static_cast<int>(
                parse_u64(need_value("thread count"), "invalid thread count"));
        } else if (arg == "--a-min") {
            options.a_min = parse_u64(need_value("A minimum"), "invalid A minimum");
        } else if (arg == "--a-max") {
            options.a_max = parse_u64(need_value("A maximum"), "invalid A maximum");
        } else if (arg == "--output") {
            options.output = need_value("output filename");
        } else if (arg == "--self-test-only") {
            options.self_test_only = true;
        } else if (arg == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    if (options.height < 2 || options.height > kMaximumHeight)
        throw std::invalid_argument("HEIGHT must lie in [2,1000000]");
    if (options.threads < 1)
        throw std::invalid_argument("thread count must be positive");
    if (options.a_max == 0) options.a_max = options.height - 1;
    options.a_min = std::max<std::uint64_t>(1, options.a_min);
    options.a_max = std::min(options.a_max, options.height - 1);
    if (options.a_min > options.a_max)
        throw std::invalid_argument("empty canonical A range");
    return options;
}

void write_header(std::ostream& out) {
    out << "A\tB\tC\ty\tz\ts_num\ts_den\tm_num\tm_den\tn_num\tn_den"
           "\tcanonical_min_AC\tsurface_orbit_key\n";
}

void write_hit(std::ostream& out, const Hit& h) {
    out << h.A << '\t' << h.B << '\t' << h.C << '\t' << h.y << '\t'
        << h.z << '\t' << to_string_i128(h.s.numerator) << '\t'
        << to_string_i128(h.s.denominator) << '\t'
        << to_string_i128(h.m.numerator) << '\t'
        << to_string_i128(h.m.denominator) << '\t'
        << to_string_i128(h.n.numerator) << '\t'
        << to_string_i128(h.n.denominator) << '\t' << h.canonical_A << '\t'
        << h.orbit.left[0] << ',' << h.orbit.left[1] << ','
        << h.orbit.left[2] << '|' << h.orbit.right[0] << ','
        << h.orbit.right[1] << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parse_options(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        usage(argv[0]);
        return 2;
    }

    if (!startup_self_test()) return 1;
    if (options.self_test_only) {
        std::cerr << "SELF_TEST ok\n";
        return 0;
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    const auto started = std::chrono::steady_clock::now();

    std::unique_ptr<std::ofstream> output_file;
    std::ostream* output = &std::cout;
    if (options.output != "-") {
        output_file = std::make_unique<std::ofstream>(options.output);
        if (!*output_file) {
            std::cerr << "error: cannot open output file " << options.output << '\n';
            return 1;
        }
        output = output_file.get();
    }
    write_header(*output);

    std::vector<std::uint64_t> squares(options.height + 1);
    for (std::uint64_t v = 0; v <= options.height; ++v) squares[v] = v * v;

    // Descending large primes are empirically the strongest first masks.
    const std::array<int, 6> block_moduli{{79, 73, 71, 59, 43, 31}};
    std::vector<BlockPrime> block_primes;
    block_primes.reserve(block_moduli.size());
    for (int p : block_moduli) block_primes.emplace_back(p, options.threads);

    // Four further odd primes plus a 2-adic condition.  These are evaluated
    // only at the very sparse survivors of all six block masks.
    const std::array<int, 5> scalar_moduli{{109, 97, 89, 83, 64}};
    std::vector<ScalarFilter> scalar_filters;
    scalar_filters.reserve(scalar_moduli.size());
    for (int m : scalar_moduli) scalar_filters.emplace_back(m);

    SearchData data{options.height, &squares, &block_primes, &scalar_filters};
    std::vector<PaddedStats> per_thread(static_cast<std::size_t>(options.threads));
    std::vector<std::vector<Hit>> thread_hits(
        static_cast<std::size_t>(options.threads));

    u128 written = 0;
    for (std::uint64_t A = options.a_min; A <= options.a_max; ++A) {
        for (auto& hits : thread_hits) hits.clear();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(options.threads)
#endif
        for (std::int64_t Bs = static_cast<std::int64_t>(A + 1);
             Bs <= static_cast<std::int64_t>(options.height); ++Bs) {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            scan_B(A, static_cast<std::uint64_t>(Bs), data,
                   per_thread[static_cast<std::size_t>(tid)].value,
                   thread_hits[static_cast<std::size_t>(tid)]);
        }

        std::size_t count = 0;
        for (const auto& hits : thread_hits) count += hits.size();
        std::vector<Hit> ordered;
        ordered.reserve(count);
        for (auto& hits : thread_hits) {
            ordered.insert(ordered.end(),
                           std::make_move_iterator(hits.begin()),
                           std::make_move_iterator(hits.end()));
        }
        std::sort(ordered.begin(), ordered.end());
        for (const Hit& h : ordered) write_hit(*output, h);
        written += ordered.size();
    }
    output->flush();

    Stats totals;
    for (const PaddedStats& s : per_thread) totals += s.value;
    const auto stopped = std::chrono::steady_clock::now();
    const double seconds =
        std::chrono::duration<double>(stopped - started).count();

    std::cerr << "SEARCH height=" << options.height << " canonical_A=["
              << options.a_min << ',' << options.a_max << "] threads="
              << options.threads << '\n';
    std::cerr << "POSITIVE_POINTS " << to_string_u128(totals.positive_points)
              << " BLOCKS " << to_string_u128(totals.blocks) << '\n';
    for (std::size_t k = 0; k < block_primes.size(); ++k) {
        std::cerr << "BLOCK_ZERO p=" << block_primes[k].modulus << " count="
                  << to_string_u128(totals.block_zero[k]) << '\n';
    }
    std::cerr << "BLOCK_MASK_SURVIVORS "
              << to_string_u128(totals.block_mask_survivors) << '\n';
    for (std::size_t k = 0; k < scalar_filters.size(); ++k) {
        std::cerr << "SCALAR_REJECT m=" << scalar_filters[k].modulus
                  << " count=" << to_string_u128(totals.scalar_reject[k])
                  << '\n';
    }
    std::cerr << "SCALAR_SURVIVORS " << to_string_u128(totals.scalar_survivors)
              << " IMPRIMITIVE " << to_string_u128(totals.imprimitive) << '\n';
    std::cerr << "EXACT_TESTED " << to_string_u128(totals.exact_tested)
              << " H_SQUARE " << to_string_u128(totals.h_square)
              << " H_L_SQUARE " << to_string_u128(totals.hl_square)
              << " DOUBLE_SQUARE "
              << to_string_u128(totals.exact_double_square) << '\n';
    std::cerr << "DISCRIMINANT_TESTS "
              << to_string_u128(totals.discriminant_tests)
              << " VALID_ORIENTATIONS "
              << to_string_u128(totals.discriminant_valid) << '\n';
    std::cerr << "HITS_WRITTEN " << to_string_u128(written) << " ELAPSED_SECONDS "
              << std::fixed << std::setprecision(3) << seconds << '\n';
    return 0;
}
