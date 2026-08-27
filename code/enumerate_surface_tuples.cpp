#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using i128 = __int128_t;

static i128 abs128(i128 x) { return x < 0 ? -x : x; }

static i128 gcd128(i128 a, i128 b) {
    a = abs128(a);
    b = abs128(b);
    while (b != 0) {
        i128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static std::string to_string128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static bool is_square128(i128 n, i128& root) {
    if (n < 0) return false;
    i128 lo = 0;
    i128 hi = 1;
    while (hi <= n / hi) {
        hi <<= 1;
    }
    while (lo + 1 < hi) {
        i128 mid = (lo + hi) >> 1;
        if (mid <= n / mid) lo = mid;
        else hi = mid;
    }
    root = lo;
    return lo * lo == n;
}

static bool make_curve_tuple(i128 a, i128 b, i128 c, i128 d_num, i128 d_den,
                             long long max_abs, std::vector<i128>& curve) {
    if (d_den == 0) return false;
    if (d_den < 0) {
        d_den = -d_den;
        d_num = -d_num;
    }

    i128 gden = gcd128(d_num, d_den);
    d_num /= gden;
    d_den /= gden;

    std::vector<i128> signed_tuple = {a * d_den, b * d_den, c * d_den, d_num};
    i128 g = 0;
    for (i128 x : signed_tuple) g = gcd128(g, x);
    if (g == 0) return false;

    curve.clear();
    for (i128 x : signed_tuple) {
        i128 y = abs128(x / g);
        if (y == 0 || y > max_abs) return false;
        curve.push_back(y);
    }
    std::sort(curve.begin(), curve.end());
    for (int i = 1; i < 4; ++i) {
        if (curve[i] == curve[i-1]) return false;
    }
    return true;
}

static std::string tuple_key(const std::vector<i128>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ",";
        os << to_string128(v[i]);
    }
    os << "]";
    return os.str();
}

static std::string tuple_line(const std::vector<i128>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ",";
        os << to_string128(v[i]);
    }
    os << "]\n";
    return os.str();
}

static long long mod128(i128 x, long long m) {
    i128 r = x % static_cast<i128>(m);
    if (r < 0) r += m;
    return static_cast<long long>(r);
}

static std::vector<bool> square_residue_table(long long modulus) {
    std::vector<bool> ok(static_cast<size_t>(modulus), false);
    for (long long x = 0; x < modulus; ++x) {
        ok[static_cast<size_t>((static_cast<__int128_t>(x) * x) % modulus)] = true;
    }
    return ok;
}

static std::vector<long long> parse_prime_list(const std::string& s) {
    std::vector<long long> out;
    if (s.empty()) return out;
    size_t start = 0;
    while (start <= s.size()) {
        size_t comma = s.find(',', start);
        std::string part = s.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!part.empty()) out.push_back(std::atoll(part.c_str()));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return out;
}

static long long square_mod_ll(long long a, long long p) {
    a %= p;
    if (a < 0) a += p;
    return static_cast<long long>((static_cast<__int128_t>(a) * a) % p);
}

static bool bad_reduction_at_prime(const std::vector<i128>& curve, long long p) {
    bool seen[256] = {false};
    for (i128 x : curve) {
        long long r = mod128(x, p);
        if (r == 0) return true;
        long long sq = square_mod_ll(r, p);
        if (seen[sq]) return true;
        seen[sq] = true;
    }
    return false;
}

static bool bad_reduction_at_all_required_primes(const std::vector<i128>& curve,
                                                 const std::vector<long long>& primes) {
    for (long long p : primes) {
        if (!bad_reduction_at_prime(curve, p)) return false;
    }
    return true;
}

static long long powmod_ll(long long a, long long e, long long m) {
    long long r = 1 % m;
    a = mod128(a, m);
    while (e > 0) {
        if (e & 1LL) r = static_cast<long long>((static_cast<__int128_t>(r) * a) % m);
        a = static_cast<long long>((static_cast<__int128_t>(a) * a) % m);
        e >>= 1LL;
    }
    return r;
}

static int legendre_symbol_ll(long long u, long long p) {
    u = mod128(u, p);
    if (u == 0) return 0;
    long long r = powmod_ll(u, (p - 1) / 2, p);
    return r == 1 ? 1 : -1;
}

static long long ipow_ll(long long base, int exp) {
    long long out = 1;
    for (int i = 0; i < exp; ++i) out *= base;
    return out;
}

static bool real_partition_12_34_possible_mod_pk(const std::vector<i128>& curve,
                                                  long long p, int depth) {
    long long q = ipow_ll(p, depth);
    const int pairs[4][2] = {{0,2}, {0,3}, {1,2}, {1,3}};
    bool have_class = false;
    int seen_class = 0;
    for (int k = 0; k < 4; ++k) {
        int i = pairs[k][0], j = pairs[k][1];
        i128 diff = curve[i] * curve[i] - curve[j] * curve[j];
        long long r = mod128(diff, q);
        if (r == 0) continue;
        int v = 0;
        while (v < depth && r % p == 0) {
            ++v;
            r /= p;
        }
        int leg = legendre_symbol_ll(r, p);
        int cls = (v & 1) * 2 + (leg == 1 ? 0 : 1);
        if (have_class && cls != seen_class) return false;
        have_class = true;
        seen_class = cls;
    }
    return true;
}

static bool real_partition_12_34_possible_at_all_primes(const std::vector<i128>& curve,
                                                         const std::vector<long long>& primes,
                                                         int depth) {
    for (long long p : primes) {
        if (!real_partition_12_34_possible_mod_pk(curve, p, depth)) return false;
    }
    return true;
}


static void maybe_add_curve(const std::vector<i128>& curve,
                            const std::vector<long long>& required_bad_primes,
                            int real_local_depth,
                            std::unordered_set<std::string>& seen,
                            std::vector<std::string>& rows) {
    if (!bad_reduction_at_all_required_primes(curve, required_bad_primes)) return;
    if (real_local_depth > 0 &&
        !real_partition_12_34_possible_at_all_primes(curve, required_bad_primes, real_local_depth)) return;
    std::string key = tuple_key(curve);
    if (seen.insert(key).second) rows.push_back(tuple_line(curve));
}


static long long mod_ll(long long x, long long m) {
    long long r = x % m;
    if (r < 0) r += m;
    return r;
}

static const char* BoundaryComponentNames[10] = {
    "Za", "Zb", "Zc", "Zd", "Eab", "Eac", "Ead", "Ebc", "Ebd", "Ecd"
};

struct BoundaryRequirement {
    long long p;
    unsigned short component_mask;
    std::string label;
    std::vector<unsigned short> masks;
    std::vector<unsigned char> partition_masks;
};

static std::string normalize_component_name(std::string s) {
    for (char& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (s == "any" || s == "*") return "any";
    if (s == "d" || s == "dcollision" || s == "d-collision") return "D";
    if (s == "n" || s == "nond" || s == "non-d" || s == "non_d" || s == "non-d-zero" || s == "non_d_zero") return "N";
    if (s.size() == 2 && s[0] == 'z') {
        s[0] = 'Z';
        return s;
    }
    if (s.size() == 3 && s[0] == 'e') {
        s[0] = 'E';
        return s;
    }
    return s;
}

static unsigned short all_boundary_components_mask() {
    return static_cast<unsigned short>((1u << 10) - 1);
}

static unsigned short boundary_component_mask(const std::string& raw) {
    std::string name = normalize_component_name(raw);
    if (name == "any") return all_boundary_components_mask();
    if (name == "D") {
        // Collisions involving the solved coordinate d: Ead, Ebd, Ecd.
        return static_cast<unsigned short>((1u << 6) | (1u << 8) | (1u << 9));
    }
    if (name == "N") {
        // Zero components and collisions not involving d.
        return static_cast<unsigned short>((1u << 0) | (1u << 1) | (1u << 2) | (1u << 3) |
                                           (1u << 4) | (1u << 5) | (1u << 7));
    }
    for (int i = 0; i < 10; ++i) {
        if (name == BoundaryComponentNames[i]) return static_cast<unsigned short>(1u << i);
    }
    std::cerr << "unknown boundary component/group: " << raw << "\n";
    std::exit(2);
}

static std::string boundary_component_label(const std::string& raw) {
    std::string name = normalize_component_name(raw);
    if (name == "any" || name == "D" || name == "N") return name;
    for (int i = 0; i < 10; ++i) {
        if (name == BoundaryComponentNames[i]) return name;
    }
    std::cerr << "unknown boundary component/group: " << raw << "\n";
    std::exit(2);
}

static size_t residue_index_prime(long long a, long long b, long long c, long long p) {
    return (static_cast<size_t>(a) * static_cast<size_t>(p) + static_cast<size_t>(b)) *
           static_cast<size_t>(p) + static_cast<size_t>(c);
}

static unsigned short boundary_mask_for_residues(long long a, long long b, long long c,
                                                 long long d, long long p) {
    long long x[4] = {mod_ll(a, p), mod_ll(b, p), mod_ll(c, p), mod_ll(d, p)};
    unsigned short mask = 0;
    for (int i = 0; i < 4; ++i) {
        if (x[i] == 0) mask |= static_cast<unsigned short>(1u << i);
    }
    const int pairs[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
    for (int k = 0; k < 6; ++k) {
        int i = pairs[k][0], j = pairs[k][1];
        if (square_mod_ll(x[i], p) == square_mod_ll(x[j], p)) {
            mask |= static_cast<unsigned short>(1u << (4 + k));
        }
    }
    return mask;
}

static unsigned char a2244_partition_mask_for_residues(long long a, long long b,
                                                        long long c, long long d,
                                                        long long p) {
    long long x[4] = {mod_ll(a, p), mod_ll(b, p), mod_ll(c, p), mod_ll(d, p)};
    long long sq[4] = {square_mod_ll(x[0], p), square_mod_ll(x[1], p),
                       square_mod_ll(x[2], p), square_mod_ll(x[3], p)};
    const int partitions[3][4][2] = {
        {{0,2}, {0,3}, {1,2}, {1,3}}, // 12|34
        {{0,1}, {0,3}, {2,1}, {2,3}}, // 13|24
        {{0,1}, {0,2}, {3,1}, {3,2}}  // 14|23
    };
    unsigned char mask = 0;
    for (int part = 0; part < 3; ++part) {
        bool killed = false;
        bool have_class = false;
        int seen_class = 0;
        for (int k = 0; k < 4; ++k) {
            int i = partitions[part][k][0];
            int j = partitions[part][k][1];
            long long diff = mod_ll(sq[i] - sq[j], p);
            if (diff == 0) continue;
            int cls = legendre_symbol_ll(diff, p) == 1 ? 0 : 1;
            if (have_class && cls != seen_class) {
                killed = true;
                break;
            }
            have_class = true;
            seen_class = cls;
        }
        if (!killed) mask |= static_cast<unsigned char>(1u << part);
    }
    return mask;
}

static unsigned char k3_partition_mask_for_abc(long long a, long long b, long long c,
                                                long long p, unsigned short component_mask) {
    a = mod_ll(a, p);
    b = mod_ll(b, p);
    c = mod_ll(c, p);
    long long u = (static_cast<__int128_t>(a) * b + static_cast<__int128_t>(a) * c +
                   static_cast<__int128_t>(b) * c) % p;
    long long v = (a + b + c) % p;
    long long abc = static_cast<long long>((static_cast<__int128_t>(a) * b % p) * c % p);
    long long A = square_mod_ll(v, p);
    long long B = mod_ll(2 * u * v - 4 * abc, p);
    long long C = square_mod_ll(u, p);

    unsigned char partition_mask = 0;
    for (long long d = 0; d < p; ++d) {
        if (!(A == 0 && B == 0)) {
            long long val = mod_ll(A * d * d + B * d + C, p);
            if (val != 0) continue;
        }
        unsigned short bmask = boundary_mask_for_residues(a, b, c, d, p);
        if ((bmask & component_mask) == 0) continue;
        partition_mask |= a2244_partition_mask_for_residues(a, b, c, d, p);
    }
    return partition_mask;
}

static std::vector<unsigned char> build_prime_partition_masks(long long p, unsigned short component_mask) {
    std::vector<unsigned char> masks(static_cast<size_t>(p) * static_cast<size_t>(p) * static_cast<size_t>(p), 0);
    for (long long a = 0; a < p; ++a) {
        for (long long b = 0; b < p; ++b) {
            for (long long c = 0; c < p; ++c) {
                masks[residue_index_prime(a, b, c, p)] =
                    k3_partition_mask_for_abc(a, b, c, p, component_mask);
            }
        }
    }
    return masks;
}


static unsigned short k3_boundary_mask_for_abc(long long a, long long b, long long c, long long p) {
    a = mod_ll(a, p);
    b = mod_ll(b, p);
    c = mod_ll(c, p);
    long long u = (static_cast<__int128_t>(a) * b + static_cast<__int128_t>(a) * c +
                   static_cast<__int128_t>(b) * c) % p;
    long long v = (a + b + c) % p;
    long long abc = static_cast<long long>((static_cast<__int128_t>(a) * b % p) * c % p);
    long long A = square_mod_ll(v, p);
    long long B = mod_ll(2 * u * v - 4 * abc, p);
    long long C = square_mod_ll(u, p);

    // In this degenerate affine chart, d is not determined modulo p.
    // Keep all boundary components possible so the prefilter cannot reject
    // a rational point because of a bad choice of chart.
    if (A == 0 && B == 0) return static_cast<unsigned short>((1u << 10) - 1);

    unsigned short mask = 0;
    for (long long d = 0; d < p; ++d) {
        long long val = mod_ll(A * d * d + B * d + C, p);
        if (val == 0) mask |= boundary_mask_for_residues(a, b, c, d, p);
    }
    return mask;
}

static std::vector<unsigned short> build_prime_boundary_masks(long long p) {
    std::vector<unsigned short> masks(static_cast<size_t>(p) * static_cast<size_t>(p) * static_cast<size_t>(p), 0);
    for (long long a = 0; a < p; ++a) {
        for (long long b = 0; b < p; ++b) {
            for (long long c = 0; c < p; ++c) {
                masks[residue_index_prime(a, b, c, p)] = k3_boundary_mask_for_abc(a, b, c, p);
            }
        }
    }
    return masks;
}

static std::vector<BoundaryRequirement> make_boundary_requirements(
        const std::vector<long long>& primes, const std::string& spec) {
    std::vector<BoundaryRequirement> reqs;
    for (long long p : primes) {
        BoundaryRequirement req;
        req.p = p;
        req.component_mask = all_boundary_components_mask();
        req.label = "any";
        req.masks = build_prime_boundary_masks(p);
        reqs.push_back(req);
    }

    if (!spec.empty()) {
        size_t start = 0;
        while (start <= spec.size()) {
            size_t comma = spec.find(',', start);
            std::string part = spec.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
            if (!part.empty()) {
                size_t colon = part.find(':');
                if (colon == std::string::npos) {
                    std::cerr << "bad boundary spec part: " << part << "\n";
                    std::exit(2);
                }
                long long p = std::atoll(part.substr(0, colon).c_str());
                std::string component_raw = part.substr(colon + 1);
                unsigned short component_mask = boundary_component_mask(component_raw);
                std::string component_label = boundary_component_label(component_raw);
                bool found = false;
                for (BoundaryRequirement& req : reqs) {
                    if (req.p == p) {
                        req.component_mask = component_mask;
                        req.label = component_label;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    BoundaryRequirement req;
                    req.p = p;
                    req.component_mask = component_mask;
                    req.label = component_label;
                    req.masks = build_prime_boundary_masks(p);
                    reqs.push_back(req);
                }
            }
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
    }
    for (BoundaryRequirement& req : reqs) {
        req.masks = build_prime_boundary_masks(req.p);
        req.partition_masks = build_prime_partition_masks(req.p, req.component_mask);
    }
    return reqs;
}

static std::string boundary_requirements_label(const std::vector<BoundaryRequirement>& reqs) {
    std::ostringstream os;
    for (size_t i = 0; i < reqs.size(); ++i) {
        if (i) os << ",";
        os << reqs[i].p << ":";
        os << reqs[i].label;
    }
    return os.str();
}

static size_t residue_index(long long a, long long b, long long c, long long modulus) {
    return (static_cast<size_t>(a) * static_cast<size_t>(modulus) + static_cast<size_t>(b)) *
           static_cast<size_t>(modulus) + static_cast<size_t>(c);
}

static size_t c_residue_list_index(int sc_index, long long a, long long b, long long modulus) {
    return (static_cast<size_t>(sc_index) * static_cast<size_t>(modulus) + static_cast<size_t>(a)) *
           static_cast<size_t>(modulus) + static_cast<size_t>(b);
}

static std::vector<char> build_required_bad_residue_table(
        long long modulus, const std::vector<BoundaryRequirement>& requirements,
        bool require_common_partition_modp, size_t& allowed_count) {
    allowed_count = 0;
    std::vector<char> table(static_cast<size_t>(modulus) * static_cast<size_t>(modulus) *
                            static_cast<size_t>(modulus), 0);
    for (long long ar = 0; ar < modulus; ++ar) {
        for (long long br = 0; br < modulus; ++br) {
            for (long long cr = 0; cr < modulus; ++cr) {
                bool ok = true;
                unsigned char common_partition_mask = static_cast<unsigned char>(0x7);
                for (const BoundaryRequirement& req : requirements) {
                    long long a = ar % req.p;
                    long long b = br % req.p;
                    long long c = cr % req.p;
                    size_t idx = residue_index_prime(a, b, c, req.p);
                    unsigned short mask = req.masks[idx];
                    if ((mask & req.component_mask) == 0) {
                        ok = false;
                        break;
                    }
                    if (require_common_partition_modp) {
                        common_partition_mask &= req.partition_masks[idx];
                        if (common_partition_mask == 0) {
                            ok = false;
                            break;
                        }
                    }
                }
                if (ok) {
                    table[residue_index(ar, br, cr, modulus)] = 1;
                    ++allowed_count;
                }
            }
        }
    }
    return table;
}

static std::vector<std::vector<int>> build_allowed_c_residue_lists(
        long long modulus, const std::vector<char>& residue_table, size_t& total_entries) {
    total_entries = 0;
    std::vector<std::vector<int>> lists(static_cast<size_t>(2) * static_cast<size_t>(modulus) *
                                        static_cast<size_t>(modulus));
    for (long long ar = 0; ar < modulus; ++ar) {
        for (long long br = 0; br < modulus; ++br) {
            for (int sc_index = 0; sc_index < 2; ++sc_index) {
                std::vector<int>& list = lists[c_residue_list_index(sc_index, ar, br, modulus)];
                for (long long cr_abs = 0; cr_abs < modulus; ++cr_abs) {
                    long long cr_signed = (sc_index == 0) ? mod_ll(-cr_abs, modulus) : cr_abs;
                    if (residue_table[residue_index(ar, br, cr_signed, modulus)]) {
                        list.push_back(static_cast<int>(cr_abs));
                    }
                }
                total_entries += list.size();
            }
        }
    }
    return lists;
}

static long long first_with_residue_at_least(long long lo, long long residue, long long modulus) {
    long long delta = mod_ll(residue - mod_ll(lo, modulus), modulus);
    return lo + delta;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " MAX_ABS OUTPUT_FILE [REQUIRED_BAD_PRIMES] [AX_START AX_END [BOUNDARY_SPEC [REAL_LOCAL_DEPTH]]]\n";
        return 2;
    }

    long long B = std::atoll(argv[1]);
    std::string output = argv[2];
    std::vector<long long> required_bad_primes;
    if (argc >= 4) required_bad_primes = parse_prime_list(argv[3]);
    long long ax_start = 1;
    long long ax_end = B;
    if (argc >= 6) {
        ax_start = std::max(1LL, std::atoll(argv[4]));
        ax_end = std::min(B, std::atoll(argv[5]));
    }
    std::string boundary_spec = (argc >= 7) ? argv[6] : std::string();
    int real_local_depth = (argc >= 8) ? std::atoi(argv[7]) : 0;
    std::vector<BoundaryRequirement> boundary_requirements = make_boundary_requirements(required_bad_primes, boundary_spec);

    std::ofstream out(output);
    if (!out) {
        std::cerr << "could not open output file: " << output << "\n";
        return 1;
    }

    std::unordered_set<std::string> seen;
    std::vector<std::string> rows;
    std::vector<i128> curve;

    const long long square_modulus = 3LL * 5LL * 7LL * 11LL * 13LL * 17LL * 19LL;
    const auto square_residue = square_residue_table(square_modulus);
    const int signs[2] = {-1, 1};

    long long required_bad_modulus = 1;
    for (long long p : required_bad_primes) required_bad_modulus *= p;
    std::vector<std::vector<int>> allowed_c_residue_lists;
    if (!required_bad_primes.empty()) {
        size_t allowed_count = 0;
        bool require_common_partition_modp = real_local_depth > 0;
        std::vector<char> residue_table = build_required_bad_residue_table(
            required_bad_modulus, boundary_requirements, require_common_partition_modp, allowed_count);
        size_t total_c_entries = 0;
        allowed_c_residue_lists = build_allowed_c_residue_lists(
            required_bad_modulus, residue_table, total_c_entries);
        std::cerr << "boundary requirements: " << boundary_requirements_label(boundary_requirements) << "\n";
        std::cerr << "allowed signed residue triples modulo " << required_bad_modulus
                  << ": " << allowed_count << "/"
                  << static_cast<size_t>(required_bad_modulus) * static_cast<size_t>(required_bad_modulus) *
                     static_cast<size_t>(required_bad_modulus) << "\n";
        std::cerr << "allowed positive c residue entries: " << total_c_entries << "/"
                  << static_cast<size_t>(2) * static_cast<size_t>(required_bad_modulus) *
                     static_cast<size_t>(required_bad_modulus) * static_cast<size_t>(required_bad_modulus)
                  << "\n";
        if (require_common_partition_modp) {
            std::cerr << "early common A2244 partition mod p prefilter: enabled\n";
        }
    }

    auto process_triple = [&](long long ax, long long bx, long long cx,
                              int sa, int sb, int sc) {
        i128 a = static_cast<i128>(sa) * ax;
        i128 b = static_cast<i128>(sb) * bx;
        i128 c = static_cast<i128>(sc) * cx;

        i128 u = a*b + a*c + b*c;
        i128 v = a + b + c;
        i128 abc = a*b*c;

        if (v == 0) {
            i128 den = 4 * abc;
            if (make_curve_tuple(a, b, c, u*u, den, B, curve)) {
                maybe_add_curve(curve, required_bad_primes, real_local_depth, seen, rows);
            }
            return;
        }

        i128 disc = 16 * abc * (abc - u*v);
        long long disc_mod = mod128(disc, square_modulus);
        if (!square_residue[static_cast<size_t>(disc_mod)]) return;

        i128 sqrt_disc = 0;
        if (!is_square128(disc, sqrt_disc)) return;

        i128 den = 2 * v * v;
        i128 base_num = 4 * abc - 2 * u * v;
        for (int eps : signs) {
            i128 num = base_num + static_cast<i128>(eps) * sqrt_disc;
            if (make_curve_tuple(a, b, c, num, den, B, curve)) {
                maybe_add_curve(curve, required_bad_primes, real_local_depth, seen, rows);
            }
        }
    };

    for (long long ax = ax_start; ax <= ax_end; ++ax) {
        for (long long bx = ax + 1; bx <= B; ++bx) {
            if (!required_bad_primes.empty()) {
                for (int sa_index = 0; sa_index < 2; ++sa_index) {
                    int sa = signs[sa_index];
                    long long ar = mod_ll(static_cast<long long>(sa) * ax, required_bad_modulus);
                    for (int sb_index = 0; sb_index < 2; ++sb_index) {
                        int sb = signs[sb_index];
                        long long br = mod_ll(static_cast<long long>(sb) * bx, required_bad_modulus);
                        for (int sc_index = 0; sc_index < 2; ++sc_index) {
                            int sc = signs[sc_index];
                            const std::vector<int>& c_residues =
                                allowed_c_residue_lists[c_residue_list_index(sc_index, ar, br, required_bad_modulus)];
                            for (int cr : c_residues) {
                                for (long long cx = first_with_residue_at_least(bx + 1, cr, required_bad_modulus);
                                     cx <= B; cx += required_bad_modulus) {
                                    process_triple(ax, bx, cx, sa, sb, sc);
                                }
                            }
                        }
                    }
                }
            } else {
                for (long long cx = bx + 1; cx <= B; ++cx) {
                    for (int sa : signs) {
                        for (int sb : signs) {
                            for (int sc : signs) {
                                process_triple(ax, bx, cx, sa, sb, sc);
                            }
                        }
                    }
                }
            }
        }
    }

    std::sort(rows.begin(), rows.end());
    for (const auto& row : rows) out << row;

    std::cerr << "wrote " << rows.size() << " tuples to " << output;
    if (!required_bad_primes.empty()) {
        std::cerr << " with required bad primes";
        for (long long p : required_bad_primes) std::cerr << " " << p;
    }
    if (argc >= 6) std::cerr << " ax_range " << ax_start << ".." << ax_end;
    if (!boundary_spec.empty()) std::cerr << " boundary_spec " << boundary_spec;
    if (real_local_depth > 0) std::cerr << " real_local_depth " << real_local_depth;
    std::cerr << "\n";
    return 0;
}
