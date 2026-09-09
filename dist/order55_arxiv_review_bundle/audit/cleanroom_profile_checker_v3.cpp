// Clean-room order-55 profile evaluator.  It does not call repository executables.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
constexpr int N = 55, HALF = 27;
constexpr u64 P = 2305843009213696591ULL, Q = 2305843009213697141ULL;
constexpr u128 INCUMBENT = (u128)13469409409475839533ULL * 10000000000000ULL
                          + (u128)1307111329132ULL;

u64 mul(u64 a, u64 b, u64 p) { return (u128)a * b % p; }
u64 powmod(u64 a, u64 e, u64 p) {
    u64 r = 1;
    while (e) {
        if (e & 1) r = mul(r, a, p);
        a = mul(a, a, p);
        e >>= 1;
    }
    return r;
}
u64 primitive_root(u64 p) {
    for (u64 c = 2;; ++c) {
        u64 r = powmod(c, (p - 1) / N, p);
        if (powmod(r, 5, p) != 1 && powmod(r, 11, p) != 1) return r;
    }
}
std::string decimal(u128 x) {
    if (!x) return "0";
    std::string out;
    while (x) { out.push_back(char('0' + x % 10)); x /= 10; }
    std::reverse(out.begin(), out.end());
    return out;
}
std::vector<int> parse(const std::string& line) {
    std::istringstream in(line);
    std::vector<int> row;
    int x;
    while (in >> x) row.push_back(x);
    return row;
}
std::set<std::vector<int>> read_rows(const std::string& path) {
    std::ifstream in(path);
    std::set<std::vector<int>> rows;
    std::string line;
    while (std::getline(in, line)) {
        auto row = parse(line);
        if (!row.empty()) rows.insert(std::move(row));
    }
    return rows;
}
bool unit(int x) { return std::gcd(x, N) == 1; }

// The candidate is canonical iff no unit multiplier gives a lexicographically
// smaller paired correlation vector.  This direction is deliberately checked
// independently of the production affine_canonical implementation.
bool canonical(const std::vector<int>& candidate, const std::array<int, N>& full) {
    for (int u = 1; u < N; ++u) {
        if (!unit(u)) continue;
        for (int s = 1; s <= HALF; ++s) {
            const int image = (u * s) % N;
            const int paired_image = std::min(image, N - image);
            if (full[paired_image] < candidate[s - 1]) return false;
            if (full[paired_image] > candidate[s - 1]) break;
        }
    }
    return true;
}
struct Counts { u64 ordered = 0, folded = 0, canonical_profiles = 0, above = 0; };

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const std::string cases_path = argv[1], cert = argv[2], output_path = argv[3];
    std::array<std::array<u64, HALF + 1>, HALF + 1> tp{}, tq{};
    const u64 rp = primitive_root(P), rq = primitive_root(Q);
    for (int j = 1; j <= HALF; ++j) for (int s = 1; s <= HALF; ++s) {
        const int e = j * s % N;
        tp[j][s] = (powmod(rp, e, P) + powmod(rp, (N - e) % N, P)) % P;
        tq[j][s] = (powmod(rq, e, Q) + powmod(rq, (N - e) % N, Q)) % Q;
    }
    const u64 pinv = powmod(P, Q - 2, Q);
    std::array<std::set<std::vector<int>>, 28> fold5, fold11;
    for (int k = 25; k <= 27; ++k) {
        fold5[k] = read_rows(cert + "/foldcorr_m5_k" + std::to_string(k) + ".txt");
        fold11[k] = read_rows(cert + "/foldcorr_m11_k" + std::to_string(k) + ".txt");
    }
    std::ifstream in(cases_path);
    std::ofstream out(output_path);
    if (!in || !out) return 3;
    std::array<Counts, 28> counts{};
    std::string line;
    while (std::getline(in, line)) {
        auto fields = parse(line);
        if (fields.size() != HALF + 2) return 4;
        const int case_id = fields[0], k = fields[1];
        std::vector<int> profile(fields.begin() + 2, fields.end());
        std::sort(profile.begin(), profile.end());
        do {
            ++counts[k].ordered;
            std::array<int, N> full{};
            full[0] = k;
            for (int s = 1; s <= HALF; ++s) {
                full[s] = profile[s - 1];
                full[N - s] = profile[s - 1];
            }
            std::vector<int> h5(5), h11(11);
            for (int s = 0; s < N; ++s) {
                h5[s % 5] += full[s];
                h11[s % 11] += full[s];
            }
            h5.resize(3); h11.resize(6);
            if (!fold5[k].count(h5) || !fold11[k].count(h11)) continue;
            ++counts[k].folded;
            if (!canonical(profile, full)) continue;
            ++counts[k].canonical_profiles;
            const int base = k * (k - 1) / 54;
            u64 residues[2]{};
            for (int field = 0; field < 2; ++field) {
                const u64 prime = field == 0 ? P : Q;
                const auto& table = field == 0 ? tp : tq;
                u64 product = (55 - k) % prime;
                for (int j = 1; j <= HALF; ++j) {
                    u64 spectral = (u64)(k - base) % prime;
                    for (int s = 1; s <= HALF; ++s) {
                        const int d = profile[s - 1] - base;
                        const u64 term = mul((u64)(d < 0 ? -d : d), table[j][s], prime);
                        spectral = d < 0 ? (spectral + prime - term) % prime
                                         : (spectral + term) % prime;
                    }
                    product = mul(product, spectral, prime);
                }
                residues[field] = product;
            }
            const u64 delta = (residues[1] + Q - residues[0]) % Q;
            u128 value = residues[0] + (u128)P * mul(delta, pinv, Q);
            if (value > (u128)P * Q / 2) value = (u128)P * Q - value;
            if (value >= INCUMBENT) ++counts[k].above;
            out << case_id << ' ' << k << ' ' << decimal(value);
            for (int x : profile) out << ' ' << x;
            out << '\n';
        } while (std::next_permutation(profile.begin(), profile.end()));
    }
    for (int k = 25; k <= 27; ++k)
        std::cerr << k << ' ' << counts[k].ordered << ' ' << counts[k].folded << ' '
                  << counts[k].canonical_profiles << ' ' << counts[k].above << '\n';
    return 0;
}
