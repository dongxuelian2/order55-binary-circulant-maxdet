"""Clean-room audit of the committed order-55 certificate.

The production verifier is intentionally not imported.  This checker reads
only the frozen certificate and re-derives the mathematical maps in Python:

* exact rational first/fourth-moment bounds;
* the 5/11 folded margin survivors and folded-correlation signatures;
* the sorted correlation-profile universe and multinomial counts;
* the signed two-prime product for formal correlation profiles;
* the target-to-margin-task map and every task's input/output topology.

The expensive binary joins are not silently represented as recomputed here:
their recorded per-task counters are checked for topology and consistency,
while the independent full production replay and a separate clean-room
winner verifier are reported as distinct evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter, defaultdict
from fractions import Fraction
from math import factorial, gcd, isqrt, prod
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
CERT = ROOT / "certificates" / "order55_global"
N = 55
HALF = 27
M = 134694094094758395331307111329132
P, Q = 2305843009213696591, 2305843009213697141
SCALE = 1 << 24
BOUND_SCALE = 10**80
MASK55 = (1 << 55) - 1
UNITS55 = tuple(u for u in range(1, 28) if gcd(u, 55) == 1)


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def jsonl(path: Path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def deterministic_prime(value: int) -> bool:
    """Deterministic Miller-Rabin for unsigned 64-bit inputs."""

    if value < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for divisor in small:
        if value % divisor == 0:
            return value == divisor
    exponent = value - 1
    twos = 0
    while exponent % 2 == 0:
        twos += 1
        exponent //= 2
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        x = pow(base, exponent, value)
        if x in (1, value - 1):
            continue
        for _ in range(twos - 1):
            x = x * x % value
            if x == value - 1:
                break
        else:
            return False
    return True


def sqrt_floor(x: Fraction) -> Fraction:
    x = Fraction(x)
    if x < 0:
        raise ValueError("negative radicand")
    numerator = x.numerator * BOUND_SCALE * BOUND_SCALE
    denominator = x.denominator
    return Fraction(isqrt(numerator // denominator), BOUND_SCALE)


def sqrt_ceil(x: Fraction) -> Fraction:
    x = Fraction(x)
    lo = sqrt_floor(x)
    if lo * lo == x:
        return lo
    return lo + Fraction(1, BOUND_SCALE)


def balanced(total: int, count: int) -> int:
    quotient, remainder = divmod(total, count)
    return (count - remainder) * quotient * quotient + remainder * (quotient + 1) ** 2


def product_upper(total, square_lower, count, cap=None) -> Fraction:
    """Independent exact-rational upper bound for a capped moment problem."""

    total, square_lower = Fraction(total), Fraction(square_lower)
    cap = None if cap is None else Fraction(cap)
    square_lower = max(square_lower, total * total / count)
    endpoint_counts = (range(count) if cap is not None else (0,))
    best = Fraction(0)
    for endpoints in endpoint_counts:
        remaining = count - endpoints
        residual_sum = total - endpoints * cap if cap is not None else total
        residual_sq = (square_lower - endpoints * cap * cap
                       if cap is not None else square_lower)
        if residual_sum < 0 or residual_sq < residual_sum * residual_sum / remaining:
            continue
        if cap is not None and (residual_sum > remaining * cap
                                or residual_sq > residual_sum * cap):
            continue
        mean = residual_sum / remaining
        variance = residual_sq - mean * mean * remaining
        multiplier = cap ** endpoints if cap is not None else Fraction(1)
        if variance == 0:
            best = max(best, multiplier * mean ** remaining)
            continue
        for low_count in range(1, remaining):
            low_var = variance * Fraction(remaining - low_count,
                                          remaining * low_count)
            high_var = variance * Fraction(low_count,
                                           remaining * (remaining - low_count))
            if low_var >= mean * mean:
                continue
            low = mean - sqrt_floor(low_var)
            high = mean + sqrt_ceil(high_var)
            # A lower bound on the high root exceeding the cap proves
            # infeasibility; using the lower root here is intentional.
            if cap is not None and mean + sqrt_floor(high_var) > cap:
                continue
            best = max(best, multiplier * low ** low_count
                       * high ** (remaining - low_count))
    return best


def fourth_bound(k: int, square_sum: int) -> Fraction:
    spectral_sum = Fraction(k * (55 - k), 2)
    spectral_sq = Fraction(55 * (k * k + 2 * square_sum) - k**4, 2)
    return (55 - k) * product_upper(spectral_sum, spectral_sq, HALF)


def cosine_cap(k: int, values: tuple[int, ...], cos_intervals) -> Fraction:
    base = k * (k - 1) // 54
    deviations = sorted(value - base for value in values)
    candidates = []
    for unit_class in (1, 5, 11):
        coefficients = sorted(
            cos_intervals[min(unit_class * shift % 55,
                              55 - unit_class * shift % 55)]
            for shift in range(1, 28)
        )
        total = Fraction(k - base)
        for deviation, interval in zip(deviations, coefficients):
            total += deviation * (interval[1] if deviation >= 0 else interval[0]) / SCALE
        candidates.append(total)
    return max(candidates)


def profile_key(values) -> tuple[int, ...]:
    return tuple(int(value) for value in values)


def unit_canonical(values: tuple[int, ...]) -> tuple[int, ...]:
    transforms = []
    for unit in UNITS55:
        transforms.append(tuple(values[min(unit * shift % 55,
                                          55 - unit * shift % 55) - 1]
                               for shift in range(1, 28)))
    return min(transforms)


def fold_signature(values: tuple[int, ...], k: int, modulus: int) -> tuple[int, ...]:
    result = [0] * ((modulus + 1) // 2)
    result[0] = k
    for shift, value in enumerate(values, start=1):
        residue = shift % modulus
        # A paired shift contributes to two different residue classes unless
        # it is divisible by the folding modulus.  Only the latter case puts
        # both s and -s in residue zero.
        result[min(residue, modulus - residue)] += value * (2 if residue == 0 else 1)
    return tuple(result)


def affine_canonical_small(values: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    return min(
        tuple(values[(unit * index + translation) % modulus]
              for index in range(modulus))
        for unit in range(1, modulus)
        if gcd(unit, modulus) == 1
        for translation in range(modulus)
    )


def rotations(values: tuple[int, ...]) -> set[tuple[int, ...]]:
    return {values[offset:] + values[:offset] for offset in range(len(values))}


def margin_norm(values: tuple[int, ...], modulus: int) -> int:
    x = sp.Symbol("x")
    polynomial = sum(value * x**index for index, value in enumerate(values))
    return int(sp.resultant(polynomial, sp.cyclotomic_poly(modulus, x), x))


def margin_vectors(modulus: int, weight: int, cap: int):
    """Enumerate one representative for every translation orbit, then filter."""

    vector = [0] * modulus

    def visit(position: int, remaining: int):
        if position == modulus:
            if remaining == 0:
                yield tuple(vector)
            return
        low = vector[0] if position else 0
        slots = modulus - position - 1
        for value in range(low, min(cap, remaining) + 1):
            if remaining - value > slots * cap:
                continue
            vector[position] = value
            yield from visit(position + 1, remaining - value)

    yield from visit(0, weight)


def formal_profile_product(k: int, values: tuple[int, ...]) -> int:
    """Signed CRT product for a formal half-correlation profile."""

    residues = []
    base = k * (k - 1) // 54
    for modulus in (P, Q):
        root = None
        for generator in range(2, 1000):
            candidate = pow(generator, (modulus - 1) // 55, modulus)
            if (pow(candidate, 55, modulus) == 1
                    and pow(candidate, 5, modulus) != 1
                    and pow(candidate, 11, modulus) != 1):
                root = candidate
                break
        if root is None:
            raise AssertionError("primitive root search failed")
        result = (55 - k) % modulus
        for frequency in range(1, 28):
            q = (k - base) % modulus
            for shift, value in enumerate(values, start=1):
                deviation = value - base
                coefficient = (pow(root, frequency * shift, modulus)
                               + pow(root, (-frequency * shift) % 55, modulus)) % modulus
                q = (q + deviation * coefficient) % modulus
            result = result * q % modulus
        residues.append(result)
    inverse = pow(P, -1, Q)
    delta = (residues[1] - residues[0]) % Q
    reconstructed = residues[0] + P * ((delta * inverse) % Q)
    modulus_product = P * Q
    if reconstructed > modulus_product // 2:
        reconstructed = modulus_product - reconstructed
    return reconstructed


def half_correlations(mask: int) -> tuple[int, ...]:
    return tuple((mask & (((mask << shift) | (mask >> (55 - shift))) & MASK55)).bit_count()
                 for shift in range(1, 28))


def word_fold(mask: int, modulus: int) -> tuple[int, ...]:
    return tuple(sum((mask >> index) & 1 for index in range(residue, 55, modulus))
                 for residue in range(modulus))


def verify_hashes() -> dict[str, int]:
    global_data = load(CERT / "global.json")
    manifest = load(CERT / "hash_manifest.json")
    if digest(CERT / "hash_manifest.json") != global_data["hash_manifest_sha256"]:
        raise AssertionError("global hash_manifest_sha256 mismatch")
    if manifest != global_data["components"]:
        raise AssertionError("global components do not equal hash manifest")
    for name, expected in manifest.items():
        path = CERT / name
        if path.name != name or digest(path) != expected:
            raise AssertionError(("certificate hash mismatch", name))
    source_hashes = load(CERT / "source_hashes.json")
    for name, expected in source_hashes.items():
        if digest(ROOT / name) != expected:
            raise AssertionError(("source hash mismatch", name))
    if global_data["status"] != "COMPLETE" or global_data["order"] != 55:
        raise AssertionError("incomplete or wrong-order certificate")
    if int(global_data["maximum"]) != M:
        raise AssertionError("claimed maximum changed")
    if not all(deterministic_prime(p) and (p - 1) % 55 == 0 for p in (P, Q)):
        raise AssertionError("production CRT prime check failed independently")
    return {
        "certificate_files_hashed": len(manifest),
        "source_files_hashed": len(source_hashes),
    }


def verify_weights(global_data) -> tuple[dict[int, int], dict[int, int]]:
    rows = load(CERT / "weight_screen.json")
    if [row["k"] for row in rows] != list(range(1, 28)):
        raise AssertionError("weight rows are not exactly 1..27")
    active = {}
    floors = {}
    for row in rows:
        k = row["k"]
        low = balanced(k * (k - 1) // 2, HALF)
        first_moment = Fraction(55 - k) * Fraction(k * (55 - k), 54) ** HALF
        fourth = fourth_bound(k, low)
        bound = min(first_moment, fourth)
        floor_value = int(bound)
        floors[k] = floor_value
        if floor_value != row["upper_floor"]:
            raise AssertionError(("weight bound mismatch", k, floor_value, row["upper_floor"]))
        if row["complement_weight"] != 55 - k or row["excluded"] != (floor_value < M):
            raise AssertionError(("weight row metadata mismatch", k))
        if not row["excluded"]:
            t = low
            while int(fourth_bound(k, t + 1)) >= M:
                t += 1
            if row["minimum_half_squares"] != low or row["maximum_half_squares"] != t:
                raise AssertionError(("weight limit mismatch", k))
            if row["first_excluded_upper"] != int(fourth_bound(k, t + 1)):
                raise AssertionError(("first excluded bound mismatch", k))
            active[k] = t
    if sorted(active) != [25, 26, 27] or global_data["active_weights"] != sorted(active):
        raise AssertionError(("active weights mismatch", active))
    return active, floors


def verify_cosines():
    rows = [tuple(map(int, line.split()))
            for line in (CERT / "cosine_intervals.txt").read_text().splitlines()]
    if len(rows) != 28:
        raise AssertionError("wrong cosine interval count")
    # Machin bounds and a longer independent Taylor enclosure than the
    # production construction.  Stored entries represent 2*cos(2*pi*r/55).
    def atan_bounds(q: int):
        partial = sum(Fraction((-1) ** j, (2 * j + 1) * q ** (2 * j + 1))
                      for j in range(120))
        next_term = Fraction(1, 241 * q ** 241)
        return partial, partial + next_term

    a, b = atan_bounds(5)
    c, d = atan_bounds(239)
    pi_low, pi_high = 16 * a - 4 * d, 16 * b - 4 * c
    for r, (low, high) in enumerate(rows):
        angle = (pi_low + pi_high) * r / 55
        polynomial = sum(Fraction((-1) ** j) * angle ** (2 * j) / factorial(2 * j)
                         for j in range(50))
        error = angle ** 100 / factorial(100) + (pi_high - pi_low) * r / 55
        if not (Fraction(low, SCALE) <= 2 * (polynomial - error)
                <= 2 * (polynomial + error) <= Fraction(high, SCALE)):
            raise AssertionError(("cosine interval miss", r, low, high))
    return rows


def verify_fold_profiles(active, threshold_rows):
    survivors = {}
    folded_signatures = {}
    for modulus in (5, 11):
        cap = 55 // modulus
        for k in active:
            expected = {tuple(row["profile"])
                        for row in jsonl(CERT / f"profiles_m{modulus}_k{k}.jsonl")}
            actual = set()
            for vector in margin_vectors(modulus, k, cap):
                square_sum = sum(value * value for value in vector)
                need = threshold_rows.get((modulus, k, square_sum))
                if need is None:
                    continue
                canonical = affine_canonical_small(vector, modulus)
                if canonical != vector:
                    continue
                norm = margin_norm(vector, modulus)
                if norm >= need:
                    actual.add(vector)
            if actual != expected:
                raise AssertionError(("fold profile set mismatch", modulus, k,
                                      len(actual), len(expected)))
            signatures = set()
            for vector in actual:
                for unit in range(1, modulus):
                    oriented = tuple(vector[(unit * index) % modulus]
                                     for index in range(modulus))
                    oriented = min(rotations(oriented))
                    h = tuple(sum(oriented[index] * oriented[(index + shift) % modulus]
                                  for index in range(modulus))
                              for shift in range(modulus))
                    signatures.add(h[: (modulus + 1) // 2])
            stored = {tuple(map(int, line.split()))
                      for line in (CERT / f"foldcorr_m{modulus}_k{k}.txt").read_text().splitlines()
                      if line}
            if signatures != stored:
                raise AssertionError(("folded signature mismatch", modulus, k,
                                      len(signatures), len(stored)))
            survivors[k, modulus] = actual
            folded_signatures[k, modulus] = signatures
    return survivors, folded_signatures


def enumerate_multisets(k: int, limit: int) -> set[tuple[int, ...]]:
    result = set()
    values = [0] * HALF

    def visit(position: int, minimum: int, remaining: int, square_sum: int):
        slots = HALF - position
        if slots == 0:
            if remaining == 0:
                result.add(tuple(values))
            return
        if remaining < minimum * slots:
            return
        if square_sum + balanced(remaining, slots) > limit:
            return
        for value in range(minimum, min(k, remaining // slots) + 1):
            values[position] = value
            visit(position + 1, value, remaining - value,
                  square_sum + value * value)

    visit(0, 0, k * (k - 1) // 2, 0)
    return result


def verify_correlation_profiles(active, cos_intervals, survivors, folded_signatures, global_data):
    entries = load(CERT / "correlation_profiles.json")
    by_k = defaultdict(set)
    for entry in entries:
        by_k[entry["k"]].add(tuple(entry["values"]))
    for k, limit in active.items():
        expected = enumerate_multisets(k, limit)
        if by_k[k] != expected:
            raise AssertionError(("multiset universe mismatch", k,
                                  len(expected), len(by_k[k])))

    target_rows = []
    for entry in entries:
        k = entry["k"]
        values = profile_key(entry["values"])
        if sum(value * value for value in values) != entry["squares"]:
            raise AssertionError(("profile square mismatch", k))
        expected_count = factorial(HALF)
        for multiplicity in Counter(values).values():
            expected_count //= factorial(multiplicity)
        if expected_count != entry["ordered_count"]:
            raise AssertionError(("multinomial count mismatch", k, values))
        cap = cosine_cap(k, values, cos_intervals)
        if cap != Fraction(entry["cap"]):
            raise AssertionError(("spectral cap mismatch", k, values, cap, entry["cap"]))
        sq = Fraction(55 * (k * k + 2 * entry["squares"]) - k**4, 2)
        independent_bound = (55 - k) * product_upper(
            Fraction(k * (55 - k), 2), sq, HALF, cap)
        if int(independent_bound) < entry["upper_floor"]:
            raise AssertionError(("independent capped bound is too small", k, values,
                                  int(independent_bound), entry["upper_floor"]))
        if entry["route"] == "boxed moment":
            if int(independent_bound) >= M:
                raise AssertionError(("boxed route not independently excluded", k, values))
            continue
        if entry["route"] != "exhaustive correlation enumeration":
            raise AssertionError(("unknown profile route", entry["route"]))
        source = CERT / f"corr_k{k}_part{entry['partition']}.txt"
        if tuple(map(int, source.read_text().split())) != values:
            raise AssertionError(("partition source mismatch", source))
        audit = entry["audit"]
        if audit["ordered_profiles"] != entry["ordered_count"]:
            raise AssertionError(("audit ordered count mismatch", source))
        targets = jsonl(CERT / f"corr_k{k}_part{entry['partition']}_targets.jsonl")
        if audit["above_screen"] != len(targets):
            raise AssertionError(("target line count mismatch", source))
        for raw_target in targets:
            target = dict(raw_target, source_partition=entry["partition"])
            c = profile_key(target["correlations"])
            if tuple(sorted(c)) != values or unit_canonical(c) != c:
                raise AssertionError(("target profile canonical/source mismatch", target))
            if formal_profile_product(k, c) != int(target["absolute_profile_product"]):
                raise AssertionError(("formal profile product mismatch", target))
            if fold_signature(c, k, 5) not in folded_signatures[k, 5]:
                raise AssertionError(("target failed mod-5 folded screen", target))
            if fold_signature(c, k, 11) not in folded_signatures[k, 11]:
                raise AssertionError(("target failed mod-11 folded screen", target))
            target_rows.append(target)
    target_rows.sort(key=lambda row: (-int(row["absolute_profile_product"]),
                                      row["k"], row["correlations"]))
    if len(target_rows) != global_data["correlation_targets"]:
        raise AssertionError(("target total mismatch", len(target_rows)))
    return entries, target_rows


def verify_lift_topology(active, survivors, target_rows, global_data):
    coverage = load(CERT / "lift_targets.json")
    if len(coverage) != len(target_rows):
        raise AssertionError("lift target count mismatch")
    oriented = {}
    for k in active:
        for modulus in (5, 11):
            lookup = defaultdict(set)
            for vector in survivors[k, modulus]:
                for unit in range(1, modulus):
                    transformed = tuple(vector[(unit * index) % modulus]
                                        for index in range(modulus))
                    normalized = min(rotations(transformed))
                    h = tuple(sum(normalized[index]
                                  * normalized[(index + shift) % modulus]
                                  for index in range(modulus))
                              for shift in range(modulus))
                    lookup[h].add(normalized)
            oriented[k, modulus] = lookup
    target_by_id = {}
    task_vectors = {}
    for index, row in enumerate(coverage):
        if row["id"] != index or row["target"] != target_rows[index]:
            raise AssertionError(("target ordering/id mismatch", index))
        target = row["target"]
        k = target["k"]
        c = profile_key(target["correlations"])
        sig5 = fold_signature(c, k, 5)
        full5 = tuple(sig5[i] if i < 3 else sig5[5 - i] for i in range(5))
        sig11 = fold_signature(c, k, 11)
        full11 = tuple(sig11[i] if i < 6 else sig11[11 - i] for i in range(11))
        row_profiles = sorted(oriented[k, 5].get(full5, set()))
        column_profiles = sorted(oriented[k, 11].get(full11, set()))
        if row_profiles != [tuple(v) for v in row["row_profiles"]]:
            raise AssertionError(("row margin map mismatch", index))
        if column_profiles != [tuple(v) for v in row["column_profiles"]]:
            raise AssertionError(("column margin map mismatch", index))
        expected_tasks = len(row_profiles) * len(column_profiles)
        if row["tasks"] != expected_tasks:
            raise AssertionError(("lift task multiplicity mismatch", index))
        target_by_id[index] = target
        for i, margin5 in enumerate(row_profiles):
            for j, margin11 in enumerate(column_profiles):
                task = f"{index}_{i}_{j}"
                full = list(margin5) + list(margin11) + [k] + list(c) + list(reversed(c))
                task_vectors[task] = full
    if len(task_vectors) != global_data["lift_tasks"]:
        raise AssertionError(("total task count mismatch", len(task_vectors)))

    seen = set()
    part_summaries = []
    output_words = []
    for part in range(8):
        input_path = CERT / f"lift_part{part}.txt"
        audit_path = CERT / f"lift_part{part}_audit.jsonl"
        word_path = CERT / f"lift_part{part}_words.txt"
        input_rows = input_path.read_text().splitlines()
        local_ids = set()
        for line in input_rows:
            fields = line.split()
            task = fields[0]
            vector = list(map(int, fields[1:]))
            if task not in task_vectors or task in seen or vector != task_vectors[task]:
                raise AssertionError(("lift input topology mismatch", part, task))
            seen.add(task)
            local_ids.add(task)
        audits = jsonl(audit_path)
        if {row["task"] for row in audits} != local_ids or len(audits) != len(local_ids):
            raise AssertionError(("lift audit ids mismatch", part))
        words = word_path.read_text().splitlines()
        if len(words) != sum(row["solutions"] for row in audits):
            raise AssertionError(("lift solution count mismatch", part))
        for line in words:
            task, word = line.split()
            if task not in local_ids or len(word) != 55 or set(word) - {"0", "1"}:
                raise AssertionError(("invalid lift output word", part, line))
            vector = task_vectors[task]
            mask = sum((bit == "1") << index for index, bit in enumerate(word))
            if word_fold(mask, 5) != tuple(vector[:5]):
                raise AssertionError(("mod-5 lift output mismatch", task))
            if word_fold(mask, 11) != tuple(vector[5:16]):
                raise AssertionError(("mod-11 lift output mismatch", task))
            if half_correlations(mask) != tuple(vector[17:44]):
                raise AssertionError(("correlation lift output mismatch", task))
            output_words.append((task, word))
        part_summaries.append({
            "part": part,
            "tasks": len(local_ids),
            "joined_words_recorded": sum(row["joined_words"] for row in audits),
            "solutions": len(words),
        })
    if seen != set(task_vectors):
        raise AssertionError(("missing or extra lift tasks", len(seen), len(task_vectors)))
    if sum(row["joined_words_recorded"] for row in part_summaries) != global_data["joined_words"]:
        raise AssertionError("joined-word total mismatch")
    if len(output_words) != 2:
        raise AssertionError(("unexpected normalized solution count", len(output_words)))
    if len({word for _, word in output_words}) != len(output_words):
        raise AssertionError("duplicate normalized lift word")
    return part_summaries, output_words


def random_forward_test(entries, target_rows, count: int):
    route_by_profile = {}
    limit_by_weight = {}
    for entry in entries:
        route_by_profile[entry["k"], tuple(entry["values"])] = entry["route"]
        limit_by_weight[entry["k"]] = max(limit_by_weight.get(entry["k"], 0),
                                           entry["squares"])
    target_profiles = {
        (target["k"], tuple(target["correlations"]))
        for target in target_rows
    }
    rng = random.Random(550055)
    categories = Counter()
    for _ in range(count):
        mask = rng.getrandbits(55)
        weight = mask.bit_count()
        if weight > 27:
            mask ^= MASK55
            weight = 55 - weight
        if weight == 0:
            categories["constant/singular"] += 1
            continue
        if weight <= 24:
            categories["weight-screen"] += 1
            continue
        profile = tuple(sorted(half_correlations(mask)))
        if sum(value * value for value in profile) > limit_by_weight[weight]:
            categories["fourth-moment"] += 1
            continue
        route = route_by_profile.get((weight, profile))
        if route is None:
            raise AssertionError(("random UNCLASSIFIED profile", weight, profile))
        if route == "boxed moment":
            categories["boxed-moment"] += 1
        else:
            canonical = unit_canonical(half_correlations(mask))
            if (weight, canonical) in target_profiles:
                categories["explicit-lift-target"] += 1
            else:
                categories["explicit-profile-screen"] += 1
    return categories


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--random-count", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("global_audit.json"))
    args = parser.parse_args()
    if args.random_count < 1:
        raise SystemExit("random count must be positive")

    global_data = load(CERT / "global.json")
    hash_summary = verify_hashes()
    active, floors = verify_weights(global_data)
    cos = verify_cosines()
    threshold_rows = {
        tuple(map(int, line.split()[:3])): int(line.split()[3])
        for line in (CERT / "profile_thresholds.txt").read_text().splitlines()
        if line
    }
    survivors, folded = verify_fold_profiles(active, threshold_rows)
    entries, targets = verify_correlation_profiles(
        active, cos, survivors, folded, global_data)
    parts, output_words = verify_lift_topology(active, survivors, targets, global_data)
    random_categories = random_forward_test(entries, targets, args.random_count)
    summary = {
        "status": "PASS",
        "hashes": hash_summary,
        "weight_floors": floors,
        "active_limits": active,
        "fold_survivor_counts": {
            f"k{k}_m{m}": len(survivors[k, m])
            for k in active for m in (5, 11)
        },
        "correlation_profile_entries": len(entries),
        "correlation_targets": len(targets),
        "lift_parts": parts,
        "normalized_output_words": output_words,
        "random_count": args.random_count,
        "random_categories": dict(random_categories),
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": "PASS",
        "correlation_targets": len(targets),
        "lift_tasks": global_data["lift_tasks"],
        "random_count": args.random_count,
        "random_categories": dict(random_categories),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
