#!/usr/bin/env python3
"""Shared Phase-0 audit helpers."""

from __future__ import annotations

from array import array
from collections import Counter
from typing import Any

from sympy import factorint


SMALL_PRIMORIAL_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29)


def blocker_weight(ell: int) -> int:
    """Return the Phase-0 Q1 blocker multiplicity k(ell).

    Current review convention:
    k(5)=3; k(7)=k(11)=k(13)=k(17)=2; k(ell)=1 for ell>=19.
    """

    if ell == 5:
        return 3
    if ell in {7, 11, 13, 17}:
        return 2
    if ell >= 19:
        return 1
    return 0


def primes_up_to(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    limit = int(n**0.5)
    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def omega_ge5_exact_sieve(y_start: int, y_end_exclusive: int, prime_limit: int) -> dict[str, Any]:
    """Count max omega_{>=5}(y+1) exactly for y in [y_start, y_end_exclusive).

    The sieve divides every n=y+1 by all primes 5 <= ell <= prime_limit and
    counts distinct divisors. A residual > 1 after that contributes one more
    large prime factor. For Segment A intervals, prime_limit=q covers sqrt(n).
    """

    if y_end_exclusive <= y_start:
        raise ValueError("empty y interval")

    n_start = y_start + 1
    n_stop = y_end_exclusive + 1
    length = n_stop - n_start
    residuals = array("Q", range(n_start, n_stop))
    counts = bytearray(length)
    for idx, value in enumerate(residuals):
        while value % 2 == 0:
            value //= 2
        while value % 3 == 0:
            value //= 3
        residuals[idx] = value

    for ell in primes_up_to(prime_limit):
        if ell < 5:
            continue
        first = ((n_start + ell - 1) // ell) * ell
        for n in range(first, n_stop, ell):
            idx = n - n_start
            if residuals[idx] % ell == 0:
                counts[idx] += 1
                while residuals[idx] % ell == 0:
                    residuals[idx] //= ell

    for idx, residual in enumerate(residuals):
        if residual >= 5:
            counts[idx] += 1

    hist: Counter[int] = Counter(int(c) for c in counts)
    max_count = max(hist)
    first_idx = next(i for i, c in enumerate(counts) if c == max_count)
    n = n_start + first_idx
    factors = sorted(int(p) for p in factorint(n))
    factors_ge5 = [p for p in factors if p >= 5]
    if len(factors_ge5) != max_count:
        raise AssertionError(
            f"omega_ge5 sieve mismatch for n={n}: sieve={max_count}, factorint={len(factors_ge5)}"
        )
    return {
        "mode": "exact_segmented_sieve",
        "y_start": y_start,
        "y_end_exclusive": y_end_exclusive,
        "n_start": n_start,
        "n_stop_exclusive": n_stop,
        "rows_scanned": length,
        "prime_limit": prime_limit,
        "max_omega_ge5": {
            "y": n - 1,
            "n": n,
            "omega_ge5": max_count,
            "n_factors": factors,
            "n_factors_ge5": factors_ge5,
        },
        "omega_ge5_hist": {str(k): v for k, v in sorted(hist.items())},
    }


def weighted_load_from_factors(factors_ge5: list[int]) -> float:
    return sum(blocker_weight(ell) / ell for ell in factors_ge5)


def weighted_q1_exact_sieve(y_start: int, y_end_exclusive: int, prime_limit: int) -> dict[str, Any]:
    """Compute max weighted Q1 load L(n)=sum_{ell|n, ell>=5} k(ell)/ell."""

    if y_end_exclusive <= y_start:
        raise ValueError("empty y interval")

    n_start = y_start + 1
    n_stop = y_end_exclusive + 1
    length = n_stop - n_start
    residuals = array("Q", range(n_start, n_stop))
    loads = [0.0] * length
    small_masks = bytearray(length)
    small_index = {ell: i for i, ell in enumerate(SMALL_PRIMORIAL_PRIMES)}

    for idx, value in enumerate(residuals):
        while value % 2 == 0:
            value //= 2
        while value % 3 == 0:
            value //= 3
        residuals[idx] = value

    for ell in primes_up_to(prime_limit):
        if ell < 5:
            continue
        weight = blocker_weight(ell) / ell
        first = ((n_start + ell - 1) // ell) * ell
        for n in range(first, n_stop, ell):
            idx = n - n_start
            if residuals[idx] % ell == 0:
                loads[idx] += weight
                if ell in small_index:
                    small_masks[idx] |= 1 << small_index[ell]
                while residuals[idx] % ell == 0:
                    residuals[idx] //= ell

    for idx, residual in enumerate(residuals):
        if residual >= 5:
            loads[idx] += blocker_weight(int(residual)) / residual
            if residual in small_index:
                small_masks[idx] |= 1 << small_index[int(residual)]

    max_load = max(loads)
    first_idx = next(i for i, value in enumerate(loads) if value == max_load)
    n = n_start + first_idx
    factors = sorted(int(p) for p in factorint(n))
    factors_ge5 = [p for p in factors if p >= 5]
    checked_load = weighted_load_from_factors(factors_ge5)
    if abs(checked_load - max_load) > 1e-12:
        raise AssertionError(f"weighted load mismatch for n={n}: sieve={max_load}, factorint={checked_load}")
    mask = int(small_masks[first_idx])
    small_signature = [ell for ell, i in small_index.items() if mask & (1 << i)]
    over_one_count = sum(1 for value in loads if value >= 1.0)
    return {
        "mode": "weighted_exact_segmented_sieve",
        "weight_rule": "k(5)=3; k(7,11,13,17)=2; k(ell>=19)=1",
        "y_start": y_start,
        "y_end_exclusive": y_end_exclusive,
        "n_start": n_start,
        "n_stop_exclusive": n_stop,
        "rows_scanned": length,
        "prime_limit": prime_limit,
        "max_weighted_load": {
            "y": n - 1,
            "n": n,
            "weighted_load": max_load,
            "n_factors": factors,
            "n_factors_ge5": factors_ge5,
            "small_primorial_signature": small_signature,
        },
        "rows_with_weighted_load_ge_1": over_one_count,
    }
