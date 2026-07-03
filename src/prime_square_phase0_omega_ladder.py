#!/usr/bin/env python3
"""Segmented exact omega_ge5 ladder for Phase-0 Route-A obstruction checks.

This scans n = y + 1 for y in [p^2, q^2), where q=nextprime(p) unless
provided. It computes max omega_{>=5}(n) exactly by segmented division with all
primes up to q. The runner writes resumable progress JSON after every chunk.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import threading
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from numba import njit
from sympy import factorint, nextprime

from prime_square_segment_a_phase0_utils import primes_up_to


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


class ProgressWatchdog:
    def __init__(self, *, timeout_s: float, progress_path: Path | None, enabled: bool = True) -> None:
        self.timeout_s = timeout_s
        self.progress_path = progress_path
        self.enabled = enabled and timeout_s > 0
        self._lock = threading.Lock()
        self._rows_done = 0
        self._last_progress = time.monotonic()
        self._stop = False
        self._thread: threading.Thread | None = None

    def start(self) -> None:
        if not self.enabled:
            return
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def update(self, rows_done: int) -> None:
        with self._lock:
            if rows_done != self._rows_done:
                self._rows_done = rows_done
                self._last_progress = time.monotonic()

    def stop(self) -> None:
        with self._lock:
            self._stop = True
        if self._thread is not None:
            self._thread.join(timeout=1.0)

    def _run(self) -> None:
        while True:
            time.sleep(min(1.0, max(0.1, self.timeout_s / 10.0)))
            with self._lock:
                if self._stop:
                    return
                idle_s = time.monotonic() - self._last_progress
                rows_done = self._rows_done
            if idle_s > self.timeout_s:
                payload = {
                    "runner": "prime_square_phase0_omega_ladder.v1",
                    "status": "stalled",
                    "rows_done": rows_done,
                    "idle_s": round(idle_s, 3),
                    "timeout_s": self.timeout_s,
                    "last_update_iso": now_iso(),
                    "message": "watchdog exit: rows_done did not increment within timeout",
                }
                if self.progress_path is not None:
                    write_json(self.progress_path, payload)
                print(payload["message"], flush=True)
                os._exit(2)


@njit(cache=True)
def strip_2_3(residuals: np.ndarray) -> None:
    for i in range(residuals.size):
        v = residuals[i]
        while v % 2 == 0:
            v //= 2
        while v % 3 == 0:
            v //= 3
        residuals[i] = v


@njit(cache=True)
def mark_prime_factors(n_start: int, n_stop: int, primes: np.ndarray, residuals: np.ndarray, counts: np.ndarray) -> None:
    for pi in range(primes.size):
        ell = int(primes[pi])
        if ell < 5:
            continue
        first = ((n_start + ell - 1) // ell) * ell
        for n in range(first, n_stop, ell):
            idx = n - n_start
            v = residuals[idx]
            if v % ell == 0:
                counts[idx] += 1
                while v % ell == 0:
                    v //= ell
                residuals[idx] = v


@njit(cache=True)
def add_large_residuals(residuals: np.ndarray, counts: np.ndarray) -> None:
    for i in range(residuals.size):
        if residuals[i] >= 5:
            counts[i] += 1


def chunk_scan(n_start: int, n_stop: int, primes: np.ndarray) -> tuple[int, int, dict[int, int]]:
    residuals = np.arange(n_start, n_stop, dtype=np.uint64)
    counts = np.zeros(n_stop - n_start, dtype=np.uint8)
    strip_2_3(residuals)
    mark_prime_factors(n_start, n_stop, primes, residuals, counts)
    add_large_residuals(residuals, counts)
    max_count = int(counts.max())
    first_idx = int(np.argmax(counts))
    hist_counts = np.bincount(counts.astype(np.int16), minlength=max_count + 1)
    hist = {i: int(v) for i, v in enumerate(hist_counts) if v}
    return max_count, n_start + first_idx, hist


def merge_hist(dst: Counter[int], src: dict[int, int]) -> None:
    for k, v in src.items():
        dst[int(k)] += int(v)


def load_progress(path: Path | None, p: int, q: int, n_start: int, n_stop: int) -> dict[str, Any] | None:
    if path is None or not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    if (
        payload.get("runner") != "prime_square_phase0_omega_ladder.v1"
        or payload.get("p") != p
        or payload.get("q") != q
        or payload.get("n_start") != n_start
        or payload.get("n_stop_exclusive") != n_stop
    ):
        return None
    return payload


def progress_payload(
    *,
    status: str,
    p: int,
    q: int,
    target: int,
    n_start: int,
    n_stop: int,
    next_n: int,
    chunk_size: int,
    max_record: dict[str, Any] | None,
    hist: Counter[int],
    started_at: float,
) -> dict[str, Any]:
    rows_done = next_n - n_start
    rows_total = n_stop - n_start
    return {
        "runner": "prime_square_phase0_omega_ladder.v1",
        "status": status,
        "target": target,
        "p": p,
        "q": q,
        "gap": q - p,
        "n_start": n_start,
        "n_stop_exclusive": n_stop,
        "rows_done": rows_done,
        "rows_total": rows_total,
        "percent": round(100.0 * rows_done / rows_total, 6),
        "next_n": next_n,
        "chunk_size": chunk_size,
        "last_update_iso": now_iso(),
        "wall_elapsed_s": round(time.monotonic() - started_at, 3),
        "max_omega_ge5": max_record,
        "omega_ge5_hist": {str(k): v for k, v in sorted(hist.items())},
    }


def write_md(path: Path, payload: dict[str, Any]) -> None:
    rec = payload["max_omega_ge5"]
    lines = [
        "# Phase-0 Omega Ladder Result",
        "",
        f"- generated_at: `{payload['generated_at']}`",
        f"- target: `{payload['target']}`",
        f"- p: `{payload['p']}`",
        f"- q: `{payload['q']}`",
        f"- gap: `{payload['gap']}`",
        f"- rows_scanned: `{payload['rows_total']}`",
        f"- max_omega_ge5: `{rec['omega_ge5']}`",
        f"- witness_y: `{rec['y']}`",
        f"- witness_n: `{rec['n']}`",
        f"- factors_ge5: `{rec['n_factors_ge5']}`",
        f"- logp_squared: `{payload['logp_squared']:.6f}`",
        f"- 2^omega/logp_squared: `{payload['pow2_over_logp_squared']:.6f}`",
        "",
        "## Histogram",
        "",
        "```json",
        json.dumps(payload["omega_ge5_hist"], indent=2, sort_keys=True),
        "```",
        "",
    ]
    atomic_write_text(path, "\n".join(lines))


def run(args: argparse.Namespace) -> dict[str, Any]:
    target = int(args.target)
    p = int(args.p) if args.p else int(nextprime(target))
    q = int(args.q) if args.q else int(nextprime(p))
    n_start = p * p + 1
    n_stop = q * q + 1
    rows_total = n_stop - n_start
    started_at = time.monotonic()
    progress_path = Path(args.progress_json) if args.progress_json else None

    progress = load_progress(progress_path, p, q, n_start, n_stop)
    if progress and progress.get("status") == "complete":
        return progress

    if progress:
        next_n = int(progress["next_n"])
        hist = Counter({int(k): int(v) for k, v in progress.get("omega_ge5_hist", {}).items()})
        max_record = progress.get("max_omega_ge5")
    else:
        next_n = n_start
        hist: Counter[int] = Counter()
        max_record: dict[str, Any] | None = None
    watchdog = ProgressWatchdog(
        timeout_s=args.no_progress_timeout_s,
        progress_path=progress_path,
        enabled=not args.disable_watchdog,
    )
    watchdog.update(next_n - n_start)

    print(f"Generating primes up to q={q}", flush=True)
    primes = np.array([ell for ell in primes_up_to(q) if ell >= 5], dtype=np.int64)
    print(f"Prime count >=5: {primes.size}", flush=True)
    watchdog.start()

    chunks_done = 0
    while next_n < n_stop:
        chunk_stop = min(n_stop, next_n + args.chunk_size)
        max_count, witness_n, chunk_hist = chunk_scan(next_n, chunk_stop, primes)
        merge_hist(hist, chunk_hist)
        if max_record is None or max_count > int(max_record["omega_ge5"]):
            factors = sorted(int(x) for x in factorint(witness_n))
            factors_ge5 = [x for x in factors if x >= 5]
            if len(factors_ge5) != max_count:
                raise AssertionError(f"omega mismatch at n={witness_n}: sieve={max_count}, factorint={factors_ge5}")
            max_record = {
                "y": witness_n - 1,
                "n": witness_n,
                "omega_ge5": max_count,
                "n_factors": factors,
                "n_factors_ge5": factors_ge5,
            }
        next_n = chunk_stop
        watchdog.update(next_n - n_start)
        chunks_done += 1
        payload = progress_payload(
            status="running",
            p=p,
            q=q,
            target=target,
            n_start=n_start,
            n_stop=n_stop,
            next_n=next_n,
            chunk_size=args.chunk_size,
            max_record=max_record,
            hist=hist,
            started_at=started_at,
        )
        if progress_path:
            write_json(progress_path, payload)
        print(
            f"OMEGA target={target} p={p} {payload['percent']}% rows={payload['rows_done']}/{rows_total} "
            f"max={max_record['omega_ge5'] if max_record else None}",
            flush=True,
        )
        if args.simulate_stall_after_chunks and chunks_done >= args.simulate_stall_after_chunks:
            print("Simulating stall for watchdog verification", flush=True)
            time.sleep(args.no_progress_timeout_s + 5.0)

    assert max_record is not None
    watchdog.stop()
    logp2 = math.log(p) ** 2
    omega = int(max_record["omega_ge5"])
    payload = progress_payload(
        status="complete",
        p=p,
        q=q,
        target=target,
        n_start=n_start,
        n_stop=n_stop,
        next_n=n_stop,
        chunk_size=args.chunk_size,
        max_record=max_record,
        hist=hist,
        started_at=started_at,
    )
    payload.update(
        {
            "generated_at": now_iso(),
            "logp_squared": logp2,
            "pow2_omega": 2**omega,
            "pow2_over_logp_squared": (2**omega) / logp2,
        }
    )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, required=True)
    parser.add_argument("--p", type=int)
    parser.add_argument("--q", type=int)
    parser.add_argument("--chunk-size", type=int, default=2_000_000)
    parser.add_argument("--progress-json")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-md", required=True)
    parser.add_argument("--no-progress-timeout-s", type=float, default=60.0)
    parser.add_argument("--disable-watchdog", action="store_true")
    parser.add_argument("--simulate-stall-after-chunks", type=int, default=0, help=argparse.SUPPRESS)
    args = parser.parse_args()

    payload = run(args)
    write_json(Path(args.output_json), payload)
    write_md(Path(args.output_md), payload)
    if args.progress_json:
        write_json(Path(args.progress_json), payload)
    print(f"WROTE {args.output_json}", flush=True)
    print(f"WROTE {args.output_md}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
