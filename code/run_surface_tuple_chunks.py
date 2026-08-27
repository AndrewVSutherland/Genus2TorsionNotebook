#!/usr/bin/env python3
"""Chunked driver for enumerate_surface_tuples.

Runs the C++ K3 tuple enumerator over consecutive ax-ranges, writes each
chunk separately, and merges unique tuple rows.  This is intended for long
boundary/local-filtered searches where monolithic runs are inconvenient.
"""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path


def read_rows(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def merge_chunks(chunk_paths: list[Path], output: Path) -> int:
    rows = set()
    for path in chunk_paths:
        rows.update(read_rows(path))
    ordered = sorted(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(row + "\n" for row in ordered))
    return len(ordered)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--B", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primes", default="11,23")
    parser.add_argument("--boundary", default="11:N,23:N")
    parser.add_argument("--local-depth", type=int, default=4)
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int)
    parser.add_argument("--enumerator", type=Path, default=Path("code/enumerate_surface_tuples"))
    parser.add_argument("--chunk-dir", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    end = args.end or args.B
    chunk_dir = args.chunk_dir or args.output.with_suffix(args.output.suffix + ".chunks")
    chunk_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output.with_suffix(args.output.suffix + ".log")

    chunk_paths: list[Path] = []
    t0 = time.time()
    with log_path.open("a") as log:
        log.write(
            f"run B={args.B} start={args.start} end={end} chunk_size={args.chunk_size} "
            f"primes={args.primes} boundary={args.boundary} local_depth={args.local_depth}\n"
        )
        for ax0 in range(args.start, end + 1, args.chunk_size):
            ax1 = min(end, ax0 + args.chunk_size - 1)
            chunk = chunk_dir / f"ax_{ax0}_{ax1}.txt"
            chunk_paths.append(chunk)
            if args.resume and chunk.exists():
                log.write(f"skip existing {chunk}\n")
                log.flush()
                continue
            tmp_chunk = chunk.with_name(chunk.name + ".tmp")
            if tmp_chunk.exists():
                tmp_chunk.unlink()
            cmd = [
                str(args.enumerator),
                str(args.B),
                str(tmp_chunk),
                args.primes,
                str(ax0),
                str(ax1),
                args.boundary,
                str(args.local_depth),
            ]
            ct0 = time.time()
            log.write(f"start ax={ax0}..{ax1}\n")
            log.flush()
            try:
                result = subprocess.run(cmd, text=True, capture_output=True)
            except KeyboardInterrupt:
                if tmp_chunk.exists():
                    tmp_chunk.unlink()
                log.write(f"interrupted ax={ax0}..{ax1}; removed temporary chunk\n")
                log.flush()
                raise
            elapsed = time.time() - ct0
            log.write(result.stderr)
            if result.stdout:
                log.write("stdout:\n" + result.stdout)
            if result.returncode == 0:
                tmp_chunk.replace(chunk)
            elif tmp_chunk.exists():
                tmp_chunk.unlink()
            rows = len(read_rows(chunk)) if chunk.exists() else 0
            log.write(f"finish ax={ax0}..{ax1} return={result.returncode} elapsed={elapsed:.3f}s rows={rows}\n")
            log.flush()
            print(f"ax={ax0}..{ax1} return={result.returncode} elapsed={elapsed:.1f}s rows={rows}", flush=True)
            if result.returncode != 0:
                raise SystemExit(result.returncode)

    all_chunk_paths = sorted(chunk_dir.glob("ax_*_*.txt"))
    count = merge_chunks(all_chunk_paths, args.output)
    elapsed = time.time() - t0
    with log_path.open("a") as log:
        log.write(f"merged_chunks={len(all_chunk_paths)} unique_rows={count} elapsed_total={elapsed:.3f}s output={args.output}\n")
    print(f"merged {len(all_chunk_paths)} chunks -> {count} unique rows in {args.output}")


if __name__ == "__main__":
    main()
