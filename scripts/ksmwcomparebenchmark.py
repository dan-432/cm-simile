#!/usr/bin/env python3
"""Run ksmwcompare over a pre-built CM population.

Thin wrapper: takes a multi-model CM file (as produced by
build_clan_population.py) and calls the ksmwcompare binary over it,
writing the all-vs-all comparison results directly to --out. No clan
labeling and no statistics here -- join the output against the
population's population_clans.tsv yourself downstream.

Prerequisites
-------------
- A built ksmwcompare binary:
    go build -o bin/ksmwcompare src/ksmwcompare.go
- Infernal's cmfetch and cmstat on PATH (or pass --cmfetch-bin/--cmstat-bin) --
  ksmwcompare's own Expander needs these to split --cm-file back into
  individual models internally, so they must resolve the same way they did
  when the population was built.

Usage
-----
    python scripts/ksmwcomparebenchmark.py \
        --cm-file results/clan_population/population.cm \
        --binary bin/ksmwcompare \
        --npool 100 --n-null-pools 1 --threads 8 \
        --out results/ksmwcompare/raw_results.tsv
"""
import argparse
import os
import shutil
import subprocess
import sys


def run(binary, cm_file, out_tsv, npool, n_null_pools, seed, threads,
        cmemit_bin, cmsearch_bin, cmstat_bin, cmfetch_bin,
        cmsearch_cpu, min_bitscore, disable_filters):
    cmd = [
        binary,
        "-cm-file", cm_file,
        "-npool", str(npool),
        "-n-null-pools", str(n_null_pools),
        "-seed", str(seed),
        "-threads", str(threads),
        "-out", out_tsv,
    ]
    if cmemit_bin:
        cmd += ["-cmemit-bin", cmemit_bin]
    if cmsearch_bin:
        cmd += ["-cmsearch-bin", cmsearch_bin]
    cmd += ["-cmstat-bin", cmstat_bin, "-cmfetch-bin", cmfetch_bin]
    cmd += ["-cmsearch-cpu", str(cmsearch_cpu), "-min-bitscore", str(min_bitscore)]
    if disable_filters:
        cmd += ["-disable-filters"]
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cm-file", required=True, help="multi-model CM population, e.g. from build_clan_population.py")
    p.add_argument("--binary", default="bin/ksmwcompare")
    p.add_argument("--out", required=True, help="output TSV path (written directly by the binary)")

    p.add_argument("--npool", type=int, default=100)
    p.add_argument("--n-null-pools", type=int, default=1)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--cmsearch-cpu", type=int, default=1)
    p.add_argument("--min-bitscore", type=float, default=-1000.0)
    p.add_argument("--disable-filters", action=argparse.BooleanOptionalAction, default=True,
                   help="pass -disable-filters to cmsearch so weak/shuffled sequences are scored")

    p.add_argument("--cmfetch-bin", default="cmfetch")
    p.add_argument("--cmstat-bin", default="cmstat")
    p.add_argument("--cmemit-bin", default=None)
    p.add_argument("--cmsearch-bin", default=None)
    args = p.parse_args()

    if not shutil.which(args.binary) and not os.path.isfile(args.binary):
        sys.exit(f"error: ksmwcompare binary not found at {args.binary!r} -- build it first (see README)")
    if not os.path.isfile(args.cm_file):
        sys.exit(f"error: --cm-file not found: {args.cm_file!r}")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"[ksmwcomparebenchmark] running over {args.cm_file} ...", file=sys.stderr)
    run(args.binary, args.cm_file, args.out, args.npool, args.n_null_pools,
        args.seed, args.threads, args.cmemit_bin, args.cmsearch_bin,
        args.cmstat_bin, args.cmfetch_bin, args.cmsearch_cpu, args.min_bitscore,
        args.disable_filters)
    print(f"[ksmwcomparebenchmark] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()