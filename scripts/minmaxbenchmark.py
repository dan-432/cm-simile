
"""Run minmaxcompare over a pre-built CM population.

Thin wrapper: takes a multi-model CM file (as produced by
build_clan_population.py) and calls the minmaxcompare binary over it,
writing its all-vs-all results straight to --out. No clan labeling and
no statistics here -- join the output against the population's
population_clans.tsv yourself downstream.

Prerequisites
-------------
- A built minmaxcompare binary:
    go build -o bin/minmaxcompare src/minmaxcompare.go
- Infernal's cmfetch and cmstat on PATH (or pass --cmfetch-bin/--cmstat-bin) --
  minmaxcompare's own Expander needs these to split --cm-file back into
  individual models internally, so they must resolve the same way they
  did when the population was built.

Usage
-----
    python minmaxcomparebenchmark.py \\
        --cm-file results/clan_population/population.cm \\
        --binary bin/minmaxcompare \\
        --generations 100 --pop-size 40 --threads 8 \\
        --out results/minmaxcompare/raw_results.tsv
"""
import argparse
import os
import shutil
import subprocess
import sys


def run(binary, cm_file, out_tsv, generations, pop_size, seed, threads,
        cmemit_bin, cmsearch_bin, cmstat_bin, cmfetch_bin):
    cmd = [
        binary,
        "-cm-file", cm_file,
        "-generations", str(generations),
        "-pop-size", str(pop_size),
        "-seed", str(seed),
        "-cmsearch-cpu", str(threads),
        "-out", out_tsv,
    ]
    if cmemit_bin:
        cmd += ["-cmemit-bin", cmemit_bin]
    if cmsearch_bin:
        cmd += ["-cmsearch-bin", cmsearch_bin]
    cmd += ["-cmstat-bin", cmstat_bin, "-cmfetch-bin", cmfetch_bin]
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cm-file", required=True, help="multi-model CM population, e.g. from build_clan_population.py")
    p.add_argument("--binary", default="bin/minmaxcompare")
    p.add_argument("--out", required=True, help="output TSV path (written directly by the binary)")

    p.add_argument("--generations", type=int, default=100)
    p.add_argument("--pop-size", type=int, default=40)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--threads", type=int, default=4)

    p.add_argument("--cmfetch-bin", default="cmfetch")
    p.add_argument("--cmstat-bin", default="cmstat")
    p.add_argument("--cmemit-bin", default=None)
    p.add_argument("--cmsearch-bin", default=None)
    args = p.parse_args()

    if not shutil.which(args.binary) and not os.path.isfile(args.binary):
        sys.exit(f"error: minmaxcompare binary not found at {args.binary!r} -- build it first (see README)")
    if not os.path.isfile(args.cm_file):
        sys.exit(f"error: --cm-file not found: {args.cm_file!r}")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"[run_minmaxcompare_benchmark] running over {args.cm_file} ...", file=sys.stderr)
    run(args.binary, args.cm_file, args.out, args.generations, args.pop_size, args.seed,
        args.threads, args.cmemit_bin, args.cmsearch_bin, args.cmstat_bin, args.cmfetch_bin)
    print(f"[run_minmaxcompare_benchmark] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()