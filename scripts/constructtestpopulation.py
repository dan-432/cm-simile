#!/usr/bin/env python3
"""Build a curated clan population for cm-simile benchmarks.

Extracts a small, hand-scaled population of covariance models from
Rfam.cm -- some known Rfam clans (positives) plus distractor families
(negatives) -- into one multi-model CM file, and writes a companion TSV
annotating every family's clan membership. This is the input the
benchmark runner scripts (run_bitscoresumdifference_benchmark.py,
run_minmaxcompare_benchmark.py) take; they don't know anything about
clans themselves. Join their pairwise output back against this TSV
(matching model_a/model_b against accession/name/fetch_key) to label
pairs and run whatever statistics you want.

Two distractor sources are supported, since which one you use changes
what a clean separation later demonstrates:

  singleton  -- families with no clan annotation at all. An easy
                contrast; separating these mainly confirms a statistic
                isn't degenerate, not that it's discriminating anything
                subtle.
  cross_clan -- families that ARE clan-annotated, just in clans other
                than the ones chosen for this population. A harder,
                more meaningful contrast -- these are curation-wise
                more similar to your positives than a random singleton.

By default the script draws from BOTH, so downstream you get the easy
sanity check and the harder test from the same population.

Note: cross_clan distractors are sampled independently per family, not
per-clan, so it's possible (rare at real Rfam scale, easy to hit in a
small toy population) for both members of some unchosen clan to end up
selected. The TSV's clan column reflects true clan membership from
Rfam.clanin regardless of the "role" the family was drawn for.

Prerequisites
-------------
- Infernal's cmfetch and cmstat on PATH (or pass --cmfetch-bin/--cmstat-bin).
- data/Rfam.cm and data/Rfam.clanin (see data/README.md).

Usage
-----
    python scripts/build_clan_population.py \\
        --rfam-cm data/Rfam.cm --clanin data/Rfam.clanin \\
        --n-clans 8 --n-singleton-distractors 20 --n-crossclan-distractors 20 \\
        --out results/clan_population

Or supply an exact, hand-picked family list instead of automatic random
selection (recommended once you've looked at an auto-selected set and
want to swap in specific clans/distractors relevant to your own work):

    python scripts/build_clan_population.py --families-file my_population.txt ...

(my_population.txt: one Rfam accession or name per line, comments with #)

Output
------
- <out>/population.cm        -- the extracted CM population
- <out>/population_clans.tsv -- one row per family: accession, name,
  clan (empty if none), role, fetch_key
"""
import argparse
import csv
import os
import random
import subprocess
import sys


# --------------------------------------------------------------------------
# Parsing Rfam.clanin and Rfam.cm (cmstat is the single source of truth for
# name<->accession mapping, since Rfam.clanin lists accessions but
# cmfetch's indexed key may be either depending on the Infernal build --
# see the self-test in main()).
# --------------------------------------------------------------------------

def parse_clanin(path):
    """Returns (clan_of: {accession: clan_id}, members_of: {clan_id: [accessions]})."""
    clan_of = {}
    members_of = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            clan_id, members = fields[0], fields[1:]
            members_of[clan_id] = members
            for acc in members:
                clan_of[acc] = clan_id
    return clan_of, members_of


def list_families(cmstat_bin, rfam_cm):
    """Parse cmstat's summary table into a list of (name, accession).
    Column positions are located from the header row by label where
    possible, falling back to the standard layout (name col 1, accession
    col 2, 0-indexed) if header detection fails -- verify against
    `cmstat -h` and a real Rfam.cm if this ever mis-parses.
    """
    proc = subprocess.run([cmstat_bin, rfam_cm], capture_output=True, text=True, check=True)

    name_col, acc_col = None, None
    families = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            header = line.lstrip("#").lower().split()
            for i, tok in enumerate(header):
                if "name" in tok and name_col is None:
                    name_col = i
                if "acc" in tok and acc_col is None:
                    acc_col = i
            continue
        fields = line.split()
        nc = name_col if name_col is not None else 1
        ac = acc_col if acc_col is not None else 2
        if len(fields) <= max(nc, ac):
            continue
        families.append((fields[nc], fields[ac]))
    if not families:
        raise RuntimeError(
            f"cmstat produced no parseable family rows for {rfam_cm} -- "
            "check `cmstat -h` output format against list_families()'s assumptions"
        )
    return families


# --------------------------------------------------------------------------
# CM extraction
# --------------------------------------------------------------------------

def ensure_ssi_index(cmfetch_bin, rfam_cm):
    if os.path.exists(rfam_cm + ".ssi"):
        return
    print(f"[build_clan_population] building SSI index for {rfam_cm} ...", file=sys.stderr)
    subprocess.run([cmfetch_bin, "--index", rfam_cm], check=True, capture_output=True)


def extract_population_cm(cmfetch_bin, rfam_cm, keys, out_path):
    """Extract keys from rfam_cm into one combined CM file at out_path.
    Tries `cmfetch -f <keyfile>` (the Easel-suite batch-fetch convention
    shared with hmmfetch); falls back to one cmfetch call per key if -f
    isn't supported by the installed Infernal version.
    """
    keyfile = out_path + ".keys"
    with open(keyfile, "w") as f:
        f.write("\n".join(keys) + "\n")
    try:
        with open(out_path, "w") as out:
            subprocess.run(
                [cmfetch_bin, "-f", rfam_cm, keyfile],
                stdout=out, check=True, stderr=subprocess.PIPE, text=True,
            )
    except subprocess.CalledProcessError:
        with open(out_path, "w") as out:
            for key in keys:
                r = subprocess.run([cmfetch_bin, rfam_cm, key], capture_output=True, text=True)
                if r.returncode != 0:
                    print(f"[build_clan_population] warning: cmfetch failed for {key!r}: {r.stderr.strip()}", file=sys.stderr)
                    continue
                out.write(r.stdout)
    finally:
        os.remove(keyfile)


# --------------------------------------------------------------------------
# Population construction
# --------------------------------------------------------------------------

def select_population(members_of, all_accessions, clan_of, n_clans,
                       n_singleton_distractors, n_crossclan_distractors, rng):
    """Returns (population: [(accession, role)], chosen_clans: set[clan_id])."""
    eligible_clans = [c for c, members in members_of.items() if len(members) >= 2]
    if n_clans > len(eligible_clans):
        raise ValueError(f"requested {n_clans} clans but only {len(eligible_clans)} clans have >=2 members")
    chosen_clan_ids = set(rng.sample(eligible_clans, n_clans))

    population = []
    for c in chosen_clan_ids:
        population.extend((acc, "chosen_clan") for acc in members_of[c])

    singleton_pool = [a for a in all_accessions if a not in clan_of]
    crossclan_pool = [a for a in all_accessions if a in clan_of and clan_of[a] not in chosen_clan_ids]

    if n_singleton_distractors > len(singleton_pool):
        print(f"[build_clan_population] warning: only {len(singleton_pool)} singleton families available, "
              f"using all of them instead of {n_singleton_distractors}", file=sys.stderr)
        n_singleton_distractors = len(singleton_pool)
    if n_crossclan_distractors > len(crossclan_pool):
        print(f"[build_clan_population] warning: only {len(crossclan_pool)} cross-clan families available, "
              f"using all of them instead of {n_crossclan_distractors}", file=sys.stderr)
        n_crossclan_distractors = len(crossclan_pool)

    population.extend((acc, "singleton_distractor") for acc in rng.sample(singleton_pool, n_singleton_distractors))
    population.extend((acc, "crossclan_distractor") for acc in rng.sample(crossclan_pool, n_crossclan_distractors))
    return population, chosen_clan_ids


def load_families_file(path):
    accs = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                accs.append(line)
    return [(acc, "manual") for acc in accs]


# --------------------------------------------------------------------------
# Writing the clan annotation TSV
# --------------------------------------------------------------------------

def write_clans_tsv(path, population, clan_of, name_by_acc, key_for):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "name", "clan", "role", "fetch_key"])
        for acc, role in population:
            writer.writerow([acc, name_by_acc.get(acc, ""), clan_of.get(acc, ""), role, key_for(acc)])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rfam-cm", default="data/Rfam.cm")
    p.add_argument("--clanin", default="data/Rfam.clanin")
    p.add_argument("--out", default="results/clan_population", help="output directory")

    p.add_argument("--families-file", default=None,
                    help="exact population as a file of accessions/names, one per line "
                         "(overrides --n-clans/--n-*-distractors; clan labels are still "
                         "derived automatically from --clanin)")
    p.add_argument("--n-clans", type=int, default=8)
    p.add_argument("--n-singleton-distractors", type=int, default=20)
    p.add_argument("--n-crossclan-distractors", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("--cmfetch-bin", default="cmfetch")
    p.add_argument("--cmstat-bin", default="cmstat")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print("[build_clan_population] parsing clan definitions ...", file=sys.stderr)
    clan_of, members_of = parse_clanin(args.clanin)

    print("[build_clan_population] listing Rfam families via cmstat ...", file=sys.stderr)
    families = list_families(args.cmstat_bin, args.rfam_cm)
    name_by_acc = {acc: name for name, acc in families}
    all_accessions = [acc for _, acc in families if acc]

    ensure_ssi_index(args.cmfetch_bin, args.rfam_cm)

    rng = random.Random(args.seed)

    if args.families_file:
        population = load_families_file(args.families_file)
        print(f"[build_clan_population] using {len(population)} hand-picked families from {args.families_file}", file=sys.stderr)
    else:
        population, chosen_clans = select_population(
            members_of, all_accessions, clan_of,
            args.n_clans, args.n_singleton_distractors, args.n_crossclan_distractors, rng,
        )
        print(f"[build_clan_population] selected {len(chosen_clans)} clans, "
              f"{len(population)} families total", file=sys.stderr)

    accessions = [acc for acc, _ in population]

    # Self-test: does cmfetch key on accession (as used in Rfam.clanin) or
    # does it need the family's short name instead?
    use_accession_keys = True
    test_key = accessions[0]
    r = subprocess.run([args.cmfetch_bin, args.rfam_cm, test_key], capture_output=True, text=True)
    if r.returncode != 0:
        use_accession_keys = False
        print("[build_clan_population] cmfetch didn't resolve by accession; falling back to family names", file=sys.stderr)

    def key_for(acc):
        return acc if use_accession_keys else name_by_acc.get(acc, acc)

    population_cm = os.path.join(args.out, "population.cm")
    clans_tsv = os.path.join(args.out, "population_clans.tsv")

    print(f"[build_clan_population] extracting {len(accessions)} families ...", file=sys.stderr)
    extract_population_cm(args.cmfetch_bin, args.rfam_cm, [key_for(a) for a in accessions], population_cm)

    write_clans_tsv(clans_tsv, population, clan_of, name_by_acc, key_for)

    print(f"[build_clan_population] wrote {population_cm}", file=sys.stderr)
    print(f"[build_clan_population] wrote {clans_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()