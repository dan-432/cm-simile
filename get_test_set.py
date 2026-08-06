#!/usr/bin/env python3
"""Curated clan benchmark for cm-simile.

Positive-control validation: build a small, hand-scaled population of
covariance models -- some known Rfam clans (positives: within-clan
pairs) plus distractor families (negatives: everything else) -- run one
bitscoresumdifference pass over the whole population, label every
resulting pair after the fact from clan membership, and test:

    H0: within-clan pairs and non-within-clan pairs are drawn from the
        same distribution of the similarity statistic (no separation).
    H1: within-clan pairs show systematically stronger cross-recognition.

This trades the broader survey design (random sampling across all of
Rfam) for something cheaper and more directly interpretable: a single
small all-vs-all run, not a batched sweep. It's the right first step
before spending compute on a broader validation -- if the statistic
can't even separate a curated, favourable case, a large random survey
isn't going to be more informative.

Two distractor sources are supported, since which one you use changes
what a clean separation actually demonstrates:

  singleton  -- families with no clan annotation at all. An easy
                contrast; separating these mainly confirms the
                statistic isn't degenerate, not that it's discriminating
                anything subtle.
  cross_clan -- families that ARE clan-annotated, just in clans other
                than the ones chosen for this population. A harder,
                more meaningful contrast (see the discussion in the
                conversation this script came out of) -- these are
                curation-wise more similar to your positives than a
                random singleton, so separating them from within-clan
                pairs is closer to demonstrating real discriminative
                power rather than just "related vs. utterly unrelated".

By default the script draws from BOTH, so you get the easy sanity check
and the harder test in the same run, labeled separately.

Note: cross_clan distractors are sampled independently per family, not
per-clan, so it's possible (rare at real Rfam scale, easy to hit in a
small toy population) for both members of some unchosen clan to end up
selected. Their pair is then correctly labeled within_clan -- the label
reflects true clan membership from Rfam.clanin, not whether a clan was
one you deliberately chose as a "positive" clan.

Prerequisites
-------------
- A built bitscoresumdifference binary:
    go build -o bin/bitscoresumdifference src/bitscoresumdifference.go
- Infernal's cmfetch and cmstat on PATH (or pass --cmfetch-bin/--cmstat-bin).
- data/Rfam.cm and data/Rfam.clanin (see data/README.md).

Usage
-----
    python scripts/run_clan_benchmark.py \\
        --rfam-cm data/Rfam.cm --clanin data/Rfam.clanin \\
        --binary bin/bitscoresumdifference \\
        --n-clans 8 --n-singleton-distractors 20 --n-crossclan-distractors 20 \\
        --npool 100 --nshuffles 200 --threads 8 \\
        --out results/clan_benchmark

Or supply an exact, hand-picked family list instead of automatic random
selection (recommended once you've looked at the auto-selected set and
want to swap in specific clans/distractors relevant to your own work):

    python scripts/run_clan_benchmark.py --families-file my_population.txt ...

(my_population.txt: one Rfam accession or name per line, comments with #)

Output
------
- <out>/population.cm       -- the extracted CM population, for reuse/inspection
- <out>/all_pairs.tsv        -- every pair from the run, labeled with group
- Summary statistics (pair counts, Mann-Whitney U, AUC, p-value) printed
  to stderr, comparing within_clan vs. everything else, and separately
  vs. each distractor source.
"""
import argparse
import csv
import math
import os
import random
import shutil
import subprocess
import sys
from collections import defaultdict


# --------------------------------------------------------------------------
# Parsing Rfam.clanin and Rfam.cm (same approach as discussed earlier:
# cmstat is the single source of truth for name<->accession mapping, since
# Rfam.clanin lists accessions but cmfetch's indexed key may be either
# depending on the Infernal build -- see the self-test in main()).
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
    print(f"[run_clan_benchmark] building SSI index for {rfam_cm} ...", file=sys.stderr)
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
                    print(f"[run_clan_benchmark] warning: cmfetch failed for {key!r}: {r.stderr.strip()}", file=sys.stderr)
                    continue
                out.write(r.stdout)
    finally:
        os.remove(keyfile)


# --------------------------------------------------------------------------
# Population construction
# --------------------------------------------------------------------------

def select_population(members_of, all_accessions, clan_of, n_clans,
                       n_singleton_distractors, n_crossclan_distractors, rng):
    """Returns (population: [accessions], chosen_clans: set[clan_id])."""
    eligible_clans = [c for c, members in members_of.items() if len(members) >= 2]
    if n_clans > len(eligible_clans):
        raise ValueError(f"requested {n_clans} clans but only {len(eligible_clans)} clans have >=2 members")
    chosen_clan_ids = set(rng.sample(eligible_clans, n_clans))

    population = []
    for c in chosen_clan_ids:
        population.extend(members_of[c])

    singleton_pool = [a for a in all_accessions if a not in clan_of]
    crossclan_pool = [a for a in all_accessions if a in clan_of and clan_of[a] not in chosen_clan_ids]

    if n_singleton_distractors > len(singleton_pool):
        print(f"[run_clan_benchmark] warning: only {len(singleton_pool)} singleton families available, "
              f"using all of them instead of {n_singleton_distractors}", file=sys.stderr)
        n_singleton_distractors = len(singleton_pool)
    if n_crossclan_distractors > len(crossclan_pool):
        print(f"[run_clan_benchmark] warning: only {len(crossclan_pool)} cross-clan families available, "
              f"using all of them instead of {n_crossclan_distractors}", file=sys.stderr)
        n_crossclan_distractors = len(crossclan_pool)

    population.extend(rng.sample(singleton_pool, n_singleton_distractors))
    population.extend(rng.sample(crossclan_pool, n_crossclan_distractors))
    return population, chosen_clan_ids


def load_families_file(path):
    accs = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                accs.append(line)
    return accs


# --------------------------------------------------------------------------
# Running bitscoresumdifference and labeling results
# --------------------------------------------------------------------------

def run_population(binary, population_cm, out_tsv, npool, nshuffles, seed, threads,
                    cmemit_bin, cmsearch_bin, cmstat_bin, cmfetch_bin):
    cmd = [
        binary,
        "-cm-file", population_cm,
        "-npool", str(npool),
        "-nshuffles", str(nshuffles),
        "-seed", str(seed),
        "-threads", str(threads),
        "-out", out_tsv,
    ]
    if cmemit_bin:
        cmd += ["-cmemit-bin", cmemit_bin]
    if cmsearch_bin:
        cmd += ["-cmsearch-bin", cmsearch_bin]
    # population.cm has multiple models, so bitscoresumdifference's own
    # Expander re-runs cmstat/cmfetch internally to split it back into
    # individual models -- it needs to resolve the SAME binaries this
    # script used, not whatever "cmstat"/"cmfetch" happen to be on PATH.
    cmd += ["-cmstat-bin", cmstat_bin, "-cmfetch-bin", cmfetch_bin]
    subprocess.run(cmd, check=True)


def label_pairs(out_tsv, clan_of, labeled_out):
    """Reads bitscoresumdifference's output TSV, adds group/clan columns,
    and returns the labeled rows (also written to labeled_out).

    clan_of should be keyed by the same family identifier that appears in
    the TSV output (for example the names from Rfam.clanin / the extracted
    models that cmstat/cmfetch resolve during population construction).

    group is one of: within_clan, cross_clan, background
      within_clan -- both families in the same chosen clan
      cross_clan  -- both families clan-annotated, different clans
      background  -- at least one family is a singleton (non-clan) distractor
    """
    rows = []
    with open(out_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            a, b = row["model_a"], row["model_b"]
            clan_a, clan_b = clan_of.get(a, ""), clan_of.get(b, "")
            if clan_a and clan_a == clan_b:
                group = "within_clan"
            elif clan_a and clan_b:
                group = "cross_clan"
            else:
                group = "background"
            row["group"] = group
            row["clan_a"] = clan_a
            row["clan_b"] = clan_b
            rows.append(row)

    fieldnames = ["group", "model_a", "model_b", "clan_a", "clan_b",
                  "mean_score_a_on_a", "mean_score_a_on_b",
                  "mean_score_b_on_b", "mean_score_b_on_a",
                  "stat_a_into_b", "stat_b_into_a",
                  "pvalue_a_into_b", "pvalue_b_into_a", "n_shuffles"]
    with open(labeled_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


# --------------------------------------------------------------------------
# Statistics: Mann-Whitney U / AUC, no external dependencies
# --------------------------------------------------------------------------

def mannwhitney_u_auc(group_a_values, group_b_values):
    """Two-sided Mann-Whitney U test with normal approximation (with tie
    correction), plus the AUC equivalent (AUC = U / (n1*n2), the
    probability a random group_a value exceeds a random group_b value,
    with ties counted as a half-win -- the standard Mann-Whitney/AUC
    correspondence). Implemented directly rather than via scipy so this
    script has no dependencies beyond the Python standard library.

    Returns (auc, u_statistic, z_score, p_value_two_sided).
    """
    n1, n2 = len(group_a_values), len(group_b_values)
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")

    combined = [(v, 0) for v in group_a_values] + [(v, 1) for v in group_b_values]
    combined.sort(key=lambda t: t[0])

    # Assign ranks (1-indexed), averaging ranks across ties.
    ranks = [0.0] * len(combined)
    i = 0
    tie_term = 0.0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0  # average of ranks i+1..j (1-indexed)
        for k in range(i, j):
            ranks[k] = avg_rank
        t = j - i
        tie_term += t ** 3 - t
        i = j

    rank_sum_a = sum(r for r, (_, g) in zip(ranks, combined) if g == 0)
    u1 = rank_sum_a - n1 * (n1 + 1) / 2.0
    auc = u1 / (n1 * n2)

    n = n1 + n2
    mean_u = n1 * n2 / 2.0
    # Tie-corrected variance (standard normal-approximation formula).
    var_u = (n1 * n2 / 12.0) * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0.0
    if var_u <= 0:
        return auc, u1, float("nan"), float("nan")
    z = (u1 - mean_u) / math.sqrt(var_u)
    p = 2.0 * (1.0 - _norm_cdf(abs(z)))
    return auc, u1, z, p


def _norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def summarize_comparison(label, rows_a, rows_b, statistic_key):
    values_a = [float(r[statistic_key]) for r in rows_a]
    values_b = [float(r[statistic_key]) for r in rows_b]
    if not values_a or not values_b:
        print(f"  {label}: skipped (one group is empty)", file=sys.stderr)
        return
    auc, u, z, p = mannwhitney_u_auc(values_a, values_b)
    print(f"  {label}:", file=sys.stderr)
    print(f"    n = {len(values_a)} vs {len(values_b)}", file=sys.stderr)
    print(f"    median {statistic_key}: {_median(values_a):.3f} vs {_median(values_b):.3f}", file=sys.stderr)
    print(f"    AUC = {auc:.4f}   Mann-Whitney U = {u:.1f}   z = {z:.3f}   p = {p:.4g}", file=sys.stderr)


def _median(values):
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    return (s[mid - 1] + s[mid]) / 2.0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rfam-cm", default="data/Rfam.cm")
    p.add_argument("--clanin", default="data/Rfam.clanin")
    p.add_argument("--binary", default="bin/bitscoresumdifference")
    p.add_argument("--out", default="results/clan_benchmark")

    p.add_argument("--families-file", default=None,
                    help="exact population as a file of accessions/names, one per line "
                         "(overrides --n-clans/--n-*-distractors; clan/background labels "
                         "are still derived automatically from --clanin)")
    p.add_argument("--n-clans", type=int, default=8)
    p.add_argument("--n-singleton-distractors", type=int, default=20)
    p.add_argument("--n-crossclan-distractors", type=int, default=20)

    p.add_argument("--npool", type=int, default=100)
    p.add_argument("--nshuffles", type=int, default=200)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--threads", type=int, default=4)

    p.add_argument("--cmfetch-bin", default="cmfetch")
    p.add_argument("--cmstat-bin", default="cmstat")
    p.add_argument("--cmemit-bin", default=None)
    p.add_argument("--cmsearch-bin", default=None)
    args = p.parse_args()

    if not shutil.which(args.binary) and not os.path.isfile(args.binary):
        sys.exit(f"error: bitscoresumdifference binary not found at {args.binary!r} -- build it first (see README)")

    os.makedirs(args.out, exist_ok=True)

    print("[run_clan_benchmark] parsing clan definitions ...", file=sys.stderr)
    clan_of, members_of = parse_clanin(args.clanin)

    print("[run_clan_benchmark] listing Rfam families via cmstat ...", file=sys.stderr)
    families = list_families(args.cmstat_bin, args.rfam_cm)
    name_by_acc = {acc: name for name, acc in families}
    all_accessions = [acc for _, acc in families if acc]

    ensure_ssi_index(args.cmfetch_bin, args.rfam_cm)

    rng = random.Random(args.seed)

    if args.families_file:
        population = load_families_file(args.families_file)
        print(f"[run_clan_benchmark] using {len(population)} hand-picked families from {args.families_file}", file=sys.stderr)
    else:
        population, chosen_clans = select_population(
            members_of, all_accessions, clan_of,
            args.n_clans, args.n_singleton_distractors, args.n_crossclan_distractors, rng,
        )
        print(f"[run_clan_benchmark] selected {len(chosen_clans)} clans, "
              f"{len(population)} families total", file=sys.stderr)

    # Self-test: does cmfetch key on accession (as used in Rfam.clanin) or
    # does it need the family's short name instead?
    use_accession_keys = True
    test_key = population[0]
    r = subprocess.run([args.cmfetch_bin, args.rfam_cm, test_key], capture_output=True, text=True)
    if r.returncode != 0:
        use_accession_keys = False
        print("[run_clan_benchmark] cmfetch didn't resolve by accession; falling back to family names", file=sys.stderr)

    def key_for(acc):
        return acc if use_accession_keys else name_by_acc.get(acc, acc)

    population_cm = os.path.join(args.out, "population.cm")
    raw_tsv = os.path.join(args.out, "raw_results.tsv")
    labeled_tsv = os.path.join(args.out, "all_pairs.tsv")

    print(f"[run_clan_benchmark] extracting {len(population)} families ...", file=sys.stderr)
    extract_population_cm(args.cmfetch_bin, args.rfam_cm, [key_for(a) for a in population], population_cm)

    n_pairs = len(population) * (len(population) - 1) // 2
    print(f"[run_clan_benchmark] running bitscoresumdifference over {len(population)} "
          f"families ({n_pairs} pairs) ...", file=sys.stderr)
    run_population(args.binary, population_cm, raw_tsv, args.npool, args.nshuffles,
                    args.seed, args.threads, args.cmemit_bin, args.cmsearch_bin,
                    args.cmstat_bin, args.cmfetch_bin)

    rows = label_pairs(raw_tsv, clan_of, labeled_tsv)
    print(f"[run_clan_benchmark] labeled results written to {labeled_tsv}", file=sys.stderr)

    by_group = defaultdict(list)
    for row in rows:
        by_group[row["group"]].append(row)

    print("", file=sys.stderr)
    print("Pair counts:", file=sys.stderr)
    for group in ("within_clan", "cross_clan", "background"):
        print(f"  {group:12s} {len(by_group[group])}", file=sys.stderr)

    print("", file=sys.stderr)
    print("H0: within-clan and non-within-clan pairs have the same statistic distribution.", file=sys.stderr)
    print("Comparisons (statistic = stat_a_into_b; also worth checking stat_b_into_a", file=sys.stderr)
    print("separately if the two disagree -- see the asymmetric-containment discussion", file=sys.stderr)
    print("in the top-level README):", file=sys.stderr)

    rest = by_group["cross_clan"] + by_group["background"]
    summarize_comparison("within_clan vs. everything else", by_group["within_clan"], rest, "stat_a_into_b")
    summarize_comparison("within_clan vs. cross_clan (hard negative)", by_group["within_clan"], by_group["cross_clan"], "stat_a_into_b")
    summarize_comparison("within_clan vs. background (easy negative)", by_group["within_clan"], by_group["background"], "stat_a_into_b")


if __name__ == "__main__":
    main()
