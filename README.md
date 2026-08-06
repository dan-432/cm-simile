# cm-simile

Fast, statistically defensible estimation of sequence/structure similarity between Infernal covariance models.

Two standalone Go programs, both living in `src/`, sharing the same package infrastructure (`cmemit`/`cmsearch` wrappers, multi-model CM file expansion, dinucleotide shuffling):

- **`bitscoresumdifference`** -- the primary, fast screening tool. Samples real sequences from each model, cross-scores them, and tests the observed directional mean bit-score difference against a dinucleotide-shuffled null distribution for a bootstrap p-value per pair. Batches cleanly across many models (O(N) subprocess calls, not O(N²)) -- this is the tool for screening large CM sets, e.g. novel motifs against all of Rfam.
- **`minmaxcmcompare`** -- the worst-case counterpart. Approximates CMCompare's exact MaxiMin score (the single sequence that scores as well as possible against *both* models at once) via an evolutionary search: seed a population from real `cmemit` samples, then mutate/crossover/select using `cmsearch` bit scores as the fitness function. Answers "does a single dangerous sequence exist that would confuse both models", which the distributional statistic can't. Costs meaningfully more per pair and doesn't batch the way the bootstrap tool does -- best used to confirm a handful of pairs `bitscoresumdifference` already flagged as significant, not to screen from scratch.

## Build

```
mkdir -p bin
go build -o bin/bitscoresumdifference src/bitscoresumdifference.go
go build -o bin/minmaxcmcompare src/minmaxcmcompare.go
```

Note both files are `package main` with their own `func main`, living in the same `src/` directory as the shared packages (`src/infernal`, `src/pipeline`, etc.). 
This means `go build ./...` / `go vet ./...` will fail specifically on `src` itself ("main redeclared") -- that's expected, not a bug: build and vet each program by naming its file explicitly, as above. 
The shared packages underneath are ordinary Go packages and build/vet/test normally:

```
go vet ./src/infernal/... ./src/matrix/... ./src/output/... ./src/pipeline/... ./src/seqpool/... ./src/shuffle/... ./src/stats/...
go test ./src/infernal/... ./src/matrix/... ./src/output/... ./src/pipeline/... ./src/seqpool/... ./src/shuffle/... ./src/stats/...
```

Zero external Go dependencies (stdlib only), so this builds as a static binary anywhere Go is available.

## Algorithm Approach: bitscoresumdifference

For a pair of covariance models CM_A and CM_B: does CM_B recognize sequences typical of CM_A (and vice versa) more than expected from base composition alone?

```
INPUT: CM_A, CM_B, n_pool (real samples per model), n_shuffles (bootstrap reps)

# --- Step 1: real sequence pools, one cmemit call per model ---
poolA_real = cmemit -N n_pool CM_A
poolB_real = cmemit -N n_pool CM_B

# --- Step 2: score real pools against both models (batched cmsearch, 4 calls total) ---
scoreA_on_A = cmsearch CM_A poolA_real   # --max -T <low>, so every sequence is scored, not filtered
scoreB_on_A = cmsearch CM_B poolA_real   # does B recognize A's sequences?
scoreB_on_B = cmsearch CM_B poolB_real
scoreA_on_B = cmsearch CM_A poolB_real   # does A recognize B's sequences?

# --- Step 3: directional overlap statistics ---
stat_A_into_B = mean(scoreB_on_A) - mean(scoreA_on_A)
stat_B_into_A = mean(scoreA_on_B) - mean(scoreB_on_B)
# near 0        -> near-total mutual overlap
# very negative -> little cross-recognition

# --- Step 4: null distributions via dinucleotide-preserving shuffles ---
for rep in 1..n_shuffles:
    shuffledA = [dinuc_shuffle(seq) for seq in poolA_real]   # Altschul-Erikson algorithm
    shuffledB = [dinuc_shuffle(seq) for seq in poolB_real]   # preserves dinuc freq, not naive shuffle

    null_B_on_shuffledA[rep] = mean(cmsearch CM_B shuffledA)
    null_A_on_shuffledB[rep] = mean(cmsearch CM_A shuffledB)

# --- Step 5: empirical p-values ---
p_A_into_B = fraction(null_B_on_shuffledA >= mean(scoreB_on_A))
p_B_into_A = fraction(null_A_on_shuffledB >= mean(scoreA_on_B))

# --- Step 6: interpretation ---
# both directions significant, stat_* near 0  -> models are essentially redundant
# only one direction significant              -> asymmetric containment (subset/superset)
# neither significant                          -> models capture independent signal
```

**Scoring via `cmsearch`, not `cmscore`.** `cmsearch --tblout` has a stable, documented 18-column layout (see the Infernal user guide's "tabular output" section), unlike `cmscore`'s plain-text table, whose column layout has shifted across Infernal releases. Run with `--max` (disables acceleration filters, so every sequence -- including deliberately weak dinucleotide-shuffled ones -- is fully scored by the CM rather than possibly discarded by a pre-filter) and a low `-T` raw bit-score threshold (rather than an E-value threshold, since E-value depends on pool size and would shift differently between the real and shuffled passes).

**Why dinucleotide shuffling.** CM bit scores are already partially composition-corrected internally, but not perfectly, especially for GC-skewed motifs. Dinucleotide shuffling (Altschul-Erikson, Eulerian-circuit-based) preserves local composition and dinucleotide frequency while destroying base-pairing complementarity and higher-order structure, isolating genuine shared structural/sequence signal from shared compositional bias.

**Why directional statistics.** Scoring poolA-vs-B and poolB-vs-A separately, rather than collapsing to one symmetric number, detects asymmetric containment -- one model as a specialised subset of a broader one -- which a single similarity number can't distinguish from true bidirectional redundancy.

**Cost.** Subprocess calls scale as `O(n_shuffles)` per pair (each replicate is one batched call over the full pool, not one call per sequence), and model pools/scores are built once and reused across all pairs, keeping total calls linear in the number of models rather than quadratic in the number of pairs.

**Relationship to exact/worst-case methods.** This is an approximation of the *distributional* question. For the *worst-case* question (does even one dangerous sequence exist), see `minmaxcmcompare` below, or [CMCompare](https://github.com/RaumZeit/CMCompare) itself for the exact analytic answer.

## Algorithm Approach: minmaxcmcompare

```
INPUT: CM_A, CM_B, pop_size, generations, survival_frac

# --- Phase 1: seed a population from real sequences ---
population = [consensus_A, consensus_B] + cmemit -N k CM_A + cmemit -N k CM_B

# --- Phase 2: evolutionary refinement ---
for generation in 1..generations:
    scores_A = cmsearch CM_A population     # one batched call
    scores_B = cmsearch CM_B population      # one batched call
    fitness  = min(scores_A, scores_B)

    survivors = top_fraction(population, fitness)
    best_ever = update(best_ever, max(fitness))

    population = survivors + mutate_and_crossover(survivors, pop_size)

# empirical_maximin_score ~= best_ever.fitness
# empirical_maximin_seq   ~= best_ever.sequence
```

Independent sampling struggles to find the true maximin sequence because it's often a compromise point that's somewhat atypical of *both* models individually -- a low-probability region under either model's own emission distribution. Turning it into a directed search (using `cmsearch`'s bit score as a fitness oracle) converges toward that region much faster than blind sampling, at the cost of two `cmsearch` calls per generation per pair (still batched across the whole population, just not batched across pairs the way `bitscoresumdifference` is).

The actual ceiling here is CMCompare's exact analytic solution (dynamic programming over the product of the two CMs' state graphs) -- this tool is a fast, good-enough approximation for triage, not a replacement for it.

## Experiments

The core validation question: **does either tool's similarity statistic recover known RNA family relationships?** Rfam's clan annotations provide ground truth -- clans group families believed to share a common evolutionary origin and structural/functional similarity, despite often being too divergent in primary sequence to merge into one family.

Planned evaluation:
- **Within-clan vs. between-clan scoring.** Compute pairwise statistics for all CM pairs within each Rfam clan (from `data/Rfam.clanin`), compare against a background distribution of between-clan pairs. Within-clan pairs should show systematically lower p-values / higher similarity ratios.
- **Discrimination performance.** Using clan co-membership as the binary ground-truth label, summarise separability with ROC/precision-recall curves and AUC.
- **Failure case inspection.** Clans that fail to separate cleanly are informative either way -- phylogenetic/functional relatedness without detectable CM-level signal, or a real limitation of the sampling/shuffling approach.
- **Cross-check between the two tools, and against real CMCompare**, on a manageable subset of clans, to confirm the distributional and worst-case statistics broadly agree in ranking even when their absolute scales differ.

## data/

See [`data/README.md`](data/README.md) for `Rfam.cm` / `Rfam.clanin` provenance and re-fetch instructions.

## envs/

`envs/cmsimile_env.yaml` describes the conda environment (Infernal, providing `cmemit`/`cmsearch`/`cmstat`/`cmfetch`; Go toolchain). Build with:

```
conda env create -f envs/cmsimile_env.yaml
conda activate cmsimile
```

## results/

Output TSVs from both tools land here by convention (`-out results/...`); not version-controlled by default beyond the placeholder.

## src/

- `bitscoresumdifference.go`, `minmaxcmcompare.go` -- the two standalone CLI programs (see above).
- `infernal/` -- `cmemit`/`cmsearch` subprocess wrappers, `--tblout` parsing, and multi-model CM file expansion (`cmstat`/`cmfetch`, for pointing either tool directly at `data/Rfam.cm`).
- `seqpool/` -- the tagged-sequence pool abstraction shared by both tools.
- `matrix/`, `stats/`, `pipeline/` -- the (sequence x model) score matrix, pairwise statistics, and N-way orchestration used by `bitscoresumdifference`.
- `shuffle/` -- the Altschul-Erikson dinucleotide-preserving shuffle.
- `output/` -- TSV writer.

This tool depends on and is intended to be used alongside:

- Nawrocki, E.P. and Eddy, S.R. (2013) *Infernal 1.1: 100-fold faster RNA homology searches*, Bioinformatics 29:2933-2935.
- Kalvari, I. et al. *Rfam*, Nucleic Acids Research. (Cite the release notes/publication matching the Rfam v15.1 data version in `data/` -- check rfam.org for the specific v15.1 citation.)
- Höner zu Siederdissen, C. and Hofacker, I.L. (2010) *Discriminatory power of RNA family models*, Bioinformatics 26:i453-i459. (CMCompare, prior art for exact MaxiMin comparison.)

## tests/

`tests/fixtures/` holds fake stand-ins for `cmemit`/`cmsearch`/`cmstat`/`cmfetch` (Python scripts producing plausible, correctly-shaped output) so the full subprocess-and-parsing pipeline can be smoke-tested without a real Infernal install:

```
go build -o bin/bitscoresumdifference src/bitscoresumdifference.go
./bin/bitscoresumdifference \
  -cm modelA=tests/cms/modelA.cm -cm modelB=tests/cms/modelB.cm -cm modelC=tests/cms/modelC.cm \
  -npool 30 -nshuffles 30 \
  -cmemit-bin tests/fixtures/fake_cmemit -cmsearch-bin tests/fixtures/fake_cmsearch \
  -out results/smoke_test.tsv
```

`tests/cms/` holds placeholder `.cm` files for the above (not real covariance models -- just fixtures the fake binaries key off of by filename).
