# cm-simile

Fast, statistically defensible estimation of sequence/structure similarity between Infernal covariance models.

Two standalone Go programs, both living in `src/`, sharing the same package infrastructure (`cmemit`/`cmsearch` wrappers, multi-model CM file expansion, dinucleotide shuffling):

- **`bitscoresumdifference`** -- the primary, fast screening tool. Samples real sequences from each model, cross-scores them, and tests the observed directional median bit-score difference against a dinucleotide-shuffled null distribution for a bootstrap p-value per pair. Batches cleanly across many models (O(N) subprocess calls, not O(N²)) -- this is the tool for screening large CM sets, e.g. novel motifs against all of Rfam.
- **`minmaxcmcompare`** -- the worst-case counterpart. Approximates CMCompare's exact MaxiMin score (the single sequence that scores as well as possible against *both* models at once) via an evolutionary search: seed a population from real `cmemit` samples, then mutate/crossover/select using `cmsearch` bit scores as the fitness function. Answers "does a single dangerous sequence exist that would confuse both models", which the distributional statistic can't. Costs meaningfully more per pair and doesn't batch the way the bootstrap tool does -- best used to confirm a handful of pairs `bitscoresumdifference` already flagged as significant, not to screen from scratch.
- **`ksmwcompare`** -- a simpler distributional alternative to `bitscoresumdifference`'s adaptive reservoir: one real pool and one dinucleotide-shuffled null pool per model (no growth loop, no convergence diagnostic), compared directly via two-sample Kolmogorov-Smirnov (any distributional difference, including shape -- e.g. a real pool that's a mixture of strong- and weak-signal sequences) and Mann-Whitney U (location shift, matching `bitscoresumdifference`'s own hypothesis). Trades some run-to-run stability in the exact reported numbers for a much smaller, easier-to-reason-about system; `-n-null-pools` lets you pool a few independent shuffle instantiations together as cheap insurance against an unlucky single draw.

## Build

```
mkdir -p bin
go build -o bin/bitscoresumdifference src/bitscoresumdifference.go
go build -o bin/minmaxcmcompare src/minmaxcmcompare.go
go build -o bin/ksmwcompare src/ksmwcompare.go
```

Note both files are `package main` with their own `func main`, living in the same `src/` directory as the shared packages (`src/infernal`, `src/pipeline`, etc.). This means `go build ./...` / `go vet ./...` will fail specifically on `src` itself ("main redeclared") -- that's expected, not a bug: build and vet each program by naming its file explicitly, as above. The shared packages underneath are ordinary Go packages and build/vet/test normally:

```
go vet ./src/infernal/... ./src/matrix/... ./src/output/... ./src/pipeline/... ./src/seqpool/... ./src/shuffle/... ./src/stats/... ./src/combine/...
go test ./src/infernal/... ./src/matrix/... ./src/output/... ./src/pipeline/... ./src/seqpool/... ./src/shuffle/... ./src/stats/... ./src/combine/...
```

Zero external Go dependencies (stdlib only), so this builds as a static binary anywhere Go is available.

## Algorithm Approach: bitscoresumdifference

For a pair of covariance models CM_A and CM_B: does CM_B recognize sequences typical of CM_A (and vice versa) more than expected from base composition alone?

```
INPUT: CM_A, CM_B, n_pool (real samples per model), reservoir config, n_shuffles (bootstrap resamples)

# --- Step 1: real sequence pools, one cmemit call per model ---
poolA_real = cmemit -N n_pool CM_A
poolB_real = cmemit -N n_pool CM_B

# --- Step 2: score real pools against both models (batched cmsearch, 4 calls total) ---
scoreA_on_A = cmsearch CM_A poolA_real   # --max -T <low>, so every sequence is scored, not filtered
scoreB_on_A = cmsearch CM_B poolA_real   # does B recognize A's sequences?
scoreB_on_B = cmsearch CM_B poolB_real
scoreA_on_B = cmsearch CM_A poolB_real   # does A recognize B's sequences?

# --- Step 3: directional overlap statistics ---
stat_A_into_B = median(scoreB_on_A)
stat_B_into_A = median(scoreA_on_B)
# near 0        -> near-total mutual overlap
# very negative -> little cross-recognition

# --- Step 4: build ONE null reservoir per model, ONCE, shared across every pair ---
for each model M:
    reservoir_M = []
    round = 0
    while not converged(M):
        target_size = reservoir_min * 2^round
        new_shuffled = [dinuc_shuffle(seq) for seq in poolM_real]   # repeated with fresh
                                                                      # seeds until target_size reached
        reservoir_M += new_shuffled
        score reservoir_M's NEW additions against EVERY model (batched, same trick as poolM_real)
        # convergence: bootstrap-resample the reservoir at n_pool-sized draws,
        # track the 5th-percentile of resampled means per target model;
        # stop growing once that stabilizes within tolerance, or reservoir_max is hit
        round += 1

# --- Step 5: per pair, draw null replicates from the ALREADY-SCORED reservoir -- no more cmsearch calls ---
null_B_on_shuffledA[1..n_shuffles] = bootstrap-resample means of size n_pool from reservoir_A's
                                      precomputed scores against CM_B
null_A_on_shuffledB[1..n_shuffles] = likewise, reservoir_B's scores against CM_A

# --- Step 6: empirical p-values ---
p_A_into_B = fraction(null_B_on_shuffledA >= median(scoreB_on_A))
p_B_into_A = fraction(null_A_on_shuffledB >= median(scoreA_on_B))

# --- Step 7: interpretation ---
# both directions significant, stat_* near 0  -> models are essentially redundant
# only one direction significant              -> asymmetric containment (subset/superset)
# neither significant                          -> models capture independent signal
```

**Why a reservoir instead of shuffling fresh per replicate.** The original design re-shuffled and re-scored the whole pool for every one of `n_shuffles` replicates, making replicate count directly multiply `cmsearch` call count (`O(n_shuffles x models)`). Building one large, adaptively-sized reservoir per model *once*, scoring it once (batched across all models, same trick as the real pool), and then drawing `n_shuffles` bootstrap resamples from that already-scored reservoir *in memory* collapses null-generation cost back to roughly `O(models x log(reservoir_size))` -- independent of `n_shuffles`, which is now nearly free to increase. See `-reservoir-min`/`-reservoir-max`/`-reservoir-tolerance`/`-resample-size` in `bitscoresumdifference -h`.

**Why the reservoir size is adaptive, not fixed.** A reservoir that's too small gives an unstable estimate of the null's tail behavior -- exactly the region a p-value threshold depends on -- which would quietly bias significance calls rather than erroring. Growth doubles each round and stops once a 5th-percentile bootstrap diagnostic stabilizes within tolerance (or `-reservoir-max` is hit, which prints a warning rather than silently proceeding on an under-validated reservoir).


**Why dinucleotide shuffling.** CM bit scores are already partially composition-corrected internally, but not perfectly, especially for GC-skewed motifs. Dinucleotide shuffling (Altschul-Erikson, Eulerian-circuit-based) preserves local composition and dinucleotide frequency while destroying base-pairing complementarity and higher-order structure, isolating genuine shared structural/sequence signal from shared compositional bias. Worth knowing as a real limitation: this null controls for *composition*, not for whether a sequence has *any* real secondary structure at all -- a truly unrelated but well-folded ncRNA can still outscore its own shuffled self against an unrelated CM, simply because it has a fold and the shuffled version doesn't. A structure-preserving null (randomize paired-position identity while keeping the pairing itself) would isolate family-specific signal more cleanly; not yet implemented here.

**Why directional statistics.** Scoring poolA-vs-B and poolB-vs-A separately, rather than collapsing to one symmetric number, detects asymmetric containment -- one model as a specialised subset of a broader one -- which a single similarity number can't distinguish from true bidirectional redundancy.

**Matrix output.** `-dump-matrices DIR` writes the full (sequence x model) score matrices for both the real pool and the null reservoir as long-format TSVs (`sequence_id, origin, model, score`) -- e.g. for plotting the null distribution's actual shape rather than only its summary statistics.

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

**Parallelism.** Pairs run concurrently, bounded by `-threads` (default 4). Within a single pair's refinement, each generation's scoring against model A and model B also runs concurrently -- so up to roughly `2 x threads` `cmsearch` processes can be in flight at once; reduce `-threads` if that's too much for your allocation.

**Matrix output.** `-dump-matrix PATH` writes every candidate's scoreA/scoreB from each reported pair's *final* generation (long format: `model_a, model_b, candidate_id, seq, score_a, score_b`) -- the converged population the search actually found, not the full generation-by-generation trajectory.

## Experiments

The core validation question: **does either tool's similarity statistic recover known RNA family relationships?** Rfam's clan annotations provide ground truth -- clans group families believed to share a common evolutionary origin and structural/functional similarity, despite often being too divergent in primary sequence to merge into one family.

Planned evaluation:
- **Within-clan vs. between-clan scoring.** Compute pairwise statistics for all CM pairs within each Rfam clan (from `data/Rfam.clanin`), compare against a background distribution of between-clan pairs. Within-clan pairs should show systematically lower p-values / higher similarity ratios.
- **Discrimination performance.** Using clan co-membership as the binary ground-truth label, summarise separability with ROC/precision-recall curves and AUC.

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

- `bitscoresumdifference.go`, `minmaxcmcompare.go`, `ksmwcompare.go` -- the three standalone CLI programs (see above).
- `nptest/` -- two-sample Kolmogorov-Smirnov and Mann-Whitney U test implementations, used by `ksmwcompare.go`.
- `infernal/` -- `cmemit`/`cmsearch` subprocess wrappers, `--tblout` parsing, and multi-model CM file expansion (`cmstat`/`cmfetch`, for pointing either tool directly at `data/Rfam.cm`).
- `seqpool/` -- the tagged-sequence pool abstraction shared by both tools.
- `matrix/`, `stats/`, `pipeline/` -- the (sequence x model) score matrix, pairwise statistics, and N-way orchestration used by `bitscoresumdifference`.
- `shuffle/` -- the Altschul-Erikson dinucleotide-preserving shuffle.
- `output/` -- TSV writer.

This tool depends on and is intended to be used alongside:

- Nawrocki, E.P. and Eddy, S.R. (2013) *Infernal 1.1: 100-fold faster RNA homology searches*, Bioinformatics 29:2933-2935.
- Kalvari, I. et al. *Rfam*, Nucleic Acids Research. (Cite the release notes/publication matching the Rfam v15.1 data version in `data/` -- check rfam.org for the specific v15.1 citation.)

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


## Future Direction

Implement with Hidden Markov Models (HMMs) with hmmer equivalents to cm - `hmmemit`, `hmmsearch`, `hmmstat`, `hmmfetch`.
Shouldn't be too much work to implement with existing go module dynamically taking executable paths and the command line arguments being almost identical for HMM and CM tools but would need to be tested and validated.