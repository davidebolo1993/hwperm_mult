# hwperm-mult

Standalone C++ CLI tool for fast multi-allelic autosomal Hardy-Weinberg permutation testing from a 2-line allele file.

## Origin and attribution

This tool is adapted from the algorithmic workflow used in the R package **HardyWeinberg**.
In particular, it mirrors the autosomal multi-allelic path based on:

- `AllelesToTriangular` (genotype triangular matrix construction)
- `dlevene` (Levene-Haldane probability calculation)
- `HWPerm.mult` (permutation-test principle)

References:

- CRAN package: https://CRAN.R-project.org/package=HardyWeinberg
- Package source repository (this codebase): `HardyWeinberg`
- JSS paper: Graffelman J. (2015), *Exploring Diallelic Genetic Markers: The HardyWeinberg Package*, Journal of Statistical Software, 64(3), doi:10.18637/jss.v064.i03

This repository provides a standalone C++ implementation for the same core idea in a CLI form.

## What it does

Given a file with exactly 2 lines:

- line 1: allele 1 for each individual
- line 2: allele 2 for each individual

the tool:

1. builds the lower-triangular genotype count matrix (same idea as `AllelesToTriangular`)
2. computes the Levene-Haldane sample probability (`dlevene`-style statistic)
3. runs a Monte Carlo permutation test (same principle as `HWPerm.mult` autosomal branch)
4. prints:
   - `p-value`
   - `observed_statistic`
   - `k_alleles`
   - `n_genotypes`
   - optional genotype matrix

## Algorithmic principle

Let `x` be the lower-triangular genotype count matrix.

- Allele counts are `n_i = rowSums(x) + colSums(x)`.
- For sample size `n = sum(x)`, the observed Levene-Haldane probability is:

  `P = exp( lgamma(n+1) + sum(lgamma(n_i+1)) + h*log(2) - lgamma(2n+1) - sum(lgamma(g_j+1)) )`

  where:
  - `h` is total heterozygote count (strict lower triangle),
  - `g_j` are genotype counts in the lower triangle.

Permutation test:

1. Expand allele pool from `n_i`.
2. Shuffle allele pool.
3. Pair consecutive alleles into genotypes.
4. Rebuild triangular genotype counts.
5. Recompute `P*`.
6. p-value is fraction of permutations with `P* <= P` (with near-equality tolerance `eps`).

This is a Monte Carlo estimate, so p-values vary slightly across runs unless a fixed seed is used.

## Build

### Option A: script

```bash
./scripts/build.sh
```

### Option B: direct compile

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic src/hwperm_mult.cpp -o hwperm_mult
```

## Input format

Two lines, tab-separated by default.

Example:

```text
0	1	2	0	1
0	1	2	1	1
```

Missing tokens (`NA`, `NaN`, `.`, empty) are ignored pairwise.

## Usage

```bash
./hwperm_mult --input data.tsv --nperm 20000
```

From stdin:

```bash
cat data.tsv | ./hwperm_mult --nperm 20000 --quiet
```

## CLI options

- `--input <path>`: file with 2 lines (optional if using stdin)
- `--nperm <int>`: number of permutations (default: `20000`)
- `--eps <double>`: floating-point near-equality tolerance (default: `1e-10`)
- `--sep <token>`: `tab|space|comma|semicolon|<char>` (default: `tab`)
- `--seed <int>`: fixed RNG seed for reproducible runs
- `--quiet`: prints only key-value result lines
- `--print-matrix`: prints full genotype matrix
- `--help`: help message

## Notes

- `observed_statistic` is deterministic for fixed input.
- `p-value` is stochastic unless `--seed` is set.
- Even with the same numeric seed, R and C++ can produce slightly different p-values because RNG engines differ.
