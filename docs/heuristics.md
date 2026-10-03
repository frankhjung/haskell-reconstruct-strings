# Assembly Dynamics and Parameter Heuristics

[← Back to README](../README.md)

This guide is intended for users calibrating the `-m` minimum overlap
threshold or the strand generator parameters for their specific dataset.
It provides mathematical models for collision probability, coverage bounds,
and worked parameter examples for both nucleotide (`ATGC`) and text (`A-Z`)
alphabet use cases. For parameter definitions, see the [Domain Glossary]
[glossary].

Assembler performance depends on the interaction between alphabet size,
fragment length, minimum overlap threshold, and sequencing coverage.

## Alphabet Size and Collision Probability

The probability of an incidental prefix-suffix match of length $k$ between two
independent random sequences over an alphabet $\Sigma$ scales exponentially:

$$P(\text{overlap} \ge k) = \left(\frac{1}{|\Sigma|}\right)^k$$

Comparing a 4-character nucleotide alphabet (`ATGC`) with a 26-character
alphabet (`A-Z`):

- **$k = 2$**:
  - `ATGC`: $(\frac{1}{4})^2 = \frac{1}{16}$ = 6.25%
  - `A-Z`: $(\frac{1}{26})^2 = \frac{1}{676}$ ~ 0.148%
    ($\approx 42\times$ rarer)
- **$k = 3$**:
  - `ATGC`: $(\frac{1}{4})^3 = \frac{1}{64}$ ~ 1.56%
  - `A-Z`: $(\frac{1}{26})^3$ ~ 0.0057% ($\approx 274\times$ rarer)
- **$k = 4$**:
  - `ATGC`: $(\frac{1}{4})^4 = \frac{1}{256}$ ~ 0.391%
  - `A-Z`: $(\frac{1}{26})^4$ ~ 0.000219% ($\approx 1{,}785 \times$ rarer)

## Impact on Assembly Behaviour

- **Unrelated Random Noise**: When assembling random fragments without a shared
  reference sequence, `ATGC` collapses fragments into spurious contigs due to
  frequent coincidental matches. Conversely, `A-Z` fragments rarely share
  accidental overlaps, causing greedy reduction in [`Assembler.hs`]
  [src-assembler] to halt immediately without merges. The apparent
  assembly of small alphabets is an illusion caused by chimeric joins.
- **Sensitivity to Coverage Gaps**: Over `A-Z`, 4-mers are statistically
  unique ($1$ in $456{,}976$). If physical coverage has a gap where adjacent
  fragments overlap by less than $m$, the assembler halts and outputs
  fragmented contigs. In `ATGC`, chance 4-mer matches across distant regions
  can falsely bridge coverage gaps, resulting in scrambled assemblies.
- **Repeat Ambiguity**: In `ATGC`, short sequences rapidly exhaust unique
  permutations ($4^4 = 256$), causing repeated $k$-mers that trap greedy
  heuristics in local optima. In `A-Z`, high information entropy virtually
  eliminates repeats in moderate-length sequences.

## Minimum Overlap Lower Bound

For a pool of $N$ fragments, there are $N(N - 1)$ ordered pairwise comparisons
in `findBestOverlap` within [`Assembler.hs`][src-assembler]. To ensure the
expected number of false-positive pairwise matches across the dataset is less
than 1:

$$E[\text{spurious pairs}] \approx N^2 \cdot |\Sigma|^{-m} < 1$$

$$\implies m_{\text{min}} = \lceil 2 \log_{|\Sigma|} N \rceil$$

- For $N = 1{,}000$ in `ATGC`: $m_{\text{min}} = \lceil 2 \log_4 1000 \rceil
  = 10$ bases.
- For $N = 1{,}000$ in `A-Z`: $m_{\text{min}} = \lceil 2 \log_{26} 1000
  \rceil = 5$ characters.

## Minimum Overlap Upper Bound

Under the Lander–Waterman model of sequencing, adjacent fragments in the
reference must overlap by at least $m$ to be detected. Given average fragment
length $\bar{L}$, the effective coverage $C_{\text{eff}}$ is:

$$C_{\text{eff}} = C \cdot \left(1 - \frac{m}{\bar{L}}\right)$$

where $C = \frac{N \cdot \bar{L}}{G}$ is nominal physical coverage for a target
of length $G$. As $m \to \bar{L}$, $C_{\text{eff}} \to 0$, causing exponential
fragmentation into separate contigs.

- **Rule of thumb**: Keep $m \le 0.5 \cdot \bar{L}$ to preserve at least 50% of
  nominal coverage.

## Fragment Length Requirements

- **Repeat Resolution**: To resolve repetitive elements, fragment length must
  strictly exceed the longest repeat length: $L > R_{\max}$.
- **Length Distribution**: Using variable fragment lengths (e.g. `--min` and
  `--max` in [`make-strands.sh`][make-strands]) breaks tie-breaking edge cases
  during greedy selection.

## Parameter Selection Guidelines

- **Minimum Overlap ($m$)**:
  - Set $m = \max\left(\lceil 2 \log_{|\Sigma|} N \rceil, \; 0.3 \cdot
    \bar{L}\right)$.
  - Typical range: 30% to 50% of average fragment length $\bar{L}$.
- **Fragment Length ($L$)**:
  - Ensure $L > R_{\max}$ (longer than the longest repeat).
  - Typical range: $2 \times$ to $3 \times$ the overlap threshold $m$.
- **Coverage Depth ($C$)**:
  - Target $C_{\text{eff}} = C(1 - \frac{m}{\bar{L}}) \ge 10$.
  - Nominal physical coverage: $15 \times$ to $30 \times$.

## Worked Configuration Examples

The following calibrated configurations illustrate parameter selection for a
target sequence of length $G = 1{,}000$ units and $N = 1{,}000$ fragments:

- **Genomic fragments (`ATGC`, $|\Sigma| = 4$)**:
  Incidental overlap collisions scale as $(\frac{1}{4})^k$, requiring
  higher overlap thresholds and longer fragments to prevent false joins.
  - Fragment length range: $20 \text{--} 40$ bp (mean $\bar{L} = 30$ bp).
  - Minimum overlap: $m = 10$ bases
    ($m_{\min} = \lceil 2 \log_4 1000 \rceil = 10$,
    $E[\text{spurious}] \approx 0.95 < 1$).
  - Coverage: nominal $C = 30 \times$, effective
    $C_{\text{eff}} = 30(1 - \frac{10}{30}) = 20 \times \ge 10 \times$.
  - Sample commands:

    ```bash
    # Prepare 1000 fragments from a reference genome of 1000 bp
    # with length between 20 and 40 bp
    ./scripts/make-strands.sh \
        -i reference_genome.txt -m 20 -M 40 -n 1000 > fragments.txt
    # Reconstruct the fragments with minimum overlap of 10 bp
    cabal run reconstruct-strings -- -m 10 -f fragments.txt
    ```

- **String character fragments (`A-Z`, $|\Sigma| = 26$)**:
  Higher entropy ($(\frac{1}{26})^k$) drastically reduces collision
  probability, permitting smaller thresholds and shorter fragments without
  chimera formation.
  - Fragment length range: $10 \text{--}20$ chars (mean $\bar{L} = 15$ chars).
  - Minimum overlap: $m = 5$ characters
    ($m_{\min} = \lceil 2 \log_{26} 1000 \rceil = 5$,
    $E[\text{spurious}] \approx 0.084 \ll 1$).
  - Coverage: nominal $C = 15 \times$, effective
    $C_{\text{eff}} = 15(1 - \frac{5}{15}) = 10 \times \ge 10 \times$.
  - Sample commands:

    ```bash
    # Prepare 1000 fragments from a reference text of 1000 characters
    # with length between 10 and 20 characters
    ./scripts/make-strands.sh \
        -i reference_text.txt -m 10 -M 20 -n 1000 > strands.txt
    # Reconstruct the fragments with a minimum overlap of 5 characters
    cabal run reconstruct-strings -- -m 5 -f strands.txt
    ```

- **Comparative Summary**:
  - Alphabet size: $4$ (`ATGC`) versus $26$ (`A-Z`).
  - Overlap threshold: $m = 10$ bases versus $m = 5$ characters ($E < 1$).
  - Fragment length: $20 \text{--}40$ bp versus $10 \text{--}20$ characters.
  - Nominal coverage: $30 \times$ versus $15 \times$ physical depth.
  - Effective coverage: $20 \times$ versus $10 \times$ Lander–Waterman depth.

## See Also

- [Documentation Index](README.md)
- [Domain Glossary][glossary]
- [Algorithm Specification][req-001]
- [Reconstructing DNA from Short Fragments][dna-doc]

[dna-doc]: reconstructing-complete-dna-strand-from-short-fragments.md
[glossary]: GLOSSARY.md
[make-strands]: ../scripts/make-strands.sh
[req-001]: REQ-001-functional-greedy-overlap-assembler.md
[src-assembler]: ../src/Assembler.hs
