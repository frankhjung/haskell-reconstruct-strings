# Functional Greedy Overlap Assembler Specification

## 1. Overview and Problem Statement

The assembler accepts an unordered collection of reads and constructs one or
more contigs by repeatedly merging reads whose suffixes overlap prefixes. This
is a greedy approximation of *de novo* sequence assembly designed to be
deterministic, pure, total, and mathematically rigorous.

The implementation must not depend on mutable state, imperative loops, or
in-place list modification. The algorithm is expressed as a pure functional
reduction over immutable values in Haskell.

The system produces a deterministic list of contigs from an input list of reads.
A correct assembly is one that:

- validates input parameters totally without runtime exceptions,
- deduplicates identical reads to a single representative,
- removes reads that contribute no new sequence information (containment),
- greedily merges the candidate pair with the best valid overlap at each step,
- resolves ties deterministically using a strict total order,
- eliminates reads that become contained within newly merged contigs,
- stops when no valid overlap remains above the threshold,
- and emits the remaining contigs in a canonical, permutation-invariant order.

This specification defines the required behaviour, type contracts, edge-case
handling, and output guarantees for a Haskell implementation.

## 2. Scope and Non-Goals

This document specifies:

- the accepted input model and type contracts,
- exact suffix-prefix overlap semantics,
- deterministic candidate selection and tie-breaking,
- pre-processing and dynamic containment elimination,
- the final output contract and canonical contig ordering,
- and expected computational complexity bounds.

Non-goals include:

- handling reverse complements or double-stranded DNA heuristics,
- inexact matching or sequencing error tolerance,
- scaffolding or gap filling across disjoint components,
- de Bruijn graph representations or Eulerian path algorithms.

## 3. Data Model

The domain model is defined using Haskell newtypes and sum types over strict
`Text`:

```haskell
module Assembler.Types
  ( Read(..)
  , Contig(..)
  , OverlapLength
  , OverlapCandidate(..)
  , AssemblyError(..)
  ) where

import Data.Text (Text)

newtype Read = Read { unRead :: Text }
  deriving stock (Eq, Ord, Show)

newtype Contig = Contig { unContig :: Text }
  deriving stock (Eq, Ord, Show)

type OverlapLength = Int

data OverlapCandidate = OverlapCandidate
  { prefixRead  :: !Read
  , suffixRead  :: !Read
  , matchLength :: !OverlapLength
  }
  deriving stock (Eq, Ord, Show)

data AssemblyError
  = InvalidMinOverlap !Int
  | EmptyReadEncountered
  deriving stock (Eq, Show)
```

### Required Interpretation

- `Read`: An immutable sequence of characters. Must be non-empty.
- `Contig`: A contiguous sequence produced by merging overlapping reads, or an
  unmerged singleton read.
- `OverlapLength`: A non-negative integer representing the length of an exact
  suffix-prefix match.
- A valid overlap is an exact match where a suffix of `prefixRead` matches a
  prefix of `suffixRead`. The overlap length must be strictly greater than zero,
  at least `min_overlap`, and strictly less than the length of the longer read.

## 4. Functional Requirements

### 4.1 Required Input Assumptions

The entry point `assemble` accepts:

- an input list of `Read` values: `[Read]`,
- a minimum overlap threshold: `min_overlap :: Int`.

#### Validation Rules

- If `min_overlap < 1`, the function must return
  `Left (InvalidMinOverlap min_overlap)`.
- If any read in the input has length 0, the function must return
  `Left EmptyReadEncountered`.
- The input list is treated as an unordered multiset.

### 4.2 `calculateOverlap`

`calculateOverlap r1 r2 minOverlap` computes the maximum overlap length `k` such
that:

- `k >= minOverlap`,
- the suffix of length `k` in `r1` equals the prefix of length `k` in `r2`,
- and `k < max (length r1) (length r2)`.

If no such `k` exists, the function returns `0`.

#### Required Semantics

- The overlap is directed: `calculateOverlap r1 r2 minOverlap` is not
  necessarily equal to `calculateOverlap r2 r1 minOverlap`.
- The search proceeds from `min (length r1) (length r2)` down to `minOverlap`.
- The returned value is the largest valid overlap meeting the criteria, or `0`.

```haskell
calculateOverlap :: Read -> Read -> Int -> OverlapLength
calculateOverlap (Read r1) (Read r2) minOverlap =
  let maxPossible = min (T.length r1) (T.length r2)
      candidates =
        [ len | len <- [maxPossible, maxPossible - 1 .. minOverlap]
              , T.takeEnd len r1 == T.take len r2
              , len < max (T.length r1) (T.length r2) ]
  in case candidates of
       (best : _) -> best
       []         -> 0
```

### 4.3 `mergePair`

`mergePair prefix suffix overlapLen` concatenates the prefix read with the
non-overlapping suffix remainder:

```haskell
mergePair :: Read -> Read -> OverlapLength -> Read
mergePair (Read p) (Read s) overlapLen =
  Read (p <> T.drop overlapLen s)
```

No other characters may be inserted, dropped, or modified.

### 4.4 `findBestOverlap`

`findBestOverlap pool minOverlap` inspects every ordered pair of distinct reads
`(a, b)` from `pool` where `a /= b` and computes
`calculateOverlap a b minOverlap`.

#### Output Contract

- If no valid pair satisfies `overlap >= minOverlap`, return `Nothing`.
- Otherwise, return `Just candidate` where `candidate` is the maximal candidate
  under the deterministic tie-breaking order.

#### Deterministic Tie-Breaking Order

Candidate pairs are ordered strictly by the following total ordering:

1. Greater `matchLength` takes precedence.
2. If `matchLength` is equal, lexicographically smaller `prefixRead` takes
   precedence.
3. If `prefixRead` is equal, lexicographically smaller `suffixRead` takes
   precedence.

Because the candidate pool contains pairwise distinct reads, Rules 1 to 3
guarantee a strict total order over all candidate pairs.

### 4.5 `filterContainedReads`

A read `target` is contained within another read `other` if:

- `target /= other`, and
- `unRead target` is an infix (proper substring) of `unRead other`.

Pre-processing removes exact duplicates (retaining one unique representative)
and drops any read that is a proper substring of another read in the pool:

```haskell
filterContainedReads :: [Read] -> [Read]
filterContainedReads reads =
  let uniqueReads = nub reads
  in [ r | r <- uniqueReads
         , not (any (\o -> r /= o &&
                           unRead r `T.isInfixOf` unRead o) uniqueReads) ]
```

## 5. Required Assembly Behaviour

### 5.1 High-Level Contract

The primary assembly function has the following signature:

```haskell
assemble :: [Read] -> Int -> Either AssemblyError [Contig]
```

It must guarantee:

- input validation without partial functions,
- initial duplicate deduplication and containment elimination,
- iterative greedy merging of the maximal candidate overlap,
- dynamic elimination of reads engulfed by newly formed contigs,
- and canonical, permutation-invariant ordering of the final contigs.

### 5.2 Reference Reduction Algorithm

```haskell
assemble :: [Read] -> Int -> Either AssemblyError [Contig]
assemble reads minOverlap
  | minOverlap < 1 = Left (InvalidMinOverlap minOverlap)
  | any (T.null . unRead) reads = Left EmptyReadEncountered
  | otherwise =
      let initialPool = filterContainedReads reads
          finalPool = reducePool initialPool minOverlap
          sortedContigs = sortCanonical (map toContig finalPool)
      in Right sortedContigs

reducePool :: [Read] -> Int -> [Read]
reducePool pool minOverlap
  | length pool <= 1 = pool
  | otherwise =
      case findBestOverlap pool minOverlap of
        Nothing -> pool
        Just candidate ->
          let merged = mergePair (prefixRead candidate)
                                 (suffixRead candidate)
                                 (matchLength candidate)
              remaining = [ r | r <- pool
                              , r /= prefixRead candidate
                              , r /= suffixRead candidate ]
              updatedPool = filterContainedReads (merged : remaining)
          in reducePool updatedPool minOverlap
```

### 5.3 Dynamic Containment and Invariants

- The pool is strictly immutable; every step returns a newly allocated list.
- Each merge step decreases the pool size by at least 1 (more if the merge
  engulfs another read).
- If no valid overlap remains, reduction terminates.
- `updatedPool` re-evaluates containment to eliminate any read engulfed within
  the newly merged contig.

### 5.4 Canonical Output Ordering

Because input reads constitute an unordered multiset, the output contig list
must be independent of input list permutations. Final contigs are sorted
canonically:

1. Longer `Contig` sequence length first (descending).
2. Lexicographically smaller sequence first on ties (ascending).

```haskell
sortCanonical :: [Contig] -> [Contig]
sortCanonical = sortBy compareContigs
  where
    compareContigs (Contig a) (Contig b) =
      compare (T.length b) (T.length a) <> compare a b
```

## 6. Edge Cases and Required Behaviour

### 6.1 Empty Input

If `reads` is empty, `assemble reads minOverlap` returns `Right []`.

### 6.2 Singleton Input

If `reads` contains exactly one non-empty read `[r]`, `assemble` returns
`Right [Contig (unRead r)]`.

### 6.3 All Reads Identical

If all input reads are identical (e.g. `["ACGT", "ACGT"]`), duplicates represent
redundant sequencing coverage. Pre-processing deduplicates them to a single
representative, returning `Right [Contig "ACGT"]`.

### 6.4 Contained Reads

Any read that is a proper substring of another read is discarded. This applies
both in pre-processing and dynamically during reduction when two reads merge to
span a third read.

### 6.5 No Valid Overlap

If no pair meets `min_overlap`, greedy reduction terminates immediately. The
remaining reads are wrapped as contigs and sorted in canonical order.

### 6.6 Disjoint Components

If reads originate from disconnected genomic regions, components reduce
independently. All resulting contigs are sorted canonically in the final output.

## 7. Computational Complexity and Invariants

### 7.1 Termination

The pool size strictly decreases by at least 1 per recursive iteration. For an
initial filtered pool of `N` reads, reduction executes at most `N - 1` merges.

### 7.2 Time Complexity

For `N` reads of average length `L`:

- Pairwise overlap evaluation takes $O(N^2 \cdot L)$$.
- At most $N - 1$ reduction stages occur.
- Worst-case runtime of naive reduction is $O(N^3 \cdot L)$.

### 7.3 Space Complexity

Because Haskell is non-strict by default, strict evaluation of contig text and
accumulator lists prevents space leaks. The memory footprint is $O(N \cdot L)$.

### 7.4 Invariants

Throughout all stages, the implementation guarantees:

- purity and referential transparency with total functions,
- absence of contained reads in active pools,
- strict total ordering for overlap candidates,
- and permutation invariance of the final contig list.

## 8. Acceptance Criteria

The implementation is verified if:

1. `assemble` accepts `[Read]` and `Int`, returning
   `Either AssemblyError [Contig]`.
2. Invalid overlap thresholds (`< 1`) or empty reads yield descriptive errors.
3. Identical reads collapse to a single representative contig.
4. Reads contained within other reads or newly merged contigs are eliminated.
5. Ties resolve strictly by `matchLength`, `prefixRead`, then `suffixRead`.
6. Final contig output is sorted canonically and invariant to input permutation.
7. Pure functions and immutable data structures are used exclusively.
8. Property-based tests verify determinism, totality, and containment
   invariants.

## 9. Example Test Cases

### Example A: Single Valid Overlap

```haskell
reads = [Read "ABC", Read "BCD", Read "CDE"]
minOverlap = 2
-- Expected: Right [Contig "ABCDE"]
```

### Example B: No Valid Overlap

```haskell
reads = [Read "ABC", Read "DEF"]
minOverlap = 2
-- Expected: Right [Contig "ABC", Contig "DEF"]
```

### Example C: Contained Read Removed

```haskell
reads = [Read "ACGT", Read "CGT"]
minOverlap = 2
-- Expected: Right [Contig "ACGT"]
```

### Example D: Dynamic Containment After Merge

```haskell
reads = [Read "AAATTT", Read "TTTGGG", Read "ATT"]
minOverlap = 3
-- Expected: Right [Contig "AAATTTGGG"]
-- "ATT" is engulfed by the merged sequence and eliminated.
```

### Example E: Canonical Permutation Invariance

```haskell
assemble [Read "DEF", Read "ABC"] 2 == assemble [Read "ABC", Read "DEF"] 2
-- Both return Right [Contig "ABC", Contig "DEF"]
```

## 10. Summary

This specification provides a total, pure functional blueprint for greedy
overlap sequence assembly in Haskell. By modeling domain entities as strong
newtypes, formalising dynamic containment, and enforcing a strict canonical
ordering on output contigs, the assembler guarantees deterministic,
permutation-invariant results across all valid inputs.

## 11. Architectural Decisions

The following architectural decision records document load-bearing trade-offs
settled during requirements review:

### ADR-1: Redundant Coverage Deduplication

- **Context:** An input may consist solely of identical reads (e.g.
  `["ACGT", "ACGT"]`). An earlier ambiguity suggested either discarding all
  duplicates or preserving a single copy.
- **Decision:** Deduplicate identical reads to a single representative.
  Multiple identical reads represent repeated sequencing coverage of the same
  region, not invalid data.
- **Consequence:** Output for an all-duplicate input is a single contig
  (`Right [Contig "ACGT"]`). Information is preserved without generating false
  self-overlaps.

### ADR-2: Strict Three-Tier Deterministic Tie-Breaking

- **Context:** An earlier specification included a fourth tie-breaking rule
  referencing input list indices. However, because reads in the candidate pool
  are distinct strings, rules 1 to 3 form a complete strict total order.
  Moreover, merged contigs lack input indices.
- **Decision:** Remove the fourth rule. Candidate pairs are ordered solely by:
  1. Greater `matchLength`
  2. Lexicographically smaller `prefixRead`
  3. Lexicographically smaller `suffixRead`
- **Consequence:** Eliminates unnecessary state tracking and index provenance
  across recursive reductions. Purity and totality are preserved.

### ADR-3: Dynamic Containment Elimination

- **Context:** Merging two overlapping reads can create a composite sequence
  that engulfs a third, previously uncontained read in the pool.
- **Decision:** Re-evaluate containment after each merge (`updatedPool`)
  rather than only once during pre-processing.
- **Consequence:** Contained fragments are promptly purged, preventing
  redundant sequence emissions and degenerate zero-gain merge cycles.

### ADR-4: Canonical Permutation Invariance

- **Context:** The input model is an unordered multiset of reads. However,
  naive recursive reduction prepends merged reads, which could produce
  different contig orders for permuted inputs.
- **Decision:** Sort the final contig list canonically: descending by length,
  then ascending lexicographically.
- **Consequence:** `assemble` is strictly permutation-invariant: any
  permutation of the same input list produces identical contig ordering.

### ADR-5: Total Error Handling and Strong Newtypes

- **Context:** The implementation is targetted to Haskell. Partial functions
  and stringly-typed representations introduce runtime panics and domain
  confusion.
- **Decision:** Use strict `Data.Text` wrapped in dedicated `Read` and
  `Contig` newtypes. Return `Either AssemblyError [Contig]` to handle
  invalid thresholds or empty reads as total values.
- **Consequence:** Eliminates runtime exceptions, enforces semantic
  distinctions between raw reads and assembled contigs, and guarantees
  referential transparency.
