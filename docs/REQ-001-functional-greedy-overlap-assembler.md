# Functional Greedy Overlap Assembler Specification

## 1. Overview and Problem Statement

The assembler accepts an unordered collection of fragments and constructs one or
more contigs by repeatedly merging fragments whose suffixes overlap prefixes.
This is a greedy approximation of *de novo* sequence assembly designed to be
deterministic, pure, total, and mathematically rigorous.

The implementation must not depend on mutable state, imperative loops, or
in-place list modification. The algorithm is expressed as a pure functional
reduction over immutable values in Haskell.

The system produces a deterministic list of contigs from an input list of
fragments. A correct assembly is one that:

- validates input parameters totally without runtime exceptions,
- deduplicates identical fragments to a single representative,
- removes fragments that contribute no new sequence information (containment),
- greedily merges the candidate pair with the best valid overlap at each step,
- resolves ties deterministically using a strict total order,
- eliminates fragments that become contained within newly merged contigs,
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
  ( Fragment(..)
  , Contig(..)
  , OverlapLength
  , OverlapCandidate(..)
  , AssemblyError(..)
  ) where

import Data.Text (Text)

newtype Fragment = Fragment { unFragment :: Text }
  deriving stock (Eq, Ord, Show)

newtype Contig = Contig { unContig :: Text }
  deriving stock (Eq, Ord, Show)

type OverlapLength = Int

data OverlapCandidate = OverlapCandidate
  { prefixFragment :: !Fragment
  , suffixFragment :: !Fragment
  , matchLength    :: !OverlapLength
  }
  deriving stock (Eq, Show)

instance Ord OverlapCandidate where
  compare a b =
    compare (matchLength a) (matchLength b)
      <> compare (prefixFragment b) (prefixFragment a)
      <> compare (suffixFragment b) (suffixFragment a)

data AssemblyError
  = InvalidMinOverlap !Int
  | EmptyFragmentEncountered
  deriving stock (Eq, Show)
```

### Required Interpretation

- `Fragment`: An immutable sequence of characters. Must be non-empty.
- `Contig`: A contiguous sequence produced by merging overlapping fragments, or
  an unmerged singleton fragment.
- `OverlapLength`: A non-negative integer representing the length of an exact
  suffix-prefix match.
- A valid overlap is an exact match where a suffix of `prefixFragment` matches a
  prefix of `suffixFragment`. The overlap length must be strictly greater than
  zero, at least `min_overlap`, and strictly less than the length of the longer
  fragment.

## 4. Functional Requirements

### 4.1 Required Input Assumptions

The entry point `assemble` accepts:

- an input list of `Fragment` values: `[Fragment]`,
- a minimum overlap threshold: `min_overlap :: Int`.

#### Validation Rules

- If `min_overlap < 1`, the function must return
  `Left (InvalidMinOverlap min_overlap)`.
- If any fragment in the input has length 0, the function must return
  `Left EmptyFragmentEncountered`.
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
calculateOverlap :: Fragment -> Fragment -> Int -> OverlapLength
calculateOverlap (Fragment r1) (Fragment r2) minOverlap =
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

`mergePair prefix suffix overlapLen` concatenates the prefix fragment with the
non-overlapping suffix remainder:

```haskell
mergePair :: Fragment -> Fragment -> OverlapLength -> Fragment
mergePair (Fragment p) (Fragment s) overlapLen =
  Fragment (p <> T.drop overlapLen s)
```

No other characters may be inserted, dropped, or modified.

### 4.4 `findBestOverlap`

`findBestOverlap pool minOverlap` inspects every ordered pair of distinct
fragments `(a, b)` from `pool` where `a /= b` and computes
`calculateOverlap a b minOverlap`.

#### Output Contract

- If no valid pair satisfies `overlap >= minOverlap`, return `Nothing`.
- Otherwise, return `Just candidate` where `candidate` is the maximal candidate
  under the deterministic tie-breaking order.

#### Deterministic Tie-Breaking Order

Candidate pairs are ordered strictly by the following total ordering:

1. Greater `matchLength` takes precedence.
2. If `matchLength` is equal, lexicographically smaller `prefixFragment` takes
   precedence.
3. If `prefixFragment` is equal, lexicographically smaller `suffixFragment`
   takes precedence.

Because the candidate pool contains pairwise distinct fragments, Rules 1 to 3
guarantee a strict total order over all candidate pairs.

### 4.5 `filterContainedFragments`

A fragment `target` is contained within another fragment `other` if:

- `target /= other`, and
- `unFragment target` is an infix (proper substring) of `unFragment other`.

Pre-processing removes exact duplicates (retaining one unique representative)
and drops any fragment that is a proper substring of another fragment in the
pool:

```haskell
filterContainedFragments :: [Fragment] -> [Fragment]
filterContainedFragments fragments =
  let uniqueFragments = nub fragments
  in [ r | r <- uniqueReads
         , not (any (\o -> r /= o &&
                           unFragment r `T.isInfixOf` unFragment o) uniqueReads) ]
```

## 5. Required Assembly Behaviour

### 5.1 High-Level Contract

The primary assembly function has the following signature:

```haskell
assemble :: [Fragment] -> Int -> Either AssemblyError [Contig]
```

It must guarantee:

- input validation without partial functions,
- initial duplicate deduplication and containment elimination,
- iterative greedy merging of the maximal candidate overlap,
- dynamic elimination of fragments engulfed by newly formed contigs,
- and canonical, permutation-invariant ordering of the final contigs.

### 5.2 Reference Reduction Algorithm

```haskell
assemble :: [Fragment] -> Int -> Either AssemblyError [Contig]
assemble fragments minOverlap
  | minOverlap < 1 = Left (InvalidMinOverlap minOverlap)
  | any (T.null . unFragment) fragments = Left EmptyFragmentEncountered
  | otherwise =
      let initialPool = filterContainedFragments fragments
          finalPool = reducePool initialPool minOverlap
          sortedContigs = sortCanonical (map toContig finalPool)
      in Right sortedContigs

reducePool :: [Fragment] -> Int -> [Fragment]
reducePool pool minOverlap
  | length pool <= 1 = pool
  | otherwise =
      case findBestOverlap pool minOverlap of
        Nothing -> pool
        Just candidate ->
          let merged = mergePair (prefixFragment candidate)
                                 (suffixFragment candidate)
                                 (matchLength candidate)
              remaining = [ r | r <- pool
                              , r /= prefixFragment candidate
                              , r /= suffixFragment candidate ]
              updatedPool = filterContainedFragments (merged : remaining)
          in reducePool updatedPool minOverlap
```

### 5.3 Dynamic Containment and Invariants

- The pool is strictly immutable; every step returns a newly allocated list.
- Each merge step decreases the pool size by at least 1 (more if the merge
  engulfs another fragment).
- If no valid overlap remains, reduction terminates.
- `updatedPool` re-evaluates containment to eliminate any fragment engulfed
  within the newly merged contig.

### 5.4 Canonical Output Ordering

Because input fragments constitute an unordered multiset, the output contig list
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

If `fragments` is empty, `assemble fragments minOverlap` returns `Right []`.

### 6.2 Singleton Input

If `fragments` contains exactly one non-empty fragment `[r]`, `assemble` returns
`Right [Contig (unFragment r)]`.

### 6.3 All Fragments Identical

If all input fragments are identical (e.g. `["ACGT", "ACGT"]`), duplicates
represent redundant sequencing coverage. Pre-processing deduplicates them to a
single representative, returning `Right [Contig "ACGT"]`.

### 6.4 Contained Fragments

Any fragment that is a proper substring of another fragment is discarded. This
applies both in pre-processing and dynamically during reduction when two
fragments merge to span a third fragment.

### 6.5 No Valid Overlap

If no pair meets `min_overlap`, greedy reduction terminates immediately. The
remaining fragments are wrapped as contigs and sorted in canonical order.

### 6.6 Disjoint Components

If fragments originate from disconnected genomic regions, components reduce
independently. All resulting contigs are sorted canonically in the final output.

## 7. Computational Complexity and Invariants

### 7.1 Termination

The pool size strictly decreases by at least 1 per recursive iteration. For an
initial filtered pool of `N` fragments, reduction executes at most `N - 1`
merges.

### 7.2 Time Complexity

For `N` fragments of average length `L`:

- Pairwise overlap evaluation takes $O(N^2 \cdot L)$$.
- At most $N - 1$ reduction stages occur.
- Worst-case runtime of naive reduction is $O(N^3 \cdot L)$.

### 7.3 Space Complexity

Because Haskell is non-strict by default, strict evaluation of contig text and
accumulator lists prevents space leaks. The memory footprint is $O(N \cdot L)$.

### 7.4 Invariants

Throughout all stages, the implementation guarantees:

- purity and referential transparency with total functions,
- absence of contained fragments in active pools,
- strict total ordering for overlap candidates,
- and permutation invariance of the final contig list.

## 8. Acceptance Criteria

The implementation is verified if:

1. `assemble` accepts `[Fragment]` and `Int`, returning
   `Either AssemblyError [Contig]`.
2. Invalid overlap thresholds (`< 1`) or empty fragments yield descriptive
   errors.
3. Identical fragments collapse to a single representative contig.
4. Fragments contained within other fragments or newly merged contigs are
   eliminated.
5. Ties resolve strictly by `matchLength`, `prefixFragment`, then
   `suffixFragment`.
6. Final contig output is sorted canonically and invariant to input permutation.
7. Pure functions and immutable data structures are used exclusively.
8. Property-based tests verify determinism, totality, and containment
   invariants.

## 9. Example Test Cases

### Example A: Single Valid Overlap

```haskell
fragments = [Fragment "ABC", Fragment "BCD", Fragment "CDE"]
minOverlap = 2
-- Expected: Right [Contig "ABCDE"]
```

### Example B: No Valid Overlap

```haskell
fragments = [Fragment "ABC", Fragment "DEF"]
minOverlap = 2
-- Expected: Right [Contig "ABC", Contig "DEF"]
```

### Example C: Contained Fragment Removed

```haskell
fragments = [Fragment "ACGT", Fragment "CGT"]
minOverlap = 2
-- Expected: Right [Contig "ACGT"]
```

### Example D: Dynamic Containment After Merge

```haskell
fragments = [Fragment "AAATTT", Fragment "TTTGGG", Fragment "ATT"]
minOverlap = 3
-- Expected: Right [Contig "AAATTTGGG"]
-- "ATT" is engulfed by the merged sequence and eliminated.
```

### Example E: Canonical Permutation Invariance

```haskell
assemble [Fragment "DEF", Fragment "ABC"] 2 == assemble [Fragment "ABC", Fragment "DEF"] 2
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

- **Context:** An input may consist solely of identical fragments (e.g.
  `["ACGT", "ACGT"]`). An earlier ambiguity suggested either discarding all
  duplicates or preserving a single copy.
- **Decision:** Deduplicate identical fragments to a single representative.
  Multiple identical fragments represent repeated sequencing coverage of the
  same region, not invalid data.
- **Consequence:** Output for an all-duplicate input is a single contig
  (`Right [Contig "ACGT"]`). Information is preserved without generating false
  self-overlaps.

### ADR-2: Strict Three-Tier Deterministic Tie-Breaking

- **Context:** An earlier specification included a fourth tie-breaking rule
  referencing input list indices. However, because fragments in the candidate
  pool are distinct strings, rules 1 to 3 form a complete strict total order.
  Moreover, merged contigs lack input indices.
- **Decision:** Remove the fourth rule. Candidate pairs are ordered solely by:
  1. Greater `matchLength`
  2. Lexicographically smaller `prefixFragment`
  3. Lexicographically smaller `suffixFragment`
- **Consequence:** Eliminates unnecessary state tracking and index provenance
  across recursive reductions. Purity and totality are preserved.

### ADR-3: Dynamic Containment Elimination

- **Context:** Merging two overlapping fragments can create a composite sequence
  that engulfs a third, previously uncontained fragment in the pool.
- **Decision:** Re-evaluate containment after each merge (`updatedPool`) rather
  than only once during pre-processing.
- **Consequence:** Contained fragments are promptly purged, preventing redundant
  sequence emissions and degenerate zero-gain merge cycles.

### ADR-4: Canonical Permutation Invariance

- **Context:** The input model is an unordered multiset of fragments. However,
  naive recursive reduction prepends merged fragments, which could produce
  different contig orders for permuted inputs.
- **Decision:** Sort the final contig list canonically: descending by length,
  then ascending lexicographically.
- **Consequence:** `assemble` is strictly permutation-invariant: any permutation
  of the same input list produces identical contig ordering.

### ADR-5: Total Error Handling and Strong Newtypes

- **Context:** The implementation is targetted to Haskell. Partial functions and
  stringly-typed representations introduce runtime panics and domain confusion.
- **Decision:** Use strict `Data.Text` wrapped in dedicated `Fragment` and
  `Contig` newtypes. Return `Either AssemblyError [Contig]` to handle invalid
  thresholds or empty fragments as total values.
- **Consequence:** Eliminates runtime exceptions, enforces semantic distinctions
  between raw fragments and assembled contigs, and guarantees referential
  transparency.
