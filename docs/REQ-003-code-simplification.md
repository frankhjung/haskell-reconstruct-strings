# Code Simplification Specification (REQ-003)

[← Back to Documentation Index](README.md)

## Objective

Simplify the Haskell implementation and build configuration following a code
review, while retaining every requirement in [REQ-001][req-001]. No observable
behaviour changes: the public API, command-line interface, output and error
messages are identical.

## Requirements

### R1: Total overlap calculation

`calculateOverlap` returns `Maybe OverlapLength` instead of using `0` as a
"no match" sentinel. `findBestOverlap` consumes the result by pattern matching
(`Just len <- [calculateOverlap ...]`), removing the redundant threshold
re-check (`len >= minOverlap`).

### R2: Remove the redundant proper-overlap guard

The `len < max (length prefix) (length suffix)` guard is removed. Callers pass
distinct fragments drawn from a pool that is free of duplicates and containment
(`filterContainedFragments`). A match spanning a whole fragment would require
the fragments to be identical or contained, so it cannot occur. The invariant
is documented on `calculateOverlap`.

### R3: Simplify `reducePool`

The `[]` and `[sole]` clauses are removed. A pool of fewer than two fragments
has no ordered candidate pairs, so `findBestOverlap` returns `Nothing` and the
pool is returned unchanged. The selected pair is removed from the pool with
``filter (`notElem` [prefix, suffix])``.

### R4: Consistent use of `OverlapLength`

`assemble`, `reducePool`, `findBestOverlap`, `calculateOverlap` and
`runPipeline` use the `OverlapLength` alias for the minimum overlap threshold
rather than a bare `Int`. The alias is unchanged (`type OverlapLength = Int`).

### R5: Point-free CLI pipeline

`runPipeline` is expressed as the composition `sanitizeInput`, `map Fragment`,
`assemble`, then `bimap formatAssemblyError (map unContig)`.

### R6: Simplify executable input

`Main` reads input with
`T.lines <$> maybe TIO.getContents TIO.readFile (optInputFile opts)`.
`Options` and `optionsParser` gain Haddock comments.

### R7: Trim build configuration

- Remove unused default extensions: `BangPatterns`,
  `GeneralizedNewtypeDeriving` and `LambdaCase`.
- Remove `-threaded -rtsopts -with-rtsopts=-N` from the executable and test
  suite. The assembler is single-threaded, so the threaded runtime added
  start-up cost without benefit.

## Acceptance Criteria

1. `make` (format, lint, build, test) completes with no compiler warnings and
   no HLint hints.
2. All Hspec examples pass. The QuickCheck properties use
   non-empty fragments, bounded pool sizes and valid thresholds, so they
   exercise the merge path. They now verify shuffled-permutation invariance,
   rejection of empty fragments, absence of contained contigs, preservation of
   every input fragment, and canonical ordering (REQ-001 criterion 8). The
   suite has 22 examples in total.
3. `cabal check` reports no errors or warnings.
4. [REQ-001][req-001] code samples match the implementation.

## Architectural Decisions

### ADR-1: `Maybe` over a sentinel value

- **Context:** `0` was both a valid `OverlapLength` and the "no overlap" marker,
  forcing callers to re-validate against `minOverlap`.
- **Decision:** Return `Maybe OverlapLength`.
- **Consequence:** The absent case is explicit in the type; the caller needs no
  second threshold check.

### ADR-2: Rely on the pool invariant rather than a local guard

- **Context:** The proper-overlap guard re-checked a condition that
  `filterContainedFragments` already guarantees.
- **Decision:** Drop the guard and document the precondition.
- **Consequence:** Less code per candidate pair. `calculateOverlap` is internal;
  `reducePool` only ever passes pools produced by `filterContainedFragments`,
  so the precondition always holds.

### ADR-3: Keep `OverlapLength` as a type alias

- **Context:** A `newtype` would be stricter but changes the public `assemble`
  signature and test suite.
- **Decision:** Retain the alias, applied consistently.
- **Consequence:** No API break. A `newtype` remains a possible future change.

## See Also

- [Documentation Index](README.md)
- [Assembler Specification (REQ-001)][req-001]

[req-001]: REQ-001-functional-greedy-overlap-assembler.md
