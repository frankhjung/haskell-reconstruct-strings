# Reconstruct Strings

A pure, total Haskell implementation of a greedy overlap sequence assembler
designed to reconstruct contiguous sequences (*contigs*) from short overlapping
string fragments or DNA sequencing reads.

## Overview

Reconstructing a DNA sequence from short read fragments—a foundational step in
*de novo* genome assembly—can be approximated through the Overlap-Layout-
Consensus (OLC) paradigm.

This project implements a deterministic, pure functional greedy overlap
reduction algorithm that:

* Validates input parameters totally without partial runtime exceptions.
* Deduplicates identical reads to a single representative sequence.
* Eliminates reads contained as proper substrings within longer reads.
* Merges the candidate pair with the longest valid suffix-prefix match.
* Resolves ties deterministically using a strict three-tier total order.
* Dynamically purges reads engulfed by newly formed composite contigs.
* Terminates when no pairwise overlap satisfies the minimum threshold.
* Sorts final contigs into a canonical, permutation-invariant order.

## Project Structure

* [`reconstruct-strings.cabal`][cabal-file]:
  Package configuration declaring build dependencies, compiler flags, and
  components.
* [`Makefile`][makefile]:
  Development automation targets for formatting, linting, building, testing,
  and documentation generation.
* [`src/Assembler.hs`][src-assembler]:
  Core greedy reduction logic and containment filtering.
* [`src/Assembler/Types.hs`][src-assembler-types]:
  Domain newtypes (`Read`, `Contig`), candidate records, and error types.
* [`app/Main.hs`][app-main]:
  Command-line interface with option parsing.
* [`test/Spec.hs`][test-spec]:
  Test driver with `hspec-discover`.
* [`test/AssemblerSpec.hs`][test-assembler-spec]:
  Hspec and QuickCheck test suite.
* [`docs/FunctionalGreedyOverlapAssemblerSpecification.md`][spec-doc]:
  Formal functional specification and architectural decision records (ADRs).
* [`docs/GLOSSARY.md`][glossary]:
  Domain terminology and definitions.

## Build and Development

The project uses [Cabal][cabal-url] and a `Makefile` task
runner.

### Common Make Targets

* `make` (or `make default`): Formats, lints, builds, and runs the test suite.
* `make all`: Runs the full pipeline including documentation and sample
  execution.
* `make format`: Formats source files with `stylish-haskell` and Cabal files
  with `cabal-fmt`.
* `make check`: Generates ctags and runs static analysis with `hlint` and
  `cabal check`.
* `make build`: Compiles the library, executable, and test suite with Cabal.
* `make test`: Executes unit tests and property tests.
* `make doc`: Generates Haddock documentation with hyperlinked source code.
* `make exec`: Executes sample sequence assembly runs.
* `make clean`: Cleans build artifacts and tags.

## Usage

### Command-Line Execution

Run the assembler directly using Cabal:

```bash
cabal run reconstruct-strings -- -m 2 ATGGC GGCGT CGTGCA
```

Output:

```text
ATGGCGTGCA
```

### Options

* `-m`, `--min-overlap INT`: Minimum overlap threshold (default: `2`).
* `-f`, `--file FILE`: Read fragments from a file (one read per line).
* `READ...`: Read fragments passed as positional arguments or piped via `stdin`.

## Continuous Integration

GitHub Actions pipeline configuration is defined in
[`.github/workflows/haskell.yml`][github-actions], which
automatically compiles, tests, and validates documentation builds on pushes
and pull requests.

## References

* [Functional Greedy Overlap Assembler Specification][spec-doc]
* [Reconstructing DNA from Short Fragments][dna-doc]
* [Domain Glossary][glossary]

[app-main]: app/Main.hs
[cabal-file]: reconstruct-strings.cabal
[cabal-url]: https://www.haskell.org/cabal/
[dna-doc]: docs/ReconstructingCompleteDNAStrandFromShortFragments.md
[github-actions]: .github/workflows/haskell.yml
[glossary]: docs/GLOSSARY.md
[makefile]: Makefile
[spec-doc]: docs/FunctionalGreedyOverlapAssemblerSpecification.md
[src-assembler]: src/Assembler.hs
[src-assembler-types]: src/Assembler/Types.hs
[test-assembler-spec]: test/AssemblerSpec.hs
[test-spec]: test/Spec.hs
