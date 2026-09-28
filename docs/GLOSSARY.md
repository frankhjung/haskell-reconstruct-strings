# Glossary

## AssemblyError

An explicit failure condition returned when input preconditions are violated,
such as specifying a minimum overlap threshold less than one.

## Candidate

An ordered pair of distinct reads that shares a valid suffix-prefix match
exceeding or meeting the minimum overlap threshold. Not the same as an arbitrary
pair of reads.

_Avoid_: Match

## Containment

The condition where a read is a proper substring of another sequence, rendering
the shorter read redundant as it contributes no novel information.

## Contig

A contiguous sequence produced by iteratively merging overlapping reads, or an
isolated singleton read that could not be merged. Distinct from a raw Read.

_Avoid_: Scaffold

## Overlap

An exact match where a suffix of a prefix read is identical to a prefix of a
suffix read, strictly shorter than the longer read and meeting the threshold.

_Avoid_: Alignment

## Read

A finite, immutable sequence of characters representing a single fragment of
sequenced genetic material, distinct from an assembled Contig.

_Avoid_: Fragment, K-mer
