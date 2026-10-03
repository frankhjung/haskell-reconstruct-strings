# Greedy Overlap Sequence Assembler Animation Guide

[← Back to Documentation Index](README.md)

This document describes the visual representation of the greedy overlap
sequence assembler algorithm. It is intended to guide an image generation
model in creating a step-by-step animation of the process.

For a formal sequence diagram of the same pipeline, see the
[Assembler Specification][req-001].

## Scene Setup

The scene consists of a pool of DNA sequence fragments or text strings
floating in a dynamic workspace. Each fragment is represented as a distinct,
colourful block with its sequence clearly printed on it. The background
should be a dark, sleek workspace with glowing nodes.

## Animation Steps

### 1. Initialisation and Filtering

- **Action:** The initial pool of raw fragments is displayed.
- **Visual:** Fragments that are exact duplicates or fully contained within
  longer fragments begin to fade out and disappear. The remaining fragments
  snap into an organised grid, representing the unique, uncontained initial
  pool.

### 2. Searching for Overlaps

- **Action:** The algorithm searches for the best suffix-prefix overlap.
- **Visual:** Glowing lines connect pairs of fragments, testing the end of
  one (prefix) against the start of another (suffix). As matches are found,
  the overlapping sections light up.
- **Highlight:** The longest overlap shines the brightest. If there is a tie
  in length, the algorithm resolves it by choosing the pair with the
  lexicographically smaller prefix, then suffix. The winning pair is
  highlighted with a bold, pulsating border.

### 3. Merging the Best Pair

- **Action:** The chosen prefix and suffix fragments merge.
- **Visual:** The two fragments are pulled together magnetically. The
  overlapping letters lock in place, while the redundant letters dissolve.
  They fuse into a single, longer contiguous block (a contig).

### 4. Pool Reduction and Re-filtering

- **Action:** The newly merged fragment is added back to the pool, replacing
  its two source fragments.
- **Visual:** The workspace rearranges. Any existing fragments that are now
  fully contained within this new, larger fragment fade away.

### 5. Repetition

- **Action:** Steps 2–4 repeat recursively.
- **Visual:** The process speeds up, showing continuous matching, fusing, and
  filtering until no more valid overlaps exist that meet the minimum overlap
  threshold. The pool gradually reduces to a small number of long contigs.

### 6. Canonical Sorting (Final Output)

- **Action:** The final set of assembled contigs is sorted.
- **Visual:** The remaining large blocks align in a vertical list. They are
  ordered first by length (longest at the top), and then alphabetically for
  any ties. The scene concludes with the final assembled sequences glowing
  triumphantly against the dark background.

[req-001]: REQ-001-functional-greedy-overlap-assembler.md
