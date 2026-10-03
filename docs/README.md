# Documentation Directory

[← Back to Project README](../README.md)

This directory contains specifications, theoretical background, domain guides,
and visual storyboards for the **reconstruct-strings** project.

## Document Index

### Specifications and Architecture

- [Domain Glossary][glossary]: Domain terminology and definitions covering
  fragments, contigs, overlaps, and sequencing coverage.
- [Functional Greedy Overlap Assembler Specification (REQ-001)][req-001]:
  Formal requirements, type contracts, Mermaid sequence diagram, reduction
  algorithm, and architectural decision records (ADRs) for the assembler.
- [Synthetic Strand Generator Specification (REQ-002)][req-002]: Functional
  requirements, command-line interface, and ADRs for the synthetic strand
  extraction utility.

### Domain Theory and Parameter Analysis

- [Reconstructing DNA from Short Fragments][dna-doc]: Biological background
  on *de novo* genome assembly, contrasting Overlap-Layout-Consensus (OLC)
  with de Bruijn graph (DBG) paradigms.
- [Assembly Dynamics and Parameter Heuristics][heuristics]: Mathematical
  models for collision probabilities, minimum overlap lower/upper bounds,
  Lander–Waterman coverage depth, and calibrated parameter examples.

### Visualisation and Media

- [Assembler Animation Guide][animation]: Visual storyboard detailing
  step-by-step animation of pool reduction, tie-breaking, and canonical
  contig sorting.

[animation]: animation.md
[dna-doc]: reconstructing-complete-dna-strand-from-short-fragments.md
[glossary]: GLOSSARY.md
[heuristics]: heuristics.md
[req-001]: REQ-001-functional-greedy-overlap-assembler.md
[req-002]: REQ-002-shell-script-to-make-strands.md
