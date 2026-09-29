# Reconstructing DNA from Short Fragments

Reconstructing a complete DNA strand from short fragments—a process known as
[*de novo* sequence assembly][de-novo-sequence-assembly]—amounts to solving a
massive jigsaw puzzle where the reference picture is missing and pieces overlap
with one another.

Modern sequencers cannot read a chromosome from end to end. Instead, they shred
millions of copies of the genome into short fragments called "reads" (ranging
from 150 base pairs in short-fragment sequencing to tens of thousands of base
pairs in long-fragment technologies).

The computational challenge is to identify where these fragments share identical
subsequences and stitch them together into contiguous sequences called
"contigs".

Two fundamental algorithmic paradigms govern genome assembly:

1. **Overlap-Layout-Consensus (OLC) and String Graphs:** Compares fragments
   directly to identify suffix-prefix overlaps, constructs an overlap graph, and
   threads a path through every fragment. Originally developed for Sanger
   sequencing, OLC has seen a major resurgence as the primary paradigm for
   modern long-read sequencing (e.g. PacBio HiFi and Oxford Nanopore) where
   fragments are long enough to span complex genomic repeats.
2. **[de Bruijn Graphs (DBG)][de-bruijn-graph]**: Breaks fragments into smaller
   substrings of fixed length $k$ ($k$-mers) and models assembly as an
   [Eulerian path][eulerian-path] problem. DBG is the primary workhorse for
   high-throughput short-read technologies (e.g. Illumina), where billions of
   short fragments make pairwise comparisons computationally intractable.

## 1. The Overlap Approach

At its simplest, sequence assembly is theoretically modelled as the
[Shortest Common Superstring (SCS)][scs-link] problem: given a set of strings,
find the shortest string that contains each input string as a substring.

While SCS provides an intuitive framework, finding the exact optimal superstring
is NP-hard. Pure SCS is also biologically simplified because it collapses
identical genomic repeats into a single copy. In practice, assemblers use
graph-based overlap representations that preserve repeat structure and coverage
information.

In an **overlap graph**, each fragment is represented as a vertex (node), and a
directed edge is drawn from fragment $A$ to fragment $B$ if a suffix of $A$
matches a prefix of $B$. Reconstructing the original genome by visiting each
fragment exactly once corresponds to finding a **Hamiltonian path**, which is an
NP-complete problem. A **greedy overlap heuristic** offers an effective
approximation on clean, toy-scale data.

### How the Greedy Overlap Algorithm Works

1. **Calculate Overlaps:** For every ordered pair of strings $(A, B)$, determine
   the length of the longest suffix of string $A$ that matches a prefix of
   string $B$. A minimum overlap threshold (e.g. at least 2 or 3 characters) is
   enforced to avoid coincidental matches.
2. **Select the Maximum Overlap:** Find the pair $(A, B)$ with the longest
   overlap across the entire pool.
3. **Merge the Pair:** Merge $A$ and $B$ into a composite string by appending
   the non-overlapping suffix of $B$ onto $A$.
4. **Repeat:** Replace $A$ and $B$ in the pool with the merged string. Repeat
   steps 1 to 3 until only one contig remains, or no remaining pairs overlap
   above the minimum threshold.

### Walkthrough with Toy Strings

Suppose we have three fragments from an unknown sequence:

- Fragment 1: `ATGGC`
- Fragment 2: `GGCGT`
- Fragment 3: `CGTGCA`

Assume a minimum overlap threshold of `min_overlap = 2`.

```text
Step 1: Check pairwise suffix-to-prefix overlaps
  ATGGC  -> GGCGT  : Suffix 'GGC' matches prefix 'GGC' (overlap = 3)
  GGCGT  -> CGTGCA : Suffix 'CGT' matches prefix 'CGT' (overlap = 3)
  ATGGC  -> CGTGCA : Suffix 'C' matches prefix 'C' (overlap = 1, < threshold)

Step 2 & 3: Merge highest overlap (Fragment 1 and Fragment 2 on overlap = 3)
  ATGGC
    GGCGT
  --------
  ATGGCGT  (New pool: ['ATGGCGT', 'CGTGCA'])

Step 4: Repeat for remaining strings (overlap = 3)
  ATGGCGT -> CGTGCA : Suffix 'CGT' matches prefix 'CGT' (overlap = 3)
  ATGGCGT
      CGTGCA
  ----------
  ATGGCGTGCA (Single consensus contig reconstructed)
```

## 2. The de Bruijn Graph Approach

The greedy overlap approach requires comparing fragments against one another,
scaling as $O(N^2 \cdot L)$ where $N$ is the number of fragments and $L$ is
fragment length. When handling hundreds of millions of short fragments, building
an explicit overlap graph becomes computationally intractable.

Modern short-read assemblers therefore use de Bruijn graphs instead of a direct
pairwise fragment-overlap model.

### How de Bruijn Graphs Work

1. **Deconstruct Fragments into $k$-mers:** Choose a fixed length $k$. Slide a
   window of size $k$ across every fragment to extract all constituent $k$-mers.
2. **Define Nodes and Edges:**
   - **Nodes:** Every distinct $(k - 1)$-mer (the prefix or suffix of a
     $k$-mer).
   - **Directed Edges:** Every $k$-mer forms a directed edge pointing from its
     prefix $(k - 1)$-mer to its suffix $(k - 1)$-mer.
3. **Find an Eulerian Path:** Rather than visiting every node (the NP-complete
   Hamiltonian path problem on overlap graphs), de Bruijn graphs model assembly
   as finding a path that visits every *edge* (an
   [Eulerian path][eulerian-path]). An Eulerian path can be solved in linear
   time, $O(V + E)$, using [Hierholzer's algorithm][hierholzers-algorithm].

**Note:** In real sequencing datasets with high coverage depth (e.g. 30x–100x),
$k$-mers occur multiple times. Real assemblers construct Eulerian multigraphs
where edge multiplicities and coverage weights guide the traversal, rather than
a simple single-visit path.

### Walkthrough with the Same Toy Strings

Using the exact same input fragments as Section 1 (`ATGGC`, `GGCGT`, `CGTGCA`)
with $k = 4$:

Extract all 4-mers from each fragment:

- From `ATGGC`: `ATGG`, `TGGC`
- From `GGCGT`: `GGCG`, `GCGT`
- From `CGTGCA`: `CGTG`, `GTGC`, `TGCA`

Combine the 4-mers to form the directed edges of the graph:

- `ATGG`: `ATG` -> `TGG`
- `TGGC`: `TGG` -> `GGC`
- `GGCG`: `GGC` -> `GCG`
- `GCGT`: `GCG` -> `CGT`
- `CGTG`: `CGT` -> `GTG`
- `GTGC`: `GTG` -> `TGC`
- `TGCA`: `TGC` -> `GCA`

Graph nodes are the distinct 3-mers: `ATG`, `TGG`, `GGC`, `GCG`, `CGT`, `GTG`,
`TGC`, `GCA`.

Following the directed edges from start to end yields an Eulerian path:

```text
ATG -> TGG -> GGC -> GCG -> CGT -> GTG -> TGC -> GCA
```

Reconstructed consensus sequence: `ATGGCGTGCA`.

Both the greedy overlap method and the de Bruijn graph approach yield the
identical reconstructed sequence.

## 3. Real-World Biological Complications

When transitioning from toy character strings to biological sequencing data,
assemblers must resolve multiple physical complexities:

- **Genomic Repeats:** Genomes contain extensive repetitive elements
  (transposons, segmental duplications) that are longer than individual
  fragments or $k$-mers. In a graph, repeats create tangled cycles and branching
  junctions, rendering reconstruction ambiguous.
- **Sequencing Errors:** A single misread nucleotide creates spurious $k$-mers,
  producing false dead-end branches ("tips") or alternative bubbles in the
  graph. Assemblers employ error correction and tip-clipping algorithms to prune
  these artefacts.
- **Reverse Complements:** DNA is double-stranded. Fragments originate randomly
  from either the forward or reverse strand. An assembler must match fragments
  against both forward sequences and their reverse complements (e.g. `AAGCT` and
  its reverse complement `AGCTT`).
- **Uneven Coverage:** Sequencing depth varies across the genome due to PCR bias
  and GC content. Regions with low coverage cause breaks in the graph, resulting
  in fragmented contigs rather than whole chromosomes.
- **Paradigm Selection (Short vs Long Reads):**
  - **Short Reads (Illumina):** Highly accurate (< 0.1% error) but short
    (150–300 bp). Ideal for de Bruijn graphs, but cannot span long repeats.
  - **Long Reads (PacBio HiFi, Nanopore):** Tens of thousands of base pairs in
    length, capable of spanning complex repeats. OLC and string graphs are the
    dominant paradigm because breaking long fragments into short $k$-mers would
    destroy their long-range structural information.

## 4. See Also

- [Sequence Assembly (Wikipedia)][sequence-assembly]
- [de Bruijn Graph (Wikipedia)][de-bruijn-graph]
- [Eulerian Path (Wikipedia)][eulerian-path]
- [Shortest Common Superstring (Wikipedia)][scs-link]

[de-bruijn-graph]: https://en.wikipedia.org/wiki/De_Bruijn_graph
[de-novo-sequence-assembly]: https://en.wikipedia.org/wiki/De_novo_sequence_assemblers
[eulerian-path]: https://en.wikipedia.org/wiki/Eulerian_path
[hierholzers-algorithm]: https://en.wikipedia.org/wiki/Hierholzer%27s_algorithm
[scs-link]: https://en.wikipedia.org/wiki/Shortest_common_supersequence#Shortest_common_superstring
[sequence-assembly]: https://en.wikipedia.org/wiki/Sequence_assembly
