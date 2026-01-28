# Chloroplast splicing coverage + quantification pipeline

This repo contains the scripts used to (i) convert RNA-seq BAM alignments against a concatenated chloroplast reference into per-base coverage tables, and (ii) quantify intron splicing efficiency from spliced vs unspliced junction signals.

---

## Overview

1. **Align RNA-seq reads** to a concatenated chloroplast reference (**catref.fa**) containing:
   - an **unspliced** chloroplast genome reference, and
   - a **spliced** chloroplast genome reference.

2. **Convert BAM → per-base coverage** using:
   - `spliced_coverage.jl` (CLI wrapper)
   - `countreads.jl` (coverage counting functions used by `spliced_coverage.jl`)
   - run via `splicing.sh`

   This produces count tables:
   - `At1_splicing.counts … At7_splicing.counts`

   Each `.counts` file contains read coverage for **every position** across both the unspliced and spliced reference sequences in `catref.fa`.

3. **Map intron annotations**
   - Intron annotations were taken from `AP000423.sff` and mapped onto the corresponding coordinates in the unspliced and spliced references.

4. **Quantify splicing efficiency** from the `.counts` tables:
   - **Unspliced signal:** read counts at **exon–intron boundary** positions in the unspliced reference
   - **Spliced signal:** read counts at **exon–exon junction** positions in the spliced reference
   - Splicing proportion (per intron, per sample):

     `spliced / (spliced + unspliced)`

---

## Files

### Core scripts
- `countreads.jl`  
  Implements read/pair filtering and coverage counting (used by `spliced_coverage.jl`).

- `spliced_coverage.jl`  
  Command-line tool that opens a BAM and writes a per-base forward/reverse coverage table.

### Runner
- `splicing.sh`  
  Convenience script used to run `spliced_coverage.jl` across samples.

---

## Requirements

- Julia
- samtools (for producing **name-sorted** BAMs)

### Julia packages
- ArgParse
- XAM