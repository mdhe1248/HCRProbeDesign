# HCRProbeDesign

## Overview
HCRProbeDesign is a command-line and Python package for designing HCR v3.0 split-initiator probe pairs
from a target FASTA sequence. It tiles target sequences, filters by GC content, melting temperature,
homopolymers, and hairpins, and can optionally enforce genome uniqueness with Bowtie2.

Key tools:
- `designProbes`: primary probe design CLI.
- `designProbesBatch`: batch probe design for multi-record FASTA inputs.
- `fetchMouseIndex`: download, install, and register the mm10 Bowtie2 index.
- `buildGenomeIndex`: build and register a new reference genome index.
- `listReferences`: display installed reference genomes and default parameters.

Documentation: https://www.gofflab.org/HCRProbeDesign/

## GUI

This fork adds a local Streamlit GUI as an alternative to the CLI tools above — paste, upload, or
search NCBI for a target sequence by species/gene name; tune parameters with built-in explanations
of what each one means and how it affects signal/specificity; and manage reference genomes
(one-click prebuilt index, direct URL, or manual upload; plus rename/delete) without touching the
command line.

Quick start (once the `hcrprobedesign` conda environment is set up):
```bash
conda activate hcrprobedesign
streamlit run gui/app.py
```

See [`gui/README.md`](gui/README.md) for full setup and usage instructions.

## CLI tools

This fork's `designProbes`, `designProbesBatch`, `buildGenomeIndex`, `fetchMouseIndex`, and
`listReferences` CLI tools are unchanged from upstream. For CLI installation and usage, see the
original repository: https://github.com/gofflab/HCRProbeDesign
