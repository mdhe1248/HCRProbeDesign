# HCRProbeDesign GUI

A local Streamlit front-end for the `HCRProbeDesign` package. Runs on `localhost`
for a single user — no auth, no deployment.

## Setup (one time)

Clone the repo and run everything below from inside it:

```bash
git clone https://github.com/mdhe1248/HCRProbeDesign.git
cd HCRProbeDesign
```

The package's own `environment.yaml` is outdated (Python 3.7 / macOS bowtie2 build
from 2021) and won't work here. Create a fresh conda environment instead:

```bash
conda create -n hcrprobedesign python=3.11 -y
conda install -n hcrprobedesign -c bioconda -c conda-forge bowtie2 -y
conda run -n hcrprobedesign pip install -e .
conda run -n hcrprobedesign pip install -r gui/requirements.txt
```

Verify:

```bash
conda run -n hcrprobedesign python -c "import HCRProbeDesign; print('ok')"
conda run -n hcrprobedesign which bowtie2 bowtie2-build
```

## Run

From inside the cloned repo directory:

```bash
conda activate hcrprobedesign
streamlit run gui/app.py
```

Opens at `http://localhost:8501`.

## Usage

1. Provide a target sequence, either:
   - **Paste or upload FASTA** (single record). Add `channel=B2` (etc.) to the
     header to override the sidebar channel for that record.
   - **Search by gene name**: enter a species and gene symbol, pick the matching
     gene, then pick a transcript variant from NCBI RefSeq — the sequence is
     fetched automatically and dropped into the FASTA box above. Requires an
     NCBI contact email (prefilled, editable) and internet access.
2. Adjust parameters in the sidebar (channel, tile size, GC/Gibbs/dTm ranges,
   homopolymer limits, max probes, genome off-target masking).
3. Click **Design probes**. Results appear as a table with download buttons for
   the probe TSV and an IDT ordering sheet.
4. To design against a species other than what's already registered, use the
   **Register a new reference genome** panel at the bottom. Three ways to add one
   (all a one-time step per species, not part of routine probe design):
   - **Common genome**: one-click download of a prebuilt Bowtie2 index (mouse,
     human, zebrafish, fly, *C. elegans*, rat) — no genome search, no
     `bowtie2-build`, fastest option.
   - **Genome URL**: paste a direct FASTA link (e.g. Ensembl or UCSC) — the
     server downloads it directly (no browser upload of a multi-GB file), then
     builds the index locally. Use this for any species not in the common list.
   - **Upload file(s)**: manual browser upload, for small/custom sequences only
     (a transgene, a plasmid) — real genomes are far larger than Streamlit's
     200MB upload cap.

## Not in scope for v1

- Multi-record / batch FASTA upload (the underlying package supports batch mode
  via `designProbesBatch`; the GUI only handles one record at a time for now).
- Multi-user auth, remote/networked deployment, HTTPS.
- Editing `HCRconfig.yaml` default parameters from the UI (the GUI only reads
  them to pre-fill widget defaults).
- A working repeat-mask toggle — this flag is permanently disabled upstream
  (an argparse default/action bug in the underlying package) and isn't exposed
  here.
