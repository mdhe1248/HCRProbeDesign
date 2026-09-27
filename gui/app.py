"""Streamlit GUI for HCRProbeDesign.

Wraps HCRProbeDesign's Python functions directly (no shelling out to the CLI).
Run with: streamlit run gui/app.py
"""
import io
import os
import shutil
import tempfile
import types
import urllib.error

import pandas as pd
import streamlit as st

from HCRProbeDesign import HCR, genomeFetch, geneFetch, listReferences, probeDesign, referenceGenome, sequencelib
from HCRProbeDesign._datadir import ensure_data_dir

st.set_page_config(page_title="HCR Probe Design", layout="wide")
ensure_data_dir()

BOWTIE2_AVAILABLE = shutil.which("bowtie2") is not None
BOWTIE2_BUILD_AVAILABLE = shutil.which("bowtie2-build") is not None

reference_info = listReferences.list_references()
species_names = [s["name"] for s in reference_info["species"]]
default_params = reference_info.get("default_params") or {}

st.title("HCR v3.0 Probe Design")

# ---------------------------------------------------------------------------
# 1. Target sequence input
# ---------------------------------------------------------------------------
st.header("1. Target sequence")

input_mode = st.radio("Target sequence source", ["Paste or upload FASTA", "Search by gene name"])

if input_mode == "Search by gene name":
    ncbi_email = st.text_input(
        "NCBI contact email (required by NCBI to use gene search)",
        value="leed2@hhmi.org",
        key="ncbi_email",
    )
    geneFetch.configure_entrez(ncbi_email)

    common_species = {
        "": "",
        "Human": "Homo sapiens",
        "Mouse": "Mus musculus",
        "Zebrafish": "Danio rerio",
        "Fly": "Drosophila melanogaster",
        "C. elegans": "Caenorhabditis elegans",
        "Rat": "Rattus norvegicus",
    }
    species_shortcut = st.selectbox("Common species (optional, autofills below)", options=list(common_species.keys()))
    species_input = st.text_input(
        "Species", value=common_species.get(species_shortcut, ""), placeholder="e.g. Homo sapiens or mouse"
    )
    symbol_input = st.text_input("Gene symbol", placeholder="e.g. GAPDH")

    if st.button("Search"):
        if not species_input.strip() or not symbol_input.strip():
            st.error("Enter both a species and a gene symbol.")
        else:
            try:
                gene_matches = geneFetch.search_gene(species_input.strip(), symbol_input.strip())
            except Exception as exc:
                st.error(f"Gene search failed: {exc}")
                gene_matches = None
            if gene_matches is not None:
                if not gene_matches:
                    st.warning("No gene found. Check the species name/spelling.")
                    st.session_state.pop("gene_matches", None)
                else:
                    st.session_state["gene_matches"] = gene_matches
                    st.session_state["gene_species_for_fallback"] = species_input.strip()
                    st.session_state["gene_symbol_for_fallback"] = symbol_input.strip()

    if st.session_state.get("gene_matches"):
        gene_matches = st.session_state["gene_matches"]
        gene_labels = [f"{g['symbol']} — {g['description']} ({g['organism']})" for g in gene_matches]
        chosen_idx = st.selectbox(
            "Gene match", options=range(len(gene_labels)), format_func=lambda i: gene_labels[i]
        )
        chosen_gene = gene_matches[chosen_idx]

        try:
            transcripts = geneFetch.list_transcripts(
                chosen_gene["gene_id"],
                species=st.session_state.get("gene_species_for_fallback"),
                symbol=st.session_state.get("gene_symbol_for_fallback"),
            )
        except Exception as exc:
            st.error(f"Transcript lookup failed: {exc}")
            transcripts = []

        if not transcripts:
            st.warning("No RefSeq transcripts found for this gene.")
        else:
            if not transcripts[0]["is_refseq"]:
                st.caption("No curated RefSeq transcript found — showing broader nuccore search results instead.")
            transcript_labels = [
                f"{t['accession']} — {t['title']} ({t['length']} bp, "
                f"{'mRNA' if t['is_coding'] else 'non-coding RNA'})"
                for t in transcripts
            ]
            transcript_idx = st.selectbox(
                "Transcript", options=range(len(transcript_labels)), format_func=lambda i: transcript_labels[i]
            )
            if st.button("Use this transcript"):
                try:
                    fetched_fasta = geneFetch.fetch_transcript_fasta(transcripts[transcript_idx]["accession"])
                    st.session_state["fasta_text_input"] = fetched_fasta
                    st.rerun()
                except Exception as exc:
                    st.error(f"Could not fetch sequence: {exc}")

uploaded_file = st.file_uploader("Upload a FASTA file", type=["fa", "fasta", "txt"])
pasted_fasta = st.text_area(
    "...or paste FASTA text",
    height=150,
    placeholder=(
        ">MyGene\nACGT...\n\n"
        "Add 'channel=B2' to the header to override the sidebar channel for this "
        "record, e.g. '>MyGene channel=B2'."
    ),
    key="fasta_text_input",
)
fasta_text = uploaded_file.getvalue().decode("utf-8") if uploaded_file is not None else pasted_fasta

# ---------------------------------------------------------------------------
# 2. Parameters
# ---------------------------------------------------------------------------
st.sidebar.header("Parameters")
st.sidebar.caption(
    "A separate, fixed check (not adjustable here) discards any probe predicted to fold "
    "back on itself (an intramolecular hairpin with melting temp ≥ 45°C, via primer3)."
)
channel = st.sidebar.selectbox(
    "Channel",
    options=sorted(HCR.initiators.keys()),
    index=0,
    help=(
        "Which HCR amplifier/initiator pair (B1–B5) this gene's probes will trigger. "
        "Each channel is spectrally/chemically orthogonal, so give each gene you want to "
        "image simultaneously a different channel. Choosing B1 vs. B2 etc. doesn't itself "
        "change signal strength or specificity — it only determines which fluorophore-"
        "conjugated hairpin set will bind."
    ),
)
tile_size = st.sidebar.number_input(
    "Tile size",
    min_value=10,
    value=int(default_params.get("tileSize", 52)),
    step=1,
    help=(
        "Length (nt) of each candidate target region before it's split into two ~25-nt "
        "half-probes with a 2-nt gap (25+25+2=52). This is the standard HCR v3 split-probe "
        "size, empirically validated for stable, specific hybridization at the standard "
        "37°C probe hybridization step. Shrinking it makes each half shorter than ~25 nt "
        "(weaker, less specific binding, more signal dropout); growing it makes halves "
        "longer than validated (stronger binding, but not the tested/standard design and "
        "may not match published/commercial probe sets for this gene)."
    ),
)

st.sidebar.subheader("GC content (%)")
st.sidebar.caption(
    "A coarse, reference-only guardrail against extreme outliers. Gibbs FE below is what "
    "actually predicts hybridization stability, since it's computed from the real sequence "
    "rather than just a base count — treat these as a loose sanity range, not the real filter."
)
min_gc = st.sidebar.number_input(
    "Min GC",
    value=float(default_params.get("minGC", 35.0)),
    help="Lower bound on a candidate probe's GC%. Too low weakens binding at 37°C hybridization, risking signal loss.",
)
max_gc = st.sidebar.number_input(
    "Max GC",
    value=float(default_params.get("maxGC", 65.0)),
    help="Upper bound on a candidate probe's GC%. Too high can cause overly strong/sticky, less-specific binding and secondary structure.",
)
target_gc = st.sidebar.number_input(
    "Target GC",
    value=float(default_params.get("targetGC", 50.0)),
    help="Not currently used to rank probes (ranking uses Target Gibbs FE below); kept for reference/future use.",
)

st.sidebar.subheader("Gibbs free energy (kcal/mol)")
st.sidebar.caption(
    "This measures how tightly a probe binds its target mRNA (an intermolecular RNA:DNA duplex)."
)
with st.sidebar.expander("What do -50/-60/-70 actually mean?"):
    st.caption(
        "These numbers are hard to build intuition for on their own, so here's what this "
        "tool's own math (52-nt tile, 37°C, 0.33M salt) actually computes for probes of a "
        "given GC%:\n\n"
        "| GC% | Typical Gibbs FE |\n"
        "|---|---|\n"
        "| 20% | ≈ -41 |\n"
        "| 30% | ≈ -47 |\n"
        "| 40% | ≈ -58 |\n"
        "| 45% | ≈ -62 |\n"
        "| 50% | ≈ -67 |\n"
        "| 55% | ≈ -69 |\n"
        "| 70% | ≈ -83 |\n"
        "| 90% | ≈ -102 |\n\n"
        "So the default −70 to −50 range roughly tracks the same 45–55% GC window set "
        "above — it isn't an independent, unrelated knob. What it *adds* on top of GC%: at a "
        "**fixed** 45% GC, real probes in a test batch ranged from about -51 to -81 kcal/mol "
        "depending only on how the G/C bases were arranged (clustered vs. scattered) — GC% "
        "alone can't tell those apart, but Gibbs FE can, since it's computed from the actual "
        "sequence via nearest-neighbor stacking energies.\n\n"
        "Practically: moving the range's edges by about 10 kcal/mol shifts the effective "
        "composition window by roughly 10–15 percentage points of GC. If a gene is giving you "
        "too few candidate probes (common for AT-rich transcripts, e.g. many 3'UTRs or "
        "invertebrate genes), loosen Max Gibbs FE toward -40 to admit weaker/AT-richer "
        "probes; if you're worried about off-target stickiness from an unusually GC-rich "
        "region, tighten Min Gibbs FE toward -60."
    )
min_gibbs = st.sidebar.number_input(
    "Min Gibbs FE",
    value=float(default_params.get("minGibbs", -70.0)),
    help=(
        "Most-negative (strongest-binding) probe-to-target duplex allowed, computed from "
        "nearest-neighbor RNA:DNA thermodynamics at 37°C — a more accurate stability "
        "estimate than GC% alone, since it accounts for the actual sequence, not just base "
        "composition. This is NOT the hairpin (self-folding) check, which is separate and "
        "fixed — see the note above. Values more negative than this are excluded: "
        "duplexes that stable can bind even mismatched/near-target sequences avidly, "
        "reducing discrimination against paralogs or splice variants (more off-target "
        "signal)."
    ),
)
max_gibbs = st.sidebar.number_input(
    "Max Gibbs FE",
    value=float(default_params.get("maxGibbs", -50.0)),
    help=(
        "Least-negative (weakest-binding) probe allowed. Weaker than this risks the probe "
        "not staying reliably bound through the 37°C hybridization step and subsequent "
        "washes, so it gets rinsed away with the unbound background — signal dropout for "
        "that probe rather than off-target risk."
    ),
)
target_gibbs = st.sidebar.number_input(
    "Target Gibbs FE",
    value=float(default_params.get("targetGibbs", -60.0)),
    help=(
        "The actual ranking criterion: among probes passing every other filter, the ones "
        "closest to this value are chosen first (up to Max number of probes). Aiming for a "
        "moderate value near the middle of the Min/Max Gibbs range (rather than the "
        "strongest possible binder) keeps binding strength consistent across the whole "
        "probe set, so probes contribute more evenly to signal instead of a few dominating."
    ),
)

st.sidebar.subheader("Probe-half Tm matching")
dtm_filter = st.sidebar.checkbox(
    "Filter on dTm between probe halves",
    value=False,
    help=(
        "Each tile is split into two half-probes; sequence composition means their melting "
        "temperatures (Tm) can differ. HCR relies on BOTH halves binding at once to nucleate "
        "amplification (its 'AND-gate' specificity) — a large Tm mismatch means the weaker "
        "half may not stay bound during the standard 37°C hybridization, causing that probe "
        "pair to drop out (lower signal) or, worse, letting one half bind non-specifically "
        "elsewhere while contributing less to true-target selectivity (more off-target risk). "
        "Enabling this discards imbalanced pairs, trading fewer total probes for more uniform, "
        "reliable ones."
    ),
)
dtm_max = st.sidebar.number_input(
    "Max dTm",
    value=float(default_params.get("dTmMax", 5.0)),
    disabled=not dtm_filter,
    help=(
        "Maximum allowed Tm difference (°C) between a probe pair's two halves when the "
        "filter above is enabled. 5°C is a reasonable, protocol-safe default for 37°C "
        "hybridization; tightening it (lower value) further improves half-to-half balance "
        "but rejects more candidates, so genes with unusual sequence composition may end up "
        "with fewer usable probes."
    ),
)

st.sidebar.subheader("Homopolymer runs")
max_run_length = st.sidebar.number_input(
    "Max run length",
    min_value=1,
    value=int(default_params.get("maxRunLength", 7)),
    step=1,
    help=(
        "Rejects candidate probes containing a run of this many or more consecutive C's or "
        "G's. Long G-runs especially can fold into G-quadruplex/self-structures that block "
        "the probe from binding its target (reduced signal) and are also harder to "
        "synthesize accurately (more sequence errors in the ordered oligo). Lowering this is "
        "stricter (fewer, cleaner-behaving probes); raising it is more permissive (more "
        "candidate probes/better transcript coverage, slightly higher risk of individual "
        "probes underperforming)."
    ),
)
max_run_mismatches = st.sidebar.number_input(
    "Max run mismatches",
    min_value=0,
    value=int(default_params.get("maxRunMismatches", 2)),
    step=1,
    help=(
        "How many non-matching bases are tolerated inside an otherwise-long C or G run "
        "before it no longer counts as a 'run' (e.g. with a value of 2, 'CCCCACCC' still "
        "counts as one long C-run since a single A interrupts it). Higher values catch more "
        "near-runs (stricter filtering, same rationale as Max run length above)."
    ),
)

max_probes = st.sidebar.number_input(
    "Max number of probes",
    min_value=1,
    value=int(default_params.get("maxProbes", 20)),
    step=1,
    help=(
        "Maximum number of non-overlapping probe pairs to keep for this target, chosen by "
        "closeness to Target Gibbs FE. More probes generally means brighter, more robust "
        "signal (HCR signal scales with probe pair count) at the cost of needing a longer "
        "target region and more synthesis cost; fewer probes risk weak/no visible signal if "
        "the transcript or region is short."
    ),
)

st.sidebar.subheader("Genome off-target masking")
genome_mask_on = st.sidebar.checkbox(
    "Enable genome off-target masking (requires bowtie2)",
    value=BOWTIE2_AVAILABLE,
    help=(
        "Checks each candidate probe against a registered genome index via bowtie2 "
        "and discards probes with too many genomic hits."
    ),
)
species = None
index_override = ""
num_hits_allowed = int(default_params.get("num_hits_allowed", 1))
if genome_mask_on:
    if not BOWTIE2_AVAILABLE:
        st.sidebar.warning("bowtie2 was not found on PATH — this will fail unless you provide an explicit index.")
    species_options = species_names or ["(none registered)"]
    _prefs_config = referenceGenome.load_config()
    _last_species = (_prefs_config.get("gui_preferences") or {}).get("last_species")
    _default_species_index = species_options.index(_last_species) if _last_species in species_options else 0

    def _save_species_preference():
        cfg = referenceGenome.load_config()
        cfg.setdefault("gui_preferences", {})["last_species"] = st.session_state["species_select"]
        referenceGenome.save_config(cfg)

    species = st.sidebar.selectbox(
        "Species",
        options=species_options,
        index=_default_species_index,
        key="species_select",
        on_change=_save_species_preference,
        help="Registered reference genome to align candidate probes against for off-target checking. Remembers your last choice across restarts.",
    )
    index_override = st.sidebar.text_input(
        "Explicit bowtie2 index prefix (overrides species)",
        value="",
        help="Optional: skip the species dropdown and point directly at a bowtie2 index path prefix.",
    )
    num_hits_allowed = st.sidebar.number_input(
        "Max allowed genomic hits",
        min_value=0,
        value=num_hits_allowed,
        step=1,
        help=(
            "How many places in the genome a candidate probe is allowed to align to before "
            "it's rejected as non-specific. Keep at 1 (default, matches the underlying "
            "package's own default) for standard single-copy genes — the strictest, most "
            "specific setting. Consider raising to 2–3 only if your target gene has a "
            "known close paralog or belongs to a recently-duplicated gene family, where part "
            "of the target sequence is shared with that other genomic locus (not the same as "
            "having multiple splice isoforms — those share one locus and never trigger this "
            "filter) and a fully gene-specific probe may not be achievable there; this trades "
            "stricter specificity for keeping otherwise-good probes, at the cost of tolerating "
            "some cross-reactivity with that paralogous locus."
        ),
    )

target_name = st.sidebar.text_input(
    "Target name (used for output file naming)",
    value="target",
    help="Used to name the downloaded TSV/IDT files and internal scratch files for this run — doesn't affect probe design itself.",
)

# ---------------------------------------------------------------------------
# 3. Run
# ---------------------------------------------------------------------------
st.header("2. Design probes")
run_clicked = st.button("Design probes", type="primary")

if run_clicked:
    if not fasta_text.strip():
        st.error("Please upload or paste a FASTA sequence first.")
        st.stop()

    fasta_iterator = sequencelib.FastaIterator(io.StringIO(fasta_text))
    try:
        record = next(fasta_iterator, None)
    except ValueError as exc:
        st.error(f"Could not parse FASTA input: {exc}")
        st.stop()

    if record is None:
        st.error("No FASTA record found. Input must start with a '>' header line.")
        st.stop()

    cleaned_name, channel_override = probeDesign._parse_record_channel(record["name"])
    record = {"name": cleaned_name or record["name"], "sequence": record["sequence"]}

    args = types.SimpleNamespace(
        verbose=False,
        channel=channel,
        tileSize=int(tile_size),
        targetName=target_name,
        species=species,
        minGC=float(min_gc),
        maxGC=float(max_gc),
        targetGC=float(target_gc),
        dTmMax=float(dtm_max),
        dTmFilter=bool(dtm_filter),
        # NOTE: upstream --no-genomemask is inverted from its name (default True means
        # masking RUNS). The checkbox's boolean value maps directly onto this flag.
        no_genomemask=bool(genome_mask_on),
        index=index_override or None,
        # --no-repeatmask is permanently broken upstream (argparse default/action bug
        # makes it always False), so repeat masking never actually runs. Not exposed here.
        no_repeatmask=False,
        minGibbs=float(min_gibbs),
        maxGibbs=float(max_gibbs),
        targetGibbs=float(target_gibbs),
        maxRunLength=int(max_run_length),
        maxProbes=int(max_probes),
        maxRunMismatches=int(max_run_mismatches),
        num_hits_allowed=int(num_hits_allowed),
    )

    if args.no_genomemask and not args.index and args.species not in species_names:
        st.error(
            f"Species '{args.species}' is not registered. Disable genome masking, "
            "provide an explicit index prefix, or register a species below."
        )
        st.stop()
    if args.no_genomemask and not args.index and not BOWTIE2_AVAILABLE:
        st.error(
            "bowtie2 not found on PATH. Uncheck 'Enable genome off-target masking' "
            "or install bowtie2 in this environment."
        )
        st.stop()

    best_tiles = None
    with st.spinner("Designing probes..."):
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                previous_cwd = os.getcwd()
                os.chdir(tmpdir)
                try:
                    # _design_tiles_for_record is the pipeline's real integration point;
                    # it writes genome-masking scratch files into the CWD, hence the chdir.
                    best_tiles = probeDesign._design_tiles_for_record(
                        args, record, target_name, channel_override
                    )
                finally:
                    os.chdir(previous_cwd)
        except ValueError as exc:
            st.error(f"Input error: {exc}")
        except SystemExit as exc:
            st.error(str(exc))
        except FileNotFoundError as exc:
            st.error(f"Required external tool not found: {exc}")

    if best_tiles is not None:
        if not best_tiles:
            st.warning("No probes passed all filters. Try relaxing the GC / Gibbs / dTm ranges.")
        else:
            st.session_state["best_tiles"] = best_tiles
            st.session_state["target_name"] = target_name

# ---------------------------------------------------------------------------
# 4. Results
# ---------------------------------------------------------------------------
if "best_tiles" in st.session_state:
    best_tiles = st.session_state["best_tiles"]
    result_target_name = st.session_state.get("target_name", "target")

    st.header("3. Results")

    table_handle = io.StringIO()
    probeDesign.outputTable(best_tiles, outHandle=table_handle)
    results_df = pd.read_csv(io.StringIO(table_handle.getvalue()), sep="\t")
    st.dataframe(results_df, use_container_width=True)

    # Use the channel actually baked into the designed tiles (not the current sidebar widget
    # value, which may have changed since this design was run) for a correct pool name.
    pool_name = f"{result_target_name}_{best_tiles[0].channel}_pool"

    idt_buffer = io.BytesIO()
    probeDesign.write_idt_opool_xlsx(best_tiles, pool_name, idt_buffer)

    col1, col2 = st.columns(2)
    with col1:
        st.download_button(
            "Download probes.tsv",
            data=table_handle.getvalue(),
            file_name=f"{result_target_name}_probes.tsv",
            mime="text/tab-separated-values",
        )
    with col2:
        st.download_button(
            "Download IDT oPools order sheet (.xlsx)",
            data=idt_buffer.getvalue(),
            file_name=f"{pool_name}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

    n_oligos = 2 * len(best_tiles)
    total_pmol_10 = n_oligos * 10
    stock_uM = 4.0
    resuspend_uL = total_pmol_10 / stock_uM
    working_nM = 4.0
    hyb_vol_uL = 500.0
    dilution_uL = (working_nM / 1000 * hyb_vol_uL) / stock_uM

    st.subheader("Order this pool")
    st.link_button("Order this pool at IDT oPools →", "https://www.idtdna.com/site/order/poolentry/")
    st.caption(
        f"Download the .xlsx file above and upload it in IDT's oPools Oligo Pools Entry page. "
        f"Available synthesis scales: 1, 10, or 50 pmol/oligo.\n\n"
        f"**Example resuspension → dilution** (10 pmol/oligo scale, this pool's {len(best_tiles)} "
        f"probe pairs = {n_oligos} oligos; shown for illustration, not a mandatory protocol):\n"
        f"1. Total oligo ordered: {n_oligos} oligos × 10 pmol = {total_pmol_10:.0f} pmol\n"
        f"2. Resuspend to a **{stock_uM:.0f} µM stock**: {total_pmol_10:.0f} pmol ÷ {stock_uM:.0f} µM = "
        f"**{resuspend_uL:.0f} µL** of nuclease-free water or IDTE buffer (conventional HCR "
        f"protocols typically use a 1 µM stock instead — ours is 4× more concentrated, using less "
        f"resuspension volume and lasting longer)\n"
        f"3. Dilute for routine use to the standard **{working_nM:.0f} nM working concentration** "
        f"(for a typical {hyb_vol_uL:.0f} µL hybridization volume):  \n"
        f"{working_nM:.0f} nM × {hyb_vol_uL:.0f} µL ÷ {stock_uM:.0f} µM = **{dilution_uL:.2g} µL** "
        f"of stock added to ~{hyb_vol_uL:.0f} µL of hybridization buffer  \n"
        f"(a conventional 1 µM stock would need {dilution_uL * stock_uM:.2g} µL instead for the "
        f"same dose — if that's too small a volume to pipette accurately, pre-dilute the stock "
        f"1:10 first).\n\n"
        f"Ordered a different scale (1 or 50 pmol/oligo)? Multiply step 1 by that number instead of 10."
    )

# ---------------------------------------------------------------------------
# 5. Species registration (advanced, infrequent)
# ---------------------------------------------------------------------------
st.divider()
with st.expander("Register a new reference genome (advanced)", expanded=False):
    force_overwrite = st.checkbox("Force overwrite if species/index already exists", key="register_force")
    source = st.radio(
        "Source",
        ["Common genome (fastest)", "Genome URL", "Upload file(s)"],
        key="genome_source",
    )

    def _make_progress_callback(bar, caption):
        def _callback(downloaded, total):
            if total and total > 0:
                bar.progress(min(downloaded / total, 1.0))
                caption.caption(f"{downloaded / 1e6:.0f} MB / {total / 1e6:.0f} MB")
            else:
                caption.caption(f"{downloaded / 1e6:.0f} MB downloaded")
        return _callback

    # --- Tier 1: one-click prebuilt Bowtie2 index (no bowtie2-build needed) ---
    if source == "Common genome (fastest)":
        st.caption("Downloads a prebuilt Bowtie2 index directly — no bowtie2-build required.")
        genome_label = st.selectbox("Genome", options=list(genomeFetch.PREBUILT_INDEXES.keys()))
        genome_choice = genomeFetch.PREBUILT_INDEXES[genome_label]
        species_key_input = st.text_input(
            "Species name", value=genome_choice["species"], key="prebuilt_species_name"
        )

        if st.button("Download and register", key="prebuilt_register_button"):
            species_key = species_key_input.strip()
            if not species_key:
                st.error("Enter a species name.")
            else:
                progress_bar = st.progress(0.0)
                progress_caption = st.empty()
                try:
                    index_prefix = genomeFetch.fetch_prebuilt_index(
                        genome_choice["key"],
                        species_key,
                        force=force_overwrite,
                        progress_callback=_make_progress_callback(progress_bar, progress_caption),
                    )
                    referenceGenome.register_species(
                        species=species_key, index_prefix=index_prefix, force=force_overwrite
                    )
                    st.success(f"Registered species '{species_key}'. Reload the page to see it in the Species dropdown.")
                except (FileExistsError, ValueError) as exc:
                    st.error(f"{exc} (check 'Force overwrite' to replace it.)")
                except (urllib.error.URLError, OSError) as exc:
                    st.error(f"Download failed: {exc}")

    # --- Tier 2: download a FASTA from a direct URL, then build locally ---
    elif source == "Genome URL":
        st.caption(
            "Paste a direct link to a genome FASTA (e.g. an Ensembl or UCSC .fa/.fa.gz URL). "
            "The server downloads it directly — no browser upload of a multi-GB file."
        )
        if not BOWTIE2_BUILD_AVAILABLE:
            st.warning("bowtie2-build was not found on PATH. Install bowtie2 in this environment first.")

        suggest_species = st.text_input(
            "Suggest a URL for species (optional, UCSC genome database name, e.g. 'hg38' or 'mm39')",
            key="ucsc_suggest_db",
        )
        default_url = ""
        if suggest_species.strip():
            db = suggest_species.strip()
            default_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{db}/bigZips/{db}.fa.gz"
        genome_url = st.text_input("Genome FASTA URL", value=default_url, key="genome_url")
        new_species_name = st.text_input("Species name", key="url_species_name")
        threads = st.number_input("Threads", min_value=1, value=1, step=1, key="url_register_threads")

        if st.button("Download, build index, and register", key="url_register_button"):
            species_key = new_species_name.strip()
            if not genome_url.strip():
                st.error("Enter a genome FASTA URL.")
            elif not species_key:
                st.error("Enter a species name.")
            elif not BOWTIE2_BUILD_AVAILABLE:
                st.error("bowtie2-build not found on PATH.")
            else:
                progress_bar = st.progress(0.0)
                progress_caption = st.empty()
                try:
                    with tempfile.TemporaryDirectory() as tmpdir:
                        fasta_path = genomeFetch.download_genome_fasta(
                            genome_url.strip(),
                            tmpdir,
                            progress_callback=_make_progress_callback(progress_bar, progress_caption),
                        )
                        fasta_paths = referenceGenome.collect_fasta_inputs([fasta_path])
                        with st.spinner(f"Building bowtie2 index for '{species_key}'... this may take a while."):
                            index_prefix = referenceGenome.build_bowtie2_index(
                                fasta_paths, species_key, threads=int(threads), force=force_overwrite
                            )
                        referenceGenome.register_species(
                            species=species_key, index_prefix=index_prefix, force=force_overwrite
                        )
                    st.success(f"Registered species '{species_key}'. Reload the page to see it in the Species dropdown.")
                except RuntimeError as exc:
                    st.error(str(exc))
                except (FileExistsError, ValueError) as exc:
                    st.error(f"{exc} (check 'Force overwrite' to replace it.)")
                except (urllib.error.URLError, OSError) as exc:
                    st.error(f"Download failed: {exc}")

    # --- Tier 3: manual upload (small/custom sequences only) ---
    else:
        st.caption(
            "For small or custom sequences (e.g. a transgene or plasmid), not a full genome — "
            "browser uploads are capped at 200MB, too small for most real genomes."
        )
        if not BOWTIE2_BUILD_AVAILABLE:
            st.warning("bowtie2-build was not found on PATH. Install bowtie2 in this environment first.")

        genome_files = st.file_uploader(
            "Genome FASTA file(s)",
            type=["fa", "fasta", "fna", "gz"],
            accept_multiple_files=True,
            key="genome_upload",
        )
        new_species_name = st.text_input("Species name", key="upload_species_name")
        threads = st.number_input("Threads", min_value=1, value=1, step=1, key="upload_register_threads")

        if st.button("Build index and register species", key="upload_register_button"):
            if not genome_files:
                st.error("Upload at least one genome FASTA file.")
            elif not new_species_name.strip():
                st.error("Enter a species name.")
            elif not BOWTIE2_BUILD_AVAILABLE:
                st.error("bowtie2-build not found on PATH.")
            else:
                species_key = new_species_name.strip()
                with st.spinner(f"Building bowtie2 index for '{species_key}'... this may take a while."):
                    try:
                        with tempfile.TemporaryDirectory() as tmpdir:
                            saved_paths = []
                            for uploaded_genome_file in genome_files:
                                path = os.path.join(tmpdir, uploaded_genome_file.name)
                                with open(path, "wb") as fh:
                                    fh.write(uploaded_genome_file.getvalue())
                                saved_paths.append(path)
                            fasta_paths = referenceGenome.collect_fasta_inputs(saved_paths)
                            index_prefix = referenceGenome.build_bowtie2_index(
                                fasta_paths, species_key, threads=int(threads), force=force_overwrite
                            )
                            referenceGenome.register_species(
                                species=species_key, index_prefix=index_prefix, force=force_overwrite
                            )
                        st.success(f"Registered species '{species_key}'. Reload the page to see it in the Species dropdown.")
                    except RuntimeError as exc:
                        st.error(str(exc))
                    except (FileExistsError, ValueError) as exc:
                        st.error(f"{exc} (check 'Force overwrite' to replace it.)")
                    except FileNotFoundError as exc:
                        st.error(str(exc))

# ---------------------------------------------------------------------------
# 6. Manage registered species (advanced, infrequent)
# ---------------------------------------------------------------------------
with st.expander("Manage registered species (advanced)", expanded=False):
    registered_species = listReferences.list_references()["species"]
    if not registered_species:
        st.caption("No species registered yet.")
    for sp in registered_species:
        with st.container(border=True):
            status = "installed" if sp["installed"] else ":warning: index files missing"
            st.markdown(f"**{sp['name']}** — `{sp['bowtie2_index']}` ({status})")

            rename_col, unregister_col, delete_col = st.columns([2, 1, 1])
            with rename_col:
                new_name = st.text_input(
                    "Rename to", key=f"rename_input_{sp['name']}", placeholder=sp["name"], label_visibility="collapsed"
                )
                if st.button("Rename", key=f"rename_btn_{sp['name']}"):
                    if not new_name.strip():
                        st.error("Enter a new name.")
                    else:
                        try:
                            referenceGenome.rename_species(sp["name"], new_name.strip())
                            if st.session_state.get("species_select") == sp["name"]:
                                # Widgets can't be reassigned via session_state after they've
                                # instantiated this run (raises StreamlitWidgetAlreadyInstantiatedError).
                                # Instead, persist the new name as the preferred default and clear
                                # the live selection so the widget re-picks it from that preference
                                # on rerun -- .pop() is unaffected by that restriction.
                                cfg = referenceGenome.load_config()
                                cfg.setdefault("gui_preferences", {})["last_species"] = new_name.strip()
                                referenceGenome.save_config(cfg)
                                st.session_state.pop("species_select", None)
                            st.success(f"Renamed '{sp['name']}' to '{new_name.strip()}'.")
                            st.rerun()
                        except ValueError as exc:
                            st.error(str(exc))
            with unregister_col:
                if st.button(
                    "Unregister",
                    key=f"unregister_btn_{sp['name']}",
                    help="Remove from this list but keep the index files on disk.",
                ):
                    referenceGenome.delete_species(sp["name"], delete_files=False)
                    if st.session_state.get("species_select") == sp["name"]:
                        st.session_state.pop("species_select", None)
                    st.success(f"Unregistered '{sp['name']}' (index files kept on disk).")
                    st.rerun()
            with delete_col:
                confirm_delete = st.checkbox(
                    "Confirm",
                    key=f"confirm_delete_{sp['name']}",
                    help="Check to enable permanently deleting this species' index files from disk.",
                )
                if st.button(
                    "Delete + files",
                    key=f"delete_btn_{sp['name']}",
                    disabled=not confirm_delete,
                    help="Also deletes the index files from disk to free space. Cannot be undone.",
                ):
                    deleted_dir = referenceGenome.delete_species(sp["name"], delete_files=True)
                    if st.session_state.get("species_select") == sp["name"]:
                        st.session_state.pop("species_select", None)
                    if deleted_dir:
                        st.success(f"Deleted '{sp['name']}' and removed {deleted_dir}.")
                    else:
                        st.success(
                            f"Unregistered '{sp['name']}' (index files were outside the managed "
                            "indices directory, e.g. a custom --index path, and were left in place)."
                        )
                    st.rerun()
