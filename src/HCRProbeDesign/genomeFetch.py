"""Fetch reference genomes without requiring a manual browser upload.

Two paths, both server-side downloads (never through a browser file upload):

- ``fetch_prebuilt_index`` downloads a prebuilt Bowtie2 index from the Bowtie2
  project's public index bucket, skipping ``bowtie2-build`` entirely.
- ``download_genome_fasta`` streams an arbitrary FASTA/FASTA.gz URL (e.g. an
  Ensembl or UCSC direct link) to disk, to be handed to
  ``referenceGenome.build_bowtie2_index``.
"""

import glob
import os
import shutil
import tempfile
import urllib.request
from zipfile import ZipFile

from . import index_path

PREBUILT_INDEX_BASE_URL = "https://genome-idx.s3.amazonaws.com/bt"

# Verified via HEAD request against PREBUILT_INDEX_BASE_URL -- do not add a key
# here without checking it resolves first, the naming is not derivable by pattern.
PREBUILT_INDEXES = {
    "Mouse (mm10)": {"key": "mm10", "species": "mouse"},
    "Mouse (GRCm39)": {"key": "GRCm39", "species": "mouse"},
    "Human (GRCh38, no-alt)": {"key": "GRCh38_noalt_as", "species": "human"},
    "Human (hg19)": {"key": "hg19", "species": "human"},
    "Zebrafish (GRCz11)": {"key": "GRCz11", "species": "zebrafish"},
    "Fly (BDGP6)": {"key": "BDGP6", "species": "fly"},
    "C. elegans (WBcel235)": {"key": "WBcel235", "species": "celegans"},
    "Rat (Rnor_6.0)": {"key": "Rnor_6.0", "species": "rat"},
}


def download_with_progress(url, dest_path, progress_callback=None):
    """
    Download a URL to a local path, optionally reporting progress.

    :param url: Source URL.
    :param dest_path: Local file path to write to.
    :param progress_callback: Optional callable(bytes_downloaded, total_bytes).
    :return: None.
    """
    def _reporthook(block_count, block_size, total_size):
        if progress_callback is not None:
            progress_callback(min(block_count * block_size, total_size if total_size > 0 else block_count * block_size), total_size)

    urllib.request.urlretrieve(url, dest_path, reporthook=_reporthook if progress_callback else None)


def fetch_prebuilt_index(genome_key, species, indices_dir=None, force=False, progress_callback=None):
    """
    Download and extract a prebuilt Bowtie2 index for a common genome.

    :param genome_key: Key from PREBUILT_INDEXES (e.g. "mm10").
    :param species: Species name to register the index under.
    :param indices_dir: Output directory for indices (default: user data dir).
    :param force: Overwrite an existing index for this species if True.
    :param progress_callback: Optional callable(bytes_downloaded, total_bytes).
    :return: Index prefix path.
    :raises FileExistsError: If an index already exists for this species and force is False.
    """
    indices_dir = indices_dir or index_path()
    species_dir = os.path.join(indices_dir, species)

    existing = glob.glob(os.path.join(species_dir, "*.bt2*"))
    if existing and not force:
        raise FileExistsError(f"Index already exists at {species_dir}. Use force=True to overwrite.")
    if existing and force:
        for fname in existing:
            os.remove(fname)

    os.makedirs(species_dir, exist_ok=True)
    url = f"{PREBUILT_INDEX_BASE_URL}/{genome_key}.zip"

    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = os.path.join(tmpdir, f"{genome_key}.zip")
        download_with_progress(url, zip_path, progress_callback=progress_callback)
        with ZipFile(zip_path, "r") as archive:
            archive.extractall(species_dir)

    index_prefix = os.path.join(species_dir, genome_key)
    return index_prefix


def download_genome_fasta(url, dest_dir, progress_callback=None):
    """
    Stream a FASTA (or FASTA.gz) URL directly to disk.

    :param url: Direct URL to a FASTA or FASTA.gz file.
    :param dest_dir: Directory to save the downloaded file into.
    :param progress_callback: Optional callable(bytes_downloaded, total_bytes).
    :return: Local path to the downloaded file.
    """
    filename = url.rstrip("/").rsplit("/", 1)[-1] or "genome.fa"
    dest_path = os.path.join(dest_dir, filename)
    download_with_progress(url, dest_path, progress_callback=progress_callback)
    return dest_path
