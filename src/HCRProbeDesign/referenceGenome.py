"""Build and register Bowtie2 indices for reference genomes."""

import argparse
import glob
import os
import shutil
import subprocess
import yaml

from . import index_path
from ._datadir import get_data_dir, get_config_path, get_indices_dir, ensure_data_dir

PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def load_config(config_path=None):
    """
    Load the HCRconfig.yaml file.

    :param config_path: Path to the YAML configuration file.
    :return: Parsed config dictionary (empty if missing).
    """
    if config_path is None:
        ensure_data_dir()
        config_path = get_config_path()
    if not os.path.exists(config_path):
        return {}
    with open(config_path, "r") as handle:
        return yaml.safe_load(handle) or {}


def save_config(config, config_path=None):
    """
    Write configuration data to HCRconfig.yaml.

    :param config: Configuration dictionary.
    :param config_path: Path to write the configuration.
    :return: None.
    """
    if config_path is None:
        ensure_data_dir()
        config_path = get_config_path()
    with open(config_path, "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)


def collect_fasta_inputs(paths):
    """
    Resolve FASTA inputs from files and directories.

    :param paths: List of FASTA files or directories.
    :return: Deduplicated list of FASTA file paths.
    :raises FileNotFoundError: If a path or directory has no FASTA files.
    """
    files = []
    for path in paths:
        if os.path.isdir(path):
            matches = []
            for ext in FASTA_EXTENSIONS:
                matches.extend(sorted(glob.glob(os.path.join(path, f"*{ext}"))))
            if not matches:
                raise FileNotFoundError(f"No FASTA files found in directory: {path}")
            files.extend(matches)
        else:
            if not os.path.exists(path):
                raise FileNotFoundError(f"FASTA path not found: {path}")
            files.append(path)

    seen = set()
    deduped = []
    for path in files:
        if path not in seen:
            seen.add(path)
            deduped.append(path)
    return deduped


def format_index_path(index_prefix):
    """
    Format an index prefix relative to the data directory when possible.

    :param index_prefix: Bowtie2 index prefix path.
    :return: Relative path if within data dir, otherwise absolute path.
    """
    abs_prefix = os.path.abspath(index_prefix)
    abs_root = os.path.abspath(get_data_dir())
    try:
        common = os.path.commonpath([abs_prefix, abs_root])
    except ValueError:
        common = None
    if common == abs_root:
        return os.path.relpath(abs_prefix, abs_root)
    return abs_prefix


def build_bowtie2_index(fasta_paths, species, index_name=None, indices_dir=None, threads=1, force=False, large_index=False):
    """
    Build a Bowtie2 index from the provided FASTA files.

    :param fasta_paths: List of FASTA file paths.
    :param species: Species name for the index directory.
    :param index_name: Optional index basename override.
    :param indices_dir: Output directory for indices.
    :param threads: Number of threads for bowtie2-build.
    :param force: Overwrite existing index files if True.
    :param large_index: Use --large-index for genomes > 4 billion bases.
    :return: Index prefix path.
    :raises RuntimeError: If bowtie2-build is not available.
    :raises FileExistsError: If index exists and force is False.
    """
    bowtie2_build = shutil.which("bowtie2-build")
    if not bowtie2_build:
        raise RuntimeError("bowtie2-build not found in PATH")

    indices_dir = indices_dir or index_path()
    species_dir = os.path.join(indices_dir, species)
    os.makedirs(species_dir, exist_ok=True)

    index_name = index_name or species
    index_prefix = os.path.join(species_dir, index_name)
    existing = glob.glob(f"{index_prefix}*.bt2*")
    if existing and not force:
        raise FileExistsError(f"Index already exists at {index_prefix}. Use --force to overwrite.")
    if existing and force:
        for fname in existing:
            os.remove(fname)

    cmd = [bowtie2_build]
    if large_index:
        cmd.append("--large-index")
    if threads and threads > 1:
        cmd.extend(["--threads", str(threads)])
    cmd.extend([",".join(fasta_paths), index_prefix])
    subprocess.check_call(cmd)
    return index_prefix


def register_species(config_path=None, species=None, index_prefix=None, force=False):
    """
    Register a species and its Bowtie2 index prefix in the config file.

    :param config_path: Path to HCRconfig.yaml (default: user data dir).
    :param species: Species key to register.
    :param index_prefix: Bowtie2 index prefix path.
    :param force: Overwrite an existing species entry if True.
    :return: None.
    :raises ValueError: If the species exists and force is False.
    """
    if config_path is None:
        ensure_data_dir()
        config_path = get_config_path()
    config = load_config(config_path)
    species_config = config.setdefault("species", {})
    if species in species_config and not force:
        raise ValueError(f"Species '{species}' already exists in config. Use --force to replace.")
    species_config[species] = {"bowtie2_index": format_index_path(index_prefix)}
    save_config(config, config_path)


def _resolve_absolute_index_path(index_prefix):
    """
    Resolve an index prefix to an absolute path.

    :param index_prefix: Bowtie2 index prefix (relative or absolute).
    :return: Absolute path to the index prefix.
    """
    if os.path.isabs(index_prefix):
        return index_prefix
    return os.path.join(get_data_dir(), index_prefix)


def rename_species(old_name, new_name, config_path=None, force=False):
    """
    Rename a registered species in the config file.

    Only updates the config key; does not move or rename any index files on disk.

    :param old_name: Currently registered species key.
    :param new_name: New species key.
    :param config_path: Path to HCRconfig.yaml (default: user data dir).
    :param force: Overwrite new_name's existing entry if True.
    :return: None.
    :raises ValueError: If old_name isn't registered, or new_name exists and force is False.
    """
    if config_path is None:
        ensure_data_dir()
        config_path = get_config_path()
    config = load_config(config_path)
    species_config = config.get("species", {}) or {}
    if old_name not in species_config:
        raise ValueError(f"Species '{old_name}' is not registered.")
    if new_name in species_config and not force:
        raise ValueError(f"Species '{new_name}' already exists. Use force=True to overwrite.")
    species_config[new_name] = species_config.pop(old_name)
    save_config(config, config_path)


def delete_species(name, config_path=None, delete_files=False):
    """
    Unregister a species, optionally deleting its index files from disk.

    Index files are only deleted if they live inside the managed indices
    directory (i.e. were created via build_bowtie2_index/fetch_prebuilt_index),
    to avoid removing a directory the user pointed at explicitly with --index.

    :param name: Registered species key to remove.
    :param config_path: Path to HCRconfig.yaml (default: user data dir).
    :param delete_files: Also delete the index directory on disk if True.
    :return: The deleted directory path if files were removed, else None.
    :raises ValueError: If name isn't registered.
    """
    if config_path is None:
        ensure_data_dir()
        config_path = get_config_path()
    config = load_config(config_path)
    species_config = config.get("species", {}) or {}
    if name not in species_config:
        raise ValueError(f"Species '{name}' is not registered.")

    index_prefix = species_config[name].get("bowtie2_index", "")
    abs_prefix = _resolve_absolute_index_path(index_prefix) if index_prefix else None

    del species_config[name]
    save_config(config, config_path)

    if not delete_files or not abs_prefix:
        return None

    species_dir = os.path.dirname(abs_prefix)
    indices_root = os.path.abspath(get_indices_dir())
    if os.path.isdir(species_dir) and os.path.commonpath([os.path.abspath(species_dir), indices_root]) == indices_root:
        shutil.rmtree(species_dir)
        return species_dir
    return None


def main():
    """CLI entry point for building and registering a reference genome index."""
    parser = argparse.ArgumentParser(
        description="Build a Bowtie2 index from a reference genome and register it."
    )
    parser.add_argument("--species", required=True, help="Species key to register in HCRconfig.yaml")
    parser.add_argument(
        "--fasta",
        required=True,
        action="append",
        help="FASTA file or directory (repeatable for multiple inputs)",
    )
    parser.add_argument("--index-name", help="Index basename (default: species)")
    ensure_data_dir()
    parser.add_argument("--indices-dir", help="Output directory for indices (default: user data dir)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for bowtie2-build")
    parser.add_argument("--config", default=get_config_path(), help="Path to HCRconfig.yaml")
    parser.add_argument("--force", action="store_true", help="Overwrite existing index/config entry")
    parser.add_argument("--large-index", action="store_true", help="Build a large index (for genomes > 4 billion bases)")
    args = parser.parse_args()

    fasta_paths = collect_fasta_inputs(args.fasta)
    index_prefix = build_bowtie2_index(
        fasta_paths,
        args.species,
        index_name=args.index_name,
        indices_dir=args.indices_dir,
        threads=args.threads,
        force=args.force,
        large_index=args.large_index,
    )
    register_species(args.config, args.species, index_prefix, force=args.force)
    print(f"Registered {args.species} with index {format_index_path(index_prefix)}")
