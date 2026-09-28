import os

from HCRProbeDesign import genomeMask as gm
from HCRProbeDesign import _datadir


def test_genomemask_uses_absolute_index(monkeypatch, tmp_path):
    captured = {}

    def fake_call(cmd):
        captured["cmd"] = cmd
        return 0

    monkeypatch.setattr(gm.subprocess, "call", fake_call)
    monkeypatch.chdir(tmp_path)

    gm.genomemask(">read1\nACGT\n", handleName="test", index="/abs/index")

    idx = captured["cmd"].index("-x") + 1
    assert captured["cmd"][idx] == "/abs/index"


def test_genomemask_uses_relative_index(monkeypatch, tmp_path):
    captured = {}

    def fake_call(cmd):
        captured["cmd"] = cmd
        return 0

    monkeypatch.setattr(gm.subprocess, "call", fake_call)
    data_dir = str(tmp_path / "datadir")
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", data_dir)
    monkeypatch.chdir(tmp_path)

    gm.genomemask(">read1\nACGT\n", handleName="test", index="indices/mm10/mm10")

    idx = captured["cmd"].index("-x") + 1
    assert captured["cmd"][idx] == os.path.join(data_dir, "indices/mm10/mm10")


def test_genomemask_writes_scratch_files_into_workdir(monkeypatch, tmp_path):
    # Regression test: genomemask() must place its scratch FASTA/SAM files
    # under an explicit workdir rather than the process's CWD, so that a
    # long-running multi-user server (e.g. the Streamlit GUI) can run
    # concurrent design jobs in their own directories without one job's
    # os.chdir() clobbering another's in-flight scratch files.
    def fake_call(cmd):
        return 0

    monkeypatch.setattr(gm.subprocess, "call", fake_call)
    other_cwd = tmp_path / "unrelated_cwd"
    other_cwd.mkdir()
    monkeypatch.chdir(other_cwd)

    workdir = tmp_path / "job_workdir"
    workdir.mkdir()

    gm.genomemask(">read1\nACGT\n", handleName="test", index="/abs/index", workdir=str(workdir))

    assert (workdir / "test_reads.fa").exists()
    assert not (other_cwd / "test_reads.fa").exists()
