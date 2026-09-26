"""Look up genes and fetch mRNA transcript sequences from NCBI, so a target
sequence for probe design can be found by species + gene symbol instead of
requiring the user to already have a FASTA file.
"""

from Bio import Entrez


def configure_entrez(email, tool="HCRProbeDesign-GUI"):
    """
    Set the identifying email/tool NCBI's usage policy requires on Entrez requests.

    :param email: Contact email to identify the requester to NCBI.
    :param tool: Tool name reported to NCBI.
    :return: None.
    """
    Entrez.email = email
    Entrez.tool = tool


def search_gene(species, symbol, retmax=10):
    """
    Search NCBI Gene for a symbol within a species.

    :param species: Species name (scientific or common, e.g. "Homo sapiens" or "mouse").
    :param symbol: Gene symbol (e.g. "GAPDH").
    :return: List of dicts with gene_id, symbol, description, organism.
    """
    handle = Entrez.esearch(db="gene", term=f"{symbol}[sym] AND {species}[orgn]", retmax=retmax)
    search_result = Entrez.read(handle)
    handle.close()

    gene_ids = search_result.get("IdList", [])
    if not gene_ids:
        return []

    handle = Entrez.esummary(db="gene", id=",".join(gene_ids))
    summaries = Entrez.read(handle)
    handle.close()

    results = []
    for doc in summaries["DocumentSummarySet"]["DocumentSummary"]:
        results.append({
            "gene_id": str(doc.attributes["uid"]),
            "symbol": str(doc.get("Name", "")),
            "description": str(doc.get("Description", "")),
            "organism": str(doc.get("Organism", {}).get("ScientificName", "")),
        })
    return results


def list_transcripts(gene_id, species=None, symbol=None):
    """
    List RefSeq mRNA/ncRNA transcripts linked to a gene.

    Falls back to a broader nuccore search (requires species and symbol) if the
    gene has no curated RefSeq mRNA link.

    :param gene_id: NCBI Gene ID.
    :param species: Species name, used only for the fallback search.
    :param symbol: Gene symbol, used only for the fallback search.
    :return: List of dicts with accession, title, length, is_coding, is_refseq.
    """
    handle = Entrez.elink(dbfrom="gene", db="nuccore", linkname="gene_nuccore_refseqrna", id=gene_id)
    link_result = Entrez.read(handle)
    handle.close()

    linksets = link_result[0]["LinkSetDb"]
    uids = [link["Id"] for link in linksets[0]["Link"]] if linksets else []
    is_refseq = True

    if not uids and species and symbol:
        handle = Entrez.esearch(
            db="nuccore",
            term=f"{symbol}[gene] AND {species}[orgn] AND biomol_mrna[PROP]",
            retmax=20,
        )
        fallback_result = Entrez.read(handle)
        handle.close()
        uids = fallback_result.get("IdList", [])
        is_refseq = False

    if not uids:
        return []

    handle = Entrez.esummary(db="nuccore", id=",".join(uids))
    summaries = Entrez.read(handle)
    handle.close()

    results = []
    for doc in summaries:
        title = str(doc.get("Title", ""))
        results.append({
            "accession": str(doc.get("AccessionVersion", doc.get("Caption", ""))),
            "title": title,
            "length": int(doc.get("Length", 0)),
            "is_coding": "non-coding" not in title.lower(),
            "is_refseq": is_refseq,
        })
    return results


def fetch_transcript_fasta(accession):
    """
    Fetch the FASTA sequence for a nuccore accession.

    :param accession: NCBI nucleotide accession (e.g. "NM_002046.7").
    :return: FASTA text.
    """
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
    fasta_text = handle.read()
    handle.close()
    return fasta_text
