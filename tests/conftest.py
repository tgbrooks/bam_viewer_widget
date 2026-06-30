import random

import pytest

pysam = pytest.importorskip("pysam")


@pytest.fixture(scope="session")
def bam_path(tmp_path_factory):
    """A small sorted+indexed BAM: 1500 reads on chr1, 500 on chr2."""
    d = tmp_path_factory.mktemp("data")
    path = d / "test.bam"
    random.seed(0)
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"LN": 100_000, "SN": "chr1"},
            {"LN": 50_000, "SN": "chr2"},
        ],
    }
    reads = []
    for i in range(2000):
        a = pysam.AlignedSegment()
        a.query_name = f"r{i}"
        rlen = 100
        a.query_sequence = "".join(random.choice("ACGT") for _ in range(rlen))
        a.flag = 0 if i % 2 == 0 else 16
        a.reference_id = 0 if i < 1500 else 1
        ref_len = 100_000 if i < 1500 else 50_000
        a.reference_start = random.randint(0, ref_len - rlen)
        a.mapping_quality = random.randint(0, 60)
        # Give a handful of reads a spliced (N) CIGAR to exercise block logic.
        if i % 500 == 0:
            a.cigar = [(0, 40), (3, 200), (0, 60)]
        else:
            a.cigar = [(0, rlen)]
        reads.append(a)
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(str(path))
    return str(path)


@pytest.fixture(scope="session")
def gtf_path(tmp_path_factory):
    d = tmp_path_factory.mktemp("data")
    path = d / "test.gtf"
    path.write_text(
        'chr1\ttest\tgene\t1000\t9000\t.\t+\t.\tgene_id "G1"; gene_name "GeneOne";\n'
        'chr1\ttest\ttranscript\t1000\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GeneOne";\n'
        'chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";\n'
        'chr1\ttest\texon\t5000\t6000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";\n'
        'chr1\ttest\tCDS\t5000\t6000\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\ttest\texon\t8000\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "3";\n'
    )
    return str(path)
