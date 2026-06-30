"""Generate a tiny sorted+indexed BAM and a GTF for the demo notebook.

    python examples/make_demo_data.py

Writes ``demo.bam`` (+ ``demo.bam.bai``) and ``demo.gtf`` to the current
directory. Requires pysam (``pip install pysam``).
"""

import random

import pysam


def main() -> None:
    random.seed(0)
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": 100_000, "SN": "chr1"}, {"LN": 50_000, "SN": "chr2"}],
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
        # A few spliced reads to show intron gaps.
        a.cigar = [(0, 40), (3, 200), (0, 60)] if i % 500 == 0 else [(0, rlen)]
        reads.append(a)
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile("demo.bam", "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index("demo.bam")

    with open("demo.gtf", "w") as f:
        f.write(
            'chr1\tdemo\tgene\t1000\t9000\t.\t+\t.\tgene_id "G1"; gene_name "GeneOne";\n'
            'chr1\tdemo\ttranscript\t1000\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GeneOne";\n'
            'chr1\tdemo\texon\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tdemo\texon\t5000\t6000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tdemo\tCDS\t5000\t6000\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tdemo\texon\t8000\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        )
    print("Wrote demo.bam, demo.bam.bai, demo.gtf")


if __name__ == "__main__":
    main()
