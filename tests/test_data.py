import os

import polars_bio as pb

from bam_viewer import _data
from bam_viewer.widget import BamViewer, _read_bam_contigs


def test_cigar_blocks_simple():
    assert _data.cigar_blocks(100, "100M", 199) == [[100, 199]]
    assert _data.cigar_blocks(100, None, 150) == [[100, 150]]
    assert _data.cigar_blocks(100, "*", 150) == [[100, 150]]


def test_cigar_blocks_spliced():
    # 40M 200N 60M starting at 100 -> two blocks, gap skipped.
    blocks = _data.cigar_blocks(100, "40M200N60M", 399)
    assert blocks == [[100, 139], [340, 399]]


def test_cigar_blocks_softclip_and_deletion():
    # Soft clips don't consume reference; deletions stay within the block.
    blocks = _data.cigar_blocks(100, "10S50M5D50M10S", 204)
    assert blocks == [[100, 204]]


EXONS = [[1000, 2000], [5000, 6000], [8000, 9000]]  # a 3-exon isoform


def test_read_compatible_single_exon_read():
    # Unspliced read fully inside one exon -> compatible.
    assert _data.read_compatible([[1200, 1800]], EXONS)
    # Unspliced read inside an intron -> incompatible.
    assert not _data.read_compatible([[3000, 3500]], EXONS)
    # Unspliced read straddling an exon/intron boundary -> incompatible.
    assert not _data.read_compatible([[1800, 2200]], EXONS)


def test_read_compatible_spliced():
    # Junction exactly matching annotated exon1->exon2 boundary -> compatible.
    assert _data.read_compatible([[1500, 2000], [5000, 5500]], EXONS)
    # Spanning three exons, matching both junctions -> compatible.
    assert _data.read_compatible([[1900, 2000], [5000, 6000], [8000, 8100]], EXONS)


def test_read_compatible_wrong_junction():
    # Right exons but the splice donor is off by a base -> incompatible.
    assert not _data.read_compatible([[1500, 1999], [5000, 5500]], EXONS)
    # Skips exon2 (junction exon1->exon3 is not annotated) -> incompatible.
    assert not _data.read_compatible([[1500, 2000], [8000, 8500]], EXONS)
    # An interior block must fill its exon exactly.
    assert not _data.read_compatible([[1900, 2000], [5000, 5500], [8000, 8100]], EXONS)


def test_read_compatible_beyond_transcript():
    # Extends 5' of the first exon -> some aligned bases outside the isoform.
    assert not _data.read_compatible([[900, 1500]], EXONS)
    assert not _data.read_compatible([[1500, 2000], [5000, 6000], [8000, 9500]], EXONS)


def test_transcript_exons_merges(gtf_path):
    gtf = _data.load_gtf(gtf_path)
    exons = _data.transcript_exons(gtf, "T1")
    assert exons == [[1000, 2000], [5000, 6000], [8000, 9000]]
    assert _data.transcript_exons(gtf, "nope") == []


def test_query_reads_with_selection(bam_path, gtf_path):
    exons = _data.transcript_exons(_data.load_gtf(gtf_path), "T1")
    out = _data.query_reads(bam_path, "chr1", 1000, 9000, selected_exons=exons)
    assert all("compat" in r for r in out["reads"])
    # Some reads land inside exons (compatible), some don't.
    assert any(r["compat"] for r in out["reads"])
    for r in out["reads"]:
        assert isinstance(r["compat"], bool)


def test_widget_selection_publishes_exons(bam_path, gtf_path):
    # Selecting an isoform publishes its full exon list (the frontend judges
    # compatibility from that); it must not re-query or alter the reads.
    w = BamViewer(bam_path, region="chr1:1,000-9,000", gtf_tracks={"genes": gtf_path})
    assert w._selected_exons == []
    reads_before = w._read_data["reads"]
    w.select_transcript("T1")
    assert w.selected_transcript == "T1"
    assert w._selected_exons == [[1000, 2000], [5000, 6000], [8000, 9000]]
    assert w._read_data["reads"] is reads_before  # reads untouched, no reload
    w.clear_selection()
    assert w.selected_transcript is None
    assert w._selected_exons == []


def test_transcript_exons_cds_only(tmp_path):
    # A GTF describing a transcript with CDS/UTR but no exon lines still yields
    # a usable exon model (reconstructed from CDS + UTR, merged).
    p = tmp_path / "cds.gtf"
    p.write_text(
        'chr1\tt\ttranscript\t1000\t7000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\tfive_prime_utr\t1000\t1099\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\tCDS\t1100\t1200\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\tCDS\t5000\t6900\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tt\tthree_prime_utr\t6901\t7000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    )
    exons = _data.transcript_exons(_data.load_gtf(str(p)), "T1")
    assert exons == [[1000, 1200], [5000, 7000]]


def test_parse_region():
    assert _data.parse_region("chr1") == ("chr1", None, None)
    assert _data.parse_region("chr1:1,000-9,000") == ("chr1", 1000, 9000)
    # Reversed bounds get normalised.
    assert _data.parse_region("chrX:50-10") == ("chrX", 10, 50)


def test_pack_intervals_no_overlap():
    items = [
        {"start": 1, "end": 10},
        {"start": 5, "end": 15},  # overlaps first -> new row
        {"start": 100, "end": 110},  # far from first -> shares row 0
    ]
    n = _data.pack_intervals(items, gap=1)
    assert n == 2
    assert items[0]["row"] == 0
    assert items[1]["row"] == 1
    assert items[2]["row"] == 0


def test_query_reads_region(bam_path):
    out = _data.query_reads(bam_path, "chr1", 20_000, 21_000, max_reads=5000)
    assert out["total"] == out["shown"]
    assert out["total"] > 0
    assert not out["truncated"]
    for r in out["reads"]:
        # Every returned read overlaps the requested window.
        assert r["end"] >= 20_000 and r["start"] <= 21_000
        assert r["strand"] in ("+", "-")
        assert 0 <= r["row"] < out["n_rows"]


def test_query_reads_downsample(bam_path):
    out = _data.query_reads(bam_path, "chr1", 1, 100_000, max_reads=50)
    assert out["truncated"]
    assert out["shown"] == 50
    assert out["total"] == 1500


def test_query_reads_empty(bam_path):
    out = _data.query_reads(bam_path, "chr2", 1, 100, max_reads=5000)
    # chr2 reads exist but unlikely at the very start; assert structure regardless.
    assert "reads" in out and out["n_rows"] >= 0


def test_query_features(gtf_path):
    gtf = _data.load_gtf(gtf_path)
    out = _data.query_features(gtf, "chr1", 1, 10_000)
    feats = out["features"]
    transcripts = [f for f in feats if f["exons"]]
    assert transcripts, "expected a transcript model with exons"
    t = transcripts[0]
    assert t["name"] == "GeneOne"
    assert t["strand"] == "+"
    assert t["transcript_id"] == "T1"
    assert [1000, 2000] in t["exons"]
    assert t["cds"] == [[5000, 6000]]


def test_query_features_drops_redundant_gene_span(gtf_path):
    # The GTF has both a `gene` and a `transcript` line for G1; only the
    # transcript model (with exons) should remain — not a duplicate gene bar.
    out = _data.query_features(_data.load_gtf(gtf_path), "chr1", 1, 10_000)
    assert len(out["features"]) == 1
    assert out["features"][0]["exons"]


def test_query_features_gene_only_gtf(tmp_path):
    # A file with only gene lines (no transcripts/exons) still yields spans.
    p = tmp_path / "genes_only.gtf"
    p.write_text(
        'chr1\tt\tgene\t100\t500\t.\t+\t.\tgene_id "GA"; gene_name "Alpha";\n'
        'chr1\tt\tgene\t100\t500\t.\t+\t.\tgene_id "GA"; gene_name "Alpha";\n'  # dup
        'chr1\tt\tgene\t800\t900\t.\t-\t.\tgene_id "GB"; gene_name "Beta";\n'
    )
    out = _data.query_features(_data.load_gtf(str(p)), "chr1", 1, 1000)
    names = sorted(f["name"] for f in out["features"])
    assert names == ["Alpha", "Beta"]  # deduped
    assert all(f["exons"] == [] for f in out["features"])


def test_load_gtf_accepts_raw_dataframe(gtf_path):
    # Passing the nested-attributes DataFrame from read_gtf should work: the
    # gene_id/transcript_id/gene_name are extracted from `attributes`.
    raw = pb.read_gtf(gtf_path)  # no attr_fields -> nested `attributes` column
    assert "attributes" in raw.columns
    gtf = _data.load_gtf(raw)
    assert {"gene_id", "transcript_id", "gene_name"} <= set(gtf.columns)
    out = _data.query_features(gtf, "chr1", 1, 10_000)
    assert out["features"][0]["name"] == "GeneOne"


def test_load_gtf_preloads_once(gtf_path):
    # query_features must not touch the file: filtering a preloaded frame works
    # even if the original path is gone.
    import shutil
    tmp = gtf_path + ".copy.gtf"
    shutil.copy(gtf_path, tmp)
    gtf = _data.load_gtf(tmp)
    os.remove(tmp)
    out = _data.query_features(gtf, "chr1", 1, 10_000)
    assert out["features"]


def test_query_features_region_filter(gtf_path):
    # A window past the gene should return nothing.
    out = _data.query_features(_data.load_gtf(gtf_path), "chr1", 20_000, 30_000)
    assert out["features"] == []


def test_read_bam_contigs(bam_path):
    contigs = _read_bam_contigs(bam_path)
    assert contigs == {"chr1": 100_000, "chr2": 50_000}


def test_widget_construction_and_reload(bam_path, gtf_path):
    w = BamViewer(bam_path, region="chr1:20,000-21,000", gtf_tracks={"genes": gtf_path})
    assert w.chrom == "chr1"
    assert w.start == 20_000 and w.end == 21_000
    assert w._read_data["shown"] > 0
    assert len(w._feature_data) == 1
    assert w._message == ""

    # Changing the region triggers an automatic reload via the observer.
    w._view = ["chr1", 1, 9000]
    assert w._read_data["shown"] > 0
    assert w._feature_data[0]["features"], "GTF features should load in gene region"


def test_widget_reads_window_vs_annotation_window(bam_path, gtf_path):
    # Between max_window and max_annotation_window: annotations load, reads are
    # withheld with a "zoom in" note (no global error message).
    w = BamViewer(
        bam_path, region="chr1:1-30,000", gtf_tracks={"genes": gtf_path},
        max_window=10_000, max_annotation_window=50_000,
    )
    assert w._message == ""
    assert w._read_data["reads"] == []
    assert "note" in w._read_data and "zoom in" in w._read_data["note"].lower()
    assert w._feature_data[0]["features"]  # annotations still present

    # Beyond max_annotation_window: nothing renders, a message explains why.
    w.goto("chr1:1-100,000")
    assert "too large" in w._message.lower()
    assert w._read_data["reads"] == []
    assert w._feature_data == []


def test_widget_goto(bam_path):
    w = BamViewer(bam_path, region="chr1:1-1000")
    w.goto("chr2:1,000-5,000")
    assert w.chrom == "chr2"
    assert w.start == 1000 and w.end == 5000


def test_widget_requires_index(tmp_path):
    import shutil
    fake = tmp_path / "noindex.bam"
    fake.write_bytes(b"\x1f\x8b")  # not a valid bam, but exists
    try:
        BamViewer(str(fake), region="chr1:1-100")
    except FileNotFoundError as e:
        assert "index" in str(e).lower()
    else:
        raise AssertionError("expected FileNotFoundError for missing .bai")
