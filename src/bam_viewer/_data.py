"""Data access and layout helpers for the BAM viewer widget.

All functions here are pure (no widget / UI state) so they can be unit tested
on their own. Coordinates follow polars-bio's default convention: **1-based,
closed** intervals (``start`` and ``end`` are both inclusive), which also
matches the GTF convention so reads and annotations line up.
"""

from __future__ import annotations

import os
import re
from typing import Optional, Union

import polars as pl
import polars_bio as pb

GtfSource = Union[str, "os.PathLike", pl.DataFrame, pl.LazyFrame]

# CIGAR operations that consume the reference and so advance the genomic
# position.  ``N`` (skipped region, e.g. an intron) also consumes the
# reference but is rendered as a gap rather than an aligned block.
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_REF_BLOCK_OPS = frozenset("M=XD")


def cigar_blocks(start: int, cigar: Optional[str], end: int) -> list[list[int]]:
    """Split an alignment into aligned blocks, breaking on ``N`` (introns).

    ``start``/``end`` are 1-based inclusive.  Returns a list of
    ``[block_start, block_end]`` pairs (also 1-based inclusive).  Reads without
    usable CIGAR information collapse to a single ``[start, end]`` block.
    """
    if not cigar or cigar == "*":
        return [[start, end]]

    blocks: list[list[int]] = []
    pos = start
    cur_start = start
    cur_end = start - 1
    matched = False
    for length, op in _CIGAR_RE.findall(cigar):
        matched = True
        length = int(length)
        if op in _REF_BLOCK_OPS:
            cur_end = pos + length - 1
            pos += length
        elif op == "N":
            if cur_end >= cur_start:
                blocks.append([cur_start, cur_end])
            pos += length
            cur_start = pos
            cur_end = pos - 1
        # I, S, H, P do not consume the reference.
    if cur_end >= cur_start:
        blocks.append([cur_start, cur_end])
    if not matched or not blocks:
        return [[start, end]]
    return blocks


def pack_intervals(items: list[dict], gap: int) -> int:
    """Greedily assign each item a ``row`` so rows contain no overlaps.

    ``items`` must be sorted by ``start``.  Two items share a row only if a gap
    of at least ``gap`` base pairs separates them.  Mutates each item in place
    adding an integer ``row`` key and returns the number of rows used.
    """
    row_last_end: list[int] = []
    for it in items:
        placed = False
        for i, last_end in enumerate(row_last_end):
            if it["start"] > last_end + gap:
                row_last_end[i] = it["end"]
                it["row"] = i
                placed = True
                break
        if not placed:
            it["row"] = len(row_last_end)
            row_last_end.append(it["end"])
    return len(row_last_end)


def read_compatible(blocks: list[list[int]], exons: list[list[int]]) -> bool:
    """Is a read (its aligned ``blocks``) compatible with an isoform (``exons``)?

    Compatible means every aligned base lies within an exon of the isoform and
    every splice junction in the read matches an annotated junction (i.e. the
    read walks consecutive exons, splicing exactly at their boundaries).  The
    read may start/end partway into its first/last exon.  ``exons`` must be
    sorted and non-overlapping; both coordinate systems are 1-based inclusive.
    """
    n = len(blocks)
    if n == 0 or not exons:
        return False
    # Locate the exon containing the first block's start.
    b0 = blocks[0][0]
    j = next((i for i, (es, ee) in enumerate(exons) if es <= b0 <= ee), None)
    if j is None:
        return False
    for i, (bs, be) in enumerate(blocks):
        if j >= len(exons):
            return False  # read has more spliced segments than the isoform
        es, ee = exons[j]
        if bs < es or be > ee:
            return False  # aligned base outside this exon (in an intron / flank)
        if i > 0 and bs != es:
            return False  # spliced-in edge must meet the exon's start
        if i < n - 1 and be != ee:
            return False  # spliced-out edge must meet the exon's end
        j += 1
    return True


def query_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    *,
    max_reads: int = 5000,
    selected_exons: Optional[list[list[int]]] = None,
) -> dict:
    """Load and lay out reads overlapping ``chrom:start-end`` from a BAM file.

    Only the requested region is read from disk: polars-bio uses the ``.bai``
    index for predicate pushdown, so the cost scales with the window, not the
    file.  Returns a JSON-serialisable dict describing the reads and layout.
    """
    lf = pb.scan_bam(bam_path)
    df = lf.filter(
        (pl.col("chrom") == chrom)
        & (pl.col("end") >= start)
        & (pl.col("start") <= end)
    ).select(["start", "end", "flags", "cigar", "mapping_quality"])

    df = df.collect()
    total = df.height
    truncated = total > max_reads
    if truncated:
        # Protect the browser from pathological pileups by sampling down to a
        # representative subset rather than refusing outright.
        df = df.sample(n=max_reads, seed=0)

    reads: list[dict] = []
    for row in df.iter_rows(named=True):
        r_start = int(row["start"])
        r_end = int(row["end"])
        blocks = cigar_blocks(r_start, row["cigar"], r_end)
        read = {
            "start": r_start,
            "end": r_end,
            "strand": "-" if (int(row["flags"]) & 0x10) else "+",
            "mapq": int(row["mapping_quality"]),
        }
        # Only ship blocks when the read is actually spliced; otherwise the
        # frontend draws a single rectangle from start..end.
        if len(blocks) > 1:
            read["blocks"] = blocks
        if selected_exons is not None:
            read["compat"] = read_compatible(blocks, selected_exons)
        reads.append(read)

    reads.sort(key=lambda r: r["start"])
    gap = max(1, (end - start) // 300)
    n_rows = pack_intervals(reads, gap)

    return {
        "reads": reads,
        "n_rows": n_rows,
        "total": total,
        "shown": len(reads),
        "truncated": truncated,
    }


# GTF feature types we care about for a lightweight gene model.
_EXONIC = frozenset({"exon"})
_CDS = frozenset({"CDS"})
_SPAN_TYPES = frozenset({"transcript", "mRNA", "gene"})
# Sub-features that let us reconstruct exonic extent when a GTF has no explicit
# ``exon`` lines (common in CDS-centric annotations).
_UTR = frozenset({
    "five_prime_utr", "three_prime_utr", "UTR",
    "five_prime_UTR", "three_prime_UTR", "5UTR", "3UTR",
})
_CODON = frozenset({"start_codon", "stop_codon"})
_EXONIC_FALLBACK = _CDS | _UTR | _CODON


def _merge_blocks(blocks: list[list[int]]) -> list[list[int]]:
    """Sort and merge touching/overlapping ``[start, end]`` intervals."""
    merged: list[list[int]] = []
    for s, e in sorted(blocks):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged

_GTF_CORE = ("chrom", "start", "end", "type", "strand")
_GTF_ATTR_FIELDS = ("gene_id", "transcript_id", "gene_name")


def transcript_exons(gtf: pl.DataFrame, transcript_id: str) -> list[list[int]]:
    """Return the full, merged exon list of a transcript from a preloaded GTF.

    Uses the whole in-memory annotation (not a windowed view), so read
    compatibility can be judged against the complete isoform even when only
    part of it is on screen.  Exons are 1-based inclusive, sorted, and merged
    where they touch or overlap.
    """
    sub = gtf.filter(pl.col("transcript_id") == transcript_id).select(
        ["type", "start", "end"]
    )
    if sub.height == 0:
        return []
    exon, fallback = [], []
    for t, s, e in sub.iter_rows():
        if t in _EXONIC:
            exon.append([int(s), int(e)])
        elif t in _EXONIC_FALLBACK:
            # Reconstruct exon extent from CDS/UTR when a GTF omits exon lines.
            fallback.append([int(s), int(e)])
    return _merge_blocks(exon or fallback)


def _extract_attr(tag: str) -> pl.Expr:
    """Pull a single attribute value out of polars-bio's nested ``attributes``
    column (``List(Struct{tag, value})``)."""
    return (
        pl.col("attributes")
        .list.eval(
            pl.element()
            .struct.field("value")
            .filter(pl.element().struct.field("tag") == tag)
        )
        .list.first()
        .alias(tag)
    )


def load_gtf(source: GtfSource) -> pl.DataFrame:
    """Load a GTF annotation track fully into memory, ready for repeated
    region queries.

    ``source`` may be a path (``str`` / ``os.PathLike``), or an already-loaded
    polars ``DataFrame`` / ``LazyFrame`` — e.g. the result of
    ``polars_bio.read_gtf(...)``.  GTFs are small and usually unindexed, so the
    whole annotation is held in memory and filtered per region rather than
    re-read from disk on every pan/zoom.

    The returned frame is normalised to columns
    ``chrom, start, end, type, strand, gene_id, transcript_id, gene_name``.
    """
    if isinstance(source, (str, os.PathLike)):
        df = pb.read_gtf(os.fspath(source), attr_fields=list(_GTF_ATTR_FIELDS))
    elif isinstance(source, pl.LazyFrame):
        df = source.collect()
    elif isinstance(source, pl.DataFrame):
        df = source
    else:
        raise TypeError(
            "GTF track must be a path, polars DataFrame, or LazyFrame; got "
            f"{type(source).__name__}"
        )

    missing = [c for c in _GTF_CORE if c not in df.columns]
    if missing:
        raise ValueError(
            f"GTF frame is missing required column(s) {missing}. Expected the "
            f"polars-bio GTF schema with columns {_GTF_CORE} (load it with "
            "polars_bio.read_gtf)."
        )

    have_attrs = "attributes" in df.columns
    exprs = []
    for tag in _GTF_ATTR_FIELDS:
        if tag in df.columns:
            continue
        exprs.append(_extract_attr(tag) if have_attrs
                     else pl.lit(None, dtype=pl.Utf8).alias(tag))
    if exprs:
        df = df.with_columns(exprs)
    return df.select([*_GTF_CORE, *_GTF_ATTR_FIELDS])


def query_features(
    gtf: pl.DataFrame,
    chrom: str,
    start: int,
    end: int,
    *,
    max_features: int = 2000,
) -> dict:
    """Group GTF features overlapping ``chrom:start-end`` into gene models.

    ``gtf`` is a preloaded, normalised frame from :func:`load_gtf`.  Filtering
    is an in-memory operation, so panning/zooming never re-reads the file.
    Features are grouped by ``transcript_id`` into models carrying their exon
    (and CDS) blocks so the frontend can draw a familiar exon/intron track.
    """
    df = gtf.filter(
        (pl.col("chrom") == chrom)
        & (pl.col("end") >= start)
        & (pl.col("start") <= end)
    )

    # Group rows into transcript models keyed by transcript_id. Bare gene-span
    # rows are tracked separately and only kept when no transcript represents
    # that gene, so a GENCODE-style file shows its transcripts (with exons)
    # while a gene-only file still shows something.
    models: dict[str, dict] = {}
    gene_spans: dict[str, dict] = {}
    other: list[dict] = []

    for row in df.iter_rows(named=True):
        s, e = int(row["start"]), int(row["end"])
        ftype = row["type"]
        tid = row["transcript_id"]
        gid = row["gene_id"]
        label = row["gene_name"] or gid or tid or ftype
        strand = row["strand"] or "."

        if tid is None:
            span = {"start": s, "end": e, "strand": strand, "name": label,
                    "gene_id": gid, "transcript_id": None, "exons": [], "cds": []}
            if ftype in _SPAN_TYPES:
                # Collapse duplicate gene lines for the same gene_id.
                key = gid or label
                prev = gene_spans.get(key)
                if prev is None:
                    gene_spans[key] = span
                else:
                    prev["start"] = min(prev["start"], s)
                    prev["end"] = max(prev["end"], e)
            elif ftype not in _EXONIC | _CDS:
                other.append(span)
            continue

        m = models.get(tid)
        if m is None:
            m = models[tid] = {
                "start": s, "end": e, "strand": strand, "name": label,
                "gene_id": gid, "transcript_id": tid, "exons": [], "cds": [],
            }
        m["start"] = min(m["start"], s)
        m["end"] = max(m["end"], e)
        # Prefer a real gene_name; never let a bare gene_id clobber one.
        if row["gene_name"]:
            m["name"] = row["gene_name"]
        if gid:
            m["gene_id"] = gid
        if strand != ".":
            m["strand"] = strand
        if ftype in _EXONIC:
            m["exons"].append([s, e])
        elif ftype in _CDS:
            m["cds"].append([s, e])

    transcript_gene_ids = {m["gene_id"] for m in models.values()}
    orphan_genes = [
        g for g in gene_spans.values() if g["gene_id"] not in transcript_gene_ids
    ]
    features = list(models.values()) + orphan_genes + other
    for f in features:
        f["exons"].sort()
        f["cds"].sort()
    features.sort(key=lambda f: f["start"])

    total = len(features)
    truncated = total > max_features
    if truncated:
        features = features[:max_features]

    gap = max(1, (end - start) // 200)
    n_rows = pack_intervals(features, gap)

    return {
        "features": features,
        "n_rows": n_rows,
        "total": total,
        "truncated": truncated,
    }


_REGION_RE = re.compile(
    r"^\s*([^:\s]+)\s*(?::\s*([\d,]+)\s*-\s*([\d,]+))?\s*$"
)


def parse_region(region: str) -> tuple[str, Optional[int], Optional[int]]:
    """Parse ``chrom``, ``chrom:start-end`` (commas allowed) into a tuple.

    Returns ``(chrom, start, end)`` with ``start``/``end`` as ``None`` when the
    string is just a contig name.  Raises ``ValueError`` on malformed input.
    """
    m = _REGION_RE.match(region)
    if not m:
        raise ValueError(f"Could not parse region: {region!r}")
    chrom = m.group(1)
    if m.group(2) is None:
        return chrom, None, None
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", ""))
    if end < start:
        start, end = end, start
    return chrom, start, end
