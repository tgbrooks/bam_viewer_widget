# bam_viewer_widget

A small, embeddable [marimo](https://marimo.io) widget for browsing aligned
reads from a **local** BAM file — no web server, no remote URL, no IGV. Built on
[anywidget](https://anywidget.dev) and
[polars-bio](https://github.com/biodatageeks/polars-bio).

It exists because `igv.js` needs a served URL and `igv-notebook` doesn't render
in marimo. This widget reads the BAM directly from disk and draws the reads
itself.

## Features

- **Lazy, indexed reads.** Only the visible window is read from disk, using the
  BAM `.bai` index via polars-bio's predicate pushdown. The cost scales with the
  size of the window, not the file.
- **Pan & zoom.** Drag to pan, scroll to zoom (or use the toolbar buttons), or
  type a region like `chr1:1,000-9,000`. Each move reloads just the new window.
- **GTF annotation tracks.** Add any number of gene/transcript tracks from GTF
  files (or preloaded polars DataFrames); exons, CDS, and strand are drawn as a
  familiar gene model. Annotations are loaded into memory once and filtered per
  view — GTFs are small and usually unindexed, so there's nothing to re-read.
- **Lightweight by design.** Refuses to render windows larger than
  `max_window` (default 100 kb) and samples dense pileups down to `max_reads`
  (default 5000) so the browser never chokes. There's no whole-chromosome view.
- **Spliced reads.** `N` CIGAR operations (introns) are drawn as gaps.

## Install

```bash
pip install bam-viewer-widget        # the widget
pip install "bam-viewer-widget[notebook]"   # also pulls in marimo
```

The BAM must be coordinate-sorted and indexed; the index is expected at
`<bam>.bai`:

```bash
samtools sort -o reads.sorted.bam reads.bam
samtools index reads.sorted.bam        # writes reads.sorted.bam.bai
```

## Usage

In a marimo cell:

```python
import marimo as mo
from bam_viewer import BamViewer

viewer = mo.ui.anywidget(
    BamViewer(
        "reads.sorted.bam",
        region="chr1:1,000-9,000",
        gtf_tracks={"GENCODE": "genes.gtf"},
    )
)
viewer
```

Read the current view back out from the notebook (reactively, like any marimo
UI element):

```python
viewer.region          # "chr1:1,000-9,000"
viewer.chrom, viewer.start, viewer.end
```

You can also pass an already-loaded annotation instead of a path — handy when
you want to load it once and share it across viewers, or filter/transform it
first:

```python
import polars_bio as pb

genes = pb.read_gtf("genes.gtf")          # whole annotation, in memory
viewer = mo.ui.anywidget(
    BamViewer("reads.sorted.bam", region="chr1:1,000-9,000",
              gtf_tracks={"GENCODE": genes})
)
```

Drive it from Python (attribute access is proxied to the widget):

```python
viewer.goto("chr2:500-2,500")
```

### Constructor

```python
BamViewer(
    bam_path,                 # path to a sorted, indexed BAM
    region="chr1:1000-9000",  # initial view; "chrom" alone starts at its 5' end
    gtf_tracks={"name": path},# optional annotation tracks
    max_window=100_000,       # largest window (bp) that will render
    max_reads=5000,           # pileups are sampled down to this many reads
)
```

## How it works

The Python side holds the view region as a single synced traitlet
(`[chrom, start, end]`, so each pan/zoom is one atomic reload). When the
frontend pans/zooms it writes the new region back; an observer reloads reads
(`polars_bio.scan_bam(...).filter(region)`, which uses the `.bai` index) and
filters the preloaded GTF frame(s) in memory, packs everything into
non-overlapping rows, and pushes the result to the canvas. Contig lengths are
read straight from the BAM header with the standard library (BGZF is
gzip-compatible), so **pysam is not a runtime dependency**.

## Development

```bash
pip install -e ".[dev]"
pytest
```

The test suite builds a small BAM with pysam and exercises the data layer and
the widget's reload logic.

## Limitations

- No base-level / mismatch / coverage view (kept intentionally lean).
- GTF tracks are held fully in memory (they're small and unindexed); BAM access
  is index-driven and reads only the visible window.
- No whole-chromosome overview by design — zoom in to a window of
  `max_window` bp or smaller.

## License

MIT
