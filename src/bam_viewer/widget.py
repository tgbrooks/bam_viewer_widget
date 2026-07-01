"""An embeddable anywidget for browsing a local BAM file in marimo.

Reads (and optional GTF annotation tracks) are loaded lazily with polars-bio:
only the visible window is read from disk via the ``.bai`` index, so the widget
stays responsive even on large files.  The frontend (``static/index.js``)
handles pan/zoom and asks Python to reload data whenever the view moves.
"""

from __future__ import annotations

import gzip
import pathlib
import struct
import traceback
from typing import Mapping, Optional, Union

import anywidget
import traitlets

from . import _data
from ._data import GtfSource

_STATIC = pathlib.Path(__file__).parent / "static"

_EMPTY_READS = {"reads": [], "n_rows": 0, "total": 0, "shown": 0, "truncated": False}


def _read_bam_contigs(path: str) -> dict[str, int]:
    """Return ``{contig_name: length}`` from a BAM header.

    BGZF is gzip-compatible, so the header parses with the stdlib alone — no
    pysam dependency.  Returns an empty dict if the header can't be read.
    """
    try:
        with gzip.open(path, "rb") as f:
            if f.read(4) != b"BAM\x01":
                return {}
            (l_text,) = struct.unpack("<i", f.read(4))
            f.read(l_text)
            (n_ref,) = struct.unpack("<i", f.read(4))
            contigs: dict[str, int] = {}
            for _ in range(n_ref):
                (l_name,) = struct.unpack("<i", f.read(4))
                name = f.read(l_name)[:-1].decode()
                (l_ref,) = struct.unpack("<i", f.read(4))
                contigs[name] = l_ref
            return contigs
    except (OSError, struct.error, UnicodeDecodeError):
        return {}


class BamViewer(anywidget.AnyWidget):
    """Display aligned reads from an indexed BAM file with GTF annotations.

    Parameters
    ----------
    bam_path:
        Path to a coordinate-sorted, indexed BAM file. The index is expected at
        ``<bam_path>.bai``.
    region:
        Initial view as ``"chrom:start-end"`` (1-based, commas allowed) or just
        ``"chrom"`` to start at the contig's beginning.
    gtf_tracks:
        Optional mapping of ``{track_name: source}`` rendered as gene tracks.
        Each ``source`` is a GTF path or an already-loaded polars
        ``DataFrame`` / ``LazyFrame`` (e.g. ``polars_bio.read_gtf(...)``). The
        whole annotation is preloaded into memory and filtered per view, since
        GTFs are small and usually unindexed.
    max_window:
        Largest window (in bp) that will be rendered. Zooming out past this
        shows a "zoom in" message instead of loading data.
    max_reads:
        Reads are sampled down to this many before being sent to the browser.

    Wrap the instance in ``marimo.ui.anywidget(...)`` to embed it in a notebook.
    """

    _esm = _STATIC / "index.js"
    _css = _STATIC / "style.css"

    # --- configuration (set once from Python) ---------------------------------
    bam_path = traitlets.Unicode().tag(sync=True)
    max_window = traitlets.Int(300_000).tag(sync=True)
    max_annotation_window = traitlets.Int(1_000_000).tag(sync=True)
    max_reads = traitlets.Int(5000).tag(sync=True)
    track_names = traitlets.List(traitlets.Unicode()).tag(sync=True)
    contigs = traitlets.Dict().tag(sync=True)

    # --- view region: a single [chrom, start, end] trait so a pan/zoom is one
    #     atomic update (and therefore exactly one reload), not three. ----------
    _view = traitlets.List().tag(sync=True)

    # --- selected isoform: {track, transcript_id} (empty = none). Its full exon
    #     list rides along inside _read_data (see _reload) so it reaches the
    #     frontend through the same reliably-synced channel as the reads; the
    #     frontend judges each read's compatibility against those exons. --------
    _selected = traitlets.Dict().tag(sync=True)

    # --- data pushed to the frontend (both version-stamped; see _publish_*) ---
    _read_data = traitlets.Dict().tag(sync=True)
    _feature_data = traitlets.Dict().tag(sync=True)
    _message = traitlets.Unicode("").tag(sync=True)
    _loading = traitlets.Bool(False).tag(sync=True)

    def __init__(
        self,
        bam_path: Union[str, pathlib.Path],
        region: Optional[str] = None,
        gtf_tracks: Optional[Mapping[str, GtfSource]] = None,
        *,
        max_window: int = 300_000,
        max_annotation_window: int = 1_000_000,
        max_reads: int = 5000,
        **kwargs,
    ):
        bam_path = str(bam_path)
        if not pathlib.Path(bam_path).exists():
            raise FileNotFoundError(bam_path)
        index = pathlib.Path(bam_path + ".bai")
        if not index.exists():
            raise FileNotFoundError(
                f"BAM index not found: {index}. The BAM must be sorted and "
                f"indexed (e.g. `samtools index {bam_path}`)."
            )

        # Preload each GTF track once into an in-memory frame; per-view
        # filtering is then just an in-memory operation.
        self._gtf_frames = {
            name: _data.load_gtf(src) for name, src in (gtf_tracks or {}).items()
        }
        contigs = _read_bam_contigs(bam_path)
        chrom, start, end = self._initial_region(region, contigs)

        # Suppress the reload while the initial traits are assigned; we do a
        # single explicit reload at the end of __init__ instead.
        self._ver = 0
        self._suspend_reload = True
        super().__init__(
            bam_path=bam_path,
            max_window=max_window,
            max_annotation_window=max_annotation_window,
            max_reads=max_reads,
            track_names=list(self._gtf_frames),
            contigs=contigs,
            _view=[chrom, start, end],
            **kwargs,
        )
        self._suspend_reload = False
        self._reload()

    @staticmethod
    def _initial_region(
        region: Optional[str], contigs: dict[str, int]
    ) -> tuple[str, int, int]:
        if region:
            chrom, start, end = _data.parse_region(region)
        elif contigs:
            chrom, start, end = next(iter(contigs)), None, None
        else:
            raise ValueError(
                "No region given and contigs could not be read from the BAM "
                "header; pass region='chrom:start-end'."
            )
        if start is None or end is None:
            length = contigs.get(chrom, 10_000)
            start, end = 1, min(length, 10_000)
        return chrom, int(start), int(end)

    # --- reactivity -----------------------------------------------------------
    @traitlets.observe("_view")
    def _on_view_change(self, _change):
        if not getattr(self, "_suspend_reload", False):
            self._reload()

    @traitlets.observe("_selected")
    def _on_selection_change(self, _change):
        # Re-publish the current reads with the newly selected isoform's exons
        # attached (no BAM re-query — the reads are unchanged).
        rd = dict(self._read_data)
        rd["selected_exons"] = self._compute_selected_exons()
        self._publish_reads(rd)

    def _compute_selected_exons(self):
        """Full exon list of the currently selected isoform ([] if none)."""
        sel = self._selected or {}
        track, tid = sel.get("track"), sel.get("transcript_id")
        df = self._gtf_frames.get(track)
        if df is None or not tid:
            return []
        return _data.transcript_exons(df, tid)

    def _next_version(self) -> int:
        v = getattr(self, "_ver", 0) + 1
        self._ver = v
        return v

    def _publish_reads(self, rd: dict) -> None:
        # Stamp a monotonic version so the frontend can drop stale re-syncs.
        # (Works around a marimo 0.23.12+ regression where a trait Python
        # updates in response to a frontend change is clobbered by a stale echo.)
        rd = dict(rd)
        rd["_v"] = self._next_version()
        self._read_data = rd

    def _publish_features(self, tracks: list) -> None:
        self._feature_data = {"_v": self._next_version(), "tracks": tracks}

    def _reload(self):
        """Load reads + features for the current window and push to frontend.

        Annotations render for windows up to ``max_annotation_window``; reads,
        which are far heavier, only up to the smaller ``max_window``.
        """
        chrom, start, end = self._view
        start, end = int(start), int(end)
        if end < start:
            start, end = end, start
        width = end - start + 1
        exons = self._compute_selected_exons()  # rides along in _read_data

        if width > self.max_annotation_window:
            self._message = (
                f"Window too large ({width:,} bp). Zoom in to "
                f"≤ {self.max_annotation_window:,} bp."
            )
            self._publish_reads({**_EMPTY_READS, "selected_exons": exons})
            self._publish_features([])
            return

        self._loading = True
        try:
            self._publish_features([
                {"name": name, **_data.query_features(df, chrom, start, end)}
                for name, df in self._gtf_frames.items()
            ])
            if width > self.max_window:
                # Too wide for reads, but annotations still render.
                rd = {
                    **_EMPTY_READS,
                    "note": f"Zoom in to ≤ {self.max_window:,} bp to see reads.",
                }
            else:
                rd = _data.query_reads(
                    self.bam_path, chrom, start, end, max_reads=self.max_reads
                )
            rd["selected_exons"] = exons
            self._publish_reads(rd)
            self._message = ""
        except Exception:  # surface load errors in the widget, don't crash
            self._message = "Error loading region:\n" + traceback.format_exc(limit=2)
            self._publish_reads({**_EMPTY_READS, "selected_exons": exons})
            self._publish_features([])
        finally:
            self._loading = False

    # --- convenience API ------------------------------------------------------
    @property
    def chrom(self) -> str:
        return self._view[0]

    @property
    def start(self) -> int:
        return int(self._view[1])

    @property
    def end(self) -> int:
        return int(self._view[2])

    @property
    def region(self) -> str:
        """The current view as a ``chrom:start-end`` string."""
        return f"{self.chrom}:{self.start:,}-{self.end:,}"

    def goto(self, region: str) -> None:
        """Jump the view to ``region`` (``"chrom:start-end"`` or ``"chrom"``)."""
        chrom, start, end = self._initial_region(region, dict(self.contigs))
        self._view = [chrom, start, end]  # fires the observer -> reload

    @property
    def selected_transcript(self) -> Optional[str]:
        """The transcript_id of the currently selected isoform, or ``None``."""
        return (self._selected or {}).get("transcript_id")

    def select_transcript(self, transcript_id: str, track: Optional[str] = None) -> None:
        """Select an isoform so reads incompatible with it fade out.

        ``track`` defaults to the first GTF track. Pass ``None`` transcript to
        :meth:`clear_selection`.
        """
        if track is None:
            track = next(iter(self._gtf_frames), None)
        self._selected = {"track": track, "transcript_id": transcript_id}

    def clear_selection(self) -> None:
        """Clear any selected isoform (all reads return to full opacity)."""
        self._selected = {}
