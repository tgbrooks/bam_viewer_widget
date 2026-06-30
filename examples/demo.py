import marimo

__generated_with = "0.23.11"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from bam_viewer import BamViewer

    return BamViewer, mo


@app.cell
def _(mo):
    mo.md(
        """
        # BAM viewer demo

        First generate the demo files (once): `python examples/make_demo_data.py`,
        then run this notebook from the same directory.

        Or point the paths below at your own coordinate-sorted, indexed BAM
        (`<bam>.bai` alongside it) and any GTF annotation files.

        Drag to pan, scroll to zoom, or type a region in the box.
        """
    )
    return


@app.cell
def _():
    # Edit these to your own files.
    BAM_PATH = "demo.bam"
    GTF_TRACKS = {"genes": "demo.gtf"}  # or {} for none
    REGION = "chr1:1,000-9,000"
    return BAM_PATH, GTF_TRACKS, REGION


@app.cell
def _(BAM_PATH, BamViewer, GTF_TRACKS, REGION, mo):
    viewer = mo.ui.anywidget(
        BamViewer(
            BAM_PATH,
            region=REGION,
            gtf_tracks=GTF_TRACKS,
            max_window=100_000,
            max_reads=5000,
        )
    )
    viewer
    return (viewer,)


@app.cell
def _(mo, viewer):
    # The current view is reactive: this cell re-runs whenever you pan/zoom.
    mo.md(f"**Current view:** `{viewer.region}`")
    return


if __name__ == "__main__":
    app.run()
