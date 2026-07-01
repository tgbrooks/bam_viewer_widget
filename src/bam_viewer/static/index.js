// Frontend for the BamViewer anywidget.
//
// Rendering is done on a single <canvas>; reads and annotation features carry
// genomic coordinates and are mapped to pixels against the *local* view, so
// panning/zooming feels instant. Whenever the view settles we push the new
// region to Python (debounced), which reloads only the visible window and sends
// fresh data back via the `_read_data` / `_feature_data` traits.

const RULER_H = 30;
const TRACK_HEADER_H = 15;
const FEAT_ROW_H = 16;
const READ_ROW_H = 9;
const SECTION_GAP = 6;
const VIEW_H = 480;
const MIN_WIDTH = 40; // smallest window in bp
const COMMIT_DELAY = 160; // ms debounce before asking Python to reload

const COLORS = {
  plus: "#5a7fb0",
  minus: "#b06a6a",
  intron: "#9aa0a6",
  exon: "#3a7d44",
  cds: "#2f6e3a",
  axis: "#444",
  tick: "#888",
  grid: "#f0f0f0",
  selected: "#e8a23d",
  selectedBg: "rgba(232, 162, 61, 0.18)",
};

const FADED_ALPHA = 0.12; // opacity of reads incompatible with the selection

function clamp(v, lo, hi) {
  return Math.max(lo, Math.min(hi, v));
}

// Choose a "nice" tick step (1/2/5 * 10^n) giving roughly `target` ticks.
function niceStep(span, target) {
  const raw = span / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = norm >= 5 ? 5 : norm >= 2 ? 2 : 1;
  return step * mag;
}

function fmtBp(n) {
  return Math.round(n).toLocaleString("en-US");
}

// Mirror of _data.read_compatible: are a read's aligned blocks compatible with
// an isoform's exon list? Every base must fall inside an exon and every splice
// junction must match an annotated one. Computed on the frontend so selecting
// an isoform never re-queries reads.
function readCompatible(blocks, exons) {
  const n = blocks.length;
  if (n === 0 || !exons || !exons.length) return false;
  const b0 = blocks[0][0];
  let j = -1;
  for (let i = 0; i < exons.length; i++) {
    if (exons[i][0] <= b0 && b0 <= exons[i][1]) {
      j = i;
      break;
    }
  }
  if (j === -1) return false;
  for (let i = 0; i < n; i++) {
    if (j >= exons.length) return false;
    const bs = blocks[i][0];
    const be = blocks[i][1];
    const es = exons[j][0];
    const ee = exons[j][1];
    if (bs < es || be > ee) return false;
    if (i > 0 && bs !== es) return false;
    if (i < n - 1 && be !== ee) return false;
    j++;
  }
  return true;
}

function render({ model, el }) {
  el.classList.add("bamv");
  el.innerHTML = `
    <div class="bamv-toolbar">
      <input class="bamv-region" type="text" spellcheck="false"
             title="Type a region, e.g. chr1:1,000-9,000" />
      <button class="bamv-btn bamv-go">Go</button>
      <span class="bamv-spacer"></span>
      <button class="bamv-btn bamv-zoomout" title="Zoom out">&minus;</button>
      <button class="bamv-btn bamv-zoomin" title="Zoom in">+</button>
      <button class="bamv-btn bamv-left" title="Pan left">&larr;</button>
      <button class="bamv-btn bamv-right" title="Pan right">&rarr;</button>
    </div>
    <div class="bamv-canvas-wrap">
      <canvas class="bamv-canvas"></canvas>
      <div class="bamv-overlay"></div>
      <div class="bamv-tooltip" style="display:none"></div>
    </div>
    <div class="bamv-status"></div>
  `;

  const regionInput = el.querySelector(".bamv-region");
  const canvas = el.querySelector(".bamv-canvas");
  const overlay = el.querySelector(".bamv-overlay");
  const tooltip = el.querySelector(".bamv-tooltip");
  const status = el.querySelector(".bamv-status");
  const ctx = canvas.getContext("2d");

  // Local, immediately-mutable view. Mirrors the model's `_view` ([chrom,
  // start, end]) but is updated optimistically during interaction.
  const initView = model.get("_view") || ["", 1, 1000];
  const view = { chrom: initView[0], start: initView[1], end: initView[2] };
  let hitReads = []; // drawn read rects for tooltip hit-testing
  let hitFeats = [];
  let commitTimer = null;
  // True while we are writing the region back to the model, so the resulting
  // change:* events don't re-enter and clobber our own view.
  let selfUpdating = false;
  let contentScrollY = 0; // vertical scroll offset of the content band
  let scrollbar = null; // geometry of the content scrollbar (or null)

  // Version-guarded snapshots of the data traits. marimo (>=0.23.12) can
  // re-send a *stale* value of a trait after Python has already updated it in
  // response to a frontend change; we keep the newest version we've seen and
  // ignore any older re-sync so reads/annotations never revert.
  let readState = model.get("_read_data") || { reads: [], n_rows: 0 };
  let featureState = (model.get("_feature_data") || {}).tracks || [];
  let lastReadV = readState._v ?? -1;
  let lastFeatV = (model.get("_feature_data") || {})._v ?? -1;

  const plotLeft = () => 8;
  const plotWidth = () => Math.max(10, canvas.clientWidth - 16);
  const width = () => view.end - view.start + 1;
  const bpToPx = (bp) =>
    plotLeft() + ((bp - view.start) / width()) * plotWidth();
  const pxToBp = (px) =>
    view.start + ((px - plotLeft()) / plotWidth()) * width();
  // Widest window we allow zooming out to: where even annotations stop drawing.
  const maxWidth = () =>
    Math.max(MIN_WIDTH, model.get("max_annotation_window") || 50 * model.get("max_window"));

  function contigLen() {
    const c = model.get("contigs") || {};
    return c[view.chrom] || Infinity;
  }

  // Does this feature match the current selection (for highlighting)?
  function isSelected(trackName, f) {
    const s = model.get("_selected") || {};
    return (
      !!s.transcript_id &&
      f.transcript_id === s.transcript_id &&
      s.track === trackName
    );
  }

  // Snap the view to whole-bp integers and keep it inside the contig.
  function clampView() {
    const w = clamp(Math.round(width()), MIN_WIDTH, maxWidth());
    let s = Math.round(view.start);
    const maxLen = contigLen();
    if (isFinite(maxLen)) s = clamp(s, 1, Math.max(1, maxLen - w + 1));
    else s = Math.max(1, s);
    view.start = s;
    view.end = s + w - 1;
  }

  function scheduleCommit() {
    if (commitTimer) clearTimeout(commitTimer);
    commitTimer = setTimeout(commit, COMMIT_DELAY);
  }

  function commit() {
    if (commitTimer) clearTimeout(commitTimer);
    commitTimer = null;
    clampView();
    const cur = model.get("_view") || [];
    // Nothing to do if the model already matches the view.
    if (cur[0] === view.chrom && cur[1] === view.start && cur[2] === view.end) {
      return;
    }
    selfUpdating = true;
    model.set("_view", [view.chrom, view.start, view.end]);
    model.save_changes();
    selfUpdating = false;
  }

  // ---- drawing -------------------------------------------------------------
  function setupCanvas() {
    const dpr = window.devicePixelRatio || 1;
    const w = el.querySelector(".bamv-canvas-wrap").clientWidth;
    canvas.style.width = w + "px";
    canvas.style.height = VIEW_H + "px";
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(VIEW_H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function drawRuler() {
    const w = width();
    const pxPerBp = plotWidth() / w;
    const step = niceStep(w, 8);
    ctx.fillStyle = COLORS.tick;
    ctx.strokeStyle = COLORS.grid;
    ctx.font = "10px system-ui, sans-serif";
    ctx.textBaseline = "alphabetic";
    const first = Math.ceil(view.start / step) * step;
    for (let bp = first; bp <= view.end; bp += step) {
      const x = bpToPx(bp);
      ctx.strokeStyle = COLORS.grid;
      ctx.beginPath();
      ctx.moveTo(x, RULER_H);
      ctx.lineTo(x, VIEW_H);
      ctx.stroke();
      ctx.strokeStyle = COLORS.tick;
      ctx.beginPath();
      ctx.moveTo(x, RULER_H - 6);
      ctx.lineTo(x, RULER_H);
      ctx.stroke();
      ctx.fillText(fmtBp(bp), x + 2, RULER_H - 8);
    }
    ctx.fillStyle = COLORS.axis;
    ctx.font = "11px system-ui, sans-serif";
    const label = `${view.chrom}:${fmtBp(view.start)}-${fmtBp(view.end)}  (${fmtBp(w)} bp, ${pxPerBp.toFixed(2)} px/bp)`;
    ctx.fillText(label, plotLeft(), 12);
  }

  function drawFeature(f, top, trackName) {
    const y = top + f.row * FEAT_ROW_H;
    const cy = y + FEAT_ROW_H / 2;
    const x0 = clamp(bpToPx(f.start), -5, canvas.clientWidth + 5);
    const x1 = clamp(bpToPx(f.end + 1), -5, canvas.clientWidth + 5);
    const selected = isSelected(trackName, f);
    if (selected) {
      // Highlight band behind the selected isoform's whole row.
      ctx.fillStyle = COLORS.selectedBg;
      ctx.fillRect(plotLeft(), y, plotWidth(), FEAT_ROW_H - 1);
    }
    // Intron / backbone line.
    ctx.strokeStyle = selected ? COLORS.selected : COLORS.intron;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(x0, cy);
    ctx.lineTo(x1, cy);
    ctx.stroke();
    // Strand arrows along the backbone.
    if (x1 - x0 > 20) {
      const dir = f.strand === "-" ? -1 : 1;
      for (let x = x0 + 8; x < x1 - 4; x += 18) {
        ctx.beginPath();
        ctx.moveTo(x, cy - 3);
        ctx.lineTo(x + 3 * dir, cy);
        ctx.lineTo(x, cy + 3);
        ctx.stroke();
      }
    }
    ctx.fillStyle = selected ? COLORS.selected : COLORS.exon;
    if (f.exons && f.exons.length) {
      for (const [es, ee] of f.exons) {
        const ex0 = bpToPx(es);
        const ew = Math.max(1, bpToPx(ee + 1) - ex0);
        ctx.fillRect(ex0, y + 3, ew, FEAT_ROW_H - 6);
      }
    } else {
      // Exon-less span (a bare gene line): a slim bar, not a fake exon.
      ctx.fillRect(x0, cy - 1.5, Math.max(1, x1 - x0), 3);
    }
    ctx.fillStyle = COLORS.cds;
    for (const [cs, ce] of f.cds || []) {
      const cx0 = bpToPx(cs);
      const cw = Math.max(1, bpToPx(ce + 1) - cx0);
      ctx.fillRect(cx0, y + 1, cw, FEAT_ROW_H - 2);
    }
    // Label, placed just left of the feature if there's room, else inside.
    if (f.name) {
      ctx.fillStyle = COLORS.axis;
      ctx.font = "10px system-ui, sans-serif";
      const tw = ctx.measureText(f.name).width;
      // Prefer just left of the feature; otherwise inside it, on the row's
      // centre line so the text never collides with the track header above.
      if (x0 - tw - 4 > plotLeft()) ctx.fillText(f.name, x0 - tw - 4, cy + 3);
      else ctx.fillText(f.name, Math.max(plotLeft(), x0) + 3, cy + 3);
    }
    hitFeats.push({ x0, x1, y0: y, y1: y + FEAT_ROW_H, f, track: trackName });
  }

  // Draw the read rows visible in the screen band [clipTop, clipBottom].
  // `rowsTop` is the screen y of read row 0 (already scroll-adjusted).
  function drawReads(data, rowsTop, clipTop, clipBottom) {
    const first = Math.max(0, Math.floor((clipTop - rowsTop) / READ_ROW_H) - 1);
    const last = Math.ceil((clipBottom - rowsTop) / READ_ROW_H) + 1;
    const pw = plotWidth();
    const exons = data.selected_exons || [];
    const selecting = exons.length > 0;
    for (const r of data.reads) {
      if (r.row < first || r.row > last) continue;
      const y = rowsTop + r.row * READ_ROW_H;
      const blocks = r.blocks || [[r.start, r.end]];
      const x0all = bpToPx(r.start);
      const x1all = bpToPx(r.end + 1);
      if (x1all < plotLeft() || x0all > plotLeft() + pw) continue;
      const base = r.strand === "-" ? COLORS.minus : COLORS.plus;
      // With an isoform selected, incompatible reads fade right out; otherwise
      // low mapping quality fades a little.
      ctx.globalAlpha = selecting
        ? (readCompatible(blocks, exons) ? 1 : FADED_ALPHA)
        : (r.mapq <= 0 ? 0.35 : r.mapq < 10 ? 0.6 : 1);
      // Connector line across the whole read (covers intron gaps).
      ctx.strokeStyle = COLORS.intron;
      ctx.beginPath();
      ctx.moveTo(x0all, y + READ_ROW_H / 2 - 0.5);
      ctx.lineTo(x1all, y + READ_ROW_H / 2 - 0.5);
      ctx.stroke();
      ctx.fillStyle = base;
      for (const [bs, be] of blocks) {
        const bx0 = bpToPx(bs);
        const bw = Math.max(1, bpToPx(be + 1) - bx0);
        ctx.fillRect(bx0, y + 0.5, bw, READ_ROW_H - 1.5);
      }
      // Strand arrowhead at the read's leading edge.
      if (x1all - x0all > 6) {
        ctx.beginPath();
        if (r.strand === "-") {
          ctx.moveTo(x0all, y + 0.5);
          ctx.lineTo(x0all - 3, y + READ_ROW_H / 2);
          ctx.lineTo(x0all, y + READ_ROW_H - 1);
        } else {
          ctx.moveTo(x1all, y + 0.5);
          ctx.lineTo(x1all + 3, y + READ_ROW_H / 2);
          ctx.lineTo(x1all, y + READ_ROW_H - 1);
        }
        ctx.fill();
      }
      ctx.globalAlpha = 1;
      hitReads.push({ x0: x0all, x1: x1all, y0: y, y1: y + READ_ROW_H, r });
    }
  }

  // Vertical scrollbar for the scrollable content band; also records geometry
  // for the pointer handlers. Sets the module-level `scrollbar` (or null).
  function drawScrollbar(top, height, contentH) {
    const maxScroll = Math.max(0, contentH - height);
    if (maxScroll <= 0) {
      scrollbar = null;
      return;
    }
    const x = canvas.clientWidth - 7;
    const thumbH = Math.max(24, (height / contentH) * height);
    const thumbY = top + (contentScrollY / maxScroll) * (height - thumbH);
    ctx.fillStyle = "rgba(0,0,0,0.05)";
    ctx.fillRect(x, top, 5, height);
    ctx.fillStyle = "rgba(0,0,0,0.28)";
    ctx.fillRect(x, thumbY, 5, thumbH);
    scrollbar = { x, top, avail: height, thumbY, thumbH, maxScroll };
  }

  // Everything below the ruler (annotation tracks + the alignments track) lives
  // in one vertically-scrollable band, so dense annotations can never squeeze
  // the reads out of existence — you just scroll down to them.
  function draw() {
    try {
      drawImpl();
    } catch (err) {
      // A bad frame must never wedge the widget's update stream.
      console.error("bam_viewer draw error:", err);
    }
  }

  function drawImpl() {
    setupCanvas();
    hitFeats = [];
    hitReads = [];
    scrollbar = null;
    ctx.clearRect(0, 0, canvas.clientWidth, VIEW_H);
    drawRuler();

    const readData = readState;
    const message = model.get("_message");
    if (message) {
      ctx.fillStyle = "#9a6a00";
      ctx.font = "12px system-ui, sans-serif";
      message.split("\n").forEach((line, i) => {
        ctx.fillText(line, plotLeft(), RULER_H + 20 + i * 15);
      });
      updateStatus(readData);
      return;
    }

    // Lay everything out in content space (y = 0 just below the ruler).
    const feats = featureState;
    const layout = [];
    let cy = 4;
    for (const track of feats) {
      const rows = Math.max(1, track.n_rows || 0);
      layout.push({ track, headerY: cy, rowsY: cy + TRACK_HEADER_H });
      cy += TRACK_HEADER_H + rows * FEAT_ROW_H + SECTION_GAP;
    }
    const readsHeaderY = cy;
    const readsRowsY = cy + TRACK_HEADER_H;
    const note = readData.note;
    const readsH = note ? 22 : (readData.n_rows || 0) * READ_ROW_H;
    const contentHeight = readsRowsY + readsH;

    const scrollTop = RULER_H;
    const scrollH = VIEW_H - scrollTop;
    contentScrollY = clamp(contentScrollY, 0, Math.max(0, contentHeight - scrollH));
    const off = scrollTop - contentScrollY; // content y -> screen y

    ctx.save();
    ctx.beginPath();
    ctx.rect(0, scrollTop, canvas.clientWidth, scrollH);
    ctx.clip();

    for (const { track, headerY, rowsY } of layout) {
      ctx.fillStyle = COLORS.axis;
      ctx.font = "bold 11px system-ui, sans-serif";
      const tnote = track.truncated ? " (truncated)" : "";
      ctx.fillText(`▸ ${track.name}${tnote}`, plotLeft(), off + headerY + 11);
      for (const f of track.features || []) drawFeature(f, off + rowsY, track.name);
    }

    ctx.fillStyle = COLORS.axis;
    ctx.font = "bold 11px system-ui, sans-serif";
    ctx.fillText("▸ Alignments", plotLeft(), off + readsHeaderY + 11);
    if (note) {
      ctx.fillStyle = "#9a6a00";
      ctx.font = "12px system-ui, sans-serif";
      ctx.fillText(note, plotLeft(), off + readsRowsY + 14);
    } else {
      drawReads(readData, off + readsRowsY, scrollTop, VIEW_H);
    }
    ctx.restore();

    drawScrollbar(scrollTop, scrollH, contentHeight);
    updateStatus(readData);
  }

  function updateStatus(readData) {
    if (document.activeElement !== regionInput) {
      regionInput.value = `${view.chrom}:${fmtBp(view.start)}-${fmtBp(view.end)}`;
    }
    let txt = model.get("_loading") ? "Loading… " : "";
    if (!model.get("_message")) {
      txt += `${readData.shown || 0} reads`;
      if (readData.truncated) txt += ` (sampled from ${readData.total})`;
      const sel = model.get("_selected") || {};
      const exons = readData.selected_exons || [];
      const reads = readData.reads || [];
      if (sel.transcript_id && exons.length && reads.length) {
        let compat = 0;
        for (const r of reads) {
          if (readCompatible(r.blocks || [[r.start, r.end]], exons)) compat++;
        }
        txt += ` · ${compat}/${reads.length} compatible with ${sel.transcript_id}`;
      } else if (sel.transcript_id && !exons.length) {
        txt += ` · ${sel.transcript_id} (no exon model)`;
      }
    }
    status.textContent = txt;
  }

  // ---- interaction ---------------------------------------------------------
  const DRAG_THRESHOLD = 4; // px of movement before a press counts as a drag

  function featureAt(x, y) {
    for (const h of hitFeats) {
      if (x >= h.x0 - 2 && x <= h.x1 + 2 && y >= h.y0 && y <= h.y1) return h;
    }
    return null;
  }
  function setSelection(sel) {
    model.set("_selected", sel);
    model.save_changes();
  }
  function toggleSelect(track, f) {
    const s = model.get("_selected") || {};
    if (s.transcript_id === f.transcript_id && s.track === track) setSelection({});
    else setSelection({ track, transcript_id: f.transcript_id });
  }

  let drag = null;
  let scrollDrag = null;
  overlay.addEventListener("mousedown", (e) => {
    const rect = overlay.getBoundingClientRect();
    const ox = e.clientX - rect.left;
    const oy = e.clientY - rect.top;
    // Grab the reads scrollbar if the press lands on it.
    if (scrollbar && ox >= scrollbar.x - 4 && oy >= scrollbar.top &&
        oy <= scrollbar.top + scrollbar.avail) {
      const onThumb =
        oy >= scrollbar.thumbY && oy <= scrollbar.thumbY + scrollbar.thumbH;
      if (!onThumb) {
        // Jump so the thumb centres on the click, then drag from there.
        const frac =
          (oy - scrollbar.top - scrollbar.thumbH / 2) /
          (scrollbar.avail - scrollbar.thumbH);
        contentScrollY = clamp(frac * scrollbar.maxScroll, 0, scrollbar.maxScroll);
        draw();
      }
      scrollDrag = { startOy: oy, startScroll: contentScrollY, sb: scrollbar };
      return;
    }
    drag = {
      x: e.clientX,
      downX: e.clientX,
      downY: e.clientY,
      ox,
      oy,
      start: view.start,
      end: view.end,
      moved: false,
    };
  });
  window.addEventListener("mousemove", (e) => {
    if (scrollDrag) {
      const rect = overlay.getBoundingClientRect();
      const dy = e.clientY - rect.top - scrollDrag.startOy;
      const range = scrollDrag.sb.avail - scrollDrag.sb.thumbH;
      contentScrollY = clamp(
        scrollDrag.startScroll + (dy / range) * scrollDrag.sb.maxScroll,
        0,
        scrollDrag.sb.maxScroll
      );
      draw();
      return;
    }
    if (!drag) return;
    // Stay a "click" until the pointer clearly moves, so selecting a feature
    // doesn't accidentally pan the view.
    if (
      !drag.moved &&
      Math.abs(e.clientX - drag.downX) < DRAG_THRESHOLD &&
      Math.abs(e.clientY - drag.downY) < DRAG_THRESHOLD
    ) {
      return;
    }
    if (!drag.moved) {
      drag.moved = true;
      overlay.classList.add("bamv-dragging");
    }
    const span = drag.end - drag.start + 1;
    const dbp = ((e.clientX - drag.x) / plotWidth()) * span;
    view.start = drag.start - dbp;
    view.end = view.start + span - 1;
    clampView();
    draw();
    scheduleCommit();
  });
  window.addEventListener("mouseup", () => {
    if (scrollDrag) {
      scrollDrag = null;
      return;
    }
    if (!drag) return;
    const d = drag;
    drag = null;
    overlay.classList.remove("bamv-dragging");
    if (d.moved) {
      commit();
      return;
    }
    // A click (no meaningful drag): select the isoform under the cursor, or
    // clear the selection when clicking empty space.
    const hit = featureAt(d.ox, d.oy);
    if (hit && hit.f.transcript_id) toggleSelect(hit.track, hit.f);
    else if (!hit && (model.get("_selected") || {}).transcript_id) setSelection({});
  });

  overlay.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      // Shift+wheel (or wheel while the reads overflow and Shift is held)
      // scrolls the alignments vertically instead of zooming.
      if (e.shiftKey && scrollbar) {
        contentScrollY = clamp(
          contentScrollY + e.deltaY,
          0,
          scrollbar.maxScroll
        );
        draw();
        return;
      }
      const rect = overlay.getBoundingClientRect();
      const anchorBp = pxToBp(e.clientX - rect.left);
      const factor = Math.exp(e.deltaY * 0.0015);
      const newW = clamp(width() * factor, MIN_WIDTH, maxWidth());
      const frac = (anchorBp - view.start) / width();
      view.start = anchorBp - frac * newW;
      view.end = view.start + newW - 1;
      clampView();
      draw();
      scheduleCommit();
    },
    { passive: false }
  );

  function panBy(frac) {
    const d = width() * frac;
    view.start += d;
    view.end += d;
    clampView();
    draw();
    commit();
  }
  function zoomBy(factor) {
    const center = (view.start + view.end) / 2;
    const newW = clamp(width() * factor, MIN_WIDTH, maxWidth());
    view.start = center - newW / 2;
    view.end = view.start + newW - 1;
    clampView();
    draw();
    commit();
  }
  el.querySelector(".bamv-left").onclick = () => panBy(-0.4);
  el.querySelector(".bamv-right").onclick = () => panBy(0.4);
  el.querySelector(".bamv-zoomin").onclick = () => zoomBy(0.5);
  el.querySelector(".bamv-zoomout").onclick = () => zoomBy(2);

  function applyRegionInput() {
    const m = regionInput.value.match(
      /^\s*([^:\s]+)\s*(?::\s*([\d,]+)\s*-\s*([\d,]+))?\s*$/
    );
    if (!m) {
      regionInput.classList.add("bamv-invalid");
      return;
    }
    regionInput.classList.remove("bamv-invalid");
    view.chrom = m[1];
    if (m[2]) {
      let s = parseInt(m[2].replace(/,/g, ""), 10);
      let en = parseInt(m[3].replace(/,/g, ""), 10);
      if (en < s) [s, en] = [en, s];
      view.start = s;
      view.end = en;
    }
    regionInput.blur();
    clampView();
    draw();
    commit();
  }
  el.querySelector(".bamv-go").onclick = applyRegionInput;
  regionInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") applyRegionInput();
  });

  // Tooltip hit-testing.
  overlay.addEventListener("mousemove", (e) => {
    if (drag) return;
    const rect = overlay.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const yy = e.clientY - rect.top;
    let hit = null;
    for (const h of hitReads) {
      if (x >= h.x0 - 2 && x <= h.x1 + 2 && yy >= h.y0 && yy <= h.y1) {
        const r = h.r;
        hit = `read ${view.chrom}:${fmtBp(r.start)}-${fmtBp(r.end)}<br>strand ${r.strand} · MAPQ ${r.mapq}`;
        break;
      }
    }
    if (!hit) {
      for (const h of hitFeats) {
        if (x >= h.x0 - 2 && x <= h.x1 + 2 && yy >= h.y0 && yy <= h.y1) {
          const f = h.f;
          hit = `<b>${f.name || "feature"}</b><br>${view.chrom}:${fmtBp(f.start)}-${fmtBp(f.end)} (${f.strand})`;
          if (f.transcript_id) {
            hit += `<br>transcript ${f.transcript_id}`;
            hit += isSelected(h.track, f)
              ? "<br><i>click to deselect</i>"
              : "<br><i>click to select isoform</i>";
          }
          break;
        }
      }
    }
    if (hit) {
      tooltip.innerHTML = hit;
      tooltip.style.display = "block";
      tooltip.style.left = Math.min(x + 12, overlay.clientWidth - 180) + "px";
      tooltip.style.top = yy + 12 + "px";
      overlay.style.cursor = "pointer";
    } else {
      tooltip.style.display = "none";
      overlay.style.cursor = "grab";
    }
  });
  overlay.addEventListener("mouseleave", () => {
    tooltip.style.display = "none";
  });

  // ---- model -> view sync --------------------------------------------------
  // Adopt the model's region only on an *external* change (e.g. a Python-side
  // goto()); our own commits set selfUpdating so this is skipped for them.
  function onRegionTrait() {
    if (selfUpdating || drag) return;
    const v = model.get("_view") || [];
    if (v[0] !== view.chrom || v[1] !== view.start || v[2] !== view.end) {
      // Adopt an external view change (e.g. Python goto) and redraw. When the
      // view already matches (our own commit echoing back), do nothing — the
      // reload's change:_read_data handles the redraw, avoiding a transient
      // frame drawn against half-applied data.
      view.chrom = v[0];
      view.start = v[1];
      view.end = v[2];
      draw();
    }
  }
  model.on("change:_read_data", () => {
    const rd = model.get("_read_data") || {};
    if (rd._v != null && rd._v <= lastReadV) return; // stale re-sync, ignore
    if (rd._v != null) lastReadV = rd._v;
    readState = rd;
    draw();
  });
  model.on("change:_feature_data", () => {
    const fd = model.get("_feature_data") || {};
    if (fd._v != null && fd._v <= lastFeatV) return; // stale re-sync, ignore
    if (fd._v != null) lastFeatV = fd._v;
    featureState = fd.tracks || [];
    draw();
  });
  model.on("change:_selected", draw);
  model.on("change:_message", draw);
  model.on("change:_loading", () => updateStatus(readState));
  model.on("change:_view", onRegionTrait);

  const ro = new ResizeObserver(() => draw());
  ro.observe(el.querySelector(".bamv-canvas-wrap"));

  draw();
  return () => ro.disconnect();
}

export default { render };
