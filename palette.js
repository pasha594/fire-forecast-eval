/* Data-encoding colors shared by the pages (time-of-arrival bands, hotspot age
   ramp, burn scar, perimeter, containment), editable per browser through a
   settings overlay styled by brand.css. Classic script: defines window.Palette. */
(function () {
  "use strict";
  // v1 settings were saved against the pre-reversal defaults; a stored "reversed"
  // flag would flip the new defaults straight back, so they are not carried over
  const KEY = "cornea-palette-v2";
  const TOA_EDGES = Object.freeze([12, 24, 48, 72, 96, 120, 168, 240, Infinity]);
  const HOT_STOPS_H = Object.freeze([24, 48, 72, 96, 168]);
  const HOT_MAX_H = HOT_STOPS_H[HOT_STOPS_H.length - 1];
  const DEFAULTS = Object.freeze({
    // <=12h ... >240h: Figma "Colour allocation" (plus one extrapolated purple), run
    // from purple for imminent arrival to teal for late
    toa: Object.freeze(["#B05DBD", "#A569C4", "#9A75CB", "#8F82D2", "#8691D9",
                        "#7CA2DE", "#90B9E7", "#A2CEE9", "#A9DCD6"]),
    toaReversed: false,
    // freshest (<=24h, deep crimson) ... 7 days old (salmon)
    hotspot: Object.freeze(["#A81F2D", "#CC3630", "#E85538", "#FA7A5E", "#FF9E8A"]),
    hotspotReversed: false,
    scar: "#B9A895",
    perimeter: "#A4161A",
    contained: "#D81B7A",
  });
  const LE = "≤";

  const isHex = v => typeof v === "string" && /^#[0-9a-f]{6}$/i.test(v);
  const isColors = (v, n) => Array.isArray(v) && v.length === n && v.every(isHex);
  const isBool = v => typeof v === "boolean";
  const VALID = {
    toa: v => isColors(v, TOA_EDGES.length), toaReversed: isBool,
    hotspot: v => isColors(v, HOT_STOPS_H.length), hotspotReversed: isBool,
    scar: isHex, perimeter: isHex, contained: isHex,
  };

  function copy(s) {
    return { ...s, toa: s.toa.slice(), hotspot: s.hotspot.slice() };
  }

  function merge(base, patch) {
    const out = copy(base);
    if (patch && typeof patch === "object")
      for (const k of Object.keys(VALID))
        if (VALID[k](patch[k])) out[k] = Array.isArray(patch[k]) ? patch[k].slice() : patch[k];
    return out;
  }

  function same(a, b) {
    if (Array.isArray(a)) return a.every((x, i) => same(x, b[i]));
    return typeof a === "string" ? a.toLowerCase() === b.toLowerCase() : a === b;
  }

  function storage() {
    try { return window.localStorage || null; } catch (e) { return null; }
  }

  function load() {
    try {
      const s = storage(), raw = s && s.getItem(KEY);
      return merge(DEFAULTS, raw ? JSON.parse(raw) : null);
    } catch (e) { return copy(DEFAULTS); }
  }

  function persist() {
    // only keys that differ from DEFAULTS are stored, so later changes to the
    // defaults still reach browsers that customized something else
    const diff = {};
    for (const k of Object.keys(VALID)) if (!same(current[k], DEFAULTS[k])) diff[k] = current[k];
    try {
      const s = storage();
      if (!s) return;
      if (Object.keys(diff).length) s.setItem(KEY, JSON.stringify(diff));
      else s.removeItem(KEY);
    } catch (e) { /* quota or blocked storage: the palette still applies to this page */ }
  }

  let current = load();
  let panel = null, ui = null, lastFocus = null, pending = null, timer = 0;
  const listeners = [];

  function notify() {
    refreshPanel();
    for (const fn of listeners.slice()) {
      try { fn(get()); } catch (e) { console.error(e); }
    }
  }

  function get() { return copy(current); }

  function set(partial) {
    current = merge(current, partial);
    persist();
    notify();
  }

  function reset() {
    clearTimeout(timer);
    pending = null;
    current = copy(DEFAULTS);
    try { const s = storage(); if (s) s.removeItem(KEY); } catch (e) { /* blocked storage */ }
    notify();
  }

  function onChange(fn) {
    listeners.push(fn);
    return () => { const i = listeners.indexOf(fn); if (i >= 0) listeners.splice(i, 1); };
  }

  window.addEventListener("storage", ev => {
    if (ev.key !== KEY && ev.key !== null) return;
    current = load();
    notify();
  });

  // ------------------------------------------------------------ color math

  const toaAt = i => current.toa[current.toaReversed ? current.toa.length - 1 - i : i];
  const hotAt = i => current.hotspot[current.hotspotReversed ? current.hotspot.length - 1 - i : i];

  function toaColors() {
    return TOA_EDGES.map((e, i) => toaAt(i));
  }

  // called once per raster pixel, so no allocation here
  function toaColor(hours) {
    for (let i = 0; i < TOA_EDGES.length; i++) if (hours <= TOA_EDGES[i]) return toaAt(i);
    return null;
  }

  function rgb(hex) {
    const n = parseInt(hex.slice(1), 16);
    return [n >> 16, (n >> 8) & 255, n & 255];
  }

  function mix(a, b, t) {
    const p = rgb(a), q = rgb(b);
    return "#" + p.map((v, k) => Math.round(v + (q[k] - v) * t).toString(16).padStart(2, "0")).join("");
  }

  function hotspotColor(ageH) {
    if (typeof ageH !== "number" || !(ageH >= 0 && ageH <= HOT_MAX_H)) return null;
    if (ageH <= HOT_STOPS_H[0]) return hotAt(0);
    let j = 1;
    while (ageH > HOT_STOPS_H[j]) j++;
    if (ageH === HOT_STOPS_H[j]) return hotAt(j);
    const t = (ageH - HOT_STOPS_H[j - 1]) / (HOT_STOPS_H[j] - HOT_STOPS_H[j - 1]);
    return mix(hotAt(j - 1), hotAt(j), t);
  }

  function hotspotGradient() {
    const pos = h => +((1 - h / HOT_MAX_H) * 100).toFixed(3) + "%";
    const stops = HOT_STOPS_H.map((h, j) => `${hotAt(j)} ${pos(h)}`).reverse();
    stops.push(`${hotAt(0)} 100%`);
    return `linear-gradient(to right, ${stops.join(", ")})`;
  }

  // ------------------------------------------------------- settings overlay

  const GEAR_SVG = (() => {
    let d = "";
    for (let k = 0; k < 6; k++)
      for (const [r, deg] of [[7.4, -17], [10.4, -10], [10.4, 10], [7.4, 17]]) {
        const a = (k * 60 + deg) * Math.PI / 180;
        d += `${d ? "L" : "M"}${(12 + r * Math.sin(a)).toFixed(2)} ${(12 - r * Math.cos(a)).toFixed(2)}`;
      }
    return `<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor"
      stroke-width="1.8" stroke-linejoin="round" aria-hidden="true" focusable="false">
      <path d="${d}Z"/><circle cx="12" cy="12" r="3"/></svg>`;
  })();

  function el(tag, attrs, ...kids) {
    const n = document.createElement(tag);
    for (const k in attrs || {}) n.setAttribute(k, attrs[k]);
    n.append(...kids);
    return n;
  }

  // the color pickers fire "input" on every drag tick and each set() re-renders
  // the map's raster, so picks are batched and applied 150ms after the last tick
  function stage(key, band, value) {
    pending = pending || {};
    if (band == null) pending[key] = value;
    else {
      const arr = pending[key] = pending[key] || current[key].slice();
      arr[current[key + "Reversed"] ? arr.length - 1 - band : band] = value;
    }
    clearTimeout(timer);
    timer = setTimeout(flush, 150);
  }

  function flush() {
    clearTimeout(timer);
    const p = pending;
    pending = null;
    if (p) set(p);
  }

  function colorRow(label, context, key, band, wide) {
    const input = el("input", { type: "color", "aria-label": `${context}: ${label}` });
    input.addEventListener("input", () => stage(key, band, input.value));
    input.addEventListener("change", flush);
    const row = el("label", { class: wide ? "palette-row palette-wide" : "palette-row" },
      input, el("span", {}, label));
    return { row, input };
  }

  function reverseBox(key) {
    const input = el("input", { type: "checkbox" });
    input.addEventListener("change", () => { flush(); set({ [key]: input.checked }); });
    return { row: el("label", { class: "palette-check" }, input, "Reverse order"), input };
  }

  const grid = rows => el("div", { class: "palette-grid" }, ...rows.map(r => r.row));
  const section = (title, ...kids) => el("section", { class: "palette-sec" }, el("h3", {}, title), ...kids);

  function build() {
    const toaLabels = TOA_EDGES.map((e, i) => e === Infinity ? `>${TOA_EDGES[i - 1]}h` : `${LE}${e}h`);
    const toa = toaLabels.map((lab, i) => colorRow(lab, "Time of arrival", "toa", i));
    const hot = [`${LE}24h`, "48h", "72h", "4d", "7d"]
      .map((lab, i) => colorRow(lab, "Hotspot age", "hotspot", i));
    const scar = colorRow("Burn scar (older than 7 days)", "Hotspots", "scar", null, true);
    const perimeter = colorRow("Perimeter", "Boundaries", "perimeter", null);
    const contained = colorRow("Contained", "Boundaries", "contained", null);
    const toaRev = reverseBox("toaReversed"), hotRev = reverseBox("hotspotReversed");

    const close = el("button", { type: "button", class: "palette-close",
                                 "aria-label": "Close color settings", title: "Close" }, "×");
    close.addEventListener("click", () => closeSettings(true));
    const resetBtn = el("button", { type: "button", class: "palette-btn" }, "Reset to defaults");
    resetBtn.addEventListener("click", reset);
    const copyBtn = el("button", { type: "button", class: "palette-btn palette-btn-primary" }, "Copy settings");
    copyBtn.addEventListener("click", copySettings);
    const copyBox = el("textarea", { class: "palette-copybox", rows: "9", readonly: "",
                                     spellcheck: "false", "aria-label": "Color settings JSON" });
    copyBox.hidden = true;

    ui = {
      toa: toa.map(r => r.input), hot: hot.map(r => r.input),
      scar: scar.input, perimeter: perimeter.input, contained: contained.input,
      toaRev: toaRev.input, hotRev: hotRev.input,
      toaStrip: el("div", { class: "palette-strip", "aria-hidden": "true" }, ...toa.map(() => el("span"))),
      hotBar: el("div", { class: "palette-bar", "aria-hidden": "true" }),
      close, copyBtn, copyBox, copiedTimer: 0,
    };

    panel = el("div", { class: "palette-panel", role: "dialog", "aria-labelledby": "palette-title" },
      el("div", { class: "palette-head" }, el("h2", { id: "palette-title" }, "Color settings"), close),
      section("Time of arrival", ui.toaStrip, grid(toa), toaRev.row),
      section("Hotspots (detection age)", ui.hotBar,
        el("div", { class: "palette-ends", "aria-hidden": "true" }, el("span", {}, "7d"), el("span", {}, "0h")),
        grid(hot), hotRev.row, grid([scar])),
      section("Boundaries", grid([perimeter, contained])),
      el("div", { class: "palette-foot" },
        el("div", { class: "palette-actions" }, resetBtn, copyBtn), copyBox,
        el("p", { class: "palette-note" },
          "Saved in this browser only — copy settings to share a palette.")));
    panel.hidden = true;
    document.body.append(panel);

    document.addEventListener("keydown", ev => {
      if (ev.key === "Escape" && !panel.hidden) closeSettings(true);
    });
    document.addEventListener("pointerdown", ev => {
      if (panel.hidden || panel.contains(ev.target)) return;
      if (ev.target instanceof Element && ev.target.closest(".palette-gear")) return;
      closeSettings(false);
    }, true);
  }

  function refreshPanel() {
    if (!ui) return;
    const toa = toaColors();
    toa.forEach((c, i) => { ui.toa[i].value = c; ui.toaStrip.children[i].style.background = c; });
    ui.hot.forEach((input, i) => { input.value = hotAt(i); });
    ui.hotBar.style.background = hotspotGradient();
    ui.toaRev.checked = current.toaReversed;
    ui.hotRev.checked = current.hotspotReversed;
    ui.scar.value = current.scar;
    ui.perimeter.value = current.perimeter;
    ui.contained.value = current.contained;
    ui.copyBox.hidden = true;
  }

  function copySettings() {
    flush();
    const json = JSON.stringify({ v: 1, ...get() }, null, 2);
    const fallback = () => {
      ui.copyBox.value = json;
      ui.copyBox.hidden = false;
      ui.copyBox.focus();
      ui.copyBox.select();
    };
    if (!navigator.clipboard || !navigator.clipboard.writeText) { fallback(); return; }
    navigator.clipboard.writeText(json).then(() => {
      ui.copyBtn.style.minWidth = `${ui.copyBtn.offsetWidth}px`;
      ui.copyBtn.textContent = "Copied";
      clearTimeout(ui.copiedTimer);
      ui.copiedTimer = setTimeout(() => { ui.copyBtn.textContent = "Copy settings"; }, 1500);
    }, fallback);
  }

  function openSettings() {
    if (!panel) build();
    if (!panel.hidden) return;
    lastFocus = document.activeElement;
    refreshPanel();
    panel.hidden = false;
    ui.close.focus();
  }

  function closeSettings(restoreFocus) {
    if (!panel || panel.hidden) return;
    flush();
    panel.hidden = true;
    ui.copyBox.hidden = true;
    if (restoreFocus && lastFocus && lastFocus.isConnected && lastFocus.focus) lastFocus.focus();
    lastFocus = null;
  }

  function gearButton() {
    const b = el("button", { type: "button", class: "palette-gear",
                             "aria-label": "Color settings", title: "Color settings" });
    b.innerHTML = GEAR_SVG;
    b.addEventListener("click", () => openSettings());
    return b;
  }

  window.Palette = {
    TOA_EDGES, HOT_STOPS_H, DEFAULTS,
    get, set, reset, onChange,
    toaColors, toaColor, hotspotColor, hotspotGradient,
    openSettings, gearButton,
  };
})();
