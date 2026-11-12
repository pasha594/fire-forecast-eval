/* One satellite detection reported more than once counts once, on the map and in the report
   as in the scores (analyze.py drop_repeats, HS_DUP_S). Cornea's /hotspots feed often carries
   a detection twice: one satellite pass, 0-2 minutes apart, the same reading, at locations up
   to ~300 m apart where the products it came through disagree on geolocation (~20 m on
   Skillet, 40-280 m over Dome's steep ground). Two detections from one satellite at most
   3 minutes apart are one when
    - they carry the same reading (brightness within 0.1 K, FRP within 5%) within 300 m:
      different pixels both saturated at 367 K differ in FRP, so they stay, or
    - they're closer than a sixth of the instrument's pixel (VIIRS 62 m, MODIS 167 m,
      Landsat 5 m), which two different pixels of one pass can't be.
   HotRepeats.make() -> { add(d) }: true for a new detection (now kept), false for a repeat
   of one already kept. d = { src, t (epoch ms), lon, lat, bright, frp }. */
(function () {
  const DUP_MS = 3 * 60000, COPY_M = 300;
  const pixelM = src => /landsat/i.test(src || "") ? 30 : /modis|terra|aqua/i.test(src || "") ? 1000 : 375;
  const num = v => (v === null || v === undefined || v === "" || !isFinite(+v)) ? null : +v;

  function same(a, b) {
    if (a.src !== b.src || Math.abs(a.t - b.t) > DUP_MS) return false;
    const kx = 111320 * Math.cos(a.lat * Math.PI / 180);
    const d = Math.hypot((a.lon - b.lon) * kx, (a.lat - b.lat) * 110540);
    if (d < pixelM(a.src) / 6) return true;
    if (d > COPY_M) return false;
    const ba = num(a.bright), bb = num(b.bright), fa = num(a.frp), fb = num(b.frp);
    return ba !== null && bb !== null && fa !== null && fb !== null &&
      Math.abs(ba - bb) <= 0.1 && Math.abs(fa - fb) <= Math.max(0.05, 0.05 * Math.max(fa, fb));
  }

  function make() {
    const byMin = new Map();   // "source|minute" -> the detections kept in it
    return {
      add(d) {
        // detections may arrive out of time order (the map loads days as needed): look both ways
        const m = Math.round(d.t / 60000);
        for (let k = m - 3; k <= m + 3; k++) {
          const near = byMin.get(`${d.src}|${k}`);
          if (near && near.some(o => same(o, d))) return false;
        }
        const key = `${d.src}|${m}`;
        if (!byMin.has(key)) byMin.set(key, []);
        byMin.get(key).push(d);
        return true;
      },
    };
  }

  window.HotRepeats = { make, same };
})();
