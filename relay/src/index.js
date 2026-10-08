// Pyrecast relay for the forecast map (map.html). Pyrecast's geoserver refuses
// cross-origin reads from every site but pyrecast.org, and the map has to read pixels:
// WebGL draws its tiles, and time of arrival is decoded from the raw coverage. This
// Worker passes the geoserver's responses through with CORS headers for our pages.
//
// It forwards only the two requests the map makes, rebuilt from an allowlist of
// parameters (so no other operation, format or layer can ride along): WMS GetMap tiles
// of Pyrecast's forecast layers, and WCS GetCoverage of their time-of-arrival rasters.
// It does no processing (the free plan allows ~10ms CPU a request). Responses carry
// Cache-Control for browsers; on a *.workers.dev address Cloudflare's own cache does
// nothing, so a repeat from another viewer reaches Pyrecast again.
const UPSTREAM = "https://geoserver-usw1.pyrecast.org/geoserver02/";
const LAYER = /^fire-spread-forecast_[a-z0-9-]+_\d{8}_\d{6}(:|__)elmfire_landfire_\d{1,2}_[a-z-]+$/;
// path -> the parameters it accepts (lower-cased names) and the values they must have
const ENDPOINTS = {
  wms: {
    params: ["service", "request", "version", "layers", "styles", "format", "transparent", "srs", "width",
             "height", "bbox", "time"],
    fixed: { service: "WMS", request: "GetMap", format: "image/png", width: "256", height: "256" },
    layer: "layers",
  },
  ows: {
    params: ["service", "version", "request", "coverageid", "format", "geotiff:compression"],
    fixed: { service: "WCS", request: "GetCoverage", format: "image/geotiff" },
    layer: "coverageid",
  },
};
const ORIGINS = [/^https:\/\/pasha594\.github\.io$/, /^http:\/\/(localhost|127\.0\.0\.1)(:\d+)?$/];
const HOME = "https://pasha594.github.io";
// tiles of a live run can still fill in as Pyrecast publishes its hours; a coverage is final
const MAX_AGE = { wms: 900, ows: 3600 };

function corsHeaders(req) {
  const origin = req.headers.get("Origin") || "";
  return {
    "Access-Control-Allow-Origin": ORIGINS.some(re => re.test(origin)) ? origin : HOME,
    "Access-Control-Allow-Methods": "GET",
    "Access-Control-Max-Age": "86400",
    "Vary": "Origin",
  };
}

function upstreamUrl(url) {
  // the request rebuilt from its allowlisted parameters, or null if it isn't one the map makes
  const service = url.pathname.replace(/^\/+|\/+$/g, ""), ep = ENDPOINTS[service];
  if (!ep) return null;
  const seen = new Map();
  for (const [k, v] of url.searchParams) {
    const key = k.toLowerCase();
    if (!ep.params.includes(key) || seen.has(key)) return null;   // unknown or repeated (any case)
    seen.set(key, v);
  }
  for (const [k, v] of Object.entries(ep.fixed)) {
    if ((seen.get(k) || "").toLowerCase() !== v.toLowerCase()) return null;
  }
  if (!LAYER.test(seen.get(ep.layer) || "")) return null;
  return UPSTREAM + service + "?" + new URLSearchParams([...seen]).toString();
}

export default {
  async fetch(req) {
    const cors = corsHeaders(req);
    if (req.method === "OPTIONS") return new Response(null, { status: 204, headers: cors });
    if (req.method !== "GET") return new Response("GET only", { status: 405, headers: cors });

    const url = new URL(req.url), target = upstreamUrl(url);
    if (!target) return new Response("not relayed", { status: 400, headers: cors });

    // a server-side fetch carries no Origin header, which the geoserver accepts
    const up = await fetch(target, { headers: { "User-Agent": "cornea-projector-relay (pasha594.github.io)" } });
    const res = new Response(up.body, up);
    // the geoserver reports some errors as 200 + XML: only images may be cached
    const image = up.ok && /^image\//.test(up.headers.get("Content-Type") || "");
    const service = url.pathname.replace(/^\/+|\/+$/g, "");
    res.headers.set("Cache-Control", image ? `public, max-age=${MAX_AGE[service]}` : "no-store");
    res.headers.delete("Set-Cookie");   // the geoserver's flow-control cookie is no use to a browser
    for (const [k, v] of Object.entries(cors)) res.headers.set(k, v);
    return res;
  },
};
