import { lookup } from "node:dns/promises";
import http from "node:http";
import https from "node:https";
import net from "node:net";

const DEFAULT_TIMEOUT_MS = 15_000;
const DEFAULT_MAX_BYTES = 1_000_000;
const MAX_REDIRECTS = 5;
const MAX_URL_CHARS = 4_096;
const MAX_TITLE_CHARS = 512;
const MAX_SNIPPET_CHARS = 2_000;
const SEARCH_STOP_WORDS = new Set([
  "about",
  "and",
  "are",
  "can",
  "current",
  "find",
  "for",
  "from",
  "how",
  "information",
  "official",
  "the",
  "this",
  "what",
  "with",
]);

function decodeHtmlEntities(value) {
  return value
    .replace(/&#(\d+);/g, (_match, decimal) =>
      String.fromCodePoint(Number.parseInt(decimal, 10)),
    )
    .replace(/&#x([0-9a-f]+);/gi, (_match, hexadecimal) =>
      String.fromCodePoint(Number.parseInt(hexadecimal, 16)),
    )
    .replaceAll("&amp;", "&")
    .replaceAll("&quot;", '"')
    .replaceAll("&#39;", "'")
    .replaceAll("&apos;", "'")
    .replaceAll("&lt;", "<")
    .replaceAll("&gt;", ">")
    .replaceAll("&nbsp;", " ");
}

function stripHtml(value) {
  return decodeHtmlEntities(value.replace(/<[^>]*>/g, " "))
    .replace(/\s+/g, " ")
    .trim();
}

function blockedIpv4(address) {
  const octets = address.split(".").map((part) => Number.parseInt(part, 10));
  if (
    octets.length !== 4 ||
    octets.some((part) => !Number.isInteger(part) || part < 0 || part > 255)
  ) {
    return true;
  }
  const [a, b] = octets;
  return (
    a === 0 ||
    a === 10 ||
    a === 127 ||
    (a === 100 && b >= 64 && b <= 127) ||
    (a === 169 && b === 254) ||
    (a === 172 && b >= 16 && b <= 31) ||
    (a === 192 && b === 0) ||
    (a === 192 && b === 88 && octets[2] === 99) ||
    (a === 192 && b === 168) ||
    (a === 198 && (b === 18 || b === 19)) ||
    (a === 198 && b === 51 && octets[2] === 100) ||
    (a === 203 && b === 0 && octets[2] === 113) ||
    a >= 224
  );
}

function ipv6Groups(address) {
  const pieces = address.split("::");
  if (pieces.length > 2) return null;
  const parseSide = (side) => {
    if (!side) return [];
    const groups = side.split(":");
    const last = groups.at(-1);
    if (last?.includes(".")) {
      const octets = last.split(".").map((part) => Number.parseInt(part, 10));
      if (
        octets.length !== 4 ||
        octets.some((part) => !Number.isInteger(part) || part < 0 || part > 255)
      ) {
        return null;
      }
      groups.splice(
        -1,
        1,
        ((octets[0] << 8) | octets[1]).toString(16),
        ((octets[2] << 8) | octets[3]).toString(16),
      );
    }
    const parsed = groups.map((part) => Number.parseInt(part, 16));
    return parsed.some(
      (part, index) =>
        !/^[0-9a-f]{1,4}$/i.test(groups[index]) ||
        !Number.isInteger(part) ||
        part < 0 ||
        part > 0xffff,
    )
      ? null
      : parsed;
  };
  const left = parseSide(pieces[0]);
  const right = parseSide(pieces[1] ?? "");
  if (!left || !right) return null;
  if (pieces.length === 1) return left.length === 8 ? left : null;
  const missing = 8 - left.length - right.length;
  return missing >= 1 ? [...left, ...Array(missing).fill(0), ...right] : null;
}

export function isBlockedAddress(address) {
  const normalized = address.toLowerCase().split("%")[0];
  const family = net.isIP(normalized);
  if (family === 4) {
    return blockedIpv4(normalized);
  }
  if (family !== 6) {
    return true;
  }
  const groups = ipv6Groups(normalized);
  if (!groups) return true;
  const allZero = groups.every((group) => group === 0);
  const loopback = groups.slice(0, 7).every((group) => group === 0) && groups[7] === 1;
  const ipv4Compatible = groups.slice(0, 6).every((group) => group === 0);
  const ipv4Mapped =
    groups.slice(0, 5).every((group) => group === 0) && groups[5] === 0xffff;
  const nat64WellKnown =
    groups[0] === 0x64 &&
    groups[1] === 0xff9b &&
    groups.slice(2, 6).every((group) => group === 0);
  const nat64Local = groups[0] === 0x64 && groups[1] === 0xff9b && groups[2] === 1;
  return (
    allZero ||
    loopback ||
    ipv4Compatible ||
    ipv4Mapped ||
    nat64WellKnown ||
    nat64Local ||
    (groups[0] >= 0xfc00 && groups[0] <= 0xfdff) ||
    (groups[0] & 0xffc0) === 0xfe80 ||
    (groups[0] & 0xffc0) === 0xfec0 ||
    (groups[0] & 0xff00) === 0xff00 ||
    (groups[0] === 0x100 && groups.slice(1, 4).every((group) => group === 0)) ||
    (groups[0] === 0x2001 && groups[1] <= 0x01ff) ||
    (groups[0] === 0x2001 && groups[1] === 0x0db8) ||
    groups[0] === 0x2002 ||
    (groups[0] & 0xfff0) === 0x3ff0
  );
}

export function normalizePublicUrl(rawUrl) {
  const text = String(rawUrl).trim();
  if (!text || text.length > MAX_URL_CHARS) {
    throw new Error(`URL must contain 1 to ${MAX_URL_CHARS} characters`);
  }
  let url;
  try {
    url = new URL(text);
  } catch {
    throw new Error("URL must be an absolute HTTP or HTTPS URL");
  }
  if (!['http:', 'https:'].includes(url.protocol)) {
    throw new Error("Only HTTP and HTTPS URLs are allowed");
  }
  if (url.username || url.password) {
    throw new Error("URLs containing credentials are not allowed");
  }
  const hostname = url.hostname.toLowerCase().replace(/^\[|\]$/g, "");
  if (
    !hostname ||
    hostname === "localhost" ||
    hostname.endsWith(".localhost") ||
    hostname.endsWith(".local") ||
    hostname.endsWith(".internal") ||
    hostname.endsWith(".lan")
  ) {
    throw new Error("Local and private-network hostnames are not allowed");
  }
  if (net.isIP(hostname) && isBlockedAddress(hostname)) {
    throw new Error("Private, loopback, link-local, and reserved addresses are not allowed");
  }
  url.hash = "";
  return url;
}

async function resolvePublicAddress(url) {
  const hostname = url.hostname.replace(/^\[|\]$/g, "");
  if (net.isIP(hostname)) {
    if (isBlockedAddress(hostname)) {
      throw new Error("Target address is not publicly routable");
    }
    return { address: hostname, family: net.isIP(hostname) };
  }
  const addresses = await lookup(hostname, { all: true, verbatim: true });
  if (addresses.length === 0) {
    throw new Error(`Could not resolve public host ${hostname}`);
  }
  if (addresses.some((entry) => isBlockedAddress(entry.address))) {
    throw new Error(`Host ${hostname} resolves to a blocked network address`);
  }
  return addresses[0];
}

function pinnedLookup(address, family) {
  return (_hostname, options, callback) => {
    if (typeof options === "function") {
      callback = options;
      options = {};
    }
    if (options?.all) {
      callback(null, [{ address, family }]);
    } else {
      callback(null, address, family);
    }
  };
}

async function requestOnce(url, { signal, timeoutMs, maxBytes }) {
  const { address, family } = await resolvePublicAddress(url);
  const client = url.protocol === "https:" ? https : http;
  return await new Promise((resolve, reject) => {
    let settled = false;
    const finish = (callback, value) => {
      if (settled) return;
      settled = true;
      signal?.removeEventListener("abort", abort);
      callback(value);
    };
    const request = client.request(
      url,
      {
        method: "GET",
        headers: {
          Accept: "text/html,application/xhtml+xml,application/json,text/plain,application/xml;q=0.8,*/*;q=0.1",
          "Accept-Encoding": "identity",
          "User-Agent": "GENtle-Agent-Web-Research/1.0 (+https://github.com/smoe/gentle_rs)",
        },
        lookup: pinnedLookup(address, family),
      },
      (response) => {
        const declaredLength = Number.parseInt(response.headers["content-length"] ?? "0", 10);
        if (declaredLength > maxBytes) {
          response.resume();
          finish(reject, new Error(`Response exceeds the ${maxBytes}-byte limit`));
          return;
        }
        const chunks = [];
        let byteCount = 0;
        response.on("data", (chunk) => {
          byteCount += chunk.length;
          if (byteCount > maxBytes) {
            request.destroy(new Error(`Response exceeds the ${maxBytes}-byte limit`));
            return;
          }
          chunks.push(chunk);
        });
        response.on("end", () =>
          finish(resolve, {
            statusCode: response.statusCode ?? 0,
            headers: response.headers,
            body: Buffer.concat(chunks).toString("utf8"),
            byteCount,
          }),
        );
      },
    );
    const abort = () => request.destroy(new Error("Web request cancelled"));
    signal?.addEventListener("abort", abort, { once: true });
    request.setTimeout(timeoutMs, () =>
      request.destroy(new Error(`Web request timed out after ${timeoutMs} ms`)),
    );
    request.on("error", (error) => finish(reject, error));
    request.end();
  });
}

export async function requestPublicDocument(rawUrl, options = {}) {
  const timeoutMs = Math.min(Math.max(options.timeoutMs ?? DEFAULT_TIMEOUT_MS, 1_000), 30_000);
  const maxBytes = Math.min(Math.max(options.maxBytes ?? DEFAULT_MAX_BYTES, 1_024), 2_000_000);
  let current = normalizePublicUrl(rawUrl);
  const redirects = [];
  for (let redirectCount = 0; redirectCount <= MAX_REDIRECTS; redirectCount += 1) {
    const response = await requestOnce(current, {
      signal: options.signal,
      timeoutMs,
      maxBytes,
    });
    if ([301, 302, 303, 307, 308].includes(response.statusCode)) {
      const location = response.headers.location;
      if (!location) {
        throw new Error(`HTTP ${response.statusCode} redirect omitted Location`);
      }
      if (redirectCount === MAX_REDIRECTS) {
        throw new Error(`Web request exceeded ${MAX_REDIRECTS} redirects`);
      }
      redirects.push(current.toString());
      current = normalizePublicUrl(new URL(location, current).toString());
      continue;
    }
    if (response.statusCode < 200 || response.statusCode >= 300) {
      throw new Error(`Public web request returned HTTP ${response.statusCode}`);
    }
    return {
      requestedUrl: normalizePublicUrl(rawUrl).toString(),
      finalUrl: current.toString(),
      redirects,
      contentType: String(response.headers["content-type"] ?? ""),
      body: response.body,
      byteCount: response.byteCount,
    };
  }
  throw new Error("Unreachable redirect state");
}

function rssField(item, field) {
  const match = item.match(new RegExp(`<${field}\\b[^>]*>([\\s\\S]*?)<\\/${field}>`, "i"));
  if (!match) return "";
  return stripHtml(match[1].replace(/^\s*<!\[CDATA\[/, "").replace(/\]\]>\s*$/, ""));
}

export function parseBingRssResults(xml, maxResults = 6) {
  const limit = Math.min(Math.max(maxResults, 1), 10);
  const items = [...xml.matchAll(/<item\b[^>]*>([\s\S]*?)<\/item>/gi)];
  const results = [];
  const seen = new Set();
  for (const itemMatch of items) {
    if (results.length >= limit) break;
    const item = itemMatch[1];
    const title = rssField(item, "title").slice(0, MAX_TITLE_CHARS);
    const rawUrl = rssField(item, "link");
    if (!title || !rawUrl) continue;
    let url;
    try {
      url = normalizePublicUrl(rawUrl).toString();
    } catch {
      continue;
    }
    if (seen.has(url)) continue;
    seen.add(url);
    results.push({
      title,
      url,
      snippet: rssField(item, "description").slice(0, MAX_SNIPPET_CHARS),
    });
  }
  return results;
}

export function searchResultsMatchQuery(query, results) {
  const terms = [
    ...new Set(
      String(query)
        .normalize("NFKC")
        .toLowerCase()
        .match(/[\p{L}\p{N}][\p{L}\p{N}._+-]*/gu) ?? [],
    ),
  ].filter((term) => term.length >= 3 && !SEARCH_STOP_WORDS.has(term));
  if (terms.length === 0) return results.length > 0;
  const haystack = results
    .map(({ title, url, snippet }) => `${title} ${url} ${snippet}`)
    .join("\n")
    .normalize("NFKC")
    .toLowerCase();
  return terms.some((term) => haystack.includes(term));
}

export function extractPageText(body, contentType, maxChars = 20_000) {
  const limit = Math.min(Math.max(maxChars, 2_000), 50_000);
  const lowerType = contentType.toLowerCase();
  if (
    lowerType &&
    !lowerType.includes("text/") &&
    !lowerType.includes("html") &&
    !lowerType.includes("json") &&
    !lowerType.includes("xml")
  ) {
    throw new Error(`Unsupported public page content type '${contentType}'`);
  }
  const titleMatch = body.match(/<title\b[^>]*>([\s\S]*?)<\/title>/i);
  let text;
  if (lowerType.includes("json")) {
    try {
      text = JSON.stringify(JSON.parse(body), null, 2);
    } catch {
      text = body;
    }
  } else {
    text = decodeHtmlEntities(
      body
        .replace(/<script\b[\s\S]*?<\/script>/gi, " ")
        .replace(/<style\b[\s\S]*?<\/style>/gi, " ")
        .replace(/<svg\b[\s\S]*?<\/svg>/gi, " ")
        .replace(/<\/(?:p|div|section|article|main|header|footer|li|tr|h[1-6])>/gi, "\n")
        .replace(/<br\s*\/?\s*>/gi, "\n")
        .replace(/<[^>]*>/g, " "),
    )
      .replace(/[ \t]+/g, " ")
      .replace(/\n\s*\n\s*\n+/g, "\n\n")
      .trim();
  }
  return {
    title: titleMatch ? stripHtml(titleMatch[1]).slice(0, MAX_TITLE_CHARS) : null,
    text: text.slice(0, limit),
    includedCharCount: Math.min(text.length, limit),
    truncated: text.length > limit,
  };
}
