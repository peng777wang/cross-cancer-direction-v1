/* Functional test for pleiotropy_explorer.html.
 *
 * No browser is available in this environment, so this runs the explorer's real
 * script against a minimal DOM stub and asserts the behaviour that can actually
 * break: initial render, text search, region search, class filter, numeric sort,
 * row-click detail, and CSV export.
 *
 *   node test_explorer.js
 */
const fs = require("fs");
const path = require("path");

const html = fs.readFileSync(path.join(__dirname, "pleiotropy_explorer.html"), "utf8");
const payload = html.match(
  /<script id="payload" type="application\/json">([\s\S]*?)<\/script>/)[1]
  .replace(/<\\\//g, "</");
const DATA = JSON.parse(payload);
const code = html.match(/<script>([\s\S]*?)<\/script>\s*<\/body>/)[1];

let failures = 0;
function check(name, got, want) {
  const ok = got === want;
  if (!ok) { failures++; }
  console.log(`  [${ok ? "PASS" : "FAIL"}] ${name}: got ${got}, want ${want}`);
}
function checkTrue(name, cond, detail) {
  if (!cond) { failures++; }
  console.log(`  [${cond ? "PASS" : "FAIL"}] ${name}${detail ? " - " + detail : ""}`);
}

// ---------------------------------------------------------------- DOM stub
function mkEl(tag = "div", id = "") {
  const el = {
    tagName: tag.toUpperCase(), id, dataset: {}, style: {}, children: [],
    _listeners: {}, textContent: "", value: "", parent: null, _html: "",
    classList: {
      _s: new Set(),
      add(c) { this._s.add(c); },
      remove(c) { this._s.delete(c); },
      contains(c) { return this._s.has(c); },
    },
    get innerHTML() { return this._html; },
    set innerHTML(v) { this._html = v; this.children = parseChildren(v, this); },
    appendChild(c) { c.parent = this; this.children.push(c); return c; },
    addEventListener(ev, fn) { (this._listeners[ev] = this._listeners[ev] || []).push(fn); },
    fire(ev, target) {
      (this._listeners[ev] || []).forEach((f) => f({ target: target || this }));
    },
    querySelector(sel) {
      if (sel === ".arw") {
        return this.children.find((c) => c.classList.contains("arw")) || null;
      }
      return null;
    },
    querySelectorAll() { return this.children; },
    closest(sel) {
      let n = this;
      while (n) {
        if (sel === "tr" && n.tagName === "TR") { return n; }
        if (sel === "th" && n.tagName === "TH") { return n; }
        n = n.parent;
      }
      return null;
    },
    scrollIntoView() {},
    click() { this.fire("click", this); },
  };
  return el;
}

function parseChildren(str, parent) {
  const out = [];
  let m;
  const th = /<th data-k="([^"]+)"[^>]*>([\s\S]*?)<\/th>/g;
  while ((m = th.exec(str))) {
    const e = mkEl("th");
    e.dataset.k = m[1];
    e.parent = parent;
    const arw = mkEl("span");
    arw.classList.add("arw");
    arw.parent = e;
    e.children.push(arw);
    out.push(e);
  }
  const tr = /<tr data-i="(\d+)">/g;
  while ((m = tr.exec(str))) {
    const e = mkEl("tr");
    e.dataset.i = m[1];
    e.parent = parent;
    out.push(e);
  }
  return out;
}

const registry = {};
function el(id, tag) {
  if (!registry[id]) { registry[id] = mkEl(tag || "div", id); }
  return registry[id];
}
registry.payload = mkEl("script", "payload");
registry.payload.textContent = payload;

let csvText = null;
const document = {
  getElementById: (id) => el(id),
  createElement: (tag) => {
    const e = mkEl(tag);
    if (tag === "a") {
      Object.defineProperty(e, "href", { set() {}, get() { return ""; }, configurable: true });
      e.click = function () { csvText = lastBlob; };
    }
    return e;
  },
  querySelectorAll: (sel) => {
    if (sel === "#head th") { return el("head").children; }
    if (sel === "#body tr") { return el("body").children; }
    return [];
  },
};
let lastBlob = null;
function Blob(parts) { lastBlob = parts.join(""); }
const URL = { createObjectURL: () => "blob:stub" };

// ---------------------------------------------------------------- run it
console.log("running the explorer script against the DOM stub\n");
new Function("document", "Blob", "URL", code)(document, Blob, URL);

const nLoci = DATA.loci.length;
check("initial render shows every locus", el("body").children.length, nLoci);
checkTrue("summary chip reports the locus count",
  el("chips").innerHTML.includes(">" + nLoci + "<"));
checkTrue("caveats are rendered", el("caveats").innerHTML.length > 200);
checkTrue("definitions are rendered", el("defs").innerHTML.includes("borrowed"));
checkTrue("chromosome selector was populated",
  el("chr").children.length >= 20, el("chr").children.length + " options");

// text search by rsID (use a locus that is actually in the resource)
const someSnp = DATA.loci[0].snp;
el("q").value = someSnp;
el("q").fire("input");
check("search by rsID returns one locus", el("body").children.length, 1);

// gene search
el("q").value = "IREB2";
el("q").fire("input");
check("search by gene name", el("body").children.length,
  DATA.loci.filter((l) => (l.gene || "").toUpperCase().includes("IREB2")).length);

// region search
el("q").value = "15:78000000-79500000";
el("q").fire("input");
check("search by region", el("body").children.length,
  DATA.loci.filter((l) => l.chr === 15 && l.pos >= 78000000 && l.pos <= 79500000).length);

// cytoband search. 8q24 is NOT 24 Mb: in GRCh37 it spans 117,700,000-146,364,022
// (UCSC hg19 cytoBand). The expectation below is hard-coded from that table so the
// test is independent of the explorer's own band lookup.
el("q").value = "8q24";
el("q").fire("input");
const wantBand8q24 = DATA.loci.filter((l) => l.chr === 8 &&
  l.pos >= 117700000 && l.pos <= 146364022).length;
check("search by cytoband 8q24", el("body").children.length, wantBand8q24);
checkTrue("8q24 finds at least one locus", wantBand8q24 > 0, wantBand8q24 + " loci");

el("q").value = "15q25";
el("q").fire("input");
const wantBand15q25 = DATA.loci.filter((l) => l.chr === 15 &&
  l.pos >= 78300000 && l.pos <= 89100000).length;
check("search by cytoband 15q25", el("body").children.length, wantBand15q25);

el("q").value = "8q24.21";
el("q").fire("input");
const want821 = DATA.loci.filter((l) => l.chr === 8 &&
  l.pos >= 127300000 && l.pos <= 131500000).length;
check("search by cytoband 8q24.21", el("body").children.length, want821);

// class filter
el("q").value = "";
el("cls").value = "reversed";
el("cls").fire("change");
check("reversed filter", el("body").children.length,
  DATA.loci.filter((l) => l.cls.indexOf("reversed") === 0).length);

el("cls").value = "aligned";
el("cls").fire("change");
check("aligned filter", el("body").children.length,
  DATA.loci.filter((l) => l.cls === "aligned").length);
el("cls").value = "";
el("cls").fire("change");

// k filter
el("mink").value = "5";
el("mink").fire("change");
check("k>=5 filter", el("body").children.length,
  DATA.loci.filter((l) => l.k >= 5).length);
el("mink").value = "2";
el("mink").fire("change");
check("reset of the k filter", el("body").children.length, nLoci);

// numeric sorting
const head = el("head");
const thK = head.children.find((c) => c.dataset.k === "k");
head.fire("click", thK);
head.fire("click", thK);
const rowOrder = el("body")._html;
const firstK = Number(rowOrder.match(/<td class="num">(\d+)<\/td>/)[1]);
checkTrue("descending sort puts the largest k first", firstK === 6, "k = " + firstK);

// row click -> detail
const row = el("body").children[0];
el("body").fire("click", row);
const locus = DATA.loci[Number(row.dataset.i)];
checkTrue("row click renders the detail panel",
  el("detail").innerHTML.includes(locus.snp), "looking for " + locus.snp);
checkTrue("detail lists every trait",
  DATA.meta.traits.every((t) => el("detail").innerHTML.indexOf(">" + t) >= 0));

// reset + CSV export
el("reset").click();
check("reset restores the full table", el("body").children.length, nLoci);
el("dl").click();
const lines = (csvText || "").trim().split("\n");
check("CSV export has a header plus one row per locus", lines.length, nLoci + 1);
checkTrue("CSV header is the documented one",
  lines[0].startsWith("snp,gene,chr,pos,class_mtag,DSI_mtag"));

console.log("\n" + (failures ? failures + " FAILURES" : "all checks passed"));
process.exit(failures ? 1 : 0);
