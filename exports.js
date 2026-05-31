"use strict";
/* exports.js — reproducibility & export layer for the Panel Data Multiverse Lab.
 * Generates runnable R / Stata / Python code for the current specification, exports the
 * multiverse as CSV, exports charts as SVG/PNG, and assembles a self-contained HTML lab report.
 * Plain global (window.Exports); no bundler required. */
(function (root) {

  // ---------- generic download ----------
  function download(filename, content, mime) {
    var blob = new Blob([content], { type: mime || "text/plain;charset=utf-8" });
    var url = URL.createObjectURL(blob), a = document.createElement("a");
    a.href = url; a.download = filename; document.body.appendChild(a); a.click();
    document.body.removeChild(a); setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  // ---------- shared spec → term lists ----------
  // controls present in the bundled NLSY extract as native columns vs. derived categoricals
  function rTerms(spec) {
    var c = spec.controls || {}, t = [];
    if (c.experience) t.push("exper", "expersq");
    if (c.hours) t.push("hours");
    if (c.health) t.push("poorhlth");
    if (c.race && spec.focal !== "educ") t.push("black", "hisp");
    if (c.industry) t.push("factor(industry)");
    if (c.region) t.push("factor(region)");
    if (spec.yearFE) t.push("factor(year)");
    return t;
  }
  function needsDerive(spec) { return (spec.controls && (spec.controls.industry || spec.controls.region)); }
  function outcomeExpr(spec, lang) {
    if (spec.outcome !== "wage") return "lwage";
    return lang === "r" ? "I(exp(lwage))" : lang === "py" ? "np.exp(df['lwage'])" : "wage";
  }
  function sampleComment(spec) {
    if (spec.sample === "nohealth") return "drop poor-health person-years";
    if (spec.sample === "trim") return "trim the outcome to its 1st–99th percentiles";
    return "full sample";
  }

  // ---------- R (plm) ----------
  function genR(spec) {
    var L = [], plmModel = { pooled: "pooling", between: "between", random: "random", within: "within", fd: "fd", twfe: "within", cre: "random" }[spec.estimator];
    var y = outcomeExpr(spec, "r");
    L.push("# " + label(spec) + "  —  " + sampleComment(spec));
    L.push("library(wooldridge); library(plm); library(lmtest); library(sandwich)");
    L.push("data(wagepan); d <- wagepan");
    if (needsDerive(spec)) {
      L.push("# derive the categorical controls used by the app");
      L.push('iv <- c("agric","min","construc","trad","tra","fin","bus","per","ent","manuf","pro","pub")');
      if (spec.controls.industry) L.push('d$industry <- factor(apply(d[,iv],1,function(r){i<-which(r==1);if(length(i))iv[i[1]] else "other"}))');
      if (spec.controls.region) L.push('d$region <- factor(with(d, ifelse(nrtheast==1,"northeast",ifelse(nrthcen==1,"northcentral",ifelse(south==1,"south","west")))))');
    }
    if (spec.sample === "nohealth") L.push("d <- subset(d, poorhlth != 1)");
    if (spec.sample === "trim") L.push("qq <- quantile(d$lwage, c(.01,.99)); d <- subset(d, lwage >= qq[1] & lwage <= qq[2])");
    L.push('pd <- pdata.frame(d, index = c("nr","year"))');
    var rhs = [spec.focal].concat(rTerms(spec)).join(" + ");
    if (spec.estimator === "cre") {
      L.push("# Correlated random effects (Mundlak): add person-means of the time-varying regressors");
      var tv = [spec.focal].concat(rTerms(spec).filter(function (t) { return t.indexOf("factor(") !== 0; }));
      tv.forEach(function (v) { L.push('d$' + v + '_m <- ave(d$' + v + ', d$nr)'); });
      L.push('pd <- pdata.frame(d, index = c("nr","year"))');
      rhs = rhs + " + " + tv.map(function (v) { return v + "_m"; }).join(" + ");
      L.push('m <- plm(' + y + ' ~ ' + rhs + ', data = pd, model = "random")');
    } else if (spec.estimator === "twfe") {
      L.push('m <- plm(' + y + ' ~ ' + rhs + ', data = pd, model = "within", effect = "twoways")');
    } else {
      L.push('m <- plm(' + y + ' ~ ' + rhs + ', data = pd, model = "' + plmModel + '")');
    }
    if (spec.seType === "cluster") L.push('coeftest(m, vcov = vcovHC(m, type = "HC1", cluster = "group"))   # cluster-robust by person');
    else if (spec.seType === "robust") L.push('coeftest(m, vcov = vcovHC(m, type = "HC1"))   # heteroskedasticity-robust');
    else L.push("summary(m)   # classical SEs");
    return L.join("\n");
  }

  // ---------- Stata ----------
  function genStata(spec) {
    var L = [], y = spec.outcome === "wage" ? "wage" : "lwage";
    var vce = spec.seType === "cluster" ? ", vce(cluster nr)" : spec.seType === "robust" ? ", vce(robust)" : "";
    L.push("* " + label(spec) + "  —  " + sampleComment(spec));
    L.push("* (load the wagepan extract, then:)");
    L.push("xtset nr year");
    if (needsDerive(spec)) L.push("* industry/region must be encoded as categoricals first (see app data-prep)");
    if (spec.sample === "nohealth") L.push("drop if poorhlth == 1");
    if (spec.sample === "trim") { L.push("_pctile lwage, p(1 99)"); L.push("keep if inrange(lwage, r(r1), r(r2))"); }
    var rhs = stataTerms(spec);
    if (spec.outcome === "wage") L.push("gen wage = exp(lwage)");
    if (spec.estimator === "pooled") L.push("reg " + y + " " + rhs + vce);
    else if (spec.estimator === "between") L.push("xtreg " + y + " " + rhs + ", be");
    else if (spec.estimator === "random") L.push("xtreg " + y + " " + rhs + ", re" + (vce ? " " + vce.slice(2) : ""));
    else if (spec.estimator === "within") L.push("xtreg " + y + " " + rhs + ", fe" + (vce ? " " + vce.slice(2) : ""));
    else if (spec.estimator === "twfe") L.push("xtreg " + y + " " + rhs + " i.year, fe" + (vce ? " " + vce.slice(2) : ""));
    else if (spec.estimator === "fd") L.push("reg D.(" + y + " " + rhs.replace(/i\.year/g, "") + ")" + vce + "   // first differences with drift");
    else if (spec.estimator === "cre") {
      L.push("* Mundlak correlated RE: add within-person means");
      L.push("foreach v of varlist " + rhs.replace(/i\.\w+/g, "").trim() + " { ");
      L.push("    bysort nr: egen mean_`v' = mean(`v')");
      L.push("}");
      L.push("xtreg " + y + " " + rhs + " mean_*, re" + (vce ? " " + vce.slice(2) : ""));
    }
    return L.join("\n");
  }
  function stataTerms(spec) {
    var c = spec.controls || {}, t = [spec.focal];
    if (c.experience) t.push("exper", "expersq");
    if (c.hours) t.push("hours");
    if (c.health) t.push("poorhlth");
    if (c.race && spec.focal !== "educ") t.push("black", "hisp");
    if (c.industry) t.push("i.industry");
    if (c.region) t.push("i.region");
    if (spec.yearFE) t.push("i.year");
    return t.join(" ");
  }

  // ---------- Python (linearmodels) ----------
  function genPython(spec) {
    var L = [], y = spec.outcome === "wage" ? "np.exp(df['lwage'])" : "df['lwage']";
    L.push("# " + label(spec) + "  —  " + sampleComment(spec));
    L.push("import numpy as np, pandas as pd");
    L.push("from linearmodels.panel import PooledOLS, BetweenOLS, RandomEffects, PanelOLS, FirstDifferenceOLS");
    L.push("# df: the wagepan extract with columns nr, year, lwage, union, married, educ, exper, expersq, ...");
    if (spec.sample === "nohealth") L.push("df = df[df['poorhlth'] != 1].copy()");
    if (spec.sample === "trim") L.push("lo, hi = np.quantile(df['lwage'], [.01, .99]); df = df[(df['lwage']>=lo)&(df['lwage']<=hi)].copy()");
    L.push("panel = df.set_index(['nr','year'])");
    var cols = pyCols(spec);
    L.push("y = panel['lwage']" + (spec.outcome === "wage" ? ".pipe(np.exp)" : ""));
    L.push("X = panel[" + JSON.stringify(cols) + "]");
    if (spec.controls && (spec.controls.industry || spec.controls.region)) L.push("# add pd.get_dummies(panel['industry']/['region'], drop_first=True) to X as needed");
    var addc = "import statsmodels.api as sm; X = sm.add_constant(X)";
    if (spec.estimator === "pooled") { L.push(addc); L.push("res = PooledOLS(y, X).fit(" + clusterArg(spec) + ")"); }
    else if (spec.estimator === "between") { L.push(addc); L.push("res = BetweenOLS(y, X).fit()"); }
    else if (spec.estimator === "random") { L.push(addc); L.push("res = RandomEffects(y, X).fit(" + clusterArg(spec) + ")"); }
    else if (spec.estimator === "within") L.push("res = PanelOLS(y, X, entity_effects=True, drop_absorbed=True).fit(" + clusterArg(spec) + ")");
    else if (spec.estimator === "twfe") L.push("res = PanelOLS(y, X, entity_effects=True, time_effects=True, drop_absorbed=True).fit(" + clusterArg(spec) + ")");
    else if (spec.estimator === "fd") L.push("res = FirstDifferenceOLS(y, X).fit(" + clusterArg(spec) + ")   # note: linearmodels FD omits the drift term");
    else if (spec.estimator === "cre") { L.push("for v in " + JSON.stringify(cols) + ": X[v+'_m'] = panel.groupby(level=0)[v].transform('mean')  # Mundlak means"); L.push(addc); L.push("res = RandomEffects(y, X).fit(" + clusterArg(spec) + ")"); }
    L.push("print(res)");
    return L.join("\n");
  }
  function pyCols(spec) {
    var c = spec.controls || {}, t = [spec.focal];
    if (c.experience) t.push("exper", "expersq");
    if (c.hours) t.push("hours");
    if (c.health) t.push("poorhlth");
    if (c.race && spec.focal !== "educ") t.push("black", "hisp");
    return t;
  }
  function clusterArg(spec) { return spec.seType === "cluster" ? "cov_type='clustered', cluster_entity=True" : spec.seType === "robust" ? "cov_type='robust'" : "cov_type='unadjusted'"; }

  function label(spec) {
    var est = { pooled: "Pooled OLS", between: "Between", random: "Random effects", within: "Fixed effects", fd: "First differences", cre: "Correlated RE (Mundlak)", twfe: "Two-way FE" }[spec.estimator];
    var foc = { union: "union premium", married: "marriage premium", educ: "return to schooling" }[spec.focal];
    return est + " — " + foc;
  }

  // ---------- multiverse CSV ----------
  function multiverseCSV(results) {
    var head = ["focal", "outcome", "estimator", "experience", "hours", "industry", "region", "health", "yearFE", "sample", "seType", "identified", "estimate", "std_error", "ci_low", "ci_high"];
    var lines = [head.join(",")];
    results.forEach(function (r) {
      var s = r.spec, c = s.controls || {}, lo = r.b - 1.96 * r.se, hi = r.b + 1.96 * r.se;
      lines.push([s.focal, s.outcome, r.estimator, +!!c.experience, +!!c.hours, +!!c.industry, +!!c.region, +!!c.health, +!!s.yearFE, s.sample, s.seType,
        r.identified ? 1 : 0, r.identified ? r.b.toFixed(6) : "NA", r.identified && isFinite(r.se) ? r.se.toFixed(6) : "NA",
        r.identified && isFinite(r.se) ? lo.toFixed(6) : "NA", r.identified && isFinite(r.se) ? hi.toFixed(6) : "NA"].join(","));
    });
    return lines.join("\n");
  }

  // ---------- figure export ----------
  function chartSVGString(containerId) {
    var svg = document.querySelector("#" + containerId + " svg");
    if (!svg) return null;
    var clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.style.background = "#0c0c0c";
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
  }
  function exportSVG(containerId, filename) { var s = chartSVGString(containerId); if (s) download(filename, s, "image/svg+xml"); }
  function exportPNG(containerId, filename, scale) {
    var svg = document.querySelector("#" + containerId + " svg"); if (!svg) return;
    var vb = svg.getAttribute("viewBox").split(/\s+/).map(Number), w = vb[2], h = vb[3]; scale = scale || 2;
    var img = new Image();
    img.onload = function () {
      var canvas = document.createElement("canvas"); canvas.width = w * scale; canvas.height = h * scale;
      var ctx = canvas.getContext("2d"); ctx.fillStyle = "#0c0c0c"; ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      canvas.toBlob(function (blob) {
        var url = URL.createObjectURL(blob), a = document.createElement("a");
        a.href = url; a.download = filename; document.body.appendChild(a); a.click(); document.body.removeChild(a);
        setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
      });
    };
    img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(chartSVGString(containerId))));
  }

  root.Exports = {
    download: download, genR: genR, genStata: genStata, genPython: genPython, label: label,
    multiverseCSV: multiverseCSV, exportSVG: exportSVG, exportPNG: exportPNG, chartSVGString: chartSVGString
  };
})(window);
