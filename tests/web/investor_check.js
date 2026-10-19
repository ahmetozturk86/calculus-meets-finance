// Injected into web/investor.html by tests/web_smoke.py (which also defines SAMPLE_CSV): clicks through the page like a
// visitor and writes what it finds as JSON, which web_smoke.py compares with the Python code.
(async () => {
  const res = {}, setv = (id, v) => { $(id).value = v; $(id).dispatchEvent(new Event("input")); };
  const bad = () => (document.body.innerText.match(/.{0,40}(NaN|undefined|\[object).{0,40}/) || [""])[0];
  const plan = async () => { calculate(); const done = await life.many; return { w: cur.plan.w, g: cur.plan.g, var: cur.plan.var, L: cur.p.L, done }; };
  try {
    res.default = await plan();
    const s = mixStats(cur.plan.w, cur.m, cur.p.alpha);
    res.mixStats = [s.g - cur.plan.g, s.var - cur.plan.var, s.es - cur.plan.es];
    res.sentence = $("planSentence").textContent;
    setv("r", 100 * (Math.exp(0.03) - 1)); setv("mu", 100 * (Math.exp(0.1) - 1));   // Python's continuous 3% and 10%
    res.pythonRates = await plan();
    $("resetParams").click();
    document.querySelector("[data-risk=careful]").click(); res.careful = await plan();
    document.querySelector("[data-risk=bold]").click(); res.bold = await plan();
    document.querySelector("[data-risk=balanced]").click();
    setLang("tr"); res.turkish = bad(); setLang("en");
    setv("currency", "₺"); res.lira = await plan(); setv("currency", "$"); $("resetParams").click();
    // a savings account that loses money after fees, so even all savings breaks the Careful limit
    setv("feeDep", 10); document.querySelector("[data-risk=careful]").click(); res.feeEdge = await plan();
    setv("feeDep", ""); document.querySelector("[data-risk=balanced]").click();
    // file and column names are shown as text, never run as HTML
    window.__ran = [];
    const rows = Array.from({ length: 40 }, (_, i) => `2020-${String(1 + i % 12).padStart(2, "0")}-${String(1 + Math.floor(i / 12)).padStart(2, "0")},${100 + i}`).join("\n");
    loadCsvText("Date,<img src=x onerror=\"__ran.push('column')\">\n" + rows, "a.csv", 1);
    loadCsvText("Date,Close\n" + rows, "<img src=x onerror=\"__ran.push('file')\">.csv", 2);
    await plan();
    await new Promise(r => setTimeout(r, 300));
    res.injected = window.__ran;
    $("clear1").click(); $("clear2").click();
    // the guess game: the plan rounded to whole percents (73 / 4 / 23) is a hair over the limit
    $("useGuess").checked = true; $("useGuess").dispatchEvent(new Event("input"));
    $("gStock").value = 73; $("gGold").value = 4; $("gBond").value = 23;
    await plan();
    const nv = nearValid(cur.guess.w, cur.m, cur.p.alpha, cur.p.L);
    res.guess = { var: cur.guess.var, L: cur.p.L, verdict: $("vsVerdict").textContent, game: $("lifeGame").textContent,
                  nearVar: nv ? mixStats(nv.w, cur.m, cur.p.alpha).var : null };
    $("useGuess").checked = false; $("useGuess").dispatchEvent(new Event("input"));
    // real prices
    loadCsvText(SAMPLE_CSV, "sample_prices.csv", 0);
    res.csv = { status: $("dataStatus0").textContent, plan: await plan(), history: !$("histBody").classList.contains("hidden") };
    $("clear0").click();
    // Fund Manager: five years, holding through every drop, then the end screen's "new best!"
    arcOpen(); arcStart();
    for (let y = 0; y < ARC_YEARS && arc.phase === "plan"; y++) {
      arcPlay(); renderArcade(); arcRun();
      while (arc.phase === "play" && arc.event) arcChoose("hold");
      if (arc.phase === "year") { arcNext(); renderArcade(); }
    }
    res.arcade = arc.phase;
    localStorage.removeItem(ARC_BEST);
    Object.assign(arc, { phase: "end", rank: 2, money: 1234567.7, quant: 1e6, bank: 1e6, hist: [] });
    arcSaveBest(); renderArcade();
    res.arcadeBest = /new best/i.test($("arcScreen").textContent);
    arcClose();
    res.text = bad();
  } catch (e) { res.exception = String(e.stack || e); }
  res.errors = window.__errors;
  document.getElementById("__out").textContent = JSON.stringify(res);
})();
