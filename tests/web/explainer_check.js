// Injected into web/index.html by tests/web_smoke.py: runs the page's own functions and writes the numbers as JSON,
// which web_smoke.py compares with the Python code.
(() => {
  const res = {};
  try {
    res.fBell = F_BELL;
    res.esBell = expectedShortfall(F_BELL, MU, SIG, R, 1, ALPHA);
    res.ks = Object.fromEntries([1, 12, 252, 2520].map(n => [n, ksDistanceToGbm(0.1, 0.2, 1, n)]));
    res.fCrash = maxFractionUnderVar(f => valueAtRiskJump(f, MU, SIG, R, 1, 0.1, -0.25, ALPHA), 0.15, F_MAX).f;
    res.best = { w: BEST.w, g: BEST.g };
    res.trueBest = TRUE_BEST.g;
    res.grids = Object.fromEntries(SEARCH_STEPS.map(s => [s, searchGrid(s).best.g]));
    // every slider through its range (about 40 values each): nothing may throw or write NaN
    for (const el of document.querySelectorAll("input[type=range]")) {
      const lo = +el.min, hi = +el.max, step = Math.max(el.step === "any" ? 0 : +el.step, (hi - lo) / 40);
      for (let x = lo; x <= hi + 1e-9; x += step) { el.value = x; el.dispatchEvent(new Event("input")); }
    }
    $("run").click();                                       // the bisection animation (instant in tests)
    res.text = (document.body.innerText.match(/.{0,40}(NaN|undefined|Infinity).{0,40}/) || [""])[0];
  } catch (e) { res.exception = String(e.stack || e); }
  res.errors = window.__errors;
  document.getElementById("__out").textContent = JSON.stringify(res);
})();
