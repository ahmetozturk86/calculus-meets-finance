"""Smoke test for the two web pages, in a headless browser (Chrome, Chromium or Edge).

Each page is copied with a check script injected (tests/web/*.js) that runs the page's own JavaScript and writes its
numbers as JSON. This file compares them with the Python code (src/) and checks that no JavaScript error happened.
The pages' normal CDF uses an approximation of erf, so the numbers agree to about 6 digits.

    python tests/web_smoke.py
"""
import html
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.gbm import ks_distance_to_gbm  # noqa: E402
from src.portfolio import (best_mix_on_grid, best_mix_under_var, portfolio_growth, portfolio_var,  # noqa: E402
                           refine_best_mix)
from src.risk import expected_shortfall, max_fraction_under_var, max_fraction_under_var_jump  # noqa: E402

MU, SIGMA, R, ALPHA, L = 0.10, 0.20, 0.03, 0.95, 0.15
MU3 = np.array([MU, math.log(1.05), math.log(1.045)])
SIG3 = np.array([SIGMA, 0.15, 0.06])
RHO3 = np.array([[1.0, 0.0, -0.2], [0.0, 1.0, 0.3], [-0.2, 0.3, 1.0]])
# The Investor's example numbers are ordinary yearly rates (stocks 10.5%, bank 3%), so ln(1 + x) of them
PAGE_MU, PAGE_R = np.array([math.log(1.105), math.log(1.05), math.log(1.045)]), math.log(1.03)

# Put before the page's own code: catch every error, switch off the animations, start without saved games.
HEAD = """<script>
window.__errors = [];
addEventListener("error", e => __errors.push("error: " + e.message + " (line " + e.lineno + ")"));
addEventListener("unhandledrejection", e => __errors.push("rejection: " + e.reason));
window.INVESTOR_NO_ANIM = true; window.EXPLAINER_NO_ANIM = true;
try { localStorage.clear(); } catch (e) {}
</script>"""


def find_browser():
    """A Chromium-based browser: on PATH (Linux, CI) or in its usual place (Windows, macOS)."""
    for name in ["google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "chrome", "msedge"]:
        if shutil.which(name):
            return shutil.which(name)
    for path in [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
                 r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
                 r"C:\Program Files\Google\Chrome\Application\chrome.exe",
                 "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                 "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge"]:
        if os.path.exists(path):
            return path
    return None


def run_page(browser, page, script, tmp, before=""):
    """Open web/<page> with tests/web/<script> injected and return the JSON it writes."""
    text = (ROOT / "web" / page).read_text(encoding="utf-8").replace("<head>", "<head>\n" + HEAD, 1)
    end = text.rindex("</body>")
    js = before + (ROOT / "tests" / "web" / script).read_text(encoding="utf-8")
    path = tmp / page
    path.write_text(text[:end] + "<pre id='__out'></pre>\n<script>\n" + js + "\n</script>\n" + text[end:], encoding="utf-8")
    args = [browser, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check", "--disable-extensions",
            f"--user-data-dir={tmp / ('profile-' + path.stem)}", "--virtual-time-budget=600000", "--dump-dom", path.as_uri()]
    if platform.system() == "Linux":
        args.insert(1, "--no-sandbox")                      # CI runners and containers
    dom = subprocess.run(args, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=900).stdout
    found = re.search(r"<pre id=\"__out\">(.*?)</pre>", dom, re.S)
    if not found or not found.group(1).strip():
        raise RuntimeError(f"{page}: the check script wrote nothing (did the page crash?)")
    return json.loads(html.unescape(found.group(1)))


failures = []


def check(name, ok, detail=""):
    print(("ok    " if ok else "FAIL  ") + name + (f"  ({detail})" if detail and not ok else ""))
    if not ok:
        failures.append(name)


def close(a, b, tol):
    return a is not None and abs(a - b) <= tol


def check_explainer(res):
    check("explainer: no JavaScript errors", not res.get("exception") and not res["errors"], res.get("exception") or res["errors"])
    if res.get("exception"):
        return
    f_var = max_fraction_under_var(L, MU, SIGMA, R, 1.0, ALPHA)[0]
    check("explainer: f_VaR = 0.705014", close(res["fBell"], f_var, 1e-5), res["fBell"])
    check("explainer: ES at f_VaR", close(res["esBell"], float(expected_shortfall(f_var, MU, SIGMA, R, 1.0, ALPHA)), 1e-5), res["esBell"])
    for n, v in res["ks"].items():
        check(f"explainer: KS distance, n = {n}", close(v, ks_distance_to_gbm(MU, SIGMA, 1.0, int(n)), 1e-5), v)
    f_crash = max_fraction_under_var_jump(L, MU, SIGMA, R, 1.0, 0.1, -0.25, ALPHA)[0]
    check("explainer: f_crash = 0.617447", close(res["fCrash"], f_crash, 1e-5), res["fCrash"])
    w, _, _ = best_mix_under_var(L, MU3, SIG3, RHO3, R, 1.0, ALPHA)
    check("explainer: best mix 73 / 4 / 23", np.allclose(res["best"]["w"][:3], w, atol=2e-5), res["best"]["w"])
    check("explainer: its growth", close(res["best"]["g"], portfolio_growth(w, MU3, SIG3, RHO3, R)[0], 2e-6), res["best"]["g"])
    for s, g in res["grids"].items():
        check(f"explainer: best growth on the grid, step 1/{s}", close(g, best_mix_on_grid(L, MU3, SIG3, RHO3, R, 1.0, ALPHA, int(s))[3], 2e-6), g)
    check("explainer: true best (step -> 0)", close(res["trueBest"], refine_best_mix(L, MU3, SIG3, RHO3, R, 1.0, ALPHA)[3], 2e-6), res["trueBest"])
    check("explainer: every slider, no NaN on the page", not res["text"], res["text"])


def check_investor(res):
    check("investor: no JavaScript errors", not res.get("exception") and not res["errors"], res.get("exception") or res["errors"])
    if res.get("exception"):
        return

    def same_plan(name, got, L_, mu, r):
        w, _, _ = best_mix_under_var(L_, mu, SIG3, RHO3, r, 1.0, ALPHA)
        check(f"investor: {name} plan", np.allclose(got["w"][:3], w, atol=2e-5), got["w"])
        check(f"investor: {name} growth", close(got["g"], portfolio_growth(w, mu, SIG3, RHO3, r)[0], 2e-6), got["g"])
        check(f"investor: {name}, 500 futures finished", got["done"] is True)

    same_plan("default (73 / 4 / 23)", res["default"], L, PAGE_MU, PAGE_R)
    check("investor: mixStats(plan) = plan", all(abs(x) < 1e-12 for x in res["mixStats"]), res["mixStats"])
    check("investor: the README's plan sentence", res["sentence"].startswith("Put $7,298 in stocks, $400 in gold and $2,299 in bonds"),
          res["sentence"][:80])
    same_plan("Python's rates", res["pythonRates"], L, MU3, R)
    same_plan("Careful (5%)", res["careful"], 0.05, PAGE_MU, PAGE_R)
    top, _, _ = best_mix_under_var(1.0, PAGE_MU, SIG3, RHO3, PAGE_R, 1.0, ALPHA)   # no limit: the fastest-growing mix
    check("investor: Bold's limit is the top's bad-year loss",
          close(res["bold"]["L"], max(portfolio_var(top, PAGE_MU, SIG3, RHO3, PAGE_R, 1.0, ALPHA), 0.02), 2e-5), res["bold"]["L"])
    check("investor: Bold growth", close(res["bold"]["g"], portfolio_growth(top, PAGE_MU, SIG3, RHO3, PAGE_R)[0], 2e-6), res["bold"]["g"])
    check("investor: Turkish page, no NaN or undefined", not res["turkish"], res["turkish"])
    check("investor: lira plan", res["lira"]["done"] is True)
    check("investor: a savings account that loses money after fees", res["feeEdge"]["done"] is True)
    check("investor: file and column names are not run as HTML", res["injected"] == [], res["injected"])
    g = res["guess"]
    check("investor: 73 / 4 / 23 is just over the limit", g["L"] < g["var"] < g["L"] + 0.0005, g["var"])
    check("investor: ... and says so with 3 decimals", f"{100 * g['var']:.3f}%" in g["verdict"] and "72.98% stocks" in g["verdict"], g["verdict"])
    check("investor: ... and suggests a split that keeps the limit", g["nearVar"] is not None and g["nearVar"] <= g["L"] + 1e-9, g["nearVar"])
    check("investor: ... and the game says it is just over", "Just over" in g["game"], g["game"][:120])
    check("investor: sample_prices.csv loads", res["csv"]["status"].startswith("✓ Loaded sample_prices.csv"), res["csv"]["status"])
    check("investor: plan and history check with real prices", res["csv"]["plan"]["done"] is True and res["csv"]["history"])
    check("investor: Fund Manager plays five years", res["arcade"] in ("end", "fired"), res["arcade"])
    check("investor: Fund Manager shows a new best", res["arcadeBest"] is True)
    check("investor: no NaN or undefined on the page", not res["text"], res["text"])


def main():
    sys.stdout.reconfigure(errors="replace")                # page text in a Windows console
    browser = find_browser()
    if not browser:
        sys.exit("No Chrome, Chromium or Edge found: install one, or put it on PATH.")
    print("browser:", browser)
    tmp = Path(tempfile.mkdtemp(prefix="web-smoke-"))
    try:
        check_explainer(run_page(browser, "index.html", "explainer_check.js", tmp))
        sample = (ROOT / "web" / "sample_prices.csv").read_text(encoding="utf-8")
        check_investor(run_page(browser, "investor.html", "investor_check.js", tmp, f"const SAMPLE_CSV = {json.dumps(sample)};\n"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{len(failures)} failed" if failures else "\nall web checks passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
