"""
PART 1: Validate against PROSPER (Beggs & Brill)

What is compared (Beggs & Brill in PROSPER, same inputs for both)
    * PVT       : Pb, Bo, oil density, Rs, oil viscosity, gas density, Z
    * IPR       : P-bottomhole from the IPR curve
    * VLP       : P-bottomhole at every rate in the System table, with split
                  into friction and gravity
    * Traverse  : P vs. Depth at the Gradient Traverse rate
    * Solution  : Operating point (or 'no natural flow')
    

Temperature: PROSPER predicts a diff. temperature at every rate (80F at low
rate, 125-136F at high rate). The engine gets the same wellhead temperature
at each rate (read from the System table) and the shape of the temp-vs-depth
profile from the Gradient Traverse. Bottomhole temp is the reservoir tempt.

"""

import argparse
import glob
import os
from dataclasses import replace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

import Part1_Nodal_Analysis_Engine as nae

# Conversion unit
STB, PSI, ATM, FT, LB = 6.28981, 14.5038, 14.696, 0.3048, 16.0185    # STB/sm3, psi/bar, psi, m/ft


#   === readers ===

def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _lines(path):
    return open(path, encoding="latin-1").read().replace("\r", "").split("\n")


def classify(path):
    text = "\n".join(_lines(path)[:30])
    if "SYSTEM SENSITIVITY ANALYSIS" in text:
        return "system"
    if "Gradient Traverse" in text:
        return "gradient"
    if "PVT - Calculation Results" in text:
        return "pvt"
    return None


def input_value(lines, label):
    for ln in lines:
        t = ln.split("\t")
        if t[0].strip() == label and len(t) > 1:
            return _num(t[1].strip())
    return None


def read_system(path):
    lines = _lines(path)
    rows = []
    for ln in lines:
        t = [x.strip() for x in ln.split("\t")]
        v = [_num(x) for x in t[:17]]
        if len(t) >= 17 and all(a is not None for a in v):
            rows.append(v)
    a = np.array(rows)
    sol = {}
    for key, label in (("rate", "Liquid Rate"), ("bhp", "Solution Node Pressure"),
                       ("fric", "dP Friction"),("grav", "dP Gravity")):
        after = False
        for ln in lines:
            if "Solution Point" in ln:
                after = True
            if after and ln.split("\t")[0].strip() == label:
                t = ln.split("\t")
                sol[key] = _num(t[1].strip()) if len(t) > 1 else None
                break
    return dict(wc=input_value(lines, "Water Cut"), whp=input_value(lines, "Top Node Pressure"),
                rate=a[:, 0], vlp=a[:, 4], ipr=a[:, 5], wht=a[:, 13], dpf=a[:, 15], dpg=a[:, 16],
                solution=sol)
    
    
def read_gradient(path):
    lines = _lines(path)
    rows = []
    for ln in lines:
        t = [x.strip() for x in ln.split("\t")]
        if len(t) >= 18 and _num(t[1]) is not None and _num(t[3]) is not None:
            rows.append(dict(md=_num(t[1]), p=_num(t[3]), t=_num(t[4]), hl=_num(t[6]),
                             fric=_num(t[16]), grav=_num(t[17])))
    return dict(wc=input_value(lines, "Water Cut"), rate=input_value(lines, "Liquid Rate"),
                rows=[r for r in rows if r["md"] <= 8857.0 + 1E-06])
    
    
def read_pvt(path):
    seen, rows = set(), []
    for ln in _lines(path):
        t = [x.strip() for x in ln.split("\t")]
        v = [_num(x) for x in t[:12]]
        if len(t) >= 12 and all(a is not None for a in v) and (v[0], v[1]) not in seen:
            seen.add((v[0], v[1]))
            rows.append(v)
    return np.array(rows)



#   === Engine side ===

def well_with_temp(w0, md_ft, t_f, twh_f, wc):
    """Temp profile: PROSPER's shape, scaled to the wellhead temp of this rate."""
    md_ft, t_f      = np.array(md_ft), np.array(t_f)
    shape           = (t_f[-1] - t_f) / (t_f[-1] - t_f[0])
    tz              = t_f[-1] - (t_f[-1] - twh_f) * shape
    return replace(w0, wc=wc / 100.0, t_profile=(list(md_ft * FT), [(t - 32) / 1.8 for t in tz]),
                   tvd_m=8858 * FT, n_seg=120)
    
    
def pct(a, b):
    return (np.asarray(a, float) - np.asarray(b, float)) / np.asarray(b, float) * 100.0


def main():
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("--props", default="CASE_PROPS.INC")
    ap.add_argument("--prosper-dir", default="PROSPER_exports")
    args, unknown = ap.parse_known_args()
    
    sysd, grad, pvts = {}, {}, []
    for p in sorted(glob.glob(os.path.join(args.prosper_dir, "*.[Pp][Rr][Nn]"))):
        kind = classify(p)
        if kind == "system":
            d = read_system(p); sysd[round(d["wc"])] = d
        elif kind == "gradient":
            d = read_gradient(p); grad[round(d["wc"])] = d
        elif kind == "pvt":
            pvts.append(read_pvt(p))
    print(f"Found: {len(sysd)} System, {len(grad)} Gradient Traverse, {len(pvts)} PVT export(s) in {os.path.abspath(args.prosper_dir)}")
    if not (sysd or grad or pvts):
        raise SystemExit("No PROSPER .PRN exports found: set --prosper-dir (or its default) to the folder that holds them.")
    
    w0 = nae.Well()
    f = nae.Fluid(args.props, w0.t_res_c)
    summary = {}
    
    
    # ------- PVT -------
    if pvts:
        pv = np.vstack(pvts)
        out = {k: [] for k in ("pb", "bo", "rho", "rs", "mu_hot", "mu_cold", "rhog", "z")}
        for t_f, p_psig, pb, rs, rho, mu, bo, _, rhog, _, _, z in pv:
            tc, pbar = (t_f - 32) / 1.8, (p_psig + ATM) / PSI
            r = f.rs(pbar, tc); b = f.bo(pbar, tc)
            out["pb"].append(pct(f.pb_at(tc) * PSI - ATM, pb))
            out["bo"].append(pct(b, bo))
            out["rho"].append(pct((f.rho_o_sc + r * f.rho_g_sc) / b / LB, rho))
            out["rs"].append(r * 5.614584 - rs)
            out["mu_hot" if t_f > 120 else "mu_cold"].append(pct(f.mu_o(pbar, tc), mu))
            out["rhog"].append(pct(f.rho_g_sc / f.bg(pbar, tc) / LB, rhog))
            out["z"].append(pct(f.z(pbar, tc), z))
        mx = lambda k: float(np.max(np.abs(out[k]))) if out[k] else float("nan")
        summary["pvt"] = dict(pb=mx("pb"), bo=mx("bo"), rho=mx("rho"), rs=mx("rs"),
                              rhog=mx("rhog"), z=mx("z"), mu_hot=mx("mu_hot"),
                              mu_cold=(float(np.min(out["mu_cold"])),
                              float(np.max(out["mu_cold"]))) if out["mu_cold"] else None)
        s = summary["pvt"]
        print("\nPVT (script vs. PROSPER, all temperature and pressures in the exports)")
        print(f"    bubble point {s['pb']:.2f}% | Bo {s['bo']:.2f}% | oil density {s['rho']:.2f}% | "
              f"Rs {s['rs']:.1f} scf/STB | gas density {s['rhog']:.2f}% | Z {s['z']:.2f}%")
        print(f"    oil visocity: hot (>120 F) within {s['mu_hot']:.1f}%"
              + (f", cold (<120 F) script is {s['mu_cold'][0]:+.0f}% to {s['mu_cold'][1]:.0f}% of PROSPER" if s['mu_cold'] else ""))
        
        
    
    # ------- IPR, VLP, solution, traverse -------
    cases = sorted(set(sysd) & set(grad))
    fig, ax = plt.subplots(max(len(cases), 1), 3, figsize=(16, 4.6 * max(len(cases), 1)), squeeze=False)
    for row, wc in enumerate(cases):
        S, G = sysd[wc], grad[wc]
        md = [r["md"] for r in G["rows"]]
        tf = [r["t"] for r in G["rows"]]
        bhp, fr, gr = [], [], []
        for i, q in enumerate(S["rate"]):
            w = well_with_temp(w0, md, tf, S["wht"][i], wc)
            det = nae.traverse(q / STB, w, f, detail=True)
            bhp.append(det[-1][1] * PSI - ATM); fr.append(det[-1][4] * PSI); gr.append(det[-1][5] * PSI)
        bhp, d_vlp = np.array(bhp), pct(bhp, S["vlp"])
        
        
        # IPR: above Pb (the part that decides the operating point), below Pb
        # near AOF the gauge/abs pressure basis differs, so PROSPER pressure
        # under 100 psig are skipped
        ipr = np.array([nae.pwf_ipr(q / STB, replace(w0, wc=wc / 100.0), f.pb) * PSI - ATM for q in S["rate"]])
        pb_psig = f.pb * PSI - ATM
        hi, lo = S["ipr"] > pb_psig, (S["ipr"] > 100) & (S["ipr"] <= pb_psig)
        ipr_err = float(np.max(np.abs(ipr[hi] - S["ipr"][hi])))
        ipr_err_lo = float(np.max(np.abs(ipr[lo] - S["ipr"][lo]))) if lo.any() else float("nan")
        
        
        # traverse at the gradient rate
        w           = well_with_temp(w0, md, tf, tf[0], wc)
        det         = nae.traverse(G["rate"] / STB, w, f, detail=True)
        z_ft        = np.array(md)
        p_sc        = np.interp(z_ft * FT, [d[0] for d in det], [d[1] * PSI - ATM for d in det])
        p_pr        = np.array([r["p"] for r in G["rows"]])
        trav_bhp    = pct(p_sc[-1], p_pr[-1])
        
        
        # operating point, rate dependent temp
        def bhp_at(q):
            twh = float(np.interp(q, S["rate"], S["wht"]))
            return nae.traverse(q / STB, well_with_temp(w0, md, tf, twh, wc), f) * PSI - ATM
        
        
        def ipr_at(q):
            return nae.pwf_ipr(q / STB, replace(w0, wc=wc / 100.0), f.pb) * PSI - ATM
        
        
        qs = np.linspace(S["rate"][0], S["rate"][-1] * 0.99, 200)
        dd = np.array([ipr_at(q) - bhp_at(q) for q in qs])
        idx = np.where(np.diff(np.sign(dd)) != 0)[0]
        op = None
        if len(idx):
            i   = idx[-1]
            q   = brentq(lambda x: ipr_at(x) - bhp_at(x), qs[i], qs[i + 1])
            op  = (q, ipr_at(q))
        sol = S["solution"]
        pr_op = (sol["rate"], sol["bhp"]) if sol.get("rate") else None
        
        summary[wc] = dict(vlp_min=d_vlp.min(), vlp_max=d_vlp.max(), vlp_mabs=np.abs(d_vlp).mean(),
                           trav=trav_bhp, ipr_err=ipr_err, ipr_err_lo=ipr_err_lo, fric_hi=(fr[-1]),
                           op=op, pr_op=pr_op)
        
        print(f"\n=== Water cut {wc}%   (PROSPER Beggs and Brill, WHP {S['whp']:.1f} psig) ===")
        print(" rate STB/d | PROSPER VLP | script VLP | diff % | friction PROSPER/script psi")
        for q, v, b, d, f1, f2 in zip(S["rate"], S["vlp"], bhp, d_vlp, S["dpf"], fr):
            print(f"  {q:9.1f}  | {v:9.1f}   | {b:9.1f}  | {d:+6.2f} | {f1:8.2f} / {f2:8.2f}")
        
        print(f"    VLP difference {d_vlp.min():+.2f}% to {d_vlp.max():+.2f}%, mean absolute {np.abs(d_vlp).mean():.2f}%")
        print(f"    IPR pressure: max difference {ipr_err:.1f} psi above Pb, {ipr_err_lo:.1f} psi below Pb (100 psig to Pb)")
        print(f"    Gradient traverse at {G['rate']:.0f} STB/d: BHP PROSPER {p_pr[-1]:.1f}, script {p_sc[-1]:.1f} psig ({trav_bhp:=.2f}%)")
        print(f"    Operating point: PROSPER " + (f"{pr_op[0]:.1f} STB/d, {pr_op[1]:.1f} psig" if pr_op else "none (VLP above IPR at every rate)")
              + " | script " + (f"{op[0]:.1f} STB/d, {op[1]:.1f} psig" if op else "none"))
            
            
        a0, a1, a2 = ax[row]
        a0.plot(S["rate"], S["vlp"], "ko-", label="PROSPER VLP (Beggs and Brill)")
        a0.plot(S["rate"], bhp, "r.-", label="Python VLP")
        a0.plot(S["rate"], S["ipr"], "g--", label="PROSPER IPR")
        a0.set_title(f"Water cut {wc}%: VLP and IPR"); a0.set_xlabel("Liquid rate (STB/d)")
        a0.set_ylabel("Bottom hole pressure (psig)"); a0.legend(fontsize=8); a0.grid(alpha=0.3)
        
        a1.bar(range(len(S["rate"])), d_vlp, color="tab:blue")
        a1.set_xticks(range(len(S["rate"])))
        a1.set_xticklabels([f"{r:.0f}" for r in S["rate"]], rotation=90, fontsize=7)
        a1.axhline(0, color="k", lw=0.8); a1.set_title("VLP: Python - PROSPER (%)")
        a1.set_xlabel("Liquid rate (STB/d)"); a1.grid(alpha=0.3)
        
        a2.plot(p_pr, z_ft, "ko-", ms=3, label="PROSPER"); a2.plot(p_sc, z_ft, "r--", label="Python")
        a2.invert_yaxis(); a2.set_title(f"Gradient traverse at {G['rate']:.0f} STB/d")
        a2.set_xlabel("Pressure (psig)"); a2.set_ylabel("Depth (ft)"); a2.legend(); a2.grid(alpha=0.3)
            
    fig.tight_layout()
    fig.savefig("validation_vs_prosper.png", dpi=150)
    print("\nSaved validation_vs_prosper.png")
            
            
            # ------- scorecard -------
    rows = []
    if "pvt" in summary:
        s = summary["pvt"]
        rows += [("PVT: bubble point, Bo, oil density", f"within {s['pb']:.1f}%, {s['bo']:.2f}%, {s['rho']:.2f}%", "OK"),
                 ("PVT: solution GOR", f"within {s['rs']:.1f} scf/STB (absolute)", "OK")]
        if s["mu_cold"]:
            rows.append(("PVT: oil viscosity, hot / cold", f"within {s['mu_hot']:.0f}% / script {s['mu_cold'][0]:+.0f}% to {s['mu_cold'][1]:+.0f}%", "NOTE"))
    for wc in cases:
        s = summary[wc]
        rows += [(f"VLP at 20 rates, water cut {wc}%", f"{s['vlp_min']:+.1f}% to {s['vlp_max']:+.1f}%, mean abs {s['vlp_mabs']:.2f}%", "OK"),
                 (f"Gradient traverse, water cut {wc}%", f"BHP {s['trav']:+.2f}%", "OK")]
    if cases:
        worst = max(summary[c]["ipr_err"] for c in cases)
        extra = [(f"IPR below Pb, water cut {c}%", f"max {summary[c]['ipr_err_lo']:.0f} psi (water treated differently in the Vogel part)", "NOTE")
                 for c in cases if summary[c]["ipr_err_lo"] > 10]
        rows[len(rows) - 2 * len(cases):len(rows) - 2 * len(cases)] = \
            [("IPR pressure above Pb", f"max {worst:.1f} psi", "OK")] + extra
        for wc in cases:
            s = summary[wc]
            if s["pr_op"] and s["op"]:
                rows.append((f"Operating point, water cut {wc}%", f"rate {pct(s['op'][0], s['pr_op'][0]):+.1f}%, BHP {pct(s['op'][1], s['pr_op'][1]):+.2f}%", "OK"))
            elif not s["pr_op"] and not s["op"]:
                rows.append((f"Operating point, water cut {wc}%", "both: no natural flow", "OK"))
            else:
                rows.append((f"Operating point, water cut {wc}%", "PROSPER and Python disagree", "NOTE"))
    if rows:
        fg, a = plt.subplots(figsize=(11, 0.55 * len(rows) + 1.6)); a.axis("off")
        tb = a.table(cellText=[list(r) for r in rows], colLabels=["Check", "Python engine vs PROSPER", "Status"],
                     colWidths=[0.40, 0.46, 0.10], loc="center", cellLoc="left")
        tb.auto_set_font_size(False); tb.set_fontsize(10.5); tb.scale(1, 1.8)
        for (r, c), cell in tb.get_celld().items():
            cell.set_edgecolor("#cccccc")
            if r == 0:
                cell.set_facecolor("#2f3b52"); cell.get_text().set_color("white"); cell.get_text().set_weight("bold")
            elif c == 2:
                cell.set_facecolor("#d8f0d8" if rows[r - 1][2] == "OK" else "#fff0c8"); cell.get_text().set_ha("center")
        a.set_title("Part 1 benchmark scorecard: Python nodal engine vs PROSPER (vertical well, Beggs and Brill)", fontsize=11.5, pad=8)
        fg.savefig("validation_scorecard.png", dpi=160, bbox_inches="tight")
        print("Saved validation_scorecard.png")


if __name__ == "__main__":
    main()     
            
            