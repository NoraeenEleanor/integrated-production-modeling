
"""

PART 1: Python Nodal Analysis Engine (IPR + VLP)

Single well, vertical, steady state nodal analysis for the 26x31x7 sector
model. Fluid properties are directly from the simulation deck (CASE_PROPS.INC
or the exported dataset .GRDECL: PVTW (Water), RSCONSTT, PVDO (Oil), DENSITY),
so Python, PROSPER and the simulator all start from one fluid.

1.  PVT     : deck tables (PVDO) for ABOVE Pb used Standing meanwhile Beggs-Robinson
              for BELOW Pb.
2.  IPR     : IPR Composites - PI Linear for ABOVE Pb; Vogel for BELOW Pb
3.  VLP     : Beggs & Brill (1973, original, vertical upflow) marched from
              the wellhead down to the bottomhole.
4.  Node    : intersection of IPR and VLP at the bottomhole (operating point).
5.  Output  : sensitivity plots + CSV files for the PROSPER comparison.

Each values in the ASSUMPTIONS block below is an assumption, except t_res_c (Petrel
fluid model) and pr_bar (deck EQUIL datum pressure).

"""

import argparse
import csv
import math
import re
from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G               = 9.80665       # m/s2
PSC_BAR         = 1.01325       # std pressure
TSC_K           = 288.706       # 15.56 deg C std temperature
RHO_AIR_SC      = 1.2255        # kg/m3 at std cond
RHO_WATER_SC    = 999.016       # kg/m3 fresh water at 60 deg F (for oil API)
BAR_TO_PSI      = 14.5038       
SCF_PER_SM3     = 5.614584      # (scf/stb) per (sm3/sm3)




# ASSUMPTIONS (NOT in the deck - you may changed based on your dataset)

@dataclass
class Well:
    tvd_m:          float       = 2700.0        # (8858 ft) perforation TVD. TODO: mean TVD of producer cells from the grid
    id_m:           float       = 0.1005        # 4.1/2" tubing, ID = 3.958 in
    rough_m:        float       = 1.5E-5        # commercial steel
    whp_bar:        float       = 5.0           # wellhead (flowing tubing head) pressure | same with PROSPER
    t_wh_c:         float       = 50.0          # flowing wellhead temperature
    t_res_c:        float       = 76.85         # (170 deg F) reservoir Temperature (value from PETREL)
    pr_bar:         float       = 240.0         # reservoir Pressure at TVD above (EQUIL: 240 bar @ 2700m) (value from PETREL)
    pi:             float       = 5.0           # (2.17 stb/day/psi) liquid PI (sm3/d/bar). TODO: calibrate from the simulator
    wc:             float       = 0.0           # water cut, fraction of liquid
    gor:        Optional[float] = None          # producing GOR (sm3/sm3); None -> solution GOR (NO free gas ABOVE Pb)
    sigma_o:        float       = 0.030         # GAS-OIL surface tension (N/m)
    sigma_w:        float       = 0.060         # GAS-WATER surface tension (N/m)
    n_seg:            int       = 54            # pressure traverse segments (50m each at 2700m)
    t_profile:  Optional[tuple] = None          # (depth_m list, temp_C list): overrides the linear temp profile
    
    
    
#   PVT from the deck

class Fluid:
    def __init__(self, path, t_res_c):
        txt = re.sub(r"--[^\n]*", "", open(path).read())
        
        def nums(kw):
            m = re.search(rf"(?m)^\s*{kw}\b(.*?)/", txt, flags=re.S)
            if not m:
                raise ValueError(f"{kw} not found in {path}")
            return [float(x) for x in m.group(1).split()]
        
        self.pref, self.bw_ref, self.cw, self.muw   = nums("PVTW")[:4]
        self.rsb, self.pb                           = nums("RSCONSTT")[:2]
        tab                                         = np.array(nums("PVDO")).reshape(-1, 3)
        self.p_tab, self.bo_tab, self.mu_tab        = tab[:, 0], tab[:, 1], tab[:, 2]
        self.rho_o_sc, self.rho_w_sc, self.rho_g_sc = nums("DENSITY")[:3]
        
        self.sg_o       =   self.rho_o_sc / RHO_WATER_SC
        self.api        =   141.5 / self.sg_o - 131.5
        self.sg_g       =   self.rho_g_sc / RHO_AIR_SC
        self.t_res_c    =   t_res_c
        self.bo_pb      =   float(np.interp(self.pb, self.p_tab, self.bo_tab))
        mu_pb           =   float(np.interp(self.pb, self.p_tab, self.mu_tab))
        self.k_mu       =   mu_pb / self._mu_br(self.rsb, t_res_c)
        
        
    #   ===========
    #   --- OIL ---
    #   ===========
    
    def pb_at(self, t_c):
        k = 10 ** (0.00091 * 1.8 * (t_c - self.t_res_c))
        return ((self.pb * BAR_TO_PSI / 18.2 + 1.4) * k - 1.4) * 18.2 / BAR_TO_PSI
    
    
    def rs(self, p, t_c=None):
        pb = self.pb if t_c is None else self.pb_at(t_c)
        if p >= pb:
            return self.rsb
        shape = (p * BAR_TO_PSI / 18.2 + 1.4) / (pb * BAR_TO_PSI / 18.2 + 1.4)
        return self.rsb * shape ** 1.2048
    
    
    def _bo_std(self, rs_sm3, t_c):
        t_f     = t_c * 1.8 + 32.0
        rs      = rs_sm3 * SCF_PER_SM3
        return 0.9759 + 0.00012 * (rs * math.sqrt(self.sg_g / self.sg_o) + 1.25 * t_f) ** 1.2
    
    
    def bo(self, p, t_c):
        thermal = self._bo_std(self.rs(p, t_c), t_c) / self._bo_std(self.rsb, self.t_res_c)
        return float(np.interp(max(p, self.pb), self.p_tab, self.bo_tab)) * thermal
    
    
    def _mu_br(self, rs_sm3, t_c):
        t_f         = t_c * 1.8 + 32.0
        mu_dead     = 10 ** (10 ** (3.0324 - 0.02023 * self.api) * t_f ** -1.163) - 1.0
        rs          = rs_sm3 * SCF_PER_SM3
        a           = 10.715 * (rs + 100.0) ** -0.515
        b           = 5.44 * (rs + 150.0) ** -0.338
        return a * mu_dead ** b
    
    
    def mu_o(self, p, t_c):
        mu_br = self.k_mu * self._mu_br(self.rs(p, t_c), t_c)
        if p < self.pb:
            return mu_br
        mu_deck = float(np.interp(p, self.p_tab, self.mu_tab))
        return mu_br * mu_deck / float(np.interp(self.pb, self.p_tab, self.mu_tab))
    
    
    
    #   ===========
    #   --- GAS ---
    #   ===========
    
    def z(self, p, t_c):            # Papay, Sutton pseudo criticals
        sg      = self.sg_g
        ppc     = 756.8 - 131.0 * sg - 3.6 * sg ** 2
        tpc     = 169.2 + 349.5 * sg - 74.0 * sg ** 2
        ppr     = p * BAR_TO_PSI / ppc
        tpr     = (t_c + 273.15) * 1.8 / tpc
        return 1 - 3.52 * ppr / 10 ** (0.9813 * tpr) + 0.274 * ppr ** 2 / 10 ** (0.8157 * tpr)
    
    
    def bg(self, p, t_c):           # rm3/sm3
        return self.z(p, t_c) * (t_c + 273.15) / TSC_K * PSC_BAR / p
    
    
    def mu_g(self, p, t_c):         # Lee-Gonzalez-Eakin, cP
        t_r     = (t_c + 273.15) * 1.8
        m       = 28.97 * self.sg_g
        rho     = self.rho_g_sc / self.bg(p, t_c) / 1000.0
        k       = (9.4 + 0.02 * m) * t_r ** 1.5 / (209.0 + 19.0 * m + t_r)
        x       = 3.5 + 986.0 / t_r + 0.01 * m
        y       = 2.4 - 0.2 * x
        return 1E-4 * k * math.exp(x * rho ** y)
    
    
    
    #   =============
    #   --- WATER ---
    #   =============
    
    def bw(self, p):
        return self.bw_ref * math.exp(-self.cw * (p - self.pref))
        
        
        
        
#   =======
#   Beggs & Brill (1973), vertical upflow
#   =======

_HL0    = {"seg": (0.98, 0.4846, 0.0868), "int": (0.845, 0.5351, 0.0173), "dist": (1.065, 0.5824, 0.0609)}
_CORR   = {"seg": (0.011, -3.768, 3.539, -1.614), "int": (2.96, 0.305, -0.4473, 0.0978)}


def _hl_regime(regime, lam, nfr, nlv):
    a, b, c     = _HL0[regime]
    hl0         = max(a * lam ** b / nfr ** c, lam)
    cc          = 0.0
    if regime in _CORR:
        d, e, f, g      = _CORR[regime]
        cc              = max((1 - lam) * math.log(d * lam ** e * nlv ** f * nfr ** g), 0.0)
    th          = 1.8 * math.pi / 2
    psi         = 1 + cc * (math.sin(th) - 0.333 * math.sin(th) ** 3)
    return min(hl0 * psi, 1.0)
    


def holdup(lam, nfr, nlv):
    if lam >= 0.99999:
        return 1.0
    l1 = 316 * lam ** 0.302
    l2 = 0.0009252 * lam ** -2.4684
    l3 = 0.10 * lam ** -1.4516
    l4 = 0.5 * lam ** -6.738
    
    if (lam < 0.01 and nfr < l1) or (lam >= 0.01 and nfr < l2):
        return _hl_regime("seg", lam, nfr, nlv)
    
    if lam >= 0.01 and l2 <= nfr <= l3:
        a = (l3 - nfr) / (l3 - l2)
        return a * _hl_regime("seg", lam, nfr, nlv) + (1 - a) * _hl_regime("int", lam, nfr, nlv)
    
    if (0.01 <= lam < 0.4 and l3 < nfr <= l1) or (lam >= 0.4 and l3 < nfr <= l4):
        return _hl_regime("int", lam, nfr, nlv)
    return _hl_regime("dist", lam, nfr, nlv)



def darcy_friction(re_, rel_rough):
    if re_ < 2300:
        return 64.0 / re_
    return 0.25 / math.log10(rel_rough / 3.7 + 5.74 / re_ ** 0.9) ** 2



def gradient(p, t_c, qo, qw, gor, w, f, ret_holdup=False, parts=False):
    """Pressure gradient dP/dz (bar per metre of depth) at pressure p [bar], T [C]."""
    area                = math.pi * w.id_m ** 2 / 4
    rs, bo, bw, bg      = f.rs(p, t_c), f.bo(p, t_c), f.bw(p), f.bg(p, t_c)
    q_gas_sc            = qo * max(gor - rs, 0.0)
    vsl                 = (qo * bo + qw * bw) / 86400.0 / area
    vsg                 = q_gas_sc * bg /86400.0 / area
    vm                  = vsl + vsg
    lam                 = vsl / vm
    
    rho_o       = (f.rho_o_sc + rs * f.rho_g_sc) / bo
    rho_w       = f.rho_w_sc / bw
    fo          = qo * bo / (qo * bo + qw * bw)
    rho_l       = fo * rho_o + (1 - fo) * rho_w
    mu_l        = (fo * f.mu_o(p, t_c) + (1 - fo) * f.muw) * 1E-3
    rho_g       = f.rho_g_sc / bg
    mu_g        = f.mu_g(p, t_c) * 1E-3
    sigma       = fo * w.sigma_o + (1 - fo) * w.sigma_w
    
    nfr         = vm ** 2 / (G * w.id_m)
    nlv         = vsl * (rho_l / (G * sigma)) ** 0.25
    hl          = holdup(lam, nfr, nlv)
    
    rho_s       = rho_l * hl + rho_g * (1 - hl)
    rho_ns      = rho_l * lam + rho_g * (1 - lam)
    mu_ns       = mu_l * lam + mu_g * (1 - lam)
    re_         = rho_ns * vm * w.id_m / mu_ns
    fn          = darcy_friction(re_, w.rough_m / w.id_m)
    
    
    y = lam / hl ** 2
    if 1.0 < y < 1.2:
        s = math.log(2.2 * y - 1.2)
    else:
        ly = math.log(y)
        s = ly / (-0.0523 + 3.182 * ly - 0.8725 * ly ** 2 + 0.01853 * ly ** 4)
    ftp = fn * math.exp(s)
    
    ek = rho_s * vm * vsg / (p * 1E05)
    dpdl = (rho_s * G + ftp * rho_ns * vm ** 2 / (2 * w.id_m)) / (1 - ek)
    if parts:
        return (dpdl / 1E05, hl, rho_s * G / (1 - ek) / 1E05,
                ftp * rho_ns * vm ** 2 / (2 * w.id_m) / (1 - ek) / 1E05)
    return (dpdl / 1E05, hl) if ret_holdup else dpdl / 1E05



def temp_at(z, w):
    """Flowing temperature [C] st depth z [m]."""
    if w.t_profile is not None:
        return float(np.interp(z, w.t_profile[0], w.t_profile[1]))
    return w.t_wh_c + (w.t_res_c - w.t_wh_c) * z / w.tvd_m



def traverse(q_liq, w, f, profile=False, detail=False):
    """March from wellhead to bottomhole. Returns BHP [bar] (and the profile)."""
    qo, qw      = q_liq * (1 - w.wc), q_liq * w.wc
    gor         = f.rsb if w.gor is None else w.gor
    dz          = w.tvd_m / w.n_seg
    p, z        = w.whp_bar, 0.0
    prof        = [(0.0, p)]
    det         = []
    cum_g       = cum_f = 0.0
    
    for _ in range(w.n_seg):
        t_mid = temp_at(z + dz / 2, w)
        p_end = p + 0.08 * dz
        for _ in range(4):
            g, hl, gg, gf = gradient(0.5 * (p + p_end), t_mid, qo, qw, gor, w, f, parts=True)
            p_end = p + g * dz
        cum_g, cum_f = cum_g + gg * dz, cum_f + gf * dz
        p, z = p_end, z + dz
        prof.append((z, p))
        det.append((z, p, temp_at(z, w), hl, cum_f, cum_g))
    if detail:
        return det
    return (p, prof) if profile else p




#   =======
#   IPR and node solution
#   =======

def pwf_ipr(q, w, pb):
    qb = w.pi * (w.pr_bar - pb)
    if q <= qb:
        return w.pr_bar - q / w.pi
    r = (q - qb) / (w.pi * pb / 1.8)
    if r >= 1:
        return 0.0
    x = (-0.2 + math.sqrt(0.04 + 3.2 * (1 - r))) / 1.6
    return x * pb



def aof(w, pb):
    return w.pi * (w.pr_bar - pb) + w.pi * pb / 1.8



def operating_point(w, f):
    """Largest rate intersection of IPR and VLP; None if the well cannot flow unaided."""
    qs = np.linspace(2.0, aof(w, f.pb) * 0.999, 80)
    diff = np.array([pwf_ipr(q, w, f.pb) - traverse(q, w, f) for q in qs])
    if diff[0] < 0:
        return None
    idx = np.where(np.diff(np.sign(diff)) != 0)[0]
    if len(idx) == 0:
        return (qs[-1], pwf_ipr(qs[-1], w, f.pb))
    i = idx[-1]
    q = brentq(lambda x: pwf_ipr(x, w, f.pb) - traverse(x, w, f), qs[i], qs[i + 1])
    return q, pwf_ipr(q, w, f.pb)




#   =======
#   MAIN
#   =======

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--props", default="CASE_PROPS.INC")
    args, unknown = ap.parse_known_args()      # cannot change this, gives error
    
    w0 = Well()
    f = Fluid(args.props, w0.t_res_c)
    
    print("== Fluid read from deck ==")
    print(f"    Rs = {f.rsb:.3f} sm3/sm3  Pb = {f.pb:.1f} bar  API = {f.api:.1f}  gas SG = {f.sg_g:.3f}")
    print(f"    Bo(Pb) = {f.bo(f.pb, w0.t_res_c):.4f}  mu_o(Pb) = {f.mu_o(f.pb, w0.t_res_c):.3f} cP   (deck: {f.bo_pb:.4f}, {np.interp(f.pb, f.p_tab, f.mu_tab):.3f})")
    print("     P(bar)    Rs    Bo    mu_o    Bg(rm3/sm3) (at mid well T)")
    t_mid = 0.5 * (w0.t_wh_c + w0.t_res_c)
    for p in (10, 20, 40, 80, 150, 240):
        print(f"    {p:5d}  {f.rs(p, t_mid):6.2f}  {f.bo(p, t_mid):.4f}  {f.mu_o(p, t_mid):5.2f}   {f.bg(p, t_mid):.4f}")
        
    print("\n== Base Case ==")
    op = operating_point(w0, f)
    if op is None:
        print(" VLP BHP exceeds the reservoir pressure even at near zero rate: no natural flow.")
        bhp0 = traverse(2.0, w0, f)
        print(f"    VLP BHP at 2 sm3/d = {bhp0:.1f} bar  vs  Pr = {w0.pr_bar:.1f} bar")
    else:
        print(f"    Operating point: q = {op[0]:.1f} sm3/d, BHP = {op[1]:.1f} bar")
        
    qs = np.linspace(2.0, aof(w0, f.pb), 60)
    ipr = [pwf_ipr(q, w0, f.pb) for q in qs]
    vlp = [traverse(q, w0, f) for q in qs]
    with open("vlp_curve_base.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["q_liq_sm3d", "bhp_vlp_bar", "pwf_ipr_bar"])
        wr.writerows(zip(qs, vlp, ipr))
        
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    # water cur sensitivity is run at WHP=5bar, because at WHP=20bar the base well does not flow at all.
    for name, vals, rng, base in (("WHP (bar)", "whp_bar", (5, 10, 20), w0),
                                  ("Water cut (WHP 5 bar)", "wc", (0.0, 0.5, 0.9), replace(w0, whp_bar=5.0))):
        a = ax[0] if vals == "whp_bar" else ax[1]
        a.plot(qs, ipr, "k", lw=2, label="IPR")
        print(f"\n== Sensitivity: {name} ==")
        for v in rng:
            w = replace(base, **{vals: v})
            curve = [traverse(q, w, f) for q in qs]
            a.plot(qs, curve, label=f"VLP {name} = {v}")
            o = operating_point(w, f)
            print(f"    {name} = {v}: " + ("no natural flow" if o is None else f"q = {o[0]:.1f} sm3/d, BHP = {o[1]:.1f} bar"))
        a.set_xlabel("Liquid rate (sm3/d)")
        a.set_ylabel("Bottomhole pressure (bar)")
        a.set_title(f"Nodal analysis - sensitivity to {name}")
        a.grid(alpha=0.3)
        a.legend()
    fig.tight_layout()
    fig.savefig("nodal_sensitivity.png", dpi=150)
    print("\nSaved: vlp_curve_base.csv, nodal_sensitivity.png")
    
    

if __name__ == "__main__":
    main()