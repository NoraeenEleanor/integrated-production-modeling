"""

PART 1's GIF to LinkedIn: Why a well stops flowing.

1.  VLP Surface : bottomhole pressure the tubing needs, vs. liquid rate and WC
                  (Python engine, same WHT at PROSPER at each rate)
2.  IPR Sheet   : bottomhole pressure the reservoir can deliver, vs. luqid rate
3.  Green line  : the intersection = the operating point. It moves to lower rate
                  as WC increased, then disappears: the well stops flowing.
5.  Black dots  : PROSPER VLP results at 0% and 50% WC (from the .PRN exports)


"""

import argparse
import glob
import math
import os
from dataclasses import replace

import imageio.v2 as imageio
import numpy as np
import pyvista as pv
from scipy.interpolate import RectBivariateSpline

import Part1_Nodal_Analysis_Engine as nae

XS, YS, ZS = 1.0, 28.0, 1.4      # display stretch: rate, WC, pressure
# Conversion unit
STB, PSI, ATM, FT = 6.28981, 14.5038, 14.696, 0.3048    # STB/sm3, psi/bar, psi, m/ft
RATE_MAX, WC_MAX, Z_LO, Z_HI = 1500.0, 50.0, 2700.0, 3700.0
CAMERA_DISTANCE = 2.8


#   === PROSPER exports ===

def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None
    

def _lines(path):
    return open(path, encoding="latin-1").read().replace("\r", "").split("\n")


def _input(lines, label):
    for ln in lines:
        t = ln.split("\t")
        if t[0].strip() == label and len(t) > 1:
            return _num(t[1].strip())
    return None


def read_exports(folder):
    system, gradient = {}, {}
    for p in sorted (glob.glob(os.path.join(folder, "*.[Pp][Rr][Nn]"))):
        lines = _lines(p)
        head = "\n".join(lines[:30])
        if "SYSTEM SENSITIVITY ANALYSIS" in head:
            rows = []
            for ln in lines:
                t = [x.strip() for x in ln.split("\t")]
                v = [_num(x) for x in t[:17]]
                if len(t) >= 17 and all(a is not None for a in v):
                    rows.append(v)
            a = np.array(rows)
            system[round(_input(lines, "Water Cut"))] = dict(rate=a[:, 0], vlp=a[:, 4], wht=a[:, 13])
        elif "Gradient Traverse" in head:
            rows = []
            for ln in lines:
                t = [x.strip() for x in ln.split("\t")]
                if len(t) >= 18 and _num(t[1]) is not None and _num(t[3]) is not None and _num(t[1]) <= 8857.0 + 1E-06:
                    rows.append((_num(t[1]), _num(t[4])))
            gradient[round(_input(lines, "Water Cut"))] = np.array(rows)
    return system, gradient



#   === Engine side ===

def compute_grids(f, w0, system, gradient, rates, wcs):
    md, tf = gradient[0][:, 0], gradient[0][:, 1]
    shape = (tf[-1] - tf) / (tf[-1] - tf[0])
    s0, s50 = system[0], system[50]
    
    
    def twh(q, wc):
        a = float(np.interp(q, s0["rate"], s0["wht"]))
        b = float(np.interp(q, s50["rate"], s50["wht"]))
        return (1 - wc / 50.0) * a + (wc / 50.0) * b
    
    
    def vlp(q, wc):
        tz = tf[-1] - (tf[-1] - twh(q, wc)) * shape
        w = replace(w0, wc=wc / 100.0, t_profile=(list(md * FT), [(t - 32) / 1.8 for t in tz]),
                    tvd_m=8858 * FT, n_seg=60)
        return nae.traverse(q / STB, w, f) * PSI - ATM
    
    v = np.array([[vlp(q, wc) for q in rates] for wc in wcs])
    i = np.array([nae.pwf_ipr(q / STB, w0, f.pb) * PSI - ATM for q in rates])
    return v, i


def operating_line(rates, wcs, v, i):
    """Intersection of IPR & VLP for every WC. Stable one: IPR goes from above to below"""
    pts, crit = [], None
    d_first = []
    for k, wc in enumerate(wcs):
        d = i - v[k]
        d_first.append(d[0])
        idx = np.where((d[:-1] > 0) & (d[1:] <= 0))[0]
        if len(idx):
            j = idx[-1]
            t = d[j] / (d[j] - d[j + 1])
            q = rates[j] + t * (rates[j + 1] - rates[j])
            p = i[j] + t * (i[j + 1] - i[j])
            pts.append((q, wc, p))
    d_first = np.array(d_first)
    neg = np.where(d_first <= 0)[0]
    if len(neg) and neg[0] > 0:
        k = neg[0]
        crit = wcs[k - 1] + (wcs[k] - wcs[k - 1]) * d_first[k - 1] / (d_first[k - 1] - d_first[k])
    return np.array(pts), crit



#   === axes drawn by hand ===

class Axes3D:
    """Box, floor grid and tick labels in TRUE units, placed on the two edges. Some versions of PyVista's
       axes show the stretched display coord and ignores the axes_ranges, so this steps is
       needed."""

    def __init__(self, pl, spans, ticks, titles):
        self.pl, self.span, self.ticks, self.titles = pl, np.array(spans, float), ticks, titles
        self.off = 0.045 * max(spans)
        X, Y, Z = self.span
        pl.add_mesh(pv.Box(bounds=(0, X, 0, Y, 0, Z)), style="wireframe", color="#7A7A7A", line_width=1.5)
        for v, _ in ticks[0]:
            pl.add_mesh(pv.Line((v, 0, 0), (v, Y, 0)), color="#CCCCCC", line_width=1)
        for v, _ in ticks[1]:
            pl.add_mesh(pv.Line((0, v, 0), (X, v, 0)), color="#CCCCCC", line_width=1)
            
    
    def _labels(self, pts, texts, name, size=11):
        self.pl.add_point_labels(np.array(pts, float), texts, font_size=size, text_color="black", shape=None,
                                 show_points=False, always_visible=True, name=name)
        
    
    def update(self, cam):
        X, Y, Z = self.span
        o = self.off
        corners = [(0, 0), (X, 0), (X, Y), (0, Y)]
        nx, ny = min(corners, key=lambda c: (c[0] - cam[0]) ** 2 + (c[1] - cam[1]) ** 2)
        dx, dy = (-1 if nx == 0 else 1), (-1 if ny == 0 else 1)
        
        self._labels([[v, ny + dy * o, -o * 0.5] for v, _ in self.ticks[0]], [t for _, t in self.ticks[0]], "ax_x")
        self._labels([[X / 2, ny + dy * 3.2 * o, -o * 0.5]], [self.titles[0]], "ax_xt", 13)
        self._labels([[nx + dx * o, v, -o * 0.5] for v, _ in self.ticks[1]], [t for _, t in self.ticks[1]], "ax_y")
        self._labels([[nx + dx * 3.2 * o, Y / 2, -o * 0.5]], [self.titles[1]], "ax_yt", 13)
        zx = X if nx == 0 else 0
        zdx = 1 if zx == X else -1
        self._labels([[zx + zdx * o, ny + dy * o, v] for v, _ in self.ticks[2]], [t for _, t in self.ticks[2]], "ax_z")
        self._labels([[zx + zdx * 2.5 * o, ny + dy * 2.5 * o, Z + 2.5 * o]], [self.titles[2]], "ax_zt", 13)



#   === Scene ===

def disp(q, wc, p):
    return np.column_stack([np.asarray(q, float) * XS, np.asarray(wc, float) * YS,
                            (np.asarray(p, float) - Z_LO) * ZS])
    
    
def set_view(pl, focal, radius, az_deg, el_deg, axes=None):
    az, el = math.radians(az_deg), math.radians(el_deg)
    pos = focal + radius * np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])
    pl.camera_position = [tuple(pos), tuple(focal), (0, 0, 1)]
    if axes is not None:
        axes.update(pos)
    
    
def make_scene(rates, wcs, v, i, line, crit, system, off_screen, size):
    pl = pv.Plotter(off_screen=off_screen, window_size=size)
    pl.set_background("white")
    
    
    # VLP surface (clipped to the plotted pressure window)
    qq, ww  = np.meshgrid(rates, wcs)
    vc      = np.clip(v, Z_LO, Z_HI)
    surf    = pv.StructuredGrid(qq * XS, ww * YS, (vc - Z_LO) * ZS)
    surf["VLP BHP (psig)"] = v.ravel(order="F")
    pl.add_mesh(surf, scalars="VLP BHP (psig)", cmap="YlOrRd", smooth_shading=True,
                scalar_bar_args=dict(title="VLP BHP (psig)", color="black", vertical=True, fmt="%.0f",
                                     position_x=0.90, position_y=0.40, width=0.04, height=0.42))
    
    
    # IPR surface: the same curve at every WC (water does not change the IPR above Pb)
    ic = np.clip(i, Z_LO, Z_HI)
    ipr_sheet = pv.StructuredGrid(qq * XS, ww * YS, (np.tile(ic, (len(wcs), 1)) - Z_LO) * ZS)
    pl.add_mesh(ipr_sheet, color="#2F7FE0", opacity=0.55, smooth_shading=True)
    
    
    # Operating point line
    if len(line):
        ln = disp(line[:, 0], line[:, 1], line[:, 2])
        ln[:, 2] += 6
        pl.add_mesh(pv.lines_from_points(ln), color="#0A8F3C", line_width=9, render_lines_as_tubes=True)
        lab0 = ln[[0]].copy(); lab0[0, 2] -= 55
        pl.add_point_labels(lab0, ["Operating point"], font_size=15, text_color="#0A6F2F",
                            shape=None, always_visible=True, show_points=False)
        if crit is not None:
            end = ln[[-1]].copy(); end[0, 1] = crit * YS
            pl.add_point_labels(end, [f"Flow stops at ~{crit:.0f}% water cut"], font_size=15, text_color="black",
                                shape=None, always_visible=True, show_points=False)
            
            
    # PROSPER results (black dots)
    for wc, S in system.items():
        m = S["rate"] <= RATE_MAX
        for p3 in disp(S["rate"][m], np.full(m.sum(), wc), S["vlp"][m]):
            pl.add_mesh(pv.Sphere(radius=11, center=p3 + [0, 0, 4]), color="black")
            
    ticks = [[(v * XS, f"{v:.0f}") for v in np.linspace(0, RATE_MAX, 5)],
             [(v * YS, f"{v:.0f}") for v in np.linspace(0, WC_MAX, 6)],
             [((v - Z_LO) * ZS, f"{v:.0f}") for v in np.linspace(Z_LO, Z_HI, 5)]]
    axes = Axes3D(pl, (RATE_MAX * XS, WC_MAX * YS, (Z_HI - Z_LO) * ZS), ticks,
                  ["Liquid rate (STB/d)", "Water cut (%)", "BHP (psig)"])
    pl.add_text("Why a well stops flowing", position="upper_edge", font_size=15, color="black")
    pl.add_text("VLP surface (orange) rises above the IPR sheet (blue as water cut increases.\n"
                "Black dots: PROSPER, 0% and 50% water cut", position=(18, 655), font_size=9,
                color="black")
    focal = np.array([0.5 * RATE_MAX * XS, 0.5 * WC_MAX * YS, 0.5 * (Z_HI - Z_LO) * ZS])
    radius = CAMERA_DISTANCE * math.sqrt((RATE_MAX * XS) ** 2 + (WC_MAX * YS) ** 2)
    return pl, focal, radius, axes



def main ():
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("--props", default="CASE_PROPS.INC")
    ap.add_argument("--prosper-dir", default="PROSPER_exports")
    ap.add_argument("--out", default="nodal_surface")
    ap.add_argument("--frames", type=float, default=90)
    ap.add_argument("--fps", type=float, default=12.0)
    ap.add_argument("--interactive", action="store_true")
    args, unknown = ap.parse_known_args()
    
    
    system, gradient = read_exports(args.prosper_dir)
    if 0 not in system or 50 not in system or 0 not in gradient:
        folder  = os.path.abspath(args.prosper_dir)
        found   = [os.path.basename(p) for p in glob.glob(os.path.join(folder, "*.[Pp][Rr][Nn]"))]
        print(f"Looked in: {folder}")
        print(f"    .PRN files found: {found if found else 'NONE (is --prosper-dir pointing at the right folder?)'}")
        print(f"    recognised: System at water cut {sorted(system)}, Gradient Traverse at water cut {sorted(gradient)}")
        raise SystemExit("Need the System export at water cut 0 and 50 and the Gradient Traverse at water cut 0.")
    
    w0          = nae.Well()
    f           = nae.Fluid(args.props, w0.t_res_c)
    rates_c     = np.linspace(6.0, RATE_MAX, 41)
    wcs_c       = np.linspace(0.0, WC_MAX, 21)
    
    print("Computing the VLP surface (about half a minute)...")
    
    v_c, _      = compute_grids(f, w0, system, gradient, rates_c, wcs_c)
    rates       = np.linspace(6.0, RATE_MAX, 150)
    wcs         = np.linspace(0.0, WC_MAX, 101)
    v           = RectBivariateSpline(wcs_c, rates_c, v_c, kx=3, ky=3)(wcs, rates)
    i           = np.array([nae.pwf_ipr(q / STB, w0, f.pb) * PSI - ATM for q in rates])
    line, crit  = operating_line(rates, wcs, v, i)
    
    print(f"Operating point at 0% water cut: {line[0, 0]:.0f} STB/d, {line[0, 2]:.0f} psig" if len(line) else "no operating point")
    print(f"Natural flow stops at about {crit:.0f}% water cut" if crit is not None else "")
    
    size = (1000, 750)
    if args.interactive:
        pl, focal, radius, axes = make_scene(rates, wcs, v, i, line, crit, system, False, size)
        set_view(pl, focal, radius, -55, 25, axes)
        pl.show()
        return
    pl, focal, radius, axes = make_scene(rates, wcs, v, i, line, crit, system, True, size)
    set_view(pl, focal, radius, -55, 25, axes)
    pl.show(auto_close=False)
    pl.screenshot(f"{args.out}_tilted.png")
    frames = []
    for k in range(args.frames):
        set_view(pl, focal, radius, -55 + 360.0 * k / args.frames, 25, axes)
        pl.render()
        frames.append(pl.screenshot(return_img=True))
    imageio.mimsave(f"{args.out}_orbit.gif", frames, duration=1000 / args.fps, loop=0)
    pl.close()
    print(f"Saved {args.out}_tilted.png and {args.out}_orbit.gif")
    
    
    
if __name__ == "__main__":
    main()
    