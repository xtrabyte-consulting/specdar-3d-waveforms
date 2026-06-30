# ----------------------------------------------------------------------
# ANALYTICAL near-field model of a z-oriented bowtie Robert's dipole radiating
# toward a planar conducting reflector. Standing wave = direct field +
# image field. This is an illustrative analytical model, NOT a full-wave
# solve. Physics-correct fringe spacing (lambda/2) and dipole sin(theta)
# pattern are used here for a visual approximation.
# ----------------------------------------------------------------------

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


c = 2.99792458e8
f = 1.0e9                 # 1 GHz carrier
lam = c / f               # 0.30 m
k = 2 * np.pi / lam
D = 0.45                  # reflector standoff (= 1.5 lambda)
x_img = 2 * D             # image dipole location (parallel dipole -> reversed sign)

def efield(x, y, z):
    """Complex transverse E-field magnitude in x-z type geometry.
    Direct dipole at origin (axis z), image at x=2D with reversed current."""
    eps = 1e-3
    # source 1 at origin
    r1 = np.sqrt(x**2 + y**2 + z**2) + eps
    sin1 = np.sqrt(x**2 + y**2) / r1            # angle from z-axis
    E1 = sin1 / r1 * np.exp(-1j * k * r1)
    # image source at (x_img, 0, 0), reversed sign (conducting plane)
    r2 = np.sqrt((x - x_img)**2 + y**2 + z**2) + eps
    sin2 = np.sqrt((x - x_img)**2 + y**2) / r2
    E2 = sin2 / r2 * np.exp(-1j * k * r2)
    return E1 - E2

# ---------- color / style ----------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.edgecolor": "#3a4a5a",
    "text.color": "#e8edf2",
    "axes.labelcolor": "#cfd8e3",
    "xtick.color": "#9fb0c0",
    "ytick.color": "#9fb0c0",
})
BG = "#0e141b"
cmap = "turbo"

fig = plt.figure(figsize=(15, 8.2), facecolor=BG)
gs = gridspec.GridSpec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 0.85],
                       wspace=0.18, hspace=0.30,
                       left=0.02, right=0.97, top=0.90, bottom=0.08)

# ======================================================================
# PANEL A : 3D cut-plane render
# ======================================================================
axA = fig.add_subplot(gs[:, 0], projection="3d")
axA.set_facecolor(BG)

# grid for x-z plane (y=0) and x-y plane (z=0)
nx, nz = 220, 180
xv = np.linspace(0.0, 0.60, nx)
zv = np.linspace(-0.27, 0.27, nz)
Xxz, Zxz = np.meshgrid(xv, zv)
Exz = np.abs(efield(Xxz, 0.0, Zxz))

yv = np.linspace(-0.27, 0.27, nz)
Xxy, Yxy = np.meshgrid(xv, yv)
Exy = np.abs(efield(Xxy, Yxy, 0.0))

# normalize together with robust limits
allvals = np.concatenate([Exz.ravel(), Exy.ravel()])
vmin, vmax = np.percentile(allvals, 3), np.percentile(allvals, 97)
norm = plt.Normalize(vmin, vmax)

cmA = matplotlib.colormaps[cmap]
# x-z plane at y=0
facxz = cmA(norm(Exz))
axA.plot_surface(Xxz, np.zeros_like(Xxz), Zxz, facecolors=facxz,
                 rstride=1, cstride=1, shade=False, antialiased=False,
                 linewidth=0, zorder=1)
# x-y plane at z=0
facxy = cmA(norm(Exy))
axA.plot_surface(Xxy, Yxy, np.zeros_like(Xxy), facecolors=facxy,
                 rstride=1, cstride=1, shade=False, antialiased=False,
                 linewidth=0, zorder=1)

# bowtie geometry (two triangles along z, apex at origin, in y=0 plane)
flare = 0.06
L = 0.075
tri_top = [(0, 0, 0), (flare, 0, L), (-flare, 0, L)]
tri_bot = [(0, 0, 0), (flare, 0, -L), (-flare, 0, -L)]
for tri in (tri_top, tri_bot):
    poly = Poly3DCollection([tri], facecolor="#f4d35e", edgecolor="#ffffff",
                            linewidths=1.2, alpha=0.97)
    poly.set_zorder(5)
    axA.add_collection3d(poly)

# reflector plane at x=D
ry = np.array([-0.27, 0.27]); rz = np.array([-0.27, 0.27])
RY, RZ = np.meshgrid(ry, rz)
RX = np.full_like(RY, D)
axA.plot_surface(RX, RY, RZ, color="#7fa8c9", alpha=0.18, zorder=2)
axA.text(D, 0, 0.30, "conducting\nreflector", color="#bcd4ea", fontsize=9, ha="center")

axA.set_xlabel("x  (m)  — propagation", labelpad=8)
axA.set_ylabel("y  (m)", labelpad=6)
axA.set_zlabel("z  (m)  — dipole axis", labelpad=4)
axA.set_xlim(0, 0.60); axA.set_ylim(-0.27, 0.27); axA.set_zlim(-0.27, 0.27)
axA.view_init(elev=22, azim=-58)
axA.xaxis.pane.set_facecolor(BG); axA.yaxis.pane.set_facecolor(BG); axA.zaxis.pane.set_facecolor(BG)
axA.xaxis.pane.set_alpha(1); axA.yaxis.pane.set_alpha(1); axA.zaxis.pane.set_alpha(1)
axA.grid(False)
for axis in (axA.xaxis, axA.yaxis, axA.zaxis):
    axis.line.set_color("#33424f")
axA.set_title("3D field distribution — orthogonal cut planes",
              color="#e8edf2", fontsize=12.5, pad=2)

# ======================================================================
# PANEL B : 2D x-z heat map with fringe annotation
# ======================================================================
axB = fig.add_subplot(gs[0, 1])
axB.set_facecolor(BG)
im = axB.pcolormesh(xv, zv, Exz, cmap=cmap, norm=norm, shading="auto")
axB.axvline(D, color="#bcd4ea", lw=2, ls=(0, (4, 3)))
axB.text(D + 0.005, 0.21, "reflector", color="#bcd4ea", fontsize=9, rotation=90, va="top")
# bowtie marker
axB.add_patch(Polygon([(0, 0), (flare, L), (-flare, L)], closed=True, color="#f4d35e"))
axB.add_patch(Polygon([(0, 0), (flare, -L), (-flare, -L)], closed=True, color="#f4d35e"))
# lambda/2 node spacing annotation
n0 = D - lam/2 * 0.5
axB.annotate("", xy=(D, -0.20), xytext=(D - lam/2, -0.20),
             arrowprops=dict(arrowstyle="<->", color="#ffffff", lw=1.4))
axB.text((D - lam/4), -0.235, r"$\lambda/2$ = 15 cm", color="#ffffff",
         fontsize=9.5, ha="center")
axB.set_xlim(0, 0.60); axB.set_ylim(-0.27, 0.27)
axB.set_xlabel("x  (m)"); axB.set_ylabel("z  (m)")
axB.set_title("Standing-wave |E| (x–z plane)", color="#e8edf2", fontsize=12, pad=6)
cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.02)
cb.set_label("|E|  (norm.)", color="#cfd8e3", fontsize=9)
cb.ax.yaxis.set_tick_params(color="#9fb0c0")
plt.setp(plt.getp(cb.ax, "yticklabels"), color="#9fb0c0")

# ======================================================================
# PANEL C : pulse-modulated 1 GHz waveform  (trapezoidal PRF train)
# ======================================================================
axC = fig.add_subplot(gs[1, 1])
axC.set_facecolor(BG)

# pulse timing (ns)
t_rise  = 2.0      # leading-edge transition
t_fall  = 2.0      # trailing-edge transition
t_width = 20.0     # full-amplitude flat top
t_delay = 80.0     # leading-edge start within each period
T_per   = 200.0    # period -> 5 MHz PRF
duty    = t_width / T_per * 100.0

def trap_env(t_ns):
    """Periodic trapezoidal envelope (0..1)."""
    tp = np.mod(t_ns - t_delay, T_per)     # time since leading edge
    e = np.zeros_like(tp)
    m1 = (tp >= 0) & (tp < t_rise)
    e[m1] = tp[m1] / t_rise
    m2 = (tp >= t_rise) & (tp < t_rise + t_width)
    e[m2] = 1.0
    m3 = (tp >= t_rise + t_width) & (tp < t_rise + t_width + t_fall)
    e[m3] = 1.0 - (tp[m3] - (t_rise + t_width)) / t_fall
    return e

t_ns    = np.linspace(0, 400, 40000)       # two full periods
env     = trap_env(t_ns)
sig     = env * np.cos(2 * np.pi * f * (t_ns * 1e-9))

axC.plot(t_ns, sig, color="#4cc9f0", lw=0.35, alpha=0.85)
axC.plot(t_ns, env, color="#f4d35e", lw=1.6, label="envelope")
axC.plot(t_ns, -env, color="#f4d35e", lw=1.6)
axC.set_xlim(0, 400); axC.set_ylim(-1.25, 1.5)
axC.set_xlabel("time  (ns)"); axC.set_ylabel("amplitude")
axC.set_title("Pulse-modulated 1 GHz excitation  ·  PRF 5 MHz, 10% duty",
              color="#e8edf2", fontsize=11.5, pad=6)
axC.grid(True, color="#22303c", lw=0.6)

# period (leading edge -> leading edge)
axC.annotate("", xy=(t_delay + T_per, 1.30), xytext=(t_delay, 1.30),
             arrowprops=dict(arrowstyle="<->", color="#ffffff", lw=1.3))
axC.text(t_delay + T_per/2, 1.355, "period = 200 ns", color="#ffffff",
         fontsize=8.5, ha="center", va="bottom")
# delay
axC.annotate("", xy=(t_delay, -1.08), xytext=(0, -1.08),
             arrowprops=dict(arrowstyle="<->", color="#9fb0c0", lw=1.1))
axC.text(t_delay/2, -1.20, "delay 80 ns", color="#cfd8e3",
         fontsize=8, ha="center", va="top")
# pulse width on first pulse
axC.annotate("", xy=(t_delay + t_rise + t_width, 1.10), xytext=(t_delay + t_rise, 1.10),
             arrowprops=dict(arrowstyle="<->", color="#f4d35e", lw=1.1))
axC.text(t_delay + t_rise + t_width/2, 1.14, "PW 20 ns", color="#f4d35e",
         fontsize=8, ha="center", va="bottom")
axC.legend(loc="lower right", facecolor="#16202b", edgecolor="#33424f",
           labelcolor="#cfd8e3", fontsize=8.5)

# inset: leading edge detail (2 ns rise + resolved carrier), in the
# quiet inter-pulse region so it overlaps no signal
axins = axC.inset_axes([0.40, 0.55, 0.27, 0.33])
axins.set_facecolor("#16202b")
tz = np.linspace(78, 90, 3000)
envz = trap_env(tz)
axins.plot(tz, envz * np.cos(2*np.pi*f*(tz*1e-9)), color="#4cc9f0", lw=0.6)
axins.plot(tz, envz, color="#f4d35e", lw=1.3)
axins.axvspan(80, 82, color="#f4d35e", alpha=0.12)
axins.set_xlim(78, 90); axins.set_ylim(-1.1, 1.1)
axins.tick_params(colors="#9fb0c0", labelsize=6.5)
axins.set_title("leading edge — 2 ns rise", color="#cfd8e3", fontsize=7.5, pad=2)
for s in axins.spines.values():
    s.set_color("#33424f")

# ---------- figure title & footer ----------
fig.suptitle("Pulse-Modulated 1 GHz Bowtie Dipole — Standing-Wave Field Distribution",
             color="#ffffff", fontsize=16, fontweight="bold", y=0.975)
fig.text(0.02, 0.015,
         "Analytical near-field model (direct + image dipole) · f = 1 GHz, "
         r"$\lambda$ = 30 cm, reflector standoff = 1.5$\lambda$ · "
         "pulse: 2 ns edges / 20 ns PW / 80 ns delay / 200 ns period · "
         "Illustrative — not a full-wave solve",
         color="#7d8c9a", fontsize=8.0)
fig.text(0.97, 0.015, "Vigilon Cyber", color="#5f6f7d", fontsize=8.5, ha="right")

fig.savefig("/home/claude/bowtie_standing_wave.png", dpi=300,
            facecolor=BG, bbox_inches="tight")
print("saved")