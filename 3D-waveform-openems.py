#!/usr/bin/env python3
r"""
3D-waveform-openems.py
======================
Full-wave 3-D FDTD solve (openEMS) of a *pulse-modulated 1 GHz* bowtie dipole
radiating into free space, with time-domain E/H field dumps written for ParaView.

This is the full-physics companion to ``3D-waveform-analytical.py`` (which is only
an illustrative image-dipole near-field model). Here Maxwell's equations are solved
directly on a Yee grid by openEMS; the standing/travelling-wave structure, the
finite bandwidth of the pulse, and the true radiation pattern all emerge from the
solve rather than being assumed.

Excitation (per the requested waveform)
---------------------------------------
    * carrier        : 1.0 GHz sinusoid
    * pulse width    : 20 ns full-amplitude flat top
    * rise / fall    : 2 ns linear (trapezoidal envelope)
    * delay          : 120 ns before the leading edge
A single pulse is launched (PRF/period is irrelevant for a transient solve of one
gate). The envelope is a trapezoid; multiply by the carrier to get the source.

    env(t) = 0                          t < td
             (t-td)/tr                  td      <= t < td+tr
             1                          td+tr   <= t < td+tr+tw
             (t4-t)/tf                  td+tr+tw<= t < t4
             0                          t >= t4         (t4 = td+tr+tw+tf)

    s(t)   = env(t) * sin( 2*pi*f0*(t - td) )

Geometry
--------
z-oriented planar bowtie (two triangles in the y=0 plane, narrow at the central
feed gap, flaring outward) fed by a 50 ohm lumped port across the gap. PML on all
six faces => open free-space radiation.

Outputs (in OUT_DIR)
--------------------
    Ez_xz_*.vtr   E-field time series on the y=0 plane (the classic E-plane cut)
    Ez_xy_*.vtr   E-field time series on the z=0 plane (H-plane cut)
    E_vol_*.vtr   (optional, OPENEMS_DUMP3D=1) decimated full 3-D E-field volume
plus ``excitation.png`` (waveform sanity plot) in the project directory.

Run
---
    # use the project venv that has the cp313 openEMS wheels
    C:\Users\bryce\openems-venv\Scripts\python.exe 3D-waveform-openems.py

Environment flags
-----------------
    OPENEMS_QUICK=1   small/fast smoke run (6 ns delay, coarse mesh, small box)
    OPENEMS_DUMP3D=1  also dump the full 3-D E-field volume (large!)
    OPENEMS_OUT=path  override output directory
    CSXCAD_INSTALL_PATH=path   location of the extracted openEMS DLLs
"""

import os
import sys
import numpy as np

# ----------------------------------------------------------------------
# Locate the openEMS native DLLs (Windows). The CSXCAD/openEMS extension
# modules load CSXCAD.dll, openEMS.dll, hdf5.dll, fparser.dll, ... from this
# folder. Must be set BEFORE importing CSXCAD/openEMS.
# ----------------------------------------------------------------------
OPENEMS_HOME = os.environ.get("CSXCAD_INSTALL_PATH", r"C:\Users\bryce\openEMS")
os.environ.setdefault("CSXCAD_INSTALL_PATH", OPENEMS_HOME)
if hasattr(os, "add_dll_directory") and os.path.isdir(OPENEMS_HOME):
    os.add_dll_directory(OPENEMS_HOME)

from CSXCAD import ContinuousStructure          # noqa: E402
from openEMS import openEMS                      # noqa: E402

# ======================================================================
# Run mode
# ======================================================================
QUICK   = os.environ.get("OPENEMS_QUICK", "0") == "1"
DUMP3D  = os.environ.get("OPENEMS_DUMP3D", "0") == "1"

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
# Field dumps default to a LOCAL (non-OneDrive) folder: a transient solve writes
# hundreds of small .vtr files and you do not want OneDrive syncing every one.
DEFAULT_OUT = os.path.join(os.path.expanduser("~"), "openems_runs",
                           "bowtie_pulse_quick" if QUICK else "bowtie_pulse")
OUT_DIR = os.environ.get("OPENEMS_OUT", DEFAULT_OUT)

# ======================================================================
# Physical constants & signal definition
# ======================================================================
c0 = 299792458.0                 # speed of light (m/s)
f0 = 1.0e9                       # carrier frequency (Hz)
fmax = 1.6e9                     # highest freq of interest -> drives mesh & dt

# pulse timing (seconds)
t_delay = (6.0 if QUICK else 120.0) * 1e-9
t_rise  = 2.0e-9
t_width = 20.0e-9
t_fall  = 2.0e-9

# trapezoid breakpoints
t1 = t_delay                                   # start of rise
t2 = t_delay + t_rise                          # start of flat top
t3 = t_delay + t_rise + t_width                # start of fall
t4 = t_delay + t_rise + t_width + t_fall       # end of pulse
w0 = 2.0 * np.pi * f0

# openEMS custom excitation string (parsed by the C++ function parser, var = t,
# seconds). Numeric literals only (no 'pi') to stay parser-safe.
EXC = (
    "(if(t<{t1:.9e},0,"
    "if(t<{t2:.9e},(t-{t1:.9e})/{tr:.9e},"
    "if(t<{t3:.9e},1,"
    "if(t<{t4:.9e},({t4:.9e}-t)/{tf:.9e},0)))))"
    "*sin({w0:.10e}*(t-{t1:.9e}))"
).format(t1=t1, t2=t2, t3=t3, t4=t4, tr=t_rise, tf=t_fall, w0=w0)


def envelope(t):
    """Numpy reference trapezoid envelope (for the sanity plot / validation)."""
    t = np.asarray(t, dtype=float)
    e = np.zeros_like(t)
    m = (t >= t1) & (t < t2); e[m] = (t[m] - t1) / t_rise
    m = (t >= t2) & (t < t3); e[m] = 1.0
    m = (t >= t3) & (t < t4); e[m] = (t4 - t[m]) / t_fall
    return e


def source(t):
    """Numpy reference of the full excitation s(t) = env*carrier."""
    return envelope(t) * np.sin(w0 * (np.asarray(t, float) - t1))


# ======================================================================
# Geometry (millimetres; matches the analytical bowtie)
# ======================================================================
unit  = 1e-3            # all CSX coordinates are in mm
flare = 60.0            # half-width at the flared (open) end, along x
arm_L = 75.0            # arm length along z
gap   = 2.0             # central feed gap along z
feed_R = 50.0           # lumped feed-port resistance (ohm)

# air padding around the structure, and target free-space cell size
pad      = 60.0 if QUICK else 150.0
max_res  = c0 / fmax / unit / (12.0 if QUICK else 20.0)   # cells per wavelength
fine_res = 4.0 if QUICK else 2.0                          # near-structure cell (mm)

# simulation box extents (mm)
x_lim = flare + pad
y_lim = pad
z_lim = gap / 2.0 + arm_L + pad


def build():
    """Assemble CSX geometry, mesh, ports and dump boxes; return (FDTD, CSX)."""
    FDTD = openEMS(EndCriteria=1e-4)
    FDTD.SetCustomExcite(EXC, f0, fmax)
    FDTD.SetBoundaryCond(['PML_8'] * 6)            # open boundaries (free space)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    # ---- bowtie metal (two triangles in the y=0 plane) ----
    bowtie = CSX.AddMetal('bowtie')
    top = [[0.0,  flare, -flare],                  # x coords
           [gap / 2.0, gap / 2.0 + arm_L, gap / 2.0 + arm_L]]   # z coords
    bot = [[0.0,  flare, -flare],
           [-gap / 2.0, -(gap / 2.0 + arm_L), -(gap / 2.0 + arm_L)]]
    bowtie.AddPolygon(top, 'y', 0.0, priority=10)
    bowtie.AddPolygon(bot, 'y', 0.0, priority=10)

    # ---- feed: 50 ohm lumped port across the gap, oriented along z ----
    port = FDTD.AddLumpedPort(1, feed_R,
                              [0.0, 0.0, -gap / 2.0],
                              [0.0, 0.0,  gap / 2.0],
                              'z', 1.0, priority=5)

    # ---- mesh: fixed lines on the structure, then the air box, then smooth ----
    mesh.AddLine('x', [-flare, 0.0, flare, -x_lim, x_lim])
    mesh.AddLine('y', [0.0, -y_lim, y_lim])
    mesh.AddLine('z', [-(gap / 2.0 + arm_L), -gap / 2.0, 0.0,
                       gap / 2.0, gap / 2.0 + arm_L, -z_lim, z_lim])
    # pin the engine grid to the conductor edges (thirds rule) for accuracy
    FDTD.AddEdges2Grid(dirs='all', properties=bowtie)

    for d in ('x', 'y', 'z'):
        mesh.SmoothMeshLines(d, fine_res, ratio=1.5)   # refine near fixed lines
        mesh.SmoothMeshLines(d, max_res, ratio=1.4)    # grade out to free space

    # ---- field dumps for ParaView (E-field, time-domain, VTK) ----
    # dump_mode=2 -> cell interpolation (good for visualisation); file_type=0 -> VTK
    e_xz = CSX.AddDump('Ez_xz', dump_type=0, dump_mode=2, file_type=0)
    e_xz.AddBox([-x_lim, 0.0, -z_lim], [x_lim, 0.0, z_lim])      # y=0 (E-plane)

    e_xy = CSX.AddDump('Ez_xy', dump_type=0, dump_mode=2, file_type=0)
    e_xy.AddBox([-x_lim, -y_lim, 0.0], [x_lim, y_lim, 0.0])      # z=0 (H-plane)

    if DUMP3D:
        e_vol = CSX.AddDump('E_vol', dump_type=0, dump_mode=2, file_type=0,
                            sub_sampling=[2, 2, 2])
        e_vol.AddBox([-x_lim, -y_lim, -z_lim], [x_lim, y_lim, z_lim])

    return FDTD, CSX, mesh, port


def estimate_timesteps(mesh):
    """Estimate a Courant-stable dt from the final mesh and pick NrTS so the run
    spans the whole pulse plus ring-down. EndCriteria ends it earlier if quiet."""
    dl_min = []
    for d in ('x', 'y', 'z'):
        lines = np.array(mesh.GetLines(d)) * unit       # -> metres
        dl_min.append(np.min(np.diff(lines)))
    inv2 = sum((1.0 / dl) ** 2 for dl in dl_min)
    dt = 0.95 / (c0 * np.sqrt(inv2))                     # 0.95 * Courant limit
    t_total = t4 + (20.0e-9 if not QUICK else 12.0e-9)   # pulse end + ring-down
    nts = int(np.ceil(t_total / dt * 1.10))
    return dt, nts, t_total


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- excitation sanity plot (matplotlib, headless) ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        tt = np.linspace(0, t4 + 5e-9, 60000)
        plt.figure(figsize=(11, 4), facecolor="#0e141b")
        ax = plt.gca(); ax.set_facecolor("#0e141b")
        ax.plot(tt * 1e9, source(tt), color="#4cc9f0", lw=0.5, label="s(t)")
        ax.plot(tt * 1e9,  envelope(tt), color="#f4d35e", lw=1.6, label="envelope")
        ax.plot(tt * 1e9, -envelope(tt), color="#f4d35e", lw=1.6)
        ax.set_xlabel("time (ns)", color="#cfd8e3")
        ax.set_ylabel("amplitude", color="#cfd8e3")
        ax.set_title(f"Pulse-modulated {f0/1e9:.0f} GHz excitation  "
                     f"(delay {t_delay*1e9:.0f} ns, rise/fall {t_rise*1e9:.0f} ns, "
                     f"PW {t_width*1e9:.0f} ns)", color="#e8edf2")
        ax.tick_params(colors="#9fb0c0"); ax.grid(True, color="#22303c", lw=0.5)
        ax.legend(loc="upper left", facecolor="#16202b", labelcolor="#cfd8e3")
        png = os.path.join(PROJECT_DIR, "excitation.png")
        plt.savefig(png, dpi=140, facecolor="#0e141b", bbox_inches="tight")
        plt.close()
        print(f"[plot] wrote {png}")
    except Exception as e:                               # non-fatal
        print(f"[plot] skipped ({e})")

    FDTD, CSX, mesh, port = build()
    dt, nts, t_total = estimate_timesteps(mesh)
    FDTD.SetNumberOfTimeSteps(nts)

    nx = len(mesh.GetLines('x')); ny = len(mesh.GetLines('y'))
    nz = len(mesh.GetLines('z'))
    print("=" * 64)
    print(f" mode            : {'QUICK smoke test' if QUICK else 'FULL fidelity'}")
    print(f" excite string   : {EXC[:70]}...")
    print(f" mesh            : {nx} x {ny} x {nz} = {nx*ny*nz} nodes")
    print(f" max cell        : {max_res:.2f} mm  (lambda/{c0/fmax/(max_res*unit):.0f} @ {fmax/1e9:.1f} GHz)")
    print(f" est. timestep   : {dt*1e12:.3f} ps")
    print(f" timesteps (cap) : {nts}  (~{t_total*1e9:.0f} ns window)")
    print(f" output dir      : {OUT_DIR}")
    print("=" * 64)

    # CSX/geometry export for inspection in AppCSXCAD
    CSX.Write2XML(os.path.join(OUT_DIR, "geometry.xml"))

    FDTD.Run(OUT_DIR, cleanup=True, verbose=3)

    # ---- report what was produced ----
    vtr = [f for f in os.listdir(OUT_DIR) if f.endswith(".vtr")]
    groups = {}
    for f in vtr:
        key = f.rsplit("_", 1)[0]
        groups[key] = groups.get(key, 0) + 1
    print("\n[done] VTK time-series written to", OUT_DIR)
    for k, n in sorted(groups.items()):
        print(f"   {k}_*.vtr : {n} frames")
    print("\nOpen in ParaView:  File > Open > select e.g. 'Ez_xz_..vtr' "
          "(the '..' group), then Apply and play the animation.")


if __name__ == "__main__":
    sys.exit(main())
