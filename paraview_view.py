#!/usr/bin/env pvpython
r"""
paraview_view.py
================
Load the openEMS E-field time-series produced by ``3D-waveform-openems.py`` and
set up an animated 3-D view of the pulse-modulated 1 GHz bowtie radiation.

Two ways to use it
------------------
1. Headless / batch (renders a movie + frames, no GUI needed):
       "C:\Program Files\ParaView 5.x\bin\pvpython.exe" paraview_view.py
   (or pass the run directory:  pvpython paraview_view.py C:\path\to\run )

2. Inside the ParaView GUI:
       Tools > Python Shell, then:
           exec(open(r"C:\...\paraview_view.py").read())
   This builds the pipeline in the active view so you can interact/scrub.

What it builds
--------------
* The Ez_xz (y=0, E-plane) field series, colored by the SIGNED z-component of E
  with a symmetric diverging colormap -> you literally see wave crests/troughs of
  the 1 GHz carrier propagate outward and the 20 ns pulse pass through.
* A second layer for the bowtie geometry outline (from geometry.xml if present).
* A time annotation and color bar.
* Camera looking down the +y axis onto the x-z plane.

Switch ``COLOR_BY`` below to ('E-Field','Magnitude') to see the pulse *envelope/
intensity* instead of the signed wavefronts.
"""

import os
import sys
import glob

# --------------------------------------------------------------------------
# Resolve the run directory (matches 3D-waveform-openems.py defaults)
# --------------------------------------------------------------------------
DEFAULT_OUT = os.path.join(os.path.expanduser("~"), "openems_runs", "bowtie_pulse")
RUN_DIR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("OPENEMS_OUT", DEFAULT_OUT)

PLANE   = "Ez_xz"                       # "Ez_xz" (E-plane) or "Ez_xy" (H-plane)
COLOR_BY = ("E-Field", "Z")             # signed wavefronts;  use "Magnitude" for envelope
MOVIE    = os.path.join(RUN_DIR, "bowtie_pulse.avi")
FRAMEDIR = os.path.join(RUN_DIR, "frames")

frames = sorted(glob.glob(os.path.join(RUN_DIR, f"{PLANE}_*.vtr")))
if not frames:
    raise SystemExit(f"No {PLANE}_*.vtr files in {RUN_DIR!r}. "
                     f"Run 3D-waveform-openems.py first (or pass the run dir).")

from paraview.simple import *           # noqa: E402,F403
paraview_version = GetParaViewVersion() if "GetParaViewVersion" in dir() else "?"
print(f"ParaView {paraview_version} | {len(frames)} frames | {RUN_DIR}")

# --------------------------------------------------------------------------
# Reader: a numbered .vtr list is read as one time-series source
# --------------------------------------------------------------------------
field = XMLRectilinearGridReader(registrationName=PLANE, FileName=frames)
field.PointArrayStatus = ["E-Field"]

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 900]
view.Background = [0.055, 0.078, 0.106]            # dark slate to match the project
view.OrientationAxesVisibility = 1

disp = Show(field, view)
disp.Representation = "Surface"
ColorBy(disp, ("POINTS",) + COLOR_BY)
disp.RescaleTransferFunctionToDataRange(True, False)

# symmetric diverging colormap for signed field; rescale to a robust symmetric range
lut = GetColorTransferFunction("E-Field")
if COLOR_BY[1] == "Magnitude":
    ApplyPreset("Inferno (matplotlib)", True)
else:
    ApplyPreset("Cool to Warm (Extended)", True)
    # symmetric range so 0 = white; pick a value well below the source spike
    di = field.GetDataInformation()
    arr = di.GetPointDataInformation().GetArrayInformation("E-Field")
    if arr is not None:
        comp = {"X": 0, "Y": 1, "Z": 2}.get(COLOR_BY[1], 0)
        rng = arr.GetComponentRange(comp)
        a = 0.25 * max(abs(rng[0]), abs(rng[1]), 1e-9)   # clip the feed spike
        lut.RescaleTransferFunction(-a, a)

# color bar
bar = GetScalarBar(lut, view)
bar.Title = f"E_{COLOR_BY[1]}  (V/m)"
bar.ComponentTitle = ""
disp.SetScalarBarVisibility(view, True)

# --------------------------------------------------------------------------
# Optional: overlay the bowtie geometry outline
# --------------------------------------------------------------------------
geo = os.path.join(RUN_DIR, "geometry.xml")
# (openEMS geometry.xml is a CSX file, not directly a ParaView reader; the metal
#  shows up as the zero-field notch in the dump. Left as a hook for a future
#  exported STL if desired.)

# --------------------------------------------------------------------------
# Camera: look down +y onto the x-z plane
# --------------------------------------------------------------------------
view.ResetCamera()
cam = GetActiveCamera()
cam.SetPosition(0, -1, 0)
cam.SetFocalPoint(0, 0, 0)
cam.SetViewUp(0, 0, 1)
view.ResetCamera()

# time annotation
try:
    ann = AnnotateTimeFilter(field)
    ann.Format = "t-step %.0f"
    ad = Show(ann, view); ad.Color = [0.9, 0.9, 0.95]
except Exception:
    pass

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
Render()

# --------------------------------------------------------------------------
# Batch export (only when run via pvpython, not in the GUI shell)
# --------------------------------------------------------------------------
if __name__ == "__main__" and sys.argv[0].endswith(".py"):
    try:
        os.makedirs(FRAMEDIR, exist_ok=True)
        SaveAnimation(os.path.join(FRAMEDIR, "frame.png"), view,
                      ImageResolution=[1280, 900], FrameRate=25)
        print(f"[ok] wrote PNG frames to {FRAMEDIR}")
        try:
            SaveAnimation(MOVIE, view, ImageResolution=[1280, 900], FrameRate=25)
            print(f"[ok] wrote movie {MOVIE}")
        except Exception as e:
            print(f"[movie skipped: {e}]")
    except Exception as e:
        print(f"[batch export skipped: {e}] -- pipeline is set up for interactive use.")
