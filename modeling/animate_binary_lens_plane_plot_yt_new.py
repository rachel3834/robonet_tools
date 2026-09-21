# Code by Y. Tsapras.
# Single-panel + faster artist-updates version.

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as poly
from tqdm import tqdm
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter, PillowWriter, writers

# ============================================================
# Parameters to be set by user
# ============================================================

# GM1 = mass of object 1 (as a % of total mass) e.g. 0.1
# GM2 = mass of object 2 (as a % of total mass) e.g  1-GM1 = 0.9
# D = half binary separation (between components) in Einstein radii
#     (midpoint at x=0, masses at x=-D and x=+D)
# cof1, cof2 = trajectory coefficients (trajectory of the
#              form Y = cof1 * X + cof2 [[use small values]])

GM1 = 0.1
GM2 = 1.0 - GM1
D = 0.5
cof1, cof2 = 0.0, -0.1

# Define the number of points (NPN) to use for the trajectory
# Trajectory from -XLM to +XLM R_E (sets the plot limits)
# frames is the number of frames to generate for the animation
XLM = 3
NPN = 4000
frames = 90

# Curves to display in plots: "none", "critical", "caustic", "both"
DISPLAY_CURVES = "both"

# User-set colors for the plots
COLOR_SOURCE   = "#FBB117"  # default source color
COLOR_CAUSTIC  = "#B83C08"  # default caustic color
COLOR_IMAGES   = "#FFA62F"  # default images color
COLOR_CRITICAL = "#1C1C1C"  # default critical curve color

# Define the number of source positions to use to fill the image
# angles sets how many points around the source circumference
# rnum sets how many radii to sample inside the source disk
angles = np.radians(np.arange(0, 360))
rnum = np.arange(0.0, 30.0) / 900.0 + 0.01

# Critical/caustic sampling
NR = 300000
DR = 1.0e-5

# Output filename stem
OUTSTEM = "binary_lens_plane"

# ============================================================
# Validation / normalization
# ============================================================

DISPLAY_CURVES = DISPLAY_CURVES.lower().strip()
valid_modes = {"none", "critical", "caustic", "both"}
if DISPLAY_CURVES not in valid_modes:
    raise ValueError(f"DISPLAY_CURVES must be one of {valid_modes}")

SHOW_CRITICAL = DISPLAY_CURVES in {"critical", "both"}
SHOW_CAUSTIC = DISPLAY_CURVES in {"caustic", "both"}

GM1 = float(GM1)
GM2 = float(GM2)
mtot = GM1 + GM2
if mtot <= 0:
    raise ValueError("GM1 + GM2 must be > 0")
GM1 /= mtot
GM2 /= mtot

# Precompute trig arrays for speed
COS_ANG = np.cos(angles)
SIN_ANG = np.sin(angles)

# ============================================================
# Binary lens solver
# ============================================================
def bin_ima(GM1, GM2, D, XS, YS):
    """
    Solve binary gravitational lens equation for image positions, magnifications, and parities.

    IMPORTANT COORDINATE NOTE:
    - Masses/images/critical curves belong to the lens plane.
    - Source/caustics belong to the source plane.
    In a *single-panel* visualization, we are deliberately overlaying quantities from
    two different planes.

    The lens-plane (x,y) origin is the midpoint of the binary (x-axis along the binary axis).
    The source-plane uses a parallel (u,v) system.

    Parameters:
        GM1, GM2: Mass fractions of lens components (GM1 + GM2 = 1)
        D: Half binary separation in Einstein radii
        XS, YS: Source position coordinates (units of Einstein radius)

    Returns:
        result: (5,4) array [XI, YI, AI, IP] with:
                XI, YI: Image x,y coordinates
                AI: Magnifications of images
                IP: Image parities (-1, 0, +1)
    """

    # Complex representation of lens positions and source position
    Z1 = complex(-D, 0.0)
    Z2 = complex(D, 0.0)
    ZS = complex(XS, YS)
    ZSC = complex(XS, -YS)
    ZSS = ZS * ZSC

    # Intermediate quantities for polynomial coefficients
    HSM = (GM1 + GM2) / 2.0
    HDM = (GM2 - GM1) / 2.0
    Z1_sq = Z1**2

    # Coefficients of the 5th-order polynomial (from Schneider & Weiss)
    ZC = np.zeros(6, dtype=complex)
    ZC[0] = Z1_sq * (4*HDM**2*ZS + Z1*(4*HSM*HDM + 4*HDM*ZSS +
             Z1*(2*HSM*ZSC + ZSS*ZSC - Z1*(2*HDM + Z1*ZS))))
    ZC[1] = -Z1*(8*HSM*HDM*ZS + Z1*(4*(HDM**2 + HSM**2) + 4*HSM*ZSS +
             Z1*(4*HDM*ZSC + Z1*(ZSC**2 - Z1_sq))))
    ZC[2] = 4*HSM**2*ZS + Z1*(4*HSM*HDM - 4*HDM*ZSS +
             Z1*(-2*ZSS*ZSC + Z1*(4*HDM + 2*ZS*Z1)))
    ZC[3] = 4*HSM*ZSS + Z1*(4*HDM*ZSC + Z1*(2*ZSC**2 - 2*Z1_sq))
    ZC[4] = ZSC*(ZSS - 2*HSM) - Z1*(2*HDM + ZS*Z1)
    ZC[5] = Z1_sq - ZSC**2

    # Solve polynomial for image positions
    ZI = poly.polyroots(ZC)

    # Prepare result array
    result = np.zeros((5, 4))
    result[:, 0] = ZI.real  # XI
    result[:, 1] = ZI.imag  # YI

    # Calculate lens equation and magnifications
    ZIC = np.conjugate(ZI)
    ZDC1 = ZIC - Z1
    ZDC2 = ZIC - Z2
    ZB = GM1/ZDC1 + GM2/ZDC2
    ZE = ZI - ZB
    ZD = GM1/ZDC1**2 + GM2/ZDC2**2

    AJ = 1.0 / (1.0 - np.abs(ZD)**2)
    result[:, 2] = AJ.real  # AI (magnifications)

    # Determine image parities
    EP = 1e-2
    for i in range(5):
        separation = np.abs(ZE[i] - ZS)
        if separation <= EP:
            if AJ[i] > 0:
                IP = 1
            elif AJ[i] < 0:
                IP = -1
            else:
                IP = 0
        else:
            IP = 0
        result[i, 3] = IP

    return result

# ============================================================
# Critical curves and caustics
# ============================================================

def compute_critical_and_caustic(GM1, GM2, D, NR=300000, DR=1e-5):
    """
    Compute the critical curves (lens plane) and the corresponding caustics
    (source plane) for a binary lens.

    Parameters:
        GM1 : float
            Mass fraction of the first lens component.
        GM2 : float
            Mass fraction of the second lens component.
            (GM1 + GM2 = 1)
        D : float
            Half the binary separation in units of the Einstein radius.
            The two lens masses are therefore located at x = -D and x = +D.
        NR : int, optional
            Number of radial sampling points used to estimate the critical
            curves and caustics. Larger values give smoother curves, but
            increase computation time.
        DR : float, optional
            Radial step size used in the sampling of the critical curves.
            Smaller values increase the sampling density, but also increase
            computation time.

    Returns:
        critical_x : ndarray
            x-coordinates of the critical curves in the lens plane.
        critical_y : ndarray
            y-coordinates of the critical curves in the lens plane.
        caustic_x : ndarray
            x-coordinates of the caustics in the source plane.
        caustic_y : ndarray
            y-coordinates of the caustics in the source plane.
    """
    xcrit_chunks, ycrit_chunks = [], []
    xcaus_chunks, ycaus_chunks = [], []

    m1, m2 = GM1, GM2
    ip = -1

    D2 = D * D

    # Perform repeat calculations with masses swapped over
    for _ in range(2):
        IR = np.arange(1, NR)
        R = IR * DR
        R2 = R * R
        R4 = R2 * R2

        m1s = m1 * m1
        m2s = m2 * m2
        mxm = m1 * m2

        R2P = R2 + 4.0 * D2
        R4M = R4 - m2s

        # Polynomial coeffs A, B, C as in Schneider & Weiss (eqn 9b)
        A = 16.0 * D2 * R2 * (R4M - mxm)
        B = 8.0 * R * D * (mxm * R2 - (R2 + 4.0 * D2) * R4M)
        C = (R2P * R2P) * R4M - m1s * R4 - 2.0 * mxm * R2 * (R2 + 4.0 * D2)
        C = C + 16.0 * m1 * m2 * D2 * R2

        # Determinant
        DT = B * B - 4.0 * A * C
        mask = DT >= 0.0
        if not np.any(mask):
            m1, m2 = m2, m1
            ip = -ip
            continue

        DT_sel = DT[mask]
        R_sel = R[mask]
        R2_sel = R2[mask]
        B_sel = B[mask]
        A_sel = A[mask]

        C1 = (-B_sel + np.sqrt(DT_sel)) / (2.0 * A_sel)
        C2 = (-B_sel - np.sqrt(DT_sel)) / (2.0 * A_sel)

        for Ck, Rk, R2k in ((C1, R_sel, R2_sel), (C2, R_sel, R2_sel)):
            mask_c = np.abs(Ck) <= 1.0
            if not np.any(mask_c):
                continue

            Rm = Rk[mask_c]
            R2m = R2k[mask_c]
            Ckm = Ck[mask_c]

            X = Rm * Ckm
            S = np.sqrt(1.0 - Ckm * Ckm)
            Y = Rm * S

            # Critical curve in lens plane
            Xc = X - D
            Yc = Y
            Xcrit = ip * Xc
            Ycrit = ip * Yc

            # Append top & bottom halves
            xcrit_chunks.append(np.concatenate([Xcrit, Xcrit]))
            ycrit_chunks.append(np.concatenate([Ycrit, -Ycrit]))

            # Map critical curves to caustics in source plane
            RD2 = R2m + 4.0 * D2
            UP1 = X - 2.0 * D
            UP2 = X
            DN1 = RD2 - 4.0 * D * X
            DN2 = R2m

            XC = X - m1 * (UP1 / DN1) - m2 * (UP2 / DN2)
            XC = XC - D
            YC = Y * (1.0 - m1 / DN1 - m2 / DN2)

            Xcaus = ip * XC
            Ycaus = ip * YC

            xcaus_chunks.append(np.concatenate([Xcaus, Xcaus]))
            ycaus_chunks.append(np.concatenate([Ycaus, -Ycaus]))

        # Swap the masses and repeat
        m1, m2 = m2, m1
        ip = -ip

    def cat(chunks):
        return np.concatenate(chunks) if chunks else np.array([])

    return cat(xcrit_chunks), cat(ycrit_chunks), cat(xcaus_chunks), cat(ycaus_chunks)

# ============================================================
# Trajectory (source plane)
# ============================================================

# X,Y coords of source center trajectory
XSA = np.linspace(-XLM, XLM, NPN)
YSA = cof1 * XSA + cof2

# Compute static curves only if needed
if SHOW_CRITICAL or SHOW_CAUSTIC:
    critical_x, critical_y, caustic_x, caustic_y = compute_critical_and_caustic(
        GM1, GM2, D, NR=NR, DR=DR
    )
else:
    critical_x = critical_y = caustic_x = caustic_y = np.array([])

# ============================================================
# Precompute source fill offsets (speed)
# ============================================================

# We precompute a fixed cloud of points that fills the source disk,
# centered at (0,0). For each animation frame, we just translate it
# to the current (xs, ys) location. This avoids repeatedly constructing
# circles point-by-point.
#
# Note: This is *only* for drawing the source. We still solve the lens
# equation for each sampled point to obtain the images.
dx_list = []
dy_list = []
for rs in rnum:
    dx_list.append(rs * COS_ANG)
    dy_list.append(rs * SIN_ANG)

# Flatten into (Nsrc,) arrays for convenient translation
SRC_DX = np.concatenate(dx_list)
SRC_DY = np.concatenate(dy_list)

# ============================================================
# Figure / animation (single panel)
# ============================================================

fig, ax = plt.subplots(1, 1, figsize=(6, 6))

# We will NOT clear() the axis in each frame (slow).
# Instead we draw static artists once and update only dynamic artists.

# --- Static axis styling ---
ax.set_xlim(-XLM, XLM)
ax.set_ylim(-2.0, 2.0)
ax.set_xlabel(r"$x/R_E$")
ax.set_ylabel(r"$y/R_E$")
ax.grid(True)

# --- Static elements: masses (lens plane) ---
# Note: marker size proportional to mass for a bit of intuition.
mass_sizes = 40.0 + 120.0 * np.array([GM1, GM2])
mass_artist = ax.scatter(
    [-D, D], [0.0, 0.0],
    c=COLOR_CRITICAL,
    s=mass_sizes,
    zorder=5
)

# --- Static element: source trajectory (source plane) ---
traj_artist, = ax.plot(XSA, YSA, "k-", lw=0.5, alpha=0.7, zorder=1)

# --- Static elements: critical curves / caustics (optional) ---
critical_artist = None
caustic_artist = None

if SHOW_CRITICAL and critical_x.size:
    # Use a pixel marker for speed
    critical_artist, = ax.plot(
        critical_x, critical_y,
        linestyle="none",
        marker=",",
        color=COLOR_CRITICAL,
        alpha=0.9,
        zorder=2
    )

if SHOW_CAUSTIC and caustic_x.size:
    caustic_artist, = ax.plot(
        caustic_x, caustic_y,
        linestyle="none",
        marker=",",
        color=COLOR_CAUSTIC,
        alpha=0.9,
        zorder=2
    )

# --- Dynamic artists (updated every frame) ---
# Source outline/fill (source plane) and images (lens plane)
source_artist, = ax.plot([], [], linestyle="none", marker=",", markersize=1,
                         color=COLOR_SOURCE, zorder=4)

source_center_artist, = ax.plot([], [], marker="o", ms=3, linestyle="none",
                                color=COLOR_SOURCE, zorder=6)

images_artist, = ax.plot([], [], linestyle="none", marker=",", markersize=1,
                         color=COLOR_IMAGES, zorder=3)

def init():
    # Nothing expensive here—just ensure artists start empty for blitting.
    source_artist.set_data([], [])
    source_center_artist.set_data([], [])
    images_artist.set_data([], [])
    artists = [source_artist, source_center_artist, images_artist, mass_artist, traj_artist]
    if critical_artist is not None:
        artists.append(critical_artist)
    if caustic_artist is not None:
        artists.append(caustic_artist)
    return artists

def animate(i):
    # Current source center along the trajectory
    idx = i * NPN // frames
    xs, ys = XSA[idx], YSA[idx]

    # Translate the precomputed source fill cloud to (xs, ys)
    src_x = xs + SRC_DX
    src_y = ys + SRC_DY

    # Compute images for each sampled source point
    # NOTE: this remains the main computational cost.
    img_x = []
    img_y = []

    # Tight loop: one bin_ima per sampled point.
    # Small overhead reduction by local-binding frequently used names.
    _bin_ima = bin_ima
    _GM1, _GM2, _D = GM1, GM2, D

    for cx, cy in zip(src_x, src_y):
        images = _bin_ima(_GM1, _GM2, _D, float(cx), float(cy))
        valid = images[images[:, 3] != 0]
        if valid.size:
            img_x.extend(valid[:, 0])
            img_y.extend(valid[:, 1])

    # Update dynamic artists in-place (fast; no redraw of static content)
    source_artist.set_data(src_x, src_y)
    source_center_artist.set_data([xs], [ys])
    images_artist.set_data(img_x, img_y)

    # Return artists for blitting
    artists = [source_artist, source_center_artist, images_artist]
    return artists

# Create animation
ani = animation.FuncAnimation(
    fig,
    animate,
    init_func=init,
    frames=tqdm(range(frames)),
    interval=100,
    blit=True,          # key speed-up (works best when you don't clear/redraw axes)
    repeat=False
)

# Save as MP4 if ffmpeg exists, else GIF
if writers.is_available("ffmpeg"):
    ani.save(f"{OUTSTEM}.mp4", writer=FFMpegWriter(fps=10))
    print(f"Animation complete: '{OUTSTEM}.mp4'")
else:
    ani.save(f"{OUTSTEM}.gif", writer=PillowWriter(fps=10))
    print(f"ffmpeg not found; saved '{OUTSTEM}.gif' instead")

plt.close(fig)
