# Code by Y. Tsapras.

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as poly
from tqdm import tqdm
import matplotlib.animation as animation

# Parameters to be set by user

# GM1 = mass of object 1 (as a % of total mass) e.g. 0.1
# GM2 = mass of object 2 (as a % of total mass) e.g  0.9
# D = half binary separation (between components) in Einstein radii
# cof1, cof2 = trajectory coefficients (trajectory of the 
#              form Y = cof1 * X + cof2 [[use small values]])

GM1, GM2, D, cof1, cof2 = 0.1, 0.9, 0.5, 0.0, -0.1

# Define the number of points (NPN) to use for the trajectory
# Trajectory from -XLM to +XLM R_E (sets the plot limits)
# frames is the number of frames to generate for the animation
XLM, NPN, frames = 3, 4000, 90

# Define the number of source positions to use to fill the image
angles = np.radians(np.arange(0, 360))
rnum = np.arange(0., 30.) / 900.0 + 0.01

#########################################################################
def bin_ima(GM1, GM2, D, XS, YS):
    """
    Solve binary gravitational lens equation for image positions, magnifications, and parities.

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

#########################################################################

# Pre-calculate trajectory
# X,Y coords of source and magnification ASA
XSA = np.linspace(-XLM, XLM, NPN)
YSA = -cof1 * XSA + cof2
ASA = np.array([np.sum(bin_ima(GM1, GM2, D, x, y)[:, 2] * bin_ima(GM1, GM2, D, x, y)[:, 3]) for x, y in zip(XSA, YSA)])

# Critical and caustic curves calculation
# Initialize arrays holding critical curve and caustic points
critical_points_x1, critical_points_y1 = np.array([]), np.array([])
critical_points_x2, critical_points_y2 = np.array([]), np.array([])
caustic_points_x1, caustic_points_y1  = np.array([]), np.array([])
caustic_points_x2, caustic_points_y2  = np.array([]), np.array([])

# Set plot limits
XLM = 3 # Trajectory from -XLM to +XLM R_E

IP  = -1
GM2 = 1 - GM1
RE1, ZR1 = np.sqrt(GM1), np.sqrt(GM1)
RE2, ZR2 = np.sqrt(GM2), np.sqrt(GM2)
D2 = D * D
D4 = D2 * D2

NR = 300000 # Points to use to plot critical curves and caustics
DR = 0.00001 # Used to define the sampling density of the caustics

# Estimate criticals and caustics
# Perform repeat calculations with masses swapped over
for IQ in range(0,2):
    IR = np.arange(1, NR)
    # Define some variables that are used repeatedly
    R = IR * DR
    R2 = R * R
    R4 = R2 * R2
    GM1S = GM1 * GM1
    GM2S = GM2 * GM2
    GMXM = GM1 * GM2
    R2P = R2 + 4 * D2
    R4M = R4 - GM2S
    # Polynomial coeffs A, B, C
    # as described in Schneider & Weiss 1986 (eqn 9b)
    A = 16 * D2 * R2 * (R4M - GMXM)
    B = 8 * R * D * (GMXM * R2 - (R2 + 4 * D2) * R4M)
    C = (R2P * R2P) * R4M - GM1S * R4 \
       - 2 * GMXM * R2 * (R2 + 4 * D2)
    C = C + 16 * GM1 * GM2 * D2 * R2
    # Calculate the determinant DT
    DT = B * B - 4 * A * C
    # When the determinant is >= 0 calculate the values of cos(theta) C1, C2
    DTge0 = np.compress(DT >= 0, DT)
    # Modify relevant arrays accordingly
    R_DTgt0 = np.compress(DT >= 0, R)
    R2_DTgt0 = np.compress(DT >= 0, R2)
    B_DTgt0 = np.compress(DT >= 0, B)
    A_DTgt0 = np.compress(DT >= 0, A)
    # Evaluate C1, C2
    C1 = (-B_DTgt0 + np.sqrt(DTge0))/(2 * A_DTgt0)
    C2 = (-B_DTgt0 - np.sqrt(DTge0))/(2 * A_DTgt0)
    # When abs(C1) or abs(C2) are <=1 append point to critical and caustic curves
    # Isolate the indexes of where the absolute values of C1 and C2 are <= 1
    idxC1le1 = ((np.abs(C1) <= 1)).nonzero()
    idxC2le1 = ((np.abs(C2) <= 1)).nonzero()
    # Set up subsamples from relevant arrays
    Rmod2_C1 = R_DTgt0[idxC1le1]
    R2mod2_C1 = R2_DTgt0[idxC1le1]
    Rmod2_C2 = R_DTgt0[idxC2le1]
    R2mod2_C2 = R2_DTgt0[idxC2le1]
    X1, X2 = Rmod2_C1 * C1[idxC1le1], Rmod2_C2 * C2[idxC2le1]
    S1, S2 = np.sqrt(1.0 - C1[idxC1le1] * C1[idxC1le1]), \
	    np.sqrt(1.0 - C2[idxC2le1] * C2[idxC2le1])
    Y1, Y2 = Rmod2_C1 * S1, Rmod2_C2 * S2
    X_C1, X_C2  = X1 - D, X2 - D
    Y_C1, Y_C2  = Y1, Y2
    PXCRIT_C1, PYCRIT_C1 = IP * X_C1, IP * Y_C1
    PXCRIT_C2, PYCRIT_C2 = IP * X_C2, IP * Y_C2
    # Append to critical points 
    critical_points_x1 = np.append(critical_points_x1, PXCRIT_C1)
    critical_points_x1 = np.append(critical_points_x1, PXCRIT_C2)
    critical_points_x2 = np.append(critical_points_x2, PXCRIT_C1)
    critical_points_x2 = np.append(critical_points_x2, PXCRIT_C2)
    critical_points_y1 = np.append(critical_points_y1, PYCRIT_C1)
    critical_points_y1 = np.append(critical_points_y1, PYCRIT_C2)
    critical_points_y2 = np.append(critical_points_y2, -PYCRIT_C1)
    critical_points_y2 = np.append(critical_points_y2, -PYCRIT_C2)
    # Set up arrays to use for mapping to caustics
    RD2_C1, RD2_C2 = R2mod2_C1 + 4 * D2,  R2mod2_C2 + 4 * D2
    UP1_C1, UP1_C2 = X1 - 2 * D, X2 - 2 * D
    UP2_C1, UP2_C2 = X1, X2
    DN1_C1, DN1_C2 = RD2_C1 - 4 * D * X1, RD2_C2 - 4 * D * X2
    DN2_C1, DN2_C2 = R2mod2_C1, R2mod2_C2
    # Map critical curves to caustics
    XC_C1, XC_C2 = X1 - GM1 * (UP1_C1/DN1_C1) - GM2 * (UP2_C1/DN2_C1), \
    	  X2 - GM1 * (UP1_C2/DN1_C2) - GM2 * (UP2_C2/DN2_C2)
    XC_C1, XC_C2 = XC_C1 - D, XC_C2 - D
    YC_C1, YC_C2 = Y1 * (1.0 - GM1/DN1_C1 - GM2/DN2_C1), \
		  Y2 * (1.0 - GM1/DN1_C2 - GM2/DN2_C2)
    PXCAUS_C1, PXCAUS_C2 = IP * XC_C1, IP * XC_C2
    PYCAUS_C1, PYCAUS_C2 = IP * YC_C1, IP * YC_C2
    # Append to caustic points
    caustic_points_x1 = np.append(caustic_points_x1, PXCAUS_C1)
    caustic_points_x1 = np.append(caustic_points_x1, PXCAUS_C2)
    caustic_points_x2 = np.append(caustic_points_x2, PXCAUS_C1)
    caustic_points_x2 = np.append(caustic_points_x2, PXCAUS_C2)
    caustic_points_y1 = np.append(caustic_points_y1, PYCAUS_C1)
    caustic_points_y1 = np.append(caustic_points_y1, PYCAUS_C2)
    caustic_points_y2 = np.append(caustic_points_y2, -PYCAUS_C1)
    caustic_points_y2 = np.append(caustic_points_y2, -PYCAUS_C2)
    # Swap the masses and repeat the calculation
    GM0 = GM1
    GM1 = GM2
    GM2 = GM0
    IP = -IP

# Set up figure and axes
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

# Animation initialization
def init():
    ax1.set_xlim(-XLM, XLM)
    ax1.set_ylim(-2, 2)
    ax1.grid(True)
    ax1.set_xlabel(r'$x/R_E$', fontsize=12)
    ax1.set_ylabel(r'$y/R_E$', fontsize=12)

    return []

# Animation update function
def animate(i):
    ax1.clear()

    # orange #cd730e

    # Static elements
    ax1.plot([D, -D], [0, 0], 'ko')
    #ax1.plot(critical_points_x1, critical_points_y1, c='#7a0ecd', marker='.', markersize=1, ls='none')
    #ax1.plot(caustic_points_x1, caustic_points_y1, c='#1d8505', marker='.', markersize=1, ls='none')
    ax1.plot(XSA, YSA, 'k-', lw=0.5)

    # Current source position
    idx = i * NPN // frames
    xs, ys = XSA[idx], YSA[idx]

    # Source and image positions
    src_x, src_y, img_x, img_y = [], [], [], []

    for rs in rnum:
        circle_x = xs + rs * np.cos(angles)
        circle_y = ys + rs * np.sin(angles)
        src_x.extend(circle_x)
        src_y.extend(circle_y)

        for cx, cy in zip(circle_x, circle_y):
            images = bin_ima(GM1, GM2, D, cx, cy)
            valid_images = images[images[:, 3] != 0]
            img_x.extend(valid_images[:, 0])
            img_y.extend(valid_images[:, 1])

    ax1.plot(src_x, src_y, c='#cd730e', marker=',', markersize=1, ls='none')
    ax1.plot(img_x, img_y, c='#e2a911', marker=',', markersize=1, ls='none')
    ax1.set_xlabel('$R_{\\rm{E}}$')
    ax1.set_ylabel('$R_{\\rm{E}}$')
    ax1.grid(True)
    ax1.set_xlim(-XLM, XLM)
    ax1.set_ylim(-2, 2)

    return []

# Create animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=tqdm(range(frames)), interval=100)

# Save as MP4
ani.save('binary_lens_plane.mp4', writer='ffmpeg', fps=10)
#plt.show()
plt.close(fig)
print("Animation complete: 'binary_microlensing_event2.mp4'")

