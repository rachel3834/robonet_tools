# Code by Y. Tsapras, adapted by RA. Street

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm

# Parameters
I0, umin, tau, tmax = 1.0, -0.22, 30.0, 15.0

# Time array for the light curve
x = np.linspace(-60, 90, 1500)

# Source radii and angles
radii = np.linspace(0.01, 0.1, 10)  # fewer radii for clarity
angles = np.linspace(0, 2 * np.pi, 100)

# Set up the figure and subplots
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

# Initialize the plots
def init():
    ax1.clear()

# Animation function
def animate(j):
    ax1.clear()

    # Einstein ring
    ering = plt.Circle((0, 0), radius=1.0, color='gray', alpha=0.3)
    ax1.add_patch(ering)

    # Current source position
    point = j / 15.0 - 3

    # Plot source and images
    for rs in radii:
        xs = rs * np.cos(angles) + point
        ys = rs * np.sin(angles) + umin

        beta = np.sqrt(xs**2 + ys**2)
        phi = np.arctan2(xs, ys)

        theta_p = beta / 2 + np.sqrt((beta / 2)**2 + 1)
        theta_m = beta / 2 - np.sqrt((beta / 2)**2 + 1)

        # Major image
        ximp, yimp = theta_p * np.sin(phi), theta_p * np.cos(phi)
        # Minor image
        ximm, yimm = theta_m * np.sin(phi), theta_m * np.cos(phi)

        ax1.plot(xs, ys, 'b.', markersize=1, alpha=0.5)
        ax1.plot(ximp, yimp, 'g.', markersize=1, alpha=0.5)
        ax1.plot(ximm, yimm, 'g.', markersize=1, alpha=0.5)

    # Lens
    ax1.plot(0, 0, 'ko', markersize=5)
    ax1.set(xlim=(-3, 3), ylim=(-3, 3), xlabel=r'$x/R_E$', ylabel=r'$y/R_E$')
    ax1.set_aspect('equal')
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Compute magnification curve
    u = np.sqrt(umin**2 + (2 * (x - tmax) / tau)**2)
    A = (u**2 + 2) / (u * np.sqrt(u**2 + 4))

    # Increase label sizes
    ax1.tick_params(axis='both', labelsize=12)
    ax1.xaxis.label.set_size(14)
    ax1.yaxis.label.set_size(14)

# Create animation with tqdm progress bar
frames = 90
ani = animation.FuncAnimation(fig, animate, frames=tqdm(range(frames)), interval=100, init_func=init)
#ani.save('microlensing_event.mp4', writer='ffmpeg', fps=10)

# uncomment this to view it as soon as it is done
plt.show()

