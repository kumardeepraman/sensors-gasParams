import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# Publication-quality settings
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 21,

    "figure.subplot.left": 0.12,
    "figure.subplot.right": 0.90,
    "figure.subplot.bottom": 0.15,
    "figure.subplot.top": 0.90,

    "xtick.labelsize": 21,
    "ytick.labelsize": 21,

    "legend.fontsize": 17,
})

# ============================================================
# Colors
# ============================================================

root_colors = {
    'default': '#0099ff',
    '75um':    '#009900',
    '50um':    '#cc3333',
    'UV':      '#993399'
}

# ============================================================
# Read wavelength and PDE curves
# ============================================================

wavelenghts_old = np.loadtxt(
    'PDE_old.txt',
    usecols=0,
    skiprows=2
)

# Default SiPM
PDE_old = np.loadtxt(
    'PDE_old.txt',
    usecols=1,
    skiprows=2
)

# Hamamatsu SiPMs
wavelenghts_pitch = np.loadtxt(
    'PDE_hamamatsu.txt',
    usecols=0,
    skiprows=2
)

# 50 um
PDE_UV_red = np.loadtxt(
    'PDE_hamamatsu.txt',
    usecols=1,
    skiprows=2
)

# 75 um
PDE_UV_green = np.loadtxt(
    'PDE_hamamatsu.txt',
    usecols=2,
    skiprows=2
)

# 25 um (not used)
PDE_UV_blue = np.loadtxt(
    'PDE_hamamatsu.txt',
    usecols=3,
    skiprows=2
)

# UV extended SiPM
wavelenghts_UV = np.loadtxt(
    'PDE_75um.txt',
    usecols=0,
    skiprows=2
)

PDE_UV = np.loadtxt(
    'PDE_75um.txt',
    usecols=1,
    skiprows=2
)

# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(7.5, 5), dpi=100)

plt.plot(
    wavelenghts_old,
    PDE_old,
    color=root_colors['default'],
    label='default',
    marker='.'
)

plt.plot(
    wavelenghts_pitch,
    PDE_UV_red,
    color=root_colors['50um'],
    label=r'50 $\mu$m',
    marker='.'
)

plt.plot(
    wavelenghts_pitch,
    PDE_UV_green,
    color=root_colors['75um'],
    label=r'75 $\mu$m',
    marker='.'
)

plt.plot(
    wavelenghts_UV,
    PDE_UV,
    color=root_colors['UV'],
    label='UV',
    marker='.'
)

# 280 nm reference line
plt.plot(
    [280, 280],
    [0, 1],
    linestyle='--',
    color='black',
    label='280 nm'
)

# ============================================================
# Axis labels and limits
# ============================================================

plt.xlabel(
    r"$\lambda$ [nm]",
    loc="right"
)

plt.ylabel(
    "PDE",
    loc="top",
    labelpad=-0.5
)

plt.ylim(0, 0.6)

# ============================================================
# Grid
# ============================================================

plt.grid()

# ============================================================
# Compact legend
# ============================================================

plt.legend(
    loc='best',
    fontsize=17,
    frameon=True,
    framealpha=1.0,
    facecolor='white',
    edgecolor='black',

    # Reduce size of legend box
    borderpad=0.25,
    labelspacing=0.25,
    handlelength=1.4,
    handletextpad=0.4,
    columnspacing=0.5
)

# ============================================================
# Save
# ============================================================

plt.savefig(
    "PDE_comparison.pdf",
    bbox_inches='tight'
)

plt.show()
