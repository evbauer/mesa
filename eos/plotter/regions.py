#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.style.use("mesa_eos_regions.mplstyle")


def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = data[0, :, 0]
    Xran = data[:, 0, 1]
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran


with open("eos_plotter.dat") as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()

# overwrite with fancier labels
xlabel = r"$\log(\rho/{\rm g\,cm^{-3}})$"
ylabel = r"$\log(T/{\rm K})$"

eosDT, Yran, Xran = parse("eos_plotter.dat")

apjcolwidth = 3.38
# set up plot and labels
# fig, ax = plt.subplots(figsize=(apjcolwidth,apjcolwidth*4./5.)) # for paper figures
fig, ax = plt.subplots(figsize=(5, 4))  # for website pngs

## for solar composition plot
title = r"MESA EOS Regions ($X=0.7$, $Z=0.02$)"
ax.text(4.0, 2.3, "no coverage", fontsize="small")

## for Z = 1 plot
# title = r'MESA EOS Regions ($X=0.0$, $Z=1.0$)'
# ax.text(4.8,2.3,'no coverage',fontsize='small')

ax.set_title(title)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(Xran.min(), Xran.max())
ax.set_ylim(Yran.min(), Yran.max())

# set up color map (slightly customized to make Skye blue)
# region codes from the Fortran plotter: HELM=1, OPAL/SCVH=2, FreeEOS=3, PC=4,
# Skye=5, CMS=6, ideal=7, WDM_AV=8.  The bins below lump PC+Skye and CMS+ideal,
# and give WDM_AV (the warm dense matter / additive volume component) its own bin.
accent = np.array(mpl.cm.get_cmap("Accent").colors)  # 8 colors, entries editable
my_colors = np.array(
    [
        accent[7],  # blend (code 0)
        accent[0],  # HELM (1)
        accent[1],  # OPAL/SCVH (2)
        accent[2],  # FreeEOS (3)
        accent[4],  # Skye (PC+Skye, 4-5)
        accent[5],  # ideal (CMS+ideal, 6-7)
        accent[6],  # WDM_AV (8)
    ]
)
cmap = colors.ListedColormap(my_colors)
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 5.5, 7.5, 8.5]
norm = colors.BoundaryNorm(bounds, cmap.N)
cmap.set_under("white")

pcol = ax.pcolormesh(
    Xran, Yran, eosDT[..., 2], shading="nearest", cmap=cmap, norm=norm, rasterized=True
)
cax = fig.colorbar(pcol, ticks=[0, 1, 2, 3, 4.5, 6.5, 8])
cax.set_label("")
cax.ax.tick_params(labelsize="x-small")
cax.ax.minorticks_off()
cax.ax.set_yticklabels(
    ["blend", "HELM", "OPAL/SCVH", "FreeEOS", "Skye", "ideal", "WDM_AV"]
)

# save figure
# fig.savefig('eos_regions.pdf')
fig.savefig("eos_regions.png", dpi=300)

# plt.show()
