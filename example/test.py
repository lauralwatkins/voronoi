"""
Demo by G. Brammer
"""
import numpy as np
from voronoi import bin2d
import matplotlib.pyplot as plt

# Noisy gaussian
yp, xp = np.indices((100,100))
R = np.sqrt((xp-50)**2+(yp-50)**2)
sigma = 10
g = 10*np.exp(-R**2/2/sigma**2)
s = 1
noise = np.random.normal(size=R.shape)*s

pix_bin, bin_x, bin_y, bin_sn, bin_npix, scale = bin2d.bin2d(xp.flatten(), yp.flatten(), (g+noise).flatten(), g.flatten()*0+s, 20., cvt=True, wvt=False, graphs=False, quiet=False)

# Bin stats
bad = bin_sn < 5
masked = pix_bin*1
mean_bins = pix_bin*0.
median_bins = pix_bin*0.

mea = bin_x*0.
med = bin_x*0.
bx = bin_x*0.
by = bin_y*0.

bin_ids = np.unique(pix_bin)

for i in range(len(bin_ids)):
    bin_mask = pix_bin == bin_ids[i]
    
    mea[i] = (g+noise).flatten()[bin_mask].mean()
    mean_bins[bin_mask] = mea[i]
    
    med[i] = np.median((g+noise).flatten()[bin_mask])
    median_bins[bin_mask] = med[i]
    
    bx[i] = np.sum(xp.flatten()*bin_mask)/bin_mask.sum()
    by[i] = np.sum(yp.flatten()*bin_mask)/bin_mask.sum()

for bin in np.where(bad)[0]:
    bin_mask = pix_bin == bin
    masked[bin_mask] = -99

# Plot
plt.rcParams['image.origin'] = 'lower'

fig = plt.figure(figsize=[9, 2.8])
ax = fig.add_subplot(131)
ax.imshow(pix_bin.reshape(R.shape))
ax.scatter(bin_x, bin_y, marker='.', color='k', alpha=0.1)

ax = fig.add_subplot(132)
ax.imshow(g+noise, vmin=-0.1, vmax=10, cmap='gray_r')

ax = fig.add_subplot(133)
ax.imshow(median_bins.reshape(R.shape), vmin=-0.1, vmax=10, cmap='gray_r')

for ax in fig.axes:
    ax.set_xticklabels([]); ax.set_yticklabels([])
    
fig.tight_layout(pad=0.1)
fig.savefig('test.png')
