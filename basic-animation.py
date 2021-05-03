import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import copy
from matplotlib import animation
from matplotlib import cm

example_filename = 'data/ds003059/sub-001/ses-LSD/func/sub-001_ses-LSD_task-rest_run-03_bold.nii.gz'

img = nib.load(example_filename)

print(img.shape)
print("Converting to numpy array. This might take some time...")
img_f = img.get_fdata()
print("Data is numpy now.")
#             [x, y, z,  t]
img_f2 = img_f[:, :, 45, :] # Fix the z-level (roughly in the middle)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()

cmap = copy.copy(cm.get_cmap("plasma"))
cmap.set_bad(color='black')
im = plt.imshow(img_f2[:, :, 0], cmap=cmap, interpolation='none', vmin=750, vmax=1400)
t_max = img_f2.shape[-1]

img_f2[img_f2 == 0] = np.nan # hack to make all 0 values black


# Animation Code shamelessly stolen from stackoverflow.com/questions/17212722/
# and https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/


# initialization function: plot the background of each frame
def init():
    im.set_data(img_f2[:, :, 0])
    return [im]


# animation function.  This is called sequentially
def animate(i):
    im.set_data(img_f2[:, :, i+1])
    return [im]


print("Animating...")
# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(
    fig,
    animate,
    #init_func=init,
    frames=(t_max - 1),
    interval=20,
    blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.htmlab
anim.save('sub-001_ses-LSD_task-rest_run-02_bold.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
