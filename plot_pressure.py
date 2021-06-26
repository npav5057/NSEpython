import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm
import numpy as np
import os
import sys
from glob import glob

# Enable looped gif by overwriting the PillowWriter class and adding loop=0
class myPillow(anim.PillowWriter):
    def finish(self):
        self._frames[0].save(
            self._outfile, save_all=True, append_images=self._frames[1:],
            duration=int(1000 / self.fps), loop=0)

if(len(sys.argv)<2):
    print("folder name not Entered")
    print("enter Correct  Folder name as system arguments")
    exit()

home=os.getcwd()
path=home+"/OUTPUT/"+sys.argv[1]

if os.path.exists(path):
    os.chdir(path)
    print(path)
else:
    print("Specified folder doesnt exist")
    print("enter Correct  Folder name: as system arguments")
    exit()


files=sorted(glob("*npz"))
print(len(files))
if(len(files)>3000):
    files=files[::10]

data = np.load(files[0])
p=data['pre']



fig = plt.figure(1, [5, 5])
ax = fig.gca()

# Plot colormap.
cmap = mpl.cm.get_cmap('viridis')
norm = BoundaryNorm(np.linspace(-0.1, 0.1, 21), cmap.N);
speed = ax.imshow(p, norm = norm, origin = "lower", cmap=cmap, extent = (0, 1, 0, 1))

# Plot colorbar.
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
bar = fig.colorbar(speed, cax=cax, orientation='vertical')
bar.locator = mpl.ticker.MultipleLocator(0.02)
bar.update_ticks()
plt.tight_layout()

def animate(n):
    data = np.load(files[n])
    p=data['pre']
    # print(p[1:-1,1:-1])
    speed.set_data(p)
    if(n%10==0):
        print("Frame "  + str(n))
    return speed,
writer = myPillow()
writer.fps = 2
animation = anim.FuncAnimation(fig, animate, frames=len(files), interval = 1, blit=True)
animation.save('Re1_50x50_p.gif', writer=writer)

