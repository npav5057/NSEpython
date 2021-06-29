
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm
import numpy as np
import matplotlib as mpl
from numpy.core.fromnumeric import ptp
from skimage.transform import resize
import os
import sys
from glob import glob


# Enable looped gif by overwriting the PillowWriter class and adding loop=0


if(len(sys.argv)<2):
    print("folder name not Entered")
    print("enter Correct  Folder name as system arguments")
    exit()

class myPillow(anim.PillowWriter):
    def finish(self):
        self._frames[0].save(
            self._outfile, save_all=True, append_images=self._frames[1:],
            duration=int(1000 / self.fps), loop=0)

home=os.getcwd()
path=home+"/OUTPUT/"+sys.argv[1]
if os.path.exists(path):
    os.chdir(path)
    print(path)
else:
    print("Specified folder doesnt exist")
    print("checck the path:",path)
    print("or enter Correct Folder name: as system arguments")
    exit()


files=sorted(glob("*npz"))
print(len(files))
if(len(files)>3000):
    files=files[::10]

# exit()
data = np.load(files[10])
u=data['vx']
v=data['vy']
# v=v[::10]
# u=u[::10]
del data

s = np.sqrt(u**2+ v**2)

u_red = resize(u,(11,11))
v_red = resize(v,(11,11))

print(u.shape)
print(len(files))

X,Y = np.meshgrid(np.linspace(0, 1,u.shape[0]), np.linspace(0, 1, u.shape[1]))


fig = plt.figure(1, [5, 5])
ax = fig.gca()

# Plot colormap.
cmap = mpl.cm.get_cmap('viridis')
norm = BoundaryNorm(np.linspace(0.0, 1.0, 21), cmap.N)
speed = ax.imshow(s, norm = norm, cmap=cmap, origin = "lower", extent = (0, 1, 0, 1))

# Plot colorbar.
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
bar = fig.colorbar(speed, cax=cax, orientation='vertical')
loc = mpl.ticker.MultipleLocator(0.2)
bar.locator = loc
bar.update_ticks()

# Plot vector field.
X, Y = np.mgrid[0:1:11j, 0:1:11j]
vec = ax.quiver(Y, X, u_red, v_red, scale=1.0, color="white")


plt.tight_layout()
plt.xlabel("X Axis")
plt.ylabel("Y Axis")

def animate(n):
    data = np.load(files[n])
    u=data['vx']
    v=data['vy']
    if(n%10==0):
        print("frame::",n)
  
    s = np.sqrt(u**2+v**2)
    speed.set_data(s)
    u_red = resize(u, (11, 11))
    v_red = resize(v, (11, 11))
    vec.set_UVC(u_red, v_red)
    t=n/(len(files)-1)
    ax.set_title("Velocity at time:{}".format(t))
    # plt.savefig("Vel{}.png".format(t))
    return speed, vec
writer = myPillow()
writer.fps = 1
animation = anim.FuncAnimation(fig, animate, frames=len(files), interval = 2, blit=True)
file="_{}_vel.gif".format(sys.argv[1])
animation.save(file, writer=writer)
