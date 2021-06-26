from setup import struct_mesh
from Solver import solverCLS
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import os

if __name__== "__main__":
    # print(" ||.\ || ||==== ||==== ")
    # print(" ||\.\|| !!==!! ||==== ")
    # print(" || \.|| ====|| ||====")
    file="parameters.txt"
    
    model=solverCLS(file)
    os.chdir(model.work)
    # f = open("log.txt", 'w')
    # sys.stdout = f

    


    # data = np.load("cv1000x50x50.npz")
    # vels=data["vel"]
    # prs=data["pre"]
    # print("######################")
    # print(dims)
    init=time.time()
    model.solve2()
    print("time taken :",time.time()-init)

    # file=input("enter file name::")
    # file=file+".npz"
    # np.savez_compressed(file,vx=model.vel_u,vy=model.vel_v,pre=model.sc_p)


    # yy,xx=np.mgrid[0:1:101j, 0:1:101j]
    
    # print(yy.shape,vels[0,:,:].shape)

    # for t in range(len(vels)):
    #     if(t%10!=0):
    #         continue
    #     fig,ax=plt.subplots(1,1)
    #     vel=vels[t,:,:].T
    #     cp = ax.contourf(vel )
    #     fig.colorbar(cp) # Add a colorbar to a plot
    #     plt.xlabel('x_axix')
    #     plt.ylabel('y_axix')
    #     plt.savefig('vel50x50_{}.png'.format(t))
    #     plt.close()
     
    

    # surfaceplots

    # for t in range(len(vels)):
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     ax.view_init(elev=40, azim=20)
    #     vel=vels[t,:,:]
    #     ax.plot_surface(yy,xx,vel,cmap='viridis', edgecolor='none')
    #     ax.set_title('velocity field at time= %f sec\n'%(t*model.dt), size=10,fontsize=16)
    #     # ax.text(0,-1/6,'at time= %f sec\n'%(t*model.dt), size=10,fontsize=16)
    #     plt.xlabel('x_axix')
    #     plt.ylabel('y_axix')
    #     plt.savefig('v_50x50_{}.png'.format(t))
    #     plt.close()
    # f.close()
    # os.chdir(model.home)
       
            





