import numpy as np
import os

class solverCLS():
        

    def __init__(self,file):
    
        self.readparams(file) 
        self.dt=0.1
        self.dx=self.lx/self.nx
        self.dy=self.ly/self.ny
        tt=self.Re*min(self.dx,self.dy)**2
       
        self.dt=min((self.Re/((1/self.dx)**2+(1/self.dy)**2)),tt)/2.0
        
        print(self.nx,self.ny)
        print("Chosen dt:",self.dt,end=' ')
        self.p=np.zeros((self.ny+1,self.nx+1),dtype=float)
        self.u=np.zeros((self.ny+1,self.nx+1),dtype=float)
        self.v=np.zeros((self.ny+1,self.nx+1),dtype=float)
        self.tsteps=int(np.ceil(self.t/self.dt))
        print(" timesteps: ",self.tsteps)
        print("Reynolds No:",self.Re)
        if input("want to set custum dt_(y/n):")=="y":
            print("YESS")
            self.dt=float(input())
            self.tsteps=int(np.ceil(self.t/self.dt))
            print("Chosen dt:",self.dt," timesteps: ",self.tsteps)

        self.rt=int(self.tsteps/100)


        self.home=os.getcwd()
        self.out="OUTPUT"
        work=os.path.join(self.home,self.out)
        folder=input("output Folder name: ")
        self.work=os.path.join(work,folder)
        if not os.path.exists(self.work):
            os.mkdir(self.work)
        
        
    
    def setICs(self):
        self.u[:,0]=0
        self.v[:,0]=0
        self.u[:,self.nx]=0
        self.v[:,self.nx]=0
        self.u[0,:]=0
        self.v[0,:]=0
        self.u[self.ny,:]=1.0
        self.v[self.ny,:]=0
        return self.u.copy(),self.v.copy()
    
    def solve(self):
        t=self.tsteps
        nx=self.nx
        ny=self.ny
        dx=self.dx
        dy=self.dy
        dt=self.dt
        rei=1.0/self.Re

        us=np.zeros(self.u.shape,dtype=float)
        vs=np.zeros(self.v.shape,dtype=float)
        bs=np.zeros(self.v.shape,dtype=float)
        ps=self.p.copy()
        # print(self.u.shape,us.shape)
        dy2=dy**2
        dx2=dx**2
        us,vs=self.setICs()
        ct=0

        os.chdir(self.work)
        print(os.getcwd())
        l=len(str(self.tsteps))
        file= "data_"+("0"*l)+".npz"
        np.savez_compressed(file,vx=us,vy=vs,pre=ps)
       
        for t in range(self.tsteps):
           
            us,vs=self.setICs()
            
            du_dx= (us[1:-1,2:]-us[1:-1,:-2])/(2*dx)
            dv_dy= (vs[2:,1:-1]-vs[:-2,1:-1])/(2*dy)
            # bs[1:-1,1:-1]=(du_dx+dv_dy)*(dx2*dy2/dt)
            du_dy= (us[2:,1:-1]-us[:-2,1:-1])/(2*dy)
            dv_dx= (vs[1:-1,2:]-vs[1:-1,:-2])/(2*dx)
            bs[1:-1,1:-1]=(du_dx**2+dv_dy**2+2.0*du_dy*dv_dx+(du_dx+dv_dy)/dt)*(dx2*dy2)

            for it in range(self.maxitr):
                tpi=ps[1:-1,1:-1].copy()
                ps[1:-1,1:-1]=((ps[1:-1,2:]+ps[1:-1,:-2])*dy2+(ps[2:,1:-1]-ps[:-2,1:-1])*dx2-bs[1:-1,1:-1])/(2*(dx2+dy2))
                err=np.sum((ps[1:-1,1:-1]-tpi)**2)
                err=np.sqrt(err/(nx*ny))
                if(err<1e-8):
                    break
                    # print(np.max(tpi),np.min(tpi)   
                    # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))
                ps[0,:]=ps[1,:]
                ps[ny,:]=ps[ny-1,:]
                ps[:,0]=ps[:,1]
                ps[:,nx]=ps[:,nx-1]
              
            # break
            # if(np.max(np.abs(ps[-10:-1,1:7]))>0.0):
            # print(np.max(np.abs(ps)))
            # print(us[-5:-1,-5:-1])
            # print(vs[-5:-1,-5:-1])
            # print(ps[-5:-1,-5:-1])
            
           
            # break
            d2u_dx2=(us[1:-1,2:]+us[1:-1,:-2]-2.0*us[1:-1,1:-1])*(dt/dx2)
            d2u_dy2=(us[2:,1:-1]+us[:-2,1:-1]-2.0*us[1:-1,1:-1])*(dt/dy2)
            udu_dx=us[1:-1,1:-1]*(us[1:-1,1:-1]-us[1:-1,:-2])*(dt/dx)
            # vav=vs[1:-1,1:-1]+vs[1:-1,:-2]+vs[2:,:-2]+vs[2:,1:-1]
            vdu_dy=vs[1:-1,1:-1]*(us[1:-1,1:-1]-us[:-2,1:-1])*(dt/dy)
            dp_dx=(ps[1:-1,:-2]-ps[1:-1,2:])*(dt/(2*dx))
                
            us[1:-1,1:-1]=us[1:-1,1:-1]+(d2u_dx2+d2u_dy2)/self.Re+dp_dx-udu_dx-vdu_dy
            
            d2v_dx2=(vs[1:-1,2:]+vs[1:-1,:-2]-2.0*vs[1:-1,1:-1])*(dt/dx**2)
            d2v_dy2=(vs[2:,1:-1]+vs[:-2,1:-1]-2.0*vs[1:-1,1:-1])*(dt/dy**2)
            vdv_dy=vs[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[:-2,1:-1])*(dt/dy)
            # uav=us[1:-1,1:-1]+us[1:-1,2:]+us[:-2,1:-1]+us[:-2,2:]
            udv_dx=us[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[1:-1,:-2])*(dt/dx)
            dp_dy=(ps[:-2,1:-1]-ps[2:,1:-1])*(dt/(2*dy))
            vs[1:-1,1:-1]=vs[1:-1,1:-1]+(d2v_dx2+d2v_dy2)/self.Re+dp_dy-udv_dx-vdv_dy

            
            self.u[1:-1,1:-1]=us[1:-1,1:-1]
            self.v[1:-1,1:-1]=vs[1:-1,1:-1]
            self.p=ps.copy()
            # # break
            # if(t>1):
            #     break
        
            if np.max(abs(self.u[1:-1,1:-1]))>100.0 or np.max(abs(self.v[1:-1,1:-1]))>100.0:
                print(" #########   Solution Exploded #################")
                break
             
            ct+=1
            if(self.rt==ct):
                print("time:",t*dt,"ITR:",it+1,end=' ')
                print(np.max(abs(self.u[1:-1,1:-1])),np.max(abs(self.v[1:-1,1:-1])),np.max(abs(self.p[1:-1,1:-1])))
                file= "data_"+"0"*(l-len(str(t+1)))+"{}.npz".format(t+1)
                np.savez_compressed(file,vx=self.u,vy=self.v,pre=self.p)
                ct=0
            elif ct%100==0:
                print("time:",t*dt,"ITR:",it+1,end=' ')
                print(np.max(abs(self.u[1:-1,1:-1])),np.max(abs(self.v[1:-1,1:-1])),np.max(abs(self.p[1:-1,1:-1])))

            
            self.u[1:-1,self.nx]= -np.abs(self.u[1:-1,self.nx-1])
            self.u[1:-1,0]=np.abs(self.u[1:-1,1])
            self.v[self.ny,1:-1,]= -np.abs(self.v[self.ny-1,1:-1])
            self.v[0,1:-1]= np.abs(self.v[1,1:-1])
            # break
        if(ct>0):
            file= "data_"+"0"*(l-len(str(t+1)))+"{}.npz".format(t+1)
            np.savez_compressed(file,vx=self.u,vy=self.v,pre=self.p)
        print("Output::",os.getcwd())
        os.chdir(self.home)
        

    def pressure_eq(self):
        # print(self.cvxgrid)
        # print(self.nx,self.ny)
        b=np.zeros((self.ny,self.nx),dtype=float)
        for i in range(1,self.nx):
            for j in range(1,self.ny):
                udx=(self.u[j,i+1]-self.u[j,i-1])*0.5*self.dxi
                vdy=(self.v[j+1,i]-self.v[j-1,i])*0.5*self.dyi
                b[j,i]=(udx+vdy)**2+(udx+vdy)*self.idt
        x2=self.dx**2
        y2=self.dy**2
        t1=0.5/(x2+y2) 
        tp=self.rho*x2*y2/t1
        b=tp*b
        tp=t1
        sz=self.nx*self.ny+1
        for itr in range(100):
            tP=self.pr.copy()
            for i in range(1,self.nx):
                for j in range(1,self.ny):
                    self.pr[j,i]=tp*((self.pr[j+1,i]+self.pr[j-1,i])*x2+(self.pr[j,i+1]+self.pr[j,i-1])*y2)/t1+b[j,i]
            final=np.sum((tP-self.pr)**2)/sz
            if(abs(final)<1e-40):
                # print(itr+1)
                break
        self.pr[0,:]=self.pr[1,:]
        self.pr[:,0]=self.pr[:,1]
        self.pr[self.ny,:]=self.pr[self.ny+1,:]
        self.pr[:,self.nx]=self.pr[:,self.nx+1]
                
    


           

    def next_vel(self):
        tu=self.u.copy()
        tv=self.v.copy()
        xn=self.nx
        yn=self.ny
        self.u[1:yn,1:xn]=tu[1:yn,1:xn]+self.dt*(self.mu*((tu[1:yn,2:xn+1]-2*tu[1:yn,1:xn]+tu[1:yn,0:xn-1])*(self.dxi**2)+(tu[2:yn+1,1:xn]-2*tu[1:yn,1:xn]+tu[0:yn-1,1:xn])*(self.dyi**2)) -0.5*self.dxi*(self.pr[1:yn,2:xn+1]-self.pr[1:yn,0:xn-1])/self.rho - self.dxi*tu[1:yn,1:xn]*(tu[1:yn,1:xn]-tu[1:yn,0:xn-1]) -self.dyi*tv[1:yn,1:xn]*(tu[1:yn,1:xn]-tu[0:yn-1,1:xn]))
        self.v[1:yn,1:xn]=tv[1:yn,1:xn]+self.dt*(self.mu*((tv[1:yn,2:xn+1]-2*tv[1:yn,1:xn]+tv[1:yn,0:xn-1])*(self.dxi)**2 + (tv[2:yn+1,1:xn]-2*tv[1:yn,1:xn]+tv[0:yn-1,1:xn])*(self.dyi)**2)-0.5*self.dyi*(self.pr[2:yn+1,1:xn]-self.pr[0:yn-1,1:xn])/self.rho  - self.dxi*tu[1:yn,1:xn]*(tv[1:yn,1:xn]-tv[1:yn,0:xn-1]) -self.dyi*tv[1:yn,1:xn]*(tv[1:yn,1:xn]-tv[0:yn-1,1:xn]))
        self.u[0,:]=1

        


    def readparams(self,file):
        print("Initialising Parameters:")
        ff= open(file,'r')
        
        line=ff.readline()
        self.nx=int(line.split()[0])

        line=ff.readline()
        self.ny=int(line.split()[0])

        line=ff.readline()
        self.lx=float(line.split()[0])

        line=ff.readline()
        self.ly=float(line.split()[0])

        line=ff.readline()
        self.t=float(line.split()[0])

        line=ff.readline()
        self.Re=float(line.split()[0])

        line=ff.readline()
        self.gx=float(line.split()[0])

        line=ff.readline()
        self.gy=float(line.split()[0])

        line=ff.readline()
        self.maxitr=int(line.split()[0])

        line=ff.readline()
        self.tau=float(line.split()[0])

        line=ff.readline()
        self.omega=float(line.split()[0])

        line=ff.readline()
        self.eps=float(line.split()[0])

        ff.close()
        print(" ")
    

    def solve2(self):
        print("############   SOLVER2  #################")
        t=self.tsteps
        nx=self.nx
        ny=self.ny
        dx=self.dx
        dy=self.dy
        dt=self.dt

        us=np.zeros(self.u.shape,dtype=float)
        vs=np.zeros(self.v.shape,dtype=float)
        bs=np.zeros(self.v.shape,dtype=float)
        ps=self.p.copy()
        # print(self.u.shape,us.shape)
        dy2=dy**2
        dx2=dx**2
        us,vs=self.setICs()
        ct=0

        os.chdir(self.work)
        print(os.getcwd())
        l=len(str(self.tsteps))
        file= "data_"+("0"*l)+".npz"
        np.savez_compressed(file,vx=us,vy=vs,pre=ps)
       
        for t in range(self.tsteps):
           
            us,vs=self.setICs()
            d2u_dx2=(us[1:-1,2:]+us[1:-1,:-2]-2.0*us[1:-1,1:-1])*(dt/(self.Re*dx2))
            d2u_dy2=(us[2:,1:-1]+us[:-2,1:-1]-2.0*us[1:-1,1:-1])*(dt/(self.Re*dy2))
            udu_dx=us[1:-1,1:-1]*(us[1:-1,1:-1]-us[1:-1,:-2])*(dt/dx)
            # vav=vs[1:-1,1:-1]+vs[1:-1,:-2]+vs[2:,:-2]+vs[2:,1:-1]
            vdu_dy=vs[1:-1,1:-1]*(us[1:-1,1:-1]-us[:-2,1:-1])*(dt/dy)
            # dp_dx=(ps[1:-1,:-2]-ps[1:-1,2:])*(dt/(2*dx))
                
            us[1:-1,1:-1]=us[1:-1,1:-1]+(d2u_dx2+d2u_dy2)-udu_dx-vdu_dy
            
            d2v_dx2=(vs[1:-1,2:]+vs[1:-1,:-2]-2.0*vs[1:-1,1:-1])*(dt/(self.Re*dx2))
            d2v_dy2=(vs[2:,1:-1]+vs[:-2,1:-1]-2.0*vs[1:-1,1:-1])*(dt/(self.Re*dy2))
            vdv_dy=vs[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[:-2,1:-1])*(dt/dy)
            # uav=us[1:-1,1:-1]+us[1:-1,2:]+us[:-2,1:-1]+us[:-2,2:]
            udv_dx=us[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[1:-1,:-2])*(dt/dx)
            # dp_dy=(ps[:-2,1:-1]-ps[2:,1:-1])*(dt/(2*dy))
            vs[1:-1,1:-1]=vs[1:-1,1:-1]+(d2v_dx2+d2v_dy2)-udv_dx-vdv_dy

            
            du_dx= (us[1:-1,2:]-us[1:-1,1:-1])/(dx)
            dv_dy= (vs[2:,1:-1]-vs[1:-1,1:-1])/(dy)
            bs[1:-1,1:-1]=(du_dx+dv_dy)*(dx2*dy2/dt)
            # du_dy= (us[2:,1:-1]-us[:-2,1:-1])/(2*dy)
            # dv_dx= (vs[1:-1,2:]-vs[1:-1,:-2])/(2*dx)
            # bs[1:-1,1:-1]=(du_dx**2+dv_dy**2+2.0*du_dy*dv_dx+(du_dx+dv_dy)/dt)*(dx2*dy2)

            for it in range(self.maxitr):
                tpi=ps[1:-1,1:-1].copy()
                ps[1:-1,1:-1]=((ps[1:-1,2:]+ps[1:-1,:-2])*dy2+(ps[2:,1:-1]-ps[:-2,1:-1])*dx2-bs[1:-1,1:-1])/(2*(dx2+dy2))
                err=np.sum((ps[1:-1,1:-1]-tpi)**2)
                err=np.sqrt(err/(nx*ny))
                if(err<1e-8):
                    break
                    # print(np.max(tpi),np.min(tpi)   
                    # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))
                ps[0,:]=ps[1,:]
                ps[ny,:]=ps[ny-1,:]
                ps[:,0]=ps[:,1]
                ps[:,nx]=ps[:,nx-1]
              
            # break
            # if(np.max(np.abs(ps[-10:-1,1:7]))>0.0):
            # print(np.max(np.abs(ps)))
            # print(us[-5:-1,-5:-1])
            # print(vs[-5:-1,-5:-1])
            # print(ps[-5:-1,-5:-1])
          
            dp_dx=(ps[1:-1,:-2]-ps[1:-1,2:])*(dt/(2*dx))
            us[1:-1,1:-1]=us[1:-1,1:-1]+dp_dx
           
            dp_dy=(ps[:-2,1:-1]-ps[2:,1:-1])*(dt/(2*dy))
            vs[1:-1,1:-1]=vs[1:-1,1:-1]+dp_dy

            self.u[1:-1,1:-1]=us[1:-1,1:-1]
            self.v[1:-1,1:-1]=vs[1:-1,1:-1]
            self.p=ps.copy()
            # # break
            # if(t>1):
            #     break
        
            if np.max(abs(self.u[1:-1,1:-1]))>100.0 or np.max(abs(self.v[1:-1,1:-1]))>100.0:
                print(" #########   Solution Exploded #################")
                break
             
            ct+=1
            if(self.rt==ct):
                print("time:",t*dt,"ITR:",it+1,end=' ')
                print(np.max(abs(self.u[1:-1,1:-1])),np.max(abs(self.v[1:-1,1:-1])),np.max(abs(self.p[1:-1,1:-1])))
                file= "data_"+"0"*(l-len(str(t+1)))+"{}.npz".format(t+1)
                np.savez_compressed(file,vx=self.u,vy=self.v,pre=self.p)
                ct=0
            elif ct%100==0:
                print("time:",t*dt,"ITR:",it+1,end=' ')
                print(np.max(abs(self.u[1:-1,1:-1])),np.max(abs(self.v[1:-1,1:-1])),np.max(abs(self.p[1:-1,1:-1])))

            
            self.u[1:-1,self.nx]= -np.abs(self.u[1:-1,self.nx-1])
            self.u[1:-1,0]=np.abs(self.u[1:-1,1])
            self.v[self.ny,1:-1,]= -np.abs(self.v[self.ny-1,1:-1])
            self.v[0,1:-1]= np.abs(self.v[1,1:-1])
            # break
        if(ct>0):
            file= "data_"+"0"*(l-len(str(t+1)))+"{}.npz".format(t+1)
            np.savez_compressed(file,vx=self.u,vy=self.v,pre=self.p)
        print("Output::",os.getcwd())
        os.chdir(self.home)

    # def temporalStep(self):
    #     B=np.zeros((self.nx*self.ny))
    #     tn=0
    #     tp=self.rho/self.dt
    #     for j in range(1,len(self.ygrid)):
    #         for i in range(1,len(self.xgrid)):
    #             B[tn]=tp*((self.us[i+1,j]-self.us[i,j])*self.dxi+(self.vs[i,j+1]-self.vs[i,j])*self.dyi)
    #             tn+=1
    #     ptp=np.linalg.solve(self.lap2,B)
    #     P=np.zeros((self.nx+2,self.ny+2))
    #     print(P.shape)
    #     P[1:-1,1:-1]=ptp.reshape(self.nx,self.ny)
    #     return P
    # def setVelocity(self,ui,vi,xcv,ycv):
    #     u=ui*np.ones((ycv,xcv))
    #     v=vi*np.ones((ycv,xcv))
    #     u[0,:]=1
    #     # w=vel[2]*np.ones((ycv,xcv))
    #     return u,v
    
    # def correctorP(self):
    #     tp=self.dt/self.rho
    #     # ssss = time.time()
    #     self.u[2:-1,1:-1]= self.us[2:-1,1:-1]-tp*self.dxi*(self.ssp[2:-1,1:-1]-self.ssp[1:-2,1:-1])
    #     self.v[1:-1,2:-1]= self.vs[1:-1,2:-1]-tp*self.dyi*(self.ssp[1:-1,2:-1]-self.ssp[1:-1,1:-2])
    #     # print("time:",time.time()-ssss)
    #     # print(self.v[1:-1,1:-1])
    #     # ssss = time.time()
    #     # for j in range(1,self.ny+1):
    #     #     for i in range(2,self.nx+1):
    #     #         self.u[i,j]=self.us[i,j]-tp*self.dxi*(self.ssp[i,j]-self.ssp[i-1,j])
    #     # for j in range(2,self.ny+1):
    #     #     for i in range(1,self.nx+1):
    #     #         self.v[i,j]=self.vs[i,j]-tp*self.dyi*(self.ssp[i,j]-self.ssp[i,j-1])
    #     # print("time:",time.time()-ssss)
    #     # print(self.v[1:-1,1:-1])

    # def correctorStepVel(self):
    #     tu=np.zeros(self.u.shape)
    #     for j in range(1,self.ny+1):
    #         for i in range(2,self.nx+1):
    #             tp=0.25*(self.v[i,j]+self.v[i,j-1]+self.v[i+1,j]+self.v[i+1,j-1])
    #             d2ux = (self.u[i-1,j]-2*self.u[i ,j]+self.u[i+1,j])*self.dxi**2
    #             d2uy = (self.u[i ,j-1]-2*self.u[i ,j]+self.u[i ,j+1])*self.dyi**2 
    #             udux=self.u[i,j]*(self.u[i+1,j]-self.u[i-1,j])*0.5*self.dxi 
    #             uduy =tp*(self.u[i,j+1]-self.u[i,j-1])*0.5*self.dyi
    #             tu[i,j]=self.u[i,j]+self.dt*(self.mu*(d2ux+d2uy) -udux -uduy)
    #     tv=np.zeros(self.v.shape)
    #     for j in range(2,self.ny+1):
    #         for i in range(1,self.nx+1):
    #             tp=0.25*(self.u[i,j]+self.u[i,j-1]+self.u[i+1,j]+self.u[i+1,j-1])
    #             d2vx = (self.v[i-1,j]-2*self.v[i,j]+self.v[i+1,j])*self.dxi**2 
    #             d2vy= (self.v[i ,j-1]-2*self.v[i ,j]+self.v[i ,j+1])*self.dyi**2 
    #             udvy =(self.v[i,j+1]-self.v[i,j-1])*0.5*self.dyi*self.v[i,j]
    #             udvx = (self.v[i+1,j]-self.v[i-1,j])*0.5*self.dxi*tp
    #             tv[i,j]=self.v[i,j]+self.dt*(self.mu*(d2vy+d2vx)-udvy-udvx)
    #     return tu,tv
    # def laplace(self):
    #     lp=np.zeros((self.nx*self.ny,self.nx*self.ny))
    #     ti=0
    #     for j in range(self.ny):
    #         for i in range(self.nx):
    #             lp[i+ti,i+ti]=2*(self.dxi**2+self.dyi**2)
    #             try:
    #                 for ii in range(i-1,i+2,2):
    #                     if ii>=0 and ii<self.nx:
    #                         lp[i+ti,ii+ti]=self.dxi**2
    #                     else:
    #                         lp[i+ti,i+ti]-=self.dxi**2
    #                 for jj in range(j-1,j+2,2):
    #                     if jj>=0 and jj<self.ny:
    #                         lp[i+ti,i+jj*self.nx]=self.dyi**2
    #                     else:
    #                         lp[i+ti,i+ti]-=self.dyi**2
    #             except:
    #                 print("ERRRR ",ii,ti)
    #         ti+=self.nx
    #     lp[0,:]=0
    #     lp[0,0]=1
    #     return lp



            # d2u_dx2=(self.u[1:-1,2:]+self.u[1:-1,:-2]-2.0*self.u[1:-1,1:-1])*(dt/dx2)
            # d2u_dy2=(self.u[2:,1:-1]+self.u[:-2,1:-1]-2.0*self.u[1:-1,1:-1])*(dt/dy2)
            # udu_dx=self.u[1:-1,1:-1]*(self.u[1:-1,1:-1]-self.u[1:-1,-2])*(dt/dx)
            # vav=self.v[1:-1,1:-1]+self.v[1:-1,:-2]+self.v[2:,:-2]+self.v[2:,1:-1]
            # vdu_dy=vav*(self.u[1:-1,1:-1]-self.u[:-2,1:-1,])*0.25*(dt/dy)
            
            # us[1:-1,1:-1]=self.u[1:-1,1:-1]+rei*(d2u_dx2+d2u_dy2)-udu_dx-vdu_dy
            
            # d2v_dx2=(self.v[1:-1,2:]+self.v[1:-1,:-2]-2.0*self.v[1:-1,1:-1])*(dt/dx2)
            # d2v_dy2=(self.v[2:,1:-1]+self.v[:-2,1:-1]-2.0*self.v[1:-1,1:-1])*(dt/dy2)
            # vdv_dy=self.v[1:-1,1:-1]*(self.v[1:-1,1:-1]-self.v[:-2,1:-1])*(dt/dy)
            # uav=self.u[1:-1,1:-1]+self.u[1:-1,2:]+self.u[:-2,1:-1]+self.u[:-2,2:]
            # udv_dx=uav*(self.v[1:-1,1:-1]-self.v[1:-1,:-2])*0.25*(dt/dx)
            # vs[1:-1,1:-1]=self.v[1:-1,1:-1]+rei*(d2v_dx2+d2v_dy2)-udv_dx-vdv_dy

            ##PPE Solver
            # du_dx= (self.u[1:-1,2:]-self.u[1:-1,:-2])*(0.5/dx)
            # dv_dy= (self.v[2:,1:-1]-self.v[:-2,1:-1])*(0.5/dy)
            # du_dy= (self.u[2:,1:-1]-self.u[:-2,1:-1])*(0.5/dy)
            # dv_dx= (self.v[1:-1,2:]-self.v[1:-1,:-2])*(0.5/dx)
            # bs[1:-1,1:-1]=(du_dx**2+dv_dy**2+2.0*du_dy*dv_dx)
            # # bs[1:-1,1:-1]=(du_dx+dv_dy)**2
            # for it in range(self.maxitr):
            #     tt=ps[1:-1,1:-1]
            #     ps[1:-1,1:-1]=((ps[1:-1,2:]+ps[1:-1,:-2])*dy2+(ps[2:,1:-1]-ps[:-2,1:-1])*dx2+bs[1:-1,1:-1]*dx2*dy2)*0.5/(dx2+dy2)
            #     err=np.sum((ps[1:-1,1:-1]-tt)**2)/(nx*ny)
            #     err=np.sqrt(err)
            #     if(err<1e-10):
            #         print("ITR::",it+1," Error:",err)
            #         # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))
            #         break
            #     ps[0,:]=ps[1,:]
            #     ps[ny,:]=ps[ny-1,:]
            #     ps[:,0]=ps[:,1]
            #     ps[:,nx]=ps[:,nx-1]
            
            # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))
            
            # for it in range(self.maxitr):
            #     usi=us.copy()
            #     vsi=vs.copy()
            #     d2u_dx2=(us[1:-1,2:]+us[1:-1,:-2]-2.0*us[1:-1,1:-1])*(dt/dx2)
            #     d2u_dy2=(us[2:,1:-1]+us[:-2,1:-1]-2.0*us[1:-1,1:-1])*(dt/dy2)
            #     udu_dx=us[1:-1,1:-1]*(us[1:-1,1:-1]-us[1:-1,-2])*(dt/dx)
            #     # vav=vs[1:-1,1:-1]+vs[1:-1,:-2]+vs[2:,:-2]+vs[2:,1:-1]
            #     vdu_dy=vs[1:-1,1:-1]*(us[1:-1,1:-1]-us[:-2,1:-1,])*(dt/dy)
            #     dp_dx=(ps[1:-1,:-2]-ps[1:-1,2:])*0.5*(dt/dx)
                
            #     us[1:-1,1:-1]=us[1:-1,1:-1]+(d2u_dx2+d2u_dy2)/self.Re+dp_dx-udu_dx-vdu_dy
            
            #     d2v_dx2=(vs[1:-1,2:]+vs[1:-1,:-2]-2.0*vs[1:-1,1:-1])*(dt/dx**2)
            #     d2v_dy2=(vs[2:,1:-1]+vs[:-2,1:-1]-2.0*vs[1:-1,1:-1])*(dt/dy**2)
            #     vdv_dy=vs[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[:-2,1:-1])*(dt/dy)
            #     # uav=us[1:-1,1:-1]+us[1:-1,2:]+us[:-2,1:-1]+us[:-2,2:]
            #     udv_dx=us[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[1:-1,:-2])*(dt/dx)
            #     dp_dy=(ps[:-2,1:-1]-ps[2:,1:-1])*0.5*(dt/dy)
            #     vs[1:-1,1:-1]=vs[1:-1,1:-1]+(d2v_dx2+d2v_dy2)/self.Re+dp_dy-udv_dx-vdv_dy
                
            #     err=np.sum((vs-vsi)**2)+np.sum((us-usi)**2)
            #     err/=(nx*ny)
            #     err=np.sqrt(err)
            #     if(err<1e-8):
            #         # print(np.max(us[1:-1,1:-1]),np.min(us[1:-1,1:-1]))
            #         # print(np.max(vs[1:-1,1:-1]),np.min(vs[1:-1,1:-1]))
            #         break
            # print("ITR::",it+1," Error:",err)
    
            # du_dx=(us[1:-1,2:]-us[1:-1,:-2])/dx
            # dv_dy=(vs[2:,1:-1]-vs[:-2,1:-1])/dy
            # bs[1:-1,1:-1]=(du_dx+dv_dy)*((dx2*dy2)/(2*dt))
            # for it in range(self.maxitr):
            #     tpi=ps[1:-1,1:-1].copy()
            #     ps[1:-1,1:-1]=((ps[1:-1,2:]+ps[1:-1,:-2])*dy2+(ps[2:,1:-1]-ps[:-2,1:-1])*dx2-bs[1:-1,1:-1])*0.5/(dx2+dy2)
            #     err=np.sum((ps[1:-1,1:-1]-tpi)**2)
            #     err=err/(nx*ny)
            #     err=np.sqrt(err)
            #     # break
            #     if(err<1e-8):
            #         print("ITR:",it+1," Error:",err)
            #         # print(np.max(tpi),np.min(tpi)   
            #         # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))
            #         break
            #     ps[0,:]=ps[1,:]
            #     ps[ny,:]=ps[ny-1,:]
            #     ps[:,0]=ps[:,1]
            #     ps[:,nx]=ps[:,nx-1]
            # print("ITR:",it+1," Error:",err)
            
            # dp_dx=(ps[1:-1,:-2]-ps[1:-1,2:])*0.5*(dt/dx)
            # dp_dy=(ps[:-2,1:-1]-ps[2:,1:-1])*0.5*(dt/dy)
            # self.u[1:-1,1:-1]=us[1:-1,1:-1]+dp_dx
            # self.v[1:-1,1:-1]=vs[1:-1,1:-1]+dp_dx
            # print(np.max(self.u[1:-1,1:-1]),np.min(self.u[1:-1,1:-1]))
            # self.p=ps.copy()
            # print(np.max(self.v[1:-1,1:-1]),np.min(self.v[1:-1,1:-1]))
            # print(np.max(self.p[1:-1,1:-1]),np.min(self.p[1:-1,1:-1]))
            # print(" ")


            ##PPE Solver
            # du_dx= (self.u[1:-1,2:]-self.u[1:-1,:-2])*(0.5/dx)
            # dv_dy= (self.v[2:,1:-1]-self.v[:-2,1:-1])*(0.5/dy)
            # du_dy= (self.u[2:,1:-1]-self.u[:-2,1:-1])*(0.5/dy)
            # dv_dx= (self.v[1:-1,2:]-self.v[1:-1,:-2])*(0.5/dx)
            # bs[1:-1,1:-1]=(du_dx**2+dv_dy**2+2.0*du_dy*dv_dx)
            # # bs[1:-1,1:-1]=(du_dx+dv_dy)**2
            # for it in range(self.maxitr):
            #     tt=ps[1:-1,1:-1]
            #     ps[1:-1,1:-1]=((ps[1:-1,2:]+ps[1:-1,:-2])*dy2+(ps[2:,1:-1]-ps[:-2,1:-1])*dx2+bs[1:-1,1:-1]*dx2*dy2)*0.5/(dx2+dy2)
            #     err=np.sum((ps[1:-1,1:-1]-tt)**2)/(nx*ny)
            #     err=np.sqrt(err)
            #     if(err<1e-10):
            #         print("ITR::",it+1," Error:",err)
            #         # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))
            #         break
            #     ps[0,:]=ps[1,:]
            #     ps[ny,:]=ps[ny-1,:]
            #     ps[:,0]=ps[:,1]
            #     ps[:,nx]=ps[:,nx-1]
            
            # print(np.max(ps[1:-1,1:-1]),np.min(ps[1:-1,1:-1]))

            # for it in range(self.maxitr):
            #     usi=us.copy()
            #     vsi=vs.copy()
            #     d2u_dx2=(self.u[1:-1,2:]+self.u[1:-1,:-2]-2.0*self.u[1:-1,1:-1])*(dt/dx2)
            #     d2u_dy2=(self.u[2:,1:-1]+self.u[:-2,1:-1]-2.0*self.u[1:-1,1:-1])*(dt/dy2)
            #     udu_dx=self.u[1:-1,1:-1]*(self.u[1:-1,1:-1]-self.u[1:-1,-2])*(dt/dx)
            #     vav=self.v[1:-1,1:-1]+self.v[1:-1,:-2]+self.v[2:,:-2]+self.v[2:,1:-1]
            #     vdu_dy=vav*(self.u[1:-1,1:-1]-self.u[:-2,1:-1,])*0.25*(dt/dy)
            
            #     us[1:-1,1:-1]=self.u[1:-1,1:-1]+rei*(d2u_dx2+d2u_dy2)-udu_dx-vdu_dy
            
            #     d2v_dx2=(self.v[1:-1,2:]+self.v[1:-1,:-2]-2.0*self.v[1:-1,1:-1])*(dt/dx2)
            #     d2v_dy2=(self.v[2:,1:-1]+self.v[:-2,1:-1]-2.0*self.v[1:-1,1:-1])*(dt/dy2)
            #     vdv_dy=self.v[1:-1,1:-1]*(self.v[1:-1,1:-1]-self.v[:-2,1:-1])*(dt/dy)
            #     uav=self.u[1:-1,1:-1]+self.u[1:-1,2:]+self.u[:-2,1:-1]+self.u[:-2,2:]
            #     udv_dx=uav*(self.v[1:-1,1:-1]-self.v[1:-1,:-2])*0.25*(dt/dx)
            #     vs[1:-1,1:-1]=self.v[1:-1,1:-1]+rei*(d2v_dx2+d2v_dy2)-udv_dx-vdv_dy

            #     err=np.sum((vs-vsi)**2)+np.sum((us-usi)**2)
            #     err/=(nx*ny)
            #     err=np.sqrt(err)
               
            #     if(err<1e-8):
            #         # print(np.max(us[1:-1,1:-1]),np.min(us[1:-1,1:-1]))
            #         # print(np.max(vs[1:-1,1:-1]),np.min(vs[1:-1,1:-1]))
            #         break
            # print("ITR::",it+1," Error:",err)

            
            # for it in range(self.maxitr):
            #     usi=us.copy()
            #     vsi=vs.copy()
            #     d2u_dx2=(us[1:-1,2:]+us[1:-1,:-2]-2.0*us[1:-1,1:-1])*(dt/dx2)
            #     d2u_dy2=(us[2:,1:-1]+us[:-2,1:-1]-2.0*us[1:-1,1:-1])*(dt/dy2)
            #     udu_dx=us[1:-1,1:-1]*(us[1:-1,1:-1]-us[1:-1,-2])*(dt/dx)
            #     # vav=vs[1:-1,1:-1]+vs[1:-1,:-2]+vs[2:,:-2]+vs[2:,1:-1]
            #     vdu_dy=vs[1:-1,1:-1]*(us[1:-1,1:-1]-us[:-2,1:-1,])*(dt/dy)
            #     dp_dx=(ps[1:-1,:-2]-ps[1:-1,2:])*0.5*(dt/dx)
                
            #     us[1:-1,1:-1]=us[1:-1,1:-1]+(d2u_dx2+d2u_dy2)/self.Re+dp_dx-udu_dx-vdu_dy
            
            #     d2v_dx2=(vs[1:-1,2:]+vs[1:-1,:-2]-2.0*vs[1:-1,1:-1])*(dt/dx**2)
            #     d2v_dy2=(vs[2:,1:-1]+vs[:-2,1:-1]-2.0*vs[1:-1,1:-1])*(dt/dy**2)
            #     vdv_dy=vs[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[:-2,1:-1])*(dt/dy)
            #     # uav=us[1:-1,1:-1]+us[1:-1,2:]+us[:-2,1:-1]+us[:-2,2:]
            #     udv_dx=us[1:-1,1:-1]*(vs[1:-1,1:-1]-vs[1:-1,:-2])*(dt/dx)
            #     dp_dy=(ps[:-2,1:-1]-ps[2:,1:-1])*0.5*(dt/dy)
            #     vs[1:-1,1:-1]=vs[1:-1,1:-1]+(d2v_dx2+d2v_dy2)/self.Re+dp_dy-udv_dx-vdv_dy
                
            #     err=np.sum((vs-vsi)**2)+np.sum((us-usi)**2)
            #     err/=(nx*ny)
            #     err=np.sqrt(err)
            #     if(err<1e-8):
            #         # print(np.max(us[1:-1,1:-1]),np.min(us[1:-1,1:-1]))
            #         # print(np.max(vs[1:-1,1:-1]),np.min(vs[1:-1,1:-1]))
            #         break
            # print("ITR::",it+1," Error:",err)

