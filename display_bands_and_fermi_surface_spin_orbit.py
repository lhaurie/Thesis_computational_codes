import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from numpy import cos, sin,sqrt,tanh, exp
from matplotlib.collections import LineCollection
import scipy.linalg as LA
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib import rc
import os

os.chdir(r"C:\Users\louis\Documents\Recherche\Thèse\Strong coupling COM\codes\codes spin-orbit coupling")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rcParams['figure.figsize'] = [9,15]
nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        #font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 24,
        "font.size": 24,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 24,
        "xtick.labelsize": 24,
        "ytick.labelsize": 24,
}

mpl.rcParams.update(nice_fonts)
mpl.rcParams.update(nice_fonts)
mpl.rc('font', family='serif', serif='cm10')

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r'\boldmath'

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})#plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def get_mu_NI(tbN): 
    not_found=True
    t=1
    if not_found: 
        print("looking for the optimal non interacting chemical potential")
        mu_min=np.zeros((4))
        E=np.zeros((4))
        
        for i in range(1):
            mmin=-2
            mmax=2
            S=0
            while np.abs(S-n)>0.005:
                muNI=(mmax+mmin)/2
                for nkx in range(len(kx_range)):
                    for nky in range(len(ky_range)):                    
                        kx=kx_range[nkx]
                        ky=ky_range[nky]
                        
                        al1=(cos(kx)+cos(ky))/2        
                        E[i]=-4*t*al1
                        if beta*( E[i]  -  muNI)<0:
                            S+= 1/(1+np.exp(beta*( E[i]  -  muNI)))
                            
                        if beta*( E[i]  -  muNI)>0:
                            S+=np.exp(-beta*(E[i]  -  muNI))/(1+np.exp(-beta*( E[i]  -  muNI))) 
                S=2*S/Nk**2  #2 because of spin degeneracy
                if S-n<0:  #S<n
                    mmin=muNI
                if S-n>0:
                    mmax=muNI
            mu_min[i]=muNI
    return(mu_min[0])



def Plot_(mu,dn, px, py,ex, ey, ux, uy, vxy,tx,ty,lambda_,omega,plot):
    #eigenvalues of E 
    eps=np.zeros((Nk,Nk,4))
    #Spectral functions at the fermi surface
    AxFS=np.zeros((Nk,Nk))
    AyFS=np.zeros((Nk,Nk))
    #Qparticle weights
    Qweightx=np.zeros((Nk,Nk,4))
    Qweighty=np.zeros((Nk,Nk,4))
    E=0  #electric field
    
    #orbitally resolved electron density
    nx=(n+dn)/2
    ny=(n-dn)/2

    #elements of I^{-1}
    alphax=2/(2-nx)
    betax=2/nx
    alphay=2/(2-ny)
    betay=2/ny

    for nkx in range(len(kx_range)):
            for nky in range(len(ky_range)):
                    S=np.zeros((4,4),dtype=complex)

                    kx=kx_range[nkx]
                    ky=ky_range[nky]
                    al=(cos(kx)+cos(ky))/2
                    m11x=-mu*(1-nx/2)- 4*tx*ex-4*tx*al*(1-nx+px)-lambda_*ux
                    m12x=4*tx*ex - 4*tx*al*(nx/2-px)+lambda_*ux
                    m22x=-(mu-U)*nx/2-4*tx*ex-4*tx*al*px-lambda_*ux

                    m11y=-(mu+E)*(1-ny/2)- 4*ty*ey-4*ty*al*(1-ny+py)-lambda_*uy
                    m12y=4*ty*ey - 4*ty*al*(ny/2-py)+lambda_*uy
                    m22y=-(mu-U+E)*ny/2-4*ty*ey-4*ty*al*py-lambda_*uy

                    m13=-lambda_*(1-(nx+ny)/2+vxy)
                    m14=-lambda_*(ny/2-vxy)
                    m23=-lambda_*(nx/2-vxy)
                    m24=-lambda_*vxy

                    matM=[[m11x,m12x,m13,m14],
                            [m12x,m22x,m23, m24],
                            [m13,m23,m11y,m12y],
                            [m14,m24,m12y,m22y]]

                    matI=[[1-nx/2,0,0,0],[0,nx/2,0,0],[0,0,1-ny/2,0],[0,0,0,ny/2]]
                    matIinv=[[alphax,0,0,0],[0,betax,0,0],[0,0,alphay,0],[0,0,0,betay]]
                    matE=np.matmul(matM, matIinv)
                    eps[nkx, nky, :], vp = np.linalg.eig(matE)

                    idx = eps[nkx,nky,:].argsort()  
                    eps[nkx,nky,:] = eps[nkx,nky,idx]
                    vp = vp[:,idx]

                    for i in range(4):
                            Qweightx[nkx,nky,i]=np.abs(vp[0,i])**2+np.abs(vp[1,i])**2
                            Qweighty[nkx,nky,i]=np.abs(vp[2,i])**2+np.abs(vp[3,i])**2

                    
                    vpinv=np.matmul(np.linalg.inv(vp), matI)
                    sigma=np.zeros((4,4,4))
                    for i in range(4):
                            for j in range(4):
                                    for k in range(4):
                                            sigma[k,i,j]= (vp[i,k]*vpinv[k,j])


                    etaFS=0.1
                    S[:,:]=sigma[0,:,:]/(-eps[nkx,nky,0]+1j*etaFS) + sigma[1,:,:]/(-eps[nkx,nky,1]+1j*etaFS)+ sigma[2,:,:]/(-eps[nkx,nky,2]+1j*etaFS)+ sigma[3,:,:]/(-eps[nkx,nky,3]+1j*etaFS)

                    AxFS[nkx,nky]=-1/np.pi*np.imag(S[0,0]+S[1,1]+S[0,1]+S[1,0])    
                    AyFS[nkx,nky]=-1/np.pi*np.imag(S[2,2]+S[3,3]+S[2,3]+S[3,2])
    
    N1=int(Nk/2)
    ener1=np.zeros(3*N1)
    ener2=np.zeros(3*N1)
    ener3=np.zeros(3*N1)
    ener4=np.zeros(3*N1)
    Qx=np.zeros((3*N1,4))
    Qy=np.zeros((3*N1,4))
    for i in range(N1):   #plot high symmetry points of A in func for energy
    # (0,0)->(pi,pi)
            ener1[i]=eps[N1-1+i,N1-1+i,0]
            ener2[i]=eps[N1-1+i,N1-1+i,1]
            ener3[i]=eps[N1-1+i,N1-1+i,2]
            ener4[i]=eps[N1-1+i,N1-1+i,3]
            Qx[i,:]=Qweightx[N1-1+i,N1-1+i,:]
            Qy[i,:]=Qweighty[N1-1+i,N1-1+i,:]

    for i in range(N1):
            # (pi,pi)-> (pi,0)
            ener1[i+N1]=eps[2*N1-1,2*N1-1-i,0]
            ener2[i+N1]=eps[2*N1-1,2*N1-1-i,1]
            ener3[i+N1]=eps[2*N1-1,2*N1-1-i,2]
            ener4[i+N1]=eps[2*N1-1,2*N1-1-i,3]
            Qx[i+N1,:]=Qweightx[2*N1-1,2*N1-1-i,:]
            Qy[i+N1,:]=Qweighty[2*N1-1,2*N1-1-i,:]

    for i in range(N1):
    # (pi,0)-> (0,0)
            ener1[i+2*N1]=eps[2*N1-1-i,N1-1,0]
            ener2[i+2*N1]=eps[2*N1-1-i,N1-1,1]
            ener3[i+2*N1]=eps[2*N1-1-i,N1-1,2]
            ener4[i+2*N1]=eps[2*N1-1-i,N1-1,3]
            Qx[i+2*N1,:]=Qweightx[2*N1-1-i,N1-1,:]
            Qy[i+2*N1,:]=Qweighty[2*N1-1-i,N1-1,:]

    
    kk=np.linspace(0,3*N1,3*N1)
    fig=plt.figure()
    fig.suptitle(r"Lower Hubbard bands for $n=$"+str(n)+r" $\lambda =$"+str(lambda_)+r" and $U=$"+str(U))
    ax  = plt.gca()
    x1=[0,N1,2*N1,3*N1-1]
    squad=[0,0,0,0]
    squad[0]=r'$\Gamma$'
    squad[1]=r'M'
    squad[2]='X'
    squad[3]=r'$\Gamma$'
    ax.set_xticks(x1)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    ax.set_ylim(-5,3)
    ax.set_xticklabels(squad, minor=False, rotation=0)
    ax.set_ylabel(r"$E(k)$",fontsize=26)  
    inset_ax = inset_axes(ax,width="20%", height=0.7,loc="lower center")
    inset_ax.set_xlim(0,3*N1)
    inset_ax.set_xticks(x1)
    inset_ax.spines['bottom'].set_linewidth(2)
    inset_ax.spines['top'].set_linewidth(2)
    inset_ax.spines['right'].set_linewidth(2)
    inset_ax.spines['left'].set_linewidth(2)
    inset_ax.yaxis.set_tick_params(width=2)
    inset_ax.xaxis.set_tick_params(width=2)
    inset_ax.set_ylim(U-5,U+3)
    inset_ax.set_xticklabels([], minor=False)


    ener=[ener1, ener2, ener3, ener4]
    cmap_reversed = mpl.cm.get_cmap('winter')  

    for j in range(4):
            norm = mpl.colors.Normalize(vmin=-1.75, vmax=1.75, clip=False)
            points = np.array([kk, ener[j]]).T.reshape(-1,1,2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments,cmap=cmap_reversed, norm=norm)
            dydx=np.zeros(3*N1)
            for i in range(3*N1):
                    dydx[i] = Qx[i,j]-Qy[i,j]
            lc.set_array(dydx)
            lc.set_linewidth(4)
            if j==0 or j==1:
                    line_low = ax.add_collection(lc)
            else:
                    line_up = inset_ax.add_collection(lc)
    cbar=fig.colorbar(lc, ax = ax, location = 'right', ticks=[-1.75,1.75])
    cbar.ax.set_yticklabels([r'$y$',r'$x$'])

    ax.text(2.5,1,r"$\lambda =$"+str(lambda_))
    ax.plot(kk,[0 for i in range(3*N1)],"--",c="black",linewidth=2)
    ax.set_xlim(0,3*N1)

    ###########################FS#######################################
    fig2=plt.figure()
            
    ax = plt.gca()
    fig2.suptitle(r"Fermi surface for $n=$"+str(n)+r" $\lambda =$"+str(lambda_)+r" and $U=$"+str(U))
    ax.plot(np.linspace(0, np.pi,5),np.linspace(0, np.pi,5),c='Black',linewidth=2)
    ax.arrow(np.pi/2, 0, -0.1, 0, shape='full', lw=0, length_includes_head=False, head_width=.5, color="black")
    ax.plot(np.linspace(np.pi, np.pi,5),np.linspace(np.pi, 0,5),c='Black',linewidth=2)
    ax.arrow(np.pi/2, np.pi/2, 0.1, 0.1, shape='full', lw=0, length_includes_head=False, head_width=.5, color="black")
    ax.plot(np.linspace(np.pi, 0,5),np.linspace(0, 0,5),c='Black',linewidth=2)
    ax.arrow(np.pi, np.pi/2, 0 , -0.1, shape='left', lw=0, length_includes_head=False, head_width=.5, color="black")
    ax.text(-0.5,-0.5,r'$\Gamma$',c="Black",fontsize=22)
    ax.text(np.pi-1,np.pi-0.5,r'M',c="Black",fontsize=22)
    ax.text(np.pi-0.4,-0.6,r'X',c="Black",fontsize=22)

    cmap = mpl.cm.get_cmap('copper')
    col = cmap(3/7)
    ax.contour(eps[:,:,0] , colors=[col], linewidths=4,  levels=[0],extent=[-np.pi,np.pi,-np.pi,np.pi],origin="lower")
    ax.contour(eps[:,:,1] , colors=[col], linewidths=4,  levels=[0],extent=[-np.pi,np.pi,-np.pi,np.pi],origin="lower")
    ax.contour(eps[:,:,2] , colors=[col], linewidths=4,  levels=[0],extent=[-np.pi,np.pi,-np.pi,np.pi],origin="lower")
    ax.contour(eps[:,:,3] , colors=[col], linewidths=4,  levels=[0],extent=[-np.pi,np.pi,-np.pi,np.pi],origin="lower")
    ax.set_xticks([-np.pi,0,np.pi])
    ax.set_yticks([-np.pi,0,np.pi])
    ax.set_xticklabels(["-$\pi$" , 0, "$\pi$"])
    ax.set_yticklabels(["-$\pi$" , 0, "$\pi$"])
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    plt.show()
    return(eps)


#Fourier variables
Nk=100
kx_range=np.linspace(-np.pi,np.pi,Nk)
ky_range=np.linspace(-np.pi,np.pi,Nk)
#parameters for the plot
tx=1
ty=1
#number of frequency points
Ne=100
#broadening
eta=0.08
omega=0
beta=1/10**-5
E=0
U=20
n=1.9
lambda_=0.22
xx=[1.1115242025565129       ,-9.9660601512416791E-002   ,9.0860097445313856E-002   ,8.0788997078941327E-002  ,-9.0114289709158801E-002  ,-3.6159984554853889E-004  ,-9.4451615086186217E-003  ,-1.0541428125254829E-002  ,0.21149239193632219]

mu= xx[0]
dn= xx[1]
px= xx[2]
py= xx[3]
ex= xx[4]
ey= xx[5]
ux= xx[6]
uy= xx[7]
vxy=xx[8]

Plot_(mu,dn,px, py,ex, ey, ux, uy, vxy,tx,ty,lambda_,omega,True)

