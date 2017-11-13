"""

Calculate ejection velocity distribution

Use Wiegert (2014) Formalism

Citation: https://arxiv.org/pdf/1404.2159.pdf

"""
import numpy as np
from scipy.stats import norm as gaussian
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

#############################################################
#MACROS AND BEHAVIOR
#############################################################
#np.random.seed(7)
rand=np.random.normal
norm=np.linalg.norm
verbose=0

def uniform(a,b):return a+(b-a)*np.random.rand()
def cart2sph(x,y,z):
    rho=np.sqrt(x**2+y**2)
    r=np.sqrt(rho**2+z**2)
    psi=np.pi/2-np.arccos(z/r)
    phir=np.arctan(y/x)
    if x>0:phi=phir
    elif y > 0:phi=phir+np.pi
    else:phi=phir-np.pi
    return psi*180/np.pi,phi*180/np.pi

#############################################################
#UNITS & CONSTANTS
#############################################################
AU=1.496e11 #m
MSUN=1.98e30 #kg
GCONST=6.67e-11 #m^3/(kg s^2)
RAD=180/np.pi
DEG=1/RAD

G=1.0
UL=1*AU
UM=1*MSUN
UT=np.sqrt(G*UL**3/(GCONST*UM))
UV=UL/UT
if verbose:print("Time unit:",UT)

#############################################################
#INITIAL CONDITIONS
#############################################################
#Stellar mass
Ms=1.0
ap=3.0
Mp=1e-3
Rp=7e7/UL

#Derived
mu=G*Ms
RH=ap*(Mp/(3*Ms))**(1./3)
vp=np.sqrt(mu/ap)
vesc=np.sqrt(2*mu/ap)
mup=G*Mp

#Components of velocity
vpx=-vp
vpy=0
vpz=0
vp=np.array([vpx,vpy,vpz])

if verbose:
    if verbose:print("Planetary radius:",Rp)
    if verbose:print("mu:",mu)
    if verbose:print("Hill radius:",RH)
    if verbose:print("Planetary orbital velocity: ",vp)
    if verbose:print("Planetary system escape velocity:",vesc)
    if verbose:print("mup:",mup)

#############################################################
#GENERATE RANDOM ELEMENTS
#############################################################
Npart=5000
n=0
k=0
vinfs=[]
while n<Npart:
    k+=1
    if verbose:
        print("Test particle:",n)
        #input()

    #GENERATE RANDOM ASTROCENTRIC VELOCITIES FOR TEST PARTICLES
    coswb=2
    while np.abs(coswb)>1:
        ab=rand(ap,ap/2)
        if ab<0:continue
        eb=np.random.rand()
        pb=ab*(1-eb**2)
        ib=uniform(0,90)
        hop=np.sqrt(mu/pb)
        Ob=0.0
        if np.random.rand()>0.5:Ob=180.0
        if Ob==0:
            coswb=(pb-ap)/(ap*eb)
            wpf=0.0
        else:
            coswb=(ap-pb)/(ap*eb)
            wpf=180.0
    wb=np.arccos(coswb)*RAD
    vi=np.sqrt(2*mu/ap-mu/ab)
    if verbose:print("Random elements: ab=%f,eb=%f,ib=%f,Ob=%f,wb=%f,wb+f=%f"%\
          (ab,eb,ib,Ob,wb,wpf))
    if verbose:print("Incoming astrocentric velocity:",vi)
    xdot=-hop*(np.cos(Ob*DEG)*(eb*np.sin(wb*DEG)+np.sin(wpf*DEG))+\
               np.sin(Ob*DEG)*np.cos(ib*DEG)*(eb*np.cos(wb*DEG)+np.cos(wpf*DEG)))
    ydot=-hop*(np.sin(Ob*DEG)*(eb*np.sin(wb*DEG)+np.cos(wpf*DEG))-\
               np.cos(Ob*DEG)*np.cos(ib*DEG)*(eb*np.cos(wb*DEG)+np.cos(wpf*DEG)))
    zdot=+hop*np.sin(ib*DEG)*(eb*np.cos(wb*DEG)+\
                              np.cos(wpf*DEG))
    if verbose:print("\tAstrocentric velocity (Zuluaga RF) :",xdot,ydot,zdot,
                     " (%lf)"%np.sqrt(xdot**2+ydot**2+zdot**2))

    if verbose:print("\tAstrocentric velocity (Wiegert RF) :",-ydot,+xdot,zdot,
                     " (%lf)"%np.sqrt(xdot**2+ydot**2+zdot**2))

    #PLANETOCENTRIC VELOCITY (ROTATED FOR COMPLAIN WEIGERT 2014)
    xdotrel=(-ydot)-vpx
    ydotrel=(+xdot)-vpy
    zdotrel=(+zdot)-vpz
    vrho=np.sqrt(xdotrel**2+ydotrel**2+zdotrel**2)
    vr=np.sqrt(xdotrel**2+ydotrel**2+zdotrel**2)
    if verbose:print("\tPlanetocentric velocity :",xdotrel,ydotrel,zdotrel,
                     " (%lf)"%vr)
    Vi=np.array([xdotrel,ydotrel,zdotrel])
    if verbose:print("\tIncoming velocity: ",Vi)

    #PLANETOCENTRIC DIRECTION
    psi,phi=cart2sph(xdotrel,ydotrel,zdotrel)
    if verbose:print("\tIncoming direction phi=%f,psi=%f"%(phi,psi))
    
    #GENERATE RANDOM INCOMING IMPACT PARAMETER
    fh=0.1
    xtp=uniform(-fh*RH,fh*RH)
    ytp=uniform(-fh*RH,fh*RH)
    beta=np.arcsin(xtp/RH)*RAD;csi=np.arcsin(ytp/RH)*RAD
    if verbose:print("\tIncoming impact parameter beta=%f,csi=%f"%(beta,csi))

    #INCOMING POSITION
    Ri=fh*RH*np.array([-np.sin((phi+beta)*DEG)*np.cos((psi+csi)*DEG),
                       +np.cos((phi+beta)*DEG)*np.cos((psi+csi)*DEG),
                       +np.sin((psi+csi)*DEG)])
    if verbose:print("\tIncoming position: ",Ri)

    #COMPUTE U
    U=np.cross(Ri,Vi)
    if verbose:print("\tIncoming pole: ",U)
    u=U/norm(U)    

    #IS THE PLANET APPROACHING?
    qap=np.dot(Ri,Vi)
    if qap>0:
        if verbose:print("\t\tParticle is receeding")
        continue

    #IMPACT PARAMETER
    B=np.linalg.norm(U)/np.linalg.norm(Vi)
    if verbose:print("\tImpact parameter: ",B)
    if B<1.1*Rp:
        if verbose:print("\t\tObject collided")

    #GAMMA
    Vin=norm(Vi)
    gamma=2*np.arctan(mup/(B*Vin*Vin))*RAD
    if verbose:print("\tGamma: ",gamma)

    #PLANETOCENTRIC ECCENTRICITY
    e=1/np.sin(gamma*DEG/2)
    if verbose:print("\te (planetocentric): ",e)
    
    #PERICENTER
    q=mup*(e-1)/(Vin*Vin)
    if verbose:print("\tPericenter: ",q)

    if q<1.1*Rp:
        if verbose:print("\t\tObject collided")
        continue

    #ROTATION OF INCOMING VELOCITY VECTOR
    c=np.cos(gamma*DEG);s=np.sin(gamma*DEG)
    ux=u[0];uy=u[1];uz=u[2]
    M=np.array([[c+(1-c)*ux**2,(1-c)*uy*ux-s*uz,(1-c)*uz*ux+s*uy],
                [(1-c)*ux*uy+s*uz,c+(1-c)*uy**2,(1-c)*ux*uy-s*ux],
                [(1-c)*ux*uz-s*uy,(1-c)*uy*uz+s*ux,c+(1-c)*uz**2]])
    Vf=np.dot(M,Vi)
    if verbose:print("\tOutbound velocity: ",Vf)

    #ASTROCENTRIC OUTBOUND VELOCITY
    vf=Vf+vp
    vfn=norm(vf)
    if verbose:print("\tOutbound astrocentric velocity: ",vf)

    #CHECK IF OBJECT IS BOUND
    if vfn<vesc:
        if verbose:print("\t\tObject is still bound (vfn = %e, vesc = %e)\n"%(vfn,vesc))
        continue
        
    #INFINITE VELOCITY
    vinf=np.sqrt(vfn**2-vesc**2)
    if verbose:print("\tVelocity at infinite (km/s): ",vinf*UV/1e3)
    
    n+=1
    vinfs+=[vinf*UV/1e3]
    #if verbose:input()
    print(n,vinf*UV/1e3)
    #break

print("Efficiency:",n/(1.*k))
vinfs=np.array(vinfs)
vmean=vinfs.mean()
vstd=vinfs.std()
print("Ratio vstd/vmean :",vstd/vmean)
print("Average velocity (km/s):",vmean)
print("Velocity dispersion:",vstd)

#HISTOGRAM
nbins=10
hs,vs=np.histogram(vinfs,nbins)
vm=(vs[1:]+vs[:-1])/2
hs=np.concatenate(([0],hs,[0]))
vm=np.concatenate(([0],vm,[vs[-1]]))

#PLOT
fig=plt.figure()
ax=fig.gca()
ax.hist(vinfs,nbins)
ax.plot(vm,hs)
ax.plot(vm,20000*gaussian.pdf(vm,vmean,vstd))
fig.savefig("vinfs.png")
