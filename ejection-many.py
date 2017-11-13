#!/usr/bin/env python
"""

Calculate ejection velocity distribution

Use Wiegert (2014) Formalism

Citation: https://arxiv.org/pdf/1404.2159.pdf

"""
import numpy as np
from scipy.stats import norm as gaussian

#############################################################
#MACROS AND BEHAVIOR
#############################################################
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

TOLMEAN=1e-5
TOLSTD=1e-3
TOLRATIO=1e-3

#############################################################
#INITIAL CONDITIONS
#############################################################
def velocityDistribution(Ms=1.0,ap=3.0,Mp=1e-3,Rp=7e7/UL):
    global G
    #Derived
    mu=G*Ms
    fh=1.0
    RH=fh*ap*(Mp/(3*Ms))**(1./3)
    vp=np.sqrt(mu/ap)
    vesc=np.sqrt(2*mu/ap)
    mup=G*Mp

    #Components of velocity
    vpx=-vp
    vpy=0
    vpz=0
    vp=np.array([vpx,vpy,vpz])

    #############################################################
    #GENERATE RANDOM ELEMENTS
    #############################################################
    Npart=5000
    n=0
    vinfs=[]

    vmean_old=0
    vstd_old=0
    ratio_old=0
    while n<Npart:
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #GENERATE RANDOM ASTROCENTRIC VELOCITIES FOR TEST PARTICLES
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        xdot=-hop*(np.cos(Ob*DEG)*(eb*np.sin(wb*DEG)+np.sin(wpf*DEG))+\
                   np.sin(Ob*DEG)*np.cos(ib*DEG)*(eb*np.cos(wb*DEG)+np.cos(wpf*DEG)))
        ydot=-hop*(np.sin(Ob*DEG)*(eb*np.sin(wb*DEG)+np.cos(wpf*DEG))-\
                   np.cos(Ob*DEG)*np.cos(ib*DEG)*(eb*np.cos(wb*DEG)+np.cos(wpf*DEG)))
        zdot=+hop*np.sin(ib*DEG)*(eb*np.cos(wb*DEG)+\
                                  np.cos(wpf*DEG))

        #PLANETOCENTRIC VELOCITY (ROTATED FOR COMPLAIN WEIGERT 2014)
        xdotrel=(-ydot)-vpx
        ydotrel=(+xdot)-vpy
        zdotrel=(+zdot)-vpz
        vrho=np.sqrt(xdotrel**2+ydotrel**2+zdotrel**2)
        vr=np.sqrt(xdotrel**2+ydotrel**2+zdotrel**2)
        Vi=np.array([xdotrel,ydotrel,zdotrel])

        #PLANETOCENTRIC DIRECTION
        psi,phi=cart2sph(xdotrel,ydotrel,zdotrel)

        #GENERATE RANDOM INCOMING IMPACT PARAMETER
        xtp=uniform(-RH,RH)
        ytp=uniform(-RH,RH)
        beta=np.arcsin(xtp/RH)*RAD;csi=np.arcsin(ytp/RH)*RAD

        #INCOMING POSITION
        Ri=np.array([-RH*np.sin((phi+beta)*DEG)*np.cos((psi+csi)*DEG),
                     +RH*np.cos((phi+beta)*DEG)*np.cos((psi+csi)*DEG),
                     +RH*np.sin((psi+csi)*DEG)])

        #COMPUTE U
        U=np.cross(Ri,Vi)
        u=U/norm(U)    

        #IS THE PLANET APPROACHING?
        qap=np.dot(Ri,Vi)
        if qap>0:continue

        #IMPACT PARAMETER
        B=np.linalg.norm(U)/np.linalg.norm(Vi)
        if B<1.1*Rp:continue

        #GAMMA
        Vin=norm(Vi)
        gamma=2*np.arctan(mup/(B*Vin*Vin))*RAD

        #PLANETOCENTRIC ECCENTRICITY
        e=1/np.sin(gamma*DEG/2)

        #PERICENTER
        q=mup*(e-1)/(Vin*Vin)

        if q<1.1*Rp:continue

        #ROTATION OF INCOMING VELOCITY VECTOR
        c=np.cos(gamma*DEG);s=np.sin(gamma*DEG)
        ux=u[0];uy=u[1];uz=u[2]
        M=np.array([[c+(1-c)*ux**2,(1-c)*uy*ux-s*uz,(1-c)*uz*ux+s*uy],
                    [(1-c)*ux*uy+s*uz,c+(1-c)*uy**2,(1-c)*ux*uy-s*ux],
                    [(1-c)*ux*uz-s*uy,(1-c)*uy*uz+s*ux,c+(1-c)*uz**2]])
        Vf=np.dot(M,Vi)

        #ASTROCENTRIC OUTBOUND VELOCITY
        vf=Vf+vp
        vfn=norm(vf)

        #CHECK IF OBJECT IS BOUND
        if vfn<vesc:continue

        #INFINITE VELOCITY
        vinf=np.sqrt(vfn**2-vesc**2)

        n+=1
        vinfs+=[vinf*UV/1e3]

        vmean=np.mean(vinfs)
        vstd=np.std(vinfs)
        ratio=vstd/vmean

        if np.abs(vmean-vmean_old)/vmean<TOLMEAN and \
           np.abs(vstd-vstd_old)/vstd<TOLSTD and \
           np.abs(ratio-ratio_old)/ratio<TOLRATIO:
            break

        #print(vstd,np.abs(vstd-vstd_old)/vstd)
        #print(ratio,np.abs(ratio-ratio_old)/ratio)

        vmean_old=vmean
        vstd_old=vstd
        ratio_old=ratio
        print(n,vmean,vstd,ratio)

    return n,vmean,vstd,ratio

table=[]
k=1

n,vmean,vstd,ratio=velocityDistribution(Ms=1.0,ap=1.0)
print(n,vmean,vstd,ratio)

"""
aps=np.arange(0.5,5.0,0.1)
Mps=np.logspace(np.log10(3e-4),np.log10(1e-2),10)
ntot=len(aps)*len(Mps)
print("Number of simulations ",ntot)
input()

for ap in aps:
    for Mp in Mps:
        print("Running set %d with parameters:"%k,ap,Mp)
        n,vmean,vstd,ratio=velocityDistribution(Mp=Mp,ap=ap)
        print("\tResults:",n,vmean,vstd,ratio)
        table+=[[ap,Mp,vmean,vstd,ratio]]
        k+=1

table=np.array(table)
np.savetxt("ejection.data",table)
#print(k)
#
#print(n,vmean,vstd,ratio)
"""
