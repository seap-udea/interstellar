import matplotlib.pyplot as plt
import numpy as np

def plotCloud():
    ###################################################
    # READ DATA
    ###################################################
    data=np.loadtxt("cloud.data")

    ###################################################
    # PREPARE DATA
    ###################################################
    """
    Data map:
    0:i
    1:tdb (terminal)
    2:tdb (future)
    3-8:Position Ecliptic J2000
    9-14:Position J2000
    15-20:Position Galactic J2000
    21:RA(h) (terminal)
    22:DEC(deg)
    23:l(deg)
    24:b(deg)
    25:d(AU)
    26-33:Asymptotic elements, q,e,i,W,w,Mo,to,mu
    34-39:Future Position Ecliptic J2000
    40-45:Future Position Galactic 
    46:RA(h) (future)
    47:DEC(deg)
    48:l(deg)
    49:b(deg)
    50:d(pc)
    """
    RAs=data[:,21]
    DECs=data[:,22]

    ###################################################
    # PLOT CLOUD IN THE SKY
    ###################################################
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(RAs,DECs,'ko')
    fig.savefig("cloud-sky.png")

plotCloud()
