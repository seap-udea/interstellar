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
    1:tdb
    2-7:Position Ecliptic J2000
    8-13:Position J2000
    14-19:Position Galactic J2000
    20:RA(h) (terminal)
    21:DEC(deg)
    22:l(deg)
    23:b(deg)
    24:d(AU)
    25-32:Asymptotic elements, q,e,i,W,w,Mo,to,mu
    """

    RAs=data[]

    ###################################################
    # PLOT CLOUD IN THE SKY
    ###################################################
    fig=plt.figure()
    ax=fig.gca()

    ax.plot(RAs,DECs,'ko')

    fig.savefig("cloud-sky.png")
