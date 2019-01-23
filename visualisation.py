import numpy as np
import matplotlib.pyplot as plt


def analyze():
    """
    Plots output of MPI model in part 3 by loading in f.dat, the output file
    from p3.f90.
    """
    C = np.loadtxt('contamination.dat')
    n = len(C)-2
    r = np.linspace(1,1+np.pi,n+2)
    t = np.linspace(0,np.pi,n+2) #theta

    plt.figure(1, figsize=(20,16))
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.contourf(t,r,C,cmap=plt.cm.jet)
    ax.set_rlim(0,5)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    plt.title('Final concentration field')

    return None


if __name__=='__main__':
    analyze()
