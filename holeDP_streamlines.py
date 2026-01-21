'''
Python code to plot schematic of hole flow + DP field
Warren and Sear 2025/2026
'''
# %%
# -*- coding: utf-8 -*-
import numpy as np
from numpy import pi as π
import matplotlib.pyplot as plt
import final_flowfields
#
print(π)

drift = final_flowfields.drift_hole



Γ = 150.0
Ds = 1600.0
R1 = 1.0
α = 0.3



rcutoff = 1.0#0.0

def streamline_calc():
    global width
    '''
        first find if there is a fixed point
    '''
    egg=False
    root=0.0
    neggpts=100
    xp_egg=np.zeros(neggpts)
    z_egg=np.zeros(neggpts)
    aquad=k*Γ/Ds-3.0
    bquad=-3.0*k*λhalf
    cquad=k*Γ/Ds
    discquad=bquad**2-4.0*aquad*cquad
    print('discrim ',discquad,'and ',aquad)
    if(discquad > 0.0):
        egg=True
        root_approx=3*k*Q/(2*π*(k*Γ-3*Ds))
        root=(-bquad+np.sqrt(bquad**2-4.0*aquad*cquad))/(2.0*aquad)
        print('fixed point along z at z=',root,root_approx)
    # compute egg
        r_egg=np.linspace(1.0e-3,root*0.9999,neggpts)
        for i in range(0,neggpts):
            z_egg[i]=np.sqrt( (k*Γ/(3.0*Ds)) * r_egg[i]**3/(r_egg[i]+k*λhalf) )
            xp_egg[i]=np.sqrt(r_egg[i]**2-z_egg[i]**2)
#
#    kwargs = {"Q": Q, "k" : k,"Γ" : Γ, "λhalf" : λhalf, "rcutoff" : rcutoff}
#
    width=4#00
    print('width of grid for streamplot in um',width)
# 100 j means 100 points
    x, z = np.mgrid[-width:width:200j, -width*0.01:width:200j]
    print(x.shape)
    ux=np.zeros((200,200))
    uz=np.zeros((200,200))
    r1=np.zeros(3)
    for iz in range(0,200):
        for ix in range(0,200):
            r1[:]=x[iz,ix], 0.0, z[iz,ix]
            ux[iz,ix],dummy,uz[iz,ix]=drift(r1,Q,R1,k,Γ,λhalf,rcutoff)
#
    return egg,z_egg,xp_egg,root,x,z,ux,uz
#
#
n_subplots=2
fig, ax = plt.subplots(1,n_subplots,figsize=(3*n_subplots, 4),sharex=True,sharey=True)
for i in range(0,n_subplots):
## Q in pl/s
    if(i==2):
        Q=1
    else:
        Q = 8.0e-1
    Q = 1e3 * Q # convert to um^3/s
    if(i==0):
        k=0
    else:
        k = 200
    print(i,Q,k)
    v1 = Q / (π*R1**2) # flow speed (definition)
    λhalf = Q / (2*π*Ds) # definition
    print('lambda half ',λhalf)
#
    egg,z_egg,xp_egg,root,x,z,ux,uz=streamline_calc()
    if(egg):
#       ax[i].plot(z_egg,xp_egg,lw=3,label='sep',color='red',zorder=19,ls='dashed')
#       ax[i].plot(z_egg,-xp_egg,lw=3,color='red',zorder=19,ls='dashed')
       ax[i].scatter(root,np.array([0.0]),s=80,color='magenta',lw=3,marker='+',zorder=20)
    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2,density=1.5)

# Add labels and title
#ax.tick_params(axis='both', which='major', labelsize=18)
    ax[i].set_xlabel(r'$z$ ($\mathrm{\mu}$m)',fontsize=18)
    ax[i].set_aspect('equal', adjustable='box')
    ax[i].set_xlim(-width*0.01, width)
    ax[i].set_ylim(-width, width)
    if(i==0): 
        ax[i].set_ylabel(r'$x$ ($\mathrm{\mu}$m)',fontsize=18)
# Add text
        ax[i].text(-0.45*width,0.9*width,'(a)',fontsize=18,bbox=dict(facecolor='w', alpha=1.0,edgecolor='w'))

    else:
        ax[i].text(-0.35*width,0.9*width,'(b)',fontsize=18,bbox=dict(facecolor='w', alpha=1.0,edgecolor='w'))
#
plt.tight_layout(pad=0.5, w_pad=-5.00)
plt.savefig('wallstreamlines.pdf')
plt.show()