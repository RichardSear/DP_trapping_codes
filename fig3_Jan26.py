# %%
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from numpy import pi as π
import matplotlib.pyplot as plt
from scipy.integrate import quad
import final_flowfields
#
print(π)

drift = final_flowfields.drift_pipette_spherical

wallpore=True
wallpore=False
if(wallpore):
#    drift = drift_spherical_wallpore
    drift=hole_drift



k = 1000.0
Γ = 100.0
Ds = 1600.0
Dp = 2.0
R1 = 1.0
α = 0.3
# pore in wall
h2=R1**2



# Q in pl/s
Q = 2.0

Q = 1e3 * Q # convert to um^3/s
#
rcutoff = 1.0#0.0
sqrt2Dp = np.sqrt(2*Dp) # for convenience
print(Dp,sqrt2Dp)
rstar = π*R1/α # where stokeslet and radial outflow match

#
def func_quad(z1):
    r = np.array([0, 0, z1])
    u = drift(r,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#    print(u)
    uz1=u[2]
    return -uz1

def calc_barrier_height():
    r_stable,r_saddle=final_flowfields.calc_fixedpts(Q,k,Γ,Ds,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#    print('stable and unstable fps ',r_stable,r_saddle)
    I = quad(func_quad, r_stable,r_saddle)#,args=(Qbarrier))
    return I

#data_in=np.loadtxt('BDout.txt')
#Qvary_pLs,kλstarvary,rmsvary0,err0,rmsvary1,err1,rmsvary2,err2= \
#    data_in[:,0],data_in[:,1],data_in[:,2],data_in[:,3], \
#        data_in[:,4],data_in[:,5],data_in[:,6],data_in[:,7]
#
Q_BIF=4.0*π*Ds*rstar*(np.sqrt(k*Γ/Ds)-1.0)**2/k
print('Q at BIFurcation ',Q_BIF,' um^3/s')
print('Q at BIFurcation ',Q_BIF*1.0e-3,' pL/s')
n_plotpts=800
Q_fp_plot=np.geomspace(1.0e-5,Q_BIF*0.9999e-3,n_plotpts)
r_fp_m=np.zeros(n_plotpts)
r_fp_p=np.zeros(n_plotpts)
barrier=np.zeros(n_plotpts)
for iQ in range(0,n_plotpts):
    Q=Q_fp_plot[iQ]*1.0e3
    λ = Q / (4*π*Ds) # definition
    λstar = λ / rstar # should be the same as Pe (salt)
    v1 = Q / (π*R1**2) # flow speed (definition)
    Pbyη = α*R1*v1 # from Secchi et al, should also = Q / r*
    Pe = Pbyη / (4*π*Ds) # definition from Secchi et al
#
    r_fp_m[iQ],r_fp_p[iQ]=final_flowfields.calc_fixedpts(Q,k,Γ,Ds,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#    print('fixed points ',r_fp_m[iQ],r_fp_p[iQ])
#
    barrier[iQ],mean_speed=calc_barrier_height()
#    print('barrier ',barrier[iQ])





# Plot the data
num_subplots=2
fig, axs = plt.subplots(num_subplots,1, \
    figsize=(6,3*num_subplots),constrained_layout=True, gridspec_kw={'height_ratios': [2, 1.2]})
# Ensure axes is iterable
if num_subplots == 1:
    axs = [axs]

#
axs[0].plot(Q_fp_plot,r_fp_m,lw=3,c='cyan')#,label='stable fp')
axs[0].plot(Q_fp_plot,r_fp_p,lw=3,c='magenta')#,label='saddle')
axs[0].plot(Q_fp_plot[-1],r_fp_p[-1],lw=0,marker='o',color='k',markersize=10)
#
#
#axs[0].plot(Qvary_pLs,kλstarvary,lw=3,label='kλ*',alpha=0.75,color='green',zorder=-19)
#axs[0].scatter(Qvary_pLs,rmsvary0,s=40,label='$r_c=0.1$',alpha=0.85,lw=4,marker='v',zorder=23,color='darkcyan')
#axs[0].scatter(Qvary_pLs,rmsvary1,s=40,label='$r_c=1$',alpha=0.85,lw=6,marker='2',zorder=23,color='crimson')
#axs[0].scatter(Qvary_pLs,rmsvary2,s=40,label='$r_c=3$',alpha=0.85,zorder=19,color='purple',marker='2')
#axs[0].errorbar(Qvary_pLs,rmsvary0,yerr=err0,fmt='o', capsize=5,label='$D_p=2$',alpha=0.85,lw=3,marker='v',zorder=23,color='blue')
#axs[0].errorbar(Qvary_pLs,rmsvary1,yerr=err1,fmt='o', capsize=5,label='$D_p=20$',alpha=0.85,lw=3,marker='2',zorder=23,color='crimson')
#axs[0].errorbar(Qvary_pLs,rmsvary2,yerr=err2,fmt='o', capsize=5,label='$D_p=50$',alpha=0.85,zorder=19,lw=3,color='purple',marker='2')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlim(1.0e-4, 1.0e2)
axs[0].set_ylim(1.0e0, 5.0e4)
# Add labels and title
axs[0].tick_params(axis='both', which='major', labelsize=18)
axs[0].set_xlabel('$Q$ (pL/s)',fontsize=16)
axs[0].set_ylabel('RMS ($\mathrm{\mu}$m)',fontsize=16)
# Add a legend
#axs[0].legend(fontsize=14,loc='center')
axs[0].legend(fontsize=14,loc='lower left')
axs[0].text(0.5e-5,5.0e4,'(a)',fontsize=18,bbox=dict(facecolor='w', alpha=1.0,edgecolor='w'))
y = np.linspace(1.0e-3, 1.0e5, 100)
x1 = np.ones(100) * 1.5e-2  # Vertical line position
x2 = np.ones(100) * 9  # Vertical line position
axs[0].fill_betweenx(y, x1, x2, color='darkcyan', alpha=0.3)
#
if(num_subplots==2):
    axs[1].loglog(Q_fp_plot,barrier,lw=4,label='$1$',color='red')
#axs[1].loglog(Q_fp_plot*1.0e15,barrier*1.0e11,lw=4,label='$10$',color='blue')
    axs[1].set_ylim([0.1,1000])
    axs[1].set_xlim([1.0e-4,100])
    axs[1].set_xlabel('$Q$ (pL/s)', fontsize=16)
    axs[1].set_ylabel('$\Delta\mathfrak{S}$  ($\mathrm{\mu}$m$^2$/s)', fontsize=16)
    axs[1].tick_params(axis='both', which='major', labelsize=16)
    axs[1].tick_params(axis='both', which='minor', labelsize=16)
    axs[1].text(0.5e-5,1.0e3,'(b)',fontsize=18,bbox=dict(facecolor='w', alpha=1.0,edgecolor='w'))

#
#plt.tight_layout()
plt.savefig('varDp.pdf')
plt.show()
