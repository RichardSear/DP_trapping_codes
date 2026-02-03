# %%
# -*- coding: utf-8 -*-
import numpy as np
from numpy import pi as π
import matplotlib.pyplot as plt
import final_flowfields
#
print(π)


drift = final_flowfields.drift_pipette_spherical

def reverse_drift(r,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff):
    minus_drift=(-1.0)*drift(r,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#    print(drift,minus_drift,r)
    return minus_drift


def f2D_drift_reverse(tt,zx, **kwargs):
#    print(kwargs)
#    Q=**kwargs[]
#    Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff = **kwargs
    r1=np.array([zx[1], 0.0, zx[0]])
    vxx,dummy,vzz=reverse_drift(r1,**kwargs)#Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#    print('xx',zx[0],zx[1],vzz,vxx)
    return np.array([vzz,vxx])

def f2D_drift(tt,zx, **kwargs):
#    print(kwargs)
#    print('xx',zx[0],zx[1])
    r1=np.array([zx[1], 0.0, zx[0]])
    vxx,dummy,vzz=drift(r1,**kwargs)#Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
    print('fwd drift xx',zx[0],zx[1],vzz,vxx)
    return np.array([vzz,vxx])

k = 200.0
Γ = 150.0
Ds = 1600.0
R1 = 1.0
α = 0.3



rcutoff = 1.0#0.0

'''
adaptive Euler from MS Copilot, modified to end when withon rcutoff or beyond a cutoff
in z (Y[0])
so NB Y[0] = z coordinate, Y[1] = x coordinate 
'''

def euler_adaptive(f, Y0, t0, z_end, h0, tol=1e-3, **kwargs):
    t, Y, h = t0, np.array(Y0, dtype=float), h0
    ts, Ys = [t], [Y.copy()]

#    while t < t_end:
    while( Y[0]**2+Y[1]**2 > rcutoff**2 and np.abs(Y[0]) < z_end ):
#        if t + h > t_end:
#            h = t_end - t

        # One full step
        Y1 = Y + h * f(t, Y,  **kwargs)

        # Two half steps
        h2 = h / 2
        Y_half = Y + h2 * f(t, Y, **kwargs)
        Y2 = Y_half + h2 * f(t + h2, Y_half, **kwargs)

        # Error estimate
        error = np.linalg.norm(Y2 - Y1)

        if error < tol:
            t += h
            Y = Y2
            ts.append(t)
            Ys.append(Y.copy())
            h *= 1.5
        else:
            h *= 0.5

    return np.array(ts), np.array(Ys)



def traj_calc(reverse,z0,x0):
    '''
    blah
    '''
# initial position on-axis (cartesian method) or off-axis (spherical polars)
# the quadratic for the roots in units of r* is r^2 − (kΓbyD − kλ* − 1)r + kλ* = 0
    b = k*Γ/Ds - k*λstar - 1 # the equation is r^2 − br + c = 0
    Δ = b**2 - 4*k*λstar # discriminant of above

    if Q > 0 and Δ > 0 and b > 0: # condition for roots to exist
        root = np.array((0.5*rstar*(b-np.sqrt(Δ)), 0.5*rstar*(b+np.sqrt(Δ))))
    else:
        root = np.array((np.nan, np.nan))
    print('roots! ',root)
#
#    z0 = R1 if np.isnan(root[0]) else root[0]
#    if(cartesian):
#        r0 = np.array([0, 0, z0])
#    else:
#        r0 = np.array([1e-2, 2e-2, z0])
#    print(r0)
#
    
    kwargs = {"Q": Q, "k" : k,"Γ" : Γ,"rstar" : rstar,"v1" : v1, \
    "Pbyη" : Pbyη,"Pe" : Pe,"λ" : λ,"λstar" : λstar,"rcutoff" : rcutoff}
    if(reverse):
        z0=root[1]
        x0=1.0e-3
        print('ignoring input z0 and x0 and starting near unstable fixed point')
        print('initial position z= ',z0,' x= ',x0,' um')
        t_vals, Y_vals = euler_adaptive(
        f2D_drift_reverse,
        Y0=[z0, x0],
        t0=0.0,
        z_end=2.0e3,
        h0=1.0e3,
        tol=1e-4,
        **kwargs # parameters
        )
    else:
        print('initial position z= ',z0,' x= ',x0,' um')
        t_vals, Y_vals = euler_adaptive(
        f2D_drift,
        Y0=[z0, x0],
        t0=0.0,
        z_end=2.0e3,
        h0=1.0e2,
        tol=1e-4,
        **kwargs # parameters
        )
# now put into arrays for return
    z_traj=Y_vals[:,0]
    x_traj=Y_vals[:,1]
    t_traj=t_vals
    print('z traj ',z_traj)
#    
    return root,z_traj,x_traj,t_traj


#
width=100
print('width of grid for streamplot in um',width)
#
fig, ax = plt.subplots(1,2,figsize=(6, 3.2),sharex=True,sharey=True)

#ax.tick_params(axis='both', which='minor', labelsize=8)
for i in range(0,2):
    if(i==0):
    # Q in pl/s
        Q = 10.0e0
    else:
        Q=1000.0e0
    Q = 1e3 * Q # convert to um^3/s
    print(Q)
#
    rstar = π*R1/α # where stokeslet and radial outflow match
    v1 = Q / (π*R1**2) # flow speed (definition)
    Pbyη = α*R1*v1 # from Secchi et al, should also = Q / r*
    Pe = Pbyη / (4*π*Ds) # definition from Secchi et al
    λ = Q / (4*π*Ds) # definition
    λstar = λ / rstar # should be the same as Pe (salt)
    if(i==0):
        '''
 compute the separatix
        '''
# these values will be ignored....
        z_start=0.0
        x_start=0.0
        fps,z_sep,x_sep,t_sep=traj_calc(reverse=True,z0=z_start,x0=x_start)
        ax[i].plot(z_sep,x_sep,lw=3,label='sep',color='red',zorder=19,ls='dashed')
        ax[i].plot(z_sep,-x_sep,lw=3,color='red',zorder=19,ls='dashed')
        ax[i].scatter(fps[0],np.array([0.0]),s=120,color='magenta',lw=3,marker='+',zorder=99)
        ax[i].scatter(fps[1],np.array([0.0]),s=80,color='red',lw=4,marker='o',zorder=99)
# 100 j means 100 points
    x, z = np.mgrid[-width:width:100j, -width:width:100j]
    print(x.shape)
    ux=np.zeros((100,100))
    uz=np.zeros((100,100))
    r1=np.zeros(3)
    for iz in range(0,100):
        for ix in range(0,100):
            r1[:]=x[iz,ix], 0.0, z[iz,ix]
            ux[iz,ix],dummy,uz[iz,ix]=drift(r1,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#        
    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2,density=1.0)
    ax[i].set_xlim(-width, width)
    ax[i].set_ylim(-width, width)
    # represent pipette
    ax[i].plot([-width,0],[0.0,0.0],lw=4,c='k')
# Add labels and title
    ax[i].set_aspect('equal')
#ax.tick_params(axis='both', which='major', labelsize=18)
    ax[i].set_xlabel(r'$z$ $/\mathrm{\mu}$m',fontsize=16)
    if(i==0): ax[i].set_ylabel(r'$x$ /$\mathrm{\mu}$m',fontsize=16)
#    ax[i].set_ylabel('$x$ ($\mathrm{\mu}$m)',fontsize=18)
    if(i==1): ax[i].set_xticks(ticks=[-50,0,50,100])
#
ax[0].text(-1.75*width,0.9*width,'(a)',fontsize=14,bbox=dict(facecolor='w', alpha=1.0,edgecolor='w'))
ax[1].text(-1.3*width,0.9*width,'(b)',fontsize=14,bbox=dict(facecolor='w', alpha=1.0,edgecolor='w'))
# Set ticks to point inward for both axes
# We change the fontsize of minor ticks label 
ax[0].tick_params(axis='both', which='major', labelsize=16,direction='in')
ax[1].tick_params(axis='both', which='major', labelsize=16,direction='in')
ax[0].tick_params(axis='both', which='minor', labelsize=16,direction='in')
ax[1].tick_params(axis='both', which='minor', labelsize=16,direction='in')
# Add a legend
#ax.legend(fontsize=14,loc='upper left')
plt.tight_layout(pad=0.5, w_pad=-1.00)
plt.savefig('pipette_streamlines.pdf')
plt.show()
