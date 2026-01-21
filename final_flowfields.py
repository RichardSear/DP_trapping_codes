import numpy as np
from numpy import pi as π

'''
Expressions for drift field in LS jet and hole in disc problems
Warren and Sear 2025/2026
'''
def drift_pipette_spherical(r,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff):
    x, y, z = r[:] # z is normal distance = r cosθ
    xy = np.sqrt(x**2 + y**2) # in-plane distance = r sinθ
    rr = np.sqrt(x**2 + y**2 + z**2) # radial distance from origin
    cosθ, sinθ = z/rr, xy/rr # polar angle, as cos and sin
    if rr < rcutoff: # crude regularisation
        u = np.zeros_like(r)
    else:
        ur = Q/(4*π*rr**2) + Pbyη*cosθ/(4*π*rr) - k*λ*Γ/(rr*(rr+k*λ))
#        ur = 0.0/(4*π*rr**2) + Pbyη*cosθ/(4*π*rr) - k*λ*Γ/(rr*(rr+k*λ))
        uθ = - Pbyη*sinθ/(8*π*rr)
        uxy = ur*sinθ + uθ*cosθ # parallel xy-component of velocity
        uz = ur*cosθ - uθ*sinθ # normal z-component of velocity
        ux, uy = uxy*x/xy, uxy*y/xy # resolved x, y components
        u = np.array((ux, uy, uz))
    return u


def oblate_spheroidal_calc(x,y,z):
    rho=np.sqrt(x**2+y**2)
    b=1.0-z**2-rho**2
    lambda_sq=0.5*(-b+np.sqrt(b**2+4*z**2))
    OS_lambda=np.sqrt(lambda_sq)
    OS_zeta=z/OS_lambda
    return OS_zeta,OS_lambda

def Sampsonff_calc(x,y,z,Q,R1):
    OS_zeta,OS_lambda=oblate_spheroidal_calc(x,y,z)
    xy = np.sqrt(x**2 + y**2) # in-plane distance = r sinθ
    prefactor=3.0*Q/(2.0*np.pi*R1**3)
    v_z=prefactor*OS_zeta**3/(OS_lambda**2+OS_zeta**2)
    term1=OS_lambda*OS_zeta**2/(OS_lambda**2+OS_zeta**2)
    term2=np.sqrt((1.0-OS_zeta**2)/(1.0+OS_lambda**2))
    v_rho=prefactor*term1*term2
    v_x=v_rho*x/xy
    v_y=v_rho*y/xy
    return v_x,v_y,v_z

def drift_hole(r,Q,R1,k,Γ,λhalf,rcutoff): # return drift velocity at a given position (hole in wall case)
    x, y, z = r[:] # z will be normal distance = r cosθ
    xy = np.sqrt(x**2 + y**2) # in-plane distance = r sinθ
    r = np.sqrt(xy**2 + z**2) # radial distance from origin
    cosθ, sinθ = z/r, xy/r # polar angle, as cos and sin
    if r < rcutoff: # crude regularisation
        ur, uθ = 0, 0
    else: # note that k --> 2k in this problem because salt diffusion into half-space
#        ur = 3*Q*cosθ**2/(2*π*r**2) - k*λhalf*Γ/(r*(r+k*λhalf))
        ur = - k*λhalf*Γ/(r*(r+k*λhalf))
        uθ = 0.0
    uxy = ur*sinθ + uθ*cosθ # parallel xy-component of velocity
    uz = ur*cosθ - uθ*sinθ # normal z-component of velocity
    ux, uy = uxy*x/xy, uxy*y/xy # resolved x, y components
# obtain full Sampson flow field
    v_x,v_y,v_z = Sampsonff_calc(x,y,z,Q,R1)
    ux=ux+v_x
    uy=uy+v_y
    uz=uz+v_z
    return np.array((ux, uy, uz))



def calc_fixedpts(Q,k,Γ,Ds,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff):
    min_kGamma_over_D=(1.0+np.sqrt(k*λstar))**2
    kGamma_over_D=k*Γ/Ds
#
    if(kGamma_over_D > min_kGamma_over_D ):
        b=1.0+k*λstar-kGamma_over_D
        a=1.0
        c=k*λstar
        rho_fp_minus=(-b-np.sqrt(b**2-4.0*a*c))/(2.0*a)
        rho_fp_plus=(-b+np.sqrt(b**2-4.0*a*c))/(2.0*a)
    else:
        rho_fp_minus=1.0e3
        rho_fp_plus=1.0e3
    return rho_fp_minus*rstar,rho_fp_plus*rstar

    '''
    def drift_cartesian(r,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff):
    x, y, z = r[:]
    rr = np.sqrt(x**2 + y**2 + z**2) # radial distance from origin
    if rr < rcutoff: # crude regularisation
        u = np.zeros_like(r)
    else:
        A = Pbyη/(8*π*rr**3)
        B = Q/(4*π*rr**3) - k*λ*Γ/(rr**2*(rr+k*λ))
        u = np.array((A*x*z + B*x, A*y*z + B*y, A*(x**2 + y**2 + 2*z**2) + B*z))
    return u
    '''