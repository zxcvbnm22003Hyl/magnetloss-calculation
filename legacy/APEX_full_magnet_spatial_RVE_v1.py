# -*- coding: utf-8 -*-
"""
APEX full-magnet spatial RVE / global-local FE² solver
======================================================

Current APEX baseline:
- 162 turns = 18 SP × 9 radial turns
- Ri = 125 mm, cable = 8.3 mm (r) × 16.9 mm (z)
- Ipk = 47.29 kA
- 3/1/3 ms trapezoid
- 2432 strands/cable, d = 0.2 mm
- Al RRR = 3000

What this code does
-------------------
1) Builds the entire 162-turn magnet geometry.
2) Evaluates the macro magnetic field at 4 Gauss points inside every turn
   (162 × 4 = 648 spatial RVE driving points) using distributed circular-loop
   Biot-Savart quadrature.
3) At every RVE point solves the finite magnetic diffusion problem in a
   0.2-mm round strand (m=1 transverse-field mode).
4) Couples local strand eddy loss to per-turn temperature and rho(T,B).
5) Applies the previously obtained hierarchical multistrand RVE corrections
   (19 -> 304 -> 2432) as small magnetic-interaction correction factors.
6) Integrates all 162 turns to obtain whole-magnet AC loss [J/pulse].

Important
---------
- This is the "full magnet spatial RVE" calculation, not the local 1.42 kJ/m
  hotspot cable result.
- The multistrand correction factors are fitted from the previously computed
  explicit RVE scans and are ~O(1%).
- Inter-strand electrical coupling loss is NOT included. All strands are
  assumed electrically insulated.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from numpy.polynomial.legendre import leggauss

DT = 5e-6
NR_DIFF = 48
RHO_CSV = "APEX_Al0p2_RRR3000_rho_TB_long.csv"
mu0 = 4*np.pi*1e-7
Ipk = 47.29e3
nDP = 9
nSP = 18
nR  = 9
Ri = 125e-3
wr = 8.3e-3
wz = 16.9e-3
gr = 2.0e-3
gDPin = 0.5e-3
gDP = 3.0e-3
pitchR = wr + gr
Ns = 2432
dstrand = 0.2e-3
astr = dstrand/2
Astr = np.pi*dstrand**2/4
AAl = Ns*Astr
rho_mass = 2700.0
Tr = 3e-3
Tf = 1e-3
Td = 3e-3
Tp = Tr+Tf+Td

if not os.path.isfile(RHO_CSV):
    alt = os.path.join(os.path.dirname(os.path.abspath(__file__)), RHO_CSV)
    if os.path.isfile(alt):
        RHO_CSV = alt
    else:
        raise FileNotFoundError(RHO_CSV)

rt = pd.read_csv(RHO_CSV)
Ts = np.sort(rt["T_K"].unique())
Bs = np.sort(rt["B_T"].unique())
R = rt.pivot(index="T_K", columns="B_T", values="rho_TB_Ohm_m").loc[Ts,Bs].values
rho_interp = RegularGridInterpolator((Ts,Bs), R, bounds_error=False, fill_value=None)

def rhoTB(T,B):
    T = np.clip(np.asarray(T), Ts[0], Ts[-1])
    B = np.clip(np.asarray(B), Bs[0], Bs[-1])
    shp = np.broadcast(T,B).shape
    Tv = np.broadcast_to(T,shp).ravel()
    Bv = np.broadcast_to(B,shp).ravel()
    return rho_interp(np.column_stack([Tv,Bv])).reshape(shp)

def cpAl(T):
    T=np.maximum(np.asarray(T),4.0)
    x=np.log10(T)
    y=(46.6467 -314.292*x +866.662*x**2 -1298.3*x**3
       +1162.27*x**4 -637.795*x**5 +210.351*x**6
       -38.3094*x**7 +2.96344*x**8)
    return 10**y

Tgrid=np.linspace(4.2,293.15,20000)
hgrid=np.concatenate([[0.0], cumulative_trapezoid(cpAl(Tgrid),Tgrid)])
def T_from_h(h):
    return np.interp(h,hgrid,Tgrid)

def waveform(t):
    t=np.asarray(t)
    f=np.zeros_like(t,dtype=float)
    m=(t>=0)&(t<Tr)
    f[m]=t[m]/Tr
    m=(t>=Tr)&(t<Tr+Tf)
    f[m]=1.0
    m=(t>=Tr+Tf)&(t<=Tp)
    f[m]=1-(t[m]-(Tr+Tf))/Td
    return f

Lz=nDP*(2*wz+gDPin)+(nDP-1)*gDP
zcent=[]
zcur=-Lz/2
for idp in range(nDP):
    zcent.append(zcur+wz/2); zcur += wz+gDPin
    zcent.append(zcur+wz/2); zcur += wz
    if idp<nDP-1:
        zcur += gDP
zcent=np.asarray(zcent)
rcent=Ri+wr/2+np.arange(nR)*pitchR
length=2*np.pi*np.broadcast_to(rcent[None,:],(nSP,nR))
mass=rho_mass*AAl*length

def loop_fields_vectorized(a,zs,sw,r,z):
    dz=z-zs
    m=4*a*r/((a+r)**2+dz**2)
    m=np.clip(m,0.0,1.0-1e-12)
    K=ellipk(m)
    E=ellipe(m)
    den=np.sqrt((a+r)**2+dz**2)
    delta=(a-r)**2+dz**2
    Br=mu0*dz/(2*np.pi*r*den)*(-K+(a*a+r*r+dz*dz)/delta*E)
    Bz=mu0/(2*np.pi*den)*(K+(a*a-r*r-dz*dz)/delta*E)
    return np.sum(sw*Br), np.sum(sw*Bz)

gs,ws=leggauss(4)
src=[]
for zc in zcent:
    for rc in rcent:
        for i,xi in enumerate(gs):
            rs=rc+0.5*wr*xi
            for j,eta in enumerate(gs):
                zs=zc+0.5*wz*eta
                src.append((rs,zs,(ws[i]/2)*(ws[j]/2)))
src=np.asarray(src)
a_src,z_src,w_src=src[:,0],src[:,1],src[:,2]

gt,wt=leggauss(2)
targets=[]
meta=[]
for iz,zc in enumerate(zcent):
    for ir,rc in enumerate(rcent):
        for xi in gt:
            for eta in gt:
                targets.append((rc+0.5*wr*xi, zc+0.5*wz*eta))
                meta.append((iz,ir))
targets=np.asarray(targets)

Bpk=np.zeros((len(targets),2))
for k,(r,z) in enumerate(targets):
    Br1,Bz1=loop_fields_vectorized(a_src,z_src,w_src,r,z)
    Bpk[k]=Ipk*np.array([Br1,Bz1])

Bpk4=Bpk.reshape(nSP,nR,4,2)
Bmagpk=np.linalg.norm(Bpk4,axis=-1)

def C1(B):
    return np.interp(B,[0,5,10,20.6],[1.0100,1.0070,1.0055,1.004560835])

def C2(B):
    return np.interp(B,[0,10,20.5],[1.00831,1.00188,1.00115])

def C3r(B):
    return np.interp(B,[0,5,10,15,20.5],[1.0065,1.000640,1.000444,1.000296,1.000160])

def C3z(B):
    return np.interp(B,[0,5,10,15,20.5],[0.9900,0.997271,0.997597,0.997849,0.998082])

Bflat=np.linalg.norm(Bpk,axis=1)
Br=Bpk[:,0]
Bz=Bpk[:,1]
wr_dir=Br**2/(Bflat**2+1e-30)
wz_dir=1.0-wr_dir
Cmulti=C1(Bflat)*C2(Bflat)*(wr_dir*C3r(Bflat)+wz_dir*C3z(Bflat))

def batched_tridiag(lower,main,upper,rhs):
    M,N=main.shape
    cp=np.empty((M,N-1))
    dp=np.empty((M,N))
    den=main[:,0]
    cp[:,0]=upper[:,0]/den
    dp[:,0]=rhs[:,0]/den
    for i in range(1,N-1):
        den=main[:,i]-lower[:,i-1]*cp[:,i-1]
        cp[:,i]=upper[:,i]/den
        dp[:,i]=(rhs[:,i]-lower[:,i-1]*dp[:,i-1])/den
    den=main[:,-1]-lower[:,-1]*cp[:,-1]
    dp[:,-1]=(rhs[:,-1]-lower[:,-1]*dp[:,-2])/den
    x=np.empty_like(rhs)
    x[:,-1]=dp[:,-1]
    for i in range(N-2,-1,-1):
        x[:,i]=dp[:,i]-cp[:,i]*x[:,i+1]
    return x

Nr=NR_DIFF
faces=np.linspace(0,1,Nr+1)
dx=1/Nr
V=(faces[1:]**4-faces[:-1]**4)/4
lo=np.zeros(Nr); di=np.zeros(Nr); up=np.zeros(Nr); bc=np.zeros(Nr)
for i in range(Nr):
    xl=faces[i]; xh=faces[i+1]; Vi=V[i]
    if i>0:
        q=xl**3/(dx*Vi)
        lo[i]=q
        di[i]-=q
    if i<Nr-1:
        q=xh**3/(dx*Vi)
        up[i]=q
        di[i]-=q
    else:
        q=xh**3/((dx/2)*Vi)
        di[i]-=q
        bc[i]=q
W3=(faces[1:]**4-faces[:-1]**4)/4

t=np.arange(0,Tp+DT/2,DT)
f=waveform(t)
u=np.zeros((len(targets),Nr))
hturn=np.zeros((nSP,nR))
Qeddy=0.0
Qj=0.0
qRVE=np.zeros(len(targets))
hist=[]
for k in range(len(t)-1):
    b0=f[k]
    b1=f[k+1]
    bm=0.5*(b0+b1)
    Tturn=T_from_h(hturn)
    Ttar=np.repeat(Tturn.reshape(-1),4)
    Bmid=Bflat*bm
    rho=rhoTB(Ttar,Bmid)
    alpha=rho/(mu0*astr**2)
    h=0.5*DT*alpha
    rhs=(1+h[:,None]*di[None,:])*u
    rhs[:,1:]+=h[:,None]*lo[None,1:]*u[:,:-1]
    rhs[:,:-1]+=h[:,None]*up[None,:-1]*u[:,1:]
    rhs+=h[:,None]*bc[None,:]*(b0+b1)
    lower=-h[:,None]*lo[None,1:]
    main=1-h[:,None]*di[None,:]
    upper=-h[:,None]*up[None,:-1]
    unew=batched_tridiag(lower,main,upper,rhs)
    du=(unew-u)/DT
    Pnorm=np.pi*astr**4/rho*np.sum(W3[None,:]*du**2,axis=1)
    Pstrand=Pnorm*(Bflat**2)*Cmulti
    Peddy_turn=Ns*Pstrand.reshape(nSP,nR,4).mean(axis=2)*length
    rhoavg=rho.reshape(nSP,nR,4).mean(axis=2)
    I=Ipk*bm
    Pj_turn=I**2*rhoavg*length/AAl
    Qeddy += Peddy_turn.sum()*DT
    Qj += Pj_turn.sum()*DT
    qRVE += Pstrand*DT
    hturn += (Peddy_turn+Pj_turn)*DT/mass
    u=unew
    Tnew=T_from_h(hturn)
    if (k % max(1,int(20e-6/DT))) == 0:
        hist.append([t[k+1]*1e3,Peddy_turn.sum()/1e6,Pj_turn.sum()/1e6,np.average(Tnew,weights=mass),Tnew.max()])
qturn=qRVE.reshape(nSP,nR,4).mean(axis=2)*Ns*length
Bturn=Bmagpk.mean(axis=2)
Tend=T_from_h(hturn)
turn_rows=[]
for iz in range(nSP):
    for ir in range(nR):
        turn_rows.append([iz+1,ir+1,rcent[ir]*1e3,zcent[iz]*1e3,Bturn[iz,ir],qturn[iz,ir],Tend[iz,ir]])
turn_df=pd.DataFrame(turn_rows,columns=["SP","radial_turn","r_mm","z_mm","Bpeak_mean_T","Qeddy_RVE_J_per_pulse","Tend_K"])
turn_df.to_csv("APEX_full_magnet_spatial_RVE_turns.csv",index=False)
hist_df=pd.DataFrame(hist,columns=["t_ms","Peddy_MW","Pj_MW","Tmean_K","Tmax_K"])
hist_df.to_csv("APEX_full_magnet_spatial_RVE_timehistory.csv",index=False)
summary=pd.DataFrame([{"N_turns":nSP*nR,"N_local_RVE_points":len(targets),"Bpeak_local_min_T":Bflat.min(),"Bpeak_local_max_T":Bflat.max(),"Qeddy_RVE_kJ":Qeddy/1e3,"Qj_kJ":Qj/1e3,"Qtotal_kJ":(Qeddy+Qj)/1e3,"Tmean_end_K":np.average(Tend,weights=mass),"Tmax_end_K":Tend.max(),"DT_us":DT*1e6,"Nr_diff":NR_DIFF}])
summary.to_csv("APEX_full_magnet_spatial_RVE_summary.csv",index=False)
print(summary.to_string(index=False))
