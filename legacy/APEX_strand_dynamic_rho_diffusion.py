"""
APEX 0.2 mm insulated-Al strand: time-domain finite magnetic diffusion.
Independent validation of the complete-penetration strand-loss formula.
Uses the current 3/1/3 ms trapezoidal field and rho(4.2 K,B) magnetoresistance.
The Phase-I design itself intentionally retains the complete-penetration formula
for conservatism; this script is validation-only.
"""
import numpy as np
import pandas as pd
from scipy.linalg import solve_banded

mu0=4*np.pi*1e-7
d=0.2e-3; a=d/2
Bpk=20.619192
Tr,Tf,Td=3e-3,1e-3,3e-3
Tpulse=Tr+Tf+Td
rho0=9.832581476e-12
rho293=2.94978e-8
Krho=rho293/rho0

def rho_B(B):
    H=0.01*np.abs(B)*Krho
    MR=H**2*(1+0.00177*H)/(1.8+1.6*H+0.53*H**2)
    return rho0*(1+MR)

def Btrap(t):
    if t<0: return 0.0
    if t<Tr: return Bpk*t/Tr
    if t<Tr+Tf: return Bpk
    if t<=Tpulse: return Bpk*(Tpulse-t)/Td
    return 0.0

def build_operator(N):
    dr=a/N
    r=np.arange(N)*dr
    lo=np.zeros(N); dg=np.zeros(N); up=np.zeros(N); b=np.zeros(N)
    dg[0]=-8/dr**2; up[0]=8/dr**2
    q=a/(2*dr); den=2+3*q
    alphaB=2/den; c1=4*q/den; c2=-q/den
    for i in range(1,N):
        ri=r[i]
        cm=1/dr**2-3/(2*ri*dr)
        c0=-2/dr**2
        cp=1/dr**2+3/(2*ri*dr)
        lo[i]+=cm; dg[i]+=c0
        if i<N-1: up[i]+=cp
        else:
            dg[i]+=cp*c1; lo[i]+=cp*c2; b[i]+=cp*alphaB
    return r,dr,lo,dg,up,b

def tri_mv(lo,dg,up,x):
    y=dg*x
    y[1:]+=lo[1:]*x[:-1]
    y[:-1]+=up[:-1]*x[1:]
    return y

def solve_case(N=160,dt=1e-6,tend=30e-3):
    r,dr,lo,dg,up,b=build_operator(N)
    g=np.zeros(N)
    Qpulse=Qtail=Pmax=0.0; tPmax=0.0
    for n in range(int(round(tend/dt))):
        tm=(n+0.5)*dt; Bm=Btrap(tm); rho=rho_B(Bm); D=rho/mu0; c=0.5*dt*D
        rhs=g+c*tri_mv(lo,dg,up,g)+dt*D*b*Bm
        ab=np.zeros((3,N)); ab[0,1:]=-c*up[:-1]; ab[1,:]=1-c*dg; ab[2,:-1]=-c*lo[1:]
        gnew=solve_banded((1,1),ab,rhs)
        gmid=0.5*(g+gnew)
        Lg=tri_mv(lo,dg,up,gmid)+b*Bm
        Jamp=-(r/mu0)*Lg
        Jb=2*Jamp[-1]-Jamp[-2]
        rr=np.r_[r,a]; JJ=np.r_[Jamp,Jb]
        P=np.pi*rho*np.trapezoid(JJ**2*rr,rr)
        if tm<=Tpulse: Qpulse+=P*dt
        else: Qtail+=P*dt
        if P>Pmax: Pmax=P; tPmax=tm
        g=gnew
    return Qpulse,Qtail,Qpulse+Qtail,Pmax,tPmax

dt_ref=0.1e-6
tt=np.arange(dt_ref/2,Tpulse,dt_ref)
BB=np.array([Btrap(t) for t in tt])
dBdt=np.zeros_like(tt)
dBdt[tt<Tr]=Bpk/Tr
dBdt[(tt>=Tr+Tf)&(tt<Tpulse)]=-Bpk/Td
Qcp=np.sum((np.pi*a**4/(4*np.array([rho_B(B) for B in BB])))*dBdt**2*dt_ref)

rows=[]
for N,dt in [(80,5e-6),(120,2e-6),(160,1e-6)]:
    Qp,Qt,Qall,Pmax,tP=solve_case(N,dt)
    rows.append({'Nr':N,'dt_us':dt*1e6,'Q_0to7ms_J_per_m':Qp,
                 'Q_tail_after7ms_J_per_m':Qt,'Q_total_J_per_m':Qall,
                 'Qtotal_over_QCP':Qall/Qcp,'Ppeak_W_per_m':Pmax,'t_Ppeak_ms':tP*1e3})

df=pd.DataFrame(rows)
df.to_csv('APEX_strand_dynamic_rho_diffusion_convergence.csv',index=False)
best=df.iloc[-1]
summary=pd.DataFrame([{
    'Bpeak_T':Bpk,'Qcp_dynamic_rho_J_per_m':Qcp,
    'Qdiff_0to7ms_J_per_m':best.Q_0to7ms_J_per_m,
    'Qtail_7to30ms_J_per_m':best.Q_tail_after7ms_J_per_m,
    'Qdiff_total_J_per_m':best.Q_total_J_per_m,
    'Qtotal_over_Qcp':best.Q_total_J_per_m/Qcp,
    'reduction_vs_CP_pct':100*(best.Q_total_J_per_m/Qcp-1),
    'COMSOL_reference_0to7ms_J_per_m':0.5809452522,
    'Python_vs_COMSOL_pct':100*(best.Q_0to7ms_J_per_m/0.5809452522-1),
}])
summary.to_csv('APEX_strand_dynamic_rho_diffusion_summary.csv',index=False)
print(df.to_string(index=False))
print('\n',summary.to_string(index=False))
