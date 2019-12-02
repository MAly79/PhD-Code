# This is a python script that will generate a LAMMPS molecule file for use in
# Polymer Brush
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy




def read_edp(filename,sheet):
    edp_df = pd.pandas.read_csv(filename, sheet_name=sheet, header = None)
    data = edp_df.values
    Ms = data[0,1:]
    zhis = data[1,1:]
    z = data[2:,0]
    edps = data[2:,1:]

    #zhi = edp_df
    return z,edps,zhis,Ms

def brush_height(z,dp,zhi):
    z0 = 1
    s = np.sum(dp)
    t0 = 1*s/100
    t1 = 80*s/100
    cs = np.cumsum(dp)
    id0 = np.where(cs < t0)
    id1 = np.where(cs < t1)
    i0 = id0[0]
    i1 = id1[0]
    z0 = z[i0[-1]]
    z1 = z[i1[-1]]
    h = (z1 - z0) * zhi

    return h

def heights(z, edps, zhis):
    m,n = edps.shape
    hs = np.zeros(n)
    for i in range(n):
        hs[i] = brush_height(z,edps[:,i],zhis[i])
    return hs
    
def fhh(sc,a,c):
    N = 60
    t = 1.0
    bo = 0.5799
    Rg = 4.6257
    L = (N-1)*bo
    return a * (((t * bo * sc * (Rg**2.0))/(L))**(1.0/3.0)) * L + c 
    
def fit_spline(X,Y):
    xnew = np.linspace(X[0], X[-1],80)
    spl = scipy.interpolate.make_interp_spline(X,Y,k=3)
    ynew = spl(xnew)
    return xnew, ynew


def dp_ids(dp):
    s = np.sum(dp)
    t0 = 0.001*s/100
    t1 = 99.9999*s/100
    cs = np.cumsum(dp)
    id0 = np.where(cs < t0)
    id1 = np.where(cs < t1)
    i0_del = id0[0]
    i1_del = id1[0]
    i0 = i0_del[-1]
    i1 = i1_del[-1]

    return i0,i1
       
def dp_z2d(zZ,dp_z,zhi,i0,i1):
    zs = (zZ[i0:i1] - zZ[i0])*zhi
    D = (zZ[i1]-zZ[i0])*zhi
    zD = zs/D
    dp_d = dp_z[i0:i1] 
    return zD,dp_d        

    
def read_dps(path,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    if stage == 1:
        a = 'equil.csv'
        b = 'abdpe.csv'
        c = 'bbdpe.csv'
        d = 'tbdpe.csv'
    elif stage == 2:
        a = 'comp.csv'
        b = 'abdpc.csv'
        c = 'bbdpc.csv'
        d = 'tbdpc.csv'
    elif stage == 3:
        a = 'shear.csv'
        b = 'abdps.csv'
        c = 'bbdps.csv'
        d = 'tbdps.csv'
    zhi = (pd.pandas.read_csv(a).values)[-1,28]
    
    abeads = (pd.pandas.read_csv(b)).values
    zZ = abeads[1:,1]
    adp_z = abeads[1:,-1]
    i0,i1 = dp_ids(adp_z)
    zD,adp = dp_z2d(zZ,adp_z,zhi,i0,i1)
    zDadp = [zD, adp]
    
    bbeads = (pd.pandas.read_csv(c)).values
    zZ = bbeads[1:,1]
    bdp_z = bbeads[1:,-1]
    zD,bdp = dp_z2d(zZ,bdp_z,zhi,i0,i1)
    zDbdp = [zD, bdp]
    
    tbeads = (pd.pandas.read_csv(d)).values
    zZ = tbeads[1:,1]
    tdp_z = tbeads[1:,-1]
    zD,tdp = dp_z2d(zZ,tdp_z,zhi,i0,i1)
    zDtdp = [zD, tdp]
    
    #Calculate interpenetration
    I = np.trapz((bdp*tdp),zD)
    D = (zZ[i1]-zZ[i0])*zhi
    os.chdir(ocwd)
    return zDbdp,zDtdp,zDadp,I,D
    
def IvsD(path):
    ocwd = os.getcwd()
    os.chdir(path)
    dirs = filter(os.path.isdir ,os.listdir('.'))
    n = len(dirs)
    Ps = np.zeros(n)
    Is = np.zeros(n)
    Ds = np.zeros(n)
    i=0
    for dir in dirs:  
        S = dir.split("_")
        Ps[i] = float(S[0])
        c = '.\\' + dir
        x,y,z,Is[i],Ds[i] = read_dps(c,2)
        i = i + 1 
    Ps,Is,Ds = zip(*sorted(zip(Ps,Is,Ds)))
    os.chdir(ocwd)
    return Ps,Is,Ds

def read_thermo(path,cols,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    dirs = list(filter(os.path.isdir ,os.listdir('.')))
    #print(dirs)
    n = len(dirs)
    m = len(cols)
    Vs = np.zeros(n)
    Cs = np.zeros((n,m))
    if stage == 1:
        f = "\equil.csv"
    elif stage == 2:
        f = "\comp.csv"
    elif stage == 3:
        f = "\shear.csv"
    i = 0    
    for dir in dirs:
        S = dir.split("_")
        Vs[i] = float(S[0])
        #print(Vs[i])
        a = dir + f
        b = '.\\' + a
        df = pd.pandas.read_csv(b)
        j = 0
        for col in cols:
            data = df.values[:,col+1]
            Cs[i,j] = np.mean(data[-200:])
         #   print(Cs[i,j])
            j = j + 1
        i = i + 1
    A = zip(Vs,Cs)    
    B = sorted(A)
    Vs,Cs = list(zip(*B))
    Cs = np.asarray(Cs)
    os.chdir(ocwd)
    return Vs,Cs
    
def thermo_time(path,cols,stage):
    ocwd = os.getcwd()
    os.chdir(path) 
    cols = list(np.asarray(cols) + 1)
    #print(cols)
    if stage == 1:
        f = "equil.csv"
    elif stage == 2:
        f = "comp.csv"
    elif stage == 3:
        f = "shear.csv"
    df = pd.pandas.read_csv(f)
    dt = df.iloc[:,1]
    ts = dt.values
    vdf = df.iloc[:,cols]
    Vt = vdf.values
    os.chdir(ocwd)
    return ts,Vt

def read_tvdp(path):
    ocwd = os.getcwd()
    os.chdir(path)
    zhi = (pd.pandas.read_csv('shear.csv').values)[-1,28]
    
    abeads = (pd.pandas.read_csv('abdps.csv')).values
    zZ = abeads[1:,1]
    adpZ = abeads[1:,-1]
    i0,i1 = dp_ids(adpZ)
    zD,adp = dp_z2d(zZ,adpZ,zhi,i0,i1)
    zDDp = [zD, adp]
    
    temps = (pd.pandas.read_csv('temps.csv')).values
    zZ = temps[1:,1]
    Tdata = temps[1:,-1]
    zD,Tp = dp_z2d(zZ,Tdata,zhi,i0,i1)
    zDTp = [zD, Tp]
    
    velps = (pd.pandas.read_csv('velps.csv')).values
    zZ = velps[1:,1]
    Vdata = velps[1:,-1]
    zD,Vp = dp_z2d(zZ,Vdata,zhi,i0,i1)
    zDVp = [zD, Vp]

    
    os.chdir(ocwd)
    return zDDp,zDTp,zDVp

if __name__ == '__main__':
    

    #  0     1   2   3   4      5       6       7           8       9     10   11  12  13  14  15  16      17         18          19         20     21       22      23    24      25   26     27         28            29         30          31    32           
    #step etotal ke pe epair temp c_melTemp c_wallTemp v_Fcatom v_Pcomp press pxx pyy pzz pxy pxz pyz c_melPress c_wallPress  v_melDens v_surfcov v_aveRg v_Wall_v v_srate v_D v_bwzmax zhi c_fbwall[1] c_fbwall[3] c_ftwall[1] c_ftwall[3] v_sbot v_pbot
    
    
    
      #---T,V and density profiles---#
      
    fig, ax1 = plt.subplots(figsize=(10,8))
    ax2 = ax1.twinx()
    
    zDDp20,zDTp20,zDVp20 = read_tvdp(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 1_P\1_V')
    Tp20, = ax1.plot(zDTp20[0],zDTp20[1], 'r--')
    Dp20, = ax1.plot(zDDp20[0],zDDp20[1], 'g--')
    Vp20, = ax2.plot(zDVp20[0],zDVp20[1], 'b--')

    plt.legend([Tp20,Vp20,Dp20], [r'Temperature','Velocity','Density'], loc=2)
    ax1.set(xlim=[0,1],xlabel='$z/D$', ylabel= r'$T, \phi(z)$')
    ax2.set_ylabel(r'$V$')
    ax1.minorticks_on()
    ax2.minorticks_on()
    plt.savefig('TnVprofiles.jpg')    
#     
#     #--- Convergence check T---#
#     
#     ts,Vt = thermo_time(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\0.005_dt\1_P\MDPBB - ECS-30_M_30_N\1_V',[5,35],3)
#     
#     #print(ts)
#     #print(Vt)
#     
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     l1, = ax1.plot(ts,Vt[:,0], 'r')
#     #ax2 = ax1.twinx()
#     #l2, = ax2.plot(ts,Vt[:,1], 'b')
#     ax1.minorticks_on()
#     #ax2.minorticks_on()
#     plt.savefig('Time Evolution.jpg') 
#     
#     #---Plot of the variation of D with V---#
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     
#     Vs1,Cs1 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 1_P',[13,15,27],3)    
#     Vs3,Cs3 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 3_P',[13,15,27],3)    
#     Vs8,Cs8 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 8_P',[13,15,27],3)    
#     l1, = ax1.plot(Vs1,Cs1[:,2], 'ro')
#     l3, = ax1.plot(Vs3,Cs3[:,2], 'go')
#     l8, = ax1.plot(Vs8,Cs8[:,2], 'bo')
#     # Vs1,Cs1 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_30_N_1_P',[13,15,27],3)    
#     # Vs3,Cs3 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_30_N_3_P',[13,15,27],3)    
#     # Vs8,Cs8 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_30_N_8_P',[13,15,27],3)    
#     # l1l, = ax1.plot(Vs1,Cs1[:,2], 'rs')
#     # l3l, = ax1.plot(Vs3,Cs3[:,2], 'gs')
#     # l8l, = ax1.plot(Vs8,Cs8[:,2], 'bs')
#     plt.legend([l1,l3,l8], ['$P_{comp}=1$','$P_{comp}=3$', '$P_{comp}=8$'], loc=9)
#     ax1.set(xlabel='$V_{wall}$', ylabel= r'D')
#     plt.savefig('D vs V.jpg')
#     
#     
#     #---Shear and Normal Pressure plot vs V_Wall for different Pcomp---#
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     ax2 = ax1.twinx()
#     
#     
#     Vs1,Cs1 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 1_P',[13,15,27,26],3)
#     n1, = ax1.plot(Cs1[:,3],Cs1[:,0], 'ro')
#     s1, = ax2.plot(Cs1[:,3],Cs1[:,1], 'bo')
#     
#     Vs3,Cs3 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 3_P',[13,15,27,26],3)
#     n3, = ax1.plot(Cs3[:,3],Cs3[:,0], 'rs')
#     s3, = ax2.plot(Cs3[:,3],Cs3[:,1], 'bs')
#     
#     Vs8,Cs8 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\High P more Vs\MDPBB - ECS-30_M_60_N - 8_P',[13,15,27,26],3)
#     n8, = ax1.plot(Cs8[:,3],Cs8[:,0], 'r^')
#     s8, = ax2.plot(Cs8[:,3],Cs8[:,1], 'b^')    
#     
#     plt.legend([n1,s1], ['$P_{zz}$','$P_{xz}$'], loc=9)
#     ax1.set(xlabel=r'$\.\gamma$', ylabel= '$P_{zz}$')
#     ax2.set(ylabel = '$P_{xz}$')
#     ax1.minorticks_on()
#     ax2.minorticks_on()
#     plt.savefig('PzzPxz vs srate - vP.jpg')
# 
# 
#     #---Plot the COF vs srate---#
#     
#     mu1 = -1*Cs1[:,1]/Cs1[:,0]
#     mu3 = -1*Cs3[:,1]/Cs3[:,0]
#     mu8 = -1*Cs8[:,1]/Cs8[:,0]
#     
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     l1, = ax1.plot(Cs1[:,3],mu1, 'bo')
#     l3, = ax1.plot(Cs3[:,3],mu3, 'rs')
#     l8, = ax1.plot(Cs8[:,3],mu8, 'g^')
#     plt.legend([l1,l3,l8], ['$P_{comp}=1$, $D=11-16$','$P_{comp}=3$, $D=8-12$', '$P_{comp}=8$, $D=6-8$'], loc=1)    
#     
#     
#     ax1.minorticks_on()
#     ax1.set(xlabel=r'$\.\gamma$', ylabel= r'$\mu$' , xlim=[0,0.25])
# 
#     plt.savefig('COF vs srate - vP.jpg')
#     
#      #---Plot the COF vs srate---#
#     
#     mu1 = -1*Cs1[:,1]/Cs1[:,0]
#     mu3 = -1*Cs3[:,1]/Cs3[:,0]
#     mu8 = -1*Cs8[:,1]/Cs8[:,0]
#     
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     l1, = ax1.plot(Vs1,mu1, 'bo')
#     l3, = ax1.plot(Vs3,mu3, 'rs')
#     l8, = ax1.plot(Vs8,mu8, 'g^')
#     plt.legend([l1,l3,l8], ['$P_{comp}=1$, $D=11-16$','$P_{comp}=3$, $D=8-12$', '$P_{comp}=8$, $D=6-8$'], loc=1)    
#     
#     
#     ax1.minorticks_on()
#     ax1.set(xlabel=r'$V_{wall}$', ylabel= r'$\mu$')
# 
#     plt.savefig('COF vs Vwall - vP.jpg')
#     #print(Vs)
#     #print(np.asarray(Cs))
#     
#    
#     #---Shear and Normal Pressure plot vs V_Wall for different grafting densities---#
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     ax2 = ax1.twinx()
#     
#     
#     Vs1,Cs1 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\Grafting Densities - 3_P\MDPBB - ECS-20_M_30_N',[13,15,27,26],3)
#     n1, = ax1.plot(Cs1[:,3],Cs1[:,0], 'ro')
#     s1, = ax2.plot(Cs1[:,3],Cs1[:,1], 'bo')
#     
#     Vs3,Cs3 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\Grafting Densities - 3_P\MDPBB - ECS-30_M_30_N',[13,15,27,26],3)
#     n3, = ax1.plot(Cs3[:,3],Cs3[:,0], 'rs')
#     s3, = ax2.plot(Cs3[:,3],Cs3[:,1], 'bs')
#     
#     Vs8,Cs8 = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\Grafting Densities - 3_P\MDPBB - ECS-60_M_30_N',[13,15,27,26],3)
#     n8, = ax1.plot(Cs8[:,3],Cs8[:,0], 'r^')
#     s8, = ax2.plot(Cs8[:,3],Cs8[:,1], 'b^')    
#     
#     plt.legend([n1,s1], ['$P_{zz}$','$P_{xz}$'], loc=9)
#     ax1.set(xlabel=r'$\.\gamma$', ylabel= '$P_{zz}$')
#     ax2.set(ylabel = '$P_{xz}$')
#     ax1.minorticks_on()
#     ax2.minorticks_on()
#     plt.savefig('PzzPxz vs srate - vM.jpg')
# 
# 
#     #---Plot the COF vs srate---#
#     
#     mu1 = -1*Cs1[:,1]/Cs1[:,0]
#     mu3 = -1*Cs3[:,1]/Cs3[:,0]
#     mu8 = -1*Cs8[:,1]/Cs8[:,0]
#     
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     l1, = ax1.plot(Cs1[:,3],mu1, 'bo')
#     l3, = ax1.plot(Cs3[:,3],mu3, 'rs')
#     l8, = ax1.plot(Cs8[:,3],mu8, 'g^')
#     plt.legend([l1,l3,l8], [r'$\rho = 0.05$',r'$\rho = 0.075$', r'$\rho = 0.15$'], loc=1)     
#     
#     ax1.minorticks_on()
#     ax1.set(xlabel=r'$\.\gamma$', ylabel= r'$\mu$', xlim=[0,0.5])
# 
#     plt.savefig('COF vs srate - vM.jpg')
#     

    
  
    
    
    
    
    
