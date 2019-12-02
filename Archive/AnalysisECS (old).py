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

def fhh(sc,a,c):
    N = 60
    t = 1.0
    bo = 0.5799
    Rg = 4.6257
    L = (N-1)*bo
    return a * (((t * bo * sc * (Rg**2.0))/(L))**(1.0/3.0)) * L + c 


def heights(z, edps, zhis):
    m,n = edps.shape
    hs = np.zeros(n)
    for i in range(n):
        hs[i] = brush_height(z,edps[:,i],zhis[i])
    return hs

def fit_spline(X,Y):
    xnew = np.linspace(X[0], X[-1],80)
    spl = scipy.interpolate.make_interp_spline(X,Y,k=3)
    ynew = spl(xnew)
    return xnew, ynew

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
    
#     
# def dp_z2d2(zZ,dp_z,zhi):
#     s = np.sum(dp_z)
#     t0 = 0.0001*s/100
#     t1 = 99.9999*s/100
#     cs = np.cumsum(dp_z)
#     id0 = np.where(cs < t0)
#     id1 = np.where(cs < t1)
#     i0_del = id0[0]
#     i1_del = id1[0]
#     i0 = i0_del[-1]
#     i1 = i1_del[-1]
#     zs = (zZ[i0:i1] - zZ[i0])*zhi
#     D = (zZ[i1]-zZ[i0])*zhi
#     zD = zs/D
#     dp_d = dp_z[i0:i1] 
# 
#     return zD,dp_d 
       
def dp_z2d(zZ,dp_z,zhi,i0,i1):
    zs = (zZ[i0:i1] - zZ[i0])*zhi
    D = (zZ[i1]-zZ[i0])*zhi
    zD = zs/D
    dp_d = dp_z[i0:i1] 
    return zD,dp_d        


def find_PD(path,col):
    os.chdir(path)
    dirs = filter(os.path.isdir ,os.listdir('.'))
    n = len(dirs)
    Ps = np.zeros(n)
    Ds = np.zeros(n)
    i = 0
    for dir in dirs:  
        S = dir.split("_")
        Ps[i] = float(S[0])
        a = dir + "\comp.csv"
        b = '.\\' + a
        equil_df = pd.pandas.read_csv(b)
        data = equil_df.values[:,col]
        Ds[i] = np.mean(data[-200:])
        i = i + 1 
    Ps,Ds = zip(*sorted(zip(Ps,Ds)))
    return Ps,Ds

    
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
            data = df.values[:,col]
            Cs[i,j] = np.mean(data[-10:])
         #   print(Cs[i,j])
            j = j + 1
        i = i + 1
    A = zip(Vs,Cs)    
    B = sorted(A)
    Vs,Cs = list(zip(*B))

    os.chdir(ocwd)
    return Vs,Cs


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
    
    #Ps,Ds = find_PD(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Comp\30_N\MDPBB - ECS-30_M_30_N',27)
    
    Vs,Cs = read_thermo(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Comp\30_N\MDPBB - ECS-30_M_30_N',[27,28,29],3)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #---T,V and density profiles---#
    
    zDDp,zDTp,zDVp = read_tvdp(r'C:\Users\klay_\OneDrive - Imperial College London\TSM CDT\Major Project\My LAMMPS\MDPBB - Shear\0.005_dt\3_P\MDPBB - ECS-X_M_Y_N\1_V')
    
    print(zDTp)
   
    fig, ax1 = plt.subplots(figsize=(10,8))
    Tp, = ax1.plot(zDTp[0],zDTp[1], 'r')
    Dp, = ax1.plot(zDDp[0],zDDp[1], 'g')
    ax2 = ax1.twinx()
    Vp, = ax2.plot(zDVp[0],zDVp[1], 'b')

    #plt.legend([Tp,Vp,Dp], [r'Temperature Profile','Velocity Profile','Density Profile'], loc=1)
    ax1.set(xlim=[0,1],xlabel='$z/D$', ylabel= r'$T, \phi(z)$' , title='')
    ax2.set_ylabel(r'$V$')
    ax1.minorticks_on()
    ax2.minorticks_on()
    plt.savefig('TnVprofiles.jpg')    
 #    
 # 
 #    
    
    
  
    
    
    
    
    
