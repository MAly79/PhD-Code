# This is a python script that will generate a LAMMPS molecule file for use in
# Polymer Brush
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy


def brush_height(z,dp):
    z0 = 1
    s = np.sum(dp)
    t0 = 0.1*s/100
    t1 = 90*s/100
    cs = np.cumsum(dp)
    id0 = np.where(cs < t0)
    id1 = np.where(cs < t1)
    i0 = id0[0]
    i1 = id1[0]
    z0 = z[i0[-1]]
    z1 = z[i1[-1]]
    h = (z1 - z0)

    return h

def heights(z, edps):
    m,n = edps.shape
    hs = np.zeros(n)
    for i in range(n):
        hs[i] = brush_height(z,edps[:,i])
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
    
def read_thermoV(path,cols,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    dirs = list(filter(os.path.isdir ,os.listdir('.')))
    #print(dirs)
    n = len(dirs)
    m = len(cols)
    steps = 200
    Vs = np.zeros(n)
    Cs = np.zeros((n,m))
    Cd = np.zeros((n,m))
    if stage == 1:
        f = "\equil.csv"
    elif stage == 2:
        f = "\comp.csv"
    elif stage == 3:
        f = "\shear.csv"
    i = 0    
    for dir in dirs:
        S = dir.split("=")
        Vs[i] = float(S[1])
        #print(Vs[i])
        a = dir + f
        b = '.\\' + a
        df = pd.pandas.read_csv(b)
        j = 0
        for col in cols:
            data = df.values[:,col+1]
            Cs[i,j] = np.mean(data[-steps:])
            Cd[i,j] = np.std(data[-steps:])
         #   print(Cs[i,j])
            j = j + 1
        i = i + 1
    A = zip(Vs,Cs,Cd)    
    B = sorted(A)
    Vs,Cs,Cd = list(zip(*B))
    Cs = np.asarray(Cs)
    Cd = np.asarray(Cd)
    os.chdir(ocwd)
    return Vs,Cs,Cd
    
    
def read_thermoP(path,cols,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    dirs = list(filter(os.path.isdir ,os.listdir('.')))
    dirs.remove('LPBB-ECS')
    #print(dirs)
    n = len(dirs)
    m = len(cols)
    steps = 200
    Ps = np.zeros(n)
    Cs = np.zeros((n,m))
    #Vars = np.zeros((n,m))
    Cd = np.zeros((n,m))
    if stage == 1:
        f = "\equil.csv"
    elif stage == 2:
        f = "\comp.csv"
    elif stage == 3:
        f = "\shear.csv"
    i = 0    
    for dir in dirs:
        S = dir.split("=")
        #print(S)
        Ps[i] = float(S[1])
        #print(Vs[i])
        a = dir + '\LPBB-ECS' + f
        b = '.\\' + a
        df = pd.pandas.read_csv(b)
        j = 0
        for col in cols:
            data = df.values[:,col+1]
            Cs[i,j] = np.mean(data[-steps:])
            #Vars[i,j] = np.var(data[-steps:])
            #Cd[i,j] =((Vars[i,j])**0.5)
            Cd[i,j] = np.std(data[-steps:])
         #   print(Cs[i,j])
            j = j + 1
        i = i + 1
    A = zip(Ps,Cs,Cd)    
    B = sorted(A)
    Ps,Cs,Cd = list(zip(*B))
    Cs = np.asarray(Cs)
    Cd = np.asarray(Cd)
    os.chdir(ocwd)
    return Ps,Cs,Cd
    
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


def read_profs(path,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    if stage == 1:
        a = 'abdpe.csv'
        b = 'bbdpe.csv'
        c = 'tbdpe.csv'
        d = 'temps.csv'
        e = 'velps.csv'
        f = 'mope.csv'
        g = 'equil.csv'
    elif stage == 2:
        a = 'abdpc.csv'
        b = 'bbdpc.csv'
        c = 'tbdpc.csv'
        d = 'temps.csv'
        e = 'velps.csv'
        f = 'mopc.csv'
        g = 'comp.csv'
    elif stage == 3:
        a = 'abdps.csv'
        b = 'bbdps.csv'
        c = 'tbdps.csv'
        d = 'temps.csv'
        e = 'velps.csv'
        f = 'mops.csv'
        g = 'shear.csv'
    
    abeads = (pd.pandas.read_csv(a)).values
    zZ = abeads[1:,1]
    Da = abeads[1:,-1]
    zZDa = [zZ, Da]
    
    bbeads = (pd.pandas.read_csv(b)).values
    zZ = bbeads[1:,1]
    Db = bbeads[1:,-1]
    zZDb = [zZ, Db]
    
    tbeads = (pd.pandas.read_csv(c)).values
    zZ = tbeads[1:,1]
    Dt = tbeads[1:,-1]
    zZDt = [zZ, Dt]
    
    temps = (pd.pandas.read_csv(d)).values
    zZ = temps[1:,1]
    Tp = temps[1:,-1]
    zZTp = [zZ, Tp]
    
    velps = (pd.pandas.read_csv(e)).values
    zZ = velps[1:,1]
    Vp = velps[1:,-1]
    zZVp = [zZ, Vp]
    
    mops = (pd.pandas.read_csv(f)).values
    zZ = mops[1:,1]
    Pp = mops[1:,2:]
    zZPp = [zZ, Pp] 
    
    thermo = (pd.pandas.read_csv(g)).values
    zo = thermo[-1,33]
    D = thermo[-1,32]
    
    os.chdir(ocwd)
    return zZDa,zZDt,zZDb,zZVp,zZTp,zZPp,zo,D
    
def zZ_profs(path,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    zZDa,zZDt,zZDb,zZVp,zZTp,zZPp,zo,D = read_profs(path,stage)
    
    zZDa[0] = (zZDa[0]-zo)
    zZDt[0] = (zZDt[0]-zo)
    zZDb[0] = (zZDb[0]-zo)
    zZVp[0] = (zZVp[0]-zo)
    zZTp[0] = (zZTp[0]-zo)
    zZPp[0] = (zZPp[0]-zo)
    
    I = np.trapz((zZDb[1]*zZDt[1]),zZDb[0])

    os.chdir(ocwd)

    return zZDa,zZDt,zZDb,zZVp,zZTp,zZPp,I,D

    
def zD_profs(path,stage):
    ocwd = os.getcwd()
    os.chdir(path)
    zZDa,zZDt,zZDb,zZVp,zZTp,zZPp,zo,D = read_profs(path,stage)
    
    zZDa[0] = (zZDa[0]-zo)/D
    zZDt[0] = (zZDt[0]-zo)/D
    zZDb[0] = (zZDb[0]-zo)/D
    zZVp[0] = (zZVp[0]-zo)/D
    zZTp[0] = (zZTp[0]-zo)/D
    zZPp[0] = (zZPp[0]-zo)/D

    os.chdir(ocwd)

    return zZDa,zZDt,zZDb,zZVp,zZTp,zZPp
    
    
def read_Rg(path):
    ocwd = os.getcwd()
    os.chdir(path)
    txts = filter(lambda x: (x != 'bsmol.txt' and x != 'BSMolf.py' and x != 'FreeChain.in' and x != 'log.lammps') ,os.listdir('.'))
    n = len(txts)
    Rgs = np.zeros(n)
    Rge = np.zeros(n)
    Ns = np.zeros(n)
    i=0
    for txt in txts:
        datas = pd.read_csv(txt, delimiter=' ', comment= '#').values
        data = datas[-1000:,1]
        Rgs[i] = np.mean(data)
        Rge[i] = np.std(data)
        Ns[i] = txt.split('.')[0]
        i=i+1
    A = zip(Ns,Rgs,Rge)    
    B = sorted(A)
    Ns,Rgs,Rge = list(zip(*B))
    Ns = np.asarray(Ns)
    Rgs = np.asarray(Rgs)
    print(Ns)
    print(Rgs)
    Rge = np.asarray(Rge)
    os.chdir(ocwd)
    return Ns,Rgs,Rge
    
    
    
    
    
def read_RDF(txt):
    ocwd = os.getcwd()
    data = np.genfromtxt(txt,comments='#', dtype=float)
    print(data.shape)
    m,n = data.shape
    print(m)
    print(n)
    r = np.zeros(m)
    g =  np.zeros((m,n-1))
    
    for i in range(m):
        r[i] = data[i][0]
        g[i,:] = data[i][1:]
        
    os.chdir(ocwd)
    return r,g
    
    

if __name__ == '__main__':
    
#  0   1  2   3   4    5       6       7           8       9       10   11  12  13  14  15  16   17     18      19         20     21         22      23         24         25          26       27       28      29       30    31    32    33      34         35             36      37           38         39        40          41           
# step et ke pe epair temp c_melTemp c_wallTemp v_Fcatom v_Pcomp2 press pxx pyy pzz pxy pxz pyz c_Vir c_Vir[1] c_Vir[2] c_Vir[3] c_Vir[4] c_Vir[5] c_Vir[6] c_melPress c_wallPress v_melDens v_surfcov v_aveRg v_Vwall v_srate v_D v_bwzmax zhi c_fbwall[1] c_fbwall[3] c_ftwall[1] c_ftwall[3] c_ggbot[1] c_ggbot[3] c_ggtop[1] c_ggtop[3]
    
  
    print('Analysis')
    
#    #---RDF from Ovito---#
#    
#    r,g = read_RDF('beq')
#    
#    print(r)
#    print(g[:,1])
#    
#    plt.rcParams.update({'font.size': 15})
#    fig, ax = plt.subplots(figsize=(10,8))
#    l1, = ax.plot(r , g[:,3], 'g-')
#    plt.legend([l1], [r'$3-3$'], loc=1)
#    ax.set(xlabel='$r$', ylabel= r'$g(r)$' , title='')
#    ax.minorticks_on()
#    plt.savefig('1-RDF.jpg')      
#    
#     
#     #---FreeChain---#
#     Ns,Rgs,Rge = read_Rg(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Isolated Chain')
#     
#     deN = np.linspace(Ns[0],Ns[-1],80)
#     deR = 0.4348*(deN)**0.61
#     
#     plt.rcParams.update({'font.size': 15})
#     fig, ax = plt.subplots(figsize=(10,8))
#     ax.errorbar(Ns ,Rgs,yerr=Rge , fmt = 'ko', capsize=5, label='$Simulation~N^{0.61}$')
#     ax.plot(deN,deR,'k--', label='$deGennes~N^{0.6}$')
#     plt.legend(loc=9)
#     ax.set(xlabel='$N$', ylabel= r'$R_g$' , title='Free Chain Radii of Gryration')
#     ax.minorticks_on()
#     plt.savefig('0-FreeChain.jpg') 
#     
#     
# 

#     #---EQUILIBRIUM---#
# 
#     #---Density Profile---#
#     zZDa,zZDt,zZDb,zZVp,zZTp,zZPp,I,D = zZ_profs(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\Profiles\MDPBB-ECS1\1_V',1)
# 
#     plt.rcParams.update({'font.size': 15})
#     fig, ax = plt.subplots(figsize=(10,8))
#     l1, = ax.plot(zZDb[0] , zZDb[1], 'g--')
#     plt.legend([l1], [r'$\rho_g$'], loc=1)
#     ax.set(xlim=[0,D/3],xlabel='$z$', ylabel= r'$\phi(z)$' , title='')
#     ax.minorticks_on()
#     plt.savefig('1-EquilDP-N-30.jpg')   
#     
#     #---Brush Height---#
#     
#     h = brush_height(zZDb[0],zZDb[1])
#     print(h)
#     
    #---COMPRESSION---#
# 
    #---Compression Curve---#
    
#     #---N=30---#
#     Ps,m30l0,m30l0e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=30\M=30\L=0',[9,31],2)
#     Ps,m60l0,m60l0e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=30\M=60\L=0',[9,31],2)
#     Ps,m90l0,m90l0e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=30\M=90\L=0',[9,31],2)
#     
#     Ps,m30l1,m30l1e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=30\M=30\L=1',[9,31],2)
#     Ps,m60l1,m60l1e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=30\M=60\L=1',[9,31],2)
#     Ps,m90l1,m90l1e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=30\M=90\L=1',[9,31],2)
# 
# 
#     plt.rcParams.update({'font.size': 15})
#     fig, ax = plt.subplots(figsize=(10,8))
#     ax.set_xscale("log")
#     #ax.set_yscale("log")
#     ax.errorbar(m30l0[:,0], m30l0[:,1]/Rgs[1], yerr=m30l0e[:,1], fmt='r--o',capsize=5)
#     ax.errorbar(m30l1[:,0], m30l1[:,1]/Rgs[1], yerr=m30l1e[:,1], fmt='r-o',capsize=5, label=r'$\rho_g/{\rho_g}^* = 1.07$')
#     ax.errorbar(m60l0[:,0], m60l0[:,1]/Rgs[1], yerr=m60l0e[:,1], fmt='g--o',capsize=5)
#     ax.errorbar(m60l1[:,0], m60l1[:,1]/Rgs[1], yerr=m60l1e[:,1], fmt='g-o',capsize=5, label=r'$\rho_g/{\rho_g}^* = 2.14$')
#     ax.errorbar(m90l0[:,0], m90l0[:,1]/Rgs[1], yerr=m90l0e[:,1], fmt='b--o',capsize=5)
#     ax.errorbar(m90l1[:,0], m90l1[:,1]/Rgs[1], yerr=m90l1e[:,1], fmt='b-o',capsize=5, label=r'$\rho_g/{\rho_g}^* = 3.21$')
# 
#     
#     plt.legend( loc=1)
#     ax.set(xlim=[0,3.5],xlabel= r'$log(P)$', ylabel= '$D/{R_g}$', title='N=30')
#     ax.minorticks_on()
#     plt.savefig('3-CompCurve-N=30.jpg')
# 
#     Ps,m30l0,m30l0e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=90\M=30\L=0',[9,31],2)
#     Ps,m60l0,m60l0e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=90\M=60\L=0',[9,31],2)
#     Ps,m90l0,m90l0e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=90\M=90\L=0',[9,31],2)
#     
#     Ps,m30l1,m30l1e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=90\M=30\L=1',[9,31],2)
#     Ps,m60l1,m60l1e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=90\M=60\L=1',[9,31],2)
#     Ps,m90l1,m90l1e = read_thermoP(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\CompSims\Sims2\N=90\M=90\L=1',[9,31],2)
# 
# 
#     plt.rcParams.update({'font.size': 15})
#     fig, ax = plt.subplots(figsize=(10,8))
#     ax.set_xscale("log")
#     #ax.set_yscale("log")
#     ax.errorbar(m30l0[:,0], m30l0[:,1]/Rgs[1], yerr=m30l0e[:,1], fmt='r--o',capsize=5)
#     ax.errorbar(m30l1[:,0], m30l1[:,1]/Rgs[1], yerr=m30l1e[:,1], fmt='r-o',capsize=5, label=r'$\rho_g/{\rho_g}^* = 4.8$')
#     ax.errorbar(m60l0[:,0], m60l0[:,1]/Rgs[1], yerr=m60l0e[:,1], fmt='g--o',capsize=5)
#     ax.errorbar(m60l1[:,0], m60l1[:,1]/Rgs[1], yerr=m60l1e[:,1], fmt='g-o',capsize=5, label=r'$\rho_g/{\rho_g}^* = 9.6$')
#     ax.errorbar(m90l0[:,0], m90l0[:,1]/Rgs[1], yerr=m90l0e[:,1], fmt='b--o',capsize=5)
#     ax.errorbar(m90l1[:,0], m90l1[:,1]/Rgs[1], yerr=m90l1e[:,1], fmt='b-o',capsize=5, label=r'$\rho_g/{\rho_g}^* = 14.5$')
# 
#     
#     plt.legend( loc=1)
#     ax.set(xlim=[0,3.5],xlabel= r'$log(P)$', ylabel= '$D/{R_g}$', title='N=90')
#     ax.minorticks_on()
#     plt.savefig('3-CompCurve-N=90.jpg')
#     
#     
#     

# 
#     #---Interpenetration Profiles---#
#     
#     zDDa,zDDt,zDDb,zDVp,zDTp,zDPp = zD_profs(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\Profiles\MDPBB-ECS1\1_V',3)
# 
#     fig, ax = plt.subplots(figsize=(10,8))
#     lb1, = ax.plot(zDDb[0] , zDDb[1], 'r--')
#     lt1, = ax.plot(zDDt[0] , zDDt[1], 'b--')
#     la1, = ax.plot(zDDa[0] , zDDa[1], 'g')
#     
# #    plt.legend([la1], ['P=1'], loc=1)
#     ax.set(xlim=[0,1],xlabel='$z/D$', ylabel= r'$\phi(z)$' , title='')
#     ax.minorticks_on()
#     plt.savefig('3-Interpenetration.jpg')
#     
#     # #---Interpenetration vs D---#
#     # 
#     # 
#     # Ps,Is,Ds = IvsD(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\FiniteSize\GD-0.075\MDPBB-FS-1')
#     # 
#     # fig, ax = plt.subplots(figsize=(10,8))
#     # l1, = ax.plot(Ds,Is, 'go')
#     # 
#     # plt.legend([l1],[r'$\rho_g$'], loc=9)
#     # ax.set(xlabel='$D$', ylabel= r'$I(D)$' , title='')
#     # ax.minorticks_on()
#     # plt.savefig('4-IvsD.jpg')
#     # 
#     # 
#     
#     #---SHEARING---#
#     
#     #---T,V and Density profiles---#
#     
#     zDDa,zDDt,zDDb,zDVp,zDTp,zDPp = zD_profs(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\Profiles\MDPBB-ECS1\1_V',3)
#     
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     ax2 = ax1.twinx()   
#     T1, = ax1.plot(zDTp[0],zDTp[1], 'r-')
#     D1, = ax1.plot(zDDa[0],zDDa[1], 'g-')
#     V1, = ax2.plot(zDVp[0],zDVp[1], 'b-')
#     
#     plt.legend([T1,V1,D1], [r'Temperature','Velocity','Density'], loc=4)
#     ax1.set(xlim=[0,1],ylim=[0,1.5],xlabel='$z/D$', ylabel= r'$T, \phi(z)$')
#     ax2.set(ylim=[-0.25,0.25], ylabel= r'V')
#     ax1.minorticks_on()
#     ax2.minorticks_on()
#     plt.savefig('4-TnVprofiles.jpg') 
#     
#     #---Stress Profiles---#
#     
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     ax2 = ax1.twinx()   
#     Pzx, = ax1.plot(zDPp[0],zDPp[1][:,0], 'r-')
#     Pzz, = ax1.plot(zDPp[0],zDPp[1][:,2], 'b-')
#     D, = ax2.plot(zDDa[0],zDDa[1], 'g-')
# 
#     plt.legend([Pzx,Pzz,D], [r'$P_(zx)$',r'$P_(zz)$',r'$\phi(z)$'], loc=4)
#     ax1.set(xlim=[0,1],ylim=[0,5],xlabel='$z/D$', ylabel= r'$\phi(z)$')
#     ax2.set(ylim=[0,1.5], ylabel= r'$\phi(z)$')
#     ax1.minorticks_on()
#     ax2.minorticks_on()
#     plt.savefig('5-PnDprofiles.jpg')
#     
#     #---Time Evolution---#
#     
#     
#     ts,Vt = thermo_time(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\FiniteSize\GD-0.075\MDPBB-FS-1\1_P',[5,24],3)
# 
#     fig, ax1 = plt.subplots(figsize=(10,8))
#     l1, = ax1.plot(ts[-1000:],Vt[-1000:,0], 'ro')
#     ax2 = ax1.twinx()
#     l2, = ax2.plot(ts[-1000:],Vt[-1000:,1], 'bo')
#     ax1.minorticks_on()
#     ax2.minorticks_on()
#     plt.savefig('6-Time Evolution.jpg') 
#     
#     
  #   #---Shear and Normal stress components vs srate---#
  #   
  #   Vs1,Cs1,Ce1 = read_thermoV(r'C:\Users\klay_\OneDrive - Imperial College London\PhD\PhD-Sims\Archive\FiniteSize\P-1\MDPBB-FS-1',[13,15,23],3)    
  #   
  #   fig, ax1 = plt.subplots(figsize=(10,8))
  #   ax2 = ax1.twinx()
  #   n1, = ax1.plot(Cs1[:,2],Cs1[:,0], 'ro')
  #   s1, = ax2.plot(Cs1[:,2],Cs1[:,1], 'bo')
  #   
  #   plt.legend([n1,s1], ['$P_{zz}$','$P_{xz}$'], loc=9)
  #   ax1.set(xlabel=r'$\.\gamma$', ylabel= '$P_{zz}$')
  #   ax2.set(ylabel = '$P_{xz}$')
  #   ax1.minorticks_on()
  #   ax2.minorticks_on()
  #   plt.savefig('7-PzzPxz vs srate - vF.jpg')
  #   
  #   
  #   mu1 = -1*Cs1[:,1]/Cs1[:,0]
  # 
  #   fig, ax1 = plt.subplots(figsize=(10,8))
  #   l1, = ax1.plot(Cs1[:,2],mu1, 'bo')
  #   
  #   plt.legend([l1], [r'f=1',], loc=1)     
  #   ax1.set(xlabel=r'$\.\gamma$', ylabel= r'$\mu$', xlim=[0,0.3], ylim=[0,0.45])
  #   ax1.minorticks_on()
  #   plt.savefig('8-COF vs srate - vF.jpg')
  #   
