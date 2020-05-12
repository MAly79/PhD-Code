#!/usr/bin/python
#$: chmod 755 yourfile.py
#$: dos2unix yourfile.py
#$: ./yourfile.py

import os
import shutil
import numpy as np

def make_dirs(vals,name):
    for v in vals:
        folderName = name + str(v)
        if os.path.isdir(folderName) == False:
            os.mkdir(folderName)
    return None



def mod_main(l,line):
    f=open('main.in','r')
    lines = f.readlines()
    f.close()
    lines[l] = line
    f2 = open('main.in','w')
    f2.writelines(lines)
    f2.close()

    return None

def myRg(N):
    Rout = 0.4348*(N)**0.61
    return Rout

def myrho_c(Rg):
    return (1)/((np.pi)*Rg**2)




if __name__ == '__main__':

    # To run the sims on Home PC go to any L directory and run
    # for d in ./P*/LPBB-ECS; do (cd "d" && lmp -in main.in);
    # To run on cx1 I imagine something like
    # for d in ./N*/M*/L*/P*/LPBB-ECS; do (cd "d" && qsub lmpJS.pbs); done;

    owd = os.getcwd()
    seriesName= 'Sims2'
    if os.path.isdir(seriesName+r'/LPBB-ECS') == False:
        shutil.copytree('LPBB-ECS',seriesName+r'/LPBB-ECS')
    os.chdir(seriesName)

    Ns = [30,60,90]
    GDs = [0.5, 1, 1.5, 2, 4]
    Ms = [30, 60, 90]
    Ls = [0,0.25,0.75,1]
    Ps = [0.0001,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.6,3.2]



    make_dirs(Ns,'N=')
    Ndirs = os.listdir('.')
    print(os.getcwd())

    for Ndir in Ndirs:
        if Ndir != 'LPBB-ECS':
            if os.path.isdir(Ndir+r'/LPBB-ECS') == False:
                shutil.copytree('LPBB-ECS',Ndir+r'/LPBB-ECS')
            os.chdir(Ndir)
            make_dirs(GDs,'GD=')
            os.chdir('LPBB-ECS')
            vN = Ndir.split("=")[1]
            mod_main(3, 'variable      N      equal  '+vN+'\n')
            os.chdir('..')
            Gdirs = os.listdir('.')

            ################################

            for Gdir in Gdirs:
                if Gdir != 'LPBB-ECS':
                    if os.path.isdir(Gdir+r'/LPBB-ECS') == False:
                        shutil.copytree('LPBB-ECS',Gdir+r'/LPBB-ECS')
                    os.chdir(Gdir)
                    make_dirs(Ls,'L=')

                    os.chdir('LPBB-ECS')
                    vGD = Gdir.split("=")[1]
                    Rg = myRg(float(vN))
                    rho_c = myrho_c(Rg)
                    GD_t = float(vGD) * rho_c
                    A = 30*30
                    M_tar = A * GD_t
                    M_tar = int(round(M_tar/2))
                    vM = str(M_tar*2)
                    mod_main(2, 'variable      M      equal  ' +vM+'\n')
                    os.chdir('..')

                    Ldirs = os.listdir('.')
                    print(Ldirs)
                    #####################################################
                    for Ldir in Ldirs:
                        if Ldir != 'LPBB-ECS':
                            if os.path.isdir(Ldir+r'/LPBB-ECS') == False:
                                shutil.copytree('LPBB-ECS',Ldir+r'/LPBB-ECS')
                            os.chdir(Ldir)
                            make_dirs(Ps,'P=')
                            os.chdir('LPBB-ECS')
                            v = Ldir.split("=")[1]
                            mod_main(4, 'variable      L      equal  ' +v+ '\n')
                            os.chdir('..')
                            Pdirs = os.listdir('.')
                            print(Pdirs)

                            ###############################################
                            for Pdir in Pdirs:
                                if Pdir != 'LPBB-ECS':
                                    if os.path.isdir(Pdir+r'/LPBB-ECS') == False:
                                        shutil.copytree('LPBB-ECS',Pdir+r'/LPBB-ECS')

                                    os.chdir(Pdir)

                                    os.chdir('LPBB-ECS')
                                    v = Pdir.split("=")[1]
                                    mod_main(20, 'variable       Pcomp      index  ' +v+ '\n')
                                    os.chdir('..')

                                    os.chdir('..')

                            ###############################################
                            os.chdir('..')


                    #####################################################
                    os.chdir('..')



            #################################
            os.chdir('..')


    os.chdir(owd)
