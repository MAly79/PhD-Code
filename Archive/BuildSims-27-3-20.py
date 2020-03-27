#!/usr/bin/python
#$: chmod 755 yourfile.py
#$: dos2unix yourfile.py
#$: ./yourfile.py

import os
import shutil

def make_dirs(vals,name):
    for v in vals:
        folderName = name + str(v)
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




if __name__ == '__main__':

    # To run the sims on Home PC go to any L directory and run
    # for d in ./P*/LPBB-ECS; do (cd "d" && lmp -in main.in);
    # To run on cx1 I imagine something like
    # for d in ./N*/M*/L*/P*/LPBB-ECS; do (cd "d" && qsub lmpJS.pbs); done;

    owd = os.getcwd()
    seriesName= 'Sims'
    shutil.copytree('LPBB-ECS',seriesName+r'/LPBB-ECS')
    os.chdir(seriesName)
    Ns = [30, 90]
    Ms = [30, 60, 90]
    Ls = [0,1]
    Ps = [0.0001,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.6,3.2,6.4]



    make_dirs(Ns,'N=')
    Ndirs = os.listdir('.')
    print(os.getcwd())

    for Ndir in Ndirs:
        if Ndir != 'LPBB-ECS':
            shutil.copytree('LPBB-ECS',Ndir+r'/LPBB-ECS')

            os.chdir(Ndir)
            make_dirs(Ms,'M=')

            os.chdir('LPBB-ECS')
            v = Ndir.split("=")[1]
            mod_main(3, 'variable      N      equal  '+v+'\n')
            os.chdir('..')

            Mdirs = os.listdir('.')

            ################################

            for Mdir in Mdirs:
                if Mdir != 'LPBB-ECS':
                    shutil.copytree('LPBB-ECS',Mdir+r'/LPBB-ECS')
                    os.chdir(Mdir)
                    make_dirs(Ls,'L=')

                    os.chdir('LPBB-ECS')
                    v = Mdir.split("=")[1]
                    mod_main(2, 'variable      M      equal  ' +v+'\n')
                    os.chdir('..')

                    Ldirs = os.listdir('.')
                    print(Ldirs)
                    #####################################################
                    for Ldir in Ldirs:
                        if Ldir != 'LPBB-ECS':
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
