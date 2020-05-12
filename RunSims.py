#!/usr/bin/python
#$: chmod 755 yourfile.py
#$: dos2unix RunSims.py
#$: ./yourfile.py


import os
import shutil
import glob

def make_dirs(vals,name):
    for v in vals:
        folderName = name + str(v)
        os.mkdir(folderName)

    return None

def run_sims(b1,b2):
    owd = os.getcwd()

    os.chdir('Sims2')

    Ndirs = os.listdir('.')
    i = 1

    for Ndir in Ndirs:
        if Ndir != 'LPBB-ECS':
            os.chdir(Ndir)
            Mdirs = os.listdir('.')
            for Mdir in Mdirs:
                if Mdir != 'LPBB-ECS':
                    os.chdir(Mdir)
                    Ldirs = os.listdir('.')
                    for Ldir in Ldirs:
                        if Ldir != 'LPBB-ECS':
                            os.chdir(Ldir)
                            Pdirs = os.listdir('.')
                            for Pdir in Pdirs:
                                if Pdir != 'LPBB-ECS':
                                    os.chdir(Pdir)
                                    os.chdir('LPBB-ECS')
                                    status = 'E'
                                    #print(len(glob.glob('*JS.pbs.o*')))
                                    if len(glob.glob('*JS.pbs.o*')) != 0:
                                        status = 'R'
                                    if os.path.exists('shear.csv') == True:
                                        status = 'C'
                                    if i <= b2 and i >= b1 and status == 'E':
                                        status = 'Q'
                                        #print(os.getcwd())
                                        os.system('chmod 755 BSMolf.py')
                                        os.system('qsub serJS.pbs')
                                        #f= open("serJS.pbs.o982389","w+")
                                        #f.close
                                    print([i,os.getcwd(),status])
                                    i = i + 1
                                    os.chdir('..')
                                    os.chdir('..')
                            os.chdir('..')
                    os.chdir('..')
            os.chdir('..')



    os.chdir(owd)




    return None



if __name__ == '__main__':

    # To run the sims on Home PC go to any L directory and run
    # for d in ./P*/LPBB-ECS; do (cd "d" && lmp -in main.in);
    # To run on cx1 I imagine something like
    # for d in ./N*/M*/L*/P*/LPBB-ECS; do (cd "d" && qsub lmpJS.pbs);
    run_sims(201,300)
