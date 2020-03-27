#!/usr/bin/python
#$: chmod 755 yourfile.py
#$: dos2unix RunSims.py
#$: ./yourfile.py


import os
import shutil

def make_dirs(vals,name):
    for v in vals:
        folderName = name + str(v)
        os.mkdir(folderName)

    return None



if __name__ == '__main__':

    # To run the sims on Home PC go to any L directory and run
    # for d in ./P*/LPBB-ECS; do (cd "d" && lmp -in main.in);
    # To run on cx1 I imagine something like
    # for d in ./N*/M*/L*/P*/LPBB-ECS; do (cd "d" && qsub lmpJS.pbs);

    owd = os.getcwd()

    os.chdir('Sims')

    Ndirs = os.listdir('.')

    print(Ndirs)

    for Ndir in Ndirs[1]:
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
                                    #os.system('qsub lmpJS.pbs')
                                    print(os.getcwd)
                                    os.chdir('..')
                                    os.chdir('..')
                            os.chdir('..')
                    os.chdir('..')
            os.chdir('..')

    os.chdir(owd)
