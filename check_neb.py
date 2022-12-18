import os, sys
from glob import glob
from ase.io import *
import numpy as np
from  ase.io.trajectory import TrajectoryReader
import matplotlib.pyplot as plt

def vasp_energy(dir):
    init_dir=os.getcwd()
    os.chdir(dir)
    with open('OSZICAR') as infile:
        lines = infile.readlines()

    final = len(lines)-1
    for i in range(final,0,-1):
        if ' F= ' in lines[i]:
            energy=float(lines[i].split()[4])
            break
    os.chdir(init_dir)
    return(energy)

calcs = [] 
met_path = '0*/'
calcs+=sorted(glob(met_path))
energies = np.zeros(len(calcs))
images = np.zeros(len(calcs))
n=0
print('Energies (eV): ')
for calc in calcs:
    calcpath = os.path.dirname(calc)
    output = glob(os.path.join(calcpath,'OUTCAR'))
    energies[n] = vasp_energy(calcpath)
    images[n] = n
    print(energies[n])
    n+=1

print('\n'+'E_fwd: %.2f eV'%(np.amax(energies)-energies[0]))
print('E_rev: %.2f eV'%(np.amax(energies)-energies[-1]))
print('Areac: %.2f eV'%(energies[-1]-energies[0]))

energies -= energies[0]

plt.plot(images[:], energies[:])
plt.ylabel('Energy (eV)')
plt.show()
