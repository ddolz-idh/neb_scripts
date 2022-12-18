import os
import math

main_path = '/mnt/c/Users/PC/Desktop/Doctoral/NewPhases'  # update
potcars_path = '/mnt/c/Users/PC/Desktop/Doctoral/potcars'   # update

max_phases = {
    'Ti2AlC': 12990,
    'Zr2AlC': 3886,
    'Hf2InC': 22156,
    'V2AlC': 1025497,
    'Nb2AlC': 996162,
    'Ta2AlC': 1025441,
    'Cr2AlC': 9956,
    'Mo2GaC': 1079635,
    'W2AlC': 1079892,
    'Ti2AlN': 4978,
    'Zr2AlN': 4678,
    'Hf2AlN': 'missing',
    'V2AlN': 'missing',
    'Nb2AlN': 'missing',
    'Ta2AlN': 'missing',
    'Cr2AlN': 10371,
    'Mo2AlN': 'missing',
    'W2AlN': 'missing',
}

incar_tags_bulk = {
    'IBRION': 2,
    'ISIF': 3,
    'EDIFF': 1E-06,
    'EDIFFG': -0.001,
    'NSW': 500,
    'NELM': 300,
    'ISMEAR': -5,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'ENCUT': 520,
    'GGA': 'PE',
    'IVDW': 11,
    'LASPH': '.TRUE.',
    'ISPIN': 2
}

incar_tags_relaxation = {
    'IBRION': 2,
    'EDIFF': 1E-06,
    'EDIFFG': -0.01,
    'NSW': 500,
    'NELM': 300,
    'ISMEAR': 1,
    'LWAVE': '.FALSE.',
    'LCHARG' : 'FALSE.',
    'ENCUT': 415,
    'ALGO': 'Fast',
    'LREAL': 'Auto',
    'GGA': 'PE',
    'IVDW': 11,
    'LDIPOL': '.TRUE.',
    'IDIPOL': '3',
    'DIPOL': '0.5 0.5 0.5',
    'LASPH': '.TRUE.',
    'ISPIN': 2,
    'NPAR': 4  # update
}

incar_tags_nebs = {
    'IBRION': 2,
    'EDIFF': 1E-06,
    'EDIFFG': -0.01,
    'NSW': 500,
    'NELM': 300,
    'ISMEAR': 1,
    'LWAVE': '.FALSE.',
    'LCHARG' : 'FALSE.',
    'ENCUT': 415,
    'ALGO': 'Fast',
    'LREAL': 'Auto',
    'GGA': 'PE',
    'IVDW': 11,
    'LDIPOL': '.TRUE.',
    'IDIPOL': '3',
    'DIPOL': '0.5 0.5 0.5',
    'LASPH': '.TRUE.',
    'ISPIN' : 1,
    'NPAR' : 4,
    'ICHAIN' : 0,
    'IMAGES' : 8,
    'LCLIMB' : '.TRUE.',
    'SPRING' : -5,
}

incar_tags_frequency = {
    'IBRION': 5,
    'EDIFF': 1E-06,
    'EDIFFG': -0.01,
    'POTIM' : 0.03, 
    'NSW': 500,
    'NELM': 300,
    'ISMEAR': 1,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'ENCUT': 415,
    'ALGO': 'Fast',
    'LREAL': 'Auto',
    'GGA': 'PE',
    'IVDW': 11,
    'LDIPOL': '.TRUE.',
    'IDIPOL': '3',
    'DIPOL': '0.5 0.5 0.5',
    'LASPH': '.TRUE.',
    'ISPIN': 2 # update
}

incar_tags_frequencyrepeat = {
    'IBRION': 5,
    'EDIFF': 1E-06,
    'EDIFFG': -0.01,
    'POTIM' : 0.03, 
    'NSW': 500,
    'NELM': 300,
    'ISMEAR': 1,
    'LWAVE': '.FALSE.',
    'ICHARG': 1,
    'ENCUT': 415,
    'ALGO': 'Fast',
    'LREAL': 'Auto',
    'GGA': 'PE',
    'IVDW': 11,
    'LDIPOL': '.TRUE.',
    'IDIPOL': '3',
    'DIPOL': '0.5 0.5 0.5',
    'LASPH': '.TRUE.',
    'ISPIN': 2 # update
}

incar_tags_gas = {
    'IBRION': 2,
    'EDIFF': 1E-06,
    'EDIFFG': -0.01,
    'NSW': 500,
    'NELM': 300,
    'ISMEAR': 0,
    'SIGMA': 0.01,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'ENCUT': 415,
    'LREAL': 'Auto',
    'GGA': 'PE',
    'IVDW': 11,
    'LASPH': '.TRUE.',
    'ISPIN': 2
}


def create_incar(path, dict_tags):
    f = open(f"{path}/INCAR", 'w')
    for tag in dict_tags:
        f.write(f"{tag} = {dict_tags[tag]}\n")
    f.close()

def correct_positions(atoms):
    a = atoms.cell.cellpar()[0]
    b = atoms.cell.cellpar()[1]
    alpha = b = atoms.cell.cellpar()[-1]*math.pi/180
    for i in range(len(atoms)):
        if atoms.positions[i][0] >= a + atoms.positions[i][1]/math.tan(alpha)-0.5:
            atoms.positions[i][0] -= a
        if atoms.positions[i][1] >= b*math.sin(alpha)-0.5:
            atoms.positions[i][0] -= b*math.cos(alpha)
            atoms.positions[i][1] -= b*math.sin(alpha)

def create_potcar(path):
    with open(f"{path}/POSCAR") as infile:
        lines = infile.readlines()
    elements = lines[0].strip().split()
    for i in range(len(elements)):
        elements[i] = 'POTCAR-' + elements[i]
    with open(f"{path}/POTCAR", 'w') as outfile:
        for element in elements:
            with open(f"{potcars_path}/{element}") as infile:
                for line in infile:
                    outfile.write(line)


def create_kpoints(path, kpoints):
    f = open(f"{path}/KPOINTS", 'w')
    f.write('Automatic mesh\n')
    f.write('0\n')
    f.write('Gamma\n')
    f.write(f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n")
    f.write('0 0 0\n')
    f.close()


def create_submit_script_MN4(path, num_nodes, walltime_limit):
    f = open(f"{path}/vasp_submission", 'w')
    f.write('#!/bin/bash -l\n')
    f.write('#BATCH --job-name=VASP_test\n')
    f.write('#SBATCH --output=vasp_%j.out\n')
    f.write('#SBATCH --error=vasp_%j.err\n')
    f.write('#SBATCH --ntasks=96\n')
    f.write(f"#SBATCH --time={walltime_limit}:00:00\n")
    f.write('#SBATCH -D .\n')
    f.write('\n')
    f.write('\n')
    f.write('export DIR=$PWD\n')
    f.write('\n')
    f.write('mkdir -p /gpfs/scratch/ub108/ub108015/${DIR}\n')
    f.write('cp * /gpfs/scratch/ub108/ub108015/${DIR}\n')
    f.write('cd /gpfs/scratch/ub108/ub108015/${DIR}\n')
    f.write('\n')
    f.write('module purge\n')
    f.write('module load intel impi mkl vasp/5.4.4\n')
    f.write('mpirun vasp_std\n')
    f.write('cp CONTCAR OUTCAR OSZICAR PCDAT XDATCAR vasp* ${DIR}\n')
    f.write('cd ${DIR}\n')
    f.close()


