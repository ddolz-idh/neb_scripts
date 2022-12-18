import sys
import numpy as np
import shutil

from tools import *
from ase.io import read, write
from ase.build import add_adsorbate
import ase
import ase.neb

sites = {
    'topM': {'atom ids': [18], 'height': 2.2},
    'hX': {'atom ids': [10, 12, 18], 'height': 2.0},
    'br': {'atom ids': [12, 18], 'height': 2.0},
    'hM': {'atom ids': [17], 'height': 2.0},
}
sites2 = {
    'topM': {'atom ids': [24], 'height': 2.2},
    'hX': {'atom ids': [16, 18, 24], 'height': 2.0},
    'br': {'atom ids': [18, 24], 'height': 2.0},
    'hM': {'atom ids': [15], 'height': 2.0}
}

def get_site_coordinates(atoms, atom_ids):
    if len(atom_ids) == 1:
        x = atoms[atom_ids[0]].x
        y = atoms[atom_ids[0]].y
    elif len(atom_ids) == 2:
        atom1 = atoms[atom_ids[0]].position
        atom2 = atoms[atom_ids[1]].position
        x = (atom1[0] + atom2[0])/2
        y = (atom1[1] + atom2[1])/2
    elif len(atom_ids) == 3:
        atom1 = atoms[atom_ids[0]].position
        atom2 = atoms[atom_ids[1]].position
        atom3 = atoms[atom_ids[2]].position
        x, y, z = get_centroid(atom1, atom2, atom3)
    else:
        sys.exit('Invalid list of IDs')
    return np.array((x, y))


def get_centroid(p1, p2, p3):
    # Calculate middle point between 3 coordinates
    x = (p1[0] + p2[0] + p3[0])/3
    y = (p1[1] + p2[1] + p3[1])/3
    z = (p1[2] + p2[2] + p3[2])/3
    return x, y, z

for mxene in ['W2N']:

    if not os.path.exists(f"{main_path}/{mxene}/adsorption"):

        print('Creating adsorption inputs')
        os.mkdir(f"{main_path}/{mxene}/adsorption")

        for adsorbate in [ 'CO','C', 'O','H2O','H2','H','OH','CH','CH2','CH3','CH4']:

            for site in sites:

                path1 = f"{main_path}/{mxene}/adsorption/{mxene}#{adsorbate}-{site}"

                if not os.path.exists(path1):
                    os.mkdir(path1)
                    atoms_slab = read(f"{main_path}/{mxene}/CONTCAR")
                    atoms_adsorbate = read(f"{main_path}/adsorbates/POSCAR-{adsorbate}")
                    xy_coords = get_site_coordinates(atoms=atoms_slab, atom_ids=sites[site]['atom ids'])
                    add_adsorbate(slab=atoms_slab,
                                  adsorbate=atoms_adsorbate,
                                  height=sites[site]['height'],
                                  position=(xy_coords[0], xy_coords[1]))
                    atoms_slab.write(filename=f"{path1}/POSCAR", vasp5=True)
                    create_potcar(path=f"{path1}")
                    incar_tags = incar_tags_relaxation.copy()
                    incar_tags['ISPIN'] = 2
                    create_incar(path=f"{path1}", dict_tags=incar_tags)
                    create_kpoints(path=f"{path1}", kpoints=(5, 5, 1))
                    create_submit_script(path=f"{path1}", num_nodes=1, walltime_limit=5)
                else:
                    print(f'Directory {path} already exists.')

                quit()


    if not os.path.exists(f"{main_path}/{mxene}/frequency"):

        print('Creating frequency inputs')
        os.mkdir(f"{main_path}/{mxene}/frequency")

        for adsorbate in [ 'CO','C', 'O','H2O','H2','H','OH','CH','CH2','CH3','CH4','CO2-X']:

            for site in sites:

                path = f"{main_path}/{mxene}/frequency/{mxene}#{adsorbate}-{site}"

                if not os.path.exists(path):
                    os.mkdir(path)
                    atoms_slab = read(f"{main_path}/adsorption/{mxene}#{adsorbate}-{site}/CONTCAR")
                    atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)

                    afile = open(f"{path}/POSCAR","r")
                    listlines = afile.readlines()
                    listlinespi = afile.readlines(7)
                    
                    listlines[7] = 'Selective Dynamics' + '\n' +'Cartesian' + '\n'
                    for line in range(8, 35):
                        listlines[line] = listlines[line].rstrip('\n') + " F F F"+'\n'
                    if adsorbate in ['CO','H2','OH','CH']:                       
                        for line in range(35, 37):
                            listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                    if adsorbate in ['C','O','H']:                       
                        listlines[35] = listlines[35].rstrip('\n') + " T T T" +'\n'
                    if adsorbate in ['H2O','CH2','CO2-X']:                       
                        for line in range(35, 38):
                            listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                    if adsorbate in ['CH3']:                       
                        for line in range(35, 39):
                            listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                    if adsorbate in ['CH4']:                       
                        for line in range(35, 40):
                            listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                    afile = open(f"{path}/POSCAR", "w")
                    afile.writelines(listlinespi)
                    afile.writelines(listlines)
                    afile.close()

                    create_potcar(path=f"{path}")
                    incar_tags = incar_tags_frequency.copy()
                    incar_tags['ISPIN'] = 2
                    create_incar(path=f"{path}", dict_tags=incar_tags)
                    create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                    create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
                else:
                    print(f'Directory {path} already exists.')
                
        quit()

    if not os.path.exists(f"{main_path}/{mxene}/adsorptionboth"):

        print('Creating dual adsorption inputs')
        os.mkdir(f"{main_path}/{mxene}/adsorptionboth")

        for adsorbate in [ 'C']:
            for adsorbates in ['O','H']:
                for site in sites:
                    for sit in sites2:
                        path = f"{main_path}/{mxene}/adsorptionboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}"
                        # Create folder for each MAX phase
                        if not os.path.exists(path):
                            os.mkdir(path)
                            atoms_slab = read(f"{main_path}/{mxene}/CONTCAR")
                            atoms_adsorbate = read(f"{main_path}/adsorbates/POSCAR-{adsorbate}")
                            atoms_adsorbates = read(f"{main_path}/adsorbates/POSCAR-{adsorbates}")
                            xy_coords = get_site_coordinates(atoms=atoms_slab, atom_ids=sites[site]['atom ids'])
                            xy_coords2 = get_site_coordinates(atoms=atoms_slab, atom_ids=sites2[sit]['atom ids'])
                            add_adsorbate(slab=atoms_slab,
                                        adsorbate=atoms_adsorbate,
                                        height=sites[site]['height'],
                                        position=(xy_coords[0], xy_coords[1]))
                            add_adsorbate(slab=atoms_slab,
                                    adsorbate=atoms_adsorbates,
                                    height=sites2[sit]['height'],
                                    position=(xy_coords2[0], xy_coords2[1]))
                            atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)
                            create_potcar(path=f"{path}")
                            incar_tags = incar_tags_relaxation.copy()
                            incar_tags['ISPIN'] = 1
                            create_incar(path=f"{path}", dict_tags=incar_tags)
                            create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                            create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
                        else:
                            print(f'Directory {path} already exists.')


        for adsorbate in [ 'H']:
            for adsorbates in ['H','CH','CH2','CH3','O','OH']:
                for site in sites:
                    for sit in sites2:
                        path = f"{main_path}/{mxene}/adsorptionboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}"
                        # Create folder for each MAX phase
                        if not os.path.exists(path):
                            os.mkdir(path)
                            atoms_slab = read(f"{main_path}/{mxene}/CONTCAR")
                            atoms_adsorbate = read(f"{main_path}/adsorbates/POSCAR-{adsorbate}")
                            atoms_adsorbates = read(f"{main_path}/adsorbates/POSCAR-{adsorbates}")
                            xy_coords = get_site_coordinates(atoms=atoms_slab, atom_ids=sites[site]['atom ids'])
                            xy_coords2 = get_site_coordinates(atoms=atoms_slab, atom_ids=sites2[sit]['atom ids'])
                            add_adsorbate(slab=atoms_slab,
                                        adsorbate=atoms_adsorbate,
                                        height=sites[site]['height'],
                                        position=(xy_coords[0], xy_coords[1]))
                            add_adsorbate(slab=atoms_slab,
                                    adsorbate=atoms_adsorbates,
                                    height=sites2[sit]['height'],
                                    position=(xy_coords2[0], xy_coords2[1]))
                            atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)
                            create_potcar(path=f"{path}")
                            incar_tags = incar_tags_relaxation.copy()
                            incar_tags['ISPIN'] = 1
                            create_incar(path=f"{path}", dict_tags=incar_tags)
                            create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                            create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
                        else:
                            print(f'Directory {path} already exists.')
                            
        for adsorbate in [ 'CO']:
            for adsorbates in ['O']:
                for site in sites:
                    for sit in sites2:
                        path = f"{main_path}/{mxene}/adsorptionboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}"
                        # Create folder for each MAX phase
                        if not os.path.exists(path):
                            os.mkdir(path)
                            atoms_slab = read(f"{main_path}/{mxene}/CONTCAR")
                            atoms_adsorbate = read(f"{main_path}/adsorbates/POSCAR-{adsorbate}")
                            atoms_adsorbates = read(f"{main_path}/adsorbates/POSCAR-{adsorbates}")
                            xy_coords = get_site_coordinates(atoms=atoms_slab, atom_ids=sites[site]['atom ids'])
                            xy_coords2 = get_site_coordinates(atoms=atoms_slab, atom_ids=sites2[sit]['atom ids'])
                            add_adsorbate(slab=atoms_slab,
                                        adsorbate=atoms_adsorbate,
                                        height=sites[site]['height'],
                                        position=(xy_coords[0], xy_coords[1]))
                            add_adsorbate(slab=atoms_slab,
                                    adsorbate=atoms_adsorbates,
                                    height=sites2[sit]['height'],
                                    position=(xy_coords2[0], xy_coords2[1]))
                            atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)
                            create_potcar(path=f"{path}")
                            incar_tags = incar_tags_relaxation.copy()
                            incar_tags['ISPIN'] = 1
                            create_incar(path=f"{path}", dict_tags=incar_tags)
                            create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                            create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
                        else:
                            print(f'Directory {path} already exists.')
        
        quit()

        
    if not os.path.exists(f"{main_path}/{mxene}/frequencyboth"):
        print('Creating dual frequency inputs')
        os.mkdir(f"{main_path}/{mxene}/frequencyboth")

        for adsorbate in [ 'C']:
            for adsorbates in ['O','H']:
                for site in sites:
                    for sit in sites2:
                        path = f"{main_path}/{mxene}/frequencyboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}"
                        # Create folder for each MAX phase
                        if not os.path.exists(path):
                            os.mkdir(path)
                            atoms_slab = read(f"{main_path}/{mxene}/adsorptionboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}/CONTCAR")

                            atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)

                            afile = open(f"{path}/POSCAR","r")
                            listlines = afile.readlines()
                            listlinespi = afile.readlines(7)
                            
                            listlines[7] = 'Selective Dynamics' + '\n' +'Cartesian' + '\n'
                            for line in range(8, 35):
                                listlines[line] = listlines[line].rstrip('\n') + " F F F"+'\n'

                            if adsorbates in ['C','O','H']:                       
                                for line in range(35, 37):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            
                            afile = open(f"{path}/POSCAR", "w")
                            afile.writelines(listlinespi)
                            afile.writelines(listlines)
                            afile.close()


                            create_potcar(path=f"{path}")
                            incar_tags = incar_tags_frequency.copy()
                            incar_tags['ISPIN'] = 1
                            create_incar(path=f"{path}", dict_tags=incar_tags)
                            create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                            create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
                        else:
                            print(f'Directory {path} already exists.')
        for adsorbate in [ 'H']:
            for adsorbates in ['H','CH','CH2','CH3','O','OH']:
                for site in sites:
                    for sit in sites2:
                        path = f"{main_path}/{mxene}/frequencyboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}"
                        # Create folder for each MAX phase
                        if not os.path.exists(path):
                            os.mkdir(path)
                            atoms_slab = read(f"{main_path}/{mxene}/adsorptionboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}/CONTCAR")
                            
                            atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)

                            afile = open(f"{path}/POSCAR","r")
                            listlines = afile.readlines()
                            listlinespi = afile.readlines(7)
                            
                            listlines[7] = 'Selective Dynamics' + '\n' +'Cartesian' + '\n'
                            for line in range(8, 35):
                                listlines[line] = listlines[line].rstrip('\n') + " F F F"+'\n'
                            if adsorbate in ['H']:                       
                                for line in range(35, 36):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['H', 'O']:                       
                                for line in range(36, 37):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['OH','CH']:                       
                                for line in range(36, 38):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['CH2']:                       
                                for line in range(36, 39):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['CH3']:                       
                                for line in range(36, 40):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            afile = open(f"{path}/POSCAR", "w")
                            afile.writelines(listlinespi)
                            afile.writelines(listlines)
                            afile.close()

                            create_potcar(path=f"{path}")
                            incar_tags = incar_tags_frequency.copy()
                            incar_tags['ISPIN'] = 1
                            create_incar(path=f"{path}", dict_tags=incar_tags)
                            create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                            create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
                        else:
                            print(f'Directory {path} already exists.')
        for adsorbate in [ 'CO']:
            for adsorbates in ['O']:
                for site in sites:
                    for sit in sites2:
                        path = f"{main_path}/{mxene}/frequencyboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}"
                        # Create folder for each MAX phase
                        if not os.path.exists(path):
                            os.mkdir(path)
                            atoms_slab = read(f"{main_path}/{mxene}/adsorptionboth/{mxene}#{adsorbate}-{adsorbates}-{site}-{sit}/CONTCAR")
                            
                            atoms_slab.write(filename=f"{path}/POSCAR", vasp5=True)

                            afile = open(f"{path}/POSCAR","r")
                            listlines = afile.readlines()
                            listlinespi = afile.readlines(7)
                            
                            listlines[7] = 'Selective Dynamics' + '\n' +'Cartesian' + '\n'
                            for line in range(8, 35):
                                listlines[line] = listlines[line].rstrip('\n') + " F F F"+'\n'
                            if adsorbate in ['CO']:                       
                                for line in range(35, 37):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['O']:                       
                                for line in range(37, 38):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['OH','CH']:                       
                                for line in range(36, 38):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['CH2']:                       
                                for line in range(36, 39):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            if adsorbates in ['CH3']:                       
                                for line in range(36, 40):
                                    listlines[line] = listlines[line].rstrip('\n') + " T T T" +'\n'
                            afile = open(f"{path}/POSCAR", "w")
                            afile.writelines(listlinespi)
                            afile.writelines(listlines)
                            afile.close()

                            create_potcar(path=f"{path}")
                            incar_tags = incar_tags_frequency.copy()
                            incar_tags['ISPIN'] = 1
                            create_incar(path=f"{path}", dict_tags=incar_tags)
                            create_kpoints(path=f"{path}", kpoints=(5, 5, 1))
                            create_submit_script(path=f"{path}", num_nodes=1, walltime_limit=5)
        quit()


stability = []
min_list = {}

if not os.path.exists(f"{main_path}/{mxene}/nebs"):
        print('Creating neb directories inputs')
        os.mkdir(f"{main_path}/{mxene}/nebs")

for mxene in ['W2N']:

    for adsorbate in [ 'CO','C', 'O','H2O','H2','H','OH','CH','CH2','CH3','CH4', 'CO2-X']:

            for site in sites:

                path = f"{main_path}/{mxene}/adsorption/{mxene}#{adsorbate}-{site}"

                with open(f"{main_path}/{mxene}/adsorption/{mxene}#{adsorbate}-{site}/OSZICAR", "r") as file:            
                    lines = file.readlines()
                    words_1 = lines[-1].replace("  "," ")
                    words = words_1.split(" ")
                    
                    energy = float(words[5])
                
                stability.append(energy)
            
            min_in = stability.index(min(stability))

            stability.clear()

            min_list[adsorbate] = list(sites.keys())[min_in]

    react = [[mxene+'#CO-O-'+min_list['CO']+'-'+min_list['O'], mxene+'#CO2-X-'+min_list['CO2-X'],1],
            [mxene+'#C-O-'+min_list['C']+'-'+min_list['O'], mxene+'#CO-'+min_list['CO'],1],
            [mxene+'#C-H-'+min_list['C']+'-'+min_list['H'], mxene+'#CH-'+min_list['CH'],0],
            [mxene+'#H-CH-'+min_list['H']+'-'+min_list['CH'], mxene+'#CH2-'+min_list['CH2'],0],
            [mxene+'#H-CH2-'+min_list['H']+'-'+min_list['CH2'], mxene+'#CH3-'+min_list['CH3'],0],
            [mxene+'#H-CH3-'+min_list['H']+'-'+min_list['CH3'], mxene+'#CH4-'+min_list['CH4'],0],
            [mxene+'#H-H-'+min_list['H']+'-'+min_list['H'], mxene+'#H2-'+min_list['H2'],1],
            [mxene+'#H-O-'+min_list['H']+'-'+min_list['O'], mxene+'#OH-'+min_list['OH'],0],
            [mxene+'#H-OH-'+min_list['H']+'-'+min_list['OH'], mxene+'#H2O-'+min_list['H2O'],0]
            ]
    for folder in react:
        if not os.path.exists(f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}"):
            os.mkdir(f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}")
        if not os.path.exists(f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb"):
            os.mkdir(f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb")

        path1 = f"{main_path}/{mxene}/adsorption/{folder[1]}/CONTCAR"
        path1_osz = f"{main_path}/{mxene}/adsorption/{folder[1]}/OSZICAR"
        path2 = f"{main_path}/{mxene}/adsorptionboth/{folder[0]}/CONTCAR"
        path2_osz = f"{main_path}/{mxene}/adsorptionboth/{folder[0]}/OSZICAR"
 
        if folder[2]==0: 
            path1_dest = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/initial.vasp"
            path1_dest_ml = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/initial.vasp"
            path1_dest_osz = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/OSZICAR-initial"
            path2_dest = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/final.vasp"
            path2_dest_ml = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/final.vasp"
            path2_dest_osz = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/OSZICAR-final"

        if folder[2]==1:
            path1_dest = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/final.vasp"
            path1_dest_ml = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/final.vasp"
            path1_dest_osz = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/OSZICAR-final"
            path2_dest = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/initial.vasp"
            path2_dest_ml = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/initial.vasp"
            path2_dest_osz = f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/OSZICAR-initial"


        shutil.copy(path1, path1_dest)
        shutil.copy(path1, path1_dest_ml)
        shutil.copy(path1_osz, path1_dest_osz)
        shutil.copy(path2, path2_dest)
        shutil.copy(path2, path2_dest_ml)
        shutil.copy(path2_osz, path2_dest_osz)

        atoms_slab = read(path1_dest_ml)
        atoms_slab = ase.build.sort(atoms_slab)
        if folder[2]==0: 
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/initial.traj")
        if folder[2]==1:    
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/final.traj")

        atoms_slab = read(path2_dest_ml)
        atoms_slab = ase.build.sort(atoms_slab)
        if folder[2]==0: 
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/final.traj")
        if folder[2]==1:    
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/mlneb/initial.traj")

        atoms_slab = read(path1_dest)
        atoms_slab = ase.build.sort(atoms_slab)
        if folder[2]==0: 
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/initial.vasp")
        if folder[2]==1:    
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/final.vasp")

        atoms_slab = read(path2_dest)
        atoms_slab = ase.build.sort(atoms_slab)
        if folder[2]==0: 
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/final.vasp")
        if folder[2]==1:    
            atoms_slab.write(filename=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}/initial.vasp")

        Nimages = 8  # elegir 8 imagenes
        os.chdir(f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}")

        initial = read('initial.vasp')
        final = read('final.vasp')
        constraints = initial.constraints  # copy constraints of initial

        write('POSCAR', initial)
        os.system('mkdir 00')
        os.system('mv POSCAR 00/')
        os.system('mv OSZICAR-initial 00/OSZICAR')

        images = [initial]
        for i in range(Nimages):
            image = initial.copy()
            image.set_constraint(constraints)
            images.append(image)
        images.append(final)

        neb = ase.neb.NEB(images, climb=False, k=0.5)
        neb.interpolate('idpp')

        for j in range(1, Nimages + 1):
            write('POSCAR', images[j])
            os.system('mkdir 0' + str(j))
            os.system('mv POSCAR 0' + str(j))


        write('POSCAR', final)
        os.system('mkdir 0' + str(Nimages + 1))
        os.system('cp POSCAR 0' + str(Nimages + 1))
        os.system('mv OSZICAR-final 0' + str(Nimages + 1)+'/OSZICAR')
        create_potcar(path=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}")
        incar_tags = incar_tags_nebs.copy()

        create_incar(path=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}", dict_tags=incar_tags)
        create_kpoints(path=f"{main_path}/{mxene}/nebs/{folder[0]}_{folder[1]}", kpoints=(5, 5, 1))

