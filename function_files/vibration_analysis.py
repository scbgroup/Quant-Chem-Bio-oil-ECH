# %%
from pylab import *
import pandas as pd
from ase import Atoms
from ase.io import espresso
from ase.io import write as wrt

# this function will create a QE scf input file given a file location and the atoms of interest
# the main purpose of this is to standardize input data and pseudopotentials

# create pseudopotential dictionary based on atoms in structure
# this will add ultrasoft pbe pseudos for any chemicals present 
def pseudo_dict(atoms):
    pp = {}
    chems = atoms.get_chemical_symbols()
    if 'Ru' in chems:
        pp['Ru'] = 'Ru.pbe-spn-rrkjus_psl.1.0.0.UPF'
    if 'Pt' in chems:
        pp['Pt'] = 'Pt.pbe-n-rrkjus_psl.1.0.0.UPF'
    if 'Cu' in chems:
        pp['Cu'] = 'Cu.pbe-dn-rrkjus_psl.1.0.0.UPF'
    if 'C' in chems:
        pp['C'] = 'C.pbe-n-rrkjus_psl.1.0.0.UPF'
    if 'O' in chems:
        pp['O'] = 'O.pbe-n-rrkjus_psl.1.0.0.UPF'
    if 'H' in chems:
        pp['H'] = 'H.pbe-rrkjus_psl.1.0.0.UPF' 
    return pp

# write scf input file for 
def scf_input(in_file, atoms, pre, stype):
    if stype == 'gas':
        input_data = {'control' : {'calculation':'scf',
                                   'prefix': pre,
                                   'pseudo_dir':'/mnt/research/scbgroup/christina/bio-oil-ECH/pseudo',
                                  'tstress': True,
                                  'tprnfor': True}, 
                  'system': {'ecutwfc': 35,
                             'ecutrho': 279.3},
                    'electrons': {'conv_thr': 1e-12,
                                  'mixing_beta': 0.5},
                     'kpoints':(1,1,1)}
        pseudopot = pseudo_dict(atoms)
        espresso.write_espresso_in(in_file, atoms, input_data, pseudopot, kpts = (1,1,1))
    if stype == 'ads':
        input_data = {'control' : {'calculation':'scf',
                               'prefix': pre,
                               'pseudo_dir':'/mnt/research/scbgroup/christina/bio-oil-ECH/pseudo',
                              'tstress': True,
                              'tprnfor': True}, 
              'system': {'ecutwfc': 35,
                         'ecutrho': 279.3,
                         'occupations': 'smearing',
                         'degauss': 0.1,
                         'smearing':'m-p'},
                'electrons': {'conv_thr': 1e-12,
                              'mixing_beta': 0.5},
                 'kpoints':(1,1,1)}
        pseudopot = pseudo_dict(atoms)
        espresso.write_espresso_in(in_file, atoms, input_data, pseudopot, kpts = (1,1,1))

# write ph input file for gas phase
def ph_input(in_file, pre, nat_to_do, surf):
    if surf == None:
        input_data = {'inputph': {'tr2_ph': 1.0e-14,
                      'prefix': pre,
                      'ldisp': True,
                      'nq1': 1, 'nq2': 1, 'nq3': 1}}
        #qpts = [(0, 0, 0)]
        with open(in_file, 'w') as file:
            espresso.write_espresso_ph(file, input_data, qpts = (0.0, 0.0, 0.0))
    else:
        if surf == 'Pt':
            tr = 1.0e-7
        if surf == 'Ru':
            tr = 5.0e-7
        input_data = {'inputph': {'tr2_ph': tr,
                      'prefix': pre,
                      'ldisp': True,
                      'nq1': 1, 'nq2': 1, 'nq3': 1,
                      'nat_todo': len(nat_to_do)}}
        #qpts = [(0, 0, 0)]
        with open(in_file, 'w') as file:
            espresso.write_espresso_ph(file, input_data, qpts = (0.0, 0.0, 0.0),
                                       nat_todo_indices=nat_to_do)

# the input to the following functions is generally the output file from Quantum Espresso calculations

# scf_out is a function to pull the final energy and information about the geometry for scf calculation

def scf_out(output_file):
    # Open the output file and read configurations
    with open(output_file, 'r') as file:
        for atoms in espresso.read_espresso_out(file, results_required=True):
            # collect the final energy
            # The package converts the unit from Ry to eV
            energy = atoms.get_total_energy()
            #energy_list += [energy]

            # collect the trajectory through the optimization
            #traj += [atoms]
            
    return energy

# ph_out is a function to pull the vibrational freuqencies in THz
def ph_out(output_file):
    freqs = []

    # open output file and read data
    with open(output_file, 'r') as file:
        res = espresso.read_espresso_ph(file)
        
    # pull freuqnecy data for each q-point and append to one array
    for i,ii in enumerate(res):
        freqs.append(res[ii]['freqs'])
        
    return freqs[0]




