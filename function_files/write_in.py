from pylab import *
import pandas as pd
from ase.io import espresso
from ase.constraints import FixAtoms, constrained_indices


# this function will create a QE input file given a file location and the atoms of interest
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

# write input file
def QE_input(in_file, atoms):
    input_data = {'control' : {'calculation':'relax',
                           'pseudo_dir':'/mnt/research/scbgroup/christina/bio-oil-ECH/pseudo',
                           'forc_conv_thr': 0.001,
                           'nstep': 500}, 
              'system': {'ecutwfc': 35,
                         'ecutrho': 279.3,
                         'vdw_corr': 'dft-d3',
                         'occupations': 'smearing',
                         'smearing': 'm-p',
                         'degauss':  0.1},
                'electrons': {'conv_thr': 1e-4,
                              'mixing_beta': 0.5}}
    pseudopot = pseudo_dict(atoms)
    with open(in_file, 'w') as file:
        espresso.write_espresso_in(file, atoms, input_data, pseudopot)

# write a list of pseudo potentials
def pseudo_list(atoms):
    pp = ['ATOMIC_SPECIES\n']
    chems = atoms.get_chemical_symbols()
    if 'Ru' in chems:
        pp+=['Ru 101.07 Ru.pbe-spn-rrkjus_psl.1.0.0.UPF\n']
    if 'Pt' in chems:
        pp+=['Pt 195.084 Pt.pbe-n-rrkjus_psl.1.0.0.UPF\n']
    if 'Cu' in chems:
        pp+=['Cu 63.546 Cu.pbe-dn-rrkjus_psl.1.0.0.UPF\n']
    if 'C' in chems:
        pp+=['C 12.011 C.pbe-n-rrkjus_psl.1.0.0.UPF\n']
    if 'O' in chems:
        pp+=['O 15.999 O.pbe-n-rrkjus_psl.1.0.0.UPF\n']
    if 'H' in chems:
        pp+=['H 1.008 H.pbe-rrkjus_psl.1.0.0.UPF\n'] 
    return pp
    
# write NEB input file
def NEB_input(first_image, last_image, filename='neb.in'):
    # define psuedopotentials list from first image
    pp = pseudo_list(first_image)
    
    """
    Writes a Quantum ESPRESSO NEB input file from the first and last images.
    
    Parameters:
        first_image (ase.Atoms): The initial atomic configuration.
        last_image (ase.Atoms): The final atomic configuration.
        filename (str): The name of the output NEB input file.
    """
    
    with open(filename, 'w') as f:
        # write NEB controls in '&PATH' section
        neb_keys = ['BEGIN\n','BEGIN_PATH_INPUT\n','&PATH\n',
                    ' string_method = "neb"\n',
                    ' restart_mode = "restart"\n',
                    ' nstep_path = 500\n',
                    ' num_of_images = 5\n',
                    ' opt_scheme = "broyden"\n',
                    ' path_thr = 1.0d-4\n',
                    ' ds = 0.2\n','/\n','END_PATH_INPUT\n']
        f.writelines(neb_keys)
        # write SCF controls
        scf_keys = ['BEGIN_ENGINE_INPUT\n','&CONTROL\n',
                   ' calculation = "relax"\n',
                   ' nstep = 500\n',
                   ' forc_conv_thr  = 0.001\n',
                   ' pseudo_dir  = "/mnt/research/scbgroup/christina/bio-oil-ECH/pseudo"\n','/\n',
                   '&SYSTEM\n',
                   ' ecutwfc = 35\n', 
                   ' ecutrho = 279.3\n',
                   ' vdw_corr = "dft-d3"\n',
                   ' occupations = "smearing"\n',
                   ' smearing = "m-p"\n',
                   ' degauss = 0.1\n',
                   f' ntyp = {len(pp)-1}\n',
                   f' nat = {len(first_image)}\n',
                   ' ibrav = 0\n','/\n',
                   '&ELECTRONS\n',
                   ' conv_thr = 0.0001\n',
                   ' mixing_beta = 0.5\n','/\n',
                   '&IONS\n','/\n',
                   '&CELL\n','/\n',
                   '&FCP\n','/\n',
                   '&RISM\n','/\n']
        for i,p in enumerate(pp):
            scf_keys += [p]
        f.writelines(scf_keys)

        # write atomic posistions for first image
        f.writelines(['\n','BEGIN_POSITIONS\n', 'FIRST_IMAGE\n','ATOMIC_POSITIONS (angstrom)\n'])
        # get constraint indices
        consts = constrained_indices(first_image)
        for i,atom in enumerate(first_image):
            # fix atoms with constraints
            if i in consts:
                f.write(f'{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f} 0 0 0\n')
            else:
                f.write(f'{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n')
        
        # write atomic posistions for last image
        f.writelines(['LAST_IMAGE\n','ATOMIC_POSITIONS (angstrom)\n'])
        # get constraint indices
        consts = constrained_indices(last_image)
        for i,atom in enumerate(last_image):
            # fix atoms with constraints
            if i in consts:
                f.write(f'{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f} 0 0 0\n')
            else:
                f.write(f'{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n')
        f.writelines(['\n','END_POSITIONS\n'])
                
        # write kpt and cell params
        params = ['K_POINTS { gamma }\n','CELL_PARAMETERS angstrom\n']
        for i,c in enumerate(first_image.get_cell()[:]):
            para = ''
            for j,p in enumerate(c):
                para += f'{p} '
            params += [f'{para}\n']   
        f.writelines(params)

        # write end
        f.writelines(['END_ENGINE_INPUT\n','END'])