from pylab import *
import pandas as pd
from ase import Atoms
from ase.io import espresso
from ase.geometry.analysis import Analysis as geoa
from ase.io import write as wrt

# this function saves several image and coordinate data files for the start, end, and optimization process

def save_images(path, pre, traj):
    '''
    path : where the files will be saved
    pre : prefix for files
    traj: atoms trajectory
    '''
    # *.xyz file with atomic info for each step of the optimization
    wrt(path + pre + '.xyz', traj)
    
    # Starting image
    ims = traj[0]

    # top view
    wrt(path + pre + '_start_top.png', ims,
                show_unit_cell = 0, rotation = "0x,0y,0z")

    # side view
    wrt(path + pre + '_start_side.png', ims,
                show_unit_cell = 0, rotation="90x,0y,0z")
    
    # Optimized image
    imf = traj[-1]

    # top view
    wrt(path + pre + '_opt_top.png', imf,
                show_unit_cell = 0, rotation = "0x,0y,0z")

    # side view
    wrt(path + pre + '_opt_side.png', imf,
                show_unit_cell = 0, rotation="90x,0y,0z")

    return print("files saved")

# this function calculates the bond length(s) given an atoms input image and the atoms of interest

def bond_length(im,a,b):
    '''
    a : atom A (i.e.'C')
    b : atom B (i.e. 'Pt')
    im : image to be analyzed
    '''
    bonds = geoa(im).get_bonds(a,b)
    
    if bonds == [[]]:
        return None
    else:
        length = geoa(im).get_values(bonds)
        return length[0]

# this function calculates the average and standard deviation of the distance between the ring and the surface

def ring_height(im, O, C, H, surf, Nsurf):
    '''
    returns the average ring height and standard deviation 
        from the top surface layer in angstroms
        
    im : the atoms of the image for analysis
    O: number of oxygen in the ring
    C: number of carbon in the ring
    H: number of H in the ring
    surf: surface atom type (i.e. 'Pt')
    Nsurf: number of top layer surface atoms
    '''
    ringH = [atom.index for atom in im if atom.symbol == 'H'][0:H]
    ringO = [atom.index for atom in im if atom.symbol == 'O'][0:O]
    ringC = [atom.index for atom in im if atom.symbol == 'C'][0:C]
    surfAt = [atom.index for atom in im if atom.symbol == surf][-Nsurf:]

    ring = ringO + ringC

    ring_z = []
    for i in im[ring].get_positions():
        ring_z += [i[-1]]

    avg_rh = average(ring_z)

    surf_z = []

    for i in im[surfAt].get_positions():
        surf_z += [i[-1]]

    avg_sh = average(surf_z)

    ring_height = average(ring_z) - average(surf_z)
    rh_std = std(ring_z) + std(surf_z)
    return ring_height, rh_std  

# the input to the following functions is generally the output file from Quantum Espresso calculations

# ase_out is a function to pull the final energy and information about the geometry for each optimization step

def ase_out(output_file):
    energy_list = []
    traj = []
    # cell_list = []
    # coords_list = []
    # sym_list = []
    
    # Open the output file and read configurations
    with open(output_file, 'r') as file:
        for atoms in espresso.read_espresso_out(file, results_required=True):
            # collect the total energy for each opt step
            # The package converts the unit from Ry to eV
            energy = atoms.get_total_energy()
            energy_list += [energy]

            # collect the trajectory through the optimization
            traj += [atoms]
            
    return array(energy_list), traj

# function to compile the energy calculated in each scf iteration and each geometry optimization step

def energy(outpath):
    factor = 27.211396/2 # eV/2 Ry to convert energies
    itern = [] # iteration number
    scf_cycle = [] # scf cycle number
    it_energy = [] # iteration number
    conv_energy = [] # converged energy
    
    with open(outpath) as file:
        # set the scf cycle starting with cycle 1
        scf = 1 
        
        for line in file: 
            # track the scf iteration and cycle
            if 'iteration #' in line: 
                it_st = line
                itern += [int(it_st[17:20])]
                scf_cycle += [scf]
            # track the energy of the iteration    
            if 'total energy              =' in line: 
                string = line
                it_energy += [float(string[32:49])*factor]
            # track the converged energy for the scf cycle    
                if '!    total energy' in string: 
                    conv_energy += [float(string[32:49])*factor]
                    # increase the scf cycle
                    scf += 1
    # collect the scf cycle number, iteration number, and iteration energy in a data frame                
    df = pd.DataFrame({'scf cycle': array(scf_cycle),
                        'iteration': array(itern), 
                        'iteration energy': array(it_energy)})      
    
    # determine max number of iterations in each geo optimization
    maxit = max(array(itern))
    
    # determine max number of optimization steps
    maxscf = max(array(scf_cycle))

    # rearrange data so that each column represents one optimization step with the values being the energy from each iteration within that step
    dict_ = {'iteration':arange(1,maxit+1)}
    keys = [f'scf cycle {i}' for i in arange(1,maxscf+1)]

    for i in arange(0,maxscf):
        iter_energy = array(df[df["scf cycle"] == i+1 ]['iteration energy'])
        if len(iter_energy < maxit):
            iter_energy = array(append(iter_energy, [None] * (maxit-len(iter_energy))))
        key = keys[i]
        dict_[key] = iter_energy

    # return the converged energy of each optimization step, data frame with the iteration and step details
    return array(conv_energy), pd.DataFrame(dict_)

# function to compile the force outputs for each optimization step
def force(outpath):
    force = [] # total force after geometry optimization step
    
    with open(outpath) as file:
        # set the bfgs cycle starting with cycle 0
        bfgs = 0 
        
        for line in file:
            # pull the force convergence threshold
            if 'force convergence thresh. =' in line:
                f_thresh = float(line[33:])
            # pull the total force after each optimization cycle    
            if 'Total force =' in line:
                force += [float(line[23:31])]
                bfgs += 1
    bfgs_cycle = arange(0,bfgs)

    # output of this function is the optimization step and corresponding force
    return f_thresh, pd.DataFrame({'Optimization Step':bfgs_cycle, 
                         'Total Force': force})
    
# function to extract computational time
# outputs total time, average iteration time, and dataframe with SCF cycle columns and iteration time in rows
def timing(outpath):
    iter = ['init'] # iteration number, starting with the initialization info 
    scf_cycle = [0] # scf cycle number 
    it_time = [] # time up to end of iteration
    it_len = [] # actual iteration time
    
    with open(outpath) as file:
        # set the scf cycle starting with cycle 1
        scf = 0 
        
        for line in file: 
            # track the scf iteration and cycle
            if 'iteration #' in line: 
                iter += [int(line[17:20])]
                scf_cycle += [scf]
            # track the time up to the end of the scf iteration
            if 'total cpu time spent up to now is' in line:
                it_time += [float(line[40:50])]
            # check if start of new optimization step
            if 'Self-consistent Calculation' in line:
                scf += 1
            # check if end of SCF calculation to update cycle number
            # add iteration 0 to collect time for force calculation except for last SCF cycle
            if 'End of self-consistent calculation' in line:
                scf_cycle += [scf]
                iter += [0]
            
        
        # initialize actual length of time per iteration list
        init_t = it_time[0] 
        it_len = []
    
        # calculate length of iteration
        for i,t in enumerate(it_time[1:]):
            it_len += [t - it_time[i]]

        # collect timing data in a dataframe
        df = pd.DataFrame({'SCF Cycle': scf_cycle[1:-1], 
                           'Iteration': iter[1:-1], 
                           'Time spent': it_time[1:], 
                           'Iteration Time': it_len})
        
        # arrange the data by length of time in each iteration and collect them on a per SCF cycle basis
        # determine max number of iterations per cycle and max SCF cycles
        maxit = max(iter[1:])
        maxscf = max(scf_cycle)
        
        # initialize a dictionary with columns representing each SCF cycle
        dict = {'iteration':arange(0,maxit+1)}
        # create keys for each SCF cycle
        keys = [f'scf cycle {i}' for i in arange(1,maxscf+1)]
        
        for i in arange(0,maxscf):
            # collect iteration times for SCF cycle i >= 1
            iter_time = array(df[df["SCF Cycle"] == i+1 ]['Iteration Time'])
            # add 'None' to row after iterations if total iterations in cycle are less than max
            if len(iter_time) < maxit:
                iter_time = array(append(iter_time, [None] * (maxit-len(iter_time)+1)))
            # set the key for the cycle
            key = keys[i]
            # add SCF cycle i to the dictionary
            dict[key] = iter_time
          
        # determine total calculation time
        total_t = sum(it_len) + init_t
    
        # average iteration time
        avg_t = average(it_len)

        return total_t, avg_t, dict
