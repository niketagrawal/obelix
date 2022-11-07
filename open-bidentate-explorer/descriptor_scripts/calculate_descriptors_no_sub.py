from operator import concat
import numpy as np
import os
import morfeus as mf
from morfeus import BiteAngle, ConeAngle, BuriedVolume, Dispersion, SASA, read_xyz
from xtb.utils import get_method
# import shutil
import matplotlib.pyplot as plt
import pandas as pd


def boltzmann_average(values, weights):
    return np.average(values, weights=weights, axis=0)
    

def boltzmann_factor(E,T=298.15,kb=3.166811563*10**(-6)): #k_b is in Hartree/K
    return np.exp(-E/(kb*T))


def boltzmann_std(values, weights):
    avg = boltzmann_average(values,weights)
    var = np.average((values-avg)**2,weights=weights,axis=0)
    return np.sqrt(var)


def find_bidentate(xyz):
    
    with open(xyz) as file:
        listfile = []
    
        for line in file:
            listfile.append(line.strip())
    indices = [0]*3
    atoms = []
    for index, atom in enumerate(listfile):
        if atom[0] == 'P':
            if not indices[0]:
                indices[0] = index - 1
            else:
                indices[2] = index - 1
            atoms.append(atom[0])
        else:
            if atom[:2] == 'Ir':
                indices[1] = index - 1
                atoms.append(atom[:2])
    return indices


def dataframe_from_dictionary(dictionary):
    Dataframe = pd.DataFrame.from_dict({i: dictionary[i]
                           for i in dictionary.keys()}, orient='index')
    return Dataframe


def calculate_descriptors_fa_conformers(crest_conformers, charge):
    with open(crest_conformers) as file:
        nr_of_atoms_per_molecule = int(file.readline()) + 2
        for nr_of_conf , _ in enumerate(file):
            pass
        indices = np.arange(0, nr_of_conf + nr_of_atoms_per_molecule, nr_of_atoms_per_molecule)
        
        if not os.path.isdir('conformers'):    
            os.makedirs('conformers')
        os.chdir('conformers')
        with open(crest_conformers) as file:
            lines = file.readlines()

            for i in range(len(indices) - 1):
        
                with open("conformer_{}.xyz".format(i), "w") as f:
                    f.writelines(lines[indices[i]:indices[i+1]])
    conformer_file_names = os.listdir()
    
    #### change for readability
    bite_angle_per_conformer = np.zeros(len(indices) - 1)
    cone_angle_per_conformer = np.zeros(len(indices) - 1)
    burried_volume_per_conformer = np.zeros(len(indices) - 1)
    dipole_per_conformer = np.zeros(len(indices) - 1)
    ea_per_conformer = np.zeros(len(indices) - 1)
    electrophilicity_per_conformer  = np.zeros(len(indices) - 1)
    nucleophilicity_per_conformer = np.zeros(len(indices) - 1)
    electrofugality_per_conformer = np.zeros(len(indices) - 1) 
    nucleofugality_per_conformer = np.zeros(len(indices) - 1)
    homo_per_conformer = np.zeros(len(indices) - 1)
    lumo_per_conformer = np.zeros(len(indices) - 1)
    ip_per_conformer = np.zeros(len(indices) - 1)
    dispersion_per_conformer = np.zeros(len(indices) - 1)
    sasa_per_conformer = np.zeros(len(indices) - 1)
    
    for conf_index, conformer in enumerate(conformer_file_names):
        elements, coordinates = read_xyz(conformer)
        bidentate = find_bidentate(conformer)
        bite_angle_per_conformer[conf_index] = BiteAngle(coordinates, bidentate[1], bidentate[0], bidentate[2]).angle
        cone_angle_per_conformer[conf_index] = ConeAngle(elements, coordinates, bidentate[1]).cone_angle
        burried_volume_per_conformer[conf_index] = BuriedVolume(elements, coordinates, bidentate[1]).fraction_buried_volume
        dispersion_per_conformer[conf_index] = Dispersion(elements, coordinates).atom_p_int[1]
        sasa_per_conformer[conf_index] = SASA(elements, coordinates).area
        
        # Electronic descriptors
        
        instance_electronic = mf.XTB(elements, coordinates, solvent='ch2cl2')
        
        ip_per_conformer[conf_index] = instance_electronic.get_ip()
        dipole_per_conformer[conf_index] = instance_electronic.get_dipole().dot(instance_electronic.get_dipole())
        ea_per_conformer[conf_index] = instance_electronic.get_ea()
        electrofugality_per_conformer[conf_index] = instance_electronic.get_global_descriptor(variety = 'electrofugality')
        nucleofugality_per_conformer[conf_index] = instance_electronic.get_global_descriptor(variety = 'nucleofugality')
        nucleophilicity_per_conformer[conf_index] = instance_electronic.get_global_descriptor(variety = 'nucleophilicity')
        electrophilicity_per_conformer[conf_index] = instance_electronic.get_global_descriptor(variety = 'electrophilicity')

        # HOMO-LUMO
        
        homo_per_conformer[conf_index] = instance_electronic.get_homo()
        lumo_per_conformer[conf_index] = instance_electronic.get_lumo()

    os.chdir("..")
    
    if len(indices) == 2:
        energies = np.genfromtxt(os.getcwd() + '/crest.energies')[1] + float(lines[1].strip())
        bite_angle_avg = bite_angle_per_conformer[0]
        cone_angle_avg = cone_angle_per_conformer[0]
        dipole_avg = dipole_per_conformer[0]
        ea_avg = ea_per_conformer[0]
        electrofugality_avg = electrofugality_per_conformer[0]
        electrophilicity_avg = electrophilicity_per_conformer[0]     
        nucleofugality_avg = nucleofugality_per_conformer[0]
        sasa_avg = sasa_per_conformer[0]
        homo_lumo_gap_avg = homo_per_conformer[0] - lumo_per_conformer[0]
        dispersion_avg = dispersion_per_conformer[0]
        ip_avg = ip_per_conformer[0]
        burried_volume_avg = burried_volume_per_conformer[0]
        nucleophilicity_avg = nucleofugality_per_conformer[0]
    else:
        denergies = np.genfromtxt(os.getcwd() + '/crest.energies')[:, 1]
        energies = denergies + float(lines[1].strip())
        bite_angle_avg = boltzmann_average(bite_angle_per_conformer, energies)
        cone_angle_avg = boltzmann_average(cone_angle_per_conformer, energies)
        dipole_avg = boltzmann_average(dipole_per_conformer, energies)
        ea_avg = boltzmann_average(ea_per_conformer, energies)
        electrofugality_avg = boltzmann_average(electrofugality_per_conformer, energies)
        electrophilicity_avg = boltzmann_average(electrophilicity_per_conformer, energies)     
        nucleofugality_avg = boltzmann_average(nucleofugality_per_conformer, energies)
        sasa_avg = boltzmann_average(sasa_per_conformer, energies)
        homo_lumo_gap_avg = boltzmann_average(homo_per_conformer - lumo_per_conformer, energies)
        dispersion_avg = boltzmann_average(dispersion_per_conformer, energies)
        ip_avg = boltzmann_average(ip_per_conformer, energies)
        burried_volume_avg = boltzmann_average(burried_volume_per_conformer, energies)
        nucleophilicity_avg = boltzmann_average(nucleophilicity_per_conformer, energies)
        
    descriptors[crest_conformers.split('/')[-2]] = {'Bite angle'       : bite_angle_avg,
                                                    'Cone angle'       : cone_angle_avg,
                                                    'Energy'           : np.mean(energies),
                                                    'Dipole'           : dipole_avg,
                                                    'EA'               : ea_avg,
                                                    'Electrofugality'  : electrofugality_avg,
                                                    'Electrophilicity' : electrophilicity_avg,
                                                    'Nucleofugality'   : nucleofugality_avg, 
                                                    'SASA'             : sasa_avg, 
                                                    'HOMO-LUMO gap'    : homo_lumo_gap_avg, 
                                                    'Dispersion'       : dispersion_avg,
                                                    'IP'               : ip_avg, 
                                                    'Burried Volume'   : burried_volume_avg,
                                                    'Nucleophilicty'   : nucleophilicity_avg}
                                                    
    
def loop_for_descriptors(directory):
    global descriptors
    descriptors = {}

    os.chdir(directory)
    folder_names = os.listdir()

    mace_averaging_dictionary = {}
    
    for folder in folder_names:
        os.chdir(directory + '/' + folder)
        if 'SP' in folder:
            charge_complex = 3
        else:
            charge_complex = 0
        calculate_descriptors_fa_conformers(os.getcwd() + '/crest_conformers.xyz', charge_complex)

        # average over MACE    

        if folder[:-1] not in mace_averaging_dictionary.keys():
            for key in descriptors.keys():
                mace_averaging_dictionary[folder[:(len(folder) - 3)]] = descriptors[key]
                print(folder[:-1])
                print(descriptors[key])
            
            mace_averaging_dictionary[folder[:(len(folder) - 3)]]['MACE Conformers'] = 1
                    
        else:
            for column, _ in mace_averaging_dictionary[folder[:(len(folder) - 3)]].items():
                print(folder[:-1])
                if column != 'MACE Conformers':
                    mace_averaging_dictionary[folder[:(len(folder)-3)]][column] += descriptors[folder][column]
                else:
                    mace_averaging_dictionary[folder[:(len(folder)-3)]]['MACE Conformers'] += 1
        
        
    df = dataframe_from_dictionary(mace_averaging_dictionary)
    df.iloc[:, :] = df.iloc[:, :].div(df['MACE Conformers'], axis=0)
    df = df.reindex(sorted(df.columns), axis=1).drop(columns='MACE Conformers')

    # print(df)
    return df    

descriptors_dataframe_OH_no_sub = loop_for_descriptors(os.getcwd() + '/CompletedCrestCalc/OH_no_sub')
descriptors_dataframe_OH_no_sub.columns = [str(col).replace(' ', '_') + '_OH' for col in descriptors_dataframe_OH_no_sub.columns]

# return to the parent dir
os.chdir('..')
os.chdir('..')
os.chdir('..')

# descriptors_dataframe_SP = loop_for_descriptors(os.getcwd() + '/CompletedCrestCalc/SP')
# descriptors_dataframe_SP.columns = [str(col).replace(' ', '_') + '_SP' for col in descriptors_dataframe_SP.columns]

# return to the parent dir

# os.chdir('..')
# os.chdir('..')
# os.chdir('..')


# concatenated = descriptors_dataframe_SP.join(descriptors_dataframe_OH)
descriptors_dataframe_OH_no_sub.to_excel('figures/new_method/Descriptors_no_sub.xlsx')
# correlation_matrix = concatenated.corr()**2
# correlation_matrix.to_excel('figures/new_method/corr2.xlsx')

# X = np.vstack([SP_descriptor, np.ones(len(SP_descriptor))]).T
# x_min = np.min(SP_descriptor)
# x_max = np.max(SP_descriptor)
# model, resid = np.linalg.lstsq(X, OH_descriptor, rcond=None)[:2]

# R2 = 1 - resid/OH_descriptor.size/OH_descriptor.var()
# # mx + c 
# m = model[0]
# c = model[1]

# plt.figure()
# plt.plot([x_min, x_max], [model[0]*x_min + model[1], model[0]*x_max + model[1]], color = 'red')
# plt.scatter(SP_descriptor, OH_descriptor, color =  'black')

# plt.legend(['{}x + {}'.format(np.round(m, 2), np.round(c)), r'$R^2 = {}$'.format(R2)]) 

# plt.xlabel('Bite angle SP')
# plt.ylabel('Bite angle OH')
# plt.savefig('figures/SP_ba_vs_OH_ba.png')
# print(dataframe_to_dictionary(d))