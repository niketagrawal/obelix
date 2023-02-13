import re, glob
import pandas as pd
import numpy as np
import os, ast
import morfeus as mf


class DFT_descriptors:
    def __init__(self, log_file, nr_of_atoms, min_donor, max_donor, metal):
        
        self.log_file = log_file # name of log file
        self.nr_of_atoms = nr_of_atoms # nr of atoms in molecule
        self.metal = metal  # metal index
        self.min_donor = min_donor # min donor index
        self.max_donor = max_donor # max donor index
        
        print(vars(self))
        
    def donor_lone_pair_occupancy(self, atom_type_min, atom_type_max):
        print("Extracting donor atoms lone pair orbital occupancy")
        # Count nr of times " NATURAL POPULATIONS:  Natural atomic orbital occupancies" appears and take the last one
        count = 0
        count2 = 0
        with open(self.log_file) as file: # Use file to refer to the file object
            data = file.readlines()
            for line_index, line in enumerate(data):
                if "NATURAL POPULATIONS:  Natural atomic orbital occupancies" in line:    
                    count += 1
                    if count == 3:
                        counting_line = line_index + 4
                if "electrons found in the effective core potential" in line:    
                    count2 += 1
                    if count2 == 3:
                        counting_line_ = line_index - 1
                            
                # if 
                # final_line = 
        file.close()
        occupancy_raw_data = data[counting_line: counting_line_]
        
        occupancy_final_data = []
        orbital_dictionary = {}

        for occupancy in occupancy_raw_data:
            occupancy = occupancy.strip("\n")
            occupancy = occupancy.strip("\r")
            occupancy = occupancy.split()
            if occupancy != []:
                occupancy_final_data.append((occupancy))
        
        for orbital in occupancy_final_data:
            if int(orbital[2]) not in orbital_dictionary.keys():
                orbital_dictionary[int(orbital[2])] = {orbital[3] + '_' + orbital[4] + orbital[5]: float(orbital[6])}
            else:
                orbital_dictionary[int(orbital[2])][orbital[3] + '_' + orbital[4] + orbital[5]] = float(orbital[6])
        # divided by 2 since 2 electrons form the lone pair and the 3s orbital
        
        returns = []
        if atom_type_min == 'N':
            returns.append(orbital_dictionary[self.min_donor + 1]['S_Val(2S)'])
            
        else:
            returns.append(orbital_dictionary[self.min_donor + 1]['S_Val(3S)'])
            
        if atom_type_max == 'N':
            returns.append(orbital_dictionary[self.max_donor + 1]['S_Val(2S)'])
            
        else:
            returns.append(orbital_dictionary[self.max_donor + 1]['S_Val(3S)'])
            
        return returns[0], returns[1]
 
    def NBO_charge(self):
        print("Extracting donor atoms and metal natural charge")
        # Count nr of times "     Atom  No    Charge         Core      Valence    Rydberg      Total" appears and take the last one
        count = 0

        with open(self.log_file) as file: # Use file to refer to the file object
            data = file.readlines()
            for line_index, line in enumerate(data):
                if "    Atom  No    Charge         Core      Valence    Rydberg      Total" in line:    
                    count += 1
                    if count == 3:
                        counting_line = line_index + 2
                # if 
                # final_line = 
        file.close()
        charge_raw_data = data[counting_line: counting_line + self.nr_of_atoms]
        charge_final_data = []
        charge_dictionary = {}

        for charge in charge_raw_data:
            charge = charge.strip("\n")
            charge = charge.strip("\r")
            charge = charge.split()
            if charge != []:
                charge_final_data.append(charge)

        for charge in charge_final_data:
            charge_dictionary[int(charge[1])] = {float(charge[2])}
  
        # divided by 2 since 2 electrons form the lone pair and the 3s orbital
        
        return list(charge_dictionary[self.metal + 1]), list(charge_dictionary[self.min_donor + 1]), list(charge_dictionary[self.max_donor + 1])


    def mulliken_charge(self):
        print("Extracting Mulliken charge")
        with open(self.log_file) as file: # Use file to refer to the file object
            data = file.readlines()
            for line_index, line in enumerate(data):
                if "Mulliken charges" in line:        
                    Mulliken_raw_data = data[line_index + 2:line_index + self.nr_of_atoms + 2]
                    break
        file.close()
        
        mulliken_final_data = []
        
        for raw_charge in Mulliken_raw_data:
            raw_charge = raw_charge.strip("\n")
            raw_charge = raw_charge.strip("\r")
            raw_charge = raw_charge.split()
            mulliken_final_data.append(float(raw_charge[-1]))
        # print(mulliken_final_data)
        return mulliken_final_data[self.metal], mulliken_final_data[self.min_donor], mulliken_final_data[self.max_donor]


    def P_orbital_occupation(self):
        with open(self.log_file) as file: # Use file to refer to the file object
            data = file.readlines()
            count = 0
            for line_index, line in enumerate(data):
                if "NBO                        Occupancy    Energy   (geminal,vicinal,remote)" in line:        
                    count += 1
                    counting_line = line_index + 1
            
            
        raw_data_metal_min = []
        raw_data_metal_max = []
        final_data_min = []
        final_data_max = []
        
        for line_ in data[counting_line:len(data)]: 
            
            if ' (Enter /software/sse/easybuild/prefix/software/Gaussian/16.C.01-avx2-nsc1/g16/l701.exe)' in line_:
                print("Brok")
                break
            
            if "BD" in line_:      
                # remove all non numbers and dot from the lines if line contains BD or BD* 
                # (generalized as BD since BD* contains BD)
                line = re.sub('[^\d\.]', ' ', line_)          
                  
                line = line.split()
                print(line)    
                if str(self.metal + 1) in line:
                    if str(self.min_donor + 1) in line:
                        raw_data_metal_min.append(line) 
                        if len(raw_data_metal_min) == 2:
                            Rh_P_min_antibonding = float(line[4])
                        else:
                            Rh_P_min_bonding = float(line[4])
                    else:
                        Rh_P_min_antibonding = None
                        Rh_P_min_bonding = None
                    if str(self.max_donor + 1) in line:
                        raw_data_metal_max.append(line)
                        if len(raw_data_metal_max) == 2:                    
                            Rh_P_max_antibonding = float(line[4])
                        else:
                            Rh_P_max_bonding = float(line[4])
                    else:
                        Rh_P_max_antibonding = None
                        Rh_P_max_bonding = None                        
                else:
                    if str(self.min_donor + 1) in line:
                        final_data_min.append(float(line[4])) 

                    if str(self.max_donor + 1) in line:
                        final_data_max.append(float(line[4]))
        ### the return contains 
        ### final_data_min = [bonding orbital 1, bonding orbital 2, bonding orbital 3, antibonding orbital 1, antibonding orbital 2, antibonding orbital 3]
        ### final_data_max = [bonding orbital 1, bonding orbital 2, bonding orbital 3, antibonding orbital 1, antibonding orbital 2, antibonding orbital 3]
        ### metal_min_bv_bonding
        ### metal_min_bv_antibonding
        ### metal_max_bv_bonding
        ### metal_max_bv_antibonding
        # print(final_data_min, final_data_max)
        
        return final_data_min, final_data_max, Rh_P_min_bonding, Rh_P_min_antibonding, Rh_P_max_bonding, Rh_P_max_antibonding

        
log_files = glob.glob(os.getcwd() + '/final_log_files/*.log')

data = pd.read_excel("mf_BD_Rh_test4.xlsx")

dataframe = pd.DataFrame(pd.DataFrame(columns=['filename','LP_occupancy_min','LP_occupancy_max','NBO_charge_Rh', 'NBO_charge_min', 
                                               'NBO_charge_max','Rh_mulliken_charge', 'min_mulliken_charge', 'max_mulliken_charge',
                                               'antibond_min_donor_1', 'antibond_min_donor_2', 'antibond_min_donor_3', 'bond_min_donor_1', 'bond_min_donor_2', 'bond_min_donor_3',
                                               'antibond_max_donor_1', 'antibond_max_donor_2', 'antibond_max_donor_3', 'bond_max_donor_1', 'bond_max_donor_2', 'bond_max_donor_3',
                                               'Rh_min_donor_bonding', 'Rh_min_donor_antibonding', 'Rh_max_donor_bonding', 'Rh_max_donor_antibonding']))
dataframe["filename"] = data['cas']
dataframe.index = data["cas"]
# for file in log_files:

for file in dataframe["filename"]:
    log_file = "final_log_files/" + file + ".log"
    print(f"doing {log_file}")
    elem, _ = mf.read_xyz("final_log_files/" + file + '.xyz')
    DFT = DFT_descriptors(log_file, len(elem), int(data[data['cas'] == file]['min_donor_atom_id']), int(data[data['cas'] == file]['max_donor_atom_id']), int(data[data['cas'] == file]['metal_id']) )
    LPO = DFT.donor_lone_pair_occupancy(data[data['cas'] == file]['min_donor_atom_type'].to_string()[-1], str(data[data['cas'] == file]['max_donor_atom_type'].to_string()[-1]))
    NBO_chrg = DFT.NBO_charge()
    mul_chrg = DFT.mulliken_charge()

    # try:
    DFT_b_ab = DFT.P_orbital_occupation()
    DFT_min_sort = np.sort(np.array(DFT_b_ab[0]))
    DFT_max_sort = np.sort(np.array(DFT_b_ab[1]))
    
    

    if len(DFT_max_sort) == 6:
        if len(DFT_min_sort) == 6:
            dataframe.loc[file] = pd.Series({'LP_occupancy_min'     : LPO[0], 
                                            'LP_occupancy_max'      : LPO[1],
                                            'NBO_charge_Rh'         : np.array(NBO_chrg)[0, 0],
                                            'NBO_charge_min'        : np.array(NBO_chrg)[1, 0],
                                            'NBO_charge_max'        : np.array(NBO_chrg)[2, 0],
                                            'Rh_mulliken_charge'    : mul_chrg[0],
                                            'min_mulliken_charge'   : mul_chrg[1],
                                            'max_mulliken_charge'   : mul_chrg[2],
                                            'antibond_min_donor_1'  : DFT_min_sort[0],
                                            'antibond_min_donor_2'  : DFT_min_sort[1],
                                            'antibond_min_donor_3'  : DFT_min_sort[2],
                                            'antibond_max_donor_1'  : DFT_max_sort[0],
                                            'antibond_max_donor_2'  : DFT_max_sort[1],
                                            'antibond_max_donor_3'  : DFT_max_sort[2],
                                            'bond_min_donor_1'      : DFT_min_sort[3],
                                            'bond_min_donor_2'      : DFT_min_sort[4],
                                            'bond_min_donor_3'      : DFT_min_sort[5],
                                            'bond_max_donor_1'      : DFT_max_sort[3],
                                            'bond_max_donor_2'      : DFT_max_sort[4],
                                            'bond_max_donor_3'      : DFT_max_sort[5],
                                            'Rh_min_donor_bonding'  : DFT_b_ab[2], 
                                            'Rh_min_donor_antibonding' : DFT_b_ab[3], 
                                            'Rh_max_donor_bonding'  : DFT_b_ab[4], 
                                            'Rh_max_donor_antibonding' : DFT_b_ab[5]
                                            })
        else:
            print(log_file + '_____1')

# return final_data_min, final_data_max, Rh_P_min_bonding, Rh_P_min_antibonding, Rh_P_max_bonding, Rh_P_max_antibonding
    else:
        print(log_file + '_____2')

        dataframe.loc[file] = pd.Series({'LP_occupancy_min'     : LPO[0], 
                                         'LP_occupancy_max'      : LPO[1],
                                            'NBO_charge_Rh'         : np.array(NBO_chrg)[0, 0],
                                            'NBO_charge_min'        : np.array(NBO_chrg)[1, 0],
                                            'NBO_charge_max'        : np.array(NBO_chrg)[2, 0],
                                            'Rh_mulliken_charge'    : mul_chrg[0],
                                            'min_mulliken_charge'   : mul_chrg[1],
                                            'max_mulliken_charge'   : mul_chrg[2],
                                            'antibond_min_donor_1'  : None,
                                            'antibond_min_donor_2'  : None,
                                            'antibond_min_donor_3'  : None,
                                            'antibond_max_donor_1'  : None,
                                            'antibond_max_donor_2'  : None,
                                            'antibond_max_donor_3'  : None,
                                            'bond_min_donor_1'      : None,
                                            'bond_min_donor_2'      : None,
                                            'bond_min_donor_3'      : None,
                                            'bond_max_donor_1'      : None,
                                            'bond_max_donor_2'      : None,
                                            'bond_max_donor_3'      : None
                                            })
    # except Exception:
                
    #     dataframe.loc[file] = pd.Series({'LP_occupancy_min'     : LPO[0], 
    #                                     'LP_occupancy_max'      : LPO[1],
    #                                     'NBO_charge_Rh'         : np.array(NBO_chrg)[0, 0],
    #                                     'NBO_charge_min'        : np.array(NBO_chrg)[1, 0],
    #                                     'NBO_charge_max'        : np.array(NBO_chrg)[2, 0],
    #                                     'Rh_mulliken_charge'    : mul_chrg[0],
    #                                     'min_mulliken_charge'   : mul_chrg[1],
    #                                     'max_mulliken_charge'   : mul_chrg[2],
    #                                     'antibond_min_donor_1'  : None,
    #                                     'antibond_min_donor_2'  : None,
    #                                     'antibond_min_donor_3'  : None,
    #                                     'antibond_max_donor_1'  : None,
    #                                     'antibond_max_donor_2'  : None,
    #                                     'antibond_max_donor_3'  : None,
    #                                     'bond_min_donor_1'      : None,
    #                                     'bond_min_donor_2'      : None,
    #                                     'bond_min_donor_3'      : None,
    #                                     'bond_max_donor_1'      : None,
    #                                     'bond_max_donor_2'      : None,
    #                                     'bond_max_donor_3'      : None,
    #                                     })

dataframe.to_excel("DFT_descriptors_test_i5.xlsx")
# save_data = np.array()

# print(pd.read_excel("mf_BD_Rh_test3.xlsx"))