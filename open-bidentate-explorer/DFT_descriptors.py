import re, glob


class DFT_descriptors:
    def __init__(self, log_file, nr_of_atoms, min_donor, max_donor, metal):
        
        self.log_file = log_file # name of log file
        self.nr_of_atoms = nr_of_atoms # nr of atoms in molecule
        self.metal = metal  # metal index
        self.min_donor = min_donor # min donor index
        self.max_donor = max_donor # max donor index
        
        
    def donor_lone_pair_occupancy(self, atom_type_min, atom_type_max):
        print("Extracting donor atoms lone pair orbital occupancy:")
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
        
        return orbital_dictionary[self.min_donor + 1].values(), orbital_dictionary[self.max_donor + 1]

    
    def NBO_charge(self):
        print("Extracting donor atoms and metal natural charge:")
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
        print(charge_final_data)

        for charge in charge_final_data:
            charge_dictionary[int(charge[1])] = {float(charge[2])}
  
        # divided by 2 since 2 electrons form the lone pair and the 3s orbital
        
        return charge_dictionary[self.min_donor + 1], charge_dictionary[self.max_donor + 1]

    def mulliken_charge(self):
        print("Extracting Mulliken charge:")
        with open(self.log_file) as file: # Use file to refer to the file object
            data = file.readlines()
            for line_index, line in enumerate(data):
                if "Mulliken charges" in line:        
                    Mulliken_raw_data = data[line_index + 2:line_index + self.nr_of_atoms]
                    break
        file.close()
        
        mulliken_final_data = []
        
        for raw_charge in Mulliken_raw_data:
            raw_charge = raw_charge.strip("\n")
            raw_charge = raw_charge.strip("\r")
            raw_charge = raw_charge.split()
            mulliken_final_data.append(float(raw_charge[-1]))
    
        return mulliken_final_data[self.metal], mulliken_final_data[self.min_donor], mulliken_final_data[self.max_donor]

log_file = "158923-09-2.log"
DFT = DFT_descriptors(log_file, 80, 21, 26, 43)
NBO = DFT.NBO_charge()
print(NBO)