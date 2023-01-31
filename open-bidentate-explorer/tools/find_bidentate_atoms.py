# openbabe-based method, example reading from a log file
import os, glob
from openbabel import openbabel
import numpy as np
from morfeus.io import read_cclib, get_xyz_string
from morfeus.utils import convert_elements


def get_bonded_atoms(xyz_string, atom_index):
    """Use openbabel's methods to find the coordinates of all atoms that are bonded to a given atom
    :param source_mol_file:
    :param atom_index:
    :return: numpy array of atoms bonded to a given atom
    """
    # initalize openbabel classes
    obconversion = openbabel.OBConversion()
    obconversion.SetInFormat("xyz")
    mol = openbabel.OBMol()
    obconversion.ReadString(mol, xyz_string)

    # make atom object for atom we want
    atom = mol.GetAtom(atom_index + 1)  # for obmol get functions indexing starts from 1

    index_list = []
    for neighbour_atom in openbabel.OBAtomAtomIter(atom):
        atomic_num = neighbour_atom.GetAtomicNum()
        index = neighbour_atom.GetIdx()
        index_list.append([atomic_num, index])

    # ToDo: convert obmol to rdkit and use methods there to determine correct bidentate cycle?
    # rdkit_mol = openbabel.OBMolToMol(mol)
    # or
    # rdkit_mol = Chem.MolFromMolBlock(conversion.WriteString(mol))
    return index_list


complexes_to_calc_descriptors = glob.glob(os.path.join(os.getcwd(), '*.log'))
for complex in complexes_to_calc_descriptors:
    elements, coordinates = read_cclib(complex)
    elements = convert_elements(elements, output='symbols')
    coordinates = np.array(coordinates).reshape(-1, len(elements), 3)
    xyz_string = ""
    for coord in coordinates:
        xyz_string = get_xyz_string(elements, coord)
    print(get_bonded_atoms(xyz_string, 38))

# RDkit and mace-dependent method
import mace
from rdkit import Chem
'''Find bidentate indices in metal-ligand complex based on smallest set of rings containing at least two mapped atoms'''


def check_if_at_least_two_mapped_atoms_in_ring(list_of_mapped_idxs, list_of_ring_idxs):
    return len(set(list_of_mapped_idxs) & set(list_of_ring_idxs)) >= 2
    # return sum(x in list_of_mapped_idxs for x in list_of_ring_idxs) >= 2


test_ligand = 'CC(C)[C@H]1CC[C@@H]([P:1]1C2=CC=CC=C2[N:1]3[C@H](CC[C@@H]3C(C)C)C(C)C)C(C)C'
substrate = 'CC#[N:1]'
ligands = [substrate, test_ligand]
CA = '[Ir+3]' # SMILES of central atom
geom = 'SP'

X = mace.ComplexFromLigands(ligands, CA, geom)
donor_atoms = [atom.GetIdx() for atom in X.mol3D.GetAtoms() if atom.GetAtomMapNum()]

Xs = X.GetStereomers(regime = 'all', dropEnantiomers = True)
first_conformer = Xs[0]

# all cycles that could possibly be the metal-bidentate cycle
bidentate_cycle_idxs = None
# the final bidentate cycle indices
bidentate_idxs = []
for i, X in enumerate(Xs):
    X.AddConformers(numConfs = 5)

    # get all mapped (donor) atom indices
    donor_atoms = [atom.GetIdx() for atom in X.mol3D.GetAtoms() if atom.GetAtomMapNum()]
    # smallest set of simple rings containing at least two mapped atoms
    atoms_of_smallest_rings = [list(idx) for idx in Chem.GetSymmSSSR(X.mol3D)]
    for list_idxs in atoms_of_smallest_rings:
        if check_if_at_least_two_mapped_atoms_in_ring(donor_atoms, list_idxs):
            # these are the same for each conformer, so I just reassign the variable each time
            # technically it would be enough to do this once for the first conformer, but I'm not sure if it's possible
            # ToDo, maybe only do all of this if i == 0? but check if it's possible to have different bidentate cycles for different conformers
            bidentate_cycle_idxs = list_idxs
            break

    # find metal atom and neighbours
    atomic_number_metal_center = Chem.MolFromSmiles(CA).GetAtomWithIdx(0).GetAtomicNum()
    metal_atom = [atom for atom in X.mol3D.GetAtoms() if atom.GetAtomicNum() == atomic_number_metal_center][0]
    metal_bonds = [bond for bond in metal_atom.GetBonds()]
    for bond in metal_bonds:
        idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (idx1 == metal_atom.GetIdx() and idx2 in bidentate_cycle_idxs):
            bidentate_idxs.append(idx2)
        elif (idx1 in bidentate_cycle_idxs and idx2 == metal_atom.GetIdx()):
            bidentate_idxs.append(idx1)

    print(set(bidentate_idxs))