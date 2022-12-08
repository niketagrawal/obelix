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