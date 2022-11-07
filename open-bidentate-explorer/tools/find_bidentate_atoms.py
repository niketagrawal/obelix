import mace

test_ligand = 'CC(C)[C@H]1CC[C@@H]([P:1]1C2=CC=CC=C2[N:1]3[C@H](CC[C@@H]3C(C)C)C(C)C)C(C)C'
substrate = 'CC#[N:1]'
ligands = [test_ligand]
CA = '[Ir+3]' # SMILES of central atom
geom = 'SP'

X = mace.ComplexFromLigands(ligands, CA, geom)
Xs = X.GetStereomers(regime = 'all', dropEnantiomers = True)
print(len(Xs))

# calc = Calculator(descriptors)
first_conformer = Xs[0]
# find mapped atoms
donor_atoms = [atom.GetIdx() for atom in first_conformer.mol3D.GetAtoms() if atom.GetAtomMapNum()]
print(donor_atoms[-2:][0])


list_of_complexes = []
list_of_complexes2 = []
for i, X in enumerate(Xs):
    X.AddConformers(numConfs = 5)
    donor_atoms = [atom.GetIdx() for atom in X.mol3D.GetAtoms() if atom.GetAtomMapNum()]
    # for some reason the last 2 items are always bidentate ligand donor atoms
    print(donor_atoms)

