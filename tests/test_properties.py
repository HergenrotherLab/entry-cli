from nose.tools import *
from openbabel import openbabel
from openbabel import pybel
import os
from .context import calc_props

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def test_smiles_benzene():
    mol = calc_props.smiles_to_ob("c1ccccc1")
    assert(isinstance(mol, openbabel.OBMol))
    assert_equals(mol.NumAtoms(), 12)

def test_rb_basic():
    # DNM
    mol = calc_props.smiles_to_ob("CC(C1=CC(C(C)=CC(N2C)=O)=C2C3=C1N4CO3)=CC4=O")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 0)

    # Ribocil C
    mol = calc_props.smiles_to_ob("C1CC(CN(C1)CC2=CN(C=N2)C3=NC=CC=N3)C4=NC(=O)C=C(N4)C5=CC=CS5")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 5)

    # Triphenylphosphine
    mol = calc_props.smiles_to_ob("C1(P(C2=CC=CC=C2)C3=CC=CC=C3)=CC=CC=C1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 3)

def test_rb_alcohol():
    # n-butanol
    mol = calc_props.smiles_to_ob("CCCCO")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 2)

def test_rb_amine():
    # n-butylamine
    mol = calc_props.smiles_to_ob("CCCCN")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 2)

def test_rb_amide():
    # Ala-Ala
    mol = calc_props.smiles_to_ob("[H]N[C@H](C(N[C@H](C(O)=O)C)=O)C")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 3)

def test_rb_ketene():
    # pent-1-en-1-one
    mol = calc_props.smiles_to_ob("[H]C(CCC)=C=O")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 2)

def test_rb_allene():
    # 3-methylocta-3,4-diene
    mol = calc_props.smiles_to_ob("[H]C(CCC)=C=C(C)CC")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 3)

def test_rb_alkyne():
    # but-1-yn-1-ylbenzene
    mol = calc_props.smiles_to_ob("CCC#CC1=CC=CC=C1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 1)

def test_rb_symmetric_alkyne():
    # hex-3-yne
    mol = calc_props.smiles_to_ob("CCC#CCC")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 1)

def test_rb_cyclohexane_alkyne():
    # but-1-yn-1-ylcyclohexane
    mol = calc_props.smiles_to_ob("CCC#CC1CCCCC1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 1)

def test_rb_cyclohexene_alkyne():
    # 1-(but-1-yn-1-yl)cyclohex-1-ene
    mol = calc_props.smiles_to_ob("CCC#CC1=CCCCC1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 1)

def test_rb_alkene():
    # (E)-but-1-en-1-ylbenzene
    mol = calc_props.smiles_to_ob("CC/C=C/C1=CC=CC=C1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 2)

def test_rb_nitrile():
    # (E)-5-phenylpent-4-enenitrile
    mol = calc_props.smiles_to_ob("N#CCC/C=C/C1=CC=CC=C1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 3)

def test_rb_azide():
    # (E)-(4-azidobut-1-en-1-yl)benzene
    mol = calc_props.smiles_to_ob("[N-]=[N+]=NCC/C=C/C1=CC=CC=C1")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 4)

def test_rb_ester():
    # phenyl butyrate
    mol = calc_props.smiles_to_ob("CCCC(OC1=CC=CC=C1)=O")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 4)

def test_rb_ketone():
    # 1-phenylpentan-2-one
    mol = calc_props.smiles_to_ob("CCCC(CC1=CC=CC=C1)=O")
    pymol = pybel.Molecule(mol)
    assert_equals(calc_props.rotatable_bonds(pymol), 4)

def test_pbf():
    obmol = openbabel.OBMol()
    obConv = openbabel.OBConversion()
    obConv.SetInFormat("mol")
    obConv.ReadFile(obmol, os.path.join(THIS_DIR, "data/triphenylphosphine.mol"))
    pymol = pybel.Molecule(obmol)
    points = calc_props.get_atom_coords(pymol)
    assert_almost_equal(calc_props.calc_pbf(points), 1.0072297, 6, 1)

def test_glob():
    obmol = openbabel.OBMol()
    obConv = openbabel.OBConversion()
    obConv.SetInFormat("mol")
    obConv.ReadFile(obmol, os.path.join(THIS_DIR, "data/triphenylphosphine.mol"))
    pymol = pybel.Molecule(obmol)
    points = calc_props.get_atom_coords(pymol)
    assert_almost_equal(calc_props.calc_glob(points), 0.245503, 6, 1)

def test_glob_benzene():
    mol = calc_props.smiles_to_ob("c1ccccc1")
    properties = calc_props.average_properties(mol)
    assert_almost_equal(properties['glob'], 0, 2, 1)

def test_adamantane():
    mol = calc_props.smiles_to_ob("C1C3CC2CC(CC1C2)C3")
    properties = calc_props.average_properties(mol)
    assert_almost_equal(properties['glob'], 1, 2, 1)

def test_cipro():
    mol = calc_props.smiles_to_ob("O=C1C(C(O)=O)=CN(C2CC2)C3=CC(N4CCNCC4)=C(F)C=C31")
    properties = calc_props.average_properties(mol)
    assert_almost_equal(properties['glob'], 0.04, 2, 1)

def test_dnm():
    mol = calc_props.smiles_to_ob("CC(C1=CC(C(C)=CC(N2C)=O)=C2C3=C1N4CO3)=CC4=O")
    properties = calc_props.average_properties(mol)
    assert_almost_equal(properties['glob'], 0.024, 2, 1)

def test_primary_amine_serine():
    mol = calc_props.smiles_to_ob("C(C(C(=O)O)N)O")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), True)

def test_primary_amine_glycine():
    mol = calc_props.smiles_to_ob("C(C(=O)O)N")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), True)

def test_primary_amine_dimethylamine():
    mol = calc_props.smiles_to_ob("CNC")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), False)

def test_primary_amine_aniline():
    mol = calc_props.smiles_to_ob("Nc1ccccc1")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), False)

def test_primary_amine_imine():
    mol = calc_props.smiles_to_ob("CC1=CC(=CC=C1)N=CC2=CC=CC=C2")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), False)

def test_primary_amine_1_benzylguanidine():
    mol = calc_props.smiles_to_ob("NC(NCC1=CC=CC=C1)=N")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), False)

def test_primary_amine_methylethyl_ketone_oxime():
    mol = calc_props.smiles_to_ob("CC/C(C)=N/O")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), False)

def test_primary_amine_benzenesulfonamide():
    mol = calc_props.smiles_to_ob("O=S(C1=CC=CC=C1)(N)=O")
    pymol = pybel.Molecule(mol)
    primary_amine_smarts = calc_props.FUNCTIONAL_GROUP_TO_SMARTS["primary_amine"]
    assert_equals(calc_props.has_functional_group(pymol, primary_amine_smarts), False)
