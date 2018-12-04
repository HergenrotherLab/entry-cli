from nose.tools import assert_equals
from nose.tools import assert_almost_equal
import openbabel
import calc_props


def test_smiles_benzene():
    mol = calc_props.smiles_to_ob("c1ccccc1")
    assert(isinstance(mol, openbabel.OBMol))
    assert_equals(mol.NumAtoms(), 12)

def test_confab():
    pass

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
