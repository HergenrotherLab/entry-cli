from nose.tools import *
from openbabel import openbabel
import os
from .context import calc_props

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def test_smiles_benzene():
    mol = calc_props.smiles_to_ob("c1ccccc1")
    assert(isinstance(mol, openbabel.OBMol))
    assert_equals(mol.NumAtoms(), 12)

def test_with_empty_args():
    assert_raises(SystemExit, calc_props.parse_args, [])

def test_with_correct_smiles_options():
    try:
        calc_props.parse_args(["-s", "C1C3CC2CC(CC1C2)C3"])
    except SystemExit:
        assert False

def test_with_smiles_out_with_path():
    try:
        calc_props.parse_args(["-s", "C1C3CC2CC(CC1C2)C3", "-o", "out.csv"])
    except SystemExit:
        assert False

def test_with_correct_batch_options():
    try:
        calc_props.parse_args(["-b", os.path.join(THIS_DIR, "data/triphenylphosphine.mol")])
    except SystemExit:
        assert False

def test_with_single_and_batch():
    assert_raises(SystemExit, calc_props.parse_args, ["-s", "C1C3CC2CC(CC1C2)C3", "-b", os.path.join(THIS_DIR, "data/triphenylphosphine.mol")])