from nose.tools import *
from contextlib import contextmanager
import os
import sys
if (sys.version_info > (3, 0)):
    from io import StringIO
else:
    from io import BytesIO as StringIO
from .context import calc_props

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def build_molecule():
    properties = {
        'smiles': 'CC(C1=CC(C(C)=CC(N2C)=O)=C2C3=C1N4CO3)=CC4=O',
        'formula': 'C16H14N2O3',
        'molwt': 282.293960,
        'rb': 0,
        'glob': 0.023820,
        'pbf': 0.204747,
        'primary_amine': False
    }
    return properties

def test_write_csv():
    mols = [build_molecule()] * 10
    output = os.path.join(THIS_DIR, 'data/write_csv_test.csv')
    calc_props.write_csv(mols, output)

    # Read file that was just written and desired output
    with open(output, 'r') as result_file:
        result = result_file.read()
    with open(os.path.join(THIS_DIR, 'data/write_csv_test_example.csv'), 'r') as example_file:
        example = example_file.read()

    # Compare two files. If different, print results.
    try:
        assert_multi_line_equal(example, result)
    except:
        print(result)
        assert False
    finally:
        os.remove(output)

@contextmanager
def stdout_redirector(stream):
    old_stdout = sys.stdout
    sys.stdout = stream
    try:
        yield
    finally:
        sys.stdout = old_stdout

def test_report_properties():
    f = StringIO()
    with stdout_redirector(f):
        calc_props.report_properties(build_molecule())
    terms = ['Properties', 'Mol. Wt.', '282.29', 'Formula', 'RB', 'Glob', 'PBF', 'primary_amine']
    assert all(term in f.getvalue() for term in terms)
