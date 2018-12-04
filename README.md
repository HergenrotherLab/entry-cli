Entryway analysis
=================

[![Build Status](https://travis-ci.com/bdrown/entry-cli.svg?token=nW3s1WqrpGx8m3s2Pwu2&branch=master)](https://travis-ci.com/bdrown/entry-cli)

The eNTRy rules are a series of guidelines that can increase small-molecule
accumulation in gram-negative bacteria. A molecule is likely to accumulate if it
contains few rotatable bonds, low three dimensionality, and an ionizable
nitrogen.

The webapp [entry-way](http://www.entry-way.org) was created to aid in the application of these guidelines
by performing the necessary predictions of physiochemical properties. Although
freely available at entry-way.org, these same calculations can be performed
locally.

Dependencies
------------

Entryway relies on OpenBabel for handling chemical structures and conformer
generation and NumPy for globularity calculations. These dependencies are most
conveniently installed via [Conda](https://conda.io/docs/user-guide/install/index.html):

```bash
conda env create
conda activate entry-cli-env
```

Running calculations
--------------------

Although other file formats will be implemented shortly, molecules can most
be submitted as SMILES strings for now.

```bash
# Ampicillin
python calc_props.py -s "O=C(O)[C@@H]2N3C(=O)[C@@H](NC(=O)[C@@H](c1ccccc1)N)[C@H]3SC2(C)C"

# Deoxynybomycin
python calc_props.py -s "CC(C1=CC(C(C)=CC(N2C)=O)=C2C3=C1N4CO3)=CC4=O"
```
