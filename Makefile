init:
	conda env create

test:
	nosetests calc_props_test.py

.PHONY: init test
