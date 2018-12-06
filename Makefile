init:
	conda env create

test:
	nosetests tests

.PHONY: init test
