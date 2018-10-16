install:
	pip install .

test:
	python -m unittest -v tests.test_mhcut
