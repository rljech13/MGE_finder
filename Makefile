.PHONY: test-finder-pipeline test-atlas test-cluster test

# Override: make MGE_PYTHON=python3 …
MGE_PYTHON ?= python3

test-finder-pipeline:
	@$(MGE_PYTHON) -m pytest tests/finder_pipeline/ -q

test-atlas:
	@$(MGE_PYTHON) -m unittest tests/atlas_pipeline/test_classify.py -q

test-cluster:
	@$(MGE_PYTHON) -m unittest tests/cluster_pipeline/test_cluster_lib.py -q

test: test-finder-pipeline test-atlas test-cluster
