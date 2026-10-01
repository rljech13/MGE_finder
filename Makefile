.PHONY: test-finder-pipeline

# Override: make MGE_PYTHON=python3 …
MGE_PYTHON ?= python3

test-finder-pipeline:
	@$(MGE_PYTHON) -m pytest tests/finder_pipeline/ -q
