# Convenience targets. CI runs the same commands directly, so this is for
# local use and does not replace the workflows.
.DEFAULT_GOAL := help
PYTHON ?= python3
VENV ?= .venv
ACTIVATE = . $(VENV)/bin/activate

.PHONY: help install dev test lint shellcheck docs docs-serve build check \
        version dist clean

help:
	@echo "creTarget development targets"
	@echo
	@echo "  make install     Create .venv and install the package with dev extras"
	@echo "  make test        Run the Python test suite"
	@echo "  make lint        Run ruff"
	@echo "  make shellcheck  Run ShellCheck on the shell scripts"
	@echo "  make docs        Build the documentation site strictly"
	@echo "  make docs-serve  Serve the documentation at http://127.0.0.1:8000"
	@echo "  make version     Verify that all version pins agree"
	@echo "  make dist        Build and check the sdist and wheel"
	@echo "  make check       Run lint, tests, docs, and the version check"
	@echo "  make clean       Remove build, test, and documentation artifacts"

$(VENV)/bin/activate:
	$(PYTHON) -m venv $(VENV)
	$(ACTIVATE) && python -m pip install --upgrade pip

install: $(VENV)/bin/activate
	$(ACTIVATE) && pip install -e ".[dev,docs]"

dev: install

test:
	$(ACTIVATE) && pytest

lint:
	$(ACTIVATE) && ruff check src tests

shellcheck:
	shellcheck --severity=warning assets/*.sh bump_version.sh scripts/*.sh

docs:
	$(ACTIVATE) && mkdocs build --strict

docs-serve:
	$(ACTIVATE) && mkdocs serve

version:
	bash bump_version.sh --check

dist:
	$(ACTIVATE) && python -m build && twine check dist/*

check: lint test docs version

clean:
	rm -rf build dist site .pytest_cache .ruff_cache htmlcov .coverage
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name '*.egg-info' -prune -exec rm -rf {} +
