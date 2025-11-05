.PHONY: help install install-dev test lint format type-check clean build docs

help: ## Show this help message
	@echo "LeishDomains Development Commands"
	@echo "================================="
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install the package
	pip install -e .

install-dev: ## Install with development dependencies
	pip install -e ".[dev]"
	pre-commit install

test: ## Run tests
	pytest tests/ -v --cov=core --cov=data --cov=analysis --cov=cli --cov-report=term-missing

test-html: ## Run tests with HTML coverage report
	pytest tests/ -v --cov=core --cov=data --cov=analysis --cov=cli --cov-report=html

lint: ## Run linting
	flake8 core/ data/ analysis/ cli/ tests/ --max-line-length=127 --max-complexity=10

format: ## Format code with black
	black core/ data/ analysis/ cli/ tests/ main.py

type-check: ## Run type checking with mypy
	mypy core/ data/ analysis/ cli/ --ignore-missing-imports

check: lint type-check test ## Run all checks

pre-commit: ## Run pre-commit hooks on all files
	pre-commit run --all-files

pre-commit-install: ## Install pre-commit hooks
	pre-commit install

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean ## Build package
	python -m build

docs: ## Build documentation
	cd docs && sphinx-build -b html . _build/html

docs-serve: docs ## Build and serve documentation
	cd docs/_build/html && python -m http.server 8000

run-gui: ## Run GUI application
	python main.py

run-cli: ## Run CLI application (with help)
	python main.py --help

example: ## Run usage examples
	python examples/usage_examples.py

all: clean install-dev check build ## Run everything
