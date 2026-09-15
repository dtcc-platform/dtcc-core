UV ?= uv
COVERAGE_JSON := tests/coverage.json

.PHONY: install test coverage check-public-api verify-public-api check-public-api-strict verify-public-api-strict measure-import-time clean

install:
	$(UV) sync

test:
	cd tests && $(UV) run pytest

coverage:
	cd tests && $(UV) run pytest --maxfail=1 --disable-warnings --cov=dtcc_core \
	  --cov-report=term-missing --cov-report=json:coverage.json

check-public-api:
	$(UV) run python scripts/check_public_api_calls.py --package dtcc_core \
	  --coverage-file $(COVERAGE_JSON)

# Runs tests with coverage and verifies that all public API functions are executed
verify-public-api: coverage check-public-api

check-public-api-strict:
	$(UV) run python scripts/check_public_api_calls.py --strict --package dtcc_core \
	  --coverage-file $(COVERAGE_JSON)

verify-public-api-strict: coverage check-public-api-strict

# Reports the cost of `import dtcc_core` and the modules responsible for it
measure-import-time:
	$(PYTHON) scripts/measure_import_time.py

clean:
	rm -f $(COVERAGE_JSON)
	find . -name __pycache__ -type d -exec rm -rf {} + || true
	find . -name .pytest_cache -type d -exec rm -rf {} + || true
