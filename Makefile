.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=genome_sampler

install: all
	$(PYTHON) -m pip install -v .

dev: all
	pip install -e .

clean: distclean

distclean: ;

html:
	cd docs && q2doc autodoc .
	cd docs && jupyter book build --html
	cp -r docs/data/ docs/_build/html/data/

serve:
	npx serve docs/_build/html/ -p 4000

clean:
	rm -rf docs/_build/html/
