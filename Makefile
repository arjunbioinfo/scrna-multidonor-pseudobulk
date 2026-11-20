.PHONY: setup run ci test clean

CORES ?= 2

setup:
	pip install -r requirements.txt

run:
	snakemake --cores $(CORES) --configfile config/config.yaml

ci:
	snakemake --cores $(CORES) --configfile config/ci.yaml

test:
	pytest -v tests/

clean:
	rm -rf results/*.tsv results/*.h5ad figures/*.png data/raw/*
