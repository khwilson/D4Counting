# Code for The Number of D4 Fields Ordered by Conductor

This repository represents the code that was used to compute Tables 1 and 2 of [The number of D4 fields ordered by conductor](https://arxiv.org/pdf/1704.01729.pdf). It will be updated as the paper goes through the review process.

## Requirements

We mainly use Python 3.7 and sympy. We do use [Poetry](https://python-poetry.org/) as our dependency management system, so you'll first need to install that to package everything correctly. After that, though, you should be able to run:

```bash
poetry install
```

If you don't want to install `poetry`, then running `pip install sympy click && pip install -e .` will do the job.

## Usage

```bash
poetry run d4counting splitting    # Computes Table 1
poetry run d4counting expectation  # Computes Table 2
```

You can run the tests with

```bash
poetry run py.test
```

## How the code is arranged

The package `d4counting` contains all the relevant code. The CLI is driven by `cli.py`, whereas the details of the computations
can be found in `splitting_types.py` for Table 1 and `expected_number.py` for Table 2.

## License

MIT
