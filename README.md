# Airflow utils for biological sequences

![version_badge](https://s3.eu-central-1.amazonaws.com/pypi.covidgenomics.com/__public__/utils/badge_version.svg)

## Installation

To install this module please do:
```bash
   $ poetry install airflow-bio-utils@0.2.0
   # or
   $ pip install airflow-bio-utils==0.2.0
```

## Usage

### Airflow operators

Example usage:
```python
from airflow_bio_utils import SequenceMergeOperator

merge_fasta_files = SequenceMergeOperator(
    ['a.fasta', 'b.fasta'],
    'output.fasta',
)

```

### CLI

The package provides useful utilities for IO operations, transformations of biological data and reusable algorithms.
It also provides CLI interfaces for useful commands:
```bash
    $ poetry run create-random-sequences --help
    $ poetry run filter-sequences --help
    $ poetry run merge-sequences --help
```

## Development

To test the package please run:
```bash
    $ make test
```

Run linter with:
```bash
    $ make lint
```