# Airflow utils for biological sequences

![version_badge](https://s3.eu-central-1.amazonaws.com/pypi.covidgenomics.com/__public__/utils/badge_version.svg)

## Installation

To install this module please do:
```bash
   $ poetry install airflow-bio-utils@0.6.0
   # or
   $ pip install airflow-bio-utils==0.6.0
```

## Usage

### Airflow operators

#### Sequence merge

`SequenceMergeOperator` takes input files and merges the sequences together.
It provides efficient antiduplication mechanism to remove duplicate IDs.

Merging sequences from local paths:
```python
from airflow_bio_utils import SequenceMergeOperator

merge_fasta_files = SequenceMergeOperator(
    ['./a.fasta', 'b.fasta'],
    'output.fasta',
)
```

Merging sequences from remote locations:
```python
from airflow_bio_utils import SequenceMergeOperator

merge_fasta_files = SequenceMergeOperator(
    [
        "s3://bucket1//folder1/folder2/input1.fasta",
        "s3://bucket1//folder1/folder2/input2.fasta",
        "s3://other_bucket//some_folder/fasta_file.fasta",
    ],
    "s3://output_bucket/output.fasta",
)
```

#### Generating random sequences

`SequenceRandomOperator` can generate nucleotide or protein (translated) sequences with given lengths in given amounts.

```python
from airflow_bio_utils import SequenceRandomOperator

random_sequence_generator_task = SequenceRandomOperator(
    output="dvc://ssh@github.com/covid-genomics/data-artifacts//random.fasta",
    count=130,
    min_length=1000,
    max_length=5000,
)
```

#### Sequences metrics

`SequenceMetricsOperator` allows you to measure qualities of `fastq` metrics, trace GC content, kmers and get other statistics.
By default it generated PDF report.

```python
from airflow_bio_utils import SequenceMetricsOperator

SequenceMetricsOperator(
    input_path="some_input.fasta",
    output_path="dvc://<PAT>@github.com/covid-genomics/data-artifacts//output_metrics",
)
```

#### Sequence filtering

`SequenceFilterOperator` allows you to remove sequences based on simple criteria like length or nucleotide contents.

```python
from airflow_bio_utils import SequenceFilterOperator

SequenceFilterOperator(
    input_paths=["a.fasta", "b.fasta"],
    min_seq_len=1000,
    accepted_symbols=["A", "C", "T", "G", "N"],
    output_paths=["a_filtered.fasta", "b_filtered.fasta"],
)
```

#### Sequence for each

`SequenceForEachOperator` takes an action and executes it for each sequence.
The action can return None or sequence and if you specify `output_path` parameter,
all results will be saved into new fasta file.
This makes `SequenceForEachOperator` foreach or map/filter operator.

```python
from airflow_bio_utils import SequenceForEachOperator

task = SequenceForEachOperator(
    sequences="sequences.fasta",
    action=lambda seq: seq,
    output_path="output.fasta",
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