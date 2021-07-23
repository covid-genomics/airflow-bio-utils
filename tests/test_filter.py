import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))

import biotite.sequence.io.fasta as fasta
from typing import List

from dag_tests import execute_test_task
from airflow_bio_utils import SequenceRandomOperator, SequenceFilterOperator, FilterResultsMetadata, FilterSymbols, FilterByLength


def test_filter():
    output_file_path = "./test1.fasta"
    filtered_file_path = "./test1_filtered.fasta"
    sequences_count = 500
    min_length, max_length = 10, 100
    execute_test_task(
        SequenceRandomOperator,
        output=output_file_path,
        count=sequences_count,
        min_length=min_length-5,
        max_length=max_length+5,
    )
    filter_metadata: List[FilterResultsMetadata] = execute_test_task(
        SequenceFilterOperator,
        input_paths=output_file_path,
        output_paths=filtered_file_path,
        filters=[
            FilterSymbols(),
            FilterByLength(
                min_length=min_length,
                max_length=max_length,
            ),
        ],
    )
    fasta_file = fasta.FastaFile.read(filtered_file_path)
    seq_count = 0
    for header, seq in fasta_file.items():
        assert min_length <= len(seq) <= max_length
        assert all(nucl in ["A", "C", "T", "G", "N"] for nucl in seq)
        seq_count += 1
    assert seq_count == filter_metadata[0].accepted_sequences

if __name__ == '__main__':
    test_filter()
