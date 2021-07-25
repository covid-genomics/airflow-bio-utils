import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "helpers"))

import biotite.sequence.io.fasta as fasta
from dag_tests import execute_test_task

from airflow_bio_utils import SequenceCountOperator, SequenceRandomOperator


def test_random():
    output_file_path = "./test1.fasta"
    sequences_count = 100
    min_length, max_length = 10, 100
    execute_test_task(
        SequenceRandomOperator,
        output=output_file_path,
        count=sequences_count,
        min_length=min_length,
        max_length=max_length,
    )
    count_output = execute_test_task(
        SequenceCountOperator,
        sequences=output_file_path,
    )
    assert count_output == sequences_count

    fasta_file = fasta.FastaFile.read(output_file_path)
    seq_count = 0
    for header, seq in fasta_file.items():
        assert min_length <= len(seq) <= max_length
        assert all(nucl in ["A", "C", "T", "G", "N"] for nucl in seq)
        seq_count += 1
    assert seq_count == sequences_count
