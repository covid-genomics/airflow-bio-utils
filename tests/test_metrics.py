import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "helpers"))

from dag_tests import execute_test_task

from airflow_bio_utils import SequenceMetricsOperator, SequenceRandomOperator


def test_metrics():
    test_file_path = "./test_example.fasta"
    output_path = "./metrics_output"
    sequences_count = 100
    min_length, max_length = 10, 100
    execute_test_task(
        SequenceRandomOperator,
        output=test_file_path,
        count=sequences_count,
        min_length=min_length,
        max_length=max_length,
    )
    execute_test_task(
        SequenceMetricsOperator,
        input_path=test_file_path,
        output_path=output_path,
    )


