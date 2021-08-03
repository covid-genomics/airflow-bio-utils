import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "helpers"))

import os
from typing import Tuple

from dag_tests import execute_test_task
from dvc_fs.management.create_dvc_repo_github import \
    create_github_dvc_temporary_repo_with_s3

from airflow_bio_utils import (FilterByLength, FilterSymbols,
                               SequenceFilterOperator, SequenceMergeOperator,
                               SequenceRandomOperator)


def test_dvc_integration():
    repo = create_github_dvc_temporary_repo_with_s3(
        "covid-genomics", "temporary_dvc_repo"
    )
    with repo as fs:
        dvc_url = f"dvc://{os.environ['DVC_GITHUB_REPO_TOKEN']}@github.com/{repo.owner}/{repo.repo_name}"
        input_counts = (21, 62)

        execute_test_task(
            SequenceRandomOperator,
            output=f"{dvc_url}//input1.fasta",
            count=input_counts[0],
            min_length=9000,
            max_length=11000,
        )

        execute_test_task(
            SequenceRandomOperator,
            output=f"{dvc_url}//input2.fasta",
            count=input_counts[1],
            min_length=9000,
            max_length=11000,
        )

        merge_output: Tuple[int, int] = execute_test_task(
            SequenceMergeOperator,
            paths=[f"{dvc_url}//input1.fasta", f"{dvc_url}//input2.fasta"],
            output_path=f"{dvc_url}//out.fasta",
        )
        assert merge_output[0] == sum(input_counts)
        assert 0 <= merge_output[1] <= merge_output[0]

        execute_test_task(
            SequenceFilterOperator,
            input_paths=f"{dvc_url}//out.fasta",
            output_paths=f"{dvc_url}//out.filtered.fasta",
            filters=[
                FilterSymbols(),
                FilterByLength(
                    min_length=10000,
                    max_length=10700,
                ),
            ],
        )
