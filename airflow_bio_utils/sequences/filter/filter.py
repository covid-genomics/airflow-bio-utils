from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from Bio.SeqRecord import SeqRecord

from airflow_bio_utils.filesystem import open_url
from airflow_bio_utils.logs import LOGS
from airflow_bio_utils.sequences.filter.condition import FilterCondition
from airflow_bio_utils.sequences.transformations import \
    perform_on_each_sequence


@dataclass(frozen=True)
class FilterResultsMetadata:
    output_path: str
    input_path: str
    rejected_statistics: Dict[str, int]
    accepted_sequences: int
    input_sequences: int
    percent_of_accepted_sequences: float


def filter_sequences(
    input_paths: Sequence[str],
    filters: List[FilterCondition],
    output_paths: Optional[Sequence[str]] = None,
) -> List[FilterResultsMetadata]:
    """
    Filter file with sequences.
    Function takes list of input files and for each performs filtering of the sequences.
    A new file is generated.
    If the output_paths is not None, then it should have the same length as input_paths.
    In that case from a file input_paths[i] a file output_paths[i] is generated.
    If output_paths is None (by default), then new file will have the same name, but with ".filtered" appended.

    :param input_paths: Paths to input files
    :param output_paths: Optional list of output paths
    :return: List of all generated files
    """

    output_metadata: List[FilterResultsMetadata] = []
    for (i, path) in enumerate(input_paths):
        valid_genomes = 0
        rejected_reasons: Dict[str, int] = defaultdict(int)

        output_filename = f"{path}.filtered"
        if output_paths is not None:
            output_filename = output_paths[i]
        with open_url(output_filename, "w") as f:

            def filter_and_save_genome(record: SeqRecord):
                """
                Process single input record

                :param record: Input SeqRecord object
                """
                nonlocal valid_genomes, rejected_reasons
                reject_reason: Optional[str] = None
                for filter in filters:
                    if not filter.check(record):
                        reject_reason = filter.__class__.__name__

                if reject_reason is None:
                    valid_genomes += 1
                    print(record.format("fasta"), file=f, end="")
                else:
                    rejected_reasons[reject_reason] += 1

            all_genomes = perform_on_each_sequence(
                filter_and_save_genome, path
            )
        percent_of_valid = round(valid_genomes / all_genomes * 100)
        LOGS.filter_sequences.info(
            "The total number of valid genomes is {}, {}% of all genomes passed".format(
                valid_genomes, percent_of_valid
            )
        )
        output_metadata.append(
            FilterResultsMetadata(
                path,
                output_filename,
                rejected_reasons,
                valid_genomes,
                all_genomes,
                percent_of_valid,
            )
        )
    return output_metadata
