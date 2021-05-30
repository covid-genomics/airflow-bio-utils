from typing import List, Optional, Sequence

from Bio.SeqRecord import SeqRecord

from airflow_bio_utils.filesystem import open_url
from airflow_bio_utils.logs import LOGS
from airflow_bio_utils.sequences.transformations import \
    perform_on_each_sequence

# Default values for parameters for filter_sequences
DEFAULT_MIN_SEQ_LEN = 20000
DEFAULT_MAX_SEQ_LEN = 1500000
DEFAULT_ACCEPTED_SEQUENCE_SYMBOLS = ["A", "C", "T", "G"]


def filter_sequences(
    input_paths: Sequence[str],
    min_seq_len: int = DEFAULT_MIN_SEQ_LEN,
    max_seq_len: int = DEFAULT_MAX_SEQ_LEN,
    accepted_symbols: Optional[List[str]] = None,
    output_paths: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    Filter file with sequences.
    Function takes list of input files and for each performs filtering of the sequences.
    A new file is generated.
    If the output_paths is not None, then it should have the same length as input_paths.
    In that case from a file input_paths[i] a file output_paths[i] is generated.
    If output_paths is None (by default), then new file will have the same name, but with ".filtered" appended.

    :param input_paths: Paths to input files
    :param min_seq_len: Minimum sequence length
    :param max_seq_len: Maximum sequence length
    :param accepted_symbols: Accepted symbols in the sequence, by default (A, C, T, G)
    :param output_paths: Optional list of output paths
    :return: List of all generated files
    """
    if accepted_symbols is None:
        accepted_symbols = DEFAULT_ACCEPTED_SEQUENCE_SYMBOLS

    all_output_paths = []
    for (i, path) in enumerate(input_paths):
        valid_genomes = 0
        output_filename = f"{path}.filtered"
        if output_paths is not None:
            output_filename = output_paths[i]
        with open_url(output_filename, "w") as f:

            def filter_and_save_genome(record: SeqRecord):
                """
                Process single input record

                :param record: Input SeqRecord object
                """
                nonlocal valid_genomes
                filtered = None
                if len(record) >= min_seq_len <= max_seq_len:
                    if all((symbol in accepted_symbols) for symbol in record):
                        filtered = record.upper()
                if filtered is not None:
                    valid_genomes += 1
                    print(f'>{record.id.replace(" ", "").rstrip()}', file=f)
                    print(filtered, file=f)

            all_genomes = perform_on_each_sequence(
                filter_and_save_genome, path
            )
        LOGS.filter_sequences.info(
            "The total number of valid genomes is {}, {}% of all genomes passed".format(
                valid_genomes, round(valid_genomes / all_genomes * 100)
            )
        )
        all_output_paths.append(output_filename)
    return all_output_paths
