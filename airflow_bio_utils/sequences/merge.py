from typing import Sequence, Tuple

from sqlitedict import SqliteDict

from airflow_bio_utils.filesystem import open_url
from airflow_bio_utils.logs import LOGS


def merge_sequences(
    paths: Sequence[str],
    output_path: str,
    deduplicate: bool = True,
    large_input: bool = False,
) -> Tuple[int, int]:
    """
    Merge multiple files with sequences into one single FASTA file.
    All sequences are deduplicated using sequence IDs.

    :param paths: Input paths to merge
    :param output_path: Path to the created output FASTA file
    :param deduplicate: Should remove sequences duplicates
    :param large_input: Use when dealing with +1GB input files and limited RAM memory to prevent
        Python out of memory issues. This option worksonly if deduplicate is set to True
        and changes set implementation to sqlite interface to store keys on disc.
    :return: Tuple of integers: (number_of_all_sequences, number_of_output_sequences)
    """
    LOGS.merge_sequences.info(f"Merge sequences {paths} -> {output_path}")

    # Capture set with all of the sequence IDs
    all_sequences = set()
    if large_input:
        all_sequences = SqliteDict("./merge.keystore.sqlite", autocommit=True)
    number_of_all_sequences = 0
    number_of_output_sequences = 0
    LOGS.merge_sequences.info(f"Open output file")
    with open_url(output_path, "w") as output_file:
        for input_path in paths:
            LOGS.merge_sequences.info(f"Try to open input file {input_path}")
            add_line = False
            #
            # Perform manual FASTA parsing
            # This is done, because we don't want to load the sequence itself, just IDs and that way it's faster
            #
            with open_url(input_path, "r") as input_file:
                LOGS.merge_sequences.info(f"Opened input file {input_path}")
                for line in LOGS.merge_sequences.progress(
                    input_file, message=f"Loading file {input_path}"
                ):
                    line = line.replace("\n", "")
                    if ">" in line:
                        number_of_all_sequences = number_of_all_sequences + 1
                        add_line = False
                        if line not in all_sequences or not deduplicate:
                            add_line = True
                            number_of_output_sequences = (
                                number_of_output_sequences + 1
                            )
                            if deduplicate:
                                if large_input:
                                    all_sequences[line] = True
                                else:
                                    all_sequences.add(line)
                    if add_line:
                        output_file.write(f"{line}\n")
        LOGS.merge_sequences.info(f"Loaded all files. Saving output file.")
    LOGS.merge_sequences.info(
        f"Merged sequences to FASTA file. Input count: {number_of_all_sequences}. "
        f"Output count: {number_of_output_sequences}"
    )
    LOGS.merge_sequences.info(
        "Duplication rate: "
        f"{round((number_of_all_sequences-number_of_output_sequences)/number_of_all_sequences*100, 2)}%"
    )

    if large_input:
        all_sequences.close()

    return number_of_all_sequences, number_of_output_sequences
