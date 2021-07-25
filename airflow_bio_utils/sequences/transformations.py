import collections
from typing import Callable, Iterator, Sequence, Union, Optional, List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from airflow_bio_utils.filesystem import open_url

from airflow_bio_utils.file_operations import get_file_lines_count
from airflow_bio_utils.logs import LOGS

# Generic sequence that can be anything that is string or biopython sequence-like
GenericSequence = Union[Seq, SeqRecord, str]
# List of generic sequences or a single one
GenericSequences = Union[
    GenericSequence, Sequence[GenericSequence], Iterator[GenericSequence]
]


def normalize_single_sequence(
    sequence: GenericSequence,
    i: int,
) -> SeqRecord:
    """
    Take a single sequence and normalize it's form to Biopython SeqRecord object.
    If the sequence has no ID the function uses i parameter to assign numeric ID to the sequence record.

    Supported input sequence types are:
        - Seq object
        - SeqRecord object
        - string

    :param sequence: Input sequence of supported type
    :param i: Index of that sequence
    :return: SeqRecord object
    """
    if isinstance(sequence, SeqRecord):
        return sequence
    elif isinstance(sequence, Seq):
        return SeqRecord(sequence, str(i))
    elif isinstance(sequence, str):
        return SeqRecord(Seq(sequence), str(i))


def convert_to_biopython_sequences(
    sequences: GenericSequences,
) -> Sequence[SeqRecord]:
    """
    Take a bunch of sequences and normalize their form to list of Biopython SeqRecord objects.

    Supported input sequences types are:
        - List of Seq objects (*)
        - List of SeqRecord objects
        - List of strings (*)
        - Seq object (**)
        - SeqRecord object (**)
        - string (**)

    If the input sequences does not provide information about the sequence id (marked with asterisks *) then
    the list index (starting at 0) will be used for sequence ID.

    If the input sequences are in fact only a single object then the integer value of 0 will be used for sequence ID.
    (marked with double asterisks **)

    :param sequences: Some input of supported type
    :return: List of Biopython SeqRecord objects
    """
    if isinstance(sequences, SeqRecord):
        return [sequences]
    elif isinstance(sequences, Seq):
        return [SeqRecord(sequences, 0)]
    elif isinstance(sequences, str):
        return [SeqRecord(Seq(sequences), 0)]
    elif isinstance(sequences, collections.Iterable):
        return [
            normalize_single_sequence(seq, i)
            for i, seq in enumerate(sequences)
        ]
    raise Exception("Invalid input type for normalize_sequences() was given.")


def detect_sequences_format(
    file_path: str, default_file_format: str = "fasta"
) -> str:
    """
    For a given file path try to determine what sequences format it contains.
    For a list of BioPython supported formats please see: https://biopython.org/wiki/SeqIO

    :param file_path: Input file path
    :param default_file_format: Default format used as a fallback
    :return: Name of the format of the file (compatible with BioPython file formats)
    """
    if file_path.endswith(".fasta"):
        return "fasta"
    return default_file_format


def load_sequences(
    sequences: GenericSequences, default_file_format: str = "fasta"
) -> Sequence[SeqRecord]:
    """
    Load sequences from the file.
    - If the parameter is string, then it's treated as a path to the input file.
        Function automatically detects file format with detect_sequences_format().
    - If the parameter is not a string, then it's converted to SeqRecords using convert_to_biopython_sequences().

    :param sequences: Input sequences
    :param default_file_format: Default file format to be used
        Refer to the detect_sequences_format() for more details
    :return: Loaded SeqRecord objects
    """
    if isinstance(sequences, str):
        default_file_format = detect_sequences_format(
            sequences, default_file_format
        )
        return SeqIO.parse(open_url(sequences).local_file_path, default_file_format)
    return convert_to_biopython_sequences(sequences)


def count_sequences(
    sequences: Union[str, Sequence[SeqRecord]],
    default_file_format: str = "fasta",
) -> int:
    """
    Counts the sequences.
    - If the parameter is string, then it's treated as a path to the input file.
    - If the parameter is a sequence of SeqRecord objects, then len() is called

    :param sequences: Path to the input file or sequence of SeqRecord objects
    :param default_file_format: Default file format to be used when loading input file
        Refer to the detect_sequences_format() for more details
    :return: Number of sequences
    """
    if isinstance(sequences, str):
        default_file_format = detect_sequences_format(
            sequences, default_file_format
        )
        if default_file_format == "fasta":
            return get_file_lines_count(sequences, ">")
        else:
            raise Exception(
                f"Unsupported format for count_sequences: "
                f"Input is file path and format is {default_file_format}"
            )
    return len(sequences)


def perform_on_each_sequence(
    action: Callable[[SeqRecord], Optional[SeqRecord]],
    sequences: GenericSequences,
    default_file_format: str = "fasta",
    output_path: Optional[str] = None,
) -> int:
    """
    Perform operation on each sequence record.
    Function accepts as input various data types:
    - If sequences paremeter is a string, it's treated as a path to the input file
      The input file is then loaded
    - If sequences parameter is not a string, then it's converted to sequence of SeqRecord objects using
      convert_to_biopython_sequences() function

    :param action: Function accepting the sequence header and its contents.
    :param sequences: Path to the input file or raw input sequences
    :param default_file_format: Default file format to be used when loading a file
        Refer to the detect_sequences_format() for more details
    :param output_path: Optional path to the output file with collected results
    :return: Number of processed genomes
    """
    records_count = 0
    input_records = load_sequences(
        sequences, default_file_format=default_file_format
    )
    if LOGS.fasta_merge.is_progress_enabled():
        input_records = LOGS.fasta_merge.progress(
            input_records,
            message="Evaluating perform_on_genomes(...) on genomes",
            total=count_sequences(
                sequences, default_file_format=default_file_format
            ),
        )
    out_records: List[SeqRecord] = []
    for record in input_records:
        records_count = records_count + 1
        out = action(record)
        if out is not None and output_path is not None:
            out_records.append(out)
    if output_path is not None:
        with open_url(output_path, "w") as f:
            for record in out_records:
                print(f'>{record.id.replace(" ", "").rstrip()}', file=f)
                print(str(record), file=f)
    return records_count
