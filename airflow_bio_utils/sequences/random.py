from random import choice, randrange
from typing import Callable, List, Optional, TextIO, Tuple, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from airflow_bio_utils.filesystem import open_url

# Weight of generated nucleotides for ACTG
NUCLEOTIDES_ACTG = [("C", 10), ("G", 20), ("A", 40), ("T", 30)]

# Weight of generated nucleotides for FASTA-compat mode
NUCLEOTIDES_FULL = [
    ("C", 140),
    ("G", 180),
    ("A", 260),
    ("T", 220),
    ("R", 3),
    ("Y", 3),
    ("K", 3),
    ("M", 3),
    ("S", 3),
    ("W", 3),
    ("B", 3),
    ("D", 3),
    ("H", 3),
    ("V", 3),
    ("N", 2),
    ("-", 1),
]

SeqMetaGenerator = Optional[
    Callable[[Seq, int], Tuple[Optional[str], Optional[str], Optional[str]]]
]


def create_random_sequences(
    count: int,
    min_length: int = 10,
    max_length: int = 100,
    translate: bool = False,
    nucleotides: Optional[List[Tuple[str, int]]] = None,
) -> Union[List[Seq], str]:
    """
    Create list of random Seq objects.

    You can optionally pass list of tuples with nucleotides and their non-zero integer weights that
    determine probability of returning that nucleotide.

    :param count: Number of sequences to return
    :param min_length: Minimal sequence length
    :param max_length: Maximal sequence length
    :param translate: Translate nucleotides into codons
    :param nucleotides: Possible nucleotides list
    :return: List of generated sequences
    """
    if nucleotides is None:
        nucleotides = NUCLEOTIDES_ACTG
    choices = "".join(x * y for x, y in nucleotides)
    for seq_index in range(count):
        length = randrange(min_length, max_length)
        sequence = ""
        for nucleotide_index in range(length):
            sequence += choice(choices)
        single_sequence = Seq(sequence)
        if translate:
            single_sequence = single_sequence[
                : len(single_sequence) - len(single_sequence) % 3
            ]
            try:
                single_sequence = single_sequence.translate()
            except Exception:
                continue
        yield single_sequence


def default_str(value: Optional[str], default_value: Optional[str]) -> str:
    """
    If value is None returns default_value.
    But if default_value is also None, returns empty string.

    :param value: Optional value
    :param default_value: Optional fallback value
    :return: value if it's not None or default_value if it's not None or empty string
    """
    if value is None:
        if default_value is None:
            return ""
        return default_value
    return value


def default_seq_meta_generator(
    sequence: Seq, index: int
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Default function to generate random sequences metadata.

    :param sequence: Sequence string
    :param index: Sequence index
    :return: Tuple with sequence metadata
    """
    return (
        f"seq{index}",
        f"RandomSequence{index}",
        "Randomly generated genetic sequence record",
    )


def create_random_sequnce_records(
    count: int,
    output: Optional[TextIO] = None,
    output_mode: str = "fasta",
    min_length: int = 10,
    max_length: int = 100,
    translate: bool = False,
    nucleotides: Optional[List[Tuple[str, int]]] = None,
    meta_generator: SeqMetaGenerator = None,
) -> List[SeqRecord]:
    """
    Create list of random nucelotide sequence records.

    If you specify output parameter the sequences are written to the file handle.
    For example:

    with open_url('output.fasta', 'w') as out:
        create_random_sequnce_records(output=out, output_mode='fasta', ...)

    The meta_generator is by default default_seq_meta_generator.
    It's the function that accepts sequence string and its index and return tuple
    [id, name, description] for the single record.
    Please refer to the default_seq_meta_generator implementation for more details.

    :param count: Number of sequences to return
    :param output: Optional file handle to output generated sequences
    :param output_mode: Optional Biopython output mode (if output is specified)
    :param min_length: Minimal length of sequences
    :param max_length: Maximal length of sequences
    :param translate: Translate nucleotides into codons
    :param nucleotides: Possible nucleotides list (see create_random_sequences for more details)
    :param meta_generator: Optional function that returns sequence metadata
    :return: List of generated sequence records
    """
    if meta_generator is None:
        meta_generator = default_seq_meta_generator
    sequences = list(
        create_random_sequences(
            count,
            min_length=min_length,
            max_length=max_length,
            translate=translate,
            nucleotides=nucleotides,
        )
    )
    records = []
    for (id, seq) in enumerate(sequences):
        seq_id, name, description = meta_generator(seq, id)
        (
            default_seq_id,
            default_name,
            default_description,
        ) = default_seq_meta_generator(seq, id)
        seq_id = default_str(seq_id, default_seq_id)
        name = default_str(name, default_name)
        description = default_str(description, default_description)

        records.append(
            SeqRecord(
                seq,
                id=seq_id,
                name=name,
                description=description,
            )
        )
    if output is not None:
        SeqIO.write(records, output, output_mode)
    return records
