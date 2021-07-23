import click

from airflow_bio_utils.sequences.filter import (
    DEFAULT_ACCEPTED_SEQUENCE_SYMBOLS, filter_sequences, FilterSymbols, FilterByLength)


@click.command()
@click.argument("paths", nargs=-1)
@click.option("--min-seq-len", default=0, type=int)
@click.option("--max-seq-len", default=50000, type=int)
@click.option("--accepted-symbols", default=DEFAULT_ACCEPTED_SEQUENCE_SYMBOLS)
def filter_sequences_cli(paths, min_seq_len, max_seq_len, accepted_symbols):
    """
    Filter file with sequences.
    Function takes list of input files and for each performs filtering of the sequences.
    A new file is generated.
    If the output_paths is not None, then it should have the same length as input_paths.
    In that case from a file input_paths[i] a file output_paths[i] is generated.
    If output_paths is None (by default), then new file will have the same name, but with ".filtered" appended.
    """
    return filter_sequences(
        input_paths=paths,
        filters=[
            FilterSymbols(accepted_symbols),
            FilterByLength(
                min_length=min_seq_len,
                max_length=max_seq_len,
            ),
        ],
    )


if __name__ == "__main__":
    filter_sequences_cli()
