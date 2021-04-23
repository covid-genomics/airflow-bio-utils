import click

from airflow_bio_utils.sequences.merge import merge_sequences


@click.command()
@click.argument("paths", nargs=-1)
@click.argument("output_path", nargs=1)
def merge_sequences_cli(*args, **kwargs):
    """
    Merge multiple files with sequences into one single FASTA file.
    All sequences are deduplicated using sequence IDs.
    """
    return merge_sequences(*args, **kwargs)


if __name__ == "__main__":
    merge_sequences_cli()
