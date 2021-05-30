import click
from airflow_bio_utils.filesystem import open_url

from airflow_bio_utils.sequences.random import (NUCLEOTIDES_ACTG,
                                                NUCLEOTIDES_FULL,
                                                create_random_sequnce_records)


@click.command()
@click.argument("count")
@click.argument("output_filename")
@click.option("--output-mode", default="fasta", help="Biopython output mode")
@click.option(
    "--translate", is_flag=True, help="Translate nucleotides into codons"
)
@click.option(
    "--nucleotides",
    default="actg",
    type=click.Choice(["actg", "full"], case_sensitive=False),
    help="Specify nucleotides generation mode.",
)
@click.option("--min-length", default=10, help="Minimum sequence length")
@click.option("--max-length", default=100, help="Maximum sequence length")
def create_random_sequences_cli(
    count,
    output_filename,
    output_mode,
    translate,
    nucleotides,
    min_length,
    max_length,
):
    """
    Create list of random Seq objects.

    You can optionally pass list of tuples with nucleotides and their non-zero integer weights that
    determine probability of returning that nucleotide.
    """
    nucl = NUCLEOTIDES_ACTG
    if nucleotides == "full":
        nucl = NUCLEOTIDES_FULL

    with open_url(output_filename, "w") as output_file:
        count = len(
            create_random_sequnce_records(
                int(count),
                min_length=int(min_length),
                max_length=int(max_length),
                output=output_file,
                output_mode=output_mode,
                translate=translate,
                nucleotides=nucl,
            )
        )
        print(f"Saved {count} records to file {output_filename}.")


if __name__ == "__main__":
    create_random_sequences_cli()
