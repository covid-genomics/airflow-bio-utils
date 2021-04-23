from typing import Any, Callable, Sequence, Union

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from Bio.SeqRecord import SeqRecord

from airflow_bio_utils.sequences.transformations import \
    perform_on_each_sequence

from .utils import resolve_callable


class SequenceForEachOperator(PythonOperator):
    action: Callable[[SeqRecord], None]
    sequences: Union[
        Union[str, Sequence[SeqRecord]],
        Callable[..., Union[str, Sequence[SeqRecord]]],
    ]
    default_file_format: Union[str, Callable[..., str]] = "fasta"
    callback: Callable[[int], Any]

    @apply_defaults
    def __init__(
        self,
        action: Callable[[SeqRecord], None],
        sequences: Union[
            Union[str, Sequence[SeqRecord]],
            Callable[..., Union[str, Sequence[SeqRecord]]],
        ],
        callback: Callable[[int], Any],
        default_file_format: Union[str, Callable[..., str]] = "fasta",
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.action = action
        self.sequences = sequences
        self.default_file_format = default_file_format
        self.callback = callback

    def _execute_operator(self, *args, **kwargs):
        self.callback(
            perform_on_each_sequence(
                self.action,
                resolve_callable(self.sequences, *args, **kwargs),
                resolve_callable(self.default_file_format, *args, **kwargs),
            ),
            *args,
            **kwargs,
        )
