from typing import Any, Callable, Sequence, Union, Optional

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from Bio.SeqRecord import SeqRecord
from airflow_bio_utils.logs import LOGS

from airflow_bio_utils.sequences.transformations import \
    perform_on_each_sequence

from .utils import resolve_callable


class SequenceForEachOperator(PythonOperator):
    action: Callable[[SeqRecord], Optional[SeqRecord]]
    sequences: Union[
        Union[str, Sequence[SeqRecord]],
        Callable[..., Union[str, Sequence[SeqRecord]]],
    ]
    output_path: Union[Optional[str], Callable[..., Optional[str]]]
    default_file_format: Union[str, Callable[..., str]] = "fasta"
    callback: Callable[[int], Any]

    @apply_defaults
    def __init__(
        self,
        action: Callable[[SeqRecord], Optional[SeqRecord]],
        sequences: Union[
            Union[str, Sequence[SeqRecord]],
            Callable[..., Union[str, Sequence[SeqRecord]]],
        ],
        callback: Callable[[int], Any],
        default_file_format: Union[str, Callable[..., str]] = "fasta",
        output_path: Union[Optional[str], Callable[..., Optional[str]]] = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.action = action
        self.sequences = sequences
        self.default_file_format = default_file_format
        self.callback = callback
        self.output_path = output_path

    def _execute_operator(self, *args, **kwargs):
        try:
            output = perform_on_each_sequence(
                self.action,
                resolve_callable(self.sequences, *args, **kwargs),
                resolve_callable(self.default_file_format, *args, **kwargs),
                output_path=resolve_callable(self.output_path, *args, **kwargs),
            )

            if self.callback:
                return self.callback(
                    output,
                    *args,
                    **kwargs,
                )
            else:
                return output
        except Exception as e:
            LOGS.merge.error(str(e))
            raise e

