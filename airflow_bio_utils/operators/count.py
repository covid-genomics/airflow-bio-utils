import traceback
from typing import Any, Callable, Optional, Sequence, Union

from Bio.SeqRecord import SeqRecord

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from airflow_bio_utils.logs import LOGS
from airflow_bio_utils.sequences.transformations import count_sequences

from .utils import resolve_callable


class SequenceCountOperator(PythonOperator):
    sequences: Union[
        Union[str, Sequence[SeqRecord]],
        Callable[..., Union[str, Sequence[SeqRecord]]],
    ]
    default_file_format: Union[str, Callable[..., str]] = "fasta"
    callback: Optional[Callable[[int], Any]]

    @apply_defaults
    def __init__(
        self,
        sequences: Union[
            Union[str, Sequence[SeqRecord]],
            Callable[..., Union[str, Sequence[SeqRecord]]],
        ],
        callback: Optional[Callable[[int], Any]] = None,
        default_file_format: Union[str, Callable[..., str]] = "fasta",
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.sequences = sequences
        self.default_file_format = default_file_format
        self.callback = callback

    def _execute_operator(self, *args, **kwargs) -> int:
        try:
            seq_count = count_sequences(
                resolve_callable(self.sequences, *args, **kwargs),
                resolve_callable(self.default_file_format, *args, **kwargs),
            )
            if self.callback:
                return self.callback(
                    seq_count,
                    *args,
                    **kwargs,
                )
            else:
                return seq_count
        except Exception as e:
            LOGS.merge.error(traceback.format_exc())
            raise e
