from typing import Callable, List, Optional, Sequence, Union

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from typing import List

from airflow_bio_utils.sequences.filter import (DEFAULT_MAX_SEQ_LEN,
                                                DEFAULT_MIN_SEQ_LEN,
                                                filter_sequences)

from .utils import resolve_callable


class SequenceFilterOperator(PythonOperator):
    input_paths: Union[Sequence[str], Callable[..., Sequence[str]]]
    min_seq_len: Union[int, Callable[..., int]]
    max_seq_len: Union[int, Callable[..., int]]
    accepted_symbols: Union[
        Optional[List[str]], Callable[..., Optional[List[str]]]
    ]
    output_paths: Union[
        Optional[Sequence[str]], Callable[..., Optional[Sequence[str]]]
    ]

    @apply_defaults
    def __init__(
        self,
        input_paths: Union[Sequence[str], Callable[..., Sequence[str]]],
        min_seq_len: Union[int, Callable[..., int]] = DEFAULT_MIN_SEQ_LEN,
        max_seq_len: Union[int, Callable[..., int]] = DEFAULT_MAX_SEQ_LEN,
        accepted_symbols: Union[
            Optional[List[str]], Callable[..., Optional[List[str]]]
        ] = None,
        output_paths: Union[
            Optional[Sequence[str]], Callable[..., Optional[Sequence[str]]]
        ] = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.input_paths = input_paths
        self.min_seq_len = min_seq_len
        self.max_seq_len = max_seq_len
        self.accepted_symbols = accepted_symbols
        self.output_paths = output_paths

    def _execute_operator(self, *args, **kwargs) -> List[str]:
        return filter_sequences(
            resolve_callable(self.input_paths, *args, **kwargs),
            min_seq_len=resolve_callable(self.min_seq_len, *args, **kwargs),
            max_seq_len=resolve_callable(self.max_seq_len, *args, **kwargs),
            accepted_symbols=resolve_callable(
                self.accepted_symbols, *args, **kwargs
            ),
            output_paths=resolve_callable(self.output_paths, *args, **kwargs),
        )
