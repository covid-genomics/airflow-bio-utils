import traceback
from typing import Callable, List, Optional, Sequence, Union

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from airflow_bio_utils.logs import LOGS
from airflow_bio_utils.sequences.filter import (FilterCondition,
                                                FilterResultsMetadata,
                                                filter_sequences)

from .utils import resolve_callable


class SequenceFilterOperator(PythonOperator):
    input_paths: Union[str, Sequence[str], Callable[..., Sequence[str]]]
    filters: Union[
        Sequence[FilterCondition], Callable[..., Sequence[FilterCondition]]
    ]
    output_paths: Union[
        Optional[Sequence[str]], Callable[..., Optional[Sequence[str]]]
    ]

    @apply_defaults
    def __init__(
        self,
        input_paths: Union[str, Sequence[str], Callable[..., Sequence[str]]],
        filters: Union[
            Sequence[FilterCondition], Callable[..., Sequence[FilterCondition]]
        ],
        output_paths: Union[
            Optional[Sequence[str]], Callable[..., Optional[Sequence[str]]]
        ] = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.input_paths = input_paths
        self.output_paths = output_paths
        self.filters = filters

    def _execute_operator(
        self, *args, **kwargs
    ) -> List[FilterResultsMetadata]:
        try:
            input_paths = resolve_callable(self.input_paths, *args, **kwargs)
            if isinstance(input_paths, str):
                input_paths = [input_paths]
            output_paths = resolve_callable(self.output_paths, *args, **kwargs)
            if isinstance(output_paths, str):
                output_paths = [output_paths]
            filters = resolve_callable(self.filters, *args, **kwargs)
            return filter_sequences(
                input_paths,
                filters=filters,
                output_paths=output_paths,
            )
        except Exception as e:
            LOGS.merge.error(traceback.format_exc())
            raise e
