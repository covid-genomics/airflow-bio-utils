import traceback
from typing import Callable, List, Optional, Tuple, Union

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from airflow_bio_utils.filesystem import open_url
from airflow_bio_utils.logs import LOGS
from airflow_bio_utils.sequences.random import create_random_sequnce_records

from .utils import resolve_callable


class SequenceRandomOperator(PythonOperator):
    output: Union[Optional[str], Callable[..., Optional[str]]] = None
    count: Union[int, Callable[..., int]]
    min_length: Union[int, Callable[..., int]]
    max_length: Union[int, Callable[..., int]]
    translate: Union[bool, Callable[..., bool]]
    nucleotides: Union[
        Optional[List[Tuple[str, int]]],
        Callable[..., Optional[List[Tuple[str, int]]]],
    ]

    @apply_defaults
    def __init__(
        self,
        count: Union[int, Callable[..., int]],
        output: Union[Optional[str], Callable[..., Optional[str]]] = None,
        min_length: Union[int, Callable[..., int]] = 10,
        max_length: Union[int, Callable[..., int]] = 100,
        translate: Union[bool, Callable[..., bool]] = False,
        nucleotides: Union[
            Optional[List[Tuple[str, int]]],
            Callable[..., Optional[List[Tuple[str, int]]]],
        ] = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.count = count
        self.min_length = min_length
        self.max_length = max_length
        self.translate = translate
        self.nucleotides = nucleotides
        self.output = output

    def _execute_operator(self, *args, **kwargs):
        try:
            output = resolve_callable(self.output, *args, **kwargs)
            count = resolve_callable(self.count, *args, **kwargs)
            min_length = resolve_callable(self.min_length, *args, **kwargs)
            max_length = resolve_callable(self.max_length, *args, **kwargs)
            translate = resolve_callable(self.translate, *args, **kwargs)
            nucleotides = resolve_callable(self.nucleotides, *args, **kwargs)

            if output is not None:
                with open_url(output, "w") as out:
                    create_random_sequnce_records(
                        output=out,
                        count=count,
                        min_length=min_length,
                        max_length=max_length,
                        translate=translate,
                        nucleotides=nucleotides,
                    )
                    return output
            else:
                return create_random_sequnce_records(
                    output=output,
                    count=count,
                    min_length=min_length,
                    max_length=max_length,
                    translate=translate,
                    nucleotides=nucleotides,
                )
        except Exception as e:
            LOGS.merge.error(traceback.format_exc())
            raise e
