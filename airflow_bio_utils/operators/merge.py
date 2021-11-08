import traceback
from typing import Callable, Sequence, Tuple, Union

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from airflow_bio_utils.logs import LOGS
from airflow_bio_utils.sequences.merge import merge_sequences

from .utils import resolve_callable


class SequenceMergeOperator(PythonOperator):
    paths: Union[Sequence[str], Callable[..., Sequence[str]]]
    output_path: Union[str, Callable[..., str]]
    deduplicate: bool
    large_input: bool

    @apply_defaults
    def __init__(
        self,
        paths: Union[Sequence[str], Callable[..., Sequence[str]]],
        output_path: Union[str, Callable[..., str]],
        deduplicate: bool = True,
        large_input: bool = False ** kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.paths = paths
        self.output_path = output_path
        self.deduplicate = deduplicate
        self.large_input = large_input

    def _execute_operator(self, *args, **kwargs) -> Tuple[int, int]:
        try:
            return merge_sequences(
                resolve_callable(self.paths, *args, **kwargs),
                resolve_callable(self.output_path, *args, **kwargs),
                deduplicate=self.deduplicate,
                large_input=self.large_input,
            )
        except Exception as e:
            LOGS.merge.error(traceback.format_exc())
            raise e
