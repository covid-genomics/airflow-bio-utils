from typing import Callable, Sequence, Union

from airflow.operators.python_operator import PythonOperator
from airflow.utils.decorators import apply_defaults
from typing import Tuple, Optional, List
from airflow_bio_utils.logs import LOGS

from airflow_bio_utils.fasta_metrics import GeneFile
from .utils import resolve_callable


class SequenceMetricsOperator(PythonOperator):
    input_path: Union[str, Callable[..., str]]
    output_path: Union[str, Callable[..., str]]
    quiet: Union[bool, Callable[..., bool]]
    binsize: Union[int, Callable[..., int]]
    name: Union[str, Callable[..., str]]
    nreads: Union[int, Callable[..., int]]
    base_probs: Union[Optional[List[float]], Callable[..., Optional[List[float]]]]
    kmer: Union[int, Callable[..., int]]
    debug_output: Union[str, Callable[..., str]]
    type: Union[Optional[str], Callable[..., Optional[str]]]
    leftlimit: Union[int, Callable[..., int]]
    rightlimit: Union[int, Callable[..., int]]
    median_qual: Union[int, Callable[..., int]]
    aligned_only: Union[bool, Callable[..., bool]]
    unaligned_only: Union[bool, Callable[..., bool]]
    count_duplicates: Union[bool, Callable[..., bool]]
    output_plotly_json: Union[bool, Callable[..., bool]]
    output_plotly_charts: Union[bool, Callable[..., bool]]
    output_matplot_images: Union[bool, Callable[..., bool]]
    output_images: Union[bool, Callable[..., bool]]
    output_csv: Union[bool, Callable[..., bool]]
    output_to: Union[str, Callable[..., str]]
    examine_kmers: Union[str, Callable[..., str]]
    use_multiprocessing: Union[str, Callable[..., str]]
    read_multiprocessing_threshold: Union[int, Callable[..., int]]
    use_cache: Union[bool, Callable[..., bool]]
    calculated_data_cache_file: Union[str, Callable[..., str]]
    only_return_data: Union[bool, Callable[..., bool]]
    figure_settings: Union[Optional[dict], Callable[..., Optional[dict]]]
    create_pdf: Union[bool, Callable[..., bool]]

    @apply_defaults
    def __init__(
        self,
        input_path: Union[str, Callable[..., str]],
        output_path: Union[str, Callable[..., str]],
        quiet: Union[bool, Callable[..., bool]] = False,
        binsize: Union[int, Callable[..., int]] = 0,
        name: Union[str, Callable[..., str]] = None,
        nreads: Union[int, Callable[..., int]] = 2000000,
        base_probs: Union[Optional[List[float]], Callable[..., Optional[List[float]]]] = None,
        kmer: Union[int, Callable[..., int]] = 5,
        debug_output: Union[str, Callable[..., str]] = '-',
        type: Union[Optional[str], Callable[..., Optional[str]]] = None,
        leftlimit: Union[int, Callable[..., int]] = 1,
        rightlimit: Union[int, Callable[..., int]] = -1,
        median_qual: Union[int, Callable[..., int]] = 30,
        aligned_only: Union[bool, Callable[..., bool]] = False,
        unaligned_only: Union[bool, Callable[..., bool]] = False,
        count_duplicates: Union[bool, Callable[..., bool]] = True,
        output_plotly_json: Union[bool, Callable[..., bool]] = False,
        output_plotly_charts: Union[bool, Callable[..., bool]] = False,
        output_matplot_images: Union[bool, Callable[..., bool]] = False,
        output_images: Union[bool, Callable[..., bool]] = False,
        output_csv: Union[bool, Callable[..., bool]] = False,
        output_to: Union[str, Callable[..., str]] = 'memory',
        examine_kmers: Union[str, Callable[..., str]] = False,
        use_multiprocessing: Union[str, Callable[..., str]] = True,
        read_multiprocessing_threshold: Union[int, Callable[..., int]] = 10000,
        use_cache: Union[bool, Callable[..., bool]] = False,
        calculated_data_cache_file: Union[str, Callable[..., str]] = 'calculated_data.pickle',
        only_return_data: Union[bool, Callable[..., bool]] = False,
        figure_settings: Union[Optional[dict], Callable[..., Optional[dict]]] = None,
        create_pdf: Union[bool, Callable[..., bool]] = True,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs, python_callable=self._execute_operator)
        self.input_path = input_path
        self.output_path = output_path
        self.quiet = quiet
        self.binsize = binsize
        self.name = name
        self.nreads = nreads
        self.base_probs = base_probs
        self.kmer = kmer
        self.debug_output = debug_output
        self.type = type
        self.leftlimit = leftlimit
        self.rightlimit = rightlimit
        self.median_qual = median_qual
        self.aligned_only = aligned_only
        self.unaligned_only = unaligned_only
        self.count_duplicates = count_duplicates
        self.output_plotly_json = output_plotly_json
        self.output_plotly_charts = output_plotly_charts
        self.output_matplot_images = output_matplot_images
        self.output_images = output_images
        self.output_csv = output_csv
        self.output_to = output_to
        self.examine_kmers = examine_kmers
        self.use_multiprocessing = use_multiprocessing
        self.read_multiprocessing_threshold = read_multiprocessing_threshold
        self.use_cache = use_cache
        self.calculated_data_cache_file = calculated_data_cache_file
        self.only_return_data = only_return_data
        self.figure_settings = figure_settings
        self.create_pdf = create_pdf
        
    
    def _execute_operator(self, *args, **kwargs):
        try:
            GeneFile(resolve_callable(self.input_path, *args, **kwargs)).collect_metrics(
                output=resolve_callable(self.output_path, *args, **kwargs),
                quiet=resolve_callable(self.quiet, *args, **kwargs),
                binsize=resolve_callable(self.binsize, *args, **kwargs),
                name=resolve_callable(self.name, *args, **kwargs),
                nreads=resolve_callable(self.nreads, *args, **kwargs),
                base_probs=resolve_callable(self.base_probs, *args, **kwargs),
                kmer=resolve_callable(self.kmer, *args, **kwargs),
                debug_output=resolve_callable(self.debug_output, *args, **kwargs),
                type=resolve_callable(self.type, *args, **kwargs),
                leftlimit=resolve_callable(self.leftlimit, *args, **kwargs),
                rightlimit=resolve_callable(self.rightlimit, *args, **kwargs),
                median_qual=resolve_callable(self.median_qual, *args, **kwargs),
                aligned_only=resolve_callable(self.aligned_only, *args, **kwargs),
                unaligned_only=resolve_callable(self.unaligned_only, *args, **kwargs),
                count_duplicates=resolve_callable(self.count_duplicates, *args, **kwargs),
                output_plotly_json=resolve_callable(self.output_plotly_json, *args, **kwargs),
                output_plotly_charts=resolve_callable(self.output_plotly_charts, *args, **kwargs),
                output_matplot_images=resolve_callable(self.output_matplot_images, *args, **kwargs),
                output_images=resolve_callable(self.output_images, *args, **kwargs),
                output_csv=resolve_callable(self.output_csv, *args, **kwargs),
                output_to=resolve_callable(self.output_to, *args, **kwargs),
                examine_kmers=resolve_callable(self.examine_kmers, *args, **kwargs),
                use_multiprocessing=resolve_callable(self.use_multiprocessing, *args, **kwargs),
                read_multiprocessing_threshold=resolve_callable(self.read_multiprocessing_threshold, *args, **kwargs),
                use_cache=resolve_callable(self.use_cache, *args, **kwargs),
                calculated_data_cache_file=resolve_callable(self.calculated_data_cache_file, *args, **kwargs),
                only_return_data=resolve_callable(self.only_return_data, *args, **kwargs),
                figure_settings=resolve_callable(self.figure_settings, *args, **kwargs),
                create_pdf=resolve_callable(self.create_pdf, *args, **kwargs),
            )
        except Exception as e:
            LOGS.merge.error(str(e))
            raise e

