#!/usr/bin/python
# pylint: disable=bad-builtin
import collections
import datetime
import glob
import itertools
import math
import os
import pickle
import shlex
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from operator import mul
from subprocess import PIPE, Popen
from typing import Any, Dict, Generator, KeysView, List, Optional, Set, Tuple

import dateutil.parser as parser
from fs import open_fs
from scipy import stats
from simplesam import Reader, Sam
from six import BytesIO
from six.moves import reduce

from airflow_bio_utils.fasta_metrics.plot_exports import export_figure
from airflow_bio_utils.filesystem import open_url, url_join_path
from airflow_bio_utils.logs import LOGS

from .calc.adapter_kmers_plot import adapter_kmers_plot
from .calc.adapters import all_adapter_sequences
from .calc.cycle_spec_composition_plot import cycle_spec_composition_plot
from .calc.dates_distribution_plot import dates_distribution_plot
from .calc.gc_distribution_plot import gc_distribution_plot
from .calc.gc_plot import gc_plot
from .calc.global_properties_plot import global_properties_plot
from .calc.methylation_bias_plot import methylation_bias_plot
from .calc.mismatch_plot import mismatch_plot
from .calc.quality_distribution_plot import quality_distribution_plot
from .calc.quality_map_plot import quality_map_plot
from .calc.quality_plot import quality_plot
from .calc.sequence_length_distribution_plot import \
    sequence_length_distribution_plot
from .calc.sequence_specific_kmers_plot import sequence_specific_kmers_plot
from .calc.utils import (CollectMetricsArgs, FastaReader, Fastq, FastqReader,
                         MappedRead, bam_read_count, gc, mean, padbases,
                         percentile, window)
from .utils.fastq_generate import generate_fastq_file

DEFAULT_BASE_PROBS: List[float] = [0.25, 0.25, 0.25, 0.25, 0.1]

STATE_PICKLE_VARS: Dict[str, str] = dict(
    sequence_lengths="read_len",
    seqs_count="seqs_count",
    seqs_with_n_count="seqs_with_n",
    seqs_with_non_actgn_symbols_count="seqs_with_any_non_actgn_symbol",
    non_actgn_symbols="non_actgn_symbols",
    sequence_dates="dates",
    duplicate_sequences="duplicates",
    gc_sequence_content="cycle_gc",
    gc_content_per_position="pos_gc",
    quality_per_position="cycle_qual",
    has_quality_data="has_quality_data",
    sequence_positions="positions",
    quality_quantiles="quantiles",
    sequence_bases="bases",
    sequence_base_calls="basecalls",
    sequence_nuc="cycle_nuc",
    sequence_kmers="cycle_kmers",
    sequence_top_kmers="top_kmers",
    quality_quantile_values="quantiles",
    depths="depths",
    adapter_kmers="adapter_kmers",
    median_qual="median_qual",
    bad_kmers="bad_kmers",
)


@dataclass(frozen=True)
class PartiallyProcessedRead:
    cycle_nuc: Dict[int, Dict[str, int]]
    cycle_qual: Dict[int, Dict[int, int]]
    cycle_gc: Dict[int, int]
    cycle_kmers: Dict[int, Dict[str, int]]
    cycle_mismatch: Dict[str, Dict[str, Dict[str, int]]]
    skips: int
    read_len: Dict[int, int]
    seqs_with_n: int
    seqs_with_any_non_actgn_symbol: int
    non_actgn_symbols: Set[str]
    seqs: List[str]
    names: List[str]


def defaultdict_creator():
    return defaultdict(int)


def normalize_value(val):
    if isinstance(val, collections.abc.KeysView) or isinstance(
        val, collections.abc.ItemsView
    ):
        val = list(val)
    elif isinstance(val, type({}.keys())) or isinstance(
        val, type({}.values())
    ):
        val = list(val)
    if isinstance(val, list):
        return [normalize_value(item) for item in val]
    return val


def process_reads(reads: List[Optional[MappedRead]]) -> PartiallyProcessedRead:
    cycle_nuc = defaultdict(defaultdict_creator)
    cycle_qual = defaultdict(defaultdict_creator)
    cycle_kmers = defaultdict(defaultdict_creator)
    cycle_mismatch = {
        "C": defaultdict(defaultdict_creator),
        "G": defaultdict(defaultdict_creator),
        "A": defaultdict(defaultdict_creator),
        "T": defaultdict(defaultdict_creator),
    }
    skips = 0
    read_len = defaultdict(int)
    cycle_gc = defaultdict(int)
    seqs_with_n = 0
    seqs_with_any_non_actgn_symbol = 0
    non_actgn_symbols = set()
    reads = list(reads)
    seqs = []
    names = []

    for read in reads:
        if read is None:
            continue
        seq = read.seq
        name = read.name
        seqs.append(seq)
        names.append(name)
        qual = read.qual
        ref = read.ref
        leftlimit = read.l
        rightlimit = read.r
        kmer = read.kmer

        seq_letters = set(seq)
        if "N" in seq_letters:
            seqs_with_n += 1

        non_actgn = list(seq_letters.difference(set("ACTGN")))
        if len(non_actgn) > 0:
            seqs_with_any_non_actgn_symbol += 1
            for symbol in non_actgn:
                non_actgn_symbols.add(symbol)

        # Set up limits
        if (leftlimit == 1) and (rightlimit < 0):
            pass
        elif (leftlimit >= 1) and (rightlimit > 0):
            try:
                seq = seq[leftlimit - 1 : rightlimit]
                qual = qual[leftlimit - 1 : rightlimit]
            except IndexError:
                skips += 1
                continue
        elif (leftlimit > 1) and (rightlimit < 0):
            try:
                seq = seq[leftlimit - 1 :]
                qual = qual[leftlimit - 1 :]
            except IndexError:
                skips += 1
                continue
        if len(seq) == 0:
            skips += 1
            continue

        read_len[len(qual)] += 1
        cycle_gc[gc(seq)] += 1

        for i, (s, q) in enumerate(zip(seq, qual)):
            cycle_nuc[leftlimit + i][s] += 1
            cycle_qual[leftlimit + i][q] += 1

        if kmer is not None:
            for i, kmer in enumerate(window(seq, n=kmer)):
                cycle_kmers[leftlimit + i][kmer] += 1

        if ref is not None:
            try:
                for i, (s, r) in enumerate(zip(seq, ref)):
                    if s != r:
                        try:
                            cycle_mismatch[r][leftlimit + i][s] += 1
                        except KeyError:
                            pass
            except KeyError:
                pass

    return PartiallyProcessedRead(
        skips=skips,
        seqs=seqs,
        names=names,
        cycle_gc=cycle_gc,
        read_len=read_len,
        cycle_nuc=cycle_nuc,
        cycle_qual=cycle_qual,
        cycle_kmers=cycle_kmers,
        cycle_mismatch=cycle_mismatch,
        seqs_with_n=seqs_with_n,
        seqs_with_any_non_actgn_symbol=seqs_with_any_non_actgn_symbol,
        non_actgn_symbols=non_actgn_symbols,
    )


class CollectMetricsTask:

    args: CollectMetricsArgs
    est_counter: int
    sample_lengths: List[int]
    sample_binsizes: List[int]
    act_nlines: int
    mean_len: float
    mean_bentry: float
    seqs_with_n: int
    seqs_with_any_non_actgn_symbol: int
    non_actgn_symbols: Set[str]
    dates: Dict[datetime.datetime, int]
    has_quality_data: bool
    seqs_count: int
    read_len: Dict[int, int]
    cycle_nuc: Dict[int, Dict[str, int]]
    cycle_qual: Dict[int, Dict[int, int]]
    cycle_gc: Dict[int, int]
    cycle_kmers: Dict[int, Dict[str, int]]
    cycle_mismatch: Dict[str, Dict[str, Dict[str, int]]]
    bad_kmers: List[Tuple[str, float, float]]
    pos_gc: List[float]
    top_kmers: List[str]
    median_qual: Optional[float]
    positions: List[int]
    quantiles: List[List[int]]
    depths: List[int]
    basecalls: List[KeysView[str]]
    bases: Set[str]
    adapter_kmers: Set[str]
    quantile_values: List[float]
    sample_name: str
    duplicates: int
    percent_complete: int
    reads: int
    calculated_values: Dict[str, Any]
    reads: Generator[Fastq, None, None]

    def __init__(self, args: CollectMetricsArgs):
        self.args = args

    def run(self):
        self.section_prepare()
        self.section_load()
        self.section_generate()
        self.section_summary()
        self.upload()

    def upload(self):
        for path in self.output_fs.listdir("."):
            with open_url(url_join_path(self.args.output, path), "wb") as out:
                out.write(self.output_fs.readbytes(path))
        self.output_fs.close()

    def section_prepare(self):

        self.base_probs = self.args.base_probs
        if self.base_probs is None:
            self.base_probs = DEFAULT_BASE_PROBS

        output = sys.stdout
        if self.args.debug_output != "-":
            output = open(output, "w")

        self.bsize = None
        self.est_nlines = 0

        self.input = open_url(self.args.input).local_file_path
        if self.input == "-":
            self.input = "<stdin>"

        self.time_start = time.time()
        if self.input != "<stdin>":
            self.bsize = os.path.getsize(self.input)

        self.est_counter = int()
        self.sample_lengths = list()
        self.sample_binsizes = list()
        self.act_nlines = int()
        name, ext = os.path.splitext(self.input)
        if (self.args.leftlimit > 0) and (self.args.rightlimit > 0):
            if self.args.rightlimit < self.args.leftlimit:
                if self.args.debug_output != "-":
                    output.close()
                LOGS.fasta_metrics.fatal(
                    "Left limit must be less than right limit."
                )
        if self.args.type:
            ext = "." + self.args.type
        if (
            ext not in [".fq", ".fastq", ".sam", ".bam", ".gz", ".fasta"]
            and self.input != "<stdin>"
        ):
            if self.args.debug_output != "-":
                output.close()
            LOGS.fasta_metrics.fatal(
                f'Input file must end in either .sam, .bam, .fastq, .fasta, or .fastq.gz: Got "{self.input}" with extension {ext}'
            )

        if self.args.name:
            self.sample_name = self.args.name
        else:
            self.sample_name = self.input

        self.has_quality_data = True

        # estimate the number of lines in self.input if we can
        if ext in [".fastq", ".fq", ".fasta"]:
            with (
                FastqReader(open(self.input))
                if ext is not ".fasta"
                else FastaReader(open(self.input))
            ) as fh:
                for read in fh:
                    self.sample_lengths.append(len(read))
                    self.sample_binsizes.append(len(str(read)))
                    self.est_counter += 1
                    if self.est_counter == 10000:
                        break
                self.mean_bentry = mean(self.sample_binsizes)
                self.mean_len = mean(self.sample_lengths)
                self.est_nlines = int(self.bsize / self.mean_bentry)
                if not self.args.quiet:
                    LOGS.fasta_metrics.info(
                        "At {bytes:.0f} bytes per read of {len:.0f} length "
                        "we estimate {est:,} reads in input file.".format(
                            bytes=self.mean_bentry,
                            len=self.mean_len,
                            est=self.est_nlines,
                        )
                    )
        elif ext == ".sam":
            with Reader(open(self.input)) as fh:
                for read in fh:
                    self.sample_lengths.append(len(read))
                    self.sample_binsizes.append(len(str(read)))
                    self.est_counter += 1
                    if self.est_counter == 10000:
                        break
                self.mean_bentry = mean(self.sample_binsizes)
                self.mean_len = mean(self.sample_lengths)
                self.est_nlines = int(self.bsize / self.mean_bentry)
                if not self.args.quiet:
                    LOGS.fasta_metrics.info(
                        "At {bytes:.0f} bytes per read of {len:.0f} length "
                        "we estimate {est:,} reads in input file.".format(
                            bytes=self.mean_bentry,
                            len=self.mean_len,
                            est=self.est_nlines,
                        )
                    )
        elif ext == ".bam":
            est_nlines = sum(bam_read_count(self.input))
            if not self.args.quiet:
                LOGS.fasta_metrics.info(
                    "{est:,} reads in input file.".format(est=est_nlines)
                )
        elif ext == ".gz":
            if self.args.binsize:
                self.n = self.args.binsize
                self.est_nlines = None
                if not self.args.quiet:
                    LOGS.fasta_metrics.info(
                        "Reading from gzipped file, bin size (-s) set to {binsize:n}.".format(
                            binsize=self.n
                        )
                    )
            else:
                LOGS.fasta_metrics.info(
                    "Gzipped file detected. Reading file to determine bin size (-s)."
                )
                p1 = Popen(
                    shlex.split("gzip -dc %s" % self.input), stdout=PIPE
                )
                p2 = Popen(shlex.split("wc -l"), stdin=p1.stdout, stdout=PIPE)
                self.est_nlines, _ = p2.communicate()
                self.est_nlines = int(self.est_nlines) // 4
                if not self.args.quiet:
                    sys.stderr.write(
                        "{est:,} reads in input file.".format(
                            est=self.est_nlines
                        )
                    )
        elif name == "<stdin>":
            if self.args.binsize:
                self.n = self.args.binsize
            else:
                self.n = 1
            if not self.args.quiet:
                LOGS.fasta_metrics.info(
                    "Reading from <stdin>, bin size (-s) set to {binsize:n}.".format(
                        binsize=self.n
                    )
                )
            self.est_nlines = None
        if self.est_nlines == 0:
            if self.args.debug_output != "-":
                self.output.close()
            LOGS.fasta_metrics.fatal(
                "The input file appears empty. Please check the file for data."
            )

        if self.est_nlines is not None:
            # set up factor for sampling bin size
            if self.args.binsize:
                self.n = self.args.binsize
            else:
                nf = math.floor(self.est_nlines / self.args.nreads)
                if nf >= 1:
                    self.n = int(nf)
                else:
                    self.n = 1
            if not self.args.quiet:
                LOGS.fasta_metrics.info(
                    "Bin size (-s) set to {binsize:n}.".format(binsize=self.n)
                )

        if ext in [".sam", ".bam"]:
            infile = Reader(open(self.input))
        elif ext == ".fasta":
            self.has_quality_data = False
            infile = FastaReader(open(self.input), ext=ext)
        else:
            infile = FastqReader(open(self.input), ext=ext)

        self.seqs_with_n = 0
        self.seqs_with_any_non_actgn_symbol = 0
        self.non_actgn_symbols = set()

        self.dates = defaultdict(int)
        self.seqs_count = 0
        self.read_len = defaultdict(int)
        self.cycle_nuc = defaultdict(partial(defaultdict, int))
        self.cycle_qual = defaultdict(partial(defaultdict, int))
        self.cycle_gc = defaultdict(int)
        self.cycle_kmers = defaultdict(partial(defaultdict, int))
        self.cycle_mismatch = {
            "C": defaultdict(partial(defaultdict, int)),
            "G": defaultdict(partial(defaultdict, int)),
            "A": defaultdict(partial(defaultdict, int)),
            "T": defaultdict(partial(defaultdict, int)),
        }

        if self.args.count_duplicates:
            try:
                from pybloom import ScalableBloomFilter

                self.bloom_filter = ScalableBloomFilter(
                    mode=ScalableBloomFilter.SMALL_SET_GROWTH
                )
            except ImportError:
                if self.args.debug_output != "-":
                    output.close()
                LOGS.fasta_metrics.fatal(
                    "--count-duplicates option requires 'pybloom' package."
                )

        self.duplicates = 0
        self.percent_complete = 10
        self.reads = infile.subsample(self.n)

    def section_load(self):

        tmp_output = "fastqp_figures"
        self.fs_output_url = f"zip://{tmp_output}.zip"
        if self.args.output_to is None:
            pass
        elif self.args.output_to == "zip":
            self.fs_output_url = f"zip://{tmp_output}.zip"
        elif self.args.output_to == "directory":
            self.fs_output_url = tmp_output
        elif self.args.output_to == "memory":
            self.fs_output_url = "mem://"
        else:
            raise Exception(
                f'Invalid output mode was specified in parameter output_to="{self.args.output_to}"'
            )

        self.output_fs = open_fs(self.fs_output_url, create=True)

        if self.args.use_cache and False:
            if hasattr(self, "calculated_values"):
                for key in STATE_PICKLE_VARS:
                    setattr(
                        self,
                        STATE_PICKLE_VARS[key],
                        self.calculated_values[key],
                    )
                return None
            elif self.args.calculated_data_cache_file is not None:
                if self.output_fs.exists(self.args.calculated_data_cache_file):
                    print(
                        f"Load existing data from {self.args.calculated_data_cache_file}"
                    )
                    with open(
                        self.args.calculated_data_cache_file, "rb"
                    ) as file:
                        self.calculated_values = pickle.load(file)
                    for key in STATE_PICKLE_VARS:
                        setattr(
                            self,
                            STATE_PICKLE_VARS[key],
                            self.calculated_values[key],
                        )
                    return None

        def map_reads(read: Fastq) -> Optional[MappedRead]:
            ref: Optional[str] = None
            seq: Optional[str] = None
            qual: Optional[str] = None
            name: str = ""
            if isinstance(read, Sam):
                if self.args.aligned_only and not read.mapped:
                    return None
                elif self.args.unaligned_only and read.mapped:
                    return None
                if read.reverse:
                    seq = read.seq[::-1]
                    qual = read.qual[::-1]
                else:
                    seq = read.seq
                    qual = read.qual
                name = read.qname
            else:
                seq = read.seq
                qual = read.qual
                name = read.name

            if isinstance(read, Sam) and read.mapped:
                ref = read.parse_md()
            return MappedRead(
                ref=ref,
                seq=seq,
                qual=qual,
                name=name,
                kmer=self.args.kmer if self.args.examine_kmers else None,
                l=self.args.leftlimit,
                r=self.args.rightlimit,
            )

        use_multiprocessing = False
        if (
            self.args.use_multiprocessing
            and self.est_nlines > self.args.read_multiprocessing_threshold
        ):
            print(
                "Estimated number of reads exceeds threshold so the multiprocessing will be used"
            )
            use_multiprocessing = True

        for read in [process_reads([map_reads(read) for read in self.reads])]:
            if read is None:
                continue

            self.act_nlines += self.n * read.skips
            for k1 in read.cycle_gc:
                self.cycle_gc[k1] += read.cycle_gc[k1]

            self.seqs_count += len(read.seqs)
            for seq in read.seqs:
                if self.args.count_duplicates:
                    if seq in self.bloom_filter:
                        self.duplicates += 1
                    else:
                        self.bloom_filter.add(seq)

            for k1 in read.cycle_nuc:
                for k2 in read.cycle_nuc[k1]:
                    self.cycle_nuc[k1][k2] += read.cycle_nuc[k1][k2]
            for k1 in read.cycle_qual:
                for k2 in read.cycle_qual[k1]:
                    self.cycle_qual[k1][k2] += read.cycle_qual[k1][k2]

            for k1 in read.read_len:
                self.read_len[k1] += read.read_len[k1]

            self.seqs_with_n += read.seqs_with_n
            self.seqs_with_any_non_actgn_symbol += (
                read.seqs_with_any_non_actgn_symbol
            )

            for symbol in read.non_actgn_symbols:
                self.non_actgn_symbols.add(symbol)

            names = read.names
            for name in names:
                tokens = name.split("|")
                date = None
                for token in tokens:
                    try:
                        date = parser.parse(token)
                        break
                    except parser._parser.ParserError:
                        pass
                if date is not None:
                    self.dates[date] += 1

            if self.args.examine_kmers:
                for k1 in read.cycle_kmers:
                    for k2 in read.cycle_kmers[k1]:
                        self.cycle_kmers[k1][k2] += read.cycle_kmers[k1][k2]

            try:
                for k1 in read.cycle_mismatch:
                    for k2 in read.cycle_mismatch[k1]:
                        for k3 in read.cycle_mismatch[k1][k2]:
                            self.cycle_mismatch[k1][k2][
                                k3
                            ] += read.cycle_mismatch[k1][k2][k3]
            except KeyError:
                pass

            if self.est_nlines is not None:
                if (
                    self.act_nlines / self.est_nlines
                ) * 100 >= self.percent_complete:
                    # task.progress(self.act_nlines, total=self.est_nlines)
                    self.percent_complete += 10
            self.act_nlines += self.n

        self.positions = [k for k in sorted(self.cycle_qual.keys())]
        self.depths = [self.read_len[k] for k in sorted(self.read_len.keys())]

        self.basecalls = [
            self.cycle_nuc[k].keys() for k in sorted(self.cycle_nuc.keys())
        ]
        self.bases = set(list(itertools.chain.from_iterable(self.basecalls)))

        map(padbases(self.bases), self.cycle_nuc.values())

        self.quantile_values = [0.05, 0.25, 0.5, 0.75, 0.95]
        self.quantiles = []

        # replace ASCII quality with integer
        for _, v in sorted(self.cycle_qual.items()):
            for q in tuple(
                v.keys()
            ):  # py3 keys are iterator, so build a tuple to avoid recursion
                v[ord(str(q)) - 33] = v.pop(q)
            line = [percentile(v, p) for p in self.quantile_values]
            self.quantiles.append(line)

        # build kmer set of known adapter sequences
        self.adapter_kmers = set()
        if self.args.examine_kmers:
            for adapter in all_adapter_sequences:
                for kmer in window(adapter, n=self.args.kmer):
                    self.adapter_kmers.add(kmer)

        # test for nonuniform kmer profiles and calculate obs/exp
        observed_expected = dict()
        all_kmers = [
            self.cycle_kmers[k].keys() for k in sorted(self.cycle_kmers.keys())
        ]
        kmers = set(list(itertools.chain.from_iterable(all_kmers)))
        self.bad_kmers = []
        sequenced_bases = sum((l * n for l, n in self.read_len.items()))
        priors = tuple(self.base_probs)
        if self.args.examine_kmers:
            for kmer in kmers:
                kmer_counts = [
                    (i, self.cycle_kmers[i][kmer])
                    for i in sorted(self.cycle_kmers.keys())
                ]
                expected_fraction = reduce(
                    mul,
                    (
                        p ** kmer.count(b)
                        for b, p in zip(("A", "T", "C", "G", "N"), priors)
                    ),
                    1,
                )
                expected = expected_fraction * sequenced_bases
                observed_expected[kmer] = (
                    sum((n for _, n in kmer_counts)) / expected
                )
                slope, _, _, p_value, _ = stats.linregress(*zip(*kmer_counts))
                if abs(slope) > 2 and p_value < 0.05:
                    self.bad_kmers.append((kmer, slope, p_value))
        self.bad_kmers = sorted(self.bad_kmers, key=lambda x: x[2])[:10]

        self.pos_gc = []
        for i in self.positions:
            try:
                pg = (
                    sum([self.cycle_nuc[i]["C"], self.cycle_nuc[i]["G"]])
                    / sum(
                        [
                            self.cycle_nuc[i]["C"],
                            self.cycle_nuc[i]["G"],
                            self.cycle_nuc[i]["A"],
                            self.cycle_nuc[i]["T"],
                        ]
                    )
                    * 100
                )
            except ZeroDivisionError:
                pg = 0  # https://github.com/mdshw5/fastqp/issues/26
            self.pos_gc.append(pg)

        self.top_kmers = [fields[0] for fields in self.bad_kmers]
        self.median_qual = None

    def section_generate(self):

        all_plots = [
            global_properties_plot(
                self.args,
                self.read_len,
                self.seqs_count,
                self.seqs_with_n,
                self.seqs_with_any_non_actgn_symbol,
                self.non_actgn_symbols,
                self.dates,
                self.duplicates,
                self.pos_gc,
                self.cycle_qual,
                self.has_quality_data,
            ),
            dates_distribution_plot(self.args, self.dates),
            quality_plot(self.args, self.positions, self.quantiles)
            if self.has_quality_data
            else None,
            quality_distribution_plot(self.args, self.cycle_qual.values())
            if self.has_quality_data
            else None,
            quality_map_plot(self.args, self.cycle_qual)
            if self.has_quality_data
            else None,
            sequence_length_distribution_plot(self.args, self.read_len),
            gc_plot(self.args, self.positions, self.pos_gc),
            gc_distribution_plot(self.args, self.cycle_gc),
            cycle_spec_composition_plot(
                self.args, self.positions, self.bases, self.cycle_nuc
            ),
            sequence_specific_kmers_plot(
                self.args, self.positions, self.cycle_kmers, self.top_kmers
            )
            if self.args.examine_kmers
            else None,
            adapter_kmers_plot(
                self.args, self.positions, self.cycle_kmers, self.adapter_kmers
            )
            if self.args.examine_kmers
            else None,
        ]

        if not self.args.only_return_data:
            for item in all_plots:
                if item:
                    if list(item)[0] == "quality_distribution_plot":
                        self.median_qual = item[-1]
                    export_figure(item, self.output_fs, self.args)

        self.calculated_values = dict()
        for key in STATE_PICKLE_VARS:
            self.calculated_values[key] = normalize_value(
                getattr(self, STATE_PICKLE_VARS[key])
            )

        if (
            self.args.use_cache
            and self.args.calculated_data_cache_file is not None
        ):
            print(
                f"Save calculated data to {self.args.calculated_data_cache_file}"
            )
            with self.output_fs.open(
                self.args.calculated_data_cache_file, mode="wb"
            ) as outfile:
                pickle.dump(
                    self.calculated_values,
                    outfile,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )

        if self.args.create_pdf:
            from PIL import Image

            image_list: List[Image] = []
            for item in all_plots:
                if item is not None:
                    binary_output, output_type = export_figure(
                        item,
                        self.output_fs,
                        self.args,
                        only_return_binary=True,
                    )
                    if output_type == "png":
                        img = Image.open(binary_output)
                        img.load()  # required for png.split()
                        background = Image.new(
                            "RGB", img.size, (255, 255, 255)
                        )
                        background.paste(
                            img, mask=img.split()[3]
                        )  # 3 is the alpha channel
                        image_list.append(background)
            tf = tempfile.NamedTemporaryFile(suffix=".pdf")
            image_list[0].save(
                tf.name,
                "PDF",
                resolution=100.0,
                save_all=True,
                append_images=image_list[1:],
            )
            tf.seek(0)
            bytes_buf = BytesIO(tf.read())
            tf.close()
            with self.output_fs.open(f"metrics.pdf", mode="wb") as outfile:
                outfile.write(bytes_buf.getbuffer())

    def section_summary(self):
        time_finish = time.time()
        elapsed = time_finish - self.time_start
        if not self.args.quiet:
            LOGS.fasta_metrics.info(
                "There were {counts:,} reads in the file. Analysis finished in {sec}.".format(
                    counts=self.act_nlines,
                    sec=time.strftime("%H:%M:%S", time.gmtime(elapsed)),
                )
            )
            if len(self.bad_kmers) > 0:
                for kmer in self.bad_kmers:
                    LOGS.fasta_metrics.info(
                        "KmerWarning: kmer %s has a non-uniform profile (slope = %s, p = %s)."
                        % (kmer)
                    )
            if (
                self.has_quality_data
                and self.median_qual < self.args.median_qual
            ):
                LOGS.fasta_metrics.info(
                    "QualityWarning: median base quality score is %s."
                    % self.median_qual
                )

        if self.args.debug_output != "-":
            self.output.close()
        return self.calculated_values


class GeneFile:
    def __init__(
        self,
        input_path: str,  # Input file (one of .sam, .bam, .fq, or .fastq(.gz) or stdin (-))
    ):
        self.input_path = input_path

    def generate_example(self):
        if self.input_path != "-":
            generate_fastq_file(self.input_path)

    def collect_metrics(self, output: str = "fastqp_figures", **kwargs):
        return CollectMetricsTask(
            CollectMetricsArgs(
                **{
                    **kwargs,
                    **dict(
                        input=self.input_path,
                        output=output,
                    ),
                },
            ),
        ).run()
