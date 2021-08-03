from __future__ import annotations, division

import os
import re
import string
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import groupby, islice
from typing import Generator, List, Optional

import matplotlib as mpl
import plotly
import plotly.tools as tls
from six import PY2, PY3, BytesIO, string_types
from six.moves import range, zip

try:
    from collections import Counter
except ImportError:
    from . import Counter

from io import TextIOWrapper
from subprocess import PIPE, Popen

DEFAULT_FIGURE_SETTINGS = dict(
    figsize=(8, 6),
)


@dataclass(frozen=True)
class MappedRead:
    ref: Optional[str] = None
    seq: Optional[str] = None
    qual: Optional[str] = None
    name: str = ""
    kmer: Optional[int] = 5
    l: int = 1
    r: int = -1


@dataclass(frozen=True)
class CollectMetricsArgs:
    input: str = "-"
    quiet: bool = False
    binsize: int = 0
    name: str = None
    nreads: int = 2000000
    base_probs: Optional[List[float]] = None
    kmer: int = 5
    output: str = "fastqp_figures"
    debug_output: str = "-"
    type: Optional[str] = None
    leftlimit: int = 1
    rightlimit: int = -1
    median_qual: int = 30
    aligned_only: bool = False
    unaligned_only: bool = False
    count_duplicates: bool = True
    output_plotly_json: bool = True
    output_plotly_charts: bool = False
    output_matplot_images: bool = True
    output_images: bool = True
    output_csv: bool = True
    output_to: str = "directory"
    examine_kmers: str = True
    use_multiprocessing: str = True
    read_multiprocessing_threshold: int = 10000
    use_cache: bool = True
    calculated_data_cache_file: str = "calculated_data.pickle"
    only_return_data: bool = False
    figure_settings: Optional[dict] = None
    create_pdf: bool = True


if sys.platform is not "darwin":
    mpl.use("Agg")


def should_create_matplot(settings: CollectMetricsArgs) -> bool:
    return (
        settings.create_pdf
        or settings.output_matplot_images
        or settings.output_images
    )


def should_create_plotly(settings: CollectMetricsArgs) -> bool:
    return (
        settings.create_pdf
        or settings.output_plotly_charts
        or settings.output_plotly_json
        or settings.output_images
    )


class Gzip(object):
    """Call system gzip and maintain interface compatibility with python
    gzip module"""

    def __init__(self, filename, mode):
        self.stream, self.p = self.open(filename, mode)
        self.mode = mode
        self.filename = filename

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        return next(self.stream)

    def open(self, filename, mode):
        if "r" in mode:
            self.fh = open(filename, "rb", 0)
            p = Popen(["gzip", "-dc", filename], stdout=PIPE, stderr=PIPE)
            if "b" in mode:
                fh = p.stdout
            else:
                try:
                    fh = TextIOWrapper(p.stdout)
                except AttributeError:  # python2.7?
                    fh = p.stdout
        elif "w" in mode:
            self.fh = open(filename, "wb", 0)
            p = Popen(["gzip", "-c"], stdin=PIPE, stdout=self.fh)
            fh = p.stdin
        return (fh, p)

    def write(self, string):
        self.stream.write(string.encode("utf-8"))

    def read(self, string):
        self.stream.read(string)

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.p.communicate()
        if self.fh:
            self.fh.close()


class Fastq(object):
    """
    A class to hold features from fastq reads.
    """

    def __init__(self, name="", seq="", strand="+", qual="", conv=None):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        self.conv = conv
        self.i = int()
        assert isinstance(name, string_types)
        assert isinstance(seq, string_types)
        assert isinstance(qual, string_types)

    def __iter__(self):
        return self

    def next(self):
        if self.i < len(self):
            value, self.i = self[self.i], self.i + 1
            return value
        else:
            raise StopIteration()

    def __getitem__(self, key):
        if self.conv:
            return self.__class__(
                self.name,
                self.seq[key],
                self.strand,
                self.qual[key],
                self.conv[key],
            )
        else:
            return self.__class__(
                self.name, self.seq[key], self.strand, self.qual[key]
            )

    def __next__(self):
        return self.next()

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.name[0] != "@":
            self.name = "".join(["@", self.name])
        if self.conv:
            return (
                "\n".join(
                    [
                        "{0}:YM:Z:{1}".format(self.name, self.conv),
                        self.seq,
                        self.strand,
                        self.qual,
                    ]
                )
                + "\n"
            )
        else:
            return (
                "\n".join([self.name, self.seq, self.strand, self.qual]) + "\n"
            )

    def __len__(self):
        return len(self.seq)

    def gc(self):
        """Return the GC content of self as an int
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT', qual=' # # ##EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.gc()
        30
        """
        g = self.seq.count("G")
        c = self.seq.count("C")
        return int((g + c) / len(self) * 100)


class Sam(object):
    """ Store fields in each line of a SAM file, provided as a tuple. """

    __slots__ = [
        "qname",
        "flag",
        "rname",
        "pos",
        "mapq",
        "cigar",
        "rnext",
        "pnext",
        "tlen",
        "seq",
        "qual",
        "tags",
        "_tags",
        "_cigars",
    ]

    def __init__(self, fields):
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]
        self.tags = None
        self._tags = fields[11:]
        self._cigars = None

    def __gt__(self, other):
        if self.rname != other.rname:
            return self.rname > other.rname
        elif (self.rname == other.rname) and (self.pos != other.pos):
            return self.pos > other.pos
        else:
            return str(self) > str(other)

    def __lt__(self, other):
        if self.rname != other.rname:
            return self.rname < other.rname
        elif (self.rname == other.rname) and (self.pos != other.pos):
            return self.pos < other.pos
        else:
            return str(self) < str(other)

    def __eq__(self, other):
        if (
            (self.rname == other.rname)
            and (self.pos == other.pos)
            and (str(self) != str(other))
        ):
            return str(self) == str(other)
        else:
            return self.pos == other.pos

    def __str__(self):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return (
            "\t".join(
                (
                    self.qname,
                    str(self.flag),
                    self.rname,
                    str(self.pos),
                    str(self.mapq),
                    str(self.cigar),
                    self.rnext,
                    str(self.pnext),
                    str(self.tlen),
                    "".join(self.seq),
                    "".join(self.qual),
                )
                + tuple(
                    ":".join((tag, self.tags[tag][0], str(self.tags[tag][1])))
                    for tag in sorted(self.tags.keys())
                )
            )
            + "\n"
        )

    def __repr__(self):
        return "Sam({0}:{1}:{2})".format(self.rname, self.pos, self.qname)

    def __len__(self):
        return sum(
            c[0]
            for c in self.cigars
            if c[1] in ("M", "D", "N", "EQ", "X", "P")
        )

    def __getitem__(self, key):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return self.tags[key][1]

    def __setitem__(self, key, value):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        self.tags[key] = value

    def cigar_split(self):
        """CIGAR grouping function modified from:
        https://github.com/brentp/bwa-meth
        """
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    @property
    def conv(self):
        return self["YM"]

    @property
    def cigars(self):
        if not self._cigars:
            self._cigars = tuple(self.cigar_split())
        return self._cigars

    @property
    def mapped(self):
        return not (self.flag & 0x4)

    @property
    def secondary(self):
        return bool(self.flag & 0x100)

    @property
    def reverse(self):
        return bool(self.flag & 0x10)

    @property
    def duplicate(self):
        return bool(self.flag & 0x400)

    def gapped(self, string, gap_char="-"):
        """Return string with all deletions wrt reference
         represented as gaps '-' and all insertions wrt reference
         removed.
        i: sequence index
        """
        seq = []
        i = 0
        for n, t in self.cigars:
            if t in ("M", "N", "EQ", "X", "P"):
                seq.extend(string[i : i + n])
                i += n
            elif t in ("D",):
                seq.extend(("-",) * n)
            elif t in ("I",):
                i += n
        return "".join(seq)

    @property
    def coords(self):
        return range(self.pos, self.pos + len(self))


class FastaReader:
    """
    A class to read the name, sequence, strand and qualities from a fastq file
    """

    def __init__(self, f, ext=None):
        if ext == ".gz":
            self.file = Gzip(f.name, "r")
        else:
            self.file = f

    def __iter__(self) -> FastaReader:
        return self

    def __next__(self) -> Fastq:
        return self.next()

    def get_dummy_strand(self) -> str:
        return "+"

    def get_dummy_qual(self, seq) -> str:
        return "~" * len(seq)

    def next(self):
        try:
            name = next(self.file).strip().split()[0]  # remove whitespace
            seq = next(self.file).strip()
            strand = self.get_dummy_strand()
            qual = self.get_dummy_qual(seq)
            return Fastq(name=name, seq=seq, strand=strand, qual=qual)
        except StopIteration:
            raise StopIteration

    def subsample(self, n) -> Generator[Fastq, None, None]:
        """ Draws every nth read from self. Returns Fastq. """
        n = n * 4
        for i, line in enumerate(self.file):
            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
                yield Fastq(
                    name=name,
                    seq=seq,
                    strand=self.get_dummy_strand(),
                    qual=self.get_dummy_qual(seq),
                )

    def fileno(self):
        return self.file.fileno()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class FastqReader:
    """
    A class to read the name, sequence, strand and qualities from a fastq file
    """

    def __init__(self, f, ext=None):
        if ext == ".gz":
            self.file = Gzip(f.name, "r")
        else:
            self.file = f

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self) -> Fastq:
        try:
            name = next(self.file).strip().split()[0]  # remove whitespace
            seq = next(self.file).strip()
            strand = next(self.file).strip()
            qual = next(self.file).strip()
            if name.count(":YM:Z:") > 0:
                tag, dtype, data = name.split(":")[-3:]
                name = ":".join(name.split(":")[:-3])
                return Fastq(
                    name=name, seq=seq, strand=strand, qual=qual, conv=data
                )
            else:
                return Fastq(name=name, seq=seq, strand=strand, qual=qual)
        except StopIteration:
            raise StopIteration

    def subsample(self, n) -> Generator[Fastq, None, None]:
        """ Draws every nth read from self. Returns Fastq. """
        n = n * 4
        for i, line in enumerate(self.file):
            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
            elif i % n == 2:
                strand = line.strip()
            elif i % n == 3:
                qual = line.strip()
                if name.count(":YM:Z:") > 0:
                    tag, dtype, data = name.split(":")[-3:]
                    name = ":".join(name.split(":")[:-3])
                    yield Fastq(
                        name=name, seq=seq, strand=strand, qual=qual, conv=data
                    )
                else:
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual)

    def fileno(self):
        return self.file.fileno()

    def __enter__(self) -> FastqReader:
        return self

    def __exit__(self, *args):
        self.file.close()


class Reader(object):
    """ Read SAM/BAM format file using iterator. """

    def __init__(self, f):
        name, ext = os.path.splitext(f.name)
        if ext == ".bam":
            BamReaderSamtools.__init__(self, f)
        else:
            SamReader.__init__(self, f)

    def next(self) -> Sam:
        try:
            line = next(self.file).rstrip("\n\r")
            return Sam(tuple(line.split("\t")))
        except StopIteration:
            raise StopIteration

    def __next__(self) -> Sam:
        return self.next()

    def __iter__(self) -> Reader:
        return self

    def subsample(self, n) -> Generator[Sam, None, None]:
        """ Draws every nth read from self. Returns Sam. """
        for i, line in enumerate(self.file):
            if i % n == 0:
                yield Sam(tuple(line.rstrip().split("\t")))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class SamReader(Reader):
    """ Read SAM format file using iterator. """

    def __init__(self, f):
        self.header = []
        self.file = f
        for line in self.file:
            if line[0] == "@":
                self.header.append(line.rstrip("\n\r"))
            else:
                break


class BamReaderSamtools(Reader):
    """ Read BAM format file using iterator. """

    def __init__(self, f):
        pline = ["samtools", "view", "-H", f.name]
        try:
            p = Popen(pline, bufsize=-1, stdout=PIPE, stderr=PIPE)
        except OSError:
            sys.stderr.write(
                "Samtools must be installed for BAM file support!\n"
            )
            sys.exit(1)
        self.header = [
            line.decode("utf-8").rstrip("\n\r") for line in p.stdout
        ]
        p.wait()
        pline = ["samtools", "view", f.name]
        self.p = Popen(pline, bufsize=-1, stdout=PIPE, stderr=PIPE)
        self.file = TextIOWrapper(self.p.stdout)

    def __exit__(self, *args):
        self.p.wait()


class Stats:
    """Counter for characterization of NGS reads"""

    def __init__(self):
        self.depth = defaultdict(int)
        self.nuc = defaultdict(lambda: defaultdict(int))
        self.qual = defaultdict(lambda: defaultdict(int))
        self.gc = defaultdict(int)
        self.kmers = Counter(defaultdict(int))
        self.conv = defaultdict(lambda: defaultdict(int))

    def evaluate(self, seq: str, qual: str, conv: Optional[str] = None):
        """ Evaluate read object at each position, and fill in nuc and qual dictionaries """
        self.gc[gc(seq)] += 1
        if conv:
            cpgs = cpg_map(seq)
        for i in range(1, len(seq) + 1):
            self.depth[i] += 1
            self.nuc[i][seq[i - 1]] += 1
            self.qual[i][qual[i - 1]] += 1
            if conv:
                if cpgs[i - 1] != "N":
                    self.conv[i][conv[i - 1]] += 1

    def kmercount(self, seq: str, k: int = 5):
        """Count all kmers of k length in seq and update kmer counter."""
        for kmer in window(seq, n=k):
            self.kmers[kmer] += 1

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def bam_read_count(bamfile):
    """ Return a tuple of the number of mapped and unmapped reads in a bam file """
    p = Popen(["samtools", "idxstats", bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    return (mapped, unmapped)


def parse_sam_tags(tagfields):
    """ Return a dictionary containing the tags """
    return dict(
        (tag, (dtype, data))
        for tag, dtype, data in (decode_tag(x) for x in tagfields)
    )


def encode_tag(tag, data_type, data):
    """Write a SAM tag in the format ``TAG:TYPE:data``
    >>> encode_tag('YM', 'Z', '#""9O"1@!J')
    'YM:Z:#""9O"1@!J'
    """
    value = ":".join(list((tag.upper(), data_type.upper(), data)))
    return value


def decode_tag(tag):
    """Parse a SAM format tag to a (TAG, TYPE, data) tuple.

    TYPE in A, i, f, Z, H, B

    >>> decode_tag('YM:Z:#""9O"1@!J')
    ('YM', 'Z', '#""9O"1@!J')
    >>> decode_tag('XS:i:5')
    ('XS', 'i', 5)
    """
    values = tag.split(":")
    if len(values) != 3:
        values = (values[0], values[1], ":".join(values[2:]))
    if values[1] == "i":
        values[2] = int(values[2])
    elif values[1] == "f":
        values[2] = float(values[2])
    elif values[1] == "H":
        raise (
            NotImplementedError,
            "Hex array SAM tags are currently not parsed.",
        )
    elif values[1] == "B":
        raise (
            NotImplementedError,
            "Byte array SAM tags are currently not parsed.",
        )
    return tuple(values)


def gc(seq):
    """Return the GC content of as an int
    >>> x = tuple('TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT')
    >>> gc(x)
    30
    """
    g = seq.count("G")
    c = seq.count("C")
    return int((g + c) / len(seq) * 100)


def padbases(bases):
    """For each base call in dictionary D, add an entry base:0 if key base does not exist."""

    def inner(D):
        for key in bases:
            if key not in D:
                D[key] = 0

    return inner


def percentile(D, percent):
    """
    modified from: http://stackoverflow.com/a/2753343/717419
    Find the percentile of a list of values.

    N - is a dictionary with key=numeric value and value=frequency.
    percent - a float value from 0.0 to 1.0.

    outlier removal: http://www.purplemath.com/modules/boxwhisk3.htm

    return the percentile of the values
    """
    N = sorted(D.keys())  # dict keys
    P = [D[n] for n in N]  # dict values
    if not N:
        return None
    k = (sum(P)) * percent
    l = (sum(P)) * 0.25  # lower quartile
    u = (sum(P)) * 0.75  # upper quartile
    e = int()
    for n, p in zip(N, P):  # find percentile
        e += p
        if e >= k:
            z = n  # value at percentile
            break
    e = int()
    for n, p in zip(N, P):  # find upper quartile
        e += p
        if e >= u:
            uz = n  # value at quartile
            break
    e = int()
    for n, p in zip(N, P):  # find lower quartile
        e += p
        if e >= l:
            lz = n  # value at quartile
            break
    iqd = 1.5 * (uz - lz)  # 1.5 times the inter-quartile distance
    if (z) & (z < lz - iqd):
        return int(lz - iqd)
    elif (z) & (z > uz + iqd):
        return int(uz + iqd)
    elif z:
        return int(z)
    else:
        return N[-1]


def window(seq, n=2):
    """Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ..."""
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield "".join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield "".join(result)


def mean(s):
    return sum(s) / len(s)


def cpg_map(seq):
    """Return tuple of C/G/N.

    >>> cpg_map('CGCGTAGCCG')
    'CGCGNNNNCG'
    """
    starts = (x.start() for x in re.finditer("CG", "".join(["N", seq, "N"])))
    cpgs = ["N"] * (len(seq) + 2)
    for start in starts:
        cpgs[start] = "C"
        cpgs[start + 1] = "G"
    return "".join(cpgs[1:-1])


if __name__ == "__main__":
    import doctest

    doctest.testmod()
