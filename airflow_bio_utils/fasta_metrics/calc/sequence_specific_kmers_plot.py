import itertools
from collections import defaultdict
from typing import Dict, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler
from plotly.subplots import make_subplots

from .colors import LEGEND_MPL_BGCOLOR_PARAM_NAME
from .utils import (DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs,
                    should_create_matplot, should_create_plotly)


def sequence_specific_kmers_plot(
    settings: CollectMetricsArgs,
    positions: List[int],
    counts: Dict[int, Dict[str, int]],
    top_kmers: List[str],
):
    all_kmers = [counts[k].keys() for k in sorted(counts.keys())]
    kmers = set(list(itertools.chain.from_iterable(all_kmers)))
    kmer_len = len(tuple(kmers)[0])
    if not top_kmers:
        top_kmers = [None]
    kmer_percent = defaultdict(lambda: defaultdict(int))
    for pos, count in tuple(counts.items()):
        for kmer in kmers:
            pos_count = sum(counts[pos].values())
            if pos_count > 0:
                kmer_percent[pos][kmer] = (
                    float(count.get(kmer, 0)) / pos_count * 100
                )
            else:
                kmer_percent[pos][kmer] = 0.0

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_plotly(settings):
        df_dict = dict(positions=positions)
        for kmer in [x for x in kmers if x not in top_kmers]:
            df_dict[kmer] = [kmer_percent[pos][kmer] for pos in positions]
        df = pd.DataFrame(df_dict)
        plotly_fig = make_subplots(
            x_title="Sequence length", y_title="Kmer content (% kmer)"
        )
        for kmer in [x for x in kmers if x not in top_kmers]:
            plotly_fig.add_scatter(
                x=df["positions"], y=df[kmer], mode="lines", name=kmer
            )

    if should_create_matplot(settings):
        cmap = mpl.cm.get_cmap(name="Set1")
        colors = [cmap(i) for i in np.linspace(0, 1, len(top_kmers))]
        mpl.rcParams["axes.prop_cycle"] = cycler(color=colors)
        matplot_fig, axes = plt.subplots(
            nrows=1,
            subplot_kw={LEGEND_MPL_BGCOLOR_PARAM_NAME: "white"},
            **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
        )
        for kmer in [x for x in kmers if x not in top_kmers]:
            axes.plot(
                positions,
                [kmer_percent[pos][kmer] for pos in positions],
                color="0.4",
                linestyle="dotted",
            )
        legend_handles = []
        for kmer in top_kmers:
            legend_handles.append(
                axes.plot(
                    positions, [kmer_percent[pos][kmer] for pos in positions]
                )[0]
            )
        # Shink current axis by 20%
        box = axes.get_position()
        axes.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        axes.yaxis.grid(
            b=True, which="major", **{"color": "gray", "linestyle": ":"}
        )
        axes.set_axisbelow(True)
        axes.set_title("Kmer ({0:n}) content".format(kmer_len))
        axes.set_xlabel("Cycle")
        axes.set_ylabel("Kmer content (% kmer)")
        legend = axes.legend(
            legend_handles,
            top_kmers,
            ncol=1,
            bbox_to_anchor=(1, 0.5),
            loc="center left",
            prop={"size": 8},
        )
        frame = legend.get_frame()
        frame.set_facecolor("white")
        for label in legend.get_texts():
            label.set_color("black")

    return "sequence_specific_kmers_plot", matplot_fig, plotly_fig, df
