from collections import defaultdict
from typing import Dict, List, Set

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler
from plotly.subplots import make_subplots

from .colors import LEGEND_MPL_BGCOLOR_PARAM_NAME
from .utils import (DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs,
                    should_create_matplot, should_create_plotly)


def cycle_spec_composition_plot(
    settings: CollectMetricsArgs,
    positions: List[int],
    nucs: Set[str],
    counts: Dict[int, Dict[str, int]],
):
    nuc_order = [
        "A",
        "T",
        "C",
        "G",
        "N",
        "M",
        "R",
        "W",
        "S",
        "Y",
        "K",
        "V",
        "H",
        "D",
        "B",
    ]
    max_depth = sum(tuple(counts[1].values()))
    nuc_percent = defaultdict(lambda: defaultdict(int))
    for pos, count in tuple(counts.items()):
        max_depth = sum(tuple(count.values()))
        for nuc in nucs:
            if max_depth > 0:
                nuc_percent[pos][nuc] = float(count[nuc]) / max_depth * 100
            else:
                nuc_percent[pos][nuc] = 0.0

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_plotly(settings):
        df_dict = dict(
            sequence_positions=positions,
        )
        for nuc in nuc_order:
            if nuc in nucs:
                df_dict[nuc] = [
                    nuc_percent[pos][nuc] / 100 for pos in positions
                ]
        df = pd.DataFrame.from_dict(df_dict)

        plotly_fig = make_subplots(
            x_title="Sequence length", y_title="Base content (% basecall)"
        )
        for nuc in nuc_order:
            if nuc in nucs:
                plotly_fig.add_scatter(
                    x=df["sequence_positions"],
                    y=df[nuc],
                    mode="lines",
                    name=nuc,
                )
        plotly_fig.update_layout(yaxis_tickformat="%")

    if should_create_matplot(settings):
        cmap = mpl.cm.get_cmap(name="Set1")
        colors = [cmap(i) for i in np.linspace(0, 1, len(nuc_order))]
        mpl.rcParams["axes.prop_cycle"] = cycler(color=colors)
        matplot_fig, axes = plt.subplots(
            nrows=1,
            subplot_kw={LEGEND_MPL_BGCOLOR_PARAM_NAME: "white"},
            **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
        )
        for nuc in nuc_order:
            if nuc in nucs:
                axes.plot(
                    positions, [nuc_percent[pos][nuc] for pos in positions]
                )
        # Shink current axis by 20%
        box = axes.get_position()
        axes.set_position([box.x0, box.y0, box.width, box.height])
        axes.yaxis.grid(
            b=True, which="major", **{"color": "gray", "linestyle": ":"}
        )
        axes.set_axisbelow(True)
        axes.set_title("Base content")
        axes.set_xlabel("Cycle")
        axes.set_ylabel("Base content (% basecall)")
        legend = axes.legend(
            tuple((n for n in nuc_order if n in nucs)),
            ncol=len(nucs),
            bbox_to_anchor=(0.5, 0.25),
            loc="center",
            prop={"size": 8},
        )
        frame = legend.get_frame()
        frame.set_facecolor("white")
        for label in legend.get_texts():
            label.set_color("black")

    return "cycle_spec_composition_plot", matplot_fig, plotly_fig, df
