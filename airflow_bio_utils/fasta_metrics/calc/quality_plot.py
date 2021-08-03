from typing import List

import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots

from .utils import (DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs,
                    should_create_matplot, should_create_plotly)


def quality_plot(
    settings: CollectMetricsArgs,
    positions: List[int],
    quantiles: List[List[int]],
):
    Q0 = [q[0] for q in quantiles]
    Q1 = [q[1] for q in quantiles]
    Q2 = [q[2] for q in quantiles]
    Q3 = [q[3] for q in quantiles]
    Q4 = [q[4] for q in quantiles]

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_plotly(settings):
        Qs = dict(
            Q4="75%-100%",
            Q1="0%-25%",
            Q2="Median",
            Q0="25%-50%",
            Q3="50%-75%",
        )
        df = pd.DataFrame(
            dict(positions=positions, Q0=Q0, Q1=Q1, Q2=Q2, Q3=Q3, Q4=Q4)
        )
        plotly_fig = make_subplots(
            x_title="Sequence length", y_title="Phred score (percentile)"
        )
        for Q in Qs.keys():
            plotly_fig.add_scatter(
                x=df["positions"], y=df[Q], mode="lines", name=Qs[Q]
            )

    if should_create_matplot(settings):
        matplot_fig = plt.figure(
            **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
        )
        ax = matplot_fig.add_subplot(111)
        matplot_fig.subplots_adjust(top=0.85)
        ax.fill_between(positions, Q4, color=(1.0, 0.8, 0.0))
        ax.fill_between(positions, Q3, color=(1.0, 0.6, 0.0))
        ax.fill_between(positions, Q2, color=(1.0, 0.4, 0.0))
        ax.fill_between(positions, Q1, color=(1.0, 0.2, 0.0))
        ax.fill_between(positions, Q0, color=(1.0, 1.0, 1.0))
        ax.plot(positions, Q4, color=(1.0, 0.8, 0.0))
        ax.plot(positions, Q3, color=(1.0, 0.6, 0.0))
        ax.plot(positions, Q2, color=(1.0, 0.4, 0.0))
        ax.plot(positions, Q1, color=(1.0, 0.2, 0.0))
        ax.plot(positions, Q2, color="black")
        x1, x2, y1, y2 = ax.axis()
        ax.axis((x1, x2, 0, max(Q4)))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height])
        ax.yaxis.grid(
            b=True, which="major", **{"color": "gray", "linestyle": ":"}
        )
        ax.legend(
            ("100-75%", "75-50%", "50-25%", "25-0%", "Median"),
            bbox_to_anchor=(0, 0.25),
            loc="center left",
        )
        ax.set_title("Quality score percentiles")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("Phred score")

    return "quality_plot", matplot_fig, plotly_fig, df
