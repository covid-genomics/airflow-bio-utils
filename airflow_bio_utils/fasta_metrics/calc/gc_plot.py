from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

from .utils import (DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs,
                    should_create_matplot, should_create_plotly)


def gc_plot(
    settings: CollectMetricsArgs,
    positions: List[int],
    pos_gc: List[float],
):
    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_matplot(settings):
        matplot_fig, axes = plt.subplots(
            nrows=1, **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
        )
        axes.plot(positions, pos_gc, color=(0.1, 0.6, 0.8))
        x1, x2, y1, y2 = axes.axis()
        axes.axis((x1, x2, 0, 100))
        axes.yaxis.grid(
            b=True, which="major", **{"color": "gray", "linestyle": ":"}
        )
        axes.set_axisbelow(True)
        axes.set_title("GC content distribution")
        axes.set_xlabel("Cycle")
        axes.set_ylabel("GC% at cycle")
    if should_create_plotly(settings):
        df = pd.DataFrame.from_dict(
            dict(
                gc_content=[val / 100 for val in pos_gc],
                sequence_length=positions,
            )
        )
        plotly_fig = px.line(df, x="sequence_length", y="gc_content")
        plotly_fig.update_yaxes(title_text="GC Content (%)", range=[0, 1])
        plotly_fig.update_xaxes(title_text="Sequence Length")
        plotly_fig.update_layout(yaxis_tickformat="%")

    return "gc_plot", matplot_fig, plotly_fig, df
