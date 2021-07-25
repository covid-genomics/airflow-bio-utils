from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .utils import (DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs,
                    should_create_matplot, should_create_plotly)


def sequence_length_distribution_plot(
    settings: CollectMetricsArgs,
    lengths: Dict[int, int],
):
    max_len = max(lengths.keys())
    positions = tuple(range(max_len + 1))

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_matplot(settings):
        matplot_fig, axes = plt.subplots(
            nrows=1, **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
        )
        axes.bar(
            positions, [lengths[n] for n in positions], color=(0.1, 0.6, 0.8)
        )
        axes.yaxis.grid(
            b=True, which="major", **{"color": "gray", "linestyle": ":"}
        )
        axes.set_title("Read lengths distribution")
        axes.set_xlabel("Cycle")
        axes.set_ylabel("Number of reads at cycle")

    if should_create_plotly(settings):
        df = pd.DataFrame.from_dict(
            dict(
                number_of_sequences=[lengths[n] for n in positions],
                sequence_length=positions,
            )
        )
        # plotly_fig = px.bar(df, x="sequence_length", y="number_of_sequences")
        dfr = df.iloc[::-1]
        xmin = df[df.number_of_sequences > 0].iloc[0].sequence_length
        xmax = dfr[dfr.number_of_sequences > 0].iloc[0].sequence_length
        bin_step = (xmax - xmin) / 100
        bins = dict(start=xmin, end=xmax, size=bin_step)
        if xmax - xmin <= 50:
            bins = None
            start = max(0, min(xmin, xmax) - 50)
            end = max(xmin, xmax) + 50
            xmin = start
            xmax = end
        if bins is not None:
            plotly_fig = go.Figure()
            plotly_fig.add_trace(
                go.Histogram(
                    x=df["sequence_length"],
                    y=df["number_of_sequences"],
                    histfunc="sum",
                    xbins=bins,
                )
            )
        else:
            plotly_fig = px.bar(
                df[
                    (df.sequence_length >= xmin) & (df.sequence_length <= xmax)
                ],
                x="sequence_length",
                y="number_of_sequences",
            )
        # plotly_fig = px.histogram(df, x="sequence_length", y="number_of_sequences", histfunc="sum", xbins=bins)
        plotly_fig.update_yaxes(title_text="Number of sequences")
        plotly_fig.update_xaxes(title_text="Sequence Length")

    return "sequence_length_distribution_plot", matplot_fig, plotly_fig, df
