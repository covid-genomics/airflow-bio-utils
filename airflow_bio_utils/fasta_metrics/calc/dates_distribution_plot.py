from .utils import CollectMetricsArgs, should_create_matplot, should_create_plotly, DEFAULT_FIGURE_SETTINGS

import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from typing import Dict
import datetime
import plotly.express as px
import scipy
import numpy as np
import math


def dates_distribution_plot(settings: CollectMetricsArgs, dates: Dict[datetime.datetime, int]):

    if len(dates) == 0:
        return None

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_matplot(settings):
        matplot_fig, axes = plt.subplots(nrows=1, **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS))
        axes.bar(list(dates.keys()), list(dates.values()), color=(0.1, 0.6, 0.8))
        axes.yaxis.grid(
            b=True, which='major', **{
                'color': 'gray',
                'linestyle': ':'
            })
        axes.set_title('Sequence dates distribution')
        axes.set_xlabel('Date')
        axes.set_ylabel('Number of sequences for date')

    if should_create_plotly(settings):
        df = pd.DataFrame.from_dict(dict(
            number_of_sequences=list(dates.values()),
            dates=list(dates.keys()),
        ))

        dfr = df.iloc[::-1]
        xmin = df[df.number_of_sequences > 0].iloc[0].dates
        xmax = dfr[dfr.number_of_sequences > 0].iloc[0].dates
        bin_step = (xmax - xmin) / 100
        plotly_fig = go.Figure()
        plotly_fig.add_trace(go.Histogram(
            x=df["dates"],
            y=df["number_of_sequences"],
            histfunc="sum",
            xbins=dict(start=xmin, end=xmax, size=bin_step)
        ))
        plotly_fig.update_yaxes(title_text="Number of sequences")
        plotly_fig.update_xaxes(title_text="Sequence date")

    return "dates_distribution_plot", matplot_fig, plotly_fig, df
