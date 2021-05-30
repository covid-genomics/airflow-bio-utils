from .utils import CollectMetricsArgs, should_create_matplot, should_create_plotly, DEFAULT_FIGURE_SETTINGS

import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import scipy
import numpy as np
import math
from typing import Dict


def gc_distribution_plot(
    settings: CollectMetricsArgs,
    counts: Dict[int, int],
):
    m = int(
        sum([k * v for k, v in zip(counts.keys(), counts.values())]) / sum(
            counts.values()))
    variances = [(k - m)**2 for k in counts.keys()]
    variance = int(
        sum([k * v for k, v in zip(variances, counts.values())]) / sum(
            counts.values()))
    sigma = math.sqrt(variance)
    x = np.linspace(0, 100, 100)
    gc_content, samples_count = zip(*sorted(counts.items(), key=lambda x: x[0]))
    gc_content_ref = scipy.stats.norm.pdf(np.linspace(0, 101, 101), m, sigma)

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_matplot(settings):
        matplot_fig, axes = plt.subplots(nrows=1, **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS))
        axes.scatter(gc_content, samples_count, color=(0.1, 0.6, 0.8), marker="o")
        x1, x2, y1, y2 = axes.axis()
        axes.axis((x1, x2, 0, y2))
        axes2 = axes.twinx()
        axes2.plot(x, scipy.stats.norm.pdf(x, m, sigma), color='red')
        axes2.get_yaxis().set_visible(False)
        handles, labels = axes.get_legend_handles_labels()
        display = (0, 1, 2)
        a = plt.Line2D((0, 1), (0, 0), color=(0.1, 0.6, 0.8))
        b = plt.Line2D((0, 1), (0, 0), color='red')
        axes.legend([handle for i, handle in enumerate(handles)
                     if i in display] + [a, b],
                    [label for i, label in enumerate(labels)
                     if i in display] + ['Actual', 'Theoretical'])
        axes.yaxis.grid(
            b=True, which='major', **{
                'color': 'gray',
                'linestyle': ':'
            })
        axes.set_axisbelow(True)
        axes.set_title('Read GC content distribution')
        axes.set_xlabel('Mean read GC content (%)')
        axes.set_ylabel('Number of reads')

    if should_create_plotly(settings):
        xint = [int(val) for val in np.linspace(0, 100, 100)]
        samples_sum = sum(samples_count)

        df = pd.DataFrame.from_dict(dict(
            samples_count=[(samples_count[gc_content.index(val)] if val in gc_content else 0) for val in xint],
            gc_content=[val/100 for val in xint],
            gc_content_ref=[val*samples_sum for val in gc_content_ref[1:]],
        ))

        plotly_fig = px.scatter(df, x="gc_content", y="samples_count")
        plotly_fig.add_trace(
            go.Scatter(x=df["gc_content"], y=df["gc_content_ref"],
                       name='Reference content (normal)',
                       mode='lines',
                       line_color='red'))
        plotly_fig.update_layout(xaxis_tickformat='%')

    return "gc_distribution_plot", matplot_fig, plotly_fig, df
