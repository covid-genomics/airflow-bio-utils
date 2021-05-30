from .utils import CollectMetricsArgs, DEFAULT_FIGURE_SETTINGS, should_create_matplot, should_create_plotly, Counter

import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np
import itertools
from typing import Dict, ValuesView


def quality_distribution_plot(
    settings: CollectMetricsArgs,
    qualities: ValuesView[Dict[int, int]],
):
    values = map(Counter, qualities)
    counts = Counter()
    for value in values:
        counts += value
    phred_score, cumulative_sum = zip(*sorted(counts.items(), key=lambda x: x[0]))

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_plotly(settings):
        df = pd.DataFrame(dict(cumulative_sum=cumulative_sum, phred_score=phred_score))
        plotly_fig = px.bar(df, x='phred_score', y='cumulative_sum')
        plotly_fig.update_yaxes(title_text="Cumulative sum")
        plotly_fig.update_xaxes(title_text="Phred score")

    if should_create_matplot(settings):
        matplot_fig = plt.figure(**(settings.figure_settings or DEFAULT_FIGURE_SETTINGS))
        ax = matplot_fig.add_subplot(111)
        matplot_fig.subplots_adjust(top=0.85)
        ax.bar([_ - 1 for _ in phred_score], cumulative_sum, color=(0.1, 0.6, 0.8))
        ax.yaxis.grid(b=True, which='major', **{'color': 'gray', 'linestyle': ':'})
        ax.set_axisbelow(True)
        x1, x2, y1, y2 = ax.axis()
        ax.axis((x1, max(phred_score), 0, y2))
        ax.set_title('Quality score cumulative distribution')
        ax.set_xlabel('Phred score')
        ax.set_ylabel('Cumulative sum')

    return "quality_distribution_plot", matplot_fig, plotly_fig, df, np.median(tuple(itertools.chain(*([n] * m for n, m in counts.items()))))
