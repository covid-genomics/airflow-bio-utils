from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px

from .colors import VIRDIS_COLOR_MAP
from .utils import (DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs, Counter,
                    should_create_matplot, should_create_plotly)


def quality_map_plot(
    settings: CollectMetricsArgs,
    qualities: Dict[int, Dict[int, int]],
):
    values = map(Counter, tuple(qualities.values()))
    counts = Counter()
    for value in values:
        counts = counts + value
    max_qual = max(tuple(counts.keys()))
    max_pos = max(tuple(qualities.keys()))
    heat_map = np.zeros((max_qual, max_pos))
    for p in range(max_pos):
        for q in range(max_qual):
            try:
                heat_map[q][p] = qualities[p + 1][q + 1]
            except KeyError:
                pass

    matplot_fig = None
    plotly_fig = None
    df = None
    if should_create_plotly(settings):
        df = pd.DataFrame(np.array(heat_map))
        plotly_fig = px.imshow(
            np.array(
                [[np.array(i, dtype=np.uint8) for i in j] for j in df.values],
                dtype=np.uint8,
            ),
            labels=dict(
                x="Sequence length",
                y="Sum of phred qualities",
                color="Number of scores",
            ),
        )
        plotly_fig.update_yaxes(autorange=True)

    if should_create_matplot(settings):
        matplot_fig = plt.figure(
            **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
        )
        ax = matplot_fig.add_subplot(111)
        imax = ax.imshow(
            np.array(heat_map),
            cmap=VIRDIS_COLOR_MAP,
            origin="lower",
            interpolation="none",
            aspect="auto",
        )
        ax.axhline(y=10, linestyle=":", color="gray")
        ax.axhline(y=20, linestyle=":", color="gray")
        ax.axhline(y=30, linestyle=":", color="gray")
        cbar = matplot_fig.colorbar(imax, orientation="horizontal", shrink=0.5)
        cbar_labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
        cbar.ax.set_xticklabels(cbar_labels, rotation=45)
        cbar.ax.set_title("")
        ax.set_title("Quality score heatmap")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("Sum of Phred qualities")

    return "quality_map_plot", matplot_fig, plotly_fig, df
