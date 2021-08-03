import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from .colors import LEGEND_MPL_BGCOLOR_PARAM_NAME
from .utils import DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs


def mismatch_plot(settings: CollectMetricsArgs, positions, counts):
    cmap = mpl.cm.get_cmap(name="Set1")
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=colors)
    fig, axes = plt.subplots(
        nrows=1,
        subplot_kw={LEGEND_MPL_BGCOLOR_PARAM_NAME: "white"},
        **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
    )
    ref_alt = []
    for ref in ["C", "G", "A", "T"]:
        for alt in {"C", "G", "A", "T"} - {ref}:
            axes.plot(
                positions,
                [
                    counts[ref][pos][alt] / sum(counts[ref][pos].values())
                    for pos in positions
                ],
            )
            ref_alt.append(">".join([ref, alt]))
    # Shink current axis by 20%
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width, box.height])
    axes.yaxis.grid(
        b=True, which="major", **{"color": "gray", "linestyle": ":"}
    )
    axes.set_axisbelow(True)
    axes.set_title("Reference mismatches by cycle")
    axes.set_xlabel("Cycle")
    axes.set_ylabel("Fraction of mismatches")
    legend = axes.legend(
        ref_alt,
        ncol=int(len(ref_alt) / 3),
        bbox_to_anchor=(0.5, 0.25),
        loc="best",
        prop={"size": 8},
    )
    frame = legend.get_frame()
    frame.set_facecolor("white")
    for label in legend.get_texts():
        label.set_color("black")
    return "mismatch_plot", fig
