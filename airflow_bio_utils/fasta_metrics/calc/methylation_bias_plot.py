import matplotlib.pyplot as plt

from .utils import DEFAULT_FIGURE_SETTINGS, CollectMetricsArgs


def methylation_bias_plot(settings: CollectMetricsArgs, positions, conv_dict):
    top_strand_c = conv_dict["C"]
    bot_strand_c = conv_dict["G"]
    methyl_values = [
        (top_strand_c[pos]["C"] + bot_strand_c[pos]["G"])
        / (
            top_strand_c[pos]["Y"]
            + bot_strand_c[pos]["R"]
            + top_strand_c[pos]["C"]
            + bot_strand_c[pos]["G"]
            + 0.1
        )
        for pos in positions
    ]
    fig, axes = plt.subplots(
        nrows=1, **(settings.figure_settings or DEFAULT_FIGURE_SETTINGS)
    )
    axes.plot(positions, methyl_values, color="red")
    x1, x2, y1, y2 = axes.axis()
    axes.axis((x1, x2, 0, 1))
    axes.yaxis.grid(
        b=True, which="major", **{"color": "gray", "linestyle": ":"}
    )
    axes.set_axisbelow(True)
    axes.set_title("Methylation bias (M-Bias)")
    axes.set_xlabel("Cycle")
    axes.set_ylabel("Bias")
    return "methylation_bias_plot", fig
