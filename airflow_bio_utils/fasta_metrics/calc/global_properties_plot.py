import datetime
import statistics
from typing import Dict, List, Set

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .utils import (CollectMetricsArgs, should_create_matplot,
                    should_create_plotly)


def global_properties_plot(
    settings: CollectMetricsArgs,
    read_len: Dict[int, int],
    seqs_count: int,
    seqs_with_n: int,
    seqs_with_any_non_actgn_symbol: int,
    non_actgn_symbols: Set[str],
    dates: Dict[datetime.datetime, int],
    duplicates: int,
    pos_gc: List[float],
    cycle_qual: Dict[int, Dict[int, int]],
    has_quality_data: bool,
):
    plotly_fig = None
    df = None

    min_read = min(read_len.keys())
    max_read = max(read_len.keys())

    has_dates = len(dates) > 0

    if has_dates:
        min_date, max_date = min(dates.keys()), max(dates.keys())
    else:
        min_date, max_date = None, None

    quals = [
        qual
        for pos in cycle_qual.keys()
        for qual in cycle_qual[pos].keys()
        for i in range(cycle_qual[pos][qual])
    ]

    properties = dict(
        sequences_count=(seqs_count, False, True),
        shortest_sequence=(min_read, "count", read_len[min_read]),
        longest_sequence=(max_read, "count", read_len[max_read]),
        average_length=(
            statistics.mean(
                [
                    seq_len
                    for seq_len in read_len.keys()
                    for i in range(read_len[seq_len])
                ]
            ),
            False,
            None,
        ),
        oldest_sequence=(min_date, "count", dates[min_date])
        if has_dates
        else None,
        newest_sequence=(max_date, "count", dates[max_date])
        if has_dates
        else None,
        average_age_days=(
            int(
                statistics.mean(
                    [
                        (datetime.datetime.today() - seq_date).days
                        for seq_date in dates.keys()
                        for i in range(dates[seq_date])
                    ]
                )
            ),
            False,
            None,
        )
        if has_dates
        else None,
        sequences_with_dates_count=(sum(dates.values()), True, True),
        n_sequences_count=(seqs_with_n, True, True),
        duplicates_count=(duplicates, True, True),
        non_actgn_sequences_count=(seqs_with_any_non_actgn_symbol, True, True),
        non_actgn_symbols=(
            list(non_actgn_symbols),
            False,
            len(non_actgn_symbols),
        ),
        average_gc_content=(statistics.mean(pos_gc), False, None),
        average_quality=(statistics.mean(quals), False, None)
        if has_quality_data
        else (None, False, None),
        best_quality=(max(quals), False, None)
        if has_quality_data
        else (None, False, None),
        worst_quality=(min(quals), False, None)
        if has_quality_data
        else (None, False, None),
    )

    properties = {
        key: value for key, value in properties.items() if value is not None
    }

    if should_create_plotly(settings):
        vals = dict(
            property_name=list(properties.keys()),
            property_value=list([val[0] for val in properties.values()]),
            percentage_of_total=[
                (
                    (val[2] if val[1] is "count" else val[0])
                    / seqs_count
                    * 100
                    if val[1] is not False
                    else None
                )
                for val in properties.values()
            ],
            count=[
                (val[0] if val[2] is True else val[2])
                for val in properties.values()
            ],
        )
        df = pd.DataFrame.from_dict(vals)
        plotly_fig = go.Figure(
            data=[
                go.Table(
                    header=dict(
                        values=list(df.columns),
                        fill_color="paleturquoise",
                        align="left",
                        font_size=7,
                    ),
                    cells=dict(
                        values=[df[key] for key in vals.keys()],
                        fill_color="lavender",
                        align="left",
                        font_size=7,
                    ),
                )
            ]
        )

    return "global_properties_plot", None, plotly_fig, df
