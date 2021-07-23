from airflow_bio_utils.sequences.filter.condition import FilterCondition
from Bio.SeqRecord import SeqRecord
from typing import Optional, List, Union, Literal

DEFAULT_ACCEPTED_SEQUENCE_SYMBOLS = ["A", "C", "T", "G", "N"]


class FilterSymbols(FilterCondition):
    # Fields to apply Airflow templates
    template_fields = ["allowed_symbols"]

    def __init__(
        self,
        allowed_symbols: Optional[List[str]] = None,
    ):
        super().__init__()
        if allowed_symbols is None:
            allowed_symbols = DEFAULT_ACCEPTED_SEQUENCE_SYMBOLS
        self.allowed_symbols = allowed_symbols

    def check(self, record: SeqRecord) -> bool:
        return all((symbol in self.allowed_symbols) for symbol in record)


class FilterSymbolQuantity(FilterCondition):
    # Fields to apply Airflow templates
    template_fields = ["symbol", "quantity", "is_percentage", "compare_method"]

    symbol: str
    quantity: float
    is_percentage: bool
    compare_method: Literal["LT", "GT"]

    def __init__(
        self,
        symbol: str,
        quantity: float,
        compare_method: Literal["LT", "GT"],
        is_percentage: bool = True,
    ):
        super().__init__()
        self.symbol = symbol
        self.quantity = quantity
        self.compare_method = compare_method
        self.is_percentage = is_percentage

    def check(self, record: SeqRecord) -> bool:
        symbol_count = sum([1 for nucl in record if nucl == self.symbol])
        if self.is_percentage:
            symbol_count = symbol_count / len(record) * 100
        if self.compare_method == "LT":
            return symbol_count <= self.quantity
        if self.compare_method == "GT":
            return symbol_count >= self.quantity

