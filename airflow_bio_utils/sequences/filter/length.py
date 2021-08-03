from typing import Optional

from Bio.SeqRecord import SeqRecord

from airflow_bio_utils.sequences.filter.condition import FilterCondition


class FilterByLength(FilterCondition):
    # Fields to apply Airflow templates
    template_fields = ["min_length", "max_length"]

    def __init__(
        self,
        min_length: Optional[int],
        max_length: Optional[int],
    ):
        super().__init__()
        self.min_length = min_length
        self.max_length = max_length

        if self.max_length is None:
            self.max_length = float("inf")
        if self.min_length is None:
            self.min_length = 0

    def check(self, record: SeqRecord) -> bool:
        return self.max_length >= len(record) >= self.min_length
