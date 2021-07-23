from abc import ABCMeta, abstractmethod
from Bio.SeqRecord import SeqRecord


class FilterCondition(metaclass=ABCMeta):

    @abstractmethod
    def check(self, record: SeqRecord) -> bool:
        raise Exception("Method not implemented: FilterCondition.check(). FilterCondition is abstract class")
