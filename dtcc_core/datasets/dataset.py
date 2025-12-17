import dtcc_core
from abc import ABC, abstractmethod
from typing import Any, Dict
from pydantic import BaseModel


class DatasetDescriptor(ABC):
    """Callable, self-describing dataset."""

    name: str
    ArgsModel: BaseModel

    def __call__(self, **kwargs):
        args = self.validate(kwargs)
        return self.build(args)

    @abstractmethod
    def validate(self, kwargs):
        """Validate and normalize user arguments."""
        raise NotImplementedError

    @abstractmethod
    def build(self, validated_args):
        """Resolve the dataset and return the result."""
        raise NotImplementedError

    def show_options(self):
        return self.ArgsModel.model_json_schema()
