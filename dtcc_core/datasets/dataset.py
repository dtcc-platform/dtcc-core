import dtcc_core
from abc import ABC, abstractmethod
from typing import Any, Dict


class DatasetDescriptor(ABC):
    """Callable, self-describing dataset."""

    name: str

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

    @abstractmethod
    def show_options(self) -> Dict[str, Any]:
        """Return a structured description of valid arguments."""
        raise NotImplementedError
