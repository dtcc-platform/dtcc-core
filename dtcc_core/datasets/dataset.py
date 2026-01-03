import dtcc_core
from dtcc_core.model import Bounds
from dtcc_core.model import Object as DTCCObject
from dtcc_core.model import Geometry as DTCCGeometry

from abc import ABC, abstractmethod
from typing import Any, Sequence, Union
from pydantic import BaseModel, Field, field_validator
from pathlib import Path
import tempfile


class DatasetBaseArgs(BaseModel):
    bounds: Sequence[float] = Field(
        ...,
        description="Bounding box [minx, miny, maxx, maxy] or [minx, miny, minz, maxx, maxy, maxz]",
    )

    @field_validator("bounds")
    @classmethod
    def validate_bounds(cls, v: Sequence[float]):
        """Validate bounds length and order."""
        if len(v) not in (4, 6):
            raise ValueError("Bounds must be 4 or 6 floats")

        # Validate order
        if len(v) == 4:
            if v[0] >= v[2] or v[1] >= v[3]:
                raise ValueError("Invalid bounds: xmin < xmax, ymin < ymax")
        elif len(v) == 6:
            if v[0] >= v[3] or v[1] >= v[4] or v[2] >= v[5]:
                raise ValueError("Invalid bounds: min < max for all dimensions")
        return v


class DatasetDescriptor(ABC):
    """Callable, self-describing dataset."""

    name: str
    description: str = ""
    ArgsModel: BaseModel

    def __init_subclass__(cls, register=True, **kwargs):
        """
        Auto-register dataset subclasses when they're defined.

        Args:
            register: Whether to auto-register this dataset (default: True).
                     Set to False for abstract base classes.
            **kwargs: Additional keyword arguments passed to super().__init_subclass__
        """
        super().__init_subclass__(**kwargs)

        # Only register if:
        # - registration is enabled (register=True)
        # - class has a name attribute
        # - name is not empty
        if register and hasattr(cls, 'name') and cls.name:
            from dtcc_core.datasets.registry import _register_dataset_class
            _register_dataset_class(cls.name, cls)

    def __call__(self, **kwargs):
        args = self.validate(kwargs)
        return self.build(args)

    def validate(self, kwargs):
        args = self.ArgsModel(**kwargs)
        return args

    @abstractmethod
    def build(self, validated_args):
        """Resolve the dataset and return the result."""
        raise NotImplementedError

    def show_options(self):
        return self.ArgsModel.model_json_schema()

    @staticmethod
    def parse_bounds(bounds: Sequence[float]) -> Bounds:
        """Convert bounds list to a Bounds object.

        Args:
            bounds: [minx, miny, maxx, maxy] or
                   [minx, miny, minz, maxx, maxy, maxz]

        Returns:
            Bounds object
        """
        if len(bounds) == 4:
            return Bounds(
                xmin=bounds[0], ymin=bounds[1], xmax=bounds[2], ymax=bounds[3]
            )
        elif len(bounds) == 6:
            return Bounds(
                xmin=bounds[0],
                ymin=bounds[1],
                zmin=bounds[2],
                xmax=bounds[3],
                ymax=bounds[4],
                zmax=bounds[5],
            )
        else:
            raise ValueError(f"Bounds must be 4 or 6 floats, got {len(bounds)}")

    @staticmethod
    def export_to_bytes(
        obj: Union[DTCCObject, DTCCGeometry], format: str, as_text=False, **save_kwargs
    ) -> Union[bytes, str]:
        """Export object to bytes.

        Args:
            obj: Object with .save() method
            format: File format extension
            as_text: Return as text instead of bytes
            **save_kwargs: Passed to obj.save()

        Returns:
            File contents as bytes
        """
        with tempfile.NamedTemporaryFile(suffix=f".{format}", delete=True) as tmpfile:
            obj.save(tmpfile.name, **save_kwargs)
            if as_text:
                return Path(tmpfile.name).read_text()
            else:
                return Path(tmpfile.name).read_bytes()
