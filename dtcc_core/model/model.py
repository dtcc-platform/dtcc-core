# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

from abc import ABC, abstractmethod
import builtins
import importlib
from dataclasses import dataclass, field
from inspect import getmembers, isfunction, ismethod, ismodule
from google.protobuf.json_format import MessageToJson
from copy import deepcopy
from ..common import warning


@dataclass
class Model(ABC):
    """Base class for all DTCC Model classes."""

    @abstractmethod
    def to_proto(self):
        """
        Convert the model to its protobuf representation.

        Returns
        -------
        google.protobuf.message.Message
            Protobuf message encoding the model.
        """
        pass

    @abstractmethod
    def from_proto(self, pb):
        """
        Populate the model from a protobuf message.

        Parameters
        ----------
        pb : google.protobuf.message.Message or bytes
            Serialized or in-memory protobuf message representing the model.
        """
        pass

    def to_json(self) -> str:
        """Return a JSON representation of the object.

        Returns
        -------
            str
                A JSON string representing the object.
        """
        return MessageToJson(self.to_proto(), including_default_value_fields=True)

    def copy(self, **kwargs):
        """Return a copy of the object.

        Returns
        -------
            Model: A copy of the object.
        """
        c = deepcopy(self)
        for key, value in kwargs.items():
            setattr(c, key, value)
        return c

    @property
    def dataset_context(self):
        """Dataset v2 context attached to this object, when available."""
        return getattr(self, "_dataset_context", None)

    @dataset_context.setter
    def dataset_context(self, context):
        """Attach Dataset v2 context to this object."""
        self._dataset_context = context

    @property
    def metadata(self):
        """Dataset v2 factual metadata, when this object came from a dataset."""
        context = self.dataset_context
        return None if context is None else context.metadata

    @property
    def provenance(self):
        """Dataset v2 lineage information, when available."""
        context = self.dataset_context
        return None if context is None else context.provenance

    @property
    def presentation(self):
        """Dataset v2 presentation guidance, when available."""
        context = self.dataset_context
        return None if context is None else context.presentation

    def manifest(self):
        """Return a Dataset Manifest v2 snapshot for a dataset-produced object."""
        context = self.dataset_context
        if context is None:
            raise ValueError(
                "Cannot build a DatasetManifest: this object has no DatasetContext."
            )
        return context.manifest()

    def export(self, *args, **kwargs):
        """Export a Dataset v2 object package.

        Requires this object to carry ``DatasetContext`` from a dataset call.
        """
        from dtcc_core.datasets.package import export_model_package

        return export_model_package(self, *args, **kwargs)

    def publish(
        self,
        *,
        dataset_key: str,
        format: str | None = None,
        uploader=None,
        upload_url: str | None = None,
        token: str | None = None,
        idempotency_key: str | None = None,
    ):
        """Publish a Dataset v2 object package.

        Requires this object to carry ``DatasetContext`` from a dataset call.
        """
        from pathlib import Path
        import tempfile

        from dtcc_core.datasets.publish import DatasetUploadClient

        if self.dataset_context is None:
            raise ValueError(
                "Cannot publish Dataset v2 package: this object has no DatasetContext."
            )

        resolved_uploader = uploader or DatasetUploadClient.from_config(
            upload_url=upload_url,
            token=token,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            package = self.export(Path(tmpdir) / "dataset_package", format=format)
            return package.publish(
                dataset_key=dataset_key,
                uploader=resolved_uploader,
                idempotency_key=idempotency_key,
            )

    def info(self, print: bool = True, presentation: bool = True) -> str | None:
        """Print or return a human-readable summary of the model.

        Subclasses may override this for richer multi-line summaries. The base
        implementation intentionally follows ``str(self)`` so every model has a
        lightweight, uniform information API. Dataset-produced objects include
        their Dataset v2 presentation and metadata by default.
        """
        summary = str(self)
        if presentation and self.dataset_context is not None:
            from dtcc_core.datasets.presentation import format_dataset_context

            summary = (
                f"{summary}\n\n"
                f"{format_dataset_context(self.dataset_context, obj=self)}"
            )
        if print:
            builtins.print(summary)
            return None
        return summary

    def print_info(self, file=None) -> None:
        """Print the human-readable model summary returned by ``info()``."""
        builtins.print(self.info(print=False), file=file)

    def view(self, *args, **kwargs):
        """View the model using dtcc-viewer when available.

        The viewer is imported lazily so dtcc-core does not require graphical
        dependencies at import time. Importing ``dtcc_viewer`` registers
        object-specific ``view`` methods on DTCC model classes; this method then
        delegates to the registered implementation.
        """
        view_method = getattr(type(self), "view", None)
        if view_method is not None and view_method is not Model.view:
            return view_method(self, *args, **kwargs)

        try:
            importlib.import_module("dtcc_viewer")
        except Exception as exc:
            warning(
                f"Cannot view object: {self.__class__.__name__}. "
                "The dtcc-viewer module is not installed or graphical rendering "
                f"is not available ({exc})."
            )
            return None

        view_method = getattr(type(self), "view", None)
        if view_method is not None and view_method is not Model.view:
            return view_method(self, *args, **kwargs)

        warning(
            f"Cannot view object: {self.__class__.__name__}. "
            "No dtcc-viewer view method is registered for this model type."
        )
        return None

    @classmethod
    def add_methods(cls, module, name=None):
        """Adds methods from a module or function to the class.

        Parameters
        ----------
            module: module or function
                A function or module containing the methods to add.
            name : str
                The name of the method to add, if None use the
                function/method name (default None).

        Raises
        ------
            TypeError
             If module parameter is not a module or function.

        Returns
        -------
            None
        """
        # hack needed create a class variable for each subclass
        if not hasattr(cls, "_methods"):
            cls._methods = []

        if isfunction(module):
            if name is None:
                name = module.__name__
            _add_method(cls, module, name)
        elif ismodule(module):
            for function_name, function in getmembers(module, isfunction):
                if not function_name.startswith("_"):
                    _add_method(cls, function, function_name)
        else:
            raise TypeError("Expected a module or function")

    @classmethod
    def print_methods(cls, verbose=False):
        """
        Print the methods that have been added to the class.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print the docstring of each method (default False).

        Returns
        -------
        None
        """
        print(f"Methods for {cls.__name__}:")
        for name, parent_module, doc in cls._methods:
            print(f" - {name}: from {parent_module}")
            if verbose:
                print(f"   * {doc}")


def _add_method(cls, function, name=None):
    if name is None:
        name = function.__name__
    for idx, (function_name, _, _) in enumerate(cls._methods):
        if function_name == name:
            warning(
                f"{function} Method {function_name} already exists, replacing it."
            )
            cls._methods.pop(idx)
            break
    cls._methods.append((name, function.__module__, function.__doc__))
    setattr(cls, name, function)
