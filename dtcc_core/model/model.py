# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

from abc import ABC
import builtins
import importlib
from dataclasses import dataclass, field
from inspect import getmembers, isfunction, ismethod, ismodule
from google.protobuf.json_format import MessageToJson
from copy import deepcopy
from ..common import warning


@dataclass(repr=False)
class Model(ABC):
    """Base class for all DTCC Model classes."""

    def _summary_items(self):
        """Cheap scalar facts for compact display; subclasses may extend these."""
        return []

    def _repr_is_complete(self):
        return False

    def __repr__(self):
        from ..common._display import format_repr

        return format_repr(
            type(self).__name__,
            self._summary_items(),
            complete=self._repr_is_complete(),
        )

    def __str__(self):
        return repr(self)

    def tree(self, indent="", geometry_type=None, *, verbose=False, max_depth=None):
        """Print the model hierarchy using tree-style branches.

        The default shows objects and geometry attachments. ``verbose=True``
        adds attributes, attachment metadata and field descriptions.
        ``max_depth`` limits nesting (root is zero); ``...`` marks omitted
        descendants. Raster and Bounds attachments are leaves.
        """
        from ._display import print_tree

        if not isinstance(verbose, bool):
            raise TypeError("verbose must be a bool")
        if max_depth is not None:
            if isinstance(max_depth, bool) or not isinstance(max_depth, int):
                raise TypeError("max_depth must be a nonnegative integer or None")
            if max_depth < 0:
                raise ValueError("max_depth must be nonnegative")
        label = repr(self)
        if geometry_type is not None:
            label = f"{geometry_type}: {label}"
        print_tree(self, label, indent=indent, verbose=verbose, max_depth=max_depth)

    def _info_sections(self):
        """Detailed report content, formatted and emitted by Model.info."""
        from ..common._display import label, value_text

        rows = [(label(key), value_text(value)) for key, value in self._summary_items()]
        if self.schema_id is not None:
            rows.append(("Schema ID", self.schema_id))
        if self.schema_version is not None:
            rows.append(("Schema version", self.schema_version))
        return [("", ("Property", "Value"), rows)]

    def plot(
        self,
        ax=None,
        *,
        lod=None,
        representation=None,
        field=None,
        max_elements=20000,
        theme="dark",
        show=True,
    ):
        """Quick 3D Matplotlib preview; return the axes for further customization.

        Objects select one representation each and traverse their children.
        Use ``representation`` (attachment ID), exact ``lod`` or a ``field`` name
        to inspect a particular part of the model. Fields appear as coloured
        samples; vector fields use magnitude. Large inputs are sampled within
        ``max_elements``. This is a preview, not a full scene renderer.

        See docs/model-preview.md for selection, geometry and coordinate limits.
        """
        from ..plotting.model import _plot_model

        return _plot_model(
            self,
            ax=ax,
            lod=lod,
            representation=representation,
            field=field,
            max_elements=max_elements,
            theme=theme,
            show=show,
        )

    @property
    def schema_id(self):
        """Root semantic schema ID; None selects the bundled default.

        This is root I/O metadata, independent of Object domain-profile labels.
        It is consulted at canonical and strict CityJSON I/O boundaries.
        Only the canonical format persists the declaration.
        """
        return getattr(self, "_schema_id", None)

    @schema_id.setter
    def schema_id(self, value):
        self._schema_id = value

    @property
    def schema_version(self):
        """Root semantic schema version, paired with schema_id."""
        return getattr(self, "_schema_version", None)

    @schema_version.setter
    def schema_version(self, value):
        self._schema_version = value

    def to_proto(self, *, validate_schema=True):
        """Return the ModelFile message defined by dtcc.proto.

        Uses the same admission and default schema validation as .dtcc files.
        """
        from .exchange import _encode_model

        return _encode_model(self, validate_schema=validate_schema)

    def from_proto(self, pb, *, validate_schema=True):
        """Replace this model from a ModelFile message or its binary bytes.

        Invalid data or a different concrete root type leaves this model intact.
        """
        from .exchange import _decode_model, SUPPORTED_ROOTS

        if type(self) not in SUPPORTED_ROOTS:
            raise NotImplementedError(
                f"Protobuf model {type(self).__name__} is unsupported"
            )
        restored, _ = _decode_model(pb, validate_schema=validate_schema)
        if type(restored) is not type(self):
            raise ValueError(
                f"Expected {type(self).__name__}, received {type(restored).__name__}"
            )
        self.__dict__.clear()
        self.__dict__.update(restored.__dict__)

    def to_json(self) -> str:
        """Return a JSON representation of the object.

        Returns
        -------
            str
                A JSON string representing the object.
        """
        return MessageToJson(self.to_proto(), always_print_fields_with_no_presence=True)

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
        """Export an object package with its attached ``DatasetContext``.

        Pass ``canonical=True`` for a canonical v3 package of the supported model
        subset. ``format`` then selects an optional supplemental artifact. The
        default remains legacy v2 while producers and consumers migrate.
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
        canonical: bool = False,
    ):
        """Publish an object package; ``canonical=True`` preserves the native model.

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
            package = self.export(
                Path(tmpdir) / "dataset_package",
                format=format,
                canonical=canonical,
            )
            return package.publish(
                dataset_key=dataset_key,
                uploader=resolved_uploader,
                idempotency_key=idempotency_key,
            )

    def info(self, print: bool = True, presentation: bool = True) -> str | None:
        """Print or return a human-readable summary of the model.

        Subclasses supply report sections; compact repr/str stay independent.
        Dataset-produced objects include their metadata, presentation and
        provenance unless ``presentation=False``. ``print=False`` returns plain
        text without writing to stdout or depending on the logging level.
        """
        from ..common._display import format_info

        summary = format_info(type(self).__name__, self._info_sections(), max_rows=20)
        if presentation and self.dataset_context is not None:
            from dtcc_core.datasets.presentation import format_dataset_context

            summary = (
                f"{summary}\n\n{format_dataset_context(self.dataset_context, obj=self)}"
            )
        if print:
            builtins.print(summary)
            return None
        return summary

    def print_info(self, file=None) -> None:
        """Print the human-readable model summary returned by ``info()``."""
        builtins.print(self.info(print=False), file=file)

    def save(self, *args, **kwargs):
        """Save the model to disk.

        IO methods are imported lazily so dtcc-core does not require them
        at import time. Importing ``dtcc_core.io`` registers object-specific
        ``save`` methods on DTCC model classes; this method then delegates
        to the registered implementation.
        """
        save_method = getattr(type(self), "save", None)
        if save_method is not None and save_method is not Model.save:
            return save_method(self, *args, **kwargs)

        try:
            importlib.import_module("dtcc_core.io")
        except Exception as exc:
            raise AttributeError(
                f"Cannot save object: {self.__class__.__name__}. "
                f"Failed to load dtcc_core.io module ({exc})."
            ) from exc

        save_method = getattr(type(self), "save", None)
        if save_method is not None and save_method is not Model.save:
            return save_method(self, *args, **kwargs)

        raise AttributeError(
            f"Cannot save object: {self.__class__.__name__}. "
            "No IO save method is registered for this model type."
        )

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
            warning(f"{function} Method {function_name} already exists, replacing it.")
            cls._methods.pop(idx)
            break
    cls._methods.append((name, function.__module__, function.__doc__))
    setattr(cls, name, function)
