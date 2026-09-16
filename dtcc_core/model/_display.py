"""Small shared report sections for native model data."""

from itertools import islice


def field_section(fields):
    return ("Fields", ("Name", "Unit", "Association", "Shape", "Type", "Description"),
            [(f.name, f.unit, f.association, f.values.shape, f.values.dtype, f.description)
             for f in fields])


def sample_section(title, models):
    rows = [(i, repr(model)) for i, model in enumerate(islice(models, 20))]
    if len(models) > 20:
        title = f"{title} (20 shown, {len(models) - 20} omitted)"
    return (title, ("Index", "Summary"), rows)


def print_tree(model, label, *, indent="", verbose=False, max_depth=None):
    """Stream a hierarchy without building a second in-memory model tree."""
    from .object.object import Object
    from .geometry.geometry import Geometry

    def entries(node):
        if isinstance(node, Object):
            if verbose:
                for key, value in node.attributes.items():
                    yield f"attribute {key!r}: {value!r}", None
            for name, record in node.geometry.items():
                metadata = ""
                if verbose:
                    descriptors = [f"{key}={value!r}" for key in ("lod", "role")
                                   if (value := getattr(record, key)) is not None]
                    if descriptors:
                        metadata = f" ({', '.join(descriptors)})"
                yield f"{name!r}{metadata}: {record.geometry!r}", record.geometry
            for group in node.children.values():
                for child in group:
                    yield repr(child), child
        elif isinstance(node, Geometry) and verbose:
            for field in node.fields:
                yield f"field {field.name!r} ({field.unit}): {field.description}", None

    def walk(node, text, prefix, branch, depth):
        children = iter(entries(node))
        current = next(children, None)
        limited = max_depth is not None and depth >= max_depth and current is not None
        # Keep labels on one line so attribute strings cannot break the branches.
        text = text.replace("\r", "\\r").replace("\n", "\\n")
        print(f"{prefix}{branch}{text}{' ...' if limited else ''}")
        if limited:
            return
        child_prefix = prefix + ("    " if branch == "└── " else "│   " if branch else "")
        while current is not None:
            following = next(children, None)
            text, child = current
            walk(child, text, child_prefix, "└── " if following is None else "├── ", depth + 1)
            current = following

    walk(model, label, indent, "", 0)
