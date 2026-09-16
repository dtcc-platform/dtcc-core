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
