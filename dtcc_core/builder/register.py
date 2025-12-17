from inspect import signature
from ..model.model import Model


def register_model_method(fn):
    """
    Register a function as a DTCC model method.

    This decorator inspects the function's signature and attaches it as a
    dynamically registered method to the corresponding ``Model`` subclass.
    The function must take a DTCC ``Model`` (or subclass) instance as its
    first argument, and this parameter must be annotated with the appropriate
    type hint. The decorated function is then added to the model's method
    registry via ``Model.add_methods``.

    Parameters
    ----------
    fn : callable
        The function to register. The function must have at least one
        parameter, and its first parameter must be annotated with a type
        derived from ``Model``. The function will be attached to that
        model class under its own name.

    Returns
    -------
    callable
        The original function, marked with the ``_dtcc_model_method`` attribute
        and registered as a method on the corresponding ``Model`` class.

    Raises
    ------
    ValueError
        If the function has no parameters.
    ValueError
        If the first parameter is not annotated with a subtype of ``Model``.
        This typically indicates a missing or incorrect type hint.

    Notes
    -----
    - Registered methods are added dynamically to the model class using
      ``Model.add_methods``.
    - The decorator marks the function with ``_dtcc_model_method = True`` to
      indicate that it was registered through this mechanism.
    - This mechanism enables an extensible plugin-style system for augmenting
      DTCC model classes with user-defined operations.
    """
    params = signature(fn).parameters
    if len(params) == 0:
        raise ValueError("Model method must have at least one parameter")
    first_arg = list(params.keys())[0]
    first_type = params[first_arg].annotation
    if not issubclass(first_type, Model):
        raise ValueError(
            f"First parameter must be a DTCC Model. Did you forget to add a type hint?"
        )
    fn._dtcc_model_method = True
    first_type.add_methods(fn, fn.__name__)
    return fn
