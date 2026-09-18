"""Docstring checks for dtcc_core.

Every public class and free function is documented, docstrings only document
parameters their function accepts, and all docstrings use NumPy style. The
checks read the source with ``ast`` and never import dtcc_core.
"""

import ast
import re
import textwrap
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE = REPO_ROOT / "dtcc_core"

# Generated protobuf bindings are not written by hand.
GENERATED_FILES = {"dtcc_pb2.py"}

NON_NUMPY_SECTION = re.compile(
    r"^(Args|Arguments|Parameters|Returns|Return|Raises|Examples?|Attributes):\s*$"
    r"|^:(param|return|returns|rtype|raises)\b",
    re.MULTILINE,
)


def _source_files():
    return [
        path
        for path in sorted(PACKAGE.rglob("*.py"))
        if path.name not in GENERATED_FILES
    ]


def _is_private_module(path):
    relative = path.relative_to(PACKAGE)
    return any(
        part.startswith("_") and part != "__init__.py" for part in relative.parts
    )


def undocumented_public_names(source):
    """Return ``(lineno, name)`` of public top-level definitions without docstrings."""
    return [
        (node.lineno, node.name)
        for node in ast.parse(source).body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
        and not node.name.startswith("_")
        and not ast.get_docstring(node)
    ]


def documented_parameters(docstring):
    """Return the names listed in the NumPy parameter sections of a docstring."""
    names = []
    lines = docstring.splitlines()
    for index, line in enumerate(lines[:-1]):
        if line.strip() not in ("Parameters", "Other Parameters"):
            continue
        if set(lines[index + 1].strip()) != {"-"}:
            continue
        entries = [entry for entry in lines[index + 2:] if entry.strip()]
        if not entries:
            continue
        # Take the indentation from the entries: a misindented header must not
        # hide a parameter list from this check.
        base = len(entries[0]) - len(entries[0].lstrip())
        for entry in lines[index + 2:]:
            if not entry.strip():
                continue
            indent = len(entry) - len(entry.lstrip())
            if indent < base:
                break
            if indent > base:
                continue
            if re.fullmatch(r"[A-Z][A-Za-z ]*", entry.strip()):
                break  # The next section heading.
            head = entry.strip().split(":")[0].strip("` ")
            for part in head.split(","):
                part = part.strip().lstrip("*")
                if re.fullmatch(r"[A-Za-z_]\w*", part):
                    names.append(part)
    return names


def unknown_documented_parameters(source):
    """Return ``(lineno, name, unknown)`` for functions documenting parameters they lack.

    Functions accepting ``**kwargs`` may document extra keyword names and are
    skipped.
    """
    problems = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        docstring = ast.get_docstring(node)
        if not docstring or node.args.kwarg is not None:
            continue
        args = node.args
        accepted = {arg.arg for arg in args.posonlyargs + args.args + args.kwonlyargs}
        if args.vararg is not None:
            accepted.add(args.vararg.arg)
        accepted.update(("self", "cls"))
        unknown = [name for name in documented_parameters(docstring) if name not in accepted]
        if unknown:
            problems.append((node.lineno, node.name, unknown))
    return problems


def non_numpy_docstrings(source):
    """Return ``(lineno, name)`` of docstrings using Google or reST field sections."""
    tree = ast.parse(source)
    problems = []
    for node in ast.walk(tree):
        if not isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            continue
        docstring = ast.get_docstring(node)
        if docstring and NON_NUMPY_SECTION.search(docstring):
            problems.append((getattr(node, "lineno", 1), getattr(node, "name", "<module>")))
    return problems


def _report(check, paths):
    failures = []
    for path in paths:
        for lineno, *details in check(path.read_text(encoding="utf-8")):
            where = path.relative_to(REPO_ROOT)
            failures.append(f"{where}:{lineno} " + " ".join(str(d) for d in details))
    return failures


def test_public_classes_and_functions_have_docstrings():
    public_modules = [path for path in _source_files() if not _is_private_module(path)]
    failures = _report(undocumented_public_names, public_modules)
    assert not failures, (
        "Public classes and functions need a NumPy-style docstring:\n"
        + "\n".join(failures)
    )


def test_docstrings_only_document_accepted_parameters():
    failures = _report(unknown_documented_parameters, _source_files())
    assert not failures, (
        "Docstrings document parameters the function does not accept:\n"
        + "\n".join(failures)
    )


def test_docstrings_use_numpy_style():
    failures = _report(non_numpy_docstrings, _source_files())
    assert not failures, (
        "Docstrings must use NumPy style (Parameters / Returns with underlines):\n"
        + "\n".join(failures)
    )


def test_checks_detect_problems_in_sample_source():
    sample = textwrap.dedent(
        '''
        def documented(a):
            """Do something.

            Parameters
            ----------
            a : int
                A real parameter.
            b : int
                Not a parameter.
            """

        def undocumented():
            pass

        def _private():
            pass

        def keywords(**kwargs):
            """Accept keywords.

            Parameters
            ----------
            anything : int
                Documented keyword.
            """

        class GoogleStyle:
            """Old style.

            Args:
                x: A value.
            """
        '''
    )

    assert [name for _, name in undocumented_public_names(sample)] == ["undocumented"]
    assert unknown_documented_parameters(sample) == [(2, "documented", ["b"])]
    assert [name for _, name in non_numpy_docstrings(sample)] == ["GoogleStyle"]



def test_check_reads_parameters_under_a_misindented_header():
    """A header indented differently from its entries must still be checked."""
    sample = textwrap.dedent(
        '''
        def misindented(a):
            """Do something.

             Parameters
            ----------
            a : int
                A real parameter.
            ghost : int
                Not a parameter.
            """
        '''
    )

    assert unknown_documented_parameters(sample) == [(2, "misindented", ["ghost"])]

def test_private_modules_are_exempt_from_the_docstring_requirement():
    assert _is_private_module(PACKAGE / "model" / "_display.py")
    assert _is_private_module(PACKAGE / "_private" / "module.py")
    assert not _is_private_module(PACKAGE / "model" / "__init__.py")
    assert not _is_private_module(PACKAGE / "model" / "geometry" / "mesh.py")
