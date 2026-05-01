from dtcc_core.model import Point


def test_model_info_defaults_to_str():
    point = Point(x=1.0, y=2.0, z=3.0)

    assert point.info(print=False) == str(point)


def test_model_info_prints_by_default(capsys):
    point = Point(x=1.0, y=2.0, z=3.0)

    result = point.info()

    captured = capsys.readouterr()
    assert result is None
    assert str(point) in captured.out


def test_model_print_info(capsys):
    point = Point(x=1.0, y=2.0, z=3.0)

    point.print_info()

    captured = capsys.readouterr()
    assert str(point) in captured.out
