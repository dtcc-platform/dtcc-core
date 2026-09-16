"""Offline examples of compact display and detailed inspection.

Run from the repository: .venv/bin/python sandbox/model_display.py
"""

import numpy as np

from dtcc_core import datasets
from dtcc_core.model import Bounds, Building, City, Field, Mesh, Point


point = Point(x=1.0, y=2.0, z=3.0)
bounds = Bounds(xmin=0.0, ymin=0.0, xmax=10.0, ymax=20.0)
temperature = Field(
    name="temperature", unit="°C", description="Temperature at each vertex",
    association="vertex", values=np.array([20.0, 21.0, 22.0]),
)
mesh = Mesh(
    vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    faces=np.array([[0, 1, 2]]), fields=[temperature],
)
city = City(id="example-city")
city.add_children([Building(id="building-1"), Building(id="building-2")])

print("Compact representations (repr, str and print agree)")
for obj in (point, bounds, mesh, city, temperature, datasets.smoke):
    assert repr(obj) == str(obj)
    print(obj)

print("\nDetailed field report")
temperature.info()
print("\nDetailed mesh report")
mesh.info()
print("\nDetailed city report")
city.info()
print("\nDataset parameter help")
datasets.smoke.info()

print("\nA small value with omitted state is also a summary")
point.fields.append(Field(name="temperature", unit="°C", values=np.array([20.0])))
print(point)

# Dataset context is included by default and can be hidden explicitly.
context = datasets.smoke.create_context(
    datasets.smoke.validate({"bounds": (0.0, 0.0, 10.0, 20.0)})
)
mesh.dataset_context = context
assert "Provenance" in mesh.info(print=False)
assert "Provenance" not in mesh.info(print=False, presentation=False)
