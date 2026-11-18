"""Plot named OpenStreetMap places by category, then only food and drink."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

places = dtcc.datasets.places(bounds=bounds)
places.info()

table = places.to_geodataframe()
print(table[["name", "category", "osm_key", "osm_value"]].head(10))

places.plot()

food = dtcc.datasets.places(bounds=bounds, categories=["food_and_drink"])
food.plot(labels=20, title="Food and drink")
