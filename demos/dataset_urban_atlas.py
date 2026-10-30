"""Plot Copernicus Urban Atlas land use classes for central Gothenburg."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(318700, 6397400, 320700, 6399400)

landuse = dtcc.datasets.urban_atlas(bounds=bounds)
landuse.info()
print("Urban Atlas classes:", sorted(set(landuse.attributes["ua_class"])))
landuse.plot()
