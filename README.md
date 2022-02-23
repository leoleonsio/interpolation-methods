# interpolation-methods
Implementation of spatial interpolations methods: 
- nearest neighbor
- inverse distance weighting (IDW)
- linear TIN
- laplace

The implementation is based on the knowledge from Terrain Book authored by TU Delft professors. https://github.com/tudelft3d/terrainbook

The code runs based on the parameters in _params.json_ file, which specifies:
- input file
  - path to the space delimited text file containing point coordinates with the header line being: **x** **y** **z**
- For each interpolation method:
  - output raster file (.asc format)
  - output raster cellsize
  - interpolation specific parameters

For more details look at example _params.json_
