## Scripts
### `generatePointsPolygon.py`
Python script that samples points inside a random polygon. Dependencies in `environment.yaml`. Takes as command line arguments the width `W`, height `H`, number points on the polygon `npoints`, number of points samples inside the polygon `nsample`. Example:

```bash
python generatePointsPolygon.py 1400 1400 10 100
```

### `genplot.jl`



### Header files
`CELL.h`
`VEC2.h`
`constants.h`
`utilityFunctions.h`


## Folders
`intermediateResults`: stores data files that are linked between languages. 
`results`: stores the final simulation plots
`modulesCpp`: stores all the c++ modules.  

## Files
`environment.yaml`: stores dependencies for `generatePointsPolygon.py`.

