## Scripts
### `generatePointsPolygon.py`
Python script that samples points inside a random polygon. Dependencies in `environment.yaml`. Takes as command line arguments the width `W`, height `H`, number points on the polygon `npoints`, number of points samples inside the polygon `nsample`. Example:
```bash
python generatePointsPolygon.py 1400 1400 10 100
```

### `genplot.jl`
Contains functions to run the simulation, and generate time series plots.

### `main.cpp`
Runs the main simulation. 
```bash
./main $ncells $ntimesteps $thresholdsize
```

### `pfbound_testboundary.cpp`
Script to use and test different boundary detection methods.

### Header files
1. `boundaryFunctions.h`
2. `CELL.h`
3. `constants.h`
4. `utilityFunctions.h`
5. `VEC2.h`
6. `vecCellUtilities.h`

## Folders
`intermediateResults`: stores data files that are linked between languages. 
`results`: stores the final simulation plots
`modulesCpp`: stores all the c++ modules.  

## Files
`environment.yaml`: Python dependency files.
`Manifest.toml`, `Project.toml`: Julia dependency files.
