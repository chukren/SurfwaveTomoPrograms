Usage of phaseVelInversion:
---------------------------
1. There are numerous changes from prior simulated annealing inverse routines.  
  * First, the grid is created within the inversion routine.  
  * Second, the sensitivity kernels are also generated within the routine, but 
  the spectral files describing the periods that are effectively averaged in 
  calculating the amplitude and phase given the filters we use and the typical 
  window lengths are read in (provided as files for each period - but I haven't 
  provided the sac routines for generating the spectral files).   
  * Third, the interpolation scheme is now linear instead of gaussian, so the 
  sensitivity does not extend past the nearest node points and the velocities 
  at the node points are actual velocities in the velocity model instead of 
  coefficients that are then used to calculate what the velocities are.  
  * Fourth, there also is a strong minimum curvature constraint built in that 
  tries to minimize the covariance between velocity at a given grid point and 
  its neighboring grid points. 
  * Fifth, the search for the best wave parameters is now conducted by a grid 
  search method instead of by simulating annealing.  This approach seems to 
  give more stable results.  
  * Finally, the program outputs not only the solution at the node points, but 
  also the interpolated velocities at the typically finer grid points where 
  crustal structure had been specified on a lat/lon grid.  (There are probably 
  some other changes as well, but not as important as these.)

2. Four versions of the phase velocity inversion are provided.
  * __srchwave587__:  2-D grid of points for isotropic component of velocity plus 
  uniform 2theta azimuthal anisotropy terms
  * __srchwave589__:  2-D grid of points for isotropic component and same 2-D for 
  lateral variations in anisotropy
  * __srchwave590__:  Allows different resolution grids for isotropic and anisotropic 
  components.  __NOTE__ that there are substantial differences between this routine and 
  the others, including more files that have to be read in (specified in the general 
  input file)
  * __srchwave591__:  uses a uniform velocity starting model and no station phase 
  corrections - otherwise like 587

3. one example general input file is given, see **srchscript** for file name.  This 
file has all the file names for the individual seismograms in it as well as all the 
various input and output files with the grid info, prior predicted phase velocities, 
etc.).  NOTE that each of the inverse programs and **rdsetupsimul** (which needs to 
be run before any of the **srchwave** to do initial fourier analysis) need to be 
changed so that the programs can extract the station names from the seismogram file 
names.

4. After running the **srchwave** programs, use **combineGridC** and **combineGridAniso** 
routines to combine the output from all the periods to be used in the shear velocity 
inversion routines.  Example input files are provided.  Spectral files are provided for 
constructing the sensitivity kernels and there are some GMT mapping scripts included.
