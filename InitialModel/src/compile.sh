# sediment thickness
ifort -c read-sed.f90

# water depth
ifort -c read-watdep.f90

# JdF seafloor age from Wilson
ifort -132 -c JdF-seafloor-age.f

# Global seafloor age in research area
ifort -c readglobalcrustage.f90

# read slab depth
ifort -c read-slab.f90

# read starting model of Western U.S. 
ifort -c read-wUS.f90

# Generate initial Vs model from geotherm amd mineralogical parameters
ifort -c thermalstart.f90

# generate local models
ifort -132 -heap-arrays -O2 -g -traceback  -o  gencrustJdF gencrustJdF.f90 read-sed.o read-watdep.o JdF-seafloor-age.o readglobalcrustage.o read-slab.o read-wUS.o thermalstart.o

# testing code to check the thermal models
ifort -132 -heap-arrays -O2 -g -traceback -o testtherm  test-thermalcode.f90 thermalstart.o
# move to dir
#mv gencrustJdF ../
