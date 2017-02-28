Generate initial model (oceanic, forearc, and continental) for Cascadia Rayleigh wave tomography
------------------------------------------------------------------------------------------------
### Usage:  
   $> src/gencrustJdF < PAR_JDFCRUST  

### Input: 
  * PAR_JDFCRUST  
   beglat, endlat, latinc, beglon, endlon, loninc (36.0 54.0 0.2 -136.0 -118.0 0.2)  
   coastfile (JdF-NA-large.dat  modified so the water depth, sediments and age data all available in this region)  
   outputfile (3DcrustJdF_0.2x0.2)  
   mantlefile (mantlemodel.dat)  
   linshearfile (LinShearer.dat obsolete, not used)  
   mohofile (cas.moho.dat)  
   defaultfile (defaultcrustCascade.dat)  
   rayl8secfile (moschetti_etal_Rayleigh_c_8s.dat)  

### Output:  

  * Output model:  
   outputfile

  * Outputs for checking:  
   tempoutput  
   testbasement.dat  
   testmoho.dat  
   testregion.dat  
   testsedvs.dat  
   testslabtopo.dat  
   testwater.dat

### Notes
  * File dependence  
   defaultcrust2.dat -> defaultcrustCascade.dat  
   Crust_thick2 -> cas.moho.dat  

  * Files needed for subroutine calls:  
   readtopobin()   <--  JdFtopo.bin, JdFtopo.xyz  
   read_sed_xyz()  <--  JdFsed.xyz  
   readglobalcrustage() <-- global.age.JdF.xy  
   read_slab_xyz()	<--  cas_slab1.0_clip.xyz  
   read_wUS_xyz()  <--  wUS-SH-2010.xyz  
