Usage of shearinv:
------------------
1. Routines are provided for inverting shear velocity and anisotropic 
structure as a function of depth at all the points where crustal structure 
had been given.  

2. There are a bunch of adjustable parameters - I recommend doubling the 
standard deviations that come out of the phase velocity inversions because 
they are typically underestimated due to the fact that the wave parameters 
are fairly strongly damped in the linearized part of the inversion.  Thus, 
the full tradeoff between wave parameters and velocity structure is not 
represented in the uncertainties.  Other parameters specify in what proportion 
P-wave velocity and density will vary when S-wave velocity varies.  See 
program for more details.  

3. If you change the number of periods you are inverting, you will have to 
specify that number in a couple of places and also modify the **perstart2.d** 
file accordingly.  The inverse routine currently does not allow changes in 
crustal thickness.  There is a strong minimum curvature constraint applied, 
but this does not apply across the Moho.
