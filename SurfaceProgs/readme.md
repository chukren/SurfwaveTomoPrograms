Note about inversion strategy by Don
------------------------------------
1. I recommend an iterative approach to a finer resolution grid -- the lateral 
resolution at longer periods is not as good as at shorter periods, which effectively 
means that longer periods are more heavily damped when trying to invert on a 
fine grid.  Using the same crustal grid the whole time, I start out with a 1 or 
1.2 degree node grid initially for the phase velocity inversions.  

2. For the shear velocity inversions, I check the ranks (for velocity parameters) in 
the phase velocity inversions at different periods and then use only those periods 
that have ranks that are > 80% (arbitrary number) of the rank of the period with 
the highest rank.  I do the shear velocity inversions for 3-D model, then use the 
predicted phase velocities for that 3-D model (provided as output from the shear 
inversion routine) as the starting model for a new phase velocity inversion on a 
finer grid.  Again, I use the output from that inversion to invert for new shear 
velocity model, but again using only those periods with ranks > 80% of maximum.   
Typically on the course grid, all periods are used, but as the grid gets finer, 
only the shorter periods are retained.

