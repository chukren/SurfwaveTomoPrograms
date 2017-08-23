Usage of 2Dkernel
------------------
procedure for generating sensitivity kernels for each period.   
Originally written by Yingie Yang, 2004.

1. Run **createsac.f**
  * Note: Compile with $SACLIB/sacio.a
  * Output file name createfile020.sac

  * Look at file in sac
    decide on cut window near center of sine wave(Note: depends 
  on way data was windowed, usually ~200s for short periods, 
  may change for long periods)
     
  * Comment added by Yun Wang. How to choose cut1 and cut2 doesn't matter. 
  What really matters is the windowing length. It seems the length of 
  the first lobe depends on the windowing length. The wider of the 
  windowing the more narrow the lobe. It is independent of frequency.  

2. Run **bpcut2.m**   
  * cut1, cut2 = narrow cut
  * input filename createfile020.sac
  * ctb1,ctb2  = wide cut (just data limit)
  * output filename **createfile020.cut700t900**	

3. Look at spectrum, see where peak is, and where first lobe starts and ends
   ```bash
   fft createfile020.cut700t900    
   wsp   
   r createfile020.cut700t900.am
   ```
4. Run **getdatafromsac.f**   
   * compile with **$SACLIB/sacio.a**
   * input filename **createfile020.cut700t900.am**
   * output filename **createfile020.cut700t900.am.dat**	

   look at **f20.cut800t1000.am.dat** file   
   find highest value at spectrum peak   
   find smallest value before and after this   
   remove lines above and below this range   
   add first line with number of frequencies retained in first lobe	

5. Run **sensitivity.f**   
   (normal compiling)   
   input phase velocity   
   input spectral file: **f20.cut800t1000.am.dat**   
   > (input number of freq: #lines in modified f20.cut800t1000.am.dat)   
   > (input unsmoothed file name temp.20.sac)   
   input smoothed file name  sens20s50km.dat   
   input scalelen (50 km for GLIMPSE, ~grid spacing - this is smoothing length, not grid size)   

   * (Comment added by Yun Wang. The operation below could be done in the 
   simannerr program. That is to say, don't need to edit files like 
   **sens20s50km.dat**.)
   * (Comment by Don Forsyth. The operation below is now done in sensitivity.f)


   > edit sens20s50km.dat before running with simannerr inversions   
   > add these 2 lines at top:   
   > 301. -1500. 10.   
   > 301. -1500. 10.   

Done!
