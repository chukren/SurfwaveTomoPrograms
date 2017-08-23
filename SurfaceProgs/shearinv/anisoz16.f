c  anisoz16.f inverts for vertical distribution of azimuthal anisotropy
c  solves separately for cos2theta and sin2theta termsx.  Outputs complete covariance matrix 
c  for each inversion and optionally inputs covariance matrix from previous
c  inversion to be used as the a priori model covariance matrix 
c
c  anizoz16 < anisozinp
c
c  Doesn't do a true anisotropic inversion - just inverts as if 2theta terms were perturbations
c  in phase velocity and solving for perturbations in shear velocity that fit 
c
c  WARNING  WARNING WARNING
c  Presumes aniso and structure input are on same grid
c

c  shearz16.f optionally finds optimum smoothing parameter to minimize 
c  off-diagonal elements of a posteriori covariance matrix


c This inversion uses phase velocity data for each period
c from surface wave inversion output (e.g. simannerr1.a.f)

c Partial derivatives are simultaneously calculated(output) by s_saito.f
c for dc/dbeta, dc/dalpha, dc/drho and used in this inversion
c for the partial derivative matrix, G,  (just read in).                  

c Residuals for this inversion are obtained by calculating 
c a third set of phase velocities by running the 
c new shear wave velocities and their corresponding P velocities 
c (a layered model) through s_saito.f and comparing the
c difference in these phase velocities from the data. 

c Fix alpha or 
c calculate the new layered model (for residuals) by:
c dalpha=dbeta*sqrt(3) and
c alphanew = alphastart + dalpha 
c (betanew is calculated in this inversion similarly)

c Sequential iterations must recalculate a new starting model
c (can use the one used for residuals) from shear wave velocity 
c output of this inversion and use it's phase velocities
c as model input. Then invert here for another shear wave profile.
c Use the same phase velocity data (output from simannerr1.a.f) 
c (Again, new residuals are calculated as above)
c Also use new partial deriv's from s_saito.f

c Program has now been modified to do any number of iterations 
c desired. Still must copy desired starting model to scalifornia_refer.dat
c before running this. The file copied to scalifornia_refer.dat must be
c the same as startmod and startmod2 in input files.

c Forward model: ph.vel. predictions for an input layered model 
c (e.g. ak135 or PREM) obtained from s_saito.f 



 
      parameter (npermax=30,maxlay=100, ngridmax = 20000)  
*********************************************************************                                
      
      real*4    freq(npermax),damp,depth(maxlay)
      real*4    cmminvmmo(maxlay),tempcov,smooth
      real*4    depmin(maxlay),depmax(maxlay)
      real*4	zero,density
      real*4	xxth,xxdens,xxlong,xxdep
      real*4    depmant(100),densmant(100),almant(100),betmant(100)
      real*4    densrat(100)	
      real*4    psait(npermax),phasait(npermax),bandsait(npermax)
      real*4    depsait(npermax),xxd5,betsait(npermax)
      real*4    lat,lon,lontemp, lattemp 
c      real*4    veltemp(npermax), stdtemp(npermax)	
      real*4    Rough(maxlay,maxlay)
c      real*4    oscil2(maxlay,maxlay),oscil(maxlay,maxlay)
      real*4    penalty(maxlay,maxlay)
      real*4    thick(maxlay),dens(maxlay)
      real*4 stddvcs(npermax,ngridmax),sn2th(npermax,ngridmax)
      real*4 stddvsn(npermax,ngridmax),cs2th(npermax,ngridmax)
      real*4 velsum1(npermax),velsum2(npermax)
      real*4 sumstdev1(npermax),sumstdev2(npermax)
      real*4 anisolon(ngridmax), anisolat(ngridmax), wgt(ngridmax)
      real*4 stcs2th(ngridmax),stsn2th(ngridmax)
      real*4 lontemp2,lattemp2


      double precision xxstalpha
      double precision stalpha(maxlay),stbeta(maxlay)
      double precision crrntmodb(maxlay),crrntmoda(maxlay)
      double precision crrntmodr(maxlay)
      double precision ccobs(npermax),stbetak135(maxlay)
      double precision origmod(maxlay),eachrank(maxlay)
      double precision ccobstd(npermax),stddevdata(npermax),ddd
      double precision ccpred(npermax),change(maxlay),correct(maxlay)
      double precision ccpredout(npermax),residout(npermax)
      double precision ccpredout2(npermax)
      double precision perd(npermax),xxd1,xxd3,xxp,d(npermax)
      double precision g(npermax,maxlay),gtd(maxlay),gtdcmm(maxlay)
      double precision gtg(maxlay,maxlay),savegtg(maxlay,maxlay)
      double precision gtginv(maxlay,maxlay),stddev(maxlay)
      double precision covinv(maxlay,maxlay),cmm(maxlay,maxlay),ccc
      double precision dataimp(npermax),tempimp(maxlay),resid(npermax)
      double precision tempgtginvgt(maxlay,npermax),xxanom
      double precision fcov(maxlay,maxlay)
      double precision sumsq,sigma2,anom(maxlay)	
      double precision dca(maxlay,npermax),dcb(maxlay,npermax)
      double precision dcr(maxlay,npermax)
      double precision resmatx(maxlay,maxlay),tempmatx,resparm,rank
      
      integer*4 iter,icnt,irow,icol,indx(npermax)
      integer*4 nparam,xxd2,xxd4
      integer iwgt(ngridmax)
      
      character*70 xxd0,newmodlout,flowmantle,covmatin,covmatout
      character*70 startmod,startmod2,phdata,shvelout,shvelmatout
      character*70 stmodnext,anomalyout,absvelout,resmatrix
      character*70 label(2), anisodata, anisoout(npermax)
      character*70 prioraniso,finalaniso
      
      read(*,*) nper 
c  nper is number of periods. 
 
      read(*,*) stdmult, damp, smthmu, p2s,nit, rmvnegfac,rho2s, 
     1    inpcov, invprd,indprior
      
c  stdmult is term to multiply standard deviations by - formal estimates are often too 
c  small.. Multiplier should be large enough that period to period fluctuations or
c  lateral fluctuations are of the order of the assigned standard deviations
c
c  inpcov indicates whether prior covariance matrices are to be used - if gt
c  zero, will expect covariance input from previous run.  If > zero, then
c  damp and smthmu will have no effect as damping will be performed by 
c  a priori matrix
c
c  invprd indicates whether to invert predictions in second iteration (invprd = 1) to obtain
c  smoother model or to just stop at straight data inversion
c
c  indprior eg. 1 indicates to read in prior anisotropy model and use it as 
c  starting model for anisotropy inversion
c
c  damp is a priori model std deviation in km/s to use as damping parameter 
c  for shear velocities
c
c  smthmu is coefficient for smoothing - minimum curvature Use either 0.0 for no application of
c  minimum curvature or 1.0 for optimum damping
c
c  p2s is ratio to assume for p velocity changes compared to s changes.  0 keeps p velocity fixed
c  1.73 makes changes equivalent to Poisson's ratio of 0.25.  1.0 assumes absolute P changes are 
c  equal to S changes, which means that S changes fractionally more (such as for sediments or melt
c  effects).
c
c  nit = number of iterations

c  in this version, input nit means nothing.  It iterates once to correct for negative lobes
c  of resolution matrix.

c  rmvnegfac  -  multiplicative factor to control multiples of negative wings 
c  predicted from resolution matrix to remove from starting model
c
c  rho2s is ratio to assume for fractional density changes compared to fractional s changes.  

c        nit = 2

	read(*,*) nper2, depth1, error1,ijk1,ijk2,ijk3,topd,bottomd
	
c  These are variables controlling the search for phase velocity (roots to equations) for
c  given structure.  nper2= # periods, depth1 default depth extent, error1, criterion for 
c  ending search, ijk1 = 4 is Rayleigh spherical earth, 3 is Love spherical earth, ijk2 = 0, ijk3 = 2
c  topd and bottomd control depth range in which partial derivatives are output  

        nobs = nper

        read(*,*) startmod
	read(*,*) flowmantle
c	read(*,*) phdata
        read(*,*) anisodata
	read(*,*) newmodlout
	read(*,*) shvelout
        read(*,*) resmatrix
	read(*,*) covmatin
	read(*,*) covmatout
	read(*,*) prioraniso
c  startmod      3D crustal/mantle model down to maximum inversion depth
c  flowmantle    mantle below maximum inversion depth 
c  anisodata     output from phase velocity inversion summarized for all periods
c  newmodlout    output of this program - new 3D anisotropy cos2theta,sin2theta model
c  shvelout      details of inversions
c  resmatrix     resolution matrices
c  covmatin      apriori covariance matrix - often not used
c  covmatout     covariance matrices
c  prioraniso    3D anisotropy model to be used as starting model	
	if (indprior.eq.1) then
	   open(55, file = prioraniso)
	   read(55,*) beglat,endlat,dlat,beglon,endlon,dlon
	endif
	   
	
	open(10, file=startmod)
	open(90, file = newmodlout)
	
	read(10,*) beglat,endlat,dlat,beglon,endlon,dlon
	write(90,*) beglat,endlat,dlat,beglon,endlon,dlon

	nlat = (endlat-beglat)/dlat +1.01
        nlon = (endlon-beglon)/dlon + 1.01
        nxy = nlat*nlon

c  read in aniso observations of cos2theta and sin2theta - spacdeg is distance between
c  node points in degrees  This version presumes aniso observations on same grid as
c  shear velocity and 3D shear velocity model so naniso and spacdeg parts are
c  commented out	
	open(12, file = anisodata)	
c	read(12,*)  spacdeg
	read(12,*) nper
	naniso = nxy
c	do ianiso = 1,naniso
	do ianiso = 1,nxy
	  read(12,*) anisolon(ianiso),anisolat(ianiso)
	  do j = 1, nper
	    read(12,*) cs2th(j,ianiso),stddvcs(j,ianiso),
     1          sn2th(j,ianiso), stddvsn(j,ianiso)
            stddvcs(j,ianiso) = stdmult*stddvcs(j,ianiso)
	    stddvsn(j,ianiso) = stdmult*stddvsn(j,ianiso)
          enddo
	enddo 
	
	label(1) = ' cos2theta '
	label(2) = ' sin2theta '
	ione = 1
	do i = 1,nper
	  read(*,*) anisoout(i)
	  open(60+i, file = anisoout(i))
	  write(60+i,*) beglat,endlat,dlat,beglon,endlon,dlon
	  write(60+i,*) ione
	enddo
			
  	open(18, file = shvelout)
        open(40, file = flowmantle)
	
	open(50, file = covmatout)
	if (inpcov.gt.0) open(51, file = covmatin)

	
	read(40,*) nlowlay
        do i=1,nlowlay                                                               
	   read(40,*) depmant(i),densmant(i),almant(i),betmant(i)               
	enddo
	
	close(40)

c  read in parameters that control interval in which Saito routines search for
c  phase velocity solutions
        open(41, file = 'perstart2.d')
	do i=1,nper
	  read(41,*) psait(i),phasait(i),bandsait(i),betsait(i),depsait(i) 
	enddo
	close(41)

        open(52,  file = resmatrix)
	
	sumnormmsft2 = 0.0
c==========================================================================================================         
c  read in starting shear wave velocity model 
c       model to bottomd as starting model 
c unit=10: startmod.dat 
	do inxy = 1, nxy
          read(10,*) lontemp, lattemp
          read(10,*) nlay
	  lon = lontemp
	  lat = lattemp
           write(*,*) inxy, lon,lat
	   write(90,*) lontemp,lattemp
           write(90,*) nlay

c  Interpolate aniso grid solution onto finer a priori grid for passing along to next phase velocity inversion
c  or for finer plotting	
c	  wgtsum = 0.0
c	  do j = 1, nper
c	    velsum1(j) = 0.0
c	    velsum2(j) = 0.0
c	    sumstdev1(j) = 0.0
c	    sumstdev2(j) = 0.0
c	  enddo
c	  nwgt = 0
c	  do kk = 1,naniso
c            sclfac = cos(convdeg*anisolat(kk))
c	    xsep = abs(lon-anisolon(kk))*sclfac
c	    if ((xsep.le.spacdeg).and.
c     1          (abs(lat-anisolat(kk)).le.spacdeg))  then
c              wgt(kk) = (1.- abs(lat-anisolat(kk))/spacdeg)
c     1                       *(1.- xsep/spacdeg)
c              nwgt = nwgt+1
c              iwgt(nwgt) = kk
c              wgtsum = wgtsum + wgt(kk)
c	      do j = 1,nper
c	        velsum1(j) = velsum1(j) + wgt(kk)*cs2th(j,kk)
c	        sumstdev1(j) = sumstdev1(j) + wgt(kk)*stddvcs(j,kk)
c	        velsum2(j) = velsum2(j) + wgt(kk)*sn2th(j,kk)
c	        sumstdev2(j) = sumstdev2(j) + wgt(kk)*stddvsn(j,kk)
c	      enddo
c	     endif
c	   enddo
c	  if (nwgt.gt.0) then
c      	    do j = 1,nper
c	      velsum1(j) = velsum1(j)/wgtsum
c	      sumstdev1(j) = sumstdev1(j)/wgtsum
c	      velsum2(j) = velsum2(j)/wgtsum
c	      sumstdev2(j) = sumstdev2(j)/wgtsum
c	    enddo
c	  else
c      	    do j = 1,nper
c	      velsum1(j) = 0.0
c	      sumstdev1(j) = damp
c	      velsum2(j) = 0.0
c	      sumstdev2(j) = damp
c	    enddo
c	  endif
C  replaces above commented out code - assumes velocity points same as aniso points
	  do j = 1, nper
	    velsum1(j) = cs2th(j,inxy)
	    velsum2(j) = sn2th(j,inxy)
	    sumstdev1(j) = stddvcs(j,inxy)
	    sumstdev2(j) = stddvsn(j,inxy)
	  enddo
          

c nlayall: 2X layers, including lower mantle parts
c np=nlay: number of layers in the lithosphere and upper mantle 
c nparam: number of layers which model needs to resolve
        np = nlay
        nparam = np  
	
	if (inpcov.gt.0) then
	  do i = 1,np
	    read(51,*) dummy
	    do j=1,i
	      read(51,*) ii,jj, covinv(i,j)
	      covinv(j,i) = covinv(i,j)
	    enddo
	  enddo
	endif
        
********************************Modify here**************************      
c 54 is twice the number of layers in the lower mantle  
c        nlayall= 2*np + 54
	nlayall = 2*nlay + nlowlay
*********************************************************************
         open(30, file = 'scalifornia_refer.dat')	

         zero = 0.
	 write(30,*)'    0    0    0'
         write(30,*) '   1'
         write(30,*)'df              col # = 47'
         write(30,*) nlayall, '    6371.0                    982.0                          
     r     0 0 1'           
         write(30,*)'(F6.2,3F8.4)'
         do ii=1,np
            read(10,*) thick(ii),dens(ii),stalpha(ii),stbeta(ii)
c  guard against thin water layers that are difficult to compute response for
	    if (thick(ii).le.0.4) thick(ii) = 0.401
            write(30,210) zero,dens(ii),stalpha(ii),stbeta(ii)
            write(30,210) thick(ii),dens(ii),stalpha(ii),stbeta(ii)
         enddo
       

c mantle starting model (below our study depth)	  
	do i=1,nlowlay
           write(30,210) depmant(i),densmant(i),almant(i),betmant(i)
	enddo

c210    format (1x,f5.2,2x,f6.3,2x,f7.3,1x,f6.3)
210    format (f6.2,3f8.4)
	write(30,*) '   4   2872.3  0' 
	write(30,*) '   1'  
		
   	write(30,310) nper2, depth1, error1,ijk1,ijk2,ijk3,topd,bottomd
310     format(I5,F7.1,E12.2,1x,3I2,F6.1,F12.1)
	write(30,*) ' (6F10.6,2I2)'
*************************************************************************
         
	do i=1,nper
  	  write(30,230) psait(i),phasait(i),bandsait(i),betsait(i),
     1       depsait(i)
	enddo
         
         close(30)
         call s_saito
                 
c  for anisotropy, original model = 0 unless indprior = 1 and we consider crude approximation that
c  sensitivity is just to SV            
c           
c set up apriori damping of model parameters                                            !
                                                                                        !
c covinv(i,i) for minimum length damping only (uses diagonals of Cmm only)              !
c for non zero's only on the diagonals, matrix inverse is easy:                         !
c inverse = 1/damp**2 on individual diagonals                                           !
        if (inpcov.eq.0) then                                                                                !
	do i=1,np 
	  do j = 1,np
	    covinv(i,j) = 0.0
	  enddo                                                                      !      
	  covinv(i,i)=1./(damp**2.)                                                         !
	enddo	                                                                        !
         
c--------------------------------      
c Rough matrix 
       do irow = 1, np
       do icol = 1, np
           Rough(irow,icol) = 0.0
        enddo
       enddo	
c---------------------------------          
c                                !------ this definition of smoothness 
        Rough(2,1) = -3          !------ matrix has minimum first order 
        Rough(2,2) = 6           !------ derivative constraint on edges
        Rough(2,3) = -4          !
        Rough(2,4) = 1           !
        Rough(np-1,np-3) = 1     !
        Rough(np-1,np-2) = -4    !
        Rough(np-1,np-1) = 6     !
        Rough(np-1,np) = -3      !
c--------------------------------

	   
       do irow = 3, np-2
          icol = irow
          Rough(irow, icol-2) = 1
          Rough(irow, icol-1) = -4
          Rough(irow, icol) = 6
          Rough(irow, icol+1) = -4
          Rough(irow, icol+2) = 1
       enddo

c------------------------------------------------------------! 
c fix the water and sediment layers (by heavy damping)       !
	   do i = 1,np                                       ! 
           if(stbeta(i).lt.2.0) then                         !
             Rough(i,i) = 1.0E+4  
	     covinv(i,i) = 1.0E+4                           !
           endif                                             !
	   enddo                                             !
c------------------------------------------------------------!	
         endif
                                                                                       !
c        endif      !-------------------------------------kk------------------------------!

       
c  initialize g matrix
        do irow = 1, nobs
        do icol = 1, np
          g(irow,icol) = 0.0                                                           !
        enddo                                                                          !
        enddo                                                                          !
       
c  data vector and partial derivatives listed event by event with all
c  real data for first event followed by imaginary data, then onto next event
c  d contains misfit to starting model

c data vector is difference of observed and predicted phase velocities
c del d in inversion equation
c use apriori stdev?

c***************************************************
c  repeat separately from here for cos2theta and sin2theta inversions
         if (indprior.eq.1) then
	   read(55,*) lontemp2,lattemp2
	   read(55,*) nlay2
	   do janso = 1, nlay2
	     read(55,*) thick(janso),stcs2th(janso),stsn2th(janso)
	   enddo
	 endif
	   
         indaniso = 1
2000	 open(16, file='DERIV.DATA')
	    
	 do i=1,nobs
           read(16,*) xxd0                           !!!!unit=16:  DERIV.DATA
           read(16,*) xxd1
           read(16,*) xxd2,xxd3,perd(i),ccpred(i)
           read(16,*) xxd4
	   
 	   do j=1,np
	     read(16,*) depth(j),dcr(j,i),dca(j,i),dcb(j,i)
	     if (indaniso.eq.1) then
               g(i,j)= dcb(j,i)*thick(j)/sumstdev1(i)
	     else
	       g(i,j)= dcb(j,i)*thick(j)/sumstdev2(i)
	     endif
	   enddo
	enddo
         close (16)

        if (indprior.ne.1) then
           do i=1,np
              origmod(i)=0.0
              crrntmodb(i)=0.0
           enddo

          do i=1,nobs
           if (indaniso.eq.1) then
	     d(i)= velsum1(i)/sumstdev1(i)
	   else
	     d(i)= velsum2(i)/sumstdev2(i)
	   endif
          enddo
	else
	  if (indaniso.eq.1) then
	    do i = 1,np
	      origmod(i) = stcs2th(i)
	      crrntmodb(i) = origmod(i)
	    enddo
	  else
	    do i = 1,np
	      origmod(i) = stsn2th(i)
	      crrntmodb(i) = origmod(i)
	    enddo
	  endif
c  calculate expected effects of starting model
          do i = 1, nobs
	    gchange = 0.0
	    do jc = 1,np
	      gchange = gchange + g(i,jc)*origmod(jc)
	    enddo
            if (indaniso.eq.1) then
	      d(i)= velsum1(i)/sumstdev1(i)-gchange
c	      write(*,*) lontemp2,lattemp2,i,velsum1(i),d(i)*sumstdev1(i)
	    else
	      d(i)= velsum2(i)/sumstdev2(i)-gchange
c	      write(*,*) lontemp2,lattemp2,i,velsum2(i),d(i)*sumstdev1(i)
	    endif
	  enddo
	endif
         	    
       		
c  Calculate gtg and gtd
        do j = 1, np
          gtd(j) = 0.0
          do i = 1, nobs
            gtd(j) = gtd(j) + g(i,j)*d(i)
          enddo
           cmminvmmo(j) = 0.0
	     
c  construct gtg  
          do jj = 1,j
            gtg(jj,j) = 0.0
            do i = 1, nobs
              gtg(jj,j)= gtg(jj,j) + g(i,jj)*g(i,j)
            enddo
            gtg(j,jj) = gtg(jj,j)
            savegtg(j,jj) = gtg(j,jj)
            savegtg(jj,j) = gtg(jj,j)
          enddo
         enddo

c  Combination of minimum length and curvature
c   add smoothness constraint (minimum curvature).  Find coefficient, smthco, 
c  that minimizes, in least squares sense, the off-diagonal terms of gtg,
c   which should lead to minimal off-diagonal terms of covariance matrix
c
        
	if (inpcov.eq.0) then
	   sumsmth2 = 0.0
	   sumsmgg = 0.0
	   do ismth = 2,np
	     do jsmth = 1, ismth-1
	       sumsmth2 = sumsmth2 + rough(ismth,jsmth)**2
	       sumsmgg = sumsmgg + rough(ismth,jsmth)*
     1                        gtg(ismth,jsmth)
             enddo
	   enddo
	   smthco = -sumsmgg/sumsmth2
	   smthco = smthco*smthmu
           do j = 1,np 
             do i = 1,np   
               gtg(i,j) = gtg(i,j) +covinv(i,j)+smthco*Rough(i,j)
             enddo
c   add to gtd Tarantola term penalizing misfit to original starting model
c  For one-sided correction, do not penalize changes from original starting model in
c  the iteration
c           do i=1,np
c             tempcov = (covinv(j,i)+smthco*Rough(j,i))
c     1	               *(crrntmodb(i)-origmod(i))
c             cmminvmmo(j) = cmminvmmo(j) + tempcov
c	   enddo	  
             gtdcmm(j) = gtd(j) - cmminvmmo(j)
           enddo
	 else
           do j = 1,np 
             do i = 1,np   
               gtg(i,j) = gtg(i,j) +covinv(i,j)
             enddo
             gtdcmm(j) = gtd(j) - cmminvmmo(j)
           enddo
	 endif
	 

*****************************************************************************************

c  Invert gtg. gtg will be destroyed.  
C  Not the most efficient approach because doesn't take advantage of 
c  symmetry of gtg.  Use LU decomposition from Press et al.
        do i= 1,np
          do j = 1, np
            gtginv(i,j)= 0.0D0
          enddo
            gtginv(i,i) =1.0D0
        enddo
        
        call dludcmp(gtg,maxlay,np,nparam,indx,ddd)
        do j = 1,np
        call dlubksb(gtg,maxlay,np,nparam,indx,gtginv(1,j))
        enddo

      
c  Find change to starting model
        do i= 1, np
           change(i)=0.0
           do j = 1,np
              change(i) =change(i) + gtdcmm(j)*gtginv(i,j)
           enddo
c	   write(*,*) change(i), gtdcmm(i),gtginv(i,i)
        enddo
c	write(*,*) 'change'
c  Find rank (sum of diagonals of resolution matrix), i.e., number of
c  pieces of information or number of independent model parameters

        rank = 0.0
	
	do i=1,np
          resparm =0.0
          do j = 1,np
            resparm = resparm + gtginv(i,j)*savegtg(j,i)
          enddo
	    eachrank(i) = resparm
            rank = rank + resparm
        enddo 
	
c calculate entire resolution matrix
c can plot rows of resmatx with depth and see resolution kernals
	do k = 1,np
	   do i = 1,np
  	      tempmatx = 0.0
	      do j = 1,np
		tempmatx = tempmatx + gtginv(k,j)*savegtg(j,i)
	      enddo
	      resmatx(k,i) = tempmatx	
	   enddo
	enddo
        
c write resmatrix 
        write(52,*) lon, lat, label(indaniso)
	do i=1,np
	   write(52,152) i, (resmatx(i,j), j= 1,np)
	enddo
152     format(i2,1x,25f6.3)
c Find data importances, diagonals of g*gtginv*gt
c g*gtginv*gt : data resolution matrix
         
	do j=1,np
	   do i=1,nobs 
		tempimp(i) = 0.0
	 	do k=1,np
		   tempimp(i) = tempimp(i) + gtginv(j,k)*g(i,k)
		enddo
		tempgtginvgt(j,i) = tempimp(i)
	   enddo
	enddo

	do i=1,nobs
	   dataimp(i)=0.0
	   do j=1,np
	      dataimp(i) = dataimp(i) + g(i,j)*tempgtginvgt(j,i)
	   enddo
	enddo

        
c  Update current model  

         do ii = 1, np
          crrntmodb(ii) = crrntmodb(ii) + change(ii)
         enddo
         if (abs(stbeta(1)).lt.0.1) crrntmodb(1) = 0.0
        do j = 1, np
          stddev(j) = sqrt(gtginv(j,j))
        enddo

	
230   format (3x,f7.3,5x,f5.3,5x,f5.3,5x,f5.3,3x,f8.3)

         close(30)
        
	  write(18,*) lon,lat, label(indaniso) 
	  write(18,*) 'rank', real(rank)
	  write(18,*) 'smthco', smthco
	  write(18,*)'depth  beta   stddev  abschange change  each rank'

	do ii = 1, np
          write(18,180)depth(ii), crrntmodb(ii), stddev(ii),  
     r	  (crrntmodb(ii)-origmod(ii)), change(ii), eachrank(ii)  
c          write(19,*) depth(ii), crrntmodb(ii), stddev(ii), change(ii)
		  
        enddo

180    format(F6.2,5F8.4)  

	write(18,*) 'data importance'
	do i=1,nobs
	   write(18,*) i,real(dataimp(i))
	enddo

c calculate residuals for the new shear wave vel. profile from predicted
c  effect of changes   
        do i = 1, nobs
	  gchange = 0.0
	  do jc = 1,np
	    gchange = gchange + g(i,jc)*change(jc)
	  enddo
	  if (indaniso.eq.1) then
            residout(i) = (d(i)-gchange)*sumstdev1(i)
	    ccobs(i) = velsum1(i)
	    ccpredout(i)= ccobs(i) - residout(i)
	    stddevdata(i) = sumstdev1(i)
	  else
	    residout(i) = (d(i)-gchange)*sumstdev2(i)
	    ccobs(i) = velsum2(i)
	    ccpredout(i)= ccobs(i) - residout(i)
	    stddevdata(i) = sumstdev1(i)
	  endif
	enddo

	  if (indaniso.eq.1) then
	    write(18,*)
     1  'T    obsCos2Th    pred   datastddv  residuals  norm residuals'
           else
	    write(18,*)
     1  'T    obsSin2Th    pred   datastddv  residuals  norm residuals'
	   endif
         sumresid  = 0.0
        do i=1,nobs
            write(18,190)perd(i),ccobs(i),ccpredout(i),stddevdata(i),
     &           residout(i),residout(i)/stddevdata(i) 

            sumresid = sumresid + (residout(i)/stddevdata(i))**2
	enddo
	    RMS = sqrt(sumresid/nobs)
            standardresid = sqrt(sumresid/(nobs-rank))  !!!!nobs-rank is degree of freedom of models
	    sumnormmsft2 = sumnormmsft2 + standardresid**2
	    write(18,*)                          
	    write(18,*) '   RMS    stddv(divided by degree of freedom)'     !!!! for model
	    write(18,*) RMS,  standardresid 
	    write(18,*)  
190     format(F6.2, 2F10.4, 3F10.4) 
      
c  now if specified, invert predictions to obtain smoother model - use already generated partial
c  derivatives, etc.  

        if (invprd.eq.1) then
c  first, undo updating of model
          do ii = 1, np
           crrntmodb(ii) = crrntmodb(ii) - change(ii)
          enddo


          do i=1,nobs
	    gchange = 0.0
	    do jc = 1,np
	      gchange = gchange + g(i,jc)*origmod(jc)
	    enddo
            if (indaniso.eq.1) then
	      d(i)= ccpredout(i)/sumstdev1(i)-gchange
	    else
	      d(i)= ccpredout(i)/sumstdev2(i)-gchange
	    endif
          enddo		
c  Calculate  gtd
          do j = 1, np
            gtd(j) = 0.0
            do i = 1, nobs
              gtd(j) = gtd(j) + g(i,j)*d(i)
            enddo
            cmminvmmo(j) = 0.0
          enddo
	  if (inpcov.eq.0) then
            do j = 1,np 
c   add to gtd Tarantola term penalizing misfit to original starting model
              gtdcmm(j) = gtd(j) - cmminvmmo(j)
            enddo
	  else
            do j = 1,np 
              do i = 1,np   
                gtg(i,j) = gtg(i,j) +covinv(i,j)
              enddo
              gtdcmm(j) = gtd(j) - cmminvmmo(j)
            enddo
          endif
c  Find change to starting model
          do i= 1, np
            change(i)=0.0
            do j = 1,np
              change(i) =change(i) + gtdcmm(j)*gtginv(i,j)
            enddo
          enddo
c  Update current model  

          do ii = 1, np
            crrntmodb(ii) = crrntmodb(ii) + change(ii)
c  remember cs2theta terms for later output (put in crrntmoda)
	    if (indaniso.eq.1) crrntmoda(ii) = crrntmodb(ii)
          enddo
          if (abs(stbeta(1)).lt.0.1) crrntmodb(1) = 0.0
	  write(18,*) ' redone for smoother model fitting predictions'
	  write(18,*) lon,lat, label(indaniso) 
	  write(18,*) 'rank', real(rank)
	  write(18,*) 'smthco', smthco
	  write(18,*)'depth  beta   stddev  abschange change  each rank'

	  do ii = 1, np
            write(18,180)depth(ii), crrntmodb(ii), stddev(ii),  
     1	    (crrntmodb(ii)-origmod(ii)), change(ii), eachrank(ii)  		  
          enddo
c calculate residuals for the new shear wave vel. profile from predicted
c  effect of model   
        do i = 1, nobs
	  gchange = 0.0
	  do jc = 1,np
	    gchange = gchange + g(i,jc)*crrntmodb(jc)
	  enddo
	  if (indaniso.eq.1) then
c            residout(i) = (d(i)-gchange)*sumstdev1(i)
	    ccobs(i) = velsum1(i)
c	    ccpredout(i)= ccobs(i) - residout(i)
            ccpredout(i)= gchange*sumstdev1(i)
c  remember prediction in second array for later output
	    ccpredout2(i) = ccpredout(i)
	    residout(i) = ccobs(i) - ccpredout(i)
	    stddevdata(i) = sumstdev1(i)
	  else
c	    residout(i) = (d(i)-gchange)*sumstdev2(i)
	    ccobs(i) = velsum2(i)
            ccpredout(i)= gchange*sumstdev2(i)
	    residout(i) = ccobs(i) - ccpredout(i)
c	    ccpredout(i)= ccobs(i) - residout(i)
	    stddevdata(i) = sumstdev1(i)
	  endif
	enddo

	  if (indaniso.eq.1) then
	    write(18,*)
     1  'T    obsCos2Th    pred   datastddv  residuals  norm residuals'
           else
	    write(18,*)
     1  'T    obsSin2Th    pred   datastddv  residuals  norm residuals'
	   endif
         sumresid  = 0.0
        do i=1,nobs
            write(18,190)perd(i),ccobs(i),ccpredout(i),stddevdata(i),
     &           residout(i),residout(i)/stddevdata(i) 

            sumresid = sumresid + (residout(i)/stddevdata(i))**2
	enddo
	    RMS = sqrt(sumresid/nobs)
            standardresid = sqrt(sumresid/(nobs-rank))  !!!!nobs-rank is degree of freedom of models
	    sumnormmsft2 = sumnormmsft2 + standardresid**2
	    write(18,*)                          
	    write(18,*) '   RMS    stddv(divided by degree of freedom)'     !!!! for model
	    write(18,*) RMS,  standardresid 
	    write(18,*)  
	
	endif

        write (50,*) inxy, lon,lat,label(indaniso)  
        do i = 1,np
	  do j= 1,i
	    write(50,*) i,j,gtginv(i,j)
	  enddo
	enddo  
        indaniso = indaniso +1
	if (indaniso.eq.2) go to 2000
	
	do i = 1, nper
c	  write(60+i,*) lon,lat,velsum1(i),velsum2(i)
	  write(60+i,*) lon,lat,ccpredout2(i),ccpredout(i)
	enddo

	do ii = 1, np
	   write(90,290) thick(ii),crrntmoda(ii),crrntmodb(ii)  
	enddo	    
290     format (f6.3,2x,f5.3,2x,f6.3,2x,f6.3)


	enddo  ! close all grid points
	
	avgnormmsft2 = sumnormmsft2/(2*nxy) 
	write(50,*) avgnormmsft2, ' avgnormmsft2'
	write(18,*) avgnormmsft2, ' avgnormmsft2'
	

	close(10)
        close(unit = 12)
        close(18)
	close(50)
	close(51) 
        close(52)
	do i = 1, nper
	  close(60+i)
	enddo
	close(90)
      stop
      end

      SUBROUTINE dlubksb(a,maxlay,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(maxlay,maxlay),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END


      SUBROUTINE dludcmp(a,maxlay,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(maxlay,maxlay),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

