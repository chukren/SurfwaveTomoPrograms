!=======================================================================================
!  srchwave589.isores.nophase.JdF.f90 is fortran90 version that outputs resolution 
!  matrix for isotropic velocity parameters
!
!  srchwave589.kern.br.d.phc.kres.f  differs from 587 in solving for both aniso
!  and isotropic lateral variations 
!                                  
!  srchwave587.kern.br.d.phc.kres.f   
!  differs from 586 by correcting for first order effect of lateral velocity changes
!  on the sensitivity kernels - does not recalculate shape of kernel, but multiples
!  sensitivity by local correction factor                                             
!
!  586   differs from 584 -  does not use uniform average velocity to correct DC level
!      of starting velocity model  CJR 2/22/09
!                            DWF  08/14/09  10/01/09  11/19/09   12/5/09   3/11/10     
!                            YYR  04/29/18
!========================================================================================
!  calculates fit after linearized inversion to test whether predicted improvement is 
!  actually achieved
!
!  replaces simulated annealing search for wave parameters with grid search
!  INCLUDES wave parameters in linearized inversion
!  Differs from 534 by replacing Gaussian interpolation with linear interpolation
!  Differs from 544 by generating grid internally
!  Differs from 554 by correcting predicted starting phase velocities on fine grid
!  internally (crustal grid)
!  Differs from 564 by using minimum curvature criterion (smoothness) to
!  minimize off-diagonal elements of a posteriori model covariance matrix (actually,
!  minimize off-diagonal components of GtG matrix

!
!  Uses xmin constraint on using node sensitivities
!
!  this differs from standard 134 by damping amplitude corrections more in first
!  five iterations 
!
!  Outer grid points damped same as interior
!
!  data: observed phases and amplitudes at each station for each event
!
!  134 means to invert 2-D isotropic phase velocity ,station site responses 
!  and anisotropy in each subregions, with starting values for node coefficients
!  read in based on a priori model. phc means phase corrections,only  
!  for stations which have erroneous instrument phase response. 
!  
!  kern means using sensitivity kern to calculate the amplitude and phase 
!  at each station. The sensitivity kernel is calculated in this program, using 
!  input about frequencies involved from windowing. With the 
!  sensitivity kernel data, we can know the sensitivity at each node by 
!  interpolating the sensitivity kernel. 
!
!  anisotropy affects scattering
!     
!
!  Use minimum length criterion. Does not smooth but penalizes changes, 
!  with larger variance on outer edges, so that edges can take up
!  residual junk without propagating effects into the interior
!  Use tarantola terms to damp solution,grid of node points 
!  to interpolate velocity, and two plane waves per event
!
!  Run with piped input file like srchwave134.phc < eqlister50.
!  eqlister50 is an example file listing earthquake event at period 50s.
!  
!  This version reads in non-uniform starting model
!=============================================================================

      parameter (maxnfreq=1, maxnsta=300, maxpts=20000, nparam = 4000, &
                maxnobs = 25000, maxnodes = 2000, maxevnts = 400, &
                maxnxints = 401, maxnyints = 401,ndeg = 81 )
      
      parameter (mxfreqkern = 200)

      real*4 dph(maxevnts,maxnsta)
      real*4 staph(maxevnts,maxnsta,maxnfreq)
      real*4 staamp(maxevnts,maxnsta,maxnfreq)
      real*4 freq(maxnfreq)
      real*4 g(maxnobs,nparam)
      real*4 isores(maxnodes,maxnodes)
      real*4 stddevdata(maxevnts)
      real*4 gtd(nparam), stddev(nparam), covinv(nparam)
      real*4 stadist(maxevnts,maxnsta), staazi(maxevnts,maxnsta)
      real*4 ysta(maxevnts,maxnsta), xsta(maxevnts,maxnsta)
      real*4 streal(maxevnts,maxnsta,maxnfreq)
      real*4 stimag(maxevnts,maxnsta,maxnfreq)
      real*4 rloc(maxevnts,maxnsta), azloc(maxevnts,maxnsta)
      real*4 attnfac(maxnfreq), xmin(maxevnts,ndeg)
      real*4 amprms(maxevnts,maxnfreq)
      real*4 bazi(maxevnts,maxnsta), stadelt(maxevnts,maxnsta)
      real*4 d(maxnobs)
      real*4 beg(maxevnts,maxnsta)
      real*4 origmod(nparam),crrntmod(nparam)
      real*4 boxlat(4), boxlon(4),applat,applon
      real*4 nodelat(maxnodes),nodelon(maxnodes), nodevel(maxnodes)
      real*4 nodecos2(maxnodes), nodesin2(maxnodes)
      real*4 xbox(maxevnts,4),ybox(maxevnts,4)
      real*4 xnode(maxevnts,maxnodes),ynode(maxevnts,maxnodes)
      real*4 stalat(maxevnts,maxnsta),stalon(maxevnts,maxnsta)
      real*4 tazim, delta, tazimref, bazimnd(maxnodes)
      real*4 adistsq,wgttemp(maxnodes)
      real*4 wgtnode1(maxnsta,maxnodes,ndeg)
      real*4 ampwgtnode1(maxnsta,maxnodes,ndeg) 
      real*4 dtime(maxnsta),avslow(maxnsta)
      real*4 dtime1(maxnsta),dtime2(maxnsta)
      real*4 avslow1(maxnsta), avslow2(maxnsta)
      real*4 appvel(maxnodes),vage(maxnodes)
      real*4 cos2node(maxnodes,ndeg),sin2node(maxnodes,ndeg)
      real*4 startamp1(maxevnts),startamp2(maxevnts),stazim1(maxevnts)
      real*4 stazim2(maxevnts),stphase1(maxevnts),stphase2(maxevnts)
      real*4 pv(6),p(7,6),smsft(7),pb(6),annlscl(6),ppb(6),pppb(6)
      real*4 phcor(50)                                            
      real*4 minstd, rmsphase(maxevnts),rmsamp(maxevnts)
      real*4 sortrms(maxevnts)
      real*4 rmsdata(maxevnts)
      real*4 damp1per(maxevnts,maxnsta), damp2per(maxevnts,maxnsta)
      real*4 dphase1(maxevnts,maxnsta),dphase2(maxevnts,maxnsta)
      real*4 phase1(maxevnts,maxnsta), phase2(maxevnts,maxnsta)
      real*4 phase(maxnsta,ndeg),dphase(maxnsta,ndeg), 
      real*4 dampper(maxnsta,ndeg)
      real*4 unifvel,ampmult(maxnsta)
      real*4 dxnode,dynode,xbegkern,dxkern
      real*4 ybegkern,dykern
      real*4 residdiag(maxnodes) 
      real*4 subgtd(nparam),substddev(nparam), subcovinv(nparam)
      real*4 gchange(maxnobs)
      real*4 covar(maxnodes,maxnodes),vchange(maxnodes)
      
      real*4 b(maxnodes,maxnodes),btb(maxnodes,maxnodes)

           
      integer*4 nxkern,nykern
      integer*4 nsta(maxevnts), iref(maxevnts),idnode(maxnodes)
      integer*4 nfreq,nstapts(maxnsta),indx(nparam),nstacor,nobs
      integer*4 istanum(maxevnts,maxnsta),istacor(maxnsta)
      integer*4 nnodes, nevents, ntyp1
      integer*4 ityp1sta(maxnsta)   
      character*12 idnum(maxevnts)                            
      integer*4 nevntsta(maxnsta),istavar(maxnsta)  
      integer*4 jbox(maxevnts)
      
      double precision change(nparam), gtdcmm(nparam)
      double precision gtg(nparam,nparam),gtginv(nparam,nparam),ddd
      double precision savegtg(nparam,nparam),chmax
      double precision subgtg(nparam,nparam),subgtdcmm(nparam)
      double precision subgtginv(nparam,nparam)
      double precision subsavegtg(nparam,nparam)
      
      real*8 sensitivity(maxnxints,maxnyints)
      real*8 ampsens(maxnxints,maxnyints)
      real*4 amsens(maxnxints,maxnyints),phsens(maxnxints,maxnyints)
      real*4 kk,lamda
      real*4 amplitude(mxfreqkern), freqkern(mxfreqkern)
      
      ! Add by YYR for Q kernels 
      real*4 kern_am2qinv(maxnxints,maxnyints) 
      real*4 kern_ph2qinv(maxnxints,maxnyints)
      real*4 freq_ref, kern_Qfactor, kern_phi, kern_denom
      real*4 ph2phsv, am2phsv
      real*4 qinv, grpvel
      real*4 kern_freqdep
      
      character*2   nettemp
      character*4   statemp 
      character*97  dummy, dummy1 
      character*70  foutput,fsummary,fvariance
      character*70  finvrsnodes, fftinput
      character*70  fmaxavamp,fvelarea,ftemp, fstalist
      character*70  fvelout,fvelnodes, startvel
      character*70  sensfn
      character*70  dirpath
      character*70  fresdiag,fresmat
      character*70  fendvel
        
      character*100 inputfile
        
      character*100 buffer
      
!      character*75  dummy1
!      character*100 foutput,fn(maxevnts,maxnsta),fsummary
!  ******* changes in this area from yun version
!      character*100 foutput,fsummary

      character*4   staname(maxnsta)
      character*110 fn(maxevnts,maxnsta)

      common /residua/ sensitivity,ampsens,d,rloc,azloc,freq, &
           xsta,dtime,avslow,streal,stimag,stddevdata,phcor, &
           xbox,ybox,ysta,dxkern,dykern,dxnode,dynode, &
           unifvel,appvel,ampmult,gamma, &
           xnode, ynode,wgtnode1,ampwgtnode1,xmin, &
           phase,dphase,dampper,istavar,istanum,nxkern,nykern, &
           ityp1sta,ntyp1,iref,nsta,iev,naddat,ifreq,nnodes

      common /msft/ bazi,cos2node,sin2node,crrntmod, &
           startamp1,startamp2,stphase1,stphase2,stazim1,stazim2, &
           nevents,idnode,ideg1,ideg2,nobs,i6,iarea
     
      common /gengrid/ nodelat,nodelon,boxlat,boxlon
       
      common /update2/ covar,vchange


      pi = 3.1415928
      convdeg = 3.1415928/180.
      circ = 6371.*3.1415928/180.
      twopi = 3.1415928*2.

      ! Note we choose 50 s rather than 1 sec for the reference frequency
      ! Therefore the inversion results cannot be directly compare with 
      ! PREM model or other model which normally obtained at 1 Hz

      freq_ref = 0.02 ! Hz 
  

! Corresponding ID for abnormal stations.                   
! ntyp1: number of stations need phase correction
      ntyp1 = 0 

!        data (ityp1sta(ityp), ityp=1, 4 )/1,3,172,174/
      if (ntyp1.ne.0) then              
        ityp1sta(1) = 1
        ityp1sta(2) = 3
        ityp1sta(3) = 172
        ityp1sta(4) = 174
      endif

200   format(a75, i2)
201   format(a70, a2)
202   format(a69, a4) ! precise cut the STANM out from the path file

!  read list of files and frequencies to be analyzed and files to output results
!  Usually will pipe in data from some file like eqlistper50. 

      print *, "bp 1"  

      call getarg(1, buffer)
      read(buffer,'(a)') inputfile
      write(*,*) "read", inputfile
        
      open(99, file = inputfile)

      read(99,*) nevents
      write(*,*) nevents

      nobs = 0
      do iev = 1, nevents
        read(99,*) nsta(iev),idnum(iev)
        nobs = nobs+ 2*nsta(iev)

        do i = 1, nsta(iev)
          read(99,'(a)') fn(iev,i)
        enddo

      enddo

      read(99,*) nfreq
      read(99,*) (freq(j), j= 1, nfreq) ! Hz
      read(99,*) foutput
      read(99,*) fsummary
      read(99,*) finvrsnodes
      read(99,*) fftinput
      read(99,*) fvariance
      read(99,*) fmaxavamp      
      read(99,*) ftemp
      read(99,*) fvelnodes
      read(99,*) fstalist
      read(99,*) unifvel       ! km/s
      read(99,*) qinv, grpvel  ! reference 1/Q and group velocity
      read(99,*) iterlimit, dampvel, dampaniso, divfac
      read(99,'(a)') sensfn
      read(99,*) fvelout
      read(99,*) startvel
      read(99,*) fresdiag
      read(99,*) fendvel
      read(99,*) fresmat
      print *, "bp 1c"


!  unifvel is average velocity found from inverting data previously for a single velocity
!  parameter
!   iterlimit is maximum number of iterations to use
!  dampvel is a priori stddev for velocity terms, dampaniso for aniso coeff.

      !  read data
      open(10, file = foutput)
      print *, "bp 1d"
      open(11, file = fsummary)
      open(12, file = fftinput)
      open(13, file = ftemp)
      open(14, file = "followit12")
      open(16, file = fvariance)
      open(18, file = fstalist)
      open(30, file = fvelout)
      open(66, file = sensfn)
      open(21, file = fresdiag)
      open(17, file = fresmat)
!      open(15, file = finvrsnodes)
!      open(60, file = startvel)
!      open(66, file = dirpath//sensfn)

      print *, "bp 2"
      !  do following step to extract station name from filename,
      !  count number of events at each station and assign a station number 
      !
      !  first, read in master list of stations   

      do ista = 1, maxnsta
        nevntsta(ista) = 0
        read(18,*) staname(ista)
        if (staname(ista).eq.'nope') then
          mxnsta = ista -1
          go to 1111
        endif
      enddo
1111  continue

      print *, "bp 3"
      do iev = 1, nevents
        do ista = 1, nsta(iev)

          write(13,'(a)') fn(iev,ista)
          rewind (13)
         
          !  ******** 
          !  modify here depending on how you cut input line 
          !  to assign your station ID precisely. 
          !  Need to be improved
          !  **************
          read (13,202) dummy1, statemp

          rewind (13)
          !print *, statemp
          istanum(iev,ista) = 0

          do ista2 = 1, mxnsta
            if (statemp.eq.staname(ista2)) then
              istanum(iev,ista) = ista2
              nevntsta(ista2) = nevntsta(ista2) +1
            endif
          enddo
      
          if (istanum(iev,ista).eq.0) then
            write(*,*) 'WARNING ', statemp,' not in station list'
          endif
        
          ! write (*,*) ista, istanum(iev,ista),statemp
        enddo
      enddo

      ! count number of stations that have events and assign each 
      ! station number a variable number

      jstacnt = 0
      do ista = 1, mxnsta
        if (nevntsta(ista).gt.0) then
          jstacnt = jstacnt+1
          istavar(ista) = jstacnt
        endif
      enddo


!  read in location of velocity nodes, box for limits of tomography
!  position in decimal degrees 

!      read(15,*) dummy
!      read(15,*) nnodes
!      do i = 1,nnodes	
!        read(15,*) nodelat(i),nodelon(i)
!      enddo
!      do i = 1, 4
!      read(15, *) boxlat(i), boxlon(i)	
!      enddo
!      read(15, *) ncol
!      read(15, *) dxnode,dynode

      call genreggrid(finvrsnodes,nnodes,ncol,dxnd,dynd)
      dxnode = abs(dxnd)
      dynode = abs(dynd)  
      
!      do i = 1,nnodes	
!        write(*,*) i, nodelat(i),nodelon(i)
!      enddo
!      do i = 1, 4
!        write(*, *) boxlat(i), boxlon(i)	
!      enddo
!      write(*, *) ncol
!      write(*, *) dxnode,dynode

!     generate sensitivity kernel
!     first, read in spectral file of frequencies and relative amps involved in particular
!     window/period processing for this target frequency
      print *, "bp 4"

      radius = 6371.
      read(66,*) nxkern, xbegkern, dxkern
      read(66,*) nykern, ybegkern, dykern
      
      if ((nxkern.gt.maxnxints).or.(nykern.gt.maxnyints)) then
        write (*,*) "Too many points in sensitivity kernel"
      endif

      read(66,*) nkernfreq
      sumkernamp = 0.0

      do ifreqkern = 1, nkernfreq
        read(66,*) freqkern(ifreqkern), amplitude(ifreqkern)
      enddo

      !  add extra, interpolated frequencies	
      do ifreqkern = 1, nkernfreq-1
        amplitude(ifreqkern+nkernfreq) = 0.5*(amplitude(ifreqkern) &
                                + amplitude(ifreqkern+1))
        freqkern(ifreqkern+nkernfreq)  = 0.5*(freqkern(ifreqkern) &
                                +freqkern(ifreqkern+1))
      enddo

      nkernfreq = 2*nkernfreq -1
      
      do ifreqkern = 1, nkernfreq
        sumkernamp = sumkernamp + amplitude(ifreqkern)
      enddo
      
      do ifreqkern = 1,nkernfreq
        amplitude(ifreqkern) = amplitude(ifreqkern)/sumkernamp
      enddo   

      do ixkern = 1,nxkern
        do iykern = 1,nykern
          amsens(ixkern,iykern) = 0.0
          phsens(ixkern,iykern) = 0.0

          kern_am2qinv(ixkern,iykern) = 0.0
          kern_ph2qinv(ixkern,iykern) = 0.0

          sensitivity(ixkern,iykern) = 0.
          ampsens(ixkern,iykern) = 0.
        enddo
      enddo

      print *, "bp 4a"
      

      do ifreqkern = 1, nkernfreq
        period = 1.0/freqkern(ifreqkern) 
        lamda = unifvel*period
        kk = 2.0*pi/lamda*radius

        kern_Qfactor = 0.5*unifvel*qinv/grpvel 
        kern_freqdep = 2.0 / onepi * log(1./period / freq_ref)

        do ixkern = 1,nxkern
          xkern = (ixkern-1)*dxkern + xbegkern
          delta1 = xkern

          do iykern = 1,nykern
            ykern = (iykern-1)*dykern + ybegkern

            if(xkern.eq.0 .and. ykern.eq.0) then 
              delta2 = sqrt(dxkern**2+dykern**2)
            else
              delta2 = sqrt(xkern**2+ykern**2)
            endif

            ! 2D velocity kernels are from Zhou et al. (2004), GJI, "3-D sensitivity 
            !  kernels for surface-wave observables"

            kern_phi = kk*(delta1+delta2)/radius + pi/4.
            kern_denom = sqrt(8.*pi*kk*abs(sin(delta2/radius)))

            ph2phsv = (-2.)* kk**2. * sin(kern_phi) / kern_denom
            am2phsv = (-2.)* kk**2. * cos(kern_phi) / kern_denom

            ph2phsv = ph2phsv * amplitude(ifreqkern) 
            am2phsv = am2phsv * amplitude(ifreqkern)

            phsens(ixkern,iykern) = phsens(ixkern,iykern) + &
                      ph2phsv * ((dxkern*dykern)/radius**2) 

            amsens(ixkern,iykern) = amsens(ixkern,iykern) + &
                      am2phsv * ((dxkern*dykern)/radius**2)

            ! 2D anelasticity (1/Q) kernels are the combinations of velocity kernels
            ! Zhou (2009), GJI 
            
            kern_am2qinv(ixkern,iykern) = kern_am2qinv(ixkern,iykern) &
                      + kern_Qfactor*(kern_freqdep*ph2phsv - am2phsv)
            
            kern_am2qinv(ixkern,iykern) = kern_am2qinv(ixkern,iykern) &
                      + kern_Qfactor*(ph2phsv + kern_freqdep*am2phsv)
              
          enddo
        enddo
      enddo

      print *, "bp 5"   

!  now take these raw sensitivity kernels and find net effect for node velocities using
!  linear interpolation between nodes, essentially convolving interpolater with raw kernels
      ixoff = dxnode/dxkern
      iyoff = dynode/dykern

      do ixkern = 1, nxkern
        ixxbeg = ixkern - ixoff
        ixxend = ixkern + ixoff

        if (ixxbeg.lt.1) ixxbeg = 1
        if (ixxend.gt.nxkern) ixxend = nxkern

        xk = (ixkern-1)*dxkern + xbegkern  
        
        do iykern = 1,nykern  
            iyybeg = iykern - iyoff
            iyyend = iykern + iyoff

            if (iyybeg.lt.1) iyybeg = 1
            if (iyyend.gt.nykern) iyyend = nykern
            
            yk = (iykern-1)*dykern + ybegkern
            wgtsum = 0.0

            do ix = ixxbeg, ixxend
              xkk = (ix-1)*dxkern + xbegkern
              do iy = iyybeg,iyyend
                ykk = (iy-1)*dykern + ybegkern
                wgtkern = (1.-abs(xk-xkk)/dxnode)*(1.-abs(yk-ykk)/dynode)
                sensitivity(ixkern,iykern) = sensitivity(ixkern,iykern) &
                                  + phsens(ix,iy)*wgtkern
                ampsens(ixkern,iykern) = ampsens(ixkern,iykern) &
                                  + amsens(ix,iy)*wgtkern
                wgtsum = wgtsum+wgtkern
              enddo
            enddo

            sensitivity(ixkern,iykern) = sensitivity(ixkern,iykern)/wgtsum
            ampsens(ixkern,iykern) = ampsens(ixkern,iykern)/wgtsum
!	write(14,*) xk,yk,sensitivity(ixkern,iykern),ampsens(ixkern,iykern), &
!              phsens(ixkern,iykern),amsens(ixkern,iykern)
        enddo
      enddo

!  renormalize sensitivity kernels
      sensnorm =  (dxnode*dynode)/(dxkern*dykern)
      
      do ixkern = 1,nxkern
        do iykern = 1,nykern
         sensitivity(ixkern,iykern) = sensitivity(ixkern,iykern)*sensnorm
         ampsens(ixkern,iykern)= ampsens(ixkern,iykern)*sensnorm
        enddo
      enddo

      write (*,*) 'done with sensitivity'

! assign velocities to each grid node apriori
      write(*,*) 'bp 5a'
      write(*,*) startvel
      
      call assignstrt(startvel,nodevel,preunifvel,nnodes, &
                                 ncol,dxnode,dynode)
      
      write(*,*) 'bp 5b'

      do i = 1, nnodes
!  read in starting node velocities based on a priori crustal model
!  but modify velocities for offset between average a priori and 
!  observed average from uniform velocity inversion - lateral velocity
!  variations will be preserved, but reference velocity will be observed average

!c         nodevel(i)=unifvel

!        read(60,*) xx1,xx2, nodevel(i)
!
!  following line is where 586 differs from 584 - without this line, the second following
!  line would shift the DC level of the starting predicted phase velocity model to agree with  
!  average velocity inferred previously from a 1-D inversion
        unifvel = preunifvel
        nodevel(i) = nodevel(i) + unifvel - preunifvel
!	write(14,*) nodelat(i),nodelon(i),nodevel(i)
      enddo


! Can input gridnode file to solve anisotropy in different areas.
!	open(20, file='gridnode.aniso.id')
!	open(20, file='invnodesethp_area2_test_ridge.d')
!	read(20,*) iarea
!
!  In version 589, iarea = nnodes
!	 iarea = 1
      iarea = nnodes
      
      do i=1, nnodes
!	 read(20,*)ix, xx, xx, idnode(i)  			
!	 idnode(i) = 1 
        idnode(i) = i
      enddo

      print *, "bp 6"
      
      write(10, *) foutput
      write(11, *) foutput
      
      ! start input loop over events
      do iv = 1, nevents
        read(12,*) iev

        do ista = 1, nsta(iev)
          read(12,*) beg(iev,ista)
          read(12,*) stadist(iev,ista),staazi(iev,ista), bazi(iev,ista), &
          stadelt(iev,ista), stalat(iev,ista), stalon(iev,ista)
          read(12,*) (staamp(iev,ista,ifreq),staph(iev,ista,ifreq), &
                                     ifreq=1,nfreq)
        enddo
      enddo
    
      print *, 'bp 7'

      ! end of input loop over events

      write(10,*) nfreq, '  nfreq'


!=======================================
!  begin inversion loop over frequencies
!  YYR: Beware of this thousand-line long loop 
!       which is completely unnecessary
!=======================================

      do ifreq = 1, nfreq
        write(10,*) freq(ifreq)
        write(11,*) freq(ifreq)
        write(14,*) freq(ifreq)
        write(*,*) freq(ifreq)
!  assume a priori data covariance = 0.2  Function of this constant value is
!  to make damping based on real estimates about the right size.
        do iev = 1, nevents
          stddevdata(iev) = 0.2
        enddo

!  set up model covariances or damping to be added to diagonal of gtg
!  variables are amplitudes, directions, initial phases for each event
!  followed by isotropic velocities, anisotropic velocities for each area
!  then amplitude correction factors, phase correction factors for anomalous 
!  stations and attenuation factor gamma
        
        npnoamp = 6*nevents + 1*nnodes+ 2*iarea 
        npp = npnoamp + jstacnt 
        np = npp+ntyp1 + 1                           
        kj = nnodes/ncol
        i6 = 6*nevents

!   ***************************************
!   WARNING - the previous is specific to grid of node points - specify iages
!   even if age zones not explicitly used in inversion - needed for 
!   later smoothing
        do ii = 1, np
          change(ii) = 0.0
        enddo

        do iev = 1, nevents
          ip = (iev-1)*6

          covinv(1+ip) = 1./(0.10**2)
          covinv(2+ip) = 1./(0.10**2)
          covinv(3+ip) = 1./((3.*convdeg)**2)
          covinv(4+ip) = 1./((3.*convdeg)**2)
          covinv(5+ip) = 1./(.05**2)
          covinv(6+ip) = 1./(.05**2)
!          covinv(1+ip) = 1./(0.0010**2)
!          covinv(2+ip) = 1./(0.0010**2)
!          covinv(3+ip) = 1./((.03*convdeg)**2)
!          covinv(4+ip) = 1./((.03*convdeg)**2)
!          covinv(5+ip) = 1./(.0005**2)
!          covinv(6+ip) = 1./(.0005**2)

          !  make covinv damping much tighter for grid search than for 
          !  simulated annealing so that linearized step following grid 
          !  search for wave parameters won't go far astray.  
        enddo

        do ii= 1, nnodes
          ip = i6 + ii
          origmod(ip) = nodevel(ii)
!          origmod(ip) = unifvel
          covinv(ip) = 1./(dampvel**2)
          crrntmod(ip) = origmod(ip)
        enddo

        ! anisotropy varies in different areas	
        do i =1, iarea
          ipp = i6 + nnodes + i
          ippp = ipp + iarea 
          origmod(ipp) = 0.0
          origmod(ippp) = 0.0
          covinv(ipp) = 1./(dampaniso**2)
          covinv(ippp) = 1./(dampaniso**2)
          crrntmod(ipp) = origmod(ipp)
          crrntmod(ippp) = origmod(ippp)
        enddo  

!  initialize station correction variables and attenuation factor
!  May become singular if don't constrain average multiplication factor
!  because can multiply all by same amount and will trade-off with 
!  amplitude factors perfectly - add data point with small standard 
!  deviation

        do ii = 1, jstacnt
          ampmult(ii) = 1.0
          ip = npnoamp+ii
          origmod(ip) = ampmult(ii)
          crrntmod(ip) = origmod(ip)
!         covinv(ip) = 1.0/(0.30**2)
          covinv(ip) = 1.0/(0.10**2)          
        enddo


!  initialize phase corrections for different instruments
        if (ntyp1.ne.0) then             
          do ii=1, ntyp1      
            phcor(ii)=0.0                                               
            ip=npp + ii                                                 
            origmod(ip) = phcor(ii)                                     
            crrntmod(ip) = origmod(ip)                                  
            covinv(ip) = 1.0/(0.1**2)                                   
          enddo
        endif


!  initialize gamma factor for attenuation
        gamma = 0.0
        origmod(np) = gamma
        crrntmod(np) = origmod(np)
        covinv(np) = 1.0/(0.0002**2)


!  increase the variance for edges for velocity 
        varfac2 = 1.
        do itp = 1,1
          ityp = nnodes*(itp-1)
          
          !  right end 
          do ijk = i6+1,i6+kj
            jk = ijk+ityp
            covinv(jk) = covinv(jk)/varfac2
          enddo

          !  left end
          do ijk = i6+nnodes-kj+1,i6+nnodes
            jk = ijk+ityp
            covinv(jk) = covinv(jk)/varfac2
          enddo

          !  top
          do ijk = i6+kj+1,i6+nnodes-2*kj+1,kj
            jk = ijk+ityp
            covinv(jk) = covinv(jk)/varfac2
          enddo

          !  bottom
          do ijk = i6+2*kj,i6+nnodes-kj,kj
            jk = ijk+ityp
            covinv(jk) = covinv(jk)/varfac2
          enddo
        enddo


        !  Set up coordinate systems and normalize amplitude for each event.
        !  Begin loop over events for coordinates

        do iev = 1, nevents

!  Find reference station that has largest amplitude for each event.  
!  At this station, interfering waves should be nearly in phase and 
!  amplitude = sum of amps of individual waves

!  Also find normalizing amplitude to equalize earthquakes of different
!  size, using rms amplitude (old version used largest amplitude)

          amplarge =0.0
          amprms(iev,ifreq) = 0.0
          iref(iev) = 1

          do ista = 1, nsta(iev)
            amprms(iev,ifreq) = amprms(iev,ifreq) &
                                  + staamp(iev,ista,ifreq)**2

            if (staamp(iev,ista,ifreq).gt.amplarge) then
              amplarge = staamp(iev,ista,ifreq)
              iref(iev) = ista
            endif

          enddo
!	write(*,*) iev,iref(iev)
          amprms(iev,ifreq) = sqrt(amprms(iev,ifreq)/nsta(iev))

  
!  Use reference station for to set up local coordinate system
!  based on reference station at zero, zero.  Use distances to stations rather
!  than absolute coordinates as way of correcting for curving wavefront.  This
!  will tend to favor keeping azimuth at original azimuth along great circle.
!  + x direction is in direction of propagation along great circle path.
!  +y is 90 deg counterclockwise from x
!   staazi are measured clockwise from north
          xsta(iev,iref(iev))  = 0.0
          ysta(iev,iref(iev))  = 0.0
          rloc(iev,iref(iev))  = 0.0
          azloc(iev,iref(iev)) = 0.0

          do ista = 1, nsta(iev)
            if (ista.ne.iref(iev)) then
              xsta(iev,ista) = stadist(iev,ista) - stadist(iev,iref(iev))
              azidiff = staazi(iev,iref(iev)) - staazi(iev,ista)
              
              if (azidiff.gt.180.) azidiff = azidiff -360.
              if (azidiff.lt.-180.) azidiff = azidiff +360.
              
              ysta(iev,ista) = circ*sin(stadelt(iev,ista)*convdeg)* &
                               azidiff
              rloc(iev,ista) = sqrt(xsta(iev,ista)*xsta(iev,ista) + &
                               ysta(iev,ista)*ysta(iev,ista))
              azloc(iev,ista) = atan2(ysta(iev,ista),xsta(iev,ista))
            endif
          enddo

!  calculate apparent pole position of earthquake for spherical earth
!  instead of real location so that coordinate system for nodes will agree
!  with that for stations 
          call gohead(stalat(iev,iref(iev)),stalon(iev,iref(iev)), &
               stadelt(iev,iref(iev)),bazi(iev,iref(iev)),applat,applon)
          call disthead(applat,applon,stalat(iev,iref(iev)), &
               stalon(iev,iref(iev)),delta,tazimref)
!  now calculate x,y of each node 
!  this approach seems to work pretty well - agrees within about 0.1% with
!  station calculations based on elliptical great circle distance in terms
!  of relative position
          appcirc = stadist(iev,iref(iev))/stadelt(iev,iref(iev))
      
          do inode = 1, nnodes
            call disthead(applat,applon,nodelat(inode),nodelon(inode) &
                                             ,delta,tazim)
            call disthead(nodelat(inode),nodelon(inode),applat,applon &
                                             ,delta,bazimnd(inode))

            xnode(iev,inode) = appcirc*delta - stadist(iev,iref(iev))
            azidiff = tazimref - tazim

            if (azidiff.gt.180.) azidiff = azidiff -360.
            if (azidiff.lt.-180.) azidiff = azidiff +360.

            ynode(iev,inode) = appcirc*sin(delta*convdeg)*azidiff
          enddo

!  similarly for outlines of region of interest
          do ibox = 1, 4
            call disthead(applat,applon,boxlat(ibox),boxlon(ibox) &
                                             ,delta,tazim)
            xbox(iev,ibox) = appcirc*delta - stadist(iev,iref(iev))
            azidiff = tazimref - tazim
        
            if (azidiff.gt.180.) azidiff = azidiff -360.
            if (azidiff.lt.-180.) azidiff = azidiff +360.
        
            ybox(iev,ibox) = appcirc*sin(delta*convdeg)*azidiff
!      write(*,*) ibox, xbox(ibox),ybox(ibox),delta, tazim
          enddo

!  find closest corner along great circle path
          xboxmin = xbox(iev,1)
          jbox(iev) = 1
        
          do ibox = 2,4
            if (xbox(iev,ibox).lt.xboxmin) then
              xboxmin = xbox(iev,ibox)
              jbox(iev) = ibox
            endif
          enddo      

!  find the mininum intersection point with corner for each allowed
!  azimuthal deviation of the two plane waves
!  - do not switch corners with azimuths as this
!  may cause discontinuous changes in needed wave parameters
! xmin here is distance from reference station of wavefront along path of
!  plane waves coming in at different angles

          do 124 ideg = 1,ndeg
            azimt = ((ideg-1.) - (ndeg-1.)/2.)*convdeg
            xmin(iev,ideg) = &
                             xbox(iev,jbox(iev))*cos(azimt) &
                             + ybox(iev,jbox(iev))*sin(azimt) 

  124     enddo

              
!  generate real and imaginary components normalized by rms amplitude
!  and compared to phase at reference station
!  Phase shift relative to reference corrected for any difference in start time
                  
          do ista = 1, nsta(iev)

            dph(iev,ista) = staph(iev,ista,ifreq)-staph(iev,iref(iev), &
                 ifreq) + freq(ifreq)*(beg(iev,ista)-beg(iev,iref(iev)))
            streal(iev,ista,ifreq)=staamp(iev,ista,ifreq) &
                  *cos(dph(iev,ista)*twopi)/amprms(iev,ifreq)
            stimag(iev,ista,ifreq)=-staamp(iev,ista,ifreq) &
                  *sin(dph(iev,ista)*twopi)/amprms(iev,ifreq)
          enddo 

        enddo ! end iev


        print *, "bp 2"
!  construct form of a priori smooothness criterion matrix for velocity 
!  parameters  (minimum curvature)

!  initialize
        do ismth = 1, nnodes
          do jsmth = 1, nnodes
            b(ismth,jsmth) = 0.0
            btb(ismth,jsmth) = 0.0
          enddo
        enddo

!  curvature terms in y-direction - but skip at edges where undefined
        do ismth = 1,ncol
          do jsmth = 2, kj-1
            ism = (ismth-1)*kj + jsmth
            b(ism,ism) = 2.0 +b(ism,ism)
            b(ism,ism-1) = -1.0+ b(ism,ism-1)
            b(ism,ism+1) = -1.0+ b(ism,ism+1)
          enddo
        enddo

!  curvature terms in x-direction
        do jsmth = 1, kj
          do ismth = 2, ncol-1
            ism = (ismth-1)*kj + jsmth 
            b(ism,ism) = 2.0 +b(ism,ism)
            b(ism,ism-kj) = -1.0 + b(ism,ism-kj)
            b(ism,ism+kj) = -1.0 + b(ism,ism+kj)
          enddo
        enddo

        do ismth = 1, nnodes
          do jsmth = 1, nnodes
            write(14,*) ismth,jsmth,b(ismth,jsmth)
          enddo
        enddo

!  now construct BtransposeB	
        do j = 1, nnodes
          do jj = 1,j
            btb(jj,j) = 0.0

            do i = 1, nnodes
              btb(jj,j)= btb(jj,j) + b(i,jj)*b(i,j)
            enddo

            btb(j,jj) = btb(jj,j)
          enddo
        enddo
      


!  *************************
!  iterate from here
!  *************************
        iter = 1
        icnt = 1

!  add one data point for constraint on amplitude multipliers that has to
!  equal number of station corrections 
        nobs = nobs + 1 
        dstacnt = jstacnt

100     continue

!  increase a priori variance for station amplitude corrections after 5
!  iterations

        if (icnt.gt.5) then
          do ii = 1, jstacnt
            ip = npnoamp+ii
            covinv(ip) = 1.0/(0.30**2)
          enddo
        endif

        sumampcor = 0.0

        do ii = 1, jstacnt
          sumampcor = sumampcor + ampmult(ii)
        enddo

        d(nobs) = (dstacnt - sumampcor)/1.0e-4     
        write(*,*) iter, icnt

!  initialize g matrix
        do irow = 1, nobs
          do icol = 1, np
            g(irow,icol) = 0.0
          enddo
        enddo

        do ii = 1,jstacnt
          icol = npnoamp+ii
          g(nobs,icol) = 1.0/1.0e-4
        enddo 


!  begin loop over events for residuals and partial derivatives

        naddat = 0
 
        do iev = 1, nevents

        !  long loop calculate the sensitivity kernel for each 
        !  possible angle between -40 to 40

        do 125 ideg = 1,ndeg
          azimt = ((ideg-1.) - (ndeg-1.)/2.)*convdeg
          cs2n = cos(2.0*(convdeg*bazi(iev,iref(iev))- azimt))
          sn2n = sin(2.0*(convdeg*bazi(iev,iref(iev))- azimt))
          cs1z = cos(azimt)
          sn1z = sin(azimt)

          do ii = 1,nnodes
            ! don't really need arrays here in this version, but 
            ! used throughout so retained
            cos2node(ii,ideg)=cs2n
            sin2node(ii,ideg)=sn2n
          enddo

          !  calculate current apparent velocity at each node

          do ii = 1, nnodes
            iii  = ii+i6
            jjj  = i6+nnodes + idnode(ii)
            jjjj = jjj + iarea
            appvel(ii) = crrntmod(iii) &
                        + cos2node(ii,ideg)*crrntmod(jjj) &
                        + sin2node(ii,ideg)*crrntmod(jjjj)
          enddo



          !  long loop over stations, generating sensitivity kernels for each

          do 102 ista = 1, nsta(iev)

            do ii = 1, nnodes
              wgtnode1(ista,ii,ideg) = 0.0
              ampwgtnode1(ista,ii,ideg) = 0.0
            enddo

            xstatemp = xsta(iev,ista)*cs1z &
                       + ysta(iev,ista)*sn1z               
     
            ystatemp = -xsta(iev,ista)*sn1z &
                       + ysta(iev,ista)*cs1z     

            do ii = 1,nnodes

              xnodetemp = xnode(iev,ii)*cs1z &
                          + ynode(iev,ii)*sn1z               
     
              ynodetemp = -xnode(iev,ii)*sn1z &
                          + ynode(iev,ii)*cs1z     


              xstanode = xnodetemp - xstatemp
              ystanode = ynodetemp - ystatemp

!  find nearest point in sensitivity kernel - kernels should have been
!  calculated with smoothing so that they represent sensitivity to nodal
!  coefficient rather than velocity at point  - could interpolate, but 
! this is sufficiently accurate if kernels on fine enough scale (i.e.,
!  dxkern << smoothing length

              if (xnodetemp.ge.xmin(iev,ideg)) then
                if (xstanode.ge.0.0) then
                  ixindex = int( xstanode/dxkern +0.5) + (nxkern+1)/2
                else
                  ixindex = int( xstanode/dxkern -0.5) + (nxkern+1)/2
                endif

                if (ystanode.ge.0.0) then
                  iyindex = int( ystanode/dykern +0.5) + (nykern+1)/2
                else

                  iyindex = int( ystanode/dykern -0.5) + (nykern+1)/2
                endif

                if (ixindex .lt.1 .or. ixindex .gt. nxkern &
                    .or. iyindex .lt.1 .or. iyindex .gt. nykern) then 


                  write(*,*) 'ixindex,iyindex,iev,ista,ideg', &
                            ixindex,iyindex,iev,ista,ideg
                  wgtnode1(ista,ii,ideg)    = 0.
                  ampwgtnode1(ista,ii,ideg) = 0.  

                else

                  wgtnode1(ista,ii,ideg) = sensitivity(ixindex,iyindex)
                  ampwgtnode1(ista,ii,ideg) =     ampsens(ixindex,iyindex)
                endif 

              else 
                wgtnode1(ista,ii,ideg)    = 0.
                ampwgtnode1(ista,ii,ideg) = 0.  
              endif      
            enddo


            dphase(ista,ideg) = 0.
            dampper(ista,ideg) = 0.

!  corrections for second order effects of large velocity changes in version 587
!   old version commented out
            do inode = 1, nnodes
!          dphase(ista,ideg) = dphase(ista,ideg) &
!        + (1.0/twopi)*wgtnode1(ista,inode,ideg) &
!                 *(appvel(inode)-unifvel)/unifvel
              dphase(ista,ideg) = dphase(ista,ideg) &
                               + (1.0/twopi)*wgtnode1(ista,inode,ideg) &
                               *(appvel(inode)-unifvel)/appvel(inode) &
                               /(appvel(inode)/unifvel)

            enddo


            do inode =1, nnodes
              dampper(ista,ideg) = dampper(ista,ideg) &
                                   +ampwgtnode1(ista,inode,ideg)* &
!                                  (appvel(inode)-unifvel)/unifvel &
                                   (appvel(inode)-unifvel)/appvel(inode) &
                                   /(appvel(inode)/unifvel)
            enddo

!  end loop over stations
 102      enddo  ! ista

!  end loop over angles
 125    enddo
!         write(*,*) iev, dphase(5,5),dampper(5,5)

!  find best fitting wave parameters from grid search
        call search(pb)

!      write(*,*)iev,pb(3)/convdeg,pb(4)/convdeg,pb(1),pb(2),pb(5),pb(6)

!  use best model as event starting model for linearized inversion
        ip = (iev-1)*6
        startamp1(iev) = pb(1)
        startamp2(iev) = pb(2)
        stazim1(iev) = pb(3)
        stazim2(iev) = pb(4)
        stphase1(iev) = pb(5)
        stphase2(iev) = pb(6)
        origmod(1+ip) = startamp1(iev)
        origmod(2+ip) = startamp2(iev)
        origmod(3+ip) = stazim1(iev)
        origmod(4+ip) = stazim2(iev)
        origmod(5+ip) = stphase1(iev)
        origmod(6+ip) = stphase2(iev)
        crrntmod(1+ip) = startamp1(iev)
        crrntmod(2+ip) = startamp2(iev)
        crrntmod(3+ip) = stazim1(iev)
        crrntmod(4+ip) = stazim2(iev)
        crrntmod(5+ip) = stphase1(iev)
        crrntmod(6+ip) = stphase2(iev)


!  calculate partial derivatives for wave parameters
        ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2+0.01) + 1
        ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2+0.01) + 1

!  calculate effect of uniform velocity for first plane wave
        do 10 ista = 1, nsta(iev)
          xstatemp = xsta(iev,ista)*cos(stazim1(iev)) &
                     + ysta(iev,ista)*sin(stazim1(iev))               
     
!       ystatemp = 
!     1  -xsta(iev,ista)*sin(stazim1(iev))
!     1 + ysta(iev,ista)*cos(stazim1(iev))     

          dtime1(ista) = (xstatemp-xmin(iev,ideg1))/unifvel

  10    enddo


!       calculate for second plane wave

        do 20 ista = 1, nsta(iev)
         
          xstatemp = xsta(iev,ista)*cos(stazim2(iev)) &
                   + ysta(iev,ista)*sin(stazim2(iev))               
     
!       ystatemp = &
!        -xsta(iev,ista)*sin(stazim2(iev)) &
!       + ysta(iev,ista)*cos(stazim2(iev))     

          dtime2(ista) = (xstatemp-xmin(iev,ideg2))/unifvel

  20    enddo
   
!  convert time to phase and assign effects of lateral heterogeneities

        do ista = 1, nsta(iev)

          phase1(iev,ista) = dtime1(ista)*freq(1)
          phase2(iev,ista) = dtime2(ista)*freq(1)
          dphase1(iev,ista) =   dphase(ista,ideg1)
          dphase2(iev,ista) =   dphase(ista,ideg2)  
          damp1per(iev,ista) = dampper(ista,ideg1)
          damp2per(iev,ista) = dampper(ista,ideg2) 
        enddo


! calculate the change of ph1,ph2, damp1,damp2 with respect to azimuth at the reference station 

        aziminc1 = 2.*convdeg
        aziminc2 = 2.*convdeg
      
        if ((stazim1(iev)+aziminc1).gt.40.*convdeg) then
          aziminc1 = -aziminc1
        endif

        if ((stazim2(iev)+aziminc2).gt.40.*convdeg) then
          aziminc2 = -aziminc2
        endif
        
!	stazim1(iev) = stazim1(iev) + aziminc1
!	stazim2(iev) = stazim2(iev) + aziminc2
!
!       ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2+0.01) + 1
!       ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2+0.01) + 1
        staziminc1 = stazim1(iev) + aziminc1
        staziminc2 = stazim2(iev) + aziminc2

        ideg1 = int(staziminc1/convdeg+(ndeg-1)/2+0.01) + 1
        ideg2 = int(staziminc2/convdeg+(ndeg-1)/2+0.01) + 1


        xstatemp = xsta(iev,iref(iev))*cos(staziminc1) &
                 + ysta(iev,iref(iev))*sin(staziminc1)               
     
!       ystatemp = &
!        -xsta(iev,iref(iev))*sin(staziminc1) &
!       + ysta(iev,iref(iev))*cos(staziminc1)     

   
!  calculate delay due to uniform velocity


        dtime1ref = (xstatemp-xmin(iev,ideg1))/unifvel

!       calculate for second plane wave

        xstatemp = xsta(iev,iref(iev))*cos(staziminc2) &
                 + ysta(iev,iref(iev))*sin(staziminc2)               
     
!       ystatemp = &
!        -xsta(iev,iref(iev))*sin(staziminc2) &
!       + ysta(iev,iref(iev))*cos(staziminc2)     
      
        dtime2ref = (xstatemp-xmin(iev,ideg2))/unifvel

!    convert time to phase and assign effects of lateral heterogeneities

        phase1ref = dtime1ref*freq(1)
        phase2ref = dtime2ref*freq(1)
        dphase1ref =   dphase(iref(iev),ideg1)
        dphase2ref =   dphase(iref(iev),ideg2) 
        damp1perref = dampper(iref(iev),ideg1)
        damp2perref = dampper(iref(iev),ideg2)

!       stazim1(iev) = stazim1(iev) -aziminc1
!       stazim2(iev) = stazim2(iev) - aziminc2

!cccccccccccccccccccccccc end of for refrence station ccccccccccccccccccccccccc

        do ista = 1,nsta(iev)
!          stazim1(iev) = stazim1(iev) + aziminc1
!	   stazim2(iev) = stazim2(iev) + aziminc2
      
!          ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2) + 1
!          ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2) + 1	  

          xstatemp = xsta(iev,ista)*cos(staziminc1) &
                   + ysta(iev,ista)*sin(staziminc1)               
     
!       ystatemp = &
!        -xsta(iev,ista)*sin(staziminc1) &
!       + ysta(iev,ista)*cos(staziminc1)     

          dtime1temp = (xstatemp-xmin(iev,ideg1))/unifvel

!       calculate for second plane wave
          xstatemp = xsta(iev,ista)*cos(staziminc2) &
                   + ysta(iev,ista)*sin(staziminc2)               
     
!       ystatemp = &
!        -xsta(iev,ista)*sin(staziminc2(iev)) &
!       + ysta(iev,ista)*cos(staziminc2(iev))     

          dtime2temp = (xstatemp-xmin(iev,ideg2))/unifvel

          phase1temp = dtime1temp*freq(1)
          phase2temp = dtime2temp*freq(1)
          dphase1temp = dphase(ista,ideg1)
          dphase2temp = dphase(ista,ideg2)
          damp1pertemp = dampper(ista,ideg1)
          damp2pertemp = dampper(ista,ideg2) 
  
!	stazim1(iev) = stazim1(iev) - aziminc1
!	stazim2(iev) = stazim2(iev) - aziminc2


!  calculate partial differences with respect to azimuth by finite differences

          parph1azim = ( (phase1temp+dphase1temp)- &
                         (phase1(iev,ista)+dphase1(iev,ista)))/aziminc1 &
                         - ( (phase1ref+dphase1ref) - &
                         (phase1(iev,iref(iev))+dphase1(iev,iref(iev))))/aziminc1

          parph2azim = ( (phase2temp+dphase2temp)- &
                         (phase2(iev,ista)+dphase2(iev,ista)))/aziminc2 &
                         - ( (phase2ref+dphase2ref) - &
                         (phase2(iev,iref(iev))+dphase2(iev,iref(iev))))/aziminc2

          pardamp1azim = (damp1pertemp - damp1per(iev,ista))/aziminc1
          pardamp2azim = (damp2pertemp - damp2per(iev,ista))/aziminc2

! #################################################################

          ideg1 = int(stazim1(iev)/convdeg+(ndeg-1)/2 +0.01) + 1
          ideg2 = int(stazim2(iev)/convdeg+(ndeg-1)/2 +0.01) + 1




!########## end of parph1azim,parph2azim,pardamp1azim,pardamp2azim######## 

!        write(*,*) phase1(iev,ista),phase2(iev,ista)
          prefase1 = ( (phase1(iev,ista)+dphase1(iev,ista)) &
                     - (phase1(iev,iref(iev))+dphase1(iev,iref(iev)))) &
                     + stphase1(iev)
          prefase2 = ( (phase2(iev,ista)+dphase2(iev,ista)) &
                     - (phase2(iev,iref(iev))+dphase2(iev,iref(iev)))) &
                     + stphase2(iev)

          staamp1 = startamp1(iev)*(1.+damp1per(iev,ista))
          staamp2 = startamp2(iev)*(1.+damp2per(iev,ista))

          parph1cor1=0.                                        
          parph2cor1=0.                                        
        
          itypst2=0

          if (ntyp1.ne.0) then 
            do itypst=1, ntyp1                                    
              if (istanum(iev,ista).eq.ityp1sta(itypst)) then     
                prefase1=prefase1+phcor(itypst)                    
                prefase2=prefase2+phcor(itypst)                    
                parph1cor1=1.                                      
                parph2cor1=1.                                      
                itypst2=itypst                                     
              endif                                               
            enddo                                                
          endif

          cosph1 = cos(prefase1*twopi)
          cosph2 = cos(prefase2*twopi)
          sinph1 = sin(prefase1*twopi)
          sinph2 = sin(prefase2*twopi)
          
          prereal = staamp1*cosph1+staamp2*cosph2
!	  write(*,*) prereal  
          preimag = -1.0*(staamp1*sinph1+ staamp2*sinph2)
!           write(*,*) ampmult(istavar(istanum(iev,ista)))
          prereal = prereal*ampmult(istavar(istanum(iev,ista))) &
                    *exp(-gamma*xsta(iev,ista))
!          write(*,*) prereal
          preimag = preimag*ampmult(istavar(istanum(iev,ista))) &
                    *exp(-gamma*xsta(iev,ista))




!  data vector and partial derivatives listed event by event with all
!  real data for first event followed by imaginary data, then onto next event
!  d contains misfit to starting model


          kreal = ista + naddat
          kimag = kreal + nsta(iev)
          d(kreal) = (streal(iev,ista,ifreq) - prereal)/stddevdata(iev)
          d(kimag) = (stimag(iev,ista,ifreq) - preimag)/stddevdata(iev)


!  partial derivatives for station amplitude correction factors and attenuation
          g(kreal,npnoamp+istavar(istanum(iev,ista))) = &
            prereal/ampmult(istavar(istanum(iev,ista)))/stddevdata(iev)

          g(kimag,npnoamp+istavar(istanum(iev,ista))) = &
            preimag/ampmult(istavar(istanum(iev,ista)))/stddevdata(iev)
          g(kreal,np) = - xsta(iev,ista)*prereal/stddevdata(iev)
          g(kimag,np) = - xsta(iev,ista)*preimag/stddevdata(iev)

          ampadj = ampmult(istavar(istanum(iev,ista))) &
            *exp(-gamma*xsta(iev,ista))
        
          atte1 = 1.
          atte2 = 1.

!  partial derivatives for station phase correction  
          if (itypst2.gt.0) then                                      
            g(kreal,npp+itypst2) = parph1cor1*(-staamp1*sinph1*twopi &
                                   -staamp2*sinph2*twopi)/stddevdata(iev)*ampadj
            g(kimag,npp+itypst2) = parph1cor1*(-staamp1*cosph1*twopi &
                                   -staamp2*cosph2*twopi)/stddevdata(iev)*ampadj
          endif

          do ii = 1, nnodes

            parph1v = (1.0/twopi)*wgtnode1(ista,ii,ideg1)/unifvel &
               - (1.0/twopi)*wgtnode1(iref(iev),ii,ideg1)/unifvel

            parph2v =  (1.0/twopi)*wgtnode1(ista,ii,ideg2)/unifvel &
                - (1.0/twopi)*wgtnode1(iref(iev),ii,ideg2)/unifvel

            parph1cs = (1.0/twopi)*cos2node(ii,ideg1)* &
                       wgtnode1(ista,ii,ideg1)/unifvel &
                       - (1.0/twopi)*cos2node(ii,ideg1)* &
                       wgtnode1(iref(iev),ii,ideg1)/unifvel

            parph2cs = (1.0/twopi)*cos2node(ii,ideg2)* &
                       wgtnode1(ista,ii,ideg2)/unifvel &
                       - (1.0/twopi)*cos2node(ii,ideg2)* &
                               wgtnode1(iref(iev),ii,ideg2)/unifvel

            parph1sn = (1.0/twopi)*sin2node(ii,ideg1)* &
                       wgtnode1(ista,ii,ideg1)/unifvel &
                       - (1.0/twopi)*sin2node(ii,ideg1)* &
                       wgtnode1(iref(iev),ii,ideg1)/unifvel

            parph2sn = (1.0/twopi)*sin2node(ii,ideg2)* &
                       wgtnode1(ista,ii,ideg2)/unifvel &
                       - (1.0/twopi)*sin2node(ii,ideg2)* &
                       wgtnode1(iref(iev),ii,ideg2)/unifvel



            paramp1v = startamp1(iev)*ampwgtnode1(ista,ii,ideg1)/unifvel
            paramp2v = startamp2(iev)*ampwgtnode1(ista,ii,ideg2)/unifvel


            paramp1cs  = startamp1(iev)*cos2node(ii,ideg1)* &
                                  ampwgtnode1(ista,ii,ideg1)/unifvel

            paramp2cs  = startamp2(iev)*cos2node(ii,ideg2)* &
                                  ampwgtnode1(ista,ii,ideg2)/unifvel

            paramp1sn  = startamp1(iev)*sin2node(ii,ideg1)* &
                                  ampwgtnode1(ista,ii,ideg1)/unifvel

            paramp2sn  = startamp2(iev)*sin2node(ii,ideg2)* &
                                  ampwgtnode1(ista,ii,ideg2)/unifvel


!  partial derivatives with respect to velocity, & cos2theta & sin2theta

            jjjarea = i6 + nnodes + idnode(ii)
            jjjjarea = jjjarea + iarea

            g(kreal,i6+ii)= (-startamp1(iev)*(1.0+damp1per(iev,ista))* &
                            twopi*sinph1*parph1v*atte1 &
                            -startamp2(iev)*(1.0+damp2per(iev,ista))* &
                            twopi*sinph2*parph2v*atte2 &
                            + paramp1v*cosph1*atte1 + paramp2v*cosph2*atte2 )  &
                            *ampadj/stddevdata(iev) 

            g(kimag,i6+ii)= (-startamp1(iev)*(1.0+damp1per(iev,ista))* &
                            twopi*cosph1*parph1v*atte1 &
                            -startamp2(iev)*(1.0+damp2per(iev,ista))* &
                            twopi*cosph2*parph2v*atte2 &
                            - (paramp1v*sinph1*atte1  +  paramp2v*sinph2*atte2 ) ) &
                            *ampadj/stddevdata(iev) 


            g(kreal,jjjarea)= (-startamp1(iev)*(1.0+damp1per(iev,ista)) &
                              *twopi*sinph1*parph1cs*atte1 &
                              -startamp2(iev)*(1.0+damp2per(iev,ista)) &
                              *twopi*sinph2*parph2cs*atte2 &
                              + paramp1cs*cosph1*atte1 + paramp2cs*cosph2*atte2 ) &
                              *ampadj/stddevdata(iev)+g(kreal,jjjarea)

            g(kimag,jjjarea)= (-startamp1(iev)*(1.0+damp1per(iev,ista)) &
                              *twopi*cosph1*parph1cs*atte1 &
                              -startamp2(iev)*(1.0+damp2per(iev,ista)) &
                              *twopi*cosph2*parph2cs*atte2 &
                              - (paramp1cs*sinph1*atte1  +  paramp2cs*sinph2*atte2 )) &
                              *ampadj/stddevdata(iev)+g(kimag,jjjarea)

            g(kreal,jjjjarea)= (-startamp1(iev)*(1.0+damp1per(iev,ista)) &
                               *twopi*sinph1*parph1sn*atte1 &
                               -startamp2(iev)*(1.0+damp2per(iev,ista)) &
                               *twopi*sinph2*parph2sn*atte2 &
                               + paramp1sn*cosph1*atte1 + paramp2sn*cosph2*atte2 ) &
                               *ampadj/stddevdata(iev)+g(kreal,jjjjarea)

            g(kimag,jjjjarea)= (-startamp1(iev)*(1.0+damp1per(iev,ista)) &
                               *twopi*cosph1*parph1sn*atte1 &
                               -startamp2(iev)*(1.0+damp2per(iev,ista)) &
                               *twopi*cosph2*parph2sn*atte2 &
                               - (paramp1sn*sinph1*atte1  +  paramp2sn*sinph2*atte2 ) ) &
                               *ampadj/stddevdata(iev)+g(kimag,jjjjarea)

          enddo

!  partial derivatives in order are with respect to amplitudes,
!  azimuths, starting phases and  slowness
          ip = (iev-1)*6
          g(kreal,1+ip) = (1.0+damp1per(iev,ista))*cosph1*atte1 &
                           *ampadj/stddevdata(iev)

          g(kreal,2+ip) = (1.0+damp2per(iev,ista))*cosph2*atte2 &
                           *ampadj/stddevdata(iev)

          g(kreal,3+ip) = (-startamp1(iev)*(1.0+damp1per(iev,ista))*sinph1 &
                            *parph1azim*twopi &
                            + startamp1(iev)*pardamp1azim*cosph1)*atte1 &
                            *ampadj/stddevdata(iev)

          g(kreal,4+ip) = (-startamp2(iev)*(1.0+damp2per(iev,ista))*sinph2 &
                            *parph2azim*twopi &
                            + startamp2(iev)*pardamp2azim*cosph2)*atte2  &
                            *ampadj/stddevdata(iev)

          g(kreal,5+ip) = -startamp1(iev)*(1.0+damp1per(iev,ista))*sinph1 &
                            *twopi*atte1 &
                            *ampadj/stddevdata(iev)

          g(kreal,6+ip) = -startamp2(iev)*(1.0+damp2per(iev,ista))*sinph2 &
                            *twopi*atte2 &
                            *ampadj/stddevdata(iev)

          g(kimag,1+ip) = -(1.0+damp1per(iev,ista))*sinph1*atte1 &
                            *ampadj/stddevdata(iev)

          g(kimag,2+ip) = -(1.0+damp2per(iev,ista))*sinph2*atte2 &
                            *ampadj/stddevdata(iev)

          g(kimag,3+ip) = -(startamp1(iev)*(1.0+damp1per(iev,ista))*cosph1 &
                            *parph1azim*twopi &
                            + startamp1(iev)*pardamp1azim*sinph1)*atte1 &
                            *ampadj/stddevdata(iev)

          g(kimag,4+ip) = -(startamp2(iev)*(1.0+damp2per(iev,ista))*cosph2 &
                            *parph2azim*twopi &
                            + startamp2(iev)*pardamp2azim*sinph2)*atte2 &
                            *ampadj/stddevdata(iev)

          g(kimag,5+ip) = -startamp1(iev)*(1.0+damp1per(iev,ista))*cosph1 &
                            *twopi*atte1 &
                            *ampadj/stddevdata(iev)

          g(kimag,6+ip) = -startamp2(iev)*(1.0+damp2per(iev,ista))*cosph2 &
                             *twopi*atte2 &
                             *ampadj/stddevdata(iev)
!  normalize by stddevdata
        enddo  ! ista

        naddat = naddat + 2*nsta(iev)
      enddo

!  calculate misfit before inversion
	premisfit = 0.0
	do i = 1, nobs
	  premisfit = premisfit + d(i)*d(i)
	enddo
	premisrms = sqrt(premisfit/nobs)
	write (*,*) "premisrms ", premisrms

!        write(*,*) "end loop over events for partial derivatives"

!  Calculate gtg and gtd	
        do j = 1, np
          gtd(j) = 0.0
          do i = 1, nobs
            gtd(j) = gtd(j) + g(i,j)*d(i)
          enddo
!   add to gtd Tarantola term penalizing misfit to original starting model
!   but skip for wave parameters

          if (j.le.i6) gtdcmm(j) = gtd(j)       
          if (j.gt.i6) gtdcmm(j) = gtd(j) -  &
                              covinv(j)*(crrntmod(j)-origmod(j))
!          gtdcmm(j) = gtd(j) - covinv(j)*(crrntmod(j)-origmod(j))
!  construct gtg  
          do jj = 1,j
            gtg(jj,j) = 0.0
            do i = 1, nobs
              gtg(jj,j)= gtg(jj,j) + g(i,jj)*g(i,j)
            enddo
            gtg(j,jj) = gtg(jj,j)
            savegtg(j,jj) = gtg(j,jj)
            savegtg(jj,j) = gtg(jj,j)
          enddo
! increase damping for wave parameters
          gtg(j,j) = gtg(j,j) + covinv(j)
        enddo
!  ******************************
!   add smoothness constraint (minimum curvature).  Find coefficient, smthco, 
!  that minimizes, in least squares sense, the off-diagonal terms of gtg,
!   which should lead to minimal off-diagonal terms of covariance matrix
!
        sumsmth2 = 0.0
	sumsmgg = 0.0
	do ismth = 2,nnodes
	  do jsmth = 1, ismth-1
	    sumsmth2 = sumsmth2 + btb(ismth,jsmth)**2
	    sumsmgg = sumsmgg + btb(ismth,jsmth)* &
                               gtg(ismth+i6,jsmth+i6)
          enddo
	enddo
	smthco = -sumsmgg/sumsmth2/divfac
	write(*,*) 'smthco', smthco
	do ismth = 1,nnodes
	  do jsmth = 1,nnodes
	    gtg(ismth+i6,jsmth+i6)= gtg(ismth+i6,jsmth+i6) + &
                      btb(ismth,jsmth)*smthco
          enddo
	enddo
!  also add smoothness constraint to Tarantola term penalizing curvature in
!  solution
        do ismth = 1,nnodes
	  do jsmth = 1, nnodes
	    gtdcmm(ismth+i6) = gtdcmm(ismth+i6) -  &
          smthco*btb(ismth,jsmth)*(crrntmod(jsmth+i6)-origmod(jsmth+i6))
          enddo
	enddo
	    
!  ******************************
!  Now add smoothness constraint (minimum curvature) for cos2theta terms.  
!  Find coefficient, smthco1, 
!  that minimizes, in least squares sense, the off-diagonal terms of gtg,
!   which should lead to minimal off-diagonal terms of covariance matrix
!
        sumsmth2 = 0.0
	sumsmgg = 0.0
	i6pnod = i6 + nnodes
	do ismth = 2,nnodes
	  do jsmth = 1, ismth-1
	    sumsmth2 = sumsmth2 + btb(ismth,jsmth)**2
	    sumsmgg = sumsmgg + btb(ismth,jsmth)* &
                       gtg(ismth+i6pnod,jsmth+i6pnod)
          enddo
	enddo
	smthco1 = -sumsmgg/sumsmth2
	write(*,*) 'smthco1 ', smthco1
	do ismth = 1,nnodes
	  do jsmth = 1,nnodes
	    gtg(ismth+i6pnod,jsmth+i6pnod)= gtg(ismth+i6pnod,jsmth+i6pnod) + &
                      btb(ismth,jsmth)*smthco1
          enddo
	enddo
!  also add smoothness constraint to Tarantola term penalizing curvature in
!  solution
        do ismth = 1,nnodes
	  do jsmth = 1, nnodes
	    gtdcmm(ismth+i6pnod) = gtdcmm(ismth+i6pnod) -  &
             smthco1*btb(ismth,jsmth)* &
             (crrntmod(jsmth+i6pnod)-origmod(jsmth+i6pnod))
          enddo
	enddo
	print *, "bp 8"
!  ******************************
!  Now add smoothness constraint (minimum curvature) for sin2theta terms.  
!  Find coefficient, smthco2, 
!  that minimizes, in least squares sense, the off-diagonal terms of gtg,
!   which should lead to minimal off-diagonal terms of covariance matrix
!
        sumsmth2 = 0.0
	sumsmgg = 0.0
	i6pnod = i6 + 2*nnodes
	do ismth = 2,nnodes
	  do jsmth = 1, ismth-1
	    sumsmth2 = sumsmth2 + btb(ismth,jsmth)**2
	    sumsmgg = sumsmgg + btb(ismth,jsmth)* &
                  gtg(ismth+i6pnod,jsmth+i6pnod)
          enddo
	enddo
	smthco2 = -sumsmgg/sumsmth2
	write(*,*) 'smthco2 ', smthco2
	do ismth = 1,nnodes
	  do jsmth = 1,nnodes
	    gtg(ismth+i6pnod,jsmth+i6pnod)= gtg(ismth+i6pnod,jsmth+i6pnod) + &
                      btb(ismth,jsmth)*smthco2
          enddo
	enddo
!  also add smoothness constraint to Tarantola term penalizing curvature in
!  solution
        do ismth = 1,nnodes
	  do jsmth = 1, nnodes
	    gtdcmm(ismth+i6pnod) = gtdcmm(ismth+i6pnod) -  &
             smthco2*btb(ismth,jsmth)* &
             (crrntmod(jsmth+i6pnod)-origmod(jsmth+i6pnod))
          enddo
	enddo
	    
	    

!  ***********************
!  now eliminate wave parameters from inversion process
!  too many parameters and too unstable when small number of observations
!  per event
!  *************************
!      npsub = np - i6
!	do j = 1, npsub
!	  do jj = 1, j
!	    subgtg(jj,j) = gtg(jj+i6,j+i6)
!	    subgtg(j,jj) = subgtg(jj,j)
!            subsavegtg(j,jj) = savegtg(j+i6,jj+i6) 
!            subsavegtg(jj,j) = savegtg(jj+i6,j+i6)
!	  enddo
!	  subgtd(j) = gtd(j+i6)
!	  subgtdcmm(j) = gtdcmm(j+i6)
!	enddo
   

!  Invert gtg.  gtg will be destroyed.  
!  Not the most efficient approach because doesn't take advantage of 
!  symmetry of gtg.  Use LU decomposition from Press et al.
        do i= 1,np
          do j = 1, np
            gtginv(i,j)= 0.0D0
          enddo
          gtginv(i,i) =1.0D0
        enddo
        call dludcmp(gtg,np,nparam,indx,ddd)
        do j = 1,np
          call dlubksb(gtg,np,nparam,indx,gtginv(1,j))
        enddo
!	do i= 1,npsub
!          do j = 1, npsub
!            subgtginv(i,j)= 0.0D0
!          enddo
!          subgtginv(i,i) =1.0D0
!        enddo
!        call dludcmp(subgtg,npsub,nparam,indx,ddd)
!        do j = 1,np
!          call dlubksb(subgtg,npsub,nparam,indx,subgtginv(1,j))
!        enddo

!  Find change to starting model
        do i= 1, np
          change(i)=0.0
          do j = 1,np
            change(i) =change(i) + gtdcmm(j)*gtginv(i,j)
          enddo
        enddo
!       do i= 1, npsub
!          change(i+i6)=0.0
!          do j = 1,npsub
!            change(i+i6) =change(i+i6) + subgtdcmm(j)*subgtginv(i,j)
!          enddo
!        enddo

!  Find rank (sum of diagonals of resolution matrix), i.e., number of
!  pieces of information or number of independent model parameters
!  rank1 is contribution from source wave terms, rank2 from velocity
!  variables
        rank1 = 0.0
        rank2 = 0.0
	rank3 = 0.0
	rank4 = 0.0
	write(14,*) icnt,iter
	do i=1,np
	  write(14,*) i, crrntmod(i),change(i)
          resparm =0.0
          do j = 1,np
            resparm = resparm + gtginv(i,j)*savegtg(j,i)
          enddo
          if (i.le.i6) then
            rank1 = rank1 + resparm
	  endif
          if ((i.gt.i6).and.(i.le.i6+nnodes)) then
            rank2 = rank2 + resparm
!   Remember diagonal resolution values for isotropic velocity parameters
	    residdiag(i-i6) = resparm
	  endif
          if ((i.gt.i6+nnodes).and.(i.le.npnoamp)) then
            rank4 = rank4 + resparm
          endif
	  if (i.gt.npnoamp) then
	    rank3 = rank3 + resparm
	  endif
        enddo 
        rank = rank1 + rank2 +rank3 + rank4   
!************************************
!	do i=1,npsub
!          resparm =0.0
!          do j = 1,npsub
!            resparm = resparm + subgtginv(i,j)*subsavegtg(j,i)
!          enddo
!          if (i.le.i6) then
!            rank1 = rank1 + resparm
!          else
!          if (i.le.(npnoamp-i6)) then
!            rank2 = rank2 + resparm
!	    if (i.le.nnodes) then
!	      residdiag(i) = resparm
!	    endif
!	   else
!	     rank3 = rank3 + resparm
!	   endif
	   
!          endif
!        enddo 
!        rank = rank2 + rank3     

!  Update current model

 
        naddat = 0
        do iev = 1,nevents
          ip = (iev-1)*6
          startamp1(iev) = startamp1(iev) + change(1+ip)
          startamp2(iev) = startamp2(iev) + change(2+ip)
          stazim1(iev) = stazim1(iev) + change(3+ip)
          stazim2(iev) = stazim2(iev) + change(4+ip)

!       if(stazim1(iev) .lt. -40.*convdeg) stazim1(iev) = -40.*convdeg
!       if(stazim1(iev) .gt. 40.*convdeg)  stazim1(iev) = 40.*convdeg
!       if(stazim2(iev) .lt. -40.*convdeg) stazim2(iev) = -40.*convdeg
!       if(stazim2(iev) .gt. 40.*convdeg)  stazim2(iev) = 40.*convdeg

          stphase1(iev) = stphase1(iev) + change(5+ip)
          stphase2(iev) = stphase2(iev) + change(6+ip)
          crrntmod(1+ip) = startamp1(iev)
          crrntmod(2+ip) = startamp2(iev)
          crrntmod(3+ip) = stazim1(iev)
          crrntmod(4+ip) = stazim2(iev)
          crrntmod(5+ip) = stphase1(iev)
          crrntmod(6+ip) = stphase2(iev)

        sumsq = 0.0
        sumsqph = 0.0

          do ista = 1,nsta(iev)
            dresid1=d(ista+naddat)*stddevdata(iev)
            dresid2=d(ista+naddat+nsta(iev))*stddevdata(iev)
            predamp = amprms(iev,ifreq) &
               *sqrt((streal(iev,ista,ifreq)-dresid1)**2  &
                       + (stimag(iev,ista,ifreq)-dresid2)**2)
            prefase = atan2(-(stimag(iev,ista,ifreq)-dresid2), &
                       (streal(iev,ista,ifreq)-dresid1))/twopi
! check for prefase being off by multiple of twopi
           absmis = abs(prefase - dph(iev,ista))
	   phasemis = absmis - int(absmis)	   
           if (phasemis.gt.0.5) phasemis = phasemis-1.0

           sumsq = sumsq+(d(ista+naddat)**2+d(ista+naddat+nsta(iev))**2) &
                      *stddevdata(iev)**2
           sumsqph = sumsqph + (dph(iev,ista)-prefase)**2
          enddo

             rmsph = sqrt(sumsqph/nsta(iev))/freq(ifreq)
	     
              if ((2*nsta(iev) - rank2/nevents - 6).lt.1.0) then
!              if ((2*nsta(iev) - rank/nevents - 6).lt.1.0) then
	      sigma2=sqrt(sumsq/1.0)
	      else
	      sigma2 = sqrt(sumsq/(2*nsta(iev) - rank2/nevents - 6))
!	      sigma2 = sqrt(sumsq/(2*nsta(iev) - rank/nevents - 6))
	      endif



          naddat = naddat + 2*nsta(iev)

	write(10,*) 'iev', idnum(iev),'  sigma2 ',sigma2,  'rmsph ',rmsph
	write(10,*) stazim1(iev)/convdeg,startamp1(iev),stphase1(iev)
	write(10,*) stazim2(iev)/convdeg,startamp2(iev),stphase2(iev)

        enddo

!  update node coefficients and calculate average isotropic velocity
	avgvel = 0.0
        do ii = 1, nnodes
          crrntmod(ii+i6) = crrntmod(ii+i6) + change(ii+i6)
	  avgvel = avgvel + crrntmod(ii+i6)
        enddo
	avgvel = avgvel/nnodes
	write(*,*) "avgvel ",avgvel
        do ii = 1, iarea
          iii = ii+i6 + nnodes
          iiii = iii+iarea
          crrntmod(iii) = crrntmod(iii) + change(iii)
          crrntmod(iiii) = crrntmod(iiii) + change(iiii)
        enddo

        do ii = 1, jstacnt
          ip = ii + npnoamp
          ampmult(ii) = ampmult(ii) + change(ip)
          crrntmod(ip) = ampmult(ii)
        enddo

!      changed for no aniso
!        do ii = 1, iages
!          iii = ii+i6 + nnodes
!          iiii = iii+iages
!          crrntmod(iii) = crrntmod(iii) + change(iii)
!          crrntmod(iiii) = crrntmod(iiii) + change(iiii)
!        enddo

         if (ntyp1.ne.0) then
         do ii=1,ntyp1                                          
           ip=npp+ii                                            
           phcor(ii)=phcor(ii)+change(ip)                       
           crrntmod(ip)=phcor(ii)                               
           write(14,*)"phcor ",ii, ityp1sta(ii), phcor(ii)      
           write(*,*)"phcor ",ii, ityp1sta(ii), phcor(ii)      
        enddo                                                   
        endif
	
        gamma = gamma + change(np)
        crrntmod(np) = gamma
	
!  Calculate predicted effect of change in model parameters and compare to actual effect
        predmis = 0.0
        do ic = 1, nobs
	  gchange(ic) = 0.0
	  do jc = i6+1,np
	    gchange(ic) = gchange(ic) + g(ic,jc)*change(jc)
	  enddo
	  predmis = predmis + (d(ic)-gchange(ic))**2
	enddo
	predmisrms = sqrt(predmis/nobs)
	write (*,*) "predmisrms ",predmisrms
!	call actlmsft(actmisrms)
!  subroutine setup for integer values of azimuths only - don't call - just for debugging
!   anyway to check for non-linearity		
!	write(*,*) "predmisrms ",predmisrms," actmisrms ",actmisrms

!   Find maximum change to model (well, to velocity
!   parameters anyway) to test for convergence 
        chmax = 0.0
        do ii = 1,nnodes
          ip = i6 + ii
          chmax = dmax1(chmax, dabs(change(ip)))
        enddo
        write(*,*) icnt, chmax
        icnt = icnt +1
!  calculate misfit after inversion
!  *******************************************
!  test for convergence and finish iteration
!  ********************************************
        if ((chmax.gt.0.0005).and.(icnt.le.iterlimit)) go to 100
        write(*,*) icnt
!  The first time through, solutions are damped with data variance as an
!  estimate.  Go through again with a posteriori estimate of data variance
!  different for each event.
        if (iter.eq.1) then
           icnttot = icnt
           iter = 2
           minstd =.03
           naddat = 0
           do iev = 1, nevents
             sumsq = 0.0
             do ista = 1, 2*nsta(iev)
               sumsq = sumsq + d(ista+naddat)**2
             enddo
!  number of degrees of freedom based on amplitude misfit - only two
!  adjustable amplitude parameters per event
!             sigma2 = sumsq/(nsta(iev)-2)
!  number of degrees of freedom assuming 6 wave parameters per event and that
!  each event supplies average of about 2 pieces of information about the 
!  velocity structure
!             sigma2 = sumsq/(2*nsta(iev)-8)
!  rank2 is the rank for velocity. the for 2-D isotropy (n nodes)and uniform anisotropy 
!  the rank is not 2*n+2 because the velocity is not independent.do the if test in our case 
	     
	      if ((2*nsta(iev) - (rank2+rank3)/nevents - 6).lt.1.0) then
	      sigma2=sumsq/1.0
	      else
	      sigma2 = sumsq/(2*nsta(iev) - (rank2+rank3)/nevents - 6)
	      endif
             
	      
             stddevdata(iev) = stddevdata(iev)*sqrt(sigma2)
!  set floor on how small stddevdata is allowed to get - below .03 probably
!  unrealistic even for best data
             if (stddevdata(iev).lt.minstd) stddevdata(iev)=minstd
             naddat = naddat + 2*nsta(iev)
           enddo
           icnt = 1
           go to 100
         endif

!  By the time this converges to a solution, d should contain the normalized 
!  misfits
!  Calculate sigma**2 and standard deviation of model parameters
!  Problem with this again with too many parameters - calculate rms instead
!  or use sigma**2 = 1
        sumsq = 0.0
        naddat = 0
        totsumsqph  = 0.0
	totsumsqamp = 0.0
        do iev = 1, nevents
          write(10,*) 'event ', idnum(iev)
          sumsqph  = 0.0
	  sumsqamp =  0.
	  sumsqtemp = 0.0
          do ista = 1,nsta(iev)
            dresid1=d(ista+naddat)*stddevdata(iev)
            dresid2=d(ista+naddat+nsta(iev))*stddevdata(iev)
!  convert misfits back to original amplitude and phase
            predamp = amprms(iev,ifreq) &
               *sqrt((streal(iev,ista,ifreq)-dresid1)**2  &
                       + (stimag(iev,ista,ifreq)-dresid2)**2)
            prefase = atan2(-(stimag(iev,ista,ifreq)-dresid2), &
                       (streal(iev,ista,ifreq)-dresid1))/twopi
! check for prefase being off by twopi
           absmis = abs(prefase - dph(iev,ista))
	   phasemis = absmis - int(absmis)	   
           if (phasemis.gt.0.5) phasemis = phasemis-1.0
           write(10,*) ista, streal(iev,ista,ifreq),dresid1, &
                           stimag(iev,ista,ifreq), dresid2, &
                staamp(iev,ista,ifreq), predamp, dph(iev,ista), prefase

           sumsq = sumsq+(d(ista+naddat)**2+d(ista+naddat+nsta(iev))**2) &
                      *stddevdata(iev)**2
           sumsqtemp = sumsqtemp+(d(ista+naddat)**2+ &
                        d(ista+naddat+nsta(iev))**2) &
                      *stddevdata(iev)**2

           sumsqph = sumsqph + phasemis**2
	   sumsqamp = sumsqamp + ((predamp-staamp(iev,ista,ifreq)) &
                                  /amprms(iev,ifreq))**2

          enddo
          totsumsqph = totsumsqph + sumsqph/freq(ifreq)**2
	  totsumsqamp = totsumsqamp + sumsqamp
          rmsphase(iev) = sqrt(sumsqph/nsta(iev))/freq(ifreq)
	  rmsamp(iev)   = sqrt(sumsqamp/nsta(iev))
	  rmsdata(iev) = sqrt(sumsqtemp/(2*nsta(iev)))
          naddat = naddat + 2*nsta(iev)
        enddo
!  rmsph multiplied by 2 because only half the number of observations
        rmsph = sqrt(2.*totsumsqph/nobs)
	rmsamplitude = sqrt(2.*totsumsqamp/nobs)
        sigma2 = sqrt(sumsq/nobs)
!        do j = 1, npsub
!          stddev(j+i6) = sqrt(subgtginv(j,j))
        do j = 1, np
          stddev(j) = sqrt(gtginv(j,j))
        enddo

!  Calculate resolution matrix for isotropic velocity parameters
        do i1 = i6+1,i6+nnodes
	  do i2 = i6+1,i6+nnodes
	    resparm =0.0
	    do j = 1, np
              resparm = resparm + gtginv(i1,j)*savegtg(j,i2)
            enddo
	    isores(i1-i6,i2-i6) = resparm
	  enddo
	enddo  

        write(10,*) 'nobs',nobs, 'rank', rank, '  &
             rank from isovel params ', rank2
        write(10,*) 'rank from anisovel params ', rank4
        write(10,*) ' rank from station corrections, etc.', rank3
        write(10,*) icnttot, icnttot2,icnt, ' iterations', sigma2,  &
            ' unnormalized rms misfit', rmsph, 'rms phase misfit,  s', &
               '  rms amp mistfit   ', rmsamplitude

        write(11,*) 'nobs',nobs, 'rank', rank,  &
             ' rank from isovel params ', rank2
        write(11,*) 'rank from anisovel params ', rank4
        write(11,*) ' rank from station corrections, etc.', rank3
        write(11,*) icnttot, icnttot2,icnt, ' iterations', sigma2,  &
            ' unnormalized rms misfit', rmsph, 'rms phase misfit,  s', &
               '  rms amp mistfit   ', rmsamplitude

! find median event rms phase misfit
        imed1 = (nevents+1)/2
        imed2 = (nevents+2)/2
        do iev = 1, nevents
           sortrms(iev) = rmsphase(iev)
        enddo
        call shell(nevents, sortrms)
        rmsmedian = (sortrms(imed1)+sortrms(imed2))/2.0
        write(10,*) 'median event misfit', rmsmedian
        write(11,*) 'median event misfit', rmsmedian
        do iev = 1, nevents
          ip = (iev-1)*6
          write(10,*) 'event ', idnum(iev), rmsdata(iev),  &
         ' data std dev', rmsphase(iev), 'rms phase misfit  s', &
         '  rms amp mistfit   ',  rmsamp(iev)
          write(11,*) 'event ', idnum(iev), rmsdata(iev),  &
         ' data std dev', rmsphase(iev), 'rms phase misfit  s', &
         '  rms amp mistfit   ',  rmsamp(iev)

          wvaz1= stazim1(iev)/convdeg
          wvaz2 = stazim2(iev)/convdeg
          stdwvaz1 = stddev(ip+3)/convdeg
          stdwvaz2 = stddev(ip+4)/convdeg
          write(10,*) wvaz1, startamp1(iev), stphase1(iev)
          write(10,*) wvaz2,startamp2(iev), stphase2(iev)
          write(11,*) wvaz1, startamp1(iev), stphase1(iev)
          write(11,*) wvaz2,startamp2(iev), stphase2(iev)
        enddo

!   write parameter covariance matrix for isotropic velocity parameters
!   for use in calculating variance of velocity model.  Could be smart and
!   write only upper or lower triangle of symmetric matrix to save space.

        write(16,*) nnodes
	do ii = 1, nnodes
          do jj = 1, nnodes
!	    write(16,*) subgtginv(ii,jj)
	    write(16,*) gtginv(ii+i6,jj+i6)
	    covar(ii,jj) = gtginv(ii+i6,jj+i6)
          enddo
        enddo

!  write resolution matrix for isotropic velocity parameters
        write(17,*) nnodes
	do ii = 1, nnodes
	  write(17,*) ii
	  write(17,*) (isores(ii,jj), jj=1,nnodes)
	enddo
	  
        do ii = 1, nnodes
          write(10,*) ii, crrntmod(ii+i6), stddev(ii+i6),residdiag(ii), &
                 nodelon(ii),nodelat(ii)
          write(11,*) ii, crrntmod(ii+i6), stddev(ii+i6),residdiag(ii), &
                 nodelon(ii),nodelat(ii)
	  vchange(ii) = crrntmod(ii+i6) - nodevel(ii)
!          write(30,*) ii,nodelon(ii),nodelat(ii), 
!     1              crrntmod(ii+i6),stddev(ii+i6),vchange(ii)
        enddo
        write(30,*) nnodes
        do ii = 1, iarea
           iii= i6 + nnodes +ii
           iiii = iii+iarea
           write(10,*) crrntmod(iii),stddev(iii), crrntmod(iiii), &
                          stddev(iiii), nodelon(ii),nodelat(ii)
           write(11,*) crrntmod(iii),stddev(iii), crrntmod(iiii), &
                          stddev(iiii), nodelon(ii),nodelat(ii)
           write(30,1276) nodelon(ii),nodelat(ii),  &
               crrntmod(iii),stddev(iii), crrntmod(iiii),stddev(iiii)
	enddo
1276    format(f9.3,f8.3,4(f8.4))
        jcnt = 0
        do ii = 1,mxnsta
          if (nevntsta(ii).gt.0) then
            jcnt = jcnt +1
            ip = jcnt + npnoamp
            write(10,*) ii, ampmult(jcnt), stddev(ip),' ', staname(ii)
            write(11,*) ii, ampmult(jcnt), stddev(ip),' ', staname(ii)
          endif
        enddo
	if (ntyp1.ne.0) then
 	write (10,*) 'station phase corrections'
 	write (11,*) 'station phase corrections'
         do ii=1, ntyp1                                                      
           ip=npp+ii                                                         
           write(10,*) ii, ityp1sta(ii),phcor(ii), stddev(ip)                
           write(11,*) ii, ityp1sta(ii),phcor(ii), stddev(ip)                
        enddo 
	endif                                                               
	write(10,*) 'attenuation correction'
	write(11,*) 'attenuation correction'
        write(10,*) gamma, stddev(np)  
        write(11,*) gamma, stddev(np)

	do ii = 1, nnodes
!	  write(20,*) nodelon(ii),nodelat(ii), residrow(ii)
	  write(21,*) nodelon(ii),nodelat(ii),residdiag(ii)
	enddo
!  change predicted phase velocities from a priori model to new results on
!  same grid as a priori model
        
	call updateapri(fendvel,unifvel,preunifvel,dampvel,nnodes)
	
      enddo
!  end loop over frequencies, which is silly since nfreq should always = 1

900   close(unit = 10)
      close(unit = 11)
      close(unit = 12)
      close(unit = 13)
      close(unit = 14)
      close(unit = 15)
      close(unit = 16)
      close(unit = 17) 
      close(unit = 18)
      close(unit = 21)
      close(unit = 30) 
      close(unit = 66)  
      close(unit = 99)
      print *, "bp -1"   

      stop 
      end


      SUBROUTINE dlubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
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


      SUBROUTINE dludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX=4000,TINY=1.0e-20)
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

      subroutine disthead(slat,slon,flat,flon,delta,azim)
!  Calculates distance and azimuth on sphere from starting point s 
!  to finishing point f
      dtor= 3.1415928/180.
      slt = slat*dtor
      sln = slon*dtor
      flt = flat*dtor
      fln = flon*dtor
      delta = acos(sin(slt)*sin(flt)+cos(slt)*cos(flt)*cos(fln-sln))
      azim = atan2(sin(fln-sln)*cos(flt), &
         sin(flt)*cos(slt) - cos(fln-sln)*cos(flt)*sin(slt))
      delta = delta/dtor
      azim = azim/dtor
      return
      end

      subroutine gohead(slat,slon,delta,azim,flat,flon)
! Calculates final latitude and longitude f when starting at point s
! traveling delta degrees on sphere at initial heading azim
      dtor= 3.1415928/180.
      slt = slat*dtor
      dlta = delta*dtor
      azm = azim*dtor
      flat = asin(cos(slt)*sin(dlta)*cos(azm) + sin(slt)*cos(dlta))
      flon = atan2(sin(dlta)*sin(azm),  &
          cos(slt)*cos(dlta) - sin(slt)*sin(dlta)*cos(azm))
      flat = flat/dtor
      flon = slon + flon/dtor
      if (flon.gt.360.) flon = flon - 360.
      if (flon.lt.-360.) flon = flon + 360.
      return
      end

      subroutine shell(n,a)
!  shell sort from Press et al. Sorts a, replacing original, in ascending order
      integer n, i,j, inc
      real a(n), v
      inc = 1
1     inc = 3*inc + 1
      if (inc.le.n) go to 1
2       continue
        inc = inc/3
        do i = inc+1, n
          v = a(i)
          j = i
3         if (a(j-inc).gt.v) then
            a(j) = a(j-inc)
            j = j-inc
            if (j.le.inc) go to 4
          go to 3
          endif
4         a(j) = v
        enddo
        if(inc.gt.1) go to 2
        return
        end

      subroutine search(pb)
!  uses grid search to find estimates of best fitting wave parameters for each event,
!  with fixed velocity structure
      parameter (maxnfreq=1, maxnsta=300, maxpts=20000, nparam = 4000, &
                  maxnobs = 25000, maxnodes = 2000, maxevnts = 400, &
                  maxnxints = 401, maxnyints = 401,ndeg = 81 )

      real*4 streal(maxevnts,maxnsta,maxnfreq)
      real*4 stimag(maxevnts,maxnsta,maxnfreq)
      real*4 rloc(maxevnts,maxnsta), azloc(maxevnts,maxnsta)
      real*4 d(maxnobs), xsta(maxevnts,maxnsta),ysta(maxevnts,maxnsta)
      real*4 freq(maxnfreq)
      real*4 stddevdata(maxevnts)
      real*4 dtime(maxnsta), avslow(maxnsta)
      real*4 xbegkern,dxkern
      real*4 xbox(maxevnts,4),ybox(maxevnts,4)
      real*4 wgtnode1(maxnsta,maxnodes,ndeg)
      real*4 ampwgtnode1(maxnsta,maxnodes,ndeg)
      real*4 dxnode,dynode
      real*8 sensitivity(maxnxints,maxnyints)
      real*8 ampsens(maxnxints,maxnyints)
      real*4 damp1per(maxnsta), damp2per(maxnsta)
      real*4 dphase1(maxnsta),dphase2(maxnsta)
      real*4 phase1(maxnsta), phase2(maxnsta)
      real*4 ybegkern,dykern
      real*4 xnode(maxevnts,maxnodes),ynode(maxevnts,maxnodes)
      real*4 dtime1(maxnsta),dtime2(maxnsta),phcor(50)   
      real*4 phase(maxnsta,ndeg),dphase(maxnsta,ndeg), &
              dampper(maxnsta,ndeg)
      real*4 appvel(maxnodes),ampmult(maxnsta)
      real*4 xmin(maxevnts,ndeg), unifvel
      real*4 gamp(maxnsta,2)
      real*4 pb(6)
!	  real*4 bmisfit
        
      integer*4 nsta(maxevnts), iref(maxevnts)
      integer*4 istanum(maxevnts,maxnsta),istacor(maxnsta)
      integer*4   istavar(maxnsta)
      integer*4 nxkern,nykern, nnodes
      integer*4 ityp1sta(maxnsta),ntyp1 
                                                                        
      common /residua/ sensitivity,ampsens,d,rloc,azloc,freq, &
        xsta,dtime,avslow,streal,stimag,stddevdata,phcor,  &
        xbox,ybox,ysta,dxkern,dykern,dxnode,dynode, &
        unifvel,appvel,ampmult,gamma, &
        xnode, ynode,wgtnode1,ampwgtnode1,xmin, &
        phase,dphase,dampper,istavar,istanum,nxkern,nykern, &
        ityp1sta,ntyp1,iref,nsta,iev,naddat,ifreq,nnodes

      twopi = 3.1415928*2.
      onepi = 3.1415928
      convdeg = 3.1415928/180.
      azimlmn1 = -20.0*convdeg
      azimlmx1 = 20.01*convdeg
      azimint1 = 2.0*convdeg
      azimlmn2 = -40.0*convdeg
      azimlmx2 = 40.01*convdeg
      azimint2 = 4.0*convdeg
      phaselmn1 = -0.25
      phaselmx1 = 0.25
!      phaseint1 = 0.01
      phaseint1 = 0.02
      phaselmn2 = -0.25
      phaselmx2 = 0.25
!      phaseint2 = 0.01
      phaseint2 = 0.02
      iter = 0
      bmisfit = 1.0E+38
      funiv = freq(1)/unifvel
100   iter = iter + 1
      do srhazim1 = azimlmn1, azimlmx1,azimint1
        cs1 = cos(srhazim1)
        sn1 = sin(srhazim1)
        ideg1 = int(srhazim1/convdeg+(ndeg-1)/2 +0.01) + 1
        do srhazim2 = azimlmn2,azimlmx2,azimint2
          if (srhazim2.ne.srhazim1) then
          cs2 = cos(srhazim2)
          sn2 = sin(srhazim2)
          ideg2 = int(srhazim2/convdeg+(ndeg-1)/2 + 0.01) + 1
!	  write(*,*) ideg1,srhazim1/convdeg,ideg2,srhazim2/convdeg
          do srhphase1 = phaselmn1,phaselmx1,phaseint1
            do srhphase2 = phaselmn2,phaselmx2,phaseint2
              do 10 ista = 1, nsta(iev)
                xstatemp = xsta(iev,ista)*cs1 + ysta(iev,ista)*sn1
!		ystatemp = xsta(iev,ista)*sn1 + ysta(iev,ista)*cs1
  		phase1(ista) = (xstatemp-xmin(iev,ideg1))*funiv
		xstatemp = xsta(iev,ista)*cs2 + ysta(iev,ista)*sn2
!		ystatemp = xsta(iev,ista)*sn2 + ysta(iev,ista)*cs2
  		phase2(ista) = (xstatemp-xmin(iev,ideg2))*funiv
!  calculate effects of uniform velocity
  10 	      enddo
  	      gtg11 = 0.0
	      gtg12 = 0.0
	      gtg22 = 0.0
	      gtd1 = 0.0
	      gtd2 = 0.0
	      do ista = 1,nsta(iev)
                prefase1 = ( (phase1(ista)+dphase(ista,ideg1)) &
          	- (phase1(iref(iev))+dphase(iref(iev),ideg1))) &
                    + srhphase1
        	prefase2 = ( (phase2(ista)+dphase(ista,ideg2)) &
          	- (phase2(iref(iev))+dphase(iref(iev),ideg2))) &
                    + srhphase2
		if (ntyp1.ne.0) then
        	do itypst=1, ntyp1                                   
         	  if (istanum(iev,ista).eq.ityp1sta(itypst)) then     
          	  prefase1=prefase1+phcor(itypst)                    
          	  prefase2=prefase2+phcor(itypst)                    
	 	  endif                                               
        	enddo 
		endif
		
		ampfac  = ampmult(istavar(istanum(iev,ista))) &
       			*exp(-gamma*xsta(iev,ista))                                              
		ampfac1 = (1.+dampper(ista,ideg1))*ampfac
		ampfac2 = (1.+dampper(ista,ideg2))*ampfac
!  now do linear inversion for best amplitudes given these trial azimuths and phases
!  This is a completely linear inversion at this point, so can use zero starting model
!  for amplitudes, and for 2x2 inversion of symmetric GTG matrix can use analytical inverse.
		
!  data vector and partial derivatives listed event by event with all
!  real data for first event followed by imaginary data

               	cosph1 = cos(prefase1*twopi)
          	cosph2 = cos(prefase2*twopi)
          	sinph1 = sin(prefase1*twopi)
          	sinph2 = sin(prefase2*twopi)
          	gamp(ista,1) = ampfac1*cosph1/stddevdata(iev)
		gamp(ista,2) = ampfac2*cosph2/stddevdata(iev)
		gamp(ista+nsta(iev),1)=-ampfac1*sinph1/stddevdata(iev)
		gamp(ista+nsta(iev),2)=-ampfac2*sinph2/stddevdata(iev)
		gtg11 = gtg11+gamp(ista,1)**2+gamp(ista+nsta(iev),1)**2
		gtg22 = gtg11+gamp(ista,2)**2+gamp(ista+nsta(iev),2)**2
		gtg12 = gtg12+gamp(ista,1)*gamp(ista,2) &
       		        + gamp(ista+nsta(iev),1)*gamp(ista+nsta(iev),2)
		gtd1 = gtd1 + (gamp(ista,1)*streal(iev,ista,ifreq) &
       		  +gamp(ista+nsta(iev),1)*stimag(iev,ista,ifreq)) &
                    /stddevdata(iev)
		gtd2 = gtd2 + (gamp(ista,2)*streal(iev,ista,ifreq) &
       		  +gamp(ista+nsta(iev),2)*stimag(iev,ista,ifreq)) &
                    /stddevdata(iev)
	      enddo
	      c2ab = gtg12*gtg12-gtg11*gtg22
	      gtginv11 = -gtg22/c2ab
	      gtginv22 = -gtg11/c2ab
	      gtginv12 = gtg12/c2ab
!   stampmod1 and stampmod2 should be best estimates of amps for these angles and phases
	      stampmod1 = gtginv11*gtd1 + gtginv12*gtd2
	      stampmod2 = gtginv12*gtd1 + gtginv22*gtd2
!  now calculate residuals to starting model 
	      pmisfit = 0.0
	      do ista = 1,nsta(iev) 
	        pmisfit = pmisfit+(gamp(ista,1)*stampmod1 &
       		  +gamp(ista,2)*stampmod2  &
       		  - streal(iev,ista,ifreq)/stddevdata(iev))**2
	        pmisfit = pmisfit + (gamp(ista+nsta(iev),1)*stampmod1 &
       		  + gamp(ista+nsta(iev),2)*stampmod2  &
                   - stimag(iev,ista,ifreq)/stddevdata(iev))**2
     	      enddo
	      if (pmisfit.lt.bmisfit) then
	        bmisfit = pmisfit
	        pb(1) = stampmod1
	        pb(2) = stampmod2
	        pb(3) = srhazim1
	        pb(4) = srhazim2
	        pb(5) = srhphase1
	        pb(6) = srhphase2
	      endif
	    enddo
	  enddo
	  endif
	enddo
      enddo
! now refine search in vicinity of best model from coarse search
      if (iter.eq.1) then
        azimlmn1 = pb(3)-2.0*convdeg
        azimlmx1 = pb(3)+2.0*convdeg
        azimint1 = 1.0*convdeg
        azimlmn2 = pb(4)-4.0*convdeg
	if (azimlmn2.lt.-40.*convdeg) azimlmn2 = -40.*convdeg
        azimlmx2 = pb(4)+4.0*convdeg
	if (azimlmx2.gt.40.*convdeg) azimlmx2 = 40.*convdeg
        azimint2 = 2.0*convdeg
!        phaselmn1 = pb(5)-0.01
!        phaselmx1 = pb(5)+0.01
        phaselmn1 = pb(5)-0.02
        phaselmx1 = pb(5)+0.02
        phaseint1 = 0.0025
!        phaselmn2 = pb(6)-0.01
!        phaselmx2 = pb(6)+0.01
        phaselmn2 = pb(6)-0.02
        phaselmx2 = pb(6)+0.02
        phaseint2 = 0.0025
	go to 100
      endif
!      write(*,*) iev,pmisfit,pb(1),pb(2),pb(3),pb(4),pb(5),pb(6)
      return
      end

      subroutine actlmsft(actmisrms)
       parameter (maxnfreq=1, maxnsta=300, maxpts=20000, nparam = 4000, &
                  maxnobs = 25000, maxnodes = 2000, maxevnts = 400, &
                  maxnxints = 401, maxnyints = 401,ndeg = 81 )

      real*4 streal(maxevnts,maxnsta,maxnfreq)
      real*4 stimag(maxevnts,maxnsta,maxnfreq)
      real*4 rloc(maxevnts,maxnsta), azloc(maxevnts,maxnsta)
      real*4 d(maxnobs), xsta(maxevnts,maxnsta),ysta(maxevnts,maxnsta)
      real*4 freq(maxnfreq)
      real*4 stddevdata(maxevnts)
      real*4 dtime(maxnsta), avslow(maxnsta)
      real*4 xbegkern,dxkern
      real*4 xbox(maxevnts,4),ybox(maxevnts,4)
      real*4 wgtnode1(maxnsta,maxnodes,ndeg)
      real*4 ampwgtnode1(maxnsta,maxnodes,ndeg)
      real*4 dxnode,dynode
      real*8 sensitivity(maxnxints,maxnyints)
      real*8 ampsens(maxnxints,maxnyints)
      real*4 damp1per(maxnsta), damp2per(maxnsta)
      real*4 dphase1(maxnsta),dphase2(maxnsta)
      real*4 phase1(maxnsta), phase2(maxnsta)
      real*4 ybegkern,dykern
      real*4 xnode(maxevnts,maxnodes),ynode(maxevnts,maxnodes)
      real*4 dtime1(maxnsta),dtime2(maxnsta),phcor(50)   
      real*4 phase(maxnsta,ndeg),dphase(maxnsta,ndeg), &
              dampper(maxnsta,ndeg)
      real*4 appvel(maxnodes),ampmult(maxnsta)
      real*4 xmin(maxevnts,ndeg), unifvel
      real*4 bazi(maxevnts,maxnsta)
      real*4 cos2node(maxnodes,ndeg),sin2node(maxnodes,ndeg)
      real*4 crrntmod(nparam)      
      real*4 gamp(maxnsta,2)
      real*4 pb(6)
      real*4 startamp1(maxevnts),startamp2(maxevnts)
      real*4 stphase1(maxevnts),stphase2(maxevnts)
      real*4 stazim1(maxevnts),stazim2(maxevnts)

      integer*4 nsta(maxevnts), iref(maxevnts)
      integer*4 istanum(maxevnts,maxnsta),istacor(maxnsta)
      integer*4 idnode(maxnodes)
      integer*4   istavar(maxnsta)
      integer*4 nxkern,nykern, nnodes
      integer*4 ityp1sta(maxnsta),ntyp1 
      integer*4 nevents
                                                                                
      common /residua/ sensitivity,ampsens,d,rloc,azloc,freq, &
        xsta,dtime,avslow,streal,stimag,stddevdata,phcor,  &
        xbox,ybox,ysta,dxkern,dykern,dxnode,dynode, &
        unifvel,appvel,ampmult,gamma, &
        xnode, ynode,wgtnode1,ampwgtnode1,xmin, &
        phase,dphase,dampper,istavar,istanum,nxkern,nykern, &
        ityp1sta,ntyp1,iref,nsta,iev,naddat,ifreq,nnodes                 

      common /msft/ bazi,cos2node,sin2node,crrntmod, &
          startamp1,startamp2,stphase1,stphase2,stazim1,stazim2, &
          nevents,idnode,ideg1,ideg2,nobs,i6,iarea

        twopi = 3.1415928*2.
        onepi = 3.1415928
        convdeg = 3.1415928/180.
	actmisfit = 0.0
       do iev = 1, nevents
!   calculate the sensitivity kernel for fixed angles

        do 125 idegg = 1,2
	  if (idegg.eq.1) ideg = ideg1
	  if (idegg.eq.2) ideg = ideg2
          azimt = ((ideg-1.) - (ndeg-1.)/2.)*convdeg
          cs2n = cos(2.0*(convdeg*bazi(iev,iref(iev))- azimt))
	  sn2n = sin(2.0*(convdeg*bazi(iev,iref(iev))- azimt))
	  cs1z = cos(azimt)
	  sn1z = sin(azimt)

          do ii = 1,nnodes
! don't really need arrays here in this version, but used throughout so retained
            cos2node(ii,ideg)=cs2n
            sin2node(ii,ideg)=sn2n
          enddo

!  calculate current apparent velocity at each node

        do ii = 1, nnodes
          iii  = ii+i6
          jjj  = i6+nnodes + idnode(ii)
	  jjjj = jjj + iarea
          appvel(ii) =       crrntmod(iii) &
          + cos2node(ii,ideg)*crrntmod(jjj) &
          + sin2node(ii,ideg)*crrntmod(jjjj)
!       write(14,*) iev,ideg,ii, appvel(ii)
        enddo

!  loop over stations, generating sensitivity kernels for each
        do 102 ista = 1, nsta(iev)

          do ii = 1, nnodes
	      wgtnode1(ista,ii,ideg) = 0.0
	    ampwgtnode1(ista,ii,ideg) = 0.0
          enddo

          xstatemp =  &
             xsta(iev,ista)*cs1z &
             + ysta(iev,ista)*sn1z               
     
          ystatemp =  &
             -xsta(iev,ista)*sn1z &
             + ysta(iev,ista)*cs1z     

          do ii = 1,nnodes

           xnodetemp =  &
             xnode(iev,ii)*cs1z &
             + ynode(iev,ii)*sn1z               
     
           ynodetemp =  &
             -xnode(iev,ii)*sn1z &
             + ynode(iev,ii)*cs1z     


           xstanode = xnodetemp - xstatemp
	   ystanode = ynodetemp - ystatemp

!  find nearest point in sensitivity kernel - kernels should have been
!  calculated with smoothing so that they represent sensitivity to nodal
!  coefficient rather than velocity at point  - could interpolate, but 
! this is sufficiently accurate if kernels on fine enough scale (i.e.,
!  dxkern << smoothing length
         if (xnodetemp.ge.xmin(iev,ideg)) then
          if (xstanode.ge.0.0) then
            ixindex = int( xstanode/dxkern +0.5) + (nxkern+1)/2
	  else
	    ixindex = int( xstanode/dxkern -0.5) + (nxkern+1)/2
	  endif
          if (ystanode.ge.0.0) then
            iyindex = int( ystanode/dykern +0.5) + (nykern+1)/2
	  else
	    iyindex = int( ystanode/dykern -0.5) + (nykern+1)/2
	  endif

          if(  ixindex .lt.1 .or. ixindex .gt. nxkern &
           .or. iyindex .lt.1 .or. iyindex .gt. nxkern) then 


	    write(*,*) 'ixindex,iyindex,iev,ista,ideg',  &
                   ixindex,iyindex,iev,ista,ideg
	     wgtnode1(ista,ii,ideg)    = 0.
	     ampwgtnode1(ista,ii,ideg) = 0.  

 	   else

             wgtnode1(ista,ii,ideg) = sensitivity(ixindex,iyindex)
             ampwgtnode1(ista,ii,ideg) =     ampsens(ixindex,iyindex)
           endif 
	  else 
	     wgtnode1(ista,ii,ideg)    = 0.
	     ampwgtnode1(ista,ii,ideg) = 0.  
           endif	   
          enddo

	 dphase(ista,ideg) = 0.
         dampper(ista,ideg) = 0.
!  corrections for second order effects of large velocity changes in version 587
!   old version commented out
          do inode = 1, nnodes
!          dphase(ista,ideg) = dphase(ista,ideg) 
!     1  + (1.0/twopi)*wgtnode1(ista,inode,ideg)
!     1           *(appvel(inode)-unifvel)/unifvel
          dphase(ista,ideg) = dphase(ista,ideg)  &
         + (1.0/twopi)*wgtnode1(ista,inode,ideg) &
                  *(appvel(inode)-unifvel)/appvel(inode) &
                  /(appvel(inode)/unifvel)

          enddo


	  do inode =1, nnodes
          dampper(ista,ideg) = dampper(ista,ideg)  &
         + ampwgtnode1(ista,inode,ideg) * &
!     1   (appvel(inode)-unifvel)/unifvel &
           (appvel(inode)-unifvel)/appvel(inode) &
                  /(appvel(inode)/unifvel)

          enddo

!  end loop over stations
 102 	  enddo
!  end loop over angles
 125     enddo

!  calculate effect of uniform velocity for first plane wave
       do 10 ista = 1, nsta(iev)
       xstatemp =  &
          xsta(iev,ista)*cos(stazim1(iev)) &
        + ysta(iev,ista)*sin(stazim1(iev))               
     
   	dtime1(ista) = (xstatemp-xmin(iev,ideg1))/unifvel

  10 	  enddo

!       calculate for second plane wave
  
        do 20 ista = 1, nsta(iev)
         
       xstatemp =  &
          xsta(iev,ista)*cos(stazim2(iev)) &
        + ysta(iev,ista)*sin(stazim2(iev))               
     
   	dtime2(ista) = (xstatemp-xmin(iev,ideg2))/unifvel

	   phase1(ista) = dtime1(ista)*freq(1)
	   phase2(ista) = dtime2(ista)*freq(1)
	  dphase1(ista) =   dphase(ista,ideg1)
	  dphase2(ista) =   dphase(ista,ideg2)  
          damp1per(ista) = dampper(ista,ideg1)
	  damp2per(ista) = dampper(ista,ideg2) 

  20    enddo
   
!  convert time to phase and assign effects of lateral heterogeneities

        do ista = 1, nsta(iev)


        prefase1 = ( (phase1(ista)+dphase1(ista)) &
         - (phase1(iref(iev))+dphase1(iref(iev)))) &
                    + stphase1(iev)
        prefase2 = ( (phase2(ista)+dphase2(ista)) &
         - (phase2(iref(iev))+dphase2(iref(iev)))) &
                    + stphase2(iev)

	  staamp1 = startamp1(iev)*(1.+damp1per(ista))
	  staamp2 = startamp2(iev)*(1.+damp2per(ista))

!        parph1cor1=0.                                        
!        parph2cor1=0.                                        
	 
	 if (ntyp1.ne.0) then
!         itypst2=0 
        do itypst=1, ntyp1                                    
         if (istanum(iev,ista).eq.ityp1sta(itypst)) then     
          prefase1=prefase1+phcor(itypst)                    
          prefase2=prefase2+phcor(itypst)                    
!          parph1cor1=1.                                      
!          parph2cor1=1.                                      
!	  itypst2=itypst                                     
	 endif                                               
        enddo 
	endif                                               


          cosph1 = cos(prefase1*twopi)
          cosph2 = cos(prefase2*twopi)
          sinph1 = sin(prefase1*twopi)
          sinph2 = sin(prefase2*twopi)
          
          prereal = staamp1*cosph1+staamp2*cosph2
!	  write(*,*) prereal  
	  preimag = -1.0*(staamp1*sinph1+ staamp2*sinph2)
!           write(*,*) ampmult(istavar(istanum(iev,ista)))
          prereal = prereal*ampmult(istavar(istanum(iev,ista))) &
                 *exp(-gamma*xsta(iev,ista))
!          write(*,*) prereal
          preimag = preimag*ampmult(istavar(istanum(iev,ista))) &
                     *exp(-gamma*xsta(iev,ista))




!  data vector and partial derivatives listed event by event with all
!  real data for first event followed by imaginary data, then onto next event
!  d contains misfit to starting model


!          kreal = ista + naddat
!          kimag = kreal + nsta(iev)
          dkreal = (streal(iev,ista,ifreq) - prereal)/stddevdata(iev)
          dkimag = (stimag(iev,ista,ifreq) - preimag)/stddevdata(iev)
	  write(14,*) iev,ista,streal(iev,ista,ifreq),stimag(iev,ista,ifreq), &
               prereal,preimag
          actmisfit = actmisfit + dkreal*dkreal+dkimag*dkimag
         enddo
!	 write(*,*) iev, ' stddevdata ',stddevdata(iev),' actmisfit ', actmisfit 
        enddo
	actmisrms = sqrt(actmisfit/nobs)
!	write(*,*) 'nobs ', nobs, ' ifreq ', ifreq,' unifvel ',unifvel
	return
	end
	
	
	
      subroutine genreggrid(freggrid,ntot,nxpt,dx,dy)
!***************************************************************************
!*   Create regularly spaced grid points in an area with center point (slat,slon). *
!*   In order to avoid distortion at high latitude, change coordinates so   *
!*   center of area lies along the equator of coordinate system and there is 
!*   compensation for curvature.                                           *
!*************************************************************************** 
!  WARNING   checkout scheme by running separately first gridreg < brwestreginp.dat
       parameter (nxptmx=100, nyptmx=100)
       parameter (maxnodes = 2000)

       real*4 dazim(nxptmx),ddelt(nyptmx),glat(nxptmx,nyptmx)
       real*4 glon(nxptmx,nyptmx)
      real*4 boxlat(4), boxlon(4)
      real*4 nodelat(maxnodes),nodelon(maxnodes)
       character*70 freggrid
       common /gengrid/ nodelat,nodelon,boxlat,boxlon
       open(15, file = freggrid, status = 'unknown')
       radius = 6371.
       pi = 3.14158
       convdeg = 2.0*pi/360.
!   'Enter the center point (slat, slon):'
       read(15,*) slat,slon
       sdelta=90.
!      sazimz is azimuth from the center point that will tilt grid relative to North
       read(15,*) sazimz
!  delx is increment in degrees at equator - the true increment will vary with latitude from
!  equator of projection
       read(15,*) nxpt, delx,begx
! to follow our traditional convention of listing points from right to left and bottom to top.
!  the degree increment for distance from the pole, dely, should be negative and the
!  beginning point should be farthest in degrees from projection pole (assuming sazimz is
!  northwards), (e.g., begy = 6 (degrees south of center point) is 96 degrees from pole)
!  begx would be negative and increase with positive delx
       read(15,*) nypt,dely, begy
! read in corners for wave intercepts to pass along to main program
       do i = 1, 4
         read(15, *) boxlat(i), boxlon(i)	
       enddo
       degdist = 2.0*pi*radius/360.
       dx = delx*degdist
       dy = dely*degdist
       call gohead(slat,slon,sdelta,sazimz,plat,plon)
       call disthead(plat,plon,slat,slon,delta,azimz)
       ntot = nxpt*nypt
       do 150 j=1,nypt
         delyinc = (j-1)*dely +begy	 
         delt = delta + delyinc	 
         do 100 i=1,nxpt
	   azim = azimz + (begx + (i-1)*delx)/cos(convdeg*delyinc)
	   call gohead(plat,plon,delt,azim,glat(i,j),glon(i,j))
	   if (glon(i,j).gt.180.0) then
	     glon(i,j) = glon(i,j) - 360.0
	   endif
  100	 enddo
  150  enddo 
!  arrange points in standard order in single dimensioned arrays
       nodecnt = 0
       do i = 1, nxpt
         do j = 1, nypt
	   nodecnt = nodecnt + 1
           nodelat(nodecnt) = glat(i,j)
	   nodelon(nodecnt) = glon(i,j)
	 enddo
       enddo
       close(15)  
  200  format(f7.3,f10.3)
       return
       end  
       
       
      subroutine assignstrt(startvel,nodevel,preunifvel,ntot, &
                                  nxpt,dxnode,dynode)
      parameter (maxnodes = 2000, maxstart = 200)
      character*70 startvel
      real*4 boxlat(4), boxlon(4)
      real*4 nodelat(maxnodes),nodelon(maxnodes),nodevel(maxnodes)
      real*4 latmin,latmax,lonmin,lonmax
      real*4 prdlon(maxstart),prdlat(maxstart),prdv(maxstart,maxstart)
      common /gengrid/ nodelat,nodelon,boxlat,boxlon
      common /updatest/ prdlon,prdlat,prdv,dylat,dxlon, &
           beglat,endlat,dlat,beglon,endlon,dlon,nlat,nlon
     
      pi = 3.1415928
      convdeg = 3.1415928/180.

      deg2km = 111.194
      open(60, file = startvel)
!  assigns start values for node parameters based on a priori predictions of
!  starting model, previously generated.  Important that a priori grid extend
!  at least as far as model grid and that it is nearly as closely or more closely
!  spaced than model grid. Uses linear weighting between nodes to generate average value for
!  starting values.  A priori grid does not have same spacing and is assumed to 
!  be regularly spaced on Mercator projection, increasing in both lon and lat
!  
      read(60,*) beglat,endlat,dlat,beglon,endlon,dlon
      read(60,*) preunifvel
      nlat = (endlat-beglat)/dlat +1.01
      nlon = (endlon-beglon)/dlon + 1.01
!  file of phase velocity predictions are assumed to run through longitudes first,
!  before incrementing in latitude
      do ilat = 1, nlat
        do ilon = 1, nlon
	  read(60,*) prdlon(ilon),prdlat(ilat),prdv(ilon,ilat)
	enddo
      enddo
!  convert spacing of model back to degrees
      dylat = dynode/deg2km
      dxlon = dxnode/deg2km
      do i = 1, ntot
        sclfac = cos(convdeg*nodelat(i))
	dxxlon = dxlon/sclfac
        latmin = nodelat(i) - dylat
	latmax = nodelat(i) + dylat
	lonmin = nodelon(i) - dxxlon
	lonmax = nodelon(i) + dxxlon
	minlat = (latmin - beglat)/dlat +1
	if (minlat.lt.1) minlat = 1
	maxlat = (latmax-beglat)/dlat +1
	if (maxlat.gt.nlat) maxlat = nlat
	minlon = (lonmin-beglon)/dlon +1
	if (minlon.lt.1) minlon = 1
	maxlon = (lonmax-beglon)/dlon + 1
	if (maxlon.gt.nlon) maxlon = nlon
	wgtsum = 0.0
	sumvel = 0.0
	if ((minlat.gt.nlat).or.(maxlat.lt.1).or. &
                    (minlon.gt.nlon).or.(maxlon.lt.1)) then
           write (*,*)  &
          'PROBLEM. Node outside region of starting predictions '
	   write (*,*) i, nodelat(i),nodelon(i)
	endif
	do ii = minlon,maxlon
	  do jj = minlat,maxlat
	    wgt = (1.-abs(nodelat(i)-prdlat(jj))/dylat) &
                    *(1.-abs(nodelon(i)-prdlon(ii))/dxxlon)
            sumvel = sumvel + wgt*prdv(ii,jj)
	    wgtsum = wgtsum + wgt
	  enddo
	enddo
	nodevel(i) = sumvel/wgtsum
!	write(*,*) i, nodevel(i),nodelat(i),nodelon(i)
!	write(*,*) minlat,maxlat,minlon,maxlon
      enddo
      return
      end
 

      subroutine updateapri(fendvel,unifvel,preunifvel,dampvel,nnodes)
!  Take inversion solution of changes to starting model and add back through
!  linear interpolation to detailed, a priori starting model. Also, generate 
!  standard deviation of estimate
!
!  Slow process because nodes not on regular grid of lat&lon - have to search for
!  which nodes are within interpolation range for every point
      parameter (maxnodes = 2000, maxstart = 200)
      character*70 fendvel
      real*4 boxlat(4), boxlon(4)
      real*4 nodelat(maxnodes),nodelon(maxnodes),nodevel(maxnodes)
      real*4 wgt(maxnodes)
      real*4 latmin,latmax,lonmin,lonmax
      real*4 prdlon(maxstart),prdlat(maxstart),prdv(maxstart,maxstart)
      real*4 newvel(maxstart,maxstart),stddevap(maxstart,maxstart)
      real*4 covar(maxnodes,maxnodes),vchange(maxnodes)
      integer iwgt(maxnodes)
      common /gengrid/ nodelat,nodelon,boxlat,boxlon
      common /updatest/ prdlon,prdlat,prdv,dylat,dxlon, &
           beglat,endlat,dlat,beglon,endlon,dlon,nlat,nlon
      common /update2/ covar,vchange
     
      pi = 3.1415928
      convdeg = 3.1415928/180.

      open (70, file = fendvel)
      write(70,*) beglat,endlat,dlat,beglon,endlon,dlon
      sumall = 0.0
      do jj = 1,nlat
!  to correct for points on mercator grid being closer together away from equator
        sclfac = cos(prdlat(jj)*convdeg)
        do ii = 1,nlon
	  wgtsum = 0.0
	  velsum = 0.0
	  sumstdev = 0.0
	  nwgt = 0
	  do kk = 1,nnodes
	    xsep = abs(prdlon(ii)-nodelon(kk))*sclfac
	    if ((xsep.le.dxlon).and. &
                (abs(prdlat(jj)-nodelat(kk)).le.dylat))  then
              wgt(kk) = (1.- abs(prdlat(jj)-nodelat(kk))/dylat) &
                             *(1.- xsep/dxlon)
              nwgt = nwgt+1
              iwgt(nwgt) = kk
              wgtsum = wgtsum + wgt(kk)
	      velsum = velsum + wgt(kk)*vchange(kk)
	      sumstdev = sumstdev + wgt(kk)*sqrt(covar(kk,kk))
	    endif
	  enddo
	  if (nwgt.gt.0) then
      	    newvel(ii,jj)=prdv(ii,jj)+unifvel-preunifvel +velsum/wgtsum
	    stddevap(ii,jj) = sumstdev/wgtsum
!  Estimate uncertainty using covariance matrix - but for this application
!  we don't want to use full covariance including off-diagonal terms, because that
!  effectively gives estimates appropriate to averages over larger regions involving
!  multiple nodes.  Instead, what is desired is an estimate of the uncertainty as if
!  a node had been located at the interpolation point (although points in between really are
!  known better because the random variations from multiple nodes will tend to average out)..  	    
!	    do i = 1,nwgt
!	      do j = 1,nwgt
!	        sumcov = wgt(iwgt(i))*wgt(iwgt(j))*covar(iwgt(i),iwgt(j))
!	      enddo
!	    enddo
!	    stddevap(ii,jj) = sqrt(sumcov)/wgtsum
	   else
	     newvel(ii,jj) = prdv(ii,jj)+unifvel-preunifvel
	     stddevap(ii,jj) = dampvel
	   endif
!	  write(*,*) ii,jj,nwgt,(iwgt(ij),ij=1,nwgt),prdlon(ii),prdlat(jj)
	  sumall = sumall + newvel(ii,jj)
	enddo
      enddo
      avgvel = sumall/(nlat*nlon)
      write(70,*) avgvel
      do jj = 1,nlat
        do ii = 1,nlon
	  write(70,*) prdlon(ii),prdlat(jj),newvel(ii,jj),stddevap(ii,jj)
	enddo
      enddo
      
      close (70)
      return
      end

