c  anisopred takes given 3D model for anisotropy and predicts the cos2theta 
c  and sin2theta coefficients as function of period.  (Just uses SV
c  approximation and substitutes for anisoz16 if going from shorter periods to longer
c  where more periods are used than were done in inverting). 
c
c  anizopred < anisopredinp
c
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
 

	read(*,*) nper2, depth1, error1,ijk1,ijk2,ijk3,topd,bottomd
	
c  These are variables controlling the search for phase velocity (roots to equations) for
c  given structure.  nper2= # periods, depth1 default depth extent, error1, criterion for 
c  ending search, ijk1 = 4 is Rayleigh spherical earth, 3 is Love spherical earth, ijk2 = 0, ijk3 = 2
c  topd and bottomd control depth range in which partial derivatives are output  

        nobs = nper

        read(*,*) startmod
	read(*,*) flowmantle
	read(*,*) prioraniso
	
	   open(55, file = prioraniso)
	   read(55,*) beglat,endlat,dlat,beglon,endlon,dlon
	   
	
	open(10, file=startmod)
	
	read(10,*) beglat,endlat,dlat,beglon,endlon,dlon

	nlat = (endlat-beglat)/dlat +1.01
        nlon = (endlon-beglon)/dlon + 1.01
        nxy = nlat*nlon

	naniso = nxy
	
	ione = 1
	do i = 1,nper
	  read(*,*) anisoout(i)
	  open(60+i, file = anisoout(i))
	  write(60+i,*) beglat,endlat,dlat,beglon,endlon,dlon
	  write(60+i,*) ione
	enddo
			
        open(40, file = flowmantle)
	
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


c nlayall: 2X layers, including lower mantle parts
c np=nlay: number of layers in the lithosphere and upper mantle 
c nparam: number of layers which model needs to resolve
        np = nlay
        nparam = np  
        
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

       
c  initialize g matrix
        do irow = 1, nobs
        do icol = 1, np
          g(irow,icol) = 0.0                                                           !
        enddo                                                                          !
        enddo                                                                          !
       
c  data vector and partial derivatives listed event by event 
c***************************************************
c  repeat separately from here for cos2theta and sin2theta inversions
	   read(55,*) lontemp2,lattemp2
	   read(55,*) nlay2
	   do janso = 1, nlay2
	     read(55,*) thick(janso),stcs2th(janso),stsn2th(janso)
	   enddo
	   
         indaniso = 1
2000	 open(16, file='DERIV.DATA')
	    
	 do i=1,nobs
           read(16,*) xxd0                           !!!!unit=16:  DERIV.DATA
           read(16,*) xxd1
           read(16,*) xxd2,xxd3,perd(i),ccpred(i)
           read(16,*) xxd4
	   
 	   do j=1,np
	     read(16,*) depth(j),dcr(j,i),dca(j,i),dcb(j,i)
               g(i,j)= dcb(j,i)*thick(j)
	   enddo
	enddo
         close (16)

	  if (indaniso.eq.1) then
	    do i = 1,np
	      origmod(i) = stcs2th(i)
	    enddo
	  else
	    do i = 1,np
	      origmod(i) = stsn2th(i)
	    enddo
	  endif
c  calculate expected effects of starting model - linear
          do i = 1, nobs
	    gchange = 0.0
	    do jc = 1,np
	      gchange = gchange + g(i,jc)*origmod(jc)
	    enddo
	    if (indaniso.eq.1) then
	      ccpredout(i)= gchange
	    else
	      ccpredout2(i) = gchange
	    endif
	  enddo
         	    
       		
	
230   format (3x,f7.3,5x,f5.3,5x,f5.3,5x,f5.3,3x,f8.3)

         close(30)
        

        indaniso = indaniso +1
	if (indaniso.eq.2) go to 2000
	
	do i = 1, nper
c	  write(60+i,*) lon,lat,velsum1(i),velsum2(i)
	  write(60+i,*) lon,lat,ccpredout(i),ccpredout2(i)
	enddo


	enddo  ! close all grid points
	
	

	close(10)
	do i = 1, nper
	  close(60+i)
	enddo
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

