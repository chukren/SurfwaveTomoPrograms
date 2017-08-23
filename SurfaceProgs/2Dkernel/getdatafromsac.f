c  modified from YingJie. change maxnt=10000 to 1000, our intereted frequency
c  range if from .007 to .05. The frequency in the *.am file is from 0 to
c  1. so just reading in the subdata which contains our interested frequencies
c  to find out the first lobe.
  
	parameter( maxnt = 1000 )
	real u(maxnt)
	character *70 fn,outfn
 
    	write(*,*) 'please input input  file name'
	read(*,*) fn
	write(*,*) 'please inout output file name'
	read(*,*) outfn
	open(10,file = outfn)
             	
         call rsac1(fn,u,npts,beg,delt,maxnt,nerr)
         do i = 1, npts
	 write(10,*) beg+(i-1)*delt,u(i)	   
	 enddo
	 close(10)
	end
