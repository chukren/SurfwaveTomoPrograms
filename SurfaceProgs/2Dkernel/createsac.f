	parameter( maxnt = 100000 )
	real tt(maxnt),u(maxnt)
	character *70 fn
        write(*,*) 'please input frequency'
	read(*,*)  freq0
	write(*,*) 'please input file name'
	read(*,*) fn
	  pi = 4.* atan(1.)
	  period0 = 1./freq0	
	  nt = 3201
	  dt = 0.5

	  do it = 1,nt
	      tt(it) = (it-1)*dt
	      u(it) = cos(2*pi*tt(it)/period0)

	  enddo  
             	

	 call newhdr
         call setnhv('npts',nt,nerr)
         call setfhv('b',tt(1),nerr)
         call setfhv('e',tt(nt),nerr)
         call setfhv('delta',dt,nerr)
     
        call wsac0(fn,tt,u,nerr)
	end