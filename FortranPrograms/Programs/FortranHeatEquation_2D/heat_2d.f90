	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves heat equation in 2 dimension
	! u_t=\alpha*(u_{xx}+u_{yy}
	! using a the backward Euler method for x\in[0,2\pi]
	!
	! The boundary conditions are u(0)=u(2\pi)
	! The initial condition is u=sin(x)
	!
	! .. Parameters ..
	!  Nx	= number of modes in x - power of 2 for FFT
	!  Ny	= number of modes in y - power of 2 for FFT
	!  Nt	= number of timesteps to take
	!  Tmax	= maximum simulation time
	!  plotgap			= number of timesteps between plots
	!  FFTW_IN_PLACE 	= value for FFTW input 
	!  FFTW_MEASURE 	= value for FFTW input
	!  FFTW_EXHAUSTIVE 	= value for FFTW input
	!  FFTW_PATIENT 	= value for FFTW input    
	!  FFTW_ESTIMATE 	= value for FFTW input
	!  FFTW_FORWARD     	= value for FFTW input
	!  FFTW_BACKWARD	= value for FFTW input	
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  Lx				= width of box 
	!  Ly				= width of box 
	!  alpha			= heat conductivity
	! .. Scalars ..
	!  ix				= loop counter in x direction
	!  iy				= loop counter in y direction
	!  it				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfx			= Forward 2d fft plan in x
	!  planbx			= Backward 1d fft plan in x
	!  dt				= timestep
	! .. Arrays ..
	!  u				= approximate REAL solution
	!  v				= Fourier transform of approximate solution
	!  vna 				= temporary field
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
	!
	! REFERENCES
	!
	! ACKNOWLEDGEMENTS
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	! Check that the initial iterate is consistent with the 
	! boundary conditions for the domain specified
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
				
	PROGRAM main
				 	   
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4),	PARAMETER 	::  Nx=64,Ny=64 
	INTEGER(kind=4),	PARAMETER	::  Nt=20  
	REAL(kind=8),	PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8),	PARAMETER	::  Lx=5.0d0,Ly=5.0d0			 
	REAL(kind=8),	PARAMETER	::  alpha=0.50d0	
	REAL(kind=8)	::  dt=0.2d0/REAL(Nt,kind(0d0))		
	COMPLEX(KIND=8), DIMENSION(:),ALLOCATABLE	::  kx,ky,kx2,ky2	
	REAL(kind=8), DIMENSION(:),ALLOCATABLE	::  x,y	
	COMPLEX(KIND=8), DIMENSION(:,:),ALLOCATABLE	::  u,v	
	REAL(kind=8), DIMENSION(:),ALLOCATABLE	::  time
	COMPLEX(KIND=8), DIMENSION(:,:),ALLOCATABLE	::  vna	
	INTEGER(kind=4)	::  ix,iy,k,it
	INTEGER(kind=4)	:: start, finish, count_rate, AllocateStatus
	INTEGER(kind=4), PARAMETER	:: FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64
	INTEGER(kind=4), PARAMETER :: FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	COMPLEX(KIND=8), DIMENSION(:,:),ALLOCATABLE :: fftfxy,fftbxy
	INTEGER(kind=8)	:: planfxy,planbxy
	CHARACTER*100	:: name_config
	    		
	CALL system_clock(start,count_rate)	
	ALLOCATE(kx(1:Nx),kx2(1:Nx),x(1:Nx),ky(1:Ny),ky2(1:Nx),y(1:Ny),u(1:Nx,1:Ny),v(1:Nx,1:Ny),&
   			time(1:1+Nt),vna(1:Nx,1:Ny),fftfxy(1:Nx,1:Ny),fftbxy(1:Nx,1:Ny),&
   			stat=AllocateStatus)
   	IF (AllocateStatus .ne. 0) STOP 

		! set up ffts
	CALL dfftw_plan_dft_2d(planfxy,Nx,Ny,fftfxy(1:Nx,1:Ny),fftbxy(1:Nx,1:Ny),&
			FFTW_FORWARD,FFTW_ESTIMATE)
	CALL dfftw_plan_dft_2d(planbxy,Nx,Ny,fftbxy(1:Nx,1:Ny),fftfxy(1:Nx,1:Ny),&
			FFTW_BACKWARD,FFTW_ESTIMATE)
	
	
	PRINT *,'Setup FFTs'
		
	! setup fourier frequencies
	DO ix=1,1+Nx/2
		kx(ix)= cmplx(0.0d0,1.0d0)*REAL(ix-1,kind(0d0))/Lx  			
	END DO
	kx(1+Nx/2)=0.00d0
	DO ix = 1,Nx/2 -1
		kx(ix+1+Nx/2)=-kx(1-ix+Nx/2)
	END DO
	kx2=kx*kx
	DO ix=1,Nx
		x(ix)=(-1.00d0 + 2.00d0*REAL(ix-1,kind(0d0))/REAL(Nx,KIND(0d0)))*pi*Lx
	END DO

	DO iy=1,1+Ny/2
		ky(iy)= cmplx(0.0d0,1.0d0)*REAL(iy-1,kind(0d0))/Ly  			
	END DO
	ky(1+Ny/2)=0.00d0
	DO iy = 1,Ny/2 -1
		ky(iy+1+Ny/2)=-ky(1-iy+Ny/2)
	END DO
	ky2=ky*ky
	DO iy=1,Ny
		y(iy)=(-1.00d0 + 2.00d0*REAL(iy-1,kind(0d0))/REAL(Ny,KIND(0d0)))*pi*Ly
	END DO		
	PRINT *,'Setup grid and fourier frequencies and splitting coefficients'

	
	!u(1:Nx,1:Ny)=sin(x(1:Nx)) 
        DO ix=1,Nx
		DO iy=1,Ny
			u(ix,iy)=sin(x(ix))
		END DO
	END DO
		! transform initial data
	CALL dfftw_execute_dft_(planfxy,u,v) 
	PRINT *,'Got initial data, starting timestepping'
	time(1)=0.0d0

	vna=v
	PRINT *,'Starting timestepping'
	DO it=1,Nt
		DO ix=1,Nx
			DO iy=1,Ny
			vna(ix,iy)=vna(ix,iy)/(1-dt*(kx2(ix)+ky2(iy)))
			end do
		END DO
		PRINT *,'storing plot data ',it
		time(it+1)=time(it)+dt
		v=vna
		CALL dfftw_execute_dft_(planbxy,v,u)
		u=u/REAL(Nx*Ny,KIND(0d0))	! normalize
	END DO
	PRINT *,'Finished time stepping'
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),'for execution'

	! Write data out to disk	

	name_config = 'u.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO iy=1,Ny
		DO ix=1,Nx
			WRITE(11,*) REAL(u(ix,iy))
		END DO
	END DO
	CLOSE(11)
		
	!name_config = 'tdata.dat' 
	!OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	!REWIND(11)
	!DO j=1,1+Nt
	!	WRITE(11,*) time(j)
	!END DO
	!CLOSE(11)

	name_config = 'xcoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO ix=1,Nx
		WRITE(11,*) x(ix)
	END DO
	CLOSE(11)

	name_config = 'ycoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO iy=1,Ny
		WRITE(11,*) y(iy)
	END DO
	CLOSE(11)



	PRINT *,'Saved data'
	DEALLOCATE(kx,ky,kx2,ky2,x,u,v,&
   			time,vna,fftfxy,fftbxy,&
   			stat=AllocateStatus)
   	IF (AllocateStatus .ne. 0) STOP 
	
	CALL dfftw_destroy_plan(planbxy)
	CALL dfftw_destroy_plan(planfxy)
	CALL dfftw_cleanup()
	PRINT *,'Program execution complete'
	END PROGRAM main

