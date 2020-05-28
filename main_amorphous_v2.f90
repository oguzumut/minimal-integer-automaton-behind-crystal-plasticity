INCLUDE 'mkl_vsl.f90'
PROGRAM phase_field

!*********************************************************************************************
!*                                                                                           *
!*                                      Main  starts here                                    *
!*                                      coded by O.U. Salman and P. Zhang                    *
!*                                                                                           *
!*********************************************************************************************

USE MKL_VSL_TYPE
USE MKL_VSL
use parameters
use potential_arrays
use fourier_operators_discrete
use threshholds
IMPLICIT NONE

integer (kind=4) :: i, j, k, n_dim_calcul, ind, ind_u, no_config, N, cycle_set, disorder_method, initial_dislocation_switch
integer (kind=4) :: iter, i1, j1, k1, i2, j2, i3, j3, mark_ava_n, cycle_number, sign_load, number_stable,ovito_step, moving_distance
integer (kind=4) :: nucleation, annihilation, motion, temp_i, mark_pro, dislocation_window
integer (kind=4), dimension(:),  allocatable :: maxlocations,minlocations
integer (kind=4), dimension(:,:,:), allocatable :: work_space_property


real(kind=8) :: dx,gamma,a,b,load, temp_load, sum_av,sum_temp,sum_s,keeper,sum_d,stiff,zero,one,two, zero_fake
real(kind=8) :: load_increase, load_rate, rho_initial, load_upper, load_lower, fxx_mean, fxx_sigma, fxy_mean, fxy_sigma
real(kind=8) :: hxx_mean, hxy_mean, hxx_sigma, hxy_sigma, disorder_fraction_set, disorder_fraction_real
real(kind=8) :: disorder_size_mean, disorder_size_sigma, dist_mean, dist_sigma, fraction_threshold,f_mean, f_sigma
real(kind=8) :: energy_loading, energy_stable, energy_release, sxy_loading, sxy_stable, sxy_release
real(kind=8), dimension(:,:,:),  allocatable ::  work_space_ex, work_space_well, exy_old
real(kind=8), dimension(:),  allocatable :: work_space2


complex (kind=8), dimension(:,:,:),  allocatable :: tf_work_space

logical test_roule

character(len=20), save :: fft_name
character(len=100), save :: name_config,name_config_u
character (len=100) :: sauvegarde
character (len=5), save :: tamp_config


! ------------------------ Parameter Setting --------------------------------

!********* Parameters required to set by the users ***************
nx=256; ny=256;  nz=1; n_dim_calcul = 2  	  !  dimensions of system
stiff = 0.5                 							  !  flexibility of system
load_upper=0.4; load_lower=-0.4                    !  upper and lower bound of loading strains
a = -0.25; b= 0.25	                                            !  rang of uniform distribution
cycle_set = 10       							  !  number of cycle loading
ovito_step= 1000						        ! write ovito file
initial_dislocation_switch=0                             ! value 0: no initial dislocation  value 1: initial dislocation
rho_initial = 0.008	     						        !  initial dislocation density 

! disorder parameters
disorder_method = 1
! disorder_method=1 (internal stress at every points)
fxx_mean = 0.0 ; fxx_sigma = 0.5		!  mean and sigma of gaussian distribution
fxy_mean = 0.0 ; fxy_sigma = 0.1
! disorder_method=2 (internal stress at circle inclusions)
hxx_mean = 0.0 ; hxx_sigma = 0.15 			  !  mean and sigma of internal normal stress
hxy_mean = 0.0 ; hxy_sigma = 0.6       	      !  mean and sigma of shear stress
disorder_fraction_set = 8					  !  setting disorder fraction
disorder_size_mean = 3.0					  !  the disorder radius
disorder_size_sigma =1.0					  !  the standard deviation of disorder radius
! disorder method=3 (trapping effect, note that you may also need to set the parameters in disorder_method=1)
dist_mean=0.2  ; dist_sigma=0.1  
fraction_threshold= 1
! disorder method=4 (add themal fluctuations, note that you may also need to set the parameters in disorder_method=1)
f_mean=0.0; f_sigma=0.1


!**********************************************************

! Initial state of parameters
N=nx*ny*nz 
iter = -1
ind = 0
sign_load = 1
number_stable = -1
no_config=10000
mark_ava_n = 0
fft_name = 'Intel'
sauvegarde = 'complete'
test_roule = .TRUE.
sum_av = 0; sum_s = 0; sum_d = 0
zero = 0.0; one =1.0; two = 2.0; zero_fake=1e-9
load = 0.0; load_increase=0.0; dislocation_window=0
cycle_number = 0    
moving_distance = 0   	
nucleation=0; annihilation=0; motion=0							  
pi = 3.1415926535897932384626433832795028841971693
pi2 = 2.*3.1415926535897932384626433832795028841971693

! Open files to write informations and results	
open(1,file='test.dat',     form='formatted',status='unknown')                            			!  use for debug
open(2,file='dislocation density_eachstep.dat',     form='formatted',status='unknown')
open(3,file='moving_distance.dat',     form='formatted',status='unknown')	
open(4,file='moving_number.dat',     form='formatted',status='unknown')		
open(20,file='energy_dissip.dat',     form='formatted',status='unknown')
open(23,file='strain_stress.dat',     form='formatted',status='unknown')
open(24,file='energy_dissip_pdirection.dat',     form='formatted',status='unknown')
open(25,file='energy_dissip_ndirection.dat',     form='formatted',status='unknown')
open(30,file='information_tailles', form='formatted',status='unknown')
open(40,file='information_tailles_Eps', form='formatted',status='unknown')



write(4, 1011) 'nucleation_number', 'annihilation_number', 'moving_number', 'total_number'
write(30,1000) nx,ny,nz, sauvegarde
write(30,1002) 3, 0; 
write(40,1000) nx,ny,nz, 3; close (40)
write(20, 1009) 'update_number', 'exy','sxy_after_loading', 'sxy_after_release', 'sxy_drop', 'energy_drop'
write(24, 1009) 'update_number', 'exy','sxy_after_loading', 'sxy_after_release', 'sxy_drop', 'energy_drop'
write(25, 1009) 'update_number', 'exy','sxy_after_loading', 'sxy_after_release', 'sxy_drop', 'energy_drop'

!---------------------------- ALLOCATES TABLES  -------------------------------------------	

allocate (work_space(0:nx-1, 0:ny-1, 0:nz-1))
allocate (tf_work_space(0:nx/2, 0:ny-1, 0:nz-1))
allocate (work_space_ex(0:nx-1, 0:ny-1, 0:nz-1))
allocate (work_space_well(0:nx-1, 0:ny-1, 0:nz-1))
allocate(work_space_property(0:nx-1, 0:ny-1, 0:nz-1))
allocate (exy_old (0:nx-1, 0:ny-1, 0:nz-1))
allocate(work_space2(1))
allocate (maxlocations(1:3))
allocate (minlocations(1:3))


CALL allocate_fourier_operators_discrete  ! It initriates the fourier operators
CALL initiate_fourier_operators_discrete   ! It calculates the fourier operators
CALL allocate_fields ! It initriates strains and displacements in real and fourier space
CALL allocate_threshholds ! it initiates 2 real space tables (dist_uniform and dist_gauss)
CALL MKL_VSL_UNIFORM  (N, dist_uniform, a, b) ! Returns uniform distribution between a and  b real
!CALL MKL_VSL_GAUSSIAN(N, dist_gauss, f_mean, f_sigma) ! Returns gaussian distribution with (f_mean, f_sigma)

! Initial state of fields
d_field(:,:,:) = 0; sxy_loading=0.0; energy_loading=0.0
work_space_ex(:,:,:) = 0.0; work_space_well (:,:,:)=0.0
exx(:,:,:) =  0; exy(:,:,:) =  dist_uniform(:,:,:); threshold(:,:,:) = 0.5

! Generate initial dislocations
if (initial_dislocation_switch==1) then
	call dislocation_generator (rho_initial, d_field(0,0,0))
endif

! Generate disorders
select case(disorder_method)
	case (1)
		! Generate disorders (paper version)	
		call MKL_VSL_GAUSSIAN(N, hxx, fxx_mean, fxx_sigma) ! Returns gaussian distribution with (f_mean, f_sigma)
		call MKL_VSL_GAUSSIAN(N, hxy, fxy_mean, fxy_sigma)   
		write(30, 1008) "fxx_mean=", fxx_mean, "fxx_sigma=", fxx_sigma, "fxy_mean=", &
						   fxy_mean, "fxy_sigma=", fxy_sigma
		close (30)
	case (2)
		! Generate disorders (Umut version)	
		CALL disorder_generator (disorder_fraction_set, disorder_fraction_real, disorder_size_mean, &
									disorder_size_sigma, hxx_mean, hxx_sigma, hxy_mean, hxy_sigma, &
									hxx(0,0,0), hxy(0,0,0))
		write(30, 1006)  disorder_fraction_real
		write(30, 1007) "hxx_mean=", hxx_mean, "hxx_sigma=", hxx_sigma, "hxy_mean=", &
							   hxy_mean, "hxy_sigma=", hxy_sigma, "disorder_size_mean=", disorder_size_mean, &
							   "disorder_size_sigma=", disorder_size_sigma
		close (30)
	case (3)
		! Generate disorders (local well trapping)	
		call MKL_VSL_GAUSSIAN(N, hxx, fxx_mean, fxx_sigma) ! Returns gaussian distribution with (f_mean, f_sigma)
		call MKL_VSL_GAUSSIAN(N, hxy, fxy_mean, fxy_sigma)   
		write(30, 1008) "fxx_mean=", fxx_mean, "fxx_sigma=", fxx_sigma, "fxy_mean=", &
						   fxy_mean, "fxy_sigma=", fxy_sigma
		close (30)
		! Generate local strong well trapping
		Call threshold_generator (dist_gauss, fraction_threshold, dist_mean, dist_sigma)
		threshold(:,:,:) = 0.5+dist_gauss
	case (4) 
		! Themal fluctuations
		call MKL_VSL_GAUSSIAN(N, hxx, fxx_mean, fxx_sigma) ! Returns gaussian distribution with (f_mean, f_sigma)
		call MKL_VSL_GAUSSIAN(N, hxy, fxy_mean, fxy_sigma)   
		write(30, 1008) "fxx_mean=", fxx_mean, "fxx_sigma=", fxx_sigma, "fxy_mean=", &
						   fxy_mean, "fxy_sigma=", fxy_sigma
		close (30)
		CALL MKL_VSL_GAUSSIAN(N, dist_gauss, f_mean, f_sigma)
end select

!---------------------------- CALCULATION STARTS --------------------------------

write(*,*) "The dimension of the system", nx, ny, nz
write(*,*) "calculation starts"
CALL r2c(nx, ny, nz, fft_name, hxx(0,0,0), tf_hxx(0,0,0), n_dim_calcul)
CALL r2c(nx, ny, nz, fft_name, hxy(0,0,0), tf_hxy(0,0,0), n_dim_calcul) 

do while (cycle_number <= cycle_set) 

! Part 1: d_field -----> tf_d_field ----> tf_u_field ----> tf_exy (0) ----> tf_exy (load) ----> exy (load)

	call r2c(nx,ny,nz,fft_name,d_field(0,0,0),tf_d_field(0,0,0),n_dim_calcul) 
	tf_u_field(:,:,0) = (stiff*qyn(:,:,0)*tf_d_field(:,:,0)+stiff*(qxn(:,:,0)*tf_hxx(:,:,0)+qyn(:,:,0)*tf_hxy(:,:,0)))/&
					     (qxx(:,:,0)+stiff*qyy(:,:,0)+tol(:,:,0))  
	tf_exy(:,:,0) = (qyp(:,:,0)*tf_u_field(:,:,0))	
	tf_exy(0,0,0) = load
	call c2r(nx,ny,nz,fft_name,tf_exy(0,0,0),exy(0,0,0),n_dim_calcul)
	
	! Calculate the dislocation density and position
	tf_work_space(:,:,:) = qxp(:,:,:)*tf_d_field(:,:,:)  
	call c2r(nx,ny,nz,fft_name,tf_work_space(0,0,0),work_space(0,0,0),n_dim_calcul) 
	write(2,1004) sum(abs(work_space(:,:,:)))/(1.0*nx*ny*nz)

	
	
! Part 2:  Cycle loading  (If every atoms are inside the wells, then increase/decrease the shear strain. Otherwise, update 
!            				     d_filed untile all of the atomes are inside the energy wells.)
	
	! Impose the themal fluctuations to the strain if disorder_method is 4;  
	! Attention: you can also change the dist_gauss(:,:,:) every step by using the "CALL MKL_VSL_GAUSSIAN(N, dist_gauss, f_mean, f_sigma)"
	if (disorder_method==4) then
		exy_fluc(:,:,:)=exy(:,:,:)+dist_gauss(:,:,:)*sign_load	
	else
		exy_fluc(:,:,:) = exy(:,:,:)
	endif
	
	! work_space_well (i,j,0) > 0, the atom is exceed the threshold. Vise versa.  
	work_space_well(:,:,:) = abs(exy_fluc(:,:,:)-d_field(:,:,:))-threshold(:,:,:)
	if (maxval(work_space_well(:,:,:)) < 0.0) then    	 
			
	    ! (i): write the strain and energy at the stable configuration		
			number_stable = number_stable+1
			exy_homo(:,:,:) = exy(:,:,:)-d_field(:,:,:)       ! calculate shear strain exy_homo
			tf_exx(:,:,0) = (qxp(:,:,0)*tf_u_field(:,:,0))   ! calculate tension strain exx
			call c2r(nx,ny,nz,fft_name,tf_exx(0,0,0),exx(0,0,0),n_dim_calcul)
			sxy_stable = sum(exy_homo(:,:,:))/(1.0*nx*ny*nz) ! calculate the shear stress
			energy_stable = sum(0.5/stiff*exx(:,:,0)*exx(:,:,0))+sum(0.5*exy_homo(:,:,0)*exy_homo(:,:,0))-&
			                           sum(hxx(:,:,0)*exx(:,:,0)+hxy(:,:,0)*exy(:,:,0)) ! calculate the energy
			sxy_release = (sxy_loading-sxy_stable)*sign_load ! calculate the shear stress released during the stablization
			energy_release = energy_loading-energy_stable ! calculate the energy released during the stabilization
			! write the dissipation to three files: entire loading process, tension process, compression process
			if (number_stable > 0.5) then 
				write(20,1005)	 iter, load, sxy_loading, sxy_stable, sxy_release, energy_release
				if (sign_load > 0) then 
					write(24,1005)	 iter, load, sxy_loading, sxy_stable, sxy_release, energy_release
				else
					write(25,1005)	 iter, load, sxy_loading, sxy_stable, sxy_release, energy_release
				endif 
			endif	
			! write the stress-strain and dislocation density to file
!            dislocation_window=0
!            do i=113, 143
!                do j=113, 143
!                    dislocation_window=dislocation_window+abs(work_space(i,j,0))
!                enddo
!            enddo

            write(23,1004) load, sxy_stable, sum(abs(work_space(:,:,:)))/(1.0*nx*ny*nz)
			! write the moving dislocation density number to the file
			write(3,1010) moving_distance*1.0	
			! The dislocation position at the stable state
			work_space_property(:,:,:)=NINT(abs(work_space(:,:,:)))
			! write the number of generation, annihilation, moving dislocation
			write (4, 2000) nucleation, annihilation, motion, (nucleation+annihilation+motion)
			! set the number to zero
			nucleation=0; annihilation=0; motion=0
			
		! (ii): increase the shear strain
			work_space_ex(:,:,:) = sign_load*(exy_fluc(:,:,:)-d_field(:,:,:)-sign_load*threshold(:,:,:)) ! (sign_load= 1 or -1, it is the loading direction)
			maxlocations(:) = maxloc(work_space_ex(:,:,:))
			i1 = maxlocations(1)-1
			j1 = maxlocations(2)-1			
			load_increase = -(work_space_ex(i1,j1,0)/sign_load)+zero_fake*sign_load
			load = load + load_increase
			! judge the change of d_filed is a nucleation, annihilation or motion of dislocations
			i2=mod(nx+i1-1,nx); i3=mod(nx+i1+1,nx)
			j2=j1; j3=j1
			mark_pro=abs(d_field(i2,j2,0)-d_field(i1,j1,0))+abs(d_field(i3,j3,0)-d_field(i1,j1,0))
			d_field(i1,j1,0) = d_field(i1,j1,0)+sign_load
			mark_pro=NINT(mark_pro-abs(d_field(i2,j2,0)-d_field(i1,j1,0))-abs(d_field(i3,j3,0)-d_field(i1,j1,0)))


			select case (mark_pro)
				case (2)
					nucleation=nucleation+1
				case (-2)
					annihilation=annihilation+1
				case (0)
					if ( (abs(d_field(i2,j2,0)-d_field(i1,j1,0))) < (abs(d_field(i2,j2,0)-d_field(i1,j1,0)+sign_load)) ) then
						if ( work_space_property(i2,j2,0) /= 0 ) then
							motion=motion+1
							work_space_property(i2,j2,0)=work_space_property(i2,j2,0)-1
						endif		
					else
						if ( work_space_property(i1,j1,0) /= 0 ) then
							motion=motion+1
							work_space_property(i1,j1,0)=work_space_property(i1,j1,0)-1
						endif		
					endif	
				case default
					stop
			end select

				
			! If shear strain exceeds the setting boundaries, change loading direction
			if ((load < load_lower ) .OR. (load > load_upper)) then  
				sign_load = sign_load*(-1)
				cycle_number = cycle_number+1
			endif
			
		! (iii): calculate the energe after loading				
			exy_homo(:,:,:) = exy_homo(:,:,:) + load_increase
			exy(:,:,:) = exy(:,:,:) + load_increase
			sxy_loading = sum(exy_homo(:,:,:))/(1.0*nx*ny*nz)
			energy_loading = sum(0.5/stiff*exx(:,:,0)*exx(:,:,0))+sum(0.5*exy_homo(:,:,0)*exy_homo(:,:,0))-&
							     sum(hxx(:,:,0)*exx(:,:,0)+hxy(:,:,0)*exy(:,:,0))    ! calculate the energy
			iter= 0	
			moving_distance=1 ! set the number of moving dislocations after loading


 write(23,1004) load, sxy_loading, sum(abs(work_space(:,:,:)))/(1.0*nx*ny*nz)
			
	else
           do i=0, nx-1
               do j=0, ny-1
                   if (work_space_well(i,j,0) > 0.0) then
                       i2=mod(nx+i-1,nx); i3=mod(nx+i+1,nx)
                       j2=j; j3=j
                       mark_pro=abs(d_field(i2,j2,0)-d_field(i,j,0))+abs(d_field(i3,j3,0)-d_field(i,j,0))
                       if ((exy_fluc(i,j,0)-d_field(i,j,0)) > zero) then
                           d_field(i,j,0) = d_field(i,j,0)+1; temp_i=1
                       else
                           d_field(i,j,0) = d_field(i,j,0)-1; temp_i=-1
                       endif
                       mark_pro=NINT(mark_pro-abs(d_field(i2,j2,0)-d_field(i,j,0))-abs(d_field(i3,j3,0)-d_field(i,j,0)))
                       select case (mark_pro)
                           case (2)
                               nucleation=nucleation+1
                           case (-2)
                               annihilation=annihilation+1
                           case (0)
                               if ( (abs(d_field(i2,j2,0)-d_field(i,j,0))) < (abs(d_field(i2,j2,0)-d_field(i,j,0)+temp_i)) ) then
                                   if ( work_space_property(i2,j2,0) /= 0 ) then
                                       motion=motion+1
                                       work_space_property(i2,j2,0)=work_space_property(i2,j2,0)-1
                                   endif
                               else
                                   if ( work_space_property(i,j,0) /= 0 ) then
                                       motion=motion+1
                                       work_space_property(i,j,0)=work_space_property(i,j,0)-1
                                   endif
                               endif
                       end select
                       moving_distance=moving_distance+1
                   endif
               enddo
           enddo


 !           d_field(:,:,:) = NINT(exy_fluc(:,:,:))
			iter=iter+1
			cycle
	endif

	! Part 3: Output the ovito file when increase the load for "ovito_step" times
	 write (*,*) cycle_number, number_stable
	 if (mod(number_stable,ovito_step)==0) then
		! Path of ovito file
		name_config = 'dir_config/config_'
		name_config_u = 'dir_config_u/config_'
		ind = index(name_config,' ')-1
		ind_u = index(name_config_u,' ')-1	
		write(tamp_config,'(i5)') no_config
		name_config = name_config(1:ind)//tamp_config
		name_config_u = name_config_u(1:ind_u)//tamp_config	
		! write data
		call c2r(nx,ny,nz,fft_name,tf_u_field(0,0,0),u_field(0,0,0),n_dim_calcul)  ! calculate u_filed
		exy_homo(:,:,:) = exy_homo(:,:,:)-load_increase
		exy(:,:,:) = exy(:,:,:) -load_increase
		!exy_homo(:,:,:) = exy(:,:,:)-d_field(:,:,:)														   ! calculate real strain exy_homo      
		call write_ovito (exx,exy,exy_homo,u_field,work_space, name_config,name_config_u, load)       ! write ovito file
		no_config = no_config+1  
	end if

 enddo
 
1000 format(3(1x,i4),5x,'mode sauvegarde: ',a10)
2000 format(4(1x,i8))
1002 format(2(1x,i3))
1003 format(e14.7,(2x),e14.7)
1004 format(e14.7,(2x),e14.7, (2x),e14.7)
1005 format(i4, (2x), e14.7, (2x), e14.7, (2x), e14.7, (2x), e14.7, (2x), e14.7)
1006 format ('volume fraction:', e14.7)
1007 format (6(A, e14.7))
1008 format (4(A,e14.7))
1009 format (6(A,(2x)))
1010 format (e14.7)
1011 format (4(A,(2x)))
5000 continue

END PROGRAM phase_field




!***********************************  SUBROUTINE  *****************************************


SUBROUTINE allocate_fourier_operators_discrete()
	USE fourier_operators_discrete
	USE parameters
	IMPLICIT NONE
	! Declare local variables
	INTEGER :: AllocateStatus,i

	! Allocate storage for array A accordingly
	ALLOCATE(qxx(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(qyy(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(qzz(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(qxp(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(qyp(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(qxn(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(qyn(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	! 	ALLOCATE( qzp(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	! 	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(sqr_norme_q(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tol(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
END SUBROUTINE allocate_fourier_operators_discrete




SUBROUTINE allocate_fields()
	USE potential_arrays
	USE parameters
	IMPLICIT NONE
	! Declare local variables
	INTEGER :: AllocateStatus,i
	
    ! Allocate storage for array A accordingly
	ALLOCATE(exx(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(exy(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(exy_homo(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(exy_fluc(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(d_field(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(u_field(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
! 	ALLOCATE(stiff_field(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
! 	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(hxx(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(hxy(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tf_exx(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tf_exy(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tf_d_field(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tf_u_field(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tf_hxx(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE(tf_hxy(0:nx/2, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
END SUBROUTINE allocate_fields




SUBROUTINE initiate_fourier_operators_discrete()
	USE parameters
	USE fourier_operators_discrete
	IMPLICIT NONE
	
	INTEGER :: i,j,k
	REAL :: xx,yy,zz

	do k = 0, nz-1
		do j = 0, ny-1
			do i = 0, nx/2	 
				qxx(i,j,k) = 2*(cos(pi2*i/nx)-1.0)
				qyy(i,j,k) = 2*(cos(pi2*j/ny)-1.0)
				qzz(i,j,k) = 2*(cos(pi2*k/nz)-1.0)
				sqr_norme_q(i,j,k) = (xx + yy +zz ) 		
			end do
		end do
	end do

	do k = 0, nz-1
		do j = 0, ny-1
			do i = 0, nx/2			 
				qxp(i,j,k) = (cos(pi2*i/nx)+sin(pi2*i/nx)*cmplx(0,1)-1.0)
				qyp(i,j,k) = (cos(pi2*j/ny)+sin(pi2*j/ny)*cmplx(0,1)-1.0)	 
				qxn(i,j,k) = 1.0-(cos(pi2 *i/nx)-sin(pi2*i/nx)*cmplx(0,1))
				qyn(i,j,k) = 1.0-(cos(pi2 *j/ny)-sin(pi2*j/ny)*cmplx(0,1))
			end do
		end do
	end do
		   
	tol(:,:,:) = 0
	tol(0,0,:) = 1.
END SUBROUTINE initiate_fourier_operators_discrete



SUBROUTINE allocate_threshholds()
	USE threshholds
	USE parameters
	IMPLICIT NONE
	! Declare local variables
	INTEGER :: AllocateStatus
	
	ALLOCATE(dist_gauss(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE( dist_uniform(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	ALLOCATE( threshold(0:nx-1, 0:ny-1, 0:nz-1), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
END SUBROUTINE allocate_threshholds



SUBROUTINE write_ovito (exx,exy,exy2,u,d,name_config,name_config_u, load)         	
	USE parameters
	IMPLICIT NONE
	
	real (kind=8), intent(in), dimension(0:nx-1,0:ny-1,0:nz-1) :: exx,exy,exy2,u,d
	character*100, intent(in) :: name_config,name_config_u
	real(kind=8),   intent(in):: load
	integer (kind=4) :: np,i,j,k

	open(60,file=name_config,form='formatted',status='unknown')
	write(60,1000) nx*ny*nz
	write(60,*) " "
	do i=0,nx-1
		do j=0,ny-1	
			!write(60,1010) 1.0*i+u(i,j,0)+load*j,1.0*j, exx(i,j,0),exy(i,j,0),exy2(i,j,0), d(i,j,0)
            write(60,1010) 1.0*i+u(i,j,0),1.0*j, exx(i,j,0),exy(i,j,0),exy2(i,j,0), d(i,j,0)
		end do
	end do
	close(60)
	
	1000 format(i8)
	1010 format(e14.7,(2x),e14.7,(2x),e14.7,(2x),e14.7,(2x),e14.7,(2x),e14.7)

	RETURN
END SUBROUTINE write_ovito
  
  
  
  
SUBROUTINE MKL_VSL_UNIFORM(n,r,a,b)
      USE MKL_VSL_TYPE
      USE MKL_VSL
	  IMPLICIT NONE
		
      integer (kind=4) clock
      integer (kind=4) errcode
      integer(kind=4),save:: startu=1
      integer brng,method,seed
      integer(kind=4),  intent(in) :: n
      
	  real(kind=8),  intent(in) ::a,b
      real(kind=8),  intent(inout),  dimension(1:N) :: r

      TYPE (VSL_STREAM_STATE), save :: stream
      
       brng=VSL_BRNG_MCG31
       method=VSL_RNG_METHOD_UNIFORM_STD

		!if (start==0) then
			CALL SYSTEM_CLOCK(COUNT=clock)
			startu = startu +2 + clock
			seed = startu
			errcode=vslnewstream(stream, brng, seed)
		!end if
       errcode=vdrnguniform(method, stream, n, r, a, b )     
END SUBROUTINE MKL_VSL_UNIFORM


SUBROUTINE MKL_VSL_GAUSSIAN (n, r, f_mean, f_sigma)
      USE MKL_VSL_TYPE
      USE MKL_VSL
	  IMPLICIT NONE
		
      integer (kind=4) clock
      integer (kind=4) errcode
      integer(kind=4), save:: startg=1
      integer brng,method,seed
      integer(kind=4),  intent(in) :: n      
	  real(kind=8),  intent(in) :: f_mean, f_sigma
      real(kind=8),  intent(inout),  dimension(1:N) :: r
      TYPE (VSL_STREAM_STATE), save :: stream
      
       brng = VSL_BRNG_MCG31
       method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
		!if (start==0) then
			CALL SYSTEM_CLOCK(COUNT=clock)
			startg = startg+clock+23455677
			seed = startg
			errcode=vslnewstream(stream, brng, seed)
		!end if
       errcode=vdrnggaussian(method, stream, n, r, f_mean, f_sigma)     
END SUBROUTINE MKL_VSL_GAUSSIAN


SUBROUTINE  dislocation_generator (rho, d_field)
	USE parameters 
	IMPLICIT NONE

	real(kind=8),  intent(in) :: rho
	real(kind=8),  intent(inout),  dimension(0:nx-1, 0:ny-1, 0:nz-1) :: d_field
	
	integer (kind=4) :: N_dislo_initial, pxb, pyb, pxe, sign, temp, i, j
	real(kind=8) :: zero, one, two
	real(kind=8), dimension(0:1) :: po_ne
      real(kind=8), dimension(:),  allocatable :: dislo_pxb, dislo_pyb, dislo_len, dislo_sign
    
      zero = 0.0
      one = 1.0
      two = 2.0
 	N_dislo_initial=floor(rho/2.0*nx*ny*nz)

	allocate (dislo_pxb (0:N_dislo_initial-1))
 	allocate (dislo_pyb (0:N_dislo_initial-1))
 	allocate (dislo_len  (0:N_dislo_initial-1))
 	allocate (dislo_sign(0:N_dislo_initial-1))
 	po_ne(0) = -1.0
 	po_ne(1) =  1.0 	
 	
 	call MKL_VSL_UNIFORM(N_dislo_initial, dislo_pxb,zero,real(nx*1.0,kind=8))
 	call MKL_VSL_UNIFORM(N_dislo_initial, dislo_pyb,zero,real(ny*1.0,kind=8)) 
 	call MKL_VSL_UNIFORM(N_dislo_initial, dislo_len,one,  real(nx*1.0,kind=8))
 	call MKL_VSL_UNIFORM(N_dislo_initial, dislo_sign,zero,two)
 	
 	do i = 0, N_dislo_initial-1
 		pxb = floor(dislo_pxb(i))
 		pyb = floor(dislo_pyb(i))
		pxe = pxb+floor(dislo_len(i))
 		sign = po_ne(floor(dislo_sign(i)))	
 		do j = pxb, pxe
 			temp = mod(j, nx+1)
 			d_field (temp,pyb,0) = d_field (temp,pyb,0)+sign
 		enddo	
 	enddo
 	deallocate (dislo_pxb, dislo_pyb, dislo_len, dislo_sign)
END SUBROUTINE 


SUBROUTINE  disorder_generator (fraction_set, fraction_real, size_mean, size_sigma, &
                           hxx_mean, hxx_sigma, hxy_mean, hxy_sigma, hxx, hxy)
	USE parameters 
	IMPLICIT NONE

	real(kind=8),  intent(in) :: fraction_set, size_mean, size_sigma, hxx_mean, hxx_sigma, hxy_mean, hxy_sigma
	real(kind=8),  intent(inout) :: fraction_real
	real(kind=8),  intent(inout),  dimension(0:nx-1, 0:ny-1, 0:nz-1) :: hxx, hxy
	
	integer(kind=4) :: i, j, k, m, ii, jj, n_position, disorder_number, i_range, j_range
	real(kind=8) :: zero, one, temp,size_mean_normal, size_sigma_normal
	integer(kind=8), dimension(:),  allocatable :: ip_x, ip_y, diameter, disorder_index
      real(kind=8), dimension(:),  allocatable :: p_x, p_y, diameter_normal,  d_hxx, d_hxy
      logical :: boundary, flag
    
	fraction_real= 0.0
	zero=0.0
	one=1.0
	disorder_number = 0
	flag= .true. 
	size_mean_normal = dlog(size_mean)-0.5*dlog(1.0+size_sigma/size_mean/size_mean)
	size_sigma_normal = (dlog(1.0+size_sigma/size_mean/size_mean))**0.5
	n_position = floor (fraction_set*nx*ny/pi/(size_mean*size_mean))
	write(*,*) n_position
 	
 	allocate (p_x(0:n_position-1))		
 	allocate (p_y(0:n_position-1))
	allocate (ip_x(0:n_position-1)) 			
	allocate (ip_y(0:n_position-1))
	allocate (diameter(0:n_position-1))
	allocate (diameter_normal(0:n_position-1))
	allocate (d_hxx(0:n_position-1))	
	allocate (d_hxy(0:n_position-1))
	allocate (disorder_index(1:n_position))
 	
 	call MKL_VSL_UNIFORM (n_position, p_x(0), one, real(nx*1.0-2.0,kind=8))
 	ip_x(:)=int(p_x(:))
 	call MKL_VSL_UNIFORM (n_position, p_y(0), one, real(ny*1.0-2.0,kind=8))
 	ip_y(:)=int(p_y(:))
 	call MKL_VSL_GAUSSIAN (n_position, diameter_normal(0), size_mean_normal, size_sigma_normal)
 	diameter(:) = int(exp(diameter_normal(:)))
 	call MKL_VSL_GAUSSIAN (n_position, d_hxx(0), hxx_mean, hxx_sigma)
 	call MKL_VSL_GAUSSIAN (n_position, d_hxy(0), hxy_mean, hxy_sigma)
 	write (*,*) "d_hxx",d_hxx(:)
 	
 	do i=1, n_position
 		! (i) Check if the random disorder is over the lattice boundaries
 		boundary = ((ip_x(i)+diameter(i))>(nx-1)) .or. ((ip_x(i)-diameter(i))<0) .or. &
 		((ip_y(i)+diameter(i))>(ny-1)) .or. ((ip_y(i)-diameter(i))<0)
 		if (boundary) then
 			write (*,*) "outside"
 			cycle
 		endif
 		! (ii) Check if there are disorder overlapping with each others		
 		do j=1, disorder_number
 			ii = disorder_index (j)
 			temp = ((ip_x(i)-ip_x(ii))*(ip_x(i)-ip_x(ii)) + (ip_y(i)-ip_y(ii))*(ip_y(i)-ip_y(ii)))**0.5
 			if (temp < (diameter(i)+diameter(ii))) then
 				write (*,*) "overlape"
				flag= .false. ; exit
 			endif	
 		enddo
 		! (iii) input disorder into the lattice
 		if (flag) then
 			i_range = diameter(i)		
 			do m=-1*i_range, i_range
 				j_range=floor((i_range*i_range-m*m)**0.5)
 				do k=-1*j_range, j_range
 					ii=ip_x(i)+m; jj=ip_y(i)+k
 					hxx(ii,jj,0) = d_hxx(i)
 					hxy(ii,jj,0) = d_hxy(i)
 					fraction_real = fraction_real+1.0
 					write (*,*) fraction_real
 				enddo
 			enddo
 			disorder_number = disorder_number+1
 			write(*,*) "disorder_number", disorder_number
 			disorder_index(disorder_number) = i
 		endif	
 		flag= .true. 
 		! (iv) calculate the volume fraction
 	enddo
 	fraction_real = fraction_real/(1.0*nx)/(1.0*ny)
 	deallocate (p_x, p_y, ip_x, ip_y, diameter, diameter_normal, d_hxx, d_hxy, disorder_index)	
END SUBROUTINE 



SUBROUTINE  threshold_generator (dist_gauss, fraction, dist_mean, dist_sigma)
	USE parameters 
	IMPLICIT NONE

	real(kind=8),  intent(in) :: fraction, dist_mean, dist_sigma
	real(kind=8),  intent(inout),  dimension(0:nx-1, 0:ny-1, 0:nz-1) :: dist_gauss
	
	integer(kind=4) :: i,i1,j1, n_position
	real(kind=8) :: zero, one, dist_mean_normal, dist_sigma_normal
	integer(kind=8), dimension(:), allocatable :: ip_x, ip_y
      real(kind=8), dimension(:), allocatable :: p_x, p_y, dist_amplitude, dist_amplitude_normal
    
	zero=0.0; one=1.0; dist_gauss(:,:,:)=0.0
	
	dist_mean_normal = dlog(dist_mean)-0.5*dlog(1.0+dist_sigma/dist_mean/dist_mean)
	dist_sigma_normal = (dlog(1.0+dist_sigma/dist_mean/dist_mean))**0.5
	n_position = floor (fraction*nx*ny)
 	
 	allocate (p_x(0:n_position-1))		
 	allocate (p_y(0:n_position-1))
	allocate (ip_x(0:n_position-1)) 			
	allocate (ip_y(0:n_position-1))
	allocate (dist_amplitude(0:n_position-1))
 	allocate (dist_amplitude_normal(0:n_position-1))
 	
 	call MKL_VSL_UNIFORM (n_position, p_x(0), zero, real(nx*1.0-1.0,kind=8))
 	ip_x(:)=int(p_x(:))
 	call MKL_VSL_UNIFORM (n_position, p_y(0), zero, real(ny*1.0-1.0,kind=8))
 	ip_y(:)=int(p_y(:))
 	call MKL_VSL_GAUSSIAN (n_position, dist_amplitude_normal(0), dist_mean_normal, dist_sigma_normal)
 	dist_amplitude(:) = exp(dist_amplitude_normal(:))
 	
 	write (*,*) "dist_amplitude",dist_amplitude(:)
 	
 	do i=0, n_position-1
 		i1=ip_x(i); j1=ip_y(i)
 		dist_gauss(i1,j1,0) = dist_amplitude(i) 
 	enddo

END SUBROUTINE 












