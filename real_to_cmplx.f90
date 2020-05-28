

	subroutine r2c(lx,ly,lz,fft_name,work_real_in,work_cmplx_out,n_dim_calcul)

	use OMP_LIB
	use MKL_DFTI

!----------------------------- dimensionnements generaux ------------------------------------

    real(kind=8)       work_real_in(lx*ly*lz)
     
    complex(kind=8)    work_cmplx_out((lx/2+1)*ly*lz)    

    character*20 fft_name
        
    Integer :: status
             
    type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
    Integer(kind=4),save,       dimension(:),  allocatable ::  strides_out,strides_in,L

	real(kind=8) :: scale

    integer(kind=4) ,   intent(in)  :: n_dim_calcul

    integer(kind=4) :: compteur = 0


!********************** FFT avec "intel" *******************************************************
!

if(compteur .eq. 0)then
allocate(strides_out(n_dim_calcul+1),strides_in(n_dim_calcul+1),L(n_dim_calcul))

!--------------------------------- allocations pour "intel" (galalibier calculateur) ----------
        
        
        
      !allocate( work_cmplx(((lx/2)+1)*(ly*lz)) )


      if(n_dim_calcul .eq. 3)then
      strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = lx
      strides_in(4) = lx*ly

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = lx/2+1
	  strides_out(4) = (lx/2+1) * ly

		L(1) = lx
		L(2) = ly
		L(3) = lz
	  else if(n_dim_calcul .eq. 2)then
	  write(*,*)'here'
	  strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = lx

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = lx/2+1

		L(1) = lx
		L(2) = ly
	  else if(n_dim_calcul .eq. 1)then
	  strides_in(1) = 0
      strides_in(2) = 1

      strides_out(1) = 0
      strides_out(2) = 1

		L(1) = lx
	  end if	  
	  
	  compteur = 1
end if	  

!------------------------------------------------------------------------------

		

    Scale = 1.0/real(lx*ly, KIND=8)
	
	Status = DftiCreateDescriptor(My_Desc1_Handle,DFTI_DOUBLE,DFTI_REAL,n_dim_calcul,L)
	Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, Scale)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_INPUT_STRIDES, strides_in)
 	Status = DftiSetValue(My_Desc1_Handle, DFTI_OUTPUT_STRIDES, strides_out)
	Status = DftiCommitDescriptor(My_Desc1_Handle)	
	!tps1 = second()
	Status = DftiComputeForward(My_Desc1_Handle,work_real_in,work_cmplx_out)
	!tps2 = second()
	!write(*,*)'intel fft',tps2-tps1


	!call reshape_forward(lx,ly,lz,work_cmplx,work_cmplx_out)

	Status = DftiFreeDescriptor(My_Desc1_Handle)

!------------------------------------------------------------------------------

   ! deallocate(strides_in,strides_out,L)

!
!*******************************************************************************************

    return

  end subroutine r2c
  
subroutine c2r(lx,ly,lz,fft_name,work_cmplx_in,work_real_out,n_dim_calcul)

use OMP_LIB
use MKL_DFTI

!----------------------------- dimensionnements generaux ------------------------------------

    complex(kind=8)    work_cmplx_in((lx/2+1)*ly*lz)
    
    real(kind=8)    work_real_out(lx*ly*lz)    

    character*20 fft_name
        
    Integer :: status
             
    type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle

    Integer(kind=4),save,       dimension(:),  allocatable ::  strides_out,strides_in,L

   integer(kind=4) ,   intent(in)  :: n_dim_calcul


	integer(kind=4) :: compteur = 0


!********************** FFT avec "intel" *******************************************************
!

if(compteur .eq. 0)then
allocate(strides_out(n_dim_calcul+1),strides_in(n_dim_calcul+1),L(n_dim_calcul))

!--------------------------------- allocations pour "intel" (galalibier calculateur) ----------
        
        
        
      !allocate( work_cmplx(((lx/2)+1)*(ly*lz)) )


      if(n_dim_calcul .eq. 3)then
     
      strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = lx
      strides_in(4) = lx*ly

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = lx/2+1
	  strides_out(4) = (lx/2+1) * ly

		L(1) = lx
		L(2) = ly
		L(3) = lz
	  else if(n_dim_calcul .eq. 2)then
	  strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = lx

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = lx/2+1

		L(1) = lx
		L(2) = ly
	  else if(n_dim_calcul .eq. 1)then
	  strides_in(1) = 0
      strides_in(2) = 1

      strides_out(1) = 0
      strides_out(2) = 1

		L(1) = lx
	  end if	  
	  
	  compteur = 1
end if	  



!------------------------------------------------------------------------------

		

    !call reshape_backward(lx,ly,lz,work_cmplx_in,work_cmplx)
     ! Scale = 1.0/real(lx*ly*lz, KIND=8)

     ! Status = DftiSetValue(My_Desc1_Handle, DFTI_BACKWARD_SCALE, Scale)
   
    Status = DftiCreateDescriptor(My_Desc1_Handle,DFTI_DOUBLE,DFTI_REAL,n_dim_calcul,L)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_INPUT_STRIDES, strides_out)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_OUTPUT_STRIDES, strides_in)  
    Status = DftiCommitDescriptor(My_Desc1_Handle)

	!tps1 = second()
	Status = DftiComputeBackward(My_Desc1_Handle,work_cmplx_in,work_real_out)
	!tps2 = second()
	!write(*,*)'intel fft',tps2-tps1
	!work(:) = work_real(:)



	Status = DftiFreeDescriptor(My_Desc1_Handle)
  !  deallocate(strides_in,strides_out,L)


!------------------------------------------------------------------------------
!*******************************************************************************************

    return

  end subroutine c2r  
  
  !*****************************************************************
	subroutine r2c_spec(lx,ly,lz,fft_name,work_real_in,work_cmplx_out,n_dim_calcul)

	use OMP_LIB
	use MKL_DFTI

!----------------------------- dimensionnements generaux ------------------------------------

    real(kind=8)       work_real_in(lx*ly*lz)
     
    complex(kind=8)    work_cmplx_out((lx/2+1)*ly*lz)    

    character*20 fft_name
        
    Integer :: status
             
    type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
    Integer(kind=4),save,       dimension(:),  allocatable ::  strides_out,strides_in,L

	real(kind=8) :: scale


    integer(kind=4) :: compteur = 0


!********************** FFT avec "intel" *******************************************************
!

if(compteur .eq. 0)then
allocate(strides_out(n_dim_calcul+1),strides_in(n_dim_calcul+1),L(n_dim_calcul))

!--------------------------------- allocations pour "intel" (galalibier calculateur) ----------
        
        
        
      !allocate( work_cmplx(((lx/2)+1)*(ly*lz)) )


      if(n_dim_calcul .eq. 3)then
      strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = lx
      strides_in(4) = lx*ly

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = lx/2+1
	  strides_out(4) = (lx/2+1) * ly

		L(1) = lx
		L(2) = ly
		L(3) = lz
	  else if(n_dim_calcul .eq. 2)then
	  strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = lx

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = lx/2+1

		L(1) = lx
		L(2) = ly
	  else if(n_dim_calcul .eq. 1)then
	  strides_in(1) = 0
      strides_in(2) = 1

      strides_out(1) = 0
      strides_out(2) = 1

		L(1) = lx
	  end if	  
	  
	  compteur = 1
end if	  

!------------------------------------------------------------------------------

		

    Scale = 1.0/real((lx)*ly*lz)

	Status = DftiCreateDescriptor(My_Desc1_Handle,DFTI_DOUBLE,DFTI_REAL,n_dim_calcul,L)
	Status = DftiSetValue(My_Desc1_Handle, DFTI_FORWARD_SCALE, Scale)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
    Status = DftiSetValue(My_Desc1_Handle, DFTI_INPUT_STRIDES, strides_in)
 	Status = DftiSetValue(My_Desc1_Handle, DFTI_OUTPUT_STRIDES, strides_out)
	Status = DftiCommitDescriptor(My_Desc1_Handle)	
	!tps1 = second()
	Status = DftiComputeForward(My_Desc1_Handle,work_real_in,work_cmplx_out)
	!tps2 = second()
	!write(*,*)'intel fft',tps2-tps1


	!call reshape_forward(lx,ly,lz,work_cmplx,work_cmplx_out)

	Status = DftiFreeDescriptor(My_Desc1_Handle)

!------------------------------------------------------------------------------

   ! deallocate(strides_in,strides_out,L)

!
!*******************************************************************************************

    return

  end subroutine r2c_spec