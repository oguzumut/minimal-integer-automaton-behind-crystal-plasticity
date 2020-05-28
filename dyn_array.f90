
! MODULE Dyn_Array
! USE mkl_pardiso
! !---------------------------------------------------------------------
! !
! !  Module containing definitions needed to dynamically allocate 
! !  the values of  arrays 
! !
! !---------------------------------------------------------------------
! TYPE(MKL_PARDISO_HANDLE), save,ALLOCATABLE  :: pt(:)
! 
! END MODULE Dyn_Array
! 
! 
! MODULE Dyn_Array2
! !---------------------------------------------------------------------
! !
! !  Module containing definitions needed to dynamically allocate 
! !  the values of  arrays 
! !
! !---------------------------------------------------------------------
! INTEGER(kind=4), save,  allocatable ::   irow(:),icol(:)
! REAL(kind=8), save,  allocatable ::   nza(:)
! 
! END MODULE Dyn_Array2
! 
! MODULE Dyn_Array3
! !---------------------------------------------------------------------
! !
! !  Module containing definitions needed to dynamically allocate 
! !  the values of  arrays 
! !
! !---------------------------------------------------------------------
! INTEGER(kind=8), save,ALLOCATABLE  :: ptn(:)
! 
! END MODULE Dyn_Array3
! 
! 
! MODULE Dyn_Array4
! !---------------------------------------------------------------------
! !
! !  Module containing definitions needed to dynamically allocate 
! !  the values of  arrays 
! !
! !---------------------------------------------------------------------
! INTEGER(kind=8), save,ALLOCATABLE  :: ptn_fcc(:)
! 
! END MODULE Dyn_Array4



! MODULE phase_field_arrays
! !---------------------------------------------------------------------
! !
! !  Module containing definitions needed to dynamically allocate 
! !  the values of  arrays 
! !
! !---------------------------------------------------------------------
! 
!     real(kind=8), dimension(:,:,:),  save, allocatable   ::   theta
!     real(kind=8), dimension(:,:,:),  save, allocatable   ::   eta
! 
! 
!     complex(kind=8), dimension(:,:,:),  save, allocatable   ::   tf_theta
!     complex(kind=8), dimension(:,:,:),  save, allocatable   ::   tf_eta
! 
! 
! 
! END MODULE phase_field_arrays


MODULE potential_arrays
!---------------------------------------------------------------------
!
!  Module containing definitions needed to dynamically allocate 
!  the values of  arrays 
!
!---------------------------------------------------------------------

    real(kind=8), dimension(:,:,:),  save, allocatable   ::   exx,exy,exy_homo, exy_fluc, hxx, hxy, right, u_field_0
    real(kind=8), dimension(:,:,:),  save, allocatable   ::   d_field,u_field, work_space, k_field, g_field

    complex(kind=8), dimension(:,:,:),  save, allocatable   ::   tf_exx,tf_exy, tf_hxx, tf_hxy, tf_right
    complex(kind=8), dimension(:,:,:),  save, allocatable   ::   tf_d_field,tf_u_field, tf_u_field_0



END MODULE potential_arrays


MODULE parameters
!---------------------------------------------------------------------
!
!  Module containing parameters
!
!---------------------------------------------------------------------
  real(kind=8)    :: d_time,pi,pi2
 
  integer(kind=4) :: nx,ny,nz
  

END MODULE parameters

MODULE fourier_operators_discrete
!---------------------------------------------------------------------
!
!  Module containing definitions needed to dynamically allocate 
!  the values of  arrays 
!
!---------------------------------------------------------------------

    real(kind=8), dimension(:,:,:),  save, allocatable   ::   qxx,qyy,qzz,sqr_norme_q,tol
    complex(kind=8), dimension(:,:,:),  save, allocatable   ::   qxp,qyp
    complex(kind=8), dimension(:,:,:),  save, allocatable   ::   qxn,qyn




END MODULE fourier_operators_discrete

MODULE threshholds
!---------------------------------------------------------------------
!
!  Module containing definitions needed to dynamically allocate 
!  the values of  arrays 
!
!---------------------------------------------------------------------

    real(kind=8), dimension(:,:,:),  save, allocatable   ::   dist_gauss,dist_uniform, threshold




END MODULE threshholds



! Module potential_tree
! 
!   Type :: tree_node
!       Integer :: dnum
!       Real :: val
!       Integer :: l, u
!   End Type tree_node
! 
! 
! end Module potential_tree









