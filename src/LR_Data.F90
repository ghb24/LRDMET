module LR_Data
    use const
    implicit none
    save
    complex(dp), allocatable, target :: LinearSystem_p(:,:)
    complex(dp), allocatable, target :: LinearSystem_h(:,:)
    complex(dp), pointer :: zDirMV_Mat(:,:)
    complex(dp), pointer :: zDirMV_Mat_cmprs(:)
    integer, pointer :: zDirMV_Mat_cmprs_inds(:)
    real(dp), allocatable :: Precond_Diag(:)
!$OMP THREADPRIVATE(LinearSystem_h,LinearSystem_p)
!$OMP THREADPRIVATE(zDirMV_Mat)
!$OMP THREADPRIVATE(zDirMV_Mat_cmprs)
!$OMP THREADPRIVATE(zDirMV_Mat_cmprs_inds)
!$OMP THREADPRIVATE(Precond_Diag)

end module LR_Data
