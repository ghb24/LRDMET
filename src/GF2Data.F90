module GF2Data
    use const
    implicit none

    integer :: nMatsubara
    real(dp), allocatable :: MatsuPoints(:)

    !nImTimePoints = ScaleImTime * nMatsubara * Beta_Temp / pi
    integer :: nImTimePoints    
    real(dp), parameter :: ScaleImTime = 20.0_dp
    !The first point will be zero, and the last Beta_Temp
    real(dp), allocatable :: ImTimePoints(:)

    real(dp) :: Beta_Temp

    complex(dp), allocatable :: SE_Matsu_GV(:,:,:)
    complex(dp), allocatable :: SE_Tau_GV(:,:,:)
    complex(dp), allocatable :: GLat_Matsu_GV(:,:,:)
    complex(dp), allocatable :: GLat_Tau_GV(:,:,:)

    real(dp), allocatable :: FockMat_GV(:,:)
    real(dp), allocatable :: DensityMat_GV(:,:) !But should this be complex? 
    real(dp), allocatable :: DensityMat_MF_GV(:,:)
    real(dp) :: nElec_GV


end module GF2Data
