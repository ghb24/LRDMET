module GF2Data
    use const
    implicit none

    real(dp) :: Beta_Temp

    integer :: nMatsubara
    real(dp), allocatable :: MatsuPoints(:)

    !For fitting of tails
    logical :: tFitTails
    real(dp) :: TailStart
    integer :: iTailNeg !End point of negative frequency tail
    integer :: iTailPos !End point of positive frequency tail
    real(dp) :: c2_denom, c3_denom  !Denominators for fitting
    real(dp) :: MatsuEnergySumFac   !How much further to sum just tails for the energy

    !G0 quantities
    real(dp) :: GF0_ChemPot
    real(dp), allocatable :: C2_Coeffs_GF0(:,:)
    real(dp), allocatable :: C3_Coeffs_GF0(:,:)

    !nImTimePoints = ScaleImTime * nMatsubara * Beta_Temp / pi
    integer :: nImTimePoints    
    real(dp) :: ScaleImTime,ScaleImTimeSpline
    logical :: tSpline
    !The first point will be zero, and the last Beta_Temp
    real(dp), allocatable :: ImTimePoints(:)

    type GreensFunc
        complex(dp), allocatable :: Matsu(:,:,:)
        complex(dp), allocatable :: Tau(:,:,:)
        real(dp), allocatable :: C0_Coeffs(:,:)
        real(dp), allocatable :: C1_Coeffs(:,:)
        real(dp), allocatable :: C2_Coeffs(:,:)
        real(dp), allocatable :: C3_Coeffs(:,:)
        logical :: tGF
    end type

    type(GreensFunc) GLat_GV
    type(GreensFunc) SE_GV

    real(dp), allocatable :: FockMat_GV(:,:)
    real(dp), allocatable :: DensityMat_GV(:,:) !But should this be complex? 
    real(dp), allocatable :: DensityMat_MF_GV(:,:)
    real(dp) :: nElec_GV


end module GF2Data
