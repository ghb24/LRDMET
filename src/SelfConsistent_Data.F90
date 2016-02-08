module SC_Data
    use const
    implicit none

    real(dp), allocatable :: FreqPoints(:),Weights(:)
    logical :: tFitMatAxis               !Whether to fit to the matsubara axis 
    integer :: nIndKPnts                !Number of independent kpoints in the fitting
    real(dp) :: SCF_mu
    integer :: SCF_iCorrFn
    logical :: tReadChemPot !Do we read the chemical potential from the file too? 

    !Stretch initial hamiltonian?
    logical :: tStretchNILatticeHam
    real(dp) :: dStretchNILatticeHam

    real(dp) :: dSelfConsConv

    !Damp the QPSC Updates
    real(dp) :: PotentialUpdateDamping
    
    integer, allocatable :: KBlock_to_KInd(:)       !Map from Kpoint index to sampled kpoint index
    integer, allocatable :: KInd_to_KBlock(:)       !Map from sampled kpoint to physical kpoint index
    logical, allocatable :: KPointSampled(:)        !A logical to indicate whether the kpoint is explicitly sampled, or just reconstructed from the others
    !Same as above, but for the imposed kpoint symmetry
    integer, allocatable :: KBlock_to_KInd_Constrain(:)
    integer, allocatable :: KInd_to_KBlock_Constrain(:)
    logical, allocatable :: KPointSampled_Constrain(:)

end module SC_Data
