module SC_Data
    use const
    implicit none

    real(dp), allocatable :: FreqPoints(:),Weights(:)
    logical :: tFitMatAxis               !Whether to fit to the matsubara axis 
    integer :: nIndKPnts                !Number of independent kpoints in the fitting
    real(dp) :: SCF_mu
    integer :: SCF_iCorrFn

end module SC_Data
