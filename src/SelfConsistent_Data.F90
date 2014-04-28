module SC_Data
    use const
    implicit none

    real(dp), allocatable :: FreqPoints(:),Weights(:)
    logical :: tFitMatAxis               !Whether to fit to the matsubara axis 
    real(dp) :: SCF_mu
    integer :: SCF_iCorrFn

end module SC_Data
