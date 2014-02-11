module SC_Data
    use const
    implicit none

    real(dp), allocatable :: FreqPoints(:),Weights(:)
    logical :: tCalcRealSpectrum=.true.  !Whether to calculate the real spectrum in addition to the fit spectrum
    logical :: tFitMatAxis               !Whether to fit to the matsubara axis 
    real(dp) :: SCF_mu
    integer :: SCF_iCorrFn

end module SC_Data
