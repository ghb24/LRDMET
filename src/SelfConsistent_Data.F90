module SC_Data
    implicit none


    real(dp), allocatable :: FreqPoints(:),Weights(:)
    logical :: tCalcRealSpectrum=.true.  !Whether to calculate the real spectrum in addition to the fit spectrum
    logical :: tFitMatAxis               !Whether to fit to the matsubara axis 


end module SC_Data
