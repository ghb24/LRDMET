module SelfConsistentLR
    use const
    use errors
    use LinearResponse
    use MomSpectra
    use globals
    use Continuation
    use utils, only: get_free_unit,append_ext_real,append_ext
    use mat_tools, only: AddPeriodicImpCoupling_RealSpace 
    implicit none

    contains
    
    !Self-consistently fit a frequency *independent* lattice coupling to the impurity in order to
    !match the matsubara greens functions
    subroutine SC_FitLatticeGF_Im()
        implicit none
        real(dp) :: GFChemPot,Omega
        integer :: nESteps_Im,nESteps_Re,OmegaVal,nFitPoints
        complex(dp), allocatable :: SE_Fit(:,:,:),G_Mat_Fit(:,:,:)
        complex(dp), allocatable :: SE_Re(:,:,:),G_Mat_Re(:,:,:)
        complex(dp), allocatable :: Lattice_GF(:,:,:),Lattice_GF_Real(:,:,:),Lattice_GF_Old(:,:,:)
        real(dp), allocatable :: h_lat_fit(:,:),Couplings(:,:),FreqPoints_temp(:),Weights_temp(:)
        real(dp), allocatable :: AllDiffs(:,:),FreqPoints(:),IntWeights(:)
        complex(dp), allocatable :: DiffImpGF(:,:,:),G_Mat_Fit_Old(:,:,:)
        real(dp) :: FinalDist,LowFreq,HighFreq
        integer :: i,iter,iLatParams,j
        logical :: tFitMatAxis,tCalcRealSpectrum
        real(dp), parameter :: dDeltaImpThresh = 1.0e-4_dp 
        logical, parameter :: tUsePoswFreqPoints = .true. !For ph symmetric systems, this seems to make no difference
        character(len=*), parameter :: t_r='SC_FitLatticeGF_Im'

        if(tDiag_kSpace) then
            do i = 1,nSites
                write(6,*) i,abs(RtoK_Rot(1,i))**2
            enddo
        else
            if(tKPntSymFit) call stop_all(t_r,'Should not be trying to conserve k-point symmetry without any!') 
            do i = 1,nSites
                write(6,*) i,HFOrbs(1,i),HFOrbs(1,i)**2 
            enddo
        endif

        !Set chemical potential (This will only be right for half-filling)
        if(.not.tAnderson) then
            !Hubbard chemical potential
            GFChemPot = U/2.0_dp
        else
            !For anderson model, we want h0v to be h0 with -U/2 on impurity
            !Then, we want no chemical potential for the linear response.
            GFChemPot = zero
        endif

        !How many frequency points are there exactly?
        nESteps_Im = 0
        nESteps_Re = 0
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.true.)
            if(OmegaVal.lt.0) exit
            nESteps_Im = nESteps_Im + 1
        enddo
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.false.)
            if(OmegaVal.lt.0) exit
            nESteps_Re = nESteps_Re + 1
        enddo
        if(tFitPoints_Legendre) then
            !Calculate the Legendre polynomials for Legendre integration points and weights
            allocate(FreqPoints(nFreqPoints))
            allocate(IntWeights(nFreqPoints))
            write(6,"(A)") "Calculating frequency points..."
            call flush(6)
            LowFreq = -1.0_dp
            HighFreq = 1.0_dp
            call gauleg(LowFreq,HighFreq,FreqPoints,IntWeights,nFreqPoints)
            call ScaleFreqPointsWeights(nFreqPoints,FreqPoints,IntWeights)
            write(6,"(A)") "Frequency integration points: "
            do i = 1,nFreqPoints
                write(6,"(I9,2G20.13)") i,FreqPoints(i),IntWeights(i)
            enddo
            if(.not.tOptGF_EVals) call stop_all(t_r,'Must use evals for legendre quadrature')
            if(tUsePoswFreqPoints) then
                !Count the positive values of w
                write(6,"(A)") "Just using the positive values of omega"
                j = 0
                do i = 1,nFreqPoints
                    if(FreqPoints(i).ge.zero) j = j+1
                enddo
                allocate(FreqPoints_temp(j))
                allocate(Weights_temp(j))
                j = 0
                do i = 1,nFreqPoints
                    if(FreqPoints(i).ge.zero) then
                        j = j+1
                        FreqPoints_temp(j) = FreqPoints(i)
                        Weights_temp(j) = IntWeights(i)
                    endif
                enddo
                deallocate(FreqPoints,IntWeights)
                nFreqPoints = j
                allocate(FreqPoints(nFreqPoints))
                allocate(IntWeights(nFreqPoints))
                FreqPoints(:) = FreqPoints_temp(:)
                IntWeights(:) = Weights_temp(:)
                deallocate(FreqPoints_temp,Weights_temp)
                write(6,"(A,I7)") "New number of frequency points: ",nFreqPoints
            endif
        endif
        call flush(6)

        !Which axis do we want to do the fitting on?
        if(tFitRealFreq) then
            tFitMatAxis = .false.
            nFitPoints = nESteps_Re
            tCalcRealSpectrum = .false.
        else
            tFitMatAxis = .true.
            tCalcRealSpectrum = .true.
            if(tFitPoints_Legendre) then
                nFitPoints = nFreqPoints
            else
                nFitPoints = nESteps_Im
            endif
        endif

        !Initially, just see if we can fit the two different Matsubara spectral functions
        write(6,*) "Number of frequency points to fit across: ",nFitPoints
        if(tOptGF_EVals) then
            if(tKPntSymFit) then
                !We don't actually want to optimize all of these.
                !e(k) = e(-k), and we are always using a uniform mesh
                !If gamma-centered mesh, then have nSites/2 + 1 independent parameters (we are sampling k=0 and BZ boundary which dont pair)
                !If Shifted mesh, then we have nSites/2 independent parameters
                if(tShift_Mesh) then
                    iLatParams = nSites/2
                else
                    iLatParams = (nSites/2) + 1
                endif
            else
                iLatParams = nSites !Optimize all kpoints independently
            endif
            if(iLatticeCoups.ne.0) then
                call stop_all(t_r,'iLatticeCoups set, but expecting to optimize eigenvalues')
            endif
        else
            !Should we do including the inpurity space in these lattice couplings?
            iLatParams = iLatticeCoups
        endif

        allocate(SE_Fit(nImp,nImp,nFitPoints))
        SE_Fit(:,:,:) = zzero
        allocate(h_lat_fit(nSites,nSites))
        allocate(Couplings(iLatParams,nImp))

        !This will initialize the lattice eigenvalues/couplings
        call InitLatticeCouplings(Couplings,iLatParams,h_lat_fit)
        !This will create the lattice hamiltonian, either from the eigenvalues or the couplings
        !write(6,*) h_lat_fit(1,:)
        !call writematrix(h_lat_fit,'h_lat_fit_init',.true.)
        call AddPeriodicImpCoupling_RealSpace(h_lat_fit,nSites,iLatParams,nImp,Couplings)
        
        write(6,"(A)") "Starting *frequency independent* lattice parameters: "
        do i = 1,iLatParams
            write(6,"(I6,G20.13)") i,Couplings(i,1)
        enddo
        write(6,"(A)") "" 

        allocate(G_Mat_Fit(nImp,nImp,nFitPoints))
        allocate(G_Mat_Re(nImp,nImp,nESteps_Re))
        allocate(SE_Re(nImp,nImp,nESteps_Re))
        SE_Re(:,:,:) = zzero
        allocate(Lattice_GF(nImp,nImp,nFitPoints))
        allocate(Lattice_GF_Old(nImp,nImp,nFitPoints))
        allocate(Lattice_GF_Real(nImp,nImp,nESteps_Re))
        !Write out the initial lattice greens function
        !call FindLocalMomGF(nESteps_Im,SE_Im,Lattice_GF,tMatbrAxis=.true.,ham=h_lat_fit)
        if(tFitPoints_Legendre) then
            call FindLocalMomGF(nFitPoints,SE_Fit,Lattice_GF,tMatbrAxis=tFitMatAxis,ham=h_lat_fit,    &
                CouplingLength=iLatParams,evals=Couplings(:,1),FreqPoints=FreqPoints)
        else
            call FindLocalMomGF(nFitPoints,SE_Fit,Lattice_GF,tMatbrAxis=tFitMatAxis,ham=h_lat_fit,    &
                CouplingLength=iLatParams,evals=Couplings(:,1))
        endif

        allocate(G_Mat_Fit_Old(nImp,nImp,nFitPoints))
        allocate(DiffImpGF(nImp,nImp,nFitPoints))
        G_Mat_Fit(:,:,:) = zzero

        allocate(AllDiffs(3,0:iMaxIter_MacroFit))
        AllDiffs(:,:) = zero

        iter = 0
        do while(.not.tSkip_Lattice_Fit)
            iter = iter + 1
            
            !High level calculation on matsubara axis, with no self-energy, and just the normal lattice hamiltonian
            G_Mat_Fit_Old(:,:,:) = G_Mat_Fit(:,:,:)
            if(tFitPoints_Legendre) then
                call SchmidtGF_wSE(G_Mat_Fit,GFChemPot,SE_Fit,nFitPoints,tMatbrAxis=tFitMatAxis,ham=h_lat_fit,FreqPoints=FreqPoints)
                call writedynamicfunction(nFitPoints,G_Mat_Fit,'G_Imp_Fit',tag=iter,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
            else
                call SchmidtGF_wSE(G_Mat_Fit,GFChemPot,SE_Fit,nFitPoints,tMatbrAxis=tFitMatAxis,ham=h_lat_fit)
                call writedynamicfunction(nFitPoints,G_Mat_Fit,'G_Imp_Fit',tag=iter,tMatbrAxis=tFitMatAxis)
            endif
                
            if(tCalcRealSpectrum) then
                call SchmidtGF_wSE(G_Mat_Re,GFChemPot,SE_Re,nESteps_Re,tMatbrAxis=.false.,ham=h_lat_fit)
                call writedynamicfunction(nESteps_Re,G_Mat_Re,'G_Imp_Re',tag=iter,tMatbrAxis=.false.)
                call FindLocalMomGF(nESteps_Re,SE_Re,Lattice_GF_Real,tMatbrAxis=.false.,ham=h_lat_fit,  &
                    CouplingLength=iLatParams,evals=Couplings(:,1))
                call writedynamicfunction(nESteps_Re,Lattice_GF_Real,'G_Lat_Re',tag=iter,tMatbrAxis=.false.)
            endif
            
            !Write out lattice greens function
            if(tFitPoints_Legendre) then
                call writedynamicfunction(nFitPoints,Lattice_GF,'G_Lat_Fit',tag=iter,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
            else
                call writedynamicfunction(nFitPoints,Lattice_GF,'G_Lat_Fit',tag=iter,tMatbrAxis=tFitMatAxis)
            endif

            if(iLatticeFitType.eq.2) then
                !Invert the high-level greens function, since we are fitting the residual of the inverses
                call InvertLocalNonHermGF(nFitPoints,G_Mat_Fit)
            endif

            if((abs(Damping_SE-one).gt.1.0e-8_dp).and.(iter.gt.1)) then
                !Take an admixture of the previous two high-level calculations to damp the fit
                G_Mat_Fit(:,:,:) = (Damping_SE*G_Mat_Fit(:,:,:)) + ((one-Damping_SE)*G_Mat_Fit_Old(:,:,:))
            endif

            if(iter.eq.1) then
                !calculate the initial residual
                if(tFitPoints_Legendre) then
                    call CalcLatticeFitResidual(G_Mat_Fit,nFitPoints,Couplings,iLatParams,AllDiffs(1,0),tFitMatAxis,   & 
                        FreqPoints=FreqPoints,Weights=IntWeights)
                else
                    call CalcLatticeFitResidual(G_Mat_Fit,nFitPoints,Couplings,iLatParams,AllDiffs(1,0),tFitMatAxis)
                endif
                write(6,"(A,F20.10)") "Initial spectrum residual: ",AllDiffs(1,0)
                AllDiffs(2,0) = zero
            endif
            if(tFitPoints_Legendre) then
                call FitLatticeCouplings(G_Mat_Fit,nFitPoints,Couplings,iLatParams,FinalDist,tFitMatAxis,   &
                    FreqPoints=FreqPoints,Weights=IntWeights)
            else
                call FitLatticeCouplings(G_Mat_Fit,nFitPoints,Couplings,iLatParams,FinalDist,tFitMatAxis)
            endif
            !Add the new lattice couplings to the lattice hamiltonian
            call AddPeriodicImpCoupling_RealSpace(h_lat_fit,nSites,iLatParams,nImp,Couplings)
            !Now, calculate the local lattice greens function
            Lattice_GF_Old(:,:,:) = Lattice_GF(:,:,:)
            if(tFitPoints_Legendre) then
                call FindLocalMomGF(nFitPoints,SE_Fit,Lattice_GF,tMatbrAxis=tFitMatAxis,ham=h_lat_fit,    &
                    CouplingLength=iLatParams,evals=Couplings(:,1),FreqPoints=FreqPoints)
            else
                call FindLocalMomGF(nFitPoints,SE_Fit,Lattice_GF,tMatbrAxis=tFitMatAxis,ham=h_lat_fit,    &
                    CouplingLength=iLatParams,evals=Couplings(:,1))
            endif

            !Write out the lattice couplings
            call WriteLatticeCouplings(iLatParams,Couplings)
        
            !Calculate & write out stats (all of them)
            AllDiffs(1,iter) = FinalDist
            DiffImpGF(:,:,:) = G_Mat_Fit_Old(:,:,:) - G_Mat_Fit(:,:,:)
            DiffImpGF(:,:,:) = DiffImpGF(:,:,:) * dconjg(DiffImpGF(:,:,:))
            AllDiffs(2,iter) = sum(real(DiffImpGF(:,:,:),dp))
            !And diff for LatticeGF
            DiffImpGF(:,:,:) = Lattice_GF_Old(:,:,:) - Lattice_GF(:,:,:)
            DiffImpGF(:,:,:) = DiffImpGF(:,:,:) * dconjg(DiffImpGF(:,:,:))
            AllDiffs(3,iter) = sum(real(DiffImpGF(:,:,:),dp))

            write(6,"(A)") ""
            write(6,"(A,I7,A)") "***   COMPLETED MACROITERATION ",iter," ***"
            write(6,"(A)") "     Iter.  FitResidual      Delta_GF_Imp(iw)    Delta_GF_Lat(iw)"
            do i = 0,iter
                write(6,"(I7,3G20.13)") i,AllDiffs(1,i),AllDiffs(2,i),AllDiffs(3,i)
            enddo
            write(6,"(A)") ""
            call flush(6)
            
            if(iter.ge.iMaxIter_MacroFit) then
                write(6,"(A,I9)") "Exiting. Max iters hit of: ",iMaxIter_MacroFit
                exit
            endif

            if(AllDiffs(2,iter).lt.dDeltaImpThresh) then
                write(6,"(A)") "Success! Convergence on imaginary axis successful"
                write(6,"(A,G20.13)") "Impurity greens function changing on imaginary axis by less than: ",dDeltaImpThresh
                exit
            endif
            if(tOptGF_EVals) then
                do i = 1,iLatParams
                    if(abs(Couplings(i,1)).gt.1000.0_dp) then
                        write(6,*) Couplings(:,1)
                        call stop_all(t_r,'Lattice eigenvalues diverging. Convergence failed.')
                    endif
                enddo
            endif
        enddo

        !We should write these out to a file so that we can restart from them at a later date (possibly with more variables).
        write(6,"(A)") "Final *frequency independent* lattice parameters: "
        do i = 1,iLatParams
            write(6,"(I6,G20.13)") i,Couplings(i,1)
        enddo
        write(6,"(A)") "" 
    
        if(tOptGF_EVals.and.tDiag_kSpace) then
            !Write out the converged one-electron dispersion / bandstructure
            call WriteBandstructure(Couplings,iLatParams)
        else
            !Should fix this so that we convert to eigenvalues from couplings and write out
        endif
        if(tFitPoints_Legendre) then
            call writedynamicfunction(nFitPoints,Lattice_GF,'G_Lat_Fit_Final',tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
        else
            call writedynamicfunction(nFitPoints,Lattice_GF,'G_Lat_Fit_Final',tMatbrAxis=tFitMatAxis)
        endif

        deallocate(SE_Fit,Lattice_GF,Lattice_GF_Real,G_Mat_Fit_Old,DiffImpGF,AllDiffs,Lattice_GF_Old)

        if(tFitMatAxis) then
            allocate(Lattice_GF(nImp,nImp,nESteps_Re))
                
            write(6,"(A)") "Now calculating spectral functions on the REAL axis"
            call flush(6)
            !What does the final lattice greens function look like on the real axis?
            call FindLocalMomGF(nESteps_Re,SE_Re,Lattice_GF,tMatbrAxis=.false.,ham=h_lat_fit,   &
                CouplingLength=iLatParams,evals=Couplings(:,1))
            !call FindLocalMomGF(nESteps_Re,SE_Re,Lattice_GF,tMatbrAxis=.false.,ham=h_lat_fit)
            !Write out lattice greens function
            !Is it the same as the impurity greens function on the real axis. This would be nice.
            call writedynamicfunction(nESteps_Re,Lattice_GF,'G_Lat_Re_Final',tMatbrAxis=.false.)

            !Finally, calculate the greens function in real-frequency space.
            call SchmidtGF_wSE(G_Mat_Re,GFChemPot,SE_Re,nESteps_Re,tMatbrAxis=.false.,ham=h_lat_fit)
            call writedynamicfunction(nESteps_Re,G_Mat_Re,'G_Imp_Re_Final',tMatbrAxis=.false.)
            deallocate(Lattice_GF)
        else
            if(tFitPoints_Legendre) then
                call writedynamicfunction(nESteps_Re,G_Mat_Fit,'G_Imp_Fit_Final',tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
            else
                call writedynamicfunction(nESteps_Re,G_Mat_Fit,'G_Imp_Fit_Final',tMatbrAxis=tFitMatAxis)
            endif
        endif

        deallocate(h_lat_fit,Couplings,SE_Re,G_Mat_Re,G_Mat_Fit)
        if(tFitPoints_Legendre) deallocate(FreqPoints,IntWeights)

    end subroutine SC_FitLatticeGF_Im

    !Self-consistency in the style of DMFT, converging both an imaginary frequency self-energy and a frequency independent hamiltonian
    !The aim is to find a self-energy that accurately represents the lattice self-energy.
    !At convergence, the *lattice* greens function should equal the *impurity* greens function
    !1) Calc G_0
    !2) Find 'bath' G - i.e. without self-energy on impurity
    !3) Fit lattice
    !4) Find G_imp
    !5) Dyson equation for updated self-energy
    subroutine SC_FitLat_and_SE_Im()
        implicit none
        real(dp) :: GFChemPot,FinalDist,Omega
        integer :: nESteps_Im,nESteps_Re,OmegaVal,iter,i,iLatParams
        complex(dp), allocatable :: SE_Im(:,:,:),G_Mat_Im(:,:,:),Lattice_GF(:,:,:)
        complex(dp), allocatable :: Bath_GF(:,:,:),zero_GF(:,:,:),Cluster_GF(:,:,:)
        complex(dp), allocatable :: DeltaSE(:,:,:),DiffImpGF(:,:,:),G_Mat_Re(:,:,:)
        complex(dp), allocatable :: G_Mat_Im_Old(:,:,:)
        real(dp), allocatable :: h_lat_fit(:,:),Couplings(:,:),AllDiffs(:,:)
        real(dp), parameter :: dDeltaThresh = 1.0e-6_dp
        logical, parameter :: tCreateBathGF = .true.
        character(len=*), parameter :: t_r='SC_FitLat_and_SE_Im'

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp

        !How many frequency points are there exactly?
        nESteps_Im = 0
        nESteps_Re = 0
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.true.)
            if(OmegaVal.lt.0) exit
            nESteps_Im = nESteps_Im + 1
        enddo
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.false.)
            if(OmegaVal.lt.0) exit
            nESteps_Re = nESteps_Re + 1
        enddo

        !Initially, just see if we can fit the two different Matsubara spectral functions
        write(6,*) "Number of Matsubara frequency points: ",nESteps_Im

        !The Greens functions
        allocate(G_Mat_Im(nImp,nImp,nESteps_Im))
        allocate(G_Mat_Im_Old(nImp,nImp,nESteps_Im))
        allocate(Lattice_GF(nImp,nImp,nESteps_Im))
        allocate(Bath_GF(nImp,nImp,nESteps_Im))
        allocate(Cluster_GF(nImp,nImp,nESteps_Im))

        !Self-energies
        allocate(SE_Im(nImp,nImp,nESteps_Im))
        allocate(DeltaSE(nImp,nImp,nESteps_Im))
        
        !Lattice fitting
        allocate(h_lat_fit(nSites,nSites))
        iLatParams = iLatticeCoups
        allocate(Couplings(iLatParams,nImp))

        !Misc/convergences
        allocate(zero_GF(nImp,nImp,nESteps_Im))
        allocate(DiffImpGF(nImp,nImp,nESteps_Im))
        allocate(AllDiffs(3,1:iMaxIter_MacroFit))

        !Initialize values
        !   fit lattice & couplings
        call InitLatticeCouplings(Couplings,iLatParams,h_lat_fit)
        call AddPeriodicImpCoupling_RealSpace(h_lat_fit,nSites,iLatParams,nImp,Couplings)

        !   Self-energy
        call Init_SE(nESteps_Im,SE_Im,tMatbrAxis=.true.)
        DeltaSE(:,:,:) = zzero
        
        !   Greens functions
        G_Mat_Im(:,:,:) = zzero
        G_Mat_Im_Old(:,:,:) = zzero
        Lattice_GF(:,:,:) = zzero
        Bath_GF(:,:,:) = zzero

        !   Convergence
        zero_GF(:,:,:) = zzero
        DiffImpGF(:,:,:) = zzero
        AllDiffs(:,:) = zero
            
        write(6,"(A)") "Starting self-consistent optimization of greens function..."

        iter = 0
        do while(.not.tSkip_Lattice_Fit)
            iter = iter + 1

            write(6,"(A,I7)") "Starting iteration ",iter

            !1) Lattice calculation with self-energy
            call FindLocalMomGF(nESteps_Im,SE_Im,Lattice_GF,tMatbrAxis=.true.)
            call writedynamicfunction(nESteps_Im,Lattice_GF,'G_Lat_Im',tag=iter,tMatbrAxis=.true.)

            write(6,"(A)") "Lattice greens function calculated and written to disk."

            !2) Optionally, create a 'bath' greens function, where the self-energy is removed from the impurity
            if(tCreateBathGF) then
                !The 'bath' greens function is defined as [G_0^-1 + SE]^-1
                Bath_GF(:,:,:) = Lattice_GF(:,:,:)
                call InvertLocalNonHermGF(nESteps_Im,Bath_GF)
                Bath_GF(:,:,:) = Bath_GF(:,:,:) + SE_Im(:,:,:)
                if(iLatticeFitType.ne.2) then
                    call InvertLocalNonHermGF(nESteps_Im,Bath_GF)
                    call writedynamicfunction(nESteps_Im,Bath_GF,'G_Bath_Im',tag=iter,tMatbrAxis=.true.)
                    write(6,"(A)") "Bath greens function calculated and written to disk."
                else
                    !Don't bother inverting, since we want the inverse of the bath greens function, 
                    !since we are fitting the residual of the inverses
                    write(6,"(A)") "Bath greens function calculated (but not written)."
                endif
            else
                !Alternatively, set the 'bath' greens function to be the 'lattice' one (with SE over impurity space).
                Bath_GF(:,:,:) = Lattice_GF(:,:,:)
                if(iLatticeFitType.eq.2) then
                    !Invert the bath greens function, since we are fitting the residual of the inverses
                    call InvertLocalNonHermGF(nESteps_Im,Bath_GF)
                endif
            endif
            
            !3) Fit a lattice hamiltonian, such that it reproduces the bath greens function
            !   This will fit couplings on to h0v, with no frequency-dependent potential, to match the 'bath' GF
            write(6,"(A)") "Fitting static lattice couplings to reproduce bath greens function..."
            call FitLatticeCouplings(Bath_GF,nESteps_Im,Couplings,iLatParams,FinalDist,.true.)
            !   Construct this lattice hamiltonian. In the future, remove the embedding 
            !   potential from h0v, and just use coupling
            call AddPeriodicImpCoupling_RealSpace(h_lat_fit,nSites,iLatParams,nImp,Couplings)
            !   For interest, find the 'cluster' greens function, which should now fit the 'bath' greens function
            call FindLocalMomGF(nESteps_Im,zero_GF,Cluster_GF,tMatbrAxis=.true.,ham=h_lat_fit,CouplingLength=iLatParams)
            call writedynamicfunction(nESteps_Im,Cluster_GF,'G_Clust_Im',tag=iter,tMatbrAxis=.true.)
            !Write out the lattice couplings
            call WriteLatticeCouplings(iLatParams,Couplings)
            write(6,"(A)") "Fit successful. Couplings and cluster greens function written to disk."

            !4) Perform impurity calculation
            !   High level calculation on matsubara axis, with no self-energy, and just the fit lattice hamiltonian
            write(6,"(A)") "Starting calculation of high-level, impurity greens function"
            call SchmidtGF_wSE(G_Mat_Im,GFChemPot,zero_GF,nESteps_Im,tMatbrAxis=.true.,ham=h_lat_fit)
            call writedynamicfunction(nESteps_Im,G_Mat_Im,'G_Imp_Im',tag=iter,tMatbrAxis=.true.)
            !   2 = change in 'impurity' greens function
!            call writedynamicfunction(nESteps_Im,G_Mat_Im_Old,'G_Mat_Im_Old_temp',tag=iter,tMatbrAxis=.true.)
!            call writedynamicfunction(nESteps_Im,G_Mat_Im,'G_Mat_Im_temp',tag=iter,tMatbrAxis=.true.)
            DiffImpGF(:,:,:) = G_Mat_Im_Old(:,:,:) - G_Mat_Im(:,:,:)
!            call writedynamicfunction(nESteps_Im,DiffImpGF,'DiffImpGF',tag=iter,tMatbrAxis=.true.)
            DiffImpGF(:,:,:) = DiffImpGF(:,:,:) * dconjg(DiffImpGF(:,:,:))
            AllDiffs(2,iter) = sum(real(DiffImpGF(:,:,:),dp))
            G_Mat_Im_Old(:,:,:) = G_Mat_Im(:,:,:)

            !5) Apply Dysons equation, in order to find the self-energy
            if(iLatticeFitType.eq.2) then
                !Bath_GF is *already* the inverse of the bath greens function.
                !We don't need to invert again
                continue
            else
                !Bath_GF is the bath greens function. Invert it to apply Dysons
                call InvertLocalNonHermGF(nESteps_Im,Bath_GF)
            endif
            !   Invert the impurity greens function
            call InvertLocalNonHermGF(nESteps_Im,G_Mat_Im)  
            !   Apply Dyson
            DeltaSE(:,:,:) = SE_Im(:,:,:)   !Store old SE
            SE_Im(:,:,:) = Damping_SE * (Bath_GF(:,:,:) - G_Mat_Im(:,:,:)) + ((one-Damping_SE)*DeltaSE(:,:,:))
            DeltaSE(:,:,:) = DeltaSE(:,:,:) - SE_Im(:,:,:)  !DeltaSE is now the change to the self-energy 
            call writedynamicfunction(nESteps_Im,SE_Im,'SE_Im',tag=iter,tMatbrAxis=.true.)
            write(6,"(A)") "Self-energy updated via Dyson, and written to disk."

            !6) Calc Stats
            !   1 = residual in lattice fit
            AllDiffs(1,iter) = FinalDist    !The final residual in the lattice fit         
            !   3 = change in self-energy 
            DiffImpGF(:,:,:) = DeltaSE(:,:,:) * dconjg(DeltaSE(:,:,:))
            AllDiffs(3,iter) = sum(real(DiffImpGF(:,:,:),dp))
            
            write(6,"(A)") ""
            write(6,"(A,I7,A)") "***   COMPLETED MACROITERATION ",iter," ***"
            write(6,"(A)") "     Iter.  FitResidual        Delta_GF_Imp     Delta_SE"
            do i = 1,iter
                write(6,"(I7,3G20.13)") i,AllDiffs(1,i),AllDiffs(2,i),AllDiffs(3,i)
            enddo
            write(6,"(A)") ""
            call flush(6)

            if(iter.gt.iMaxIter_MacroFit) then
                write(6,"(A,I9)") "Exiting. Max macroiters hit of: ",iMaxIter_MacroFit
                exit
            endif

            !If the change to the self-energy is small enough, we have converged
            if(AllDiffs(3,iter).lt.dDeltaThresh) then
                write(6,"(A)") "Success! Convergence of self-energy on imaginary axis successful"
                write(6,"(A,G20.13)") "Self-energy changing on imaginary axis by less than: ",dDeltaThresh
                exit
            endif
        enddo

        !Write out final self-energy and lattice couplings
        call writedynamicfunction(nESteps_Im,SE_Im,'Converged_SE_Im',tMatbrAxis=.true.)
        write(6,"(A)") "Final *frequency independent* lattice coupling parameters for bath fit: "
        do i = 1,iLatParams
            write(6,"(I6,G20.13)") i,Couplings(i,1)
        enddo
        write(6,"(A)") "" 
        
        deallocate(SE_Im,G_Mat_Im,zero_GF,Cluster_GF,Lattice_GF,Bath_GF,DeltaSE,DiffImpGF,G_Mat_Im_Old,AllDiffs)

        write(6,"(A)") "Now calculating spectral functions on the REAL axis..."
        call flush(6)
        !First, calculate the real-frequency 'cluster' greens function
        allocate(zero_GF(nImp,nImp,nESteps_Re))
        allocate(Cluster_GF(nImp,nImp,nESteps_Re))
        zero_GF(:,:,:) = zzero
        call FindLocalMomGF(nESteps_Re,zero_GF,Cluster_GF,tMatbrAxis=.false.,ham=h_lat_fit,CouplingLength=iLatParams)
        call writedynamicfunction(nESteps_Re,Cluster_GF,'G_Clust_Re_Final',tMatbrAxis=.false.)

        !...and the real-frequency greens function
        allocate(G_Mat_Re(nImp,nImp,nESteps_Re))
        call SchmidtGF_wSE(G_Mat_Re,GFChemPot,zero_GF,nESteps_Re,tMatbrAxis=.false.,ham=h_lat_fit)
        call writedynamicfunction(nESteps_Re,G_Mat_Re,'G_Imp_Re_Final',tMatbrAxis=.false.)

        deallocate(h_lat_fit,Couplings,G_Mat_Re,Cluster_GF)

    end subroutine SC_FitLat_and_SE_Im

    !Return the residual (dist) for a set of lattice couplings to the impurity, maintaining translational symmetry, to the 
    !impurity greens function. This can come with different weighting factors.
    !This can actually be done on either the real or imaginary axis
    !There are some different options. 
    !   1) Different weighting factors. iFitGFWeighting. 0 = flat, 1 = 1/w, 2 = 1/w^2
    !   2) Minimize sum_w W(w) Tr|G_0(w) - G(w)|^2 or with the inverse functions (which renders it equivalent to fitting the hybridization in some way)
    !   3) Different lattice coupling definitions? (i.e. complex couplings? (probably not), independent bath fits for each site?)

    !NOTE: If LatticeFitType = 2 (i.e. the inverse of the functions), then G_Imp is sent in as the *inverse* of the
    !impurity greens function. Not the greens function itself.
    !Note that this can become expensive - it is a n_sites^2 calculation at best
    !If dSeperateFuncs is present, the difference between the values will be stored for each frequency point
    subroutine CalcLatticeFitResidual(G_Imp,nESteps,Couplings,iNumCoups,dist,tMatbrAxis,dSeperateFuncs, &
        dJacobian,FreqPoints,Weights,dJacobian2)
        implicit none
        integer, intent(in) :: nESteps,iNumCoups
        complex(dp), intent(in) :: G_Imp(nImp,nImp,nESteps)
        real(dp), intent(in) :: Couplings(iNumCoups,nImp)
        real(dp), intent(out) :: dist
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(out), optional :: dSeperateFuncs(nImp,nImp,nESteps)  
        real(dp), intent(out), optional :: dJacobian(nESteps,iNumCoups)  !TODO: This will be wrong for > 1 imp
        real(dp), intent(in), optional :: FreqPoints(nESteps)
        real(dp), intent(in), optional :: Weights(nESteps)
        real(dp), intent(out), optional :: dJacobian2(iNumCoups)
        real(dp), allocatable :: h_lat_fit(:,:)
        complex(dp), allocatable :: SE_Dummy(:,:,:),Lattice_GF(:,:,:),DiffMat(:,:,:)
        real(dp), allocatable :: DiffMatr(:,:,:)
        real(dp) :: Omega,LattWeight
        integer :: i,j,k
        logical :: tNonStandardGrid
        integer, parameter :: iNormPower = 2    !The power of the matrix norm for the residual
        character(len=*), parameter :: t_r='CalcLatticeFitResidual'

        if(present(FreqPoints)) then
            tNonStandardGrid = .true.
        else
            tNonStandardGrid = .false.
        endif
        
        !TODO: Fix this, so that the self-energy is an optional argument
        allocate(SE_Dummy(nImp,nImp,nESteps))
        allocate(Lattice_GF(nImp,nImp,nESteps))

        if(.not.tOptGF_Evals) then
            !We need to calculate the lattice greens function with the attached coupling
            allocate(h_lat_fit(nSites,nSites))
            h_lat_fit(:,:) = h0v(:,:)
            !Add the extended, periodized, non-local couplings to the rest of the lattice
            call AddPeriodicImpCoupling_RealSpace(h_lat_fit,nSites,iNumCoups,nImp,Couplings)

            !Now, calculate the local lattice greens function
            !call writematrix(Couplings,'Couplings',.false.)
            !call writematrix(h_lat_fit,'h_lat_fit',.true.)
            !call FindLocalMomGF(nESteps,SE_Dummy,Lattice_GF,tMatbrAxis=tMatbrAxis,ham=h_lat_fit)
            call FindLocalMomGF(nESteps,SE_Dummy,Lattice_GF,tMatbrAxis=tMatbrAxis,ham=h_lat_fit,CouplingLength=iNumCoups)
            deallocate(h_lat_fit)
        else
            if(tNonStandardGrid) then
                call FindLocalMomGF(nESteps,SE_Dummy,Lattice_GF,tMatbrAxis=tMatbrAxis,CouplingLength=iNumCoups, &
                    evals=Couplings(:,1),FreqPoints=FreqPoints)
            else
                call FindLocalMomGF(nESteps,SE_Dummy,Lattice_GF,tMatbrAxis=tMatbrAxis,CouplingLength=iNumCoups, &
                    evals=Couplings(:,1))
            endif
        endif

        if(iLatticeFitType.eq.2) then
            !If required (i.e. we are fitting the inverses of the functions), invert the lattice greens function
            call InvertLocalNonHermGF(nESteps,Lattice_GF)
        endif
        
        allocate(DiffMat(nImp,nImp,nESteps))
        allocate(DiffMatr(nImp,nImp,nESteps))

        !Now, take the difference between the functions
        DiffMat(:,:,:) = Lattice_GF(:,:,:) - G_Imp(:,:,:)

        !Take abs values 
        DiffMatr(:,:,:) = abs(DiffMat(:,:,:))
!        write(6,*) "Sep_Dists from CalcLatticeFitResidual: "
!        write(6,*) DiffMatr(:,:,:)
!        write(6,*) "enorm of sep_dists: ",enorm(nESteps,DiffMatr(1,1,:))
        if(present(dSeperateFuncs)) then
            !We want to return the individual differences
            dSeperateFuncs(:,:,:) = DiffMatr(:,:,:)
            if(tNonStandardGrid) then
                do i = 1,nESteps
                    !dSeperateFuncs(:,:,i) = Weights(i)*(dSeperateFuncs(:,:,i)**2)
                    dSeperateFuncs(:,:,i) = dSeperateFuncs(:,:,i)*sqrt(Weights(i))
                enddo
            !else
            !    do i = 1,nESteps
            !        dSeperateFuncs(:,:,i) = dSeperateFuncs(:,:,i)**2
            !    enddo
            endif
            if(iFitGFWeighting.ne.0) then
                i = 0
                do while(.true.)
                    if(tNonStandardGrid) then
                        call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis,nFreqPoints=nESteps,FreqPoints=FreqPoints)
                    else
                        call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis)
                    endif
                    if((abs(Omega).lt.1.0e-9_dp).and.(iFitGFWeighting.ne.0)) then
                        call stop_all(t_r,'Should not be sampling w=0 with a non-flat weighting function')
                    endif
                    if(i.lt.0) exit
                    if(i.gt.nESteps) call stop_all(t_r,'Too many frequency points')

                    !Now, include the appropriate weighting factor
                    do j = 1,nImp
                        do k = 1,nImp
                            if(iFitGFWeighting.eq.1) then
                                !1/w weighting
                                !We actually divide by the sqrt of Omega, since the LM optimization will automatically take the squares of the functions
                                dSeperateFuncs(k,j,i) = dSeperateFuncs(k,j,i)/sqrt(abs(Omega))
                                !dSeperateFuncs(k,j,i) = dSeperateFuncs(k,j,i)/(abs(Omega))
                            elseif(iFitGFWeighting.eq.2) then
                                !1/w^2 weighting
                                dSeperateFuncs(k,j,i) = dSeperateFuncs(k,j,i)/abs(Omega)
                                !dSeperateFuncs(k,j,i) = dSeperateFuncs(k,j,i)/(Omega**2)
                            else
                                call stop_all(t_r,'Error here')
                            endif
                        enddo
                    enddo
                enddo
            endif
            if(present(dJacobian)) then
                !Calculate the derivatives too
                if(present(FreqPoints)) then
                    call CalcJacobian(nESteps,iNumCoups,DiffMatr,DiffMat,dJacobian,Couplings(:,1),tMatbrAxis,   &
                        FreqPoints=FreqPoints,Weights=Weights)
                else
                    call CalcJacobian(nESteps,iNumCoups,DiffMatr,DiffMat,dJacobian,Couplings(:,1),tMatbrAxis)
                endif
            endif
        endif

        if(present(dJacobian2)) then
            !Calculate the derivative of dist wrt each eigenvalue
            if(present(FreqPoints)) then
                call CalcJacobian2(nESteps,iNumCoups,DiffMat,dJacobian2,Couplings(:,1),tMatbrAxis,     &
                    FreqPoints=FreqPoints,Weights=Weights)
            else
                call CalcJacobian2(nESteps,iNumCoups,DiffMat,dJacobian2,Couplings(:,1),tMatbrAxis)
            endif
        endif

        dist = zero
        LattWeight = zero
        i = 0
        do while(.true.)
            if(tNonStandardGrid) then
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis,nFreqPoints=nESteps,FreqPoints=FreqPoints)
            else
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis)
            endif
            if((abs(Omega).lt.1.0e-9_dp).and.(iFitGFWeighting.ne.0)) then
                call stop_all(t_r,'Should not be sampling w=0 with a non-flat weighting function')
            endif
            if(i.lt.0) exit
            if(i.gt.nESteps) call stop_all(t_r,'Too many frequency points')

            !Now, sum the squares, with the appropriate weighting factor
            do j = 1,nImp
                do k = 1,nImp
                    if(iFitGFWeighting.eq.0) then
                        !Flat weighting
                        if(tNonStandardGrid) then
                            dist = dist + (DiffMatr(k,j,i)**iNormPower)*Weights(i)
                        else
                            dist = dist + (DiffMatr(k,j,i)**iNormPower)
                        endif
                    elseif(iFitGFWeighting.eq.1) then
                        !1/w weighting
                        if(tNonStandardGrid) then
                            dist = dist + (DiffMatr(k,j,i)**iNormPower)*Weights(i)/abs(Omega)
                        else
                            dist = dist + (DiffMatr(k,j,i)**iNormPower)/abs(Omega)
                        endif
                    else
                        !1/w^2 weighting
                        if(tNonStandardGrid) then
                            dist = dist + (DiffMatr(k,j,i)**iNormPower)*Weights(i)/(Omega**2)
                        else
                            dist = dist + (DiffMatr(k,j,i)**iNormPower)/(Omega**2)
                        endif
                    endif
                    if(iLatticeFitType.eq.3) then
                        LattWeight = LattWeight + (abs(Lattice_GF(k,j,i)))**iNormPower
                    endif
                enddo
            enddo
        enddo

        if(iLatticeFitType.eq.3) then
            !Attempt to maximize lattice weight in the spectral window so that we don't optimize to v high frequencies
            !write(6,*) "Dist = ",dist
            !write(6,*) "Lattweight = ",LattWeight
            dist = dist/LattWeight
        endif

        deallocate(SE_Dummy,DiffMat,Lattice_GF,DiffMatr)

    end subroutine CalcLatticeFitResidual

    !Use a minimization routine to fit the greens functions by adjusting the lattice coupling
    subroutine FitLatticeCouplings(G_Imp,n,Couplings,iNumCoups,FinalErr,tMatbrAxis,FreqPoints,Weights)
        use MinAlgos
        use Levenberg_Marquardt
        use lbfgs
        implicit none
        integer, intent(in) :: n    !Number of frequency points
        integer, intent(in) :: iNumCoups    !Number of independent coupling parameters
        complex(dp), intent(in) :: G_Imp(nImp,nImp,n)   !This may be the inverse
        real(dp), intent(inout) :: Couplings(iNumCoups,nImp)    !The lattice couplings
        real(dp), intent(out) :: FinalErr
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in), optional :: FreqPoints(n),Weights(n)

        real(dp), allocatable :: step(:),vars(:),var(:),FinalVec(:),Jac(:,:)
        real(dp), allocatable :: Freqs_dum(:),Weights_dum(:),l(:),u(:),grad(:),wa(:)
        integer, allocatable :: iwork(:),nbd(:),iwa(:)
        integer :: nop,i,maxf,iprint,nloop,iquad,ierr,iRealCoupNum,nFuncs,j,ierr_tmp,corrs
        real(dp) :: stopcr,simp,rhobeg,rhoend,InitErr,dsave(29),pgtol
        logical :: tfirst,tOptEVals_,tNonStandardGrid
        character(len=60) :: task,csave
        integer :: isave(44)
        logical :: lsave(4)
        logical, parameter :: tConstrainSyms = .true.
        character(len=*), parameter :: t_r='FitLatticeCouplings'

        if(tOptGF_EVals) then
            iRealCoupNum = iNumCoups
            tOptEVals_ = .true.
        else
            if(tEveryOtherCoup) then
                !We actually only have half the number of couplings, since we are not optimizing every other coupling, and
                !just setting it to zero
                if(mod(iNumCoups,2).eq.2) then
                    iRealCoupNum = iNumCoups/2
                else
                    iRealCoupNum = (iNumCoups+1)/2
                endif
            else
                iRealCoupNum = iNumCoups
            endif
            tOptEVals_ = .false.
        endif
            
        !Initialize parameters
        tNonStandardGrid = .false.
        allocate(vars(iRealCoupNum))
        if(tOptEVals_) then
            if(tKPntSymFit) then
                if(tShift_Mesh.and.(iRealCoupNum.ne.(nSites/2))) call stop_all(t_r,'Error here')
                if((.not.tShift_Mesh).and.(iRealCoupNum.ne.((nSites/2)+1))) call stop_all(t_r,'Error here')
            else
                if(iRealCoupNum.ne.nSites) call stop_all(t_r,'Error Here')
            endif
            do i = 1,iRealCoupNum
                vars(i) = Couplings(i,1)
            enddo
            if(present(FreqPoints).and.present(Weights)) then
                tNonStandardGrid = .true.
            endif
        else
            if(tEveryOtherCoup) then
                do i = 1,iRealCoupNum
                    vars(i) = Couplings((2*i)-1,1)
                enddo
            else
                do i = 1,iRealCoupNum
                    vars(i) = Couplings(i,1)
                enddo
            endif
        endif

        if(present(FreqPoints)) then
            tNonStandardGrid = .true.
        else
            tNonStandardGrid = .false.
            allocate(Weights_dum(n))
            allocate(Freqs_dum(n))
        endif

        if(iFitAlgo.eq.4) then
            !Use L-BFGS(-constrained) to fit
            if(.not.tAnalyticDerivs) call stop_all(t_r,'Need analytic derivatives if optimizing with BFGS')

            !We can specify upper and lower bounds on the solutions. Set these to be /pm1000 
            allocate(l(iRealCoupNum))   !Lower limits
            allocate(u(iRealCoupNum))   !Upper limits
            allocate(nbd(iRealCoupNum)) !Type of limit on each variable: 0 - unbounded, 1 - lower only, 2 - both bounds, 3 - upper only
            nbd(:) = 0
            l(:) = -1000.0_dp
            u(:) = 1000.0_dp

            !rhoend = 1.0e12_dp  !Low accuracy
            !rhoend = 1.0e7_dp   !Mid accuracy
            rhoend = 10.0_dp    !High accuracy
            pgtol = 1.0e-5_dp   !Maximum gradient exit criterion

            allocate(grad(iRealCoupNum))
            
            !corrs is the maximum number of variable metric corrections
            !used to define the limited memory matrix.
            corrs = 6

            !Working arrays
            allocate(wa(2*iRealCoupNum*corrs + 5*iRealCoupNum + 11*corrs*corrs + 8*corrs))
            allocate(iwa(3*iRealCoupNum))

            iprint = 25

            task='START'

            do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START')

                call setulb(iRealCoupNum,corrs,vars,l,u,nbd,FinalErr,grad,rhoend,   &
                    pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)

                if(task(1:2).eq.'FG') then
                    !Request the function (FinalErr), and the gradient (grad) at the current parameters (vars)

                    !Compute function and gradient
                    if(tNonStandardGrid) then
                        call MinCoups(vars,FinalErr,n,iRealCoupNum,G_Imp,tMatbrAxis,tNonStandardGrid,   &
                            FreqPoints=FreqPoints,Weights=Weights,dJac=grad)
                    else
                        call MinCoups(vars,FinalErr,n,iRealCoupNum,G_Imp,tMatbrAxis,tNonStandardGrid,dJac=grad)
                    endif

                elseif(task(1:5).eq.'NEW_X') then
                    !Returned with a new iterate
                    !Have we already been through more iterations than we want?
                    if(isave(34).ge.iMaxFitMicroIter) then
                        task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
                    endif
                else
                    exit
                endif
            enddo
                
            write(6,"(A,A)") "Final task: ",task
            write(6,"(A,I8)") "Final number of bfgs iterations: ",isave(30)
            write(6,"(A,I8)") "Total number of f and g evaluations: ",isave(34)
            write(6,"(A,G20.10)") "Final norm of projected gradient: ",dsave(13)
            write(6,"(A,G20.10)") "Final residual after fit: ",FinalErr
            write(6,"(A)") "Lattice parameters: "
            write(6,*) vars(:)
            deallocate(l,u,nbd,grad,wa,iwa) 
                
        elseif(iFitAlgo.eq.3) then
            !Use a Levenberg-Marquardt algorithm for non-linear optimization with either finite-difference
            !or analytic derivatives
            if(iLatticeFitType.eq.3) then
                call stop_all(t_r,'Cannot currently do Levenberg-Marquardt optimization when dividing by norm')
            endif
            rhoend = 1.0e-7_dp  !Convergence criteria
            nFuncs = n*(nImp**2)   !Number of functions
            if(iMaxFitMicroIter.eq.0) then
                maxf = 500 * (nFuncs+1) !Max number of function evaluations
            else
                maxf = iMaxFitMicroIter !Max number of function evaluations
            endif
            write(6,"(A,I9)") "Number of functions: ", nFuncs
                
            allocate(iwork(iRealCoupNum))
            allocate(FinalVec(nFuncs))
                
            ierr_tmp = 1
            if(tNonStandardGrid) then
                call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis,tNonStandardGrid,    &
                    FreqPoints=FreqPoints,Weights=Weights)
            else
                call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis,tNonStandardGrid,    &
                    FreqPoints=Freqs_dum,Weights=Weights_dum)
            endif
!            write(6,*) "Sep_Dists from MinCoups_LM: "
!            write(6,*) FinalVec(:)
            !InitErr = (enorm(nFuncs,FinalVec))   !Because enorm also takes the sqrt
            InitErr = (enorm(nFuncs,FinalVec))**2   !Because enorm also takes the sqrt
            write(6,"(A,G20.10)") "Initial residual before fit: ",InitErr

            if(tAnalyticDerivs) then
                write(6,"(A)") "Optimizing lattice hamiltonian with analytic derivative Levenberg-Marquardt algorithm"

                !Additional memory for Jacobian matrix
                allocate(Jac(nFuncs,iRealCoupNum))

                if(.not.tNonStandardGrid) then
                    call lmder1(MinCoups_LM, nFuncs, iRealCoupNum, vars, FinalVec, Jac, rhoend, ierr, iwork, maxf, n,   &
                        G_Imp, tMatbrAxis, tNonStandardGrid, Freqs_dum, Weights_dum)
                else
                    call lmder1(MinCoups_LM, nFuncs, iRealCoupNum, vars, FinalVec, Jac, rhoend, ierr, iwork, maxf, n,   &
                        G_Imp, tMatbrAxis, tNonStandardGrid,FreqPoints,Weights)
                endif

                !Final calculation of the residual and jacobian at the end
                !ierr_tmp = 2
                !call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis,Jac)
                ierr_tmp = 1
                if(.not.tNonStandardGrid) then
                    call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis,tNonStandardGrid,    &
                        Jac,FreqPoints=Freqs_dum,Weights=Weights_dum)
                else
                    call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis,tNonStandardGrid,    &
                        Jac,FreqPoints=FreqPoints,Weights=Weights)
                endif
                !do i = 1,iRealCoupNum
                !    do j = 1,nFuncs
                !        if(abs(Jac(j,i)).gt.1e-4_dp) then
                !            write(6,*) j,i,Jac(j,i)
                !            call warning(t_r,'Jacobian not zero')
                !        endif
                !    enddo
                !enddo

                deallocate(Jac)
            else
                !Use Finite-difference jacobians
                write(6,"(A)") "Optimizing lattice hamiltonian with numerical derivative Levenberg-Marquardt algorithm"

                !Call LM algorithm
                if(.not.tNonStandardGrid) then
                    call lmdif1(MinCoups_LM , nFuncs, iRealCoupNum, vars, FinalVec, rhoend, ierr, iwork, maxf, n,   &
                        G_Imp, tMatbrAxis, tNonStandardGrid,Freqs_dum,Weights_dum)
                else
                    call lmdif1(MinCoups_LM , nFuncs, iRealCoupNum, vars, FinalVec, rhoend, ierr, iwork, maxf, n,   &
                        G_Imp, tMatbrAxis, tNonStandardGrid,FreqPoints,Weights)
                endif
            endif
                
            !FinalErr = (enorm(nFuncs,FinalVec))
            FinalErr = (enorm(nFuncs,FinalVec))**2
            write(6,"(A,G20.10)") "Final residual after fit: ",FinalErr
            write(6,"(A)") "Lattice parameters: "
            write(6,*) vars(:)
            if(InitErr.lt.FinalErr) call stop_all(t_r,'Residual actually increased during Lev-Mar optimization...?')

            select case(ierr)
                case(0)
                    call stop_all(t_r,'Improper imput parameters to lmdif1')
                case(1)
                    write(6,"(A)") "Levenberg-Marquardt with approximate gradients successful..."
                case(2)
                    write(6,"(A)") "Levenberg-Marquardt with approximate gradients successful..."
                case(3)
                    write(6,"(A)") "Levenberg-Marquardt with approximate gradients successful..."
                case(4)
                    write(6,"(A)") "Levenberg-Marquardt not successful: Objective function is orthogonal " &
                        //"to the columns of the jacobian to machine precision."
                case(5)
                    write(6,"(A)") "Levenberg-Marquardt not successful: Number of calls to fcn has " &
                        //"reached or exceeded max."
                case(6)
                    write(6,"(A)") "Levenberg-Marquardt not successful: tol is too small. No further " &
                        //"reduction in the sum of squares is possible."
                case(7)
                    write(6,"(A)") "Levenberg-Marquardt not successful: tol is too small. No further " &
                        //"improvement in the approximate solution x is possible." 
                case default
                    call stop_all(t_r,'Unknown exit flag')
            end select
            deallocate(iwork,FinalVec)

        elseif(iFitAlgo.eq.2) then
            !Use modified Powell algorithm for optimization (based on fitting to quadratics)
            write(6,"(A)") "Optimizing lattice hamiltonian with modified Powell algorithm"
            rhobeg = 0.05_dp
            rhoend = 1.0e-6_dp
            iprint = 2
            if(iMaxFitMicroIter.eq.0) then
                maxf = 20*iRealCoupNum
            else
                maxf = iMaxFitMicroIter
            endif

            !Vars is updated to be the best value to minimize the residual
            if(tNonStandardGrid) then
                call uobyqa(iRealCoupNum,vars,rhobeg,rhoend,iprint,maxf,FinalErr,MinCoups, n, G_Imp,    &
                    tMatbrAxis, tNonStandardGrid, FreqPoints, Weights)
            else
                call uobyqa(iRealCoupNum,vars,rhobeg,rhoend,iprint,maxf,FinalErr,MinCoups, n, G_Imp,    &
                    tMatbrAxis, tNonStandardGrid, Freqs_dum, Weights_dum)
            endif

        elseif(iFitAlgo.eq.1) then
            !Use simplex method for optimization without derivatives
            write(6,"(A)") "Optimizing lattice hamiltonian with Simplex algorithm"

            !Starting values are assumed to be read in
            allocate(step(iRealCoupNum))
            allocate(var(iRealCoupNum))
            step(:) = 0.05_dp
            nop = iRealCoupNum

            !Set max no of function evaluations. Default = 20*variables, print every 25
            if(iMaxFitMicroIter.eq.0) then
                maxf = 20*iRealCoupNum
            else
                maxf = iMaxFitMicroIter
            endif
            iprint = 25

            !Set value for stopping criterion. Stopping occurs when the 
            !standard deviation of the values of the objective function at
            !the points of the current simplex < stopcr
            stopcr = 1.0e-4_dp
            nloop = 2*iRealCoupNum !This is how often the stopping criterion is checked

            !Fit a quadratic surface to be sure a minimum has been found
            iquad = 1

            !As function value is being evaluated in real(dp), it should
            !be accurate to about 15 dp. If we set simp = 1.d-6, we should
            !get about 9dp accuracy in fitting the surface
            simp = 1.0e-6_dp

            !Now call minim to do the work
            tfirst = .true.
            do while(.true.)
                if(tNonStandardGrid) then
                    call minim(vars, step, iRealCoupNum, FinalErr, maxf, iprint, stopcr, nloop, &
                        iquad, simp, var, MinCoups, ierr, n, G_Imp, tMatbrAxis, tNonStandardGrid, FreqPoints,Weights)
                else
                    call minim(vars, step, iRealCoupNum, FinalErr, maxf, iprint, stopcr, nloop, &
                        iquad, simp, var, MinCoups, ierr, n, G_Imp, tMatbrAxis, tNonStandardGrid, Freqs_dum, Weights_dum)
                endif

                if(ierr.eq.0) exit
                if(.not.tFirst) exit
                tFirst = .false.    !We have found a minimum, but run again with a small number of iterations to check it is stable
                maxf = 3*iRealCoupNum
            enddo

            !On output, Err is the minimized objective function, var contains the diagonal of the inverse of the information matrix, whatever that is
            if(ierr.eq.0) then
                write(6,"(A)") "Simplex optimization successful."
                write(6,"(A,F14.7)") "Minimized residual: ",FinalErr
            elseif(ierr.eq.4) then
                call stop_all(t_r,'nloop < 1')
            elseif(ierr.eq.3) then
                call stop_all(t_r,'iRealCoupNum < 1')
            elseif(ierr.eq.2) then
                call warning(t_r,'Information matrix is not +ve semi-definite')
                write(6,"(A,F14.7)") "Final residual: ",FinalErr
            elseif(ierr.eq.1) then
                call warning(t_r,'Max number of Simplex function evaluations reached.')
                write(6,"(A,F14.7)") "Final residual: ",FinalErr
            endif

            deallocate(step,var)
        endif
            
        if(tOptGF_EVals.and.tConstrainSyms) then
            !Impose momentum inversion, and potentially ph symmetry
            !Symmetrize the resulting eigenvalues
            call Imposephsym(iRealCoupNum,vars)

            !Calculate the new final residual
            if((iFitAlgo.le.2).or.(iFitAlgo.eq.4)) then
                if(tNonStandardGrid) then
                    call MinCoups(vars,FinalErr,n,iRealCoupNum,G_Imp,tMatbrAxis, tNonStandardGrid, FreqPoints, Weights)
                else
                    call MinCoups(vars,FinalErr,n,iRealCoupNum,G_Imp,tMatbrAxis, tNonStandardGrid)
                endif
            elseif(iFitAlgo.eq.3) then
                !LM algorithm
                allocate(FinalVec(nFuncs))
                ierr_tmp = 1
                if(tNonStandardGrid) then
                    call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis, tNonStandardGrid,   &
                        FreqPoints=FreqPoints, Weights=Weights)
                else
                    call MinCoups_LM(nFuncs,iRealCoupNum,vars,FinalVec,ierr_tmp,n,G_Imp,tMatbrAxis, tNonStandardGrid,   &
                        FreqPoints=Freqs_dum, Weights=Weights_dum)
                endif
                FinalErr = (enorm(nFuncs,FinalVec))**2
                !FinalErr = (enorm(nFuncs,FinalVec))
                deallocate(FinalVec)
            endif
            write(6,"(A,G20.10)") "Final residual after symmetrization of eigenvalues: ",FinalErr
            write(6,"(A)") "Lattice after symmetrization of parameters: "
            write(6,*) vars(:)
        endif
            
        !Update couplings
        if(tOptEVals_) then
            !We need to translate these eigenvalues back into couplings?
            do i = 1,nImp
                Couplings(:,i) = vars(:)
            enddo
        else
            if(tEveryOtherCoup) then
                !Put back into the couplings array is every other slot
                Couplings(:,1) = zero
                do i = 1,iRealCoupNum
                    Couplings((2*i)-1,1) = vars(i)
                enddo
                do i = 2,nImp
                    Couplings(:,i) = Couplings(:,1)
                enddo
            else
                do i = 1,nImp
                    Couplings(:,i) = vars(:)
                enddo
            endif
        endif

        deallocate(vars)
        if(allocated(Weights_dum)) deallocate(Weights_dum)
        if(allocated(Freqs_dum)) deallocate(Freqs_dum)

    end subroutine FitLatticeCouplings

    !Wrapper function to evaluate residuals for the Leveberg-Marquardt algorithm
    subroutine MinCoups_LM(n,iNumOptCoups,vars,dists,iflag,nESteps,G,tMatbrAxis,tNonStandardGrid,Jacobian,FreqPoints,Weights)
        implicit none
        integer, intent(in) :: n, iNumOptCoups
        real(dp), intent(in) :: vars(iNumOptCoups)
        real(dp), intent(inout) :: dists(n)
        integer, intent(inout) :: iflag
        integer, intent(in) :: nESteps
        complex(dp), intent(in) :: G(nImp,nImp,nESteps)
        logical, intent(in) :: tMatbraxis
        real(dp), intent(inout), optional :: Jacobian(n,iNumOptCoups)
        logical, intent(in) :: tNonStandardGrid
        real(dp), intent(in), optional :: FreqPoints(n)
        real(dp), intent(in), optional :: Weights(n)

        real(dp), allocatable :: CoupsTemp(:,:),Sep_Dists(:,:,:)
        real(dp) :: dist,Omega,diff,dist2,NumDiff
        real(dp), allocatable :: CoupsTemp2(:,:),Sep_Dists2(:,:,:)
        logical, parameter :: tTest=.false. 
        integer :: i,j,k,ind,RealCoupsNum,ivar,ifreq,iunit2
        character(len=*), parameter :: t_r='MinCoups_LM'

        !Expand variables back to full set
        if(tOptGF_EVals) then
            RealCoupsNum = iNumOptCoups
            allocate(CoupsTemp(RealCoupsNum,nImp))
            CoupsTemp(:,:) = zero
            do i = 1,nImp
                CoupsTemp(:,i) = vars(:)
            enddo

        elseif(tEveryOtherCoup) then
            !Expand the couplings back into the full couplings (with the zeros)
            RealCoupsNum = 2*iNumOptCoups
            allocate(CoupsTemp(RealCoupsNum,nImp))
            CoupsTemp(:,:) = zero
            do i = 1,iNumOptCoups
                CoupsTemp((i*2)-1,1) = vars(i)
            enddo
            do j = 2,nImp
                CoupsTemp(:,j) = CoupsTemp(:,1)
            enddo
        else
            RealCoupsNum = iNumOptCoups
            allocate(CoupsTemp(RealCoupsNum,nImp))
            CoupsTemp(:,:) = zero
            do i = 1,nImp
                CoupsTemp(:,i) = vars(:)
            enddo
        endif

        !n is the number of function points
        if(n.ne.(nESteps*(nImp**2))) then
            call stop_all(t_r,'Wrong number of functions')
        endif

        allocate(Sep_Dists(nImp,nImp,nESteps))

        if(tNonStandardGrid) then
            if(tAnalyticDerivs.and.(iflag.eq.2)) then
                !Calculate analytic derivatives
                if(.not.present(Jacobian)) call stop_all(t_r,'Jacobian argument not present')
                call CalcLatticeFitResidual(G,nESteps,CoupsTemp,RealCoupsNum,dist,tMatbrAxis,dSeperateFuncs=Sep_Dists,  &
                    dJacobian=Jacobian,FreqPoints=FreqPoints,Weights=Weights)
            else
                call CalcLatticeFitResidual(G,nESteps,CoupsTemp,RealCoupsNum,dist,tMatbrAxis,dSeperateFuncs=Sep_Dists,  &
                    FreqPoints=FreqPoints,Weights=Weights)
            endif
        else
            if(tAnalyticDerivs.and.(iflag.eq.2)) then
                !Calculate analytic derivatives
                if(.not.present(Jacobian)) call stop_all(t_r,'Jacobian argument not present')
                call CalcLatticeFitResidual(G,nESteps,CoupsTemp,RealCoupsNum,dist,tMatbrAxis,dSeperateFuncs=Sep_Dists,  &
                    dJacobian=Jacobian)
            else
                call CalcLatticeFitResidual(G,nESteps,CoupsTemp,RealCoupsNum,dist,tMatbrAxis,dSeperateFuncs=Sep_Dists)
            endif
        endif

        if((tAnalyticDerivs.and.iflag.ne.2).or.(.not.tAnalyticDerivs)) then
            !Now, put these values of the functions back into the correct array
            !Potentially in the future, we will want to compress the function array so that there are only nImp*(nImp+1)/2 functions. Perhaps this is not a bottleneck though.
            i = 0
            ind = 0
            do while(.true.)
                if(tNonStandardGrid) then
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis,nFreqPoints=nESteps,FreqPoints=FreqPoints)
                else
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis)
                endif
                if(i.lt.0) exit
                if(i.gt.nESteps) call stop_all(t_r,'Too many frequency points')
                do j = 1,nImp
                    do k = 1,nImp
                        ind = ind + 1
                        if(ind.gt.n) call stop_all(t_r,'More functions than there should be?')
                        dists(ind) = Sep_Dists(k,j,i)
                    enddo
                enddo
            enddo
        endif

        if(tTest.and.(tAnalyticDerivs).and.(iflag.eq.2)) then
            !Test analytic derivatives

            allocate(Sep_Dists2(nImp,nImp,nESteps))
            allocate(CoupsTemp2(iNumOptCoups,nImp))
            !ivar = 1    !Variable to differentiate wrt
            !ifreq = 20   !Frequency number that is being differentiated
            iunit2 = 66

            do ifreq = 1,n  !,4
                do ivar = 1,iNumOptCoups

                    diff = 0.0001_dp
                    do while(diff.gt.1.0e-16_dp)

                        CoupsTemp2(:,:) = CoupsTemp(:,:)
                        CoupsTemp2(ivar,:) = CoupsTemp2(ivar,:) + diff
                        if(tNonStandardGrid) then
                            call CalcLatticeFitResidual(G,nESteps,CoupsTemp2,RealCoupsNum,dist2,tMatbrAxis, &
                                dSeperateFuncs=Sep_Dists2,FreqPoints=FreqPoints,Weights=Weights)
                        else
                            call CalcLatticeFitResidual(G,nESteps,CoupsTemp2,RealCoupsNum,dist2,tMatbrAxis, &
                                dSeperateFuncs=Sep_Dists2)
                        endif
                        NumDiff = ( Sep_Dists2(1,1,ifreq) - Sep_Dists(1,1,ifreq) ) / diff
                        !write(6,*) "***",Sep_Dists2(1,1,ifreq),Sep_Dists(1,1,ifreq)
                        !if(Sep_Dists2(1,1,ifreq).lt.Sep_Dists(1,1,ifreq)) write(6,*) "Negative gradient?"
                        write(iunit2,"(2I8,6G25.14)") ivar,ifreq,diff,NumDiff,Jacobian(ifreq,ivar), &
                            abs((NumDiff-Jacobian(ifreq,ivar))/NumDiff),abs(NumDiff-Jacobian(ifreq,ivar)),Sep_Dists(1,1,ifreq)
                        if(abs((NumDiff-Jacobian(ifreq,ivar))/NumDiff).gt.one) then
                            write(6,"(A,2I8,6G25.14)") "***",ivar,ifreq,diff,NumDiff,Jacobian(ifreq,ivar), &
                                abs((NumDiff-Jacobian(ifreq,ivar))/NumDiff),abs(NumDiff-Jacobian(ifreq,ivar)),Sep_Dists(1,1,ifreq)
                        endif
                        diff = diff/2.0_dp

                    enddo
                    write(iunit2,"(A)") ""
                enddo
            enddo
            deallocate(Sep_Dists2,CoupsTemp2)
            call stop_all(t_r,'end')
        endif

!        write(6,"(10G20.13)") "***",vars(:)
!        write(6,*) "Resid: ",dist
        deallocate(CoupsTemp,Sep_Dists)

    end subroutine MinCoups_LM

    !Wrapper function to evaluate residual for the simplex/Melder-Neal algorithm
    subroutine MinCoups(vars,dist,n,iNumOptCoups,G,tMatbrAxis,tNonStandardGrid,FreqPoints,Weights,dJac)
        implicit none
        integer, intent(in) :: n,iNumOptCoups
        real(dp), intent(in) :: vars(iNumOptCoups)
        real(dp), intent(out) :: dist
        complex(dp), intent(in) :: G(nImp,nImp,n)
        logical, intent(in) :: tMatbrAxis
        logical, intent(in) :: tNonStandardGrid
        real(dp), intent(in), optional :: FreqPoints(n)
        real(dp), intent(in), optional :: Weights(n)
        real(dp), intent(out), optional :: dJac(iNumOptCoups)

        real(dp) :: dist2,diff,NumDiff
        real(dp), allocatable :: CoupsTemp(:,:),CoupsTemp2(:,:)
        integer :: i,j,ivar,iunit2,RealCoupsNum
        logical, parameter :: tTest=.false. 
        character(len=*), parameter :: t_r='MinCoups'

        if(nImp.ne.1) call stop_all(t_r,'n != nfreq')

        if(tOptGF_EVals) then
            RealCoupsNum = iNumOptCoups
            allocate(CoupsTemp(RealCoupsNum,nImp))
            CoupsTemp(:,:) = zero
            do i = 1,nImp
                CoupsTemp(:,i) = vars(:)
            enddo

        elseif(tEveryOtherCoup) then
            !Expand the couplings back into the full couplings (with the zeros)
            RealCoupsNum = 2*iNumOptCoups
            allocate(CoupsTemp(RealCoupsNum,nImp))
            CoupsTemp(:,:) = zero
            do i = 1,iNumOptCoups
                CoupsTemp((i*2)-1,1) = vars(i)
            enddo
            do j = 2,nImp
                CoupsTemp(:,j) = CoupsTemp(:,1)
            enddo
        else
            RealCoupsNum = iNumOptCoups
            allocate(CoupsTemp(RealCoupsNum,nImp))
            CoupsTemp(:,:) = zero
            do i = 1,nImp
                CoupsTemp(:,i) = vars(:)
            enddo
        endif

        if(present(dJac)) then
            !Compute the gradient
            if(tNonStandardGrid) then
                if(.not.present(FreqPoints)) call stop_all(t_r,'Error here')
                call CalcLatticeFitResidual(G,n,CoupsTemp,RealCoupsNum,dist,tMatbrAxis, &
                    FreqPoints=FreqPoints,Weights=Weights,dJacobian2=dJac)
            else
                call CalcLatticeFitResidual(G,n,CoupsTemp,RealCoupsNum,dist,tMatbrAxis,dJacobian2=dJac)
            endif
            if(tTest) then
                !Test analytic derivatives
                allocate(CoupsTemp2(iNumOptCoups,nImp))
                iunit2 = 66

                do ivar = 1,iNumOptCoups

                    diff = 0.0001_dp
                    do while(diff.gt.1.0e-16_dp)

                        CoupsTemp2(:,:) = CoupsTemp(:,:)
                        CoupsTemp2(ivar,:) = CoupsTemp2(ivar,:) + diff
                        if(tNonStandardGrid) then
                            call CalcLatticeFitResidual(G,n,CoupsTemp2,RealCoupsNum,dist2,tMatbrAxis, &
                                FreqPoints=FreqPoints,Weights=Weights)
                        else
                            call CalcLatticeFitResidual(G,n,CoupsTemp2,RealCoupsNum,dist2,tMatbrAxis)
                        endif
                        NumDiff = ( dist2 - dist ) / diff
                        write(iunit2,"(I8,6G25.14)") ivar,diff,NumDiff,dJac(ivar), &
                            abs((NumDiff-dJac(ivar))/NumDiff),abs(NumDiff-dJac(ivar)),dist2
                        if(abs((NumDiff-dJac(ivar))/NumDiff).gt.one) then
                            write(6,"(A,I8,6G25.14)") "***",ivar,diff,NumDiff,dJac(ivar), &
                                abs((NumDiff-dJac(ivar))/NumDiff),abs(NumDiff-dJac(ivar)),dist2
                        endif
                        diff = diff/2.0_dp

                    enddo
                    write(iunit2,"(A)") ""
                enddo
                deallocate(CoupsTemp2)
                call stop_all(t_r,'end')
            endif
        else
            if(tNonStandardGrid) then
                if(.not.present(FreqPoints)) call stop_all(t_r,'Error here')
                call CalcLatticeFitResidual(G,n,CoupsTemp,RealCoupsNum,dist,tMatbrAxis, &
                    FreqPoints=FreqPoints,Weights=Weights)
            else
                call CalcLatticeFitResidual(G,n,CoupsTemp,RealCoupsNum,dist,tMatbrAxis)
            endif
        endif

        deallocate(CoupsTemp)
    end subroutine MinCoups

    subroutine WriteLatticeCouplings(iLatParams,Couplings)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: Couplings(iLatParams,nImp)
        integer :: iunit,i
        character(len=*), parameter :: t_r='WriteLatticeCouplings'

        if(nImp.gt.1) call stop_all(t_r,'Error here')

        iunit = get_free_unit()
        if(tOptGF_EVals) then
            open(unit=iunit,file='LatEVals',status='unknown')
        else
            open(unit=iunit,file='LatCouplings',status='unknown')
        endif
        write(iunit,*) iLatParams
        if(tDiag_kSpace) then
            do i = 1,iLatParams
                write(iunit,"(2G25.13)") KPnts(1,i),Couplings(i,1)
            enddo
        else
            do i = 1,iLatParams
                write(iunit,"(G25.13)") Couplings(i,1)
            enddo
        endif
        close(iunit)

    end subroutine WriteLatticeCouplings
        
    subroutine InitLatticeCouplings(Couplings,iLatParams,h_lat_fit)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(out) :: Couplings(iLatParams,nImp)
        real(dp), intent(out) :: h_lat_fit(nSites,nSites)
        real(dp), allocatable :: lat_evals(:)
        integer :: iunit,iloclatcoups,i
        logical :: exists
        character(len=*), parameter :: t_r='InitLatticeCouplings'

        Couplings(:,:) = zero
        !h_lat_fit(:,:) = h0(:,:)
        h_lat_fit(:,:) = h0v(:,:)

        if(tReadCouplings) then
            if(tOptGF_EVals) then
                call ReadLatticeEigenvalues(Couplings,iLatParams)
            else
                iunit = get_free_unit()
                inquire(file='LatCouplings_Read',exist=exists)
                if(.not.exists) then
                    call stop_all(t_r,'LatCouplings_Read file does not exist')
                endif
                open(unit=iunit,file='LatCouplings_Read',status='old')
                read(iunit,*) iloclatcoups  !This does not need to be the same as iLatParams
                iloclatcoups = min(iLatParams,iloclatcoups)
                do i = 1,iloclatcoups
                    read(iunit,*) Couplings(i,1) !Fill the coupslings
                enddo
                if(tEveryOtherCoup) then
                    !Set every other lattice coupling to zero
                    do i = 2,iloclatcoups,2
                        Couplings(i,1) = zero
                    enddo
                endif
                do i = 2,nImp
                    !Copy them to the other impurities
                    Couplings(:,i) = Couplings(:,1)
                enddo
                write(6,*) "Lattice couplings read in from disk, and initialised to: ",Couplings(:,1)
                close(iunit)
            endif
        else
            if(tOptGF_EVals) then
                !Here, we are optimizing the eigenvalues, rather than the lattice couplings
                if(nImp.gt.1) call stop_all(t_r,'Not working for multiple impurity sites yet')
                if(tkPntSymFit) then
                    if(tShift_Mesh.and.(iLatParams.ne.(nSites/2))) call stop_all(t_r,'Wrong # of params')
                    if((.not.tShift_Mesh).and.(iLatParams.ne.((nSites/2)+1))) call stop_all(t_r,'Wrong # of params')
                else
                    if(iLatParams.ne.nSites) call stop_all(t_r,'Wrong # of params')
                endif
                !The eigenvalues are returned in *increasing* k-point order
                allocate(lat_evals(nSites))
                call CouplingsToEVals(h_lat_fit,lat_evals)
                !call writevector(lat_evals,'lat_evals')
                Couplings(:,1) = lat_evals(1:iLatParams)
                do i = 2,nImp
                    Couplings(:,i) = Couplings(:,1)
                enddo
                deallocate(lat_evals)
            else
                !This is the normal nearest neighbour coupling in the hubbard model. We will start from this
                Couplings(1,:) = -1.0_dp    
            endif
        endif

    end subroutine InitLatticeCouplings
                
    subroutine ReadLatticeEigenvalues(Couplings,iLatParams)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(out) :: Couplings(iLatParams,nImp)
        integer :: iunit,iloclatcoups,i,ineghalf
        real(dp) :: kval,GammaBound,yp1,ypn
        real(dp), allocatable :: kmesh_read(:),evals_read(:),y2_spline(:),neghalf_evals(:)
        logical :: exists
        character(len=*), parameter :: t_r='ReadLatticeEigenvalues'

        Couplings(:,:) = zero
        if(nImp.gt.1) call stop_all(t_r,'Not working for multiple impurity sites')
        if(LatticeDim.ne.1) call stop_all(t_r,'Currently only working for 1D systems')

        inquire(file='LatEVals_Read',exist=exists)
        if(.not.exists) then
            call stop_all(t_r,'LatEVals_Read file does not exist')
        endif
        iunit = get_free_unit()
        open(unit=iunit,file='LatEVals_Read',status='old')
        read(iunit,*) iloclatcoups

        if(iloclatcoups.gt.iLatParams) then
            call stop_all(t_r,'Cannot read in from larger lattice / check inputs compatible')
        elseif(iloclatcoups.eq.iLatParams) then
            write(6,"(A)") "In lattice read, number of free parameters the same. Assuming that the "    &
                //"k-point mesh is commensurate, and just using these values."

            if(tKPntSymFit) then
                if(tShift_Mesh.and.(iloclatcoups.ne.(nSites/2))) call stop_all(t_r,'Wrong # of params read in')
                if((.not.tShift_Mesh).and.(iloclatcoups.ne.((nSites/2)+1))) call stop_all(t_r,'Wrong # of params read in')
            else
                if(iloclatcoups.ne.nSites) call stop_all(t_r,'Wrong # of params read in')
            endif
            if(tDiag_kSpace) then
                do i = 1,iloclatcoups
                    read(iunit,*) kval,Couplings(i,1)
                enddo
            else
                do i = 1,iloclatcoups
                    read(iunit,*) Couplings(i,1)
                enddo
            endif
            !Potentially shift the lattice eigenvalues (required if reading in from a previous calculation with a different chemical potential)
            do i = 1,iloclatcoups
                Couplings(i,1) = Couplings(i,1) + dShiftLatticeEvals
            enddo
        else
            !Reading in from a smaller lattice, and fitting to the kpoint mesh of a larger lattice via cubic splines
            allocate(kmesh_read(iloclatcoups))
            allocate(evals_read(iloclatcoups))

            do i = 1,iloclatcoups
                read(iunit,*), kmesh_read(i),evals_read(i)
                if(i.gt.1) then
                    !Check that the kpoints are written out in increasing order
                    if(kmesh_read(i).lt.kmesh_read(i-1)) then
                        call stop_all(t_r,'Written out eigenvalues are not in ascending k-point order')
                    endif
                endif
            enddo

            GammaBound = zero+1.0e-8
            do i = 1,iloclatcoups
                if(kmesh_read(i).gt.GammaBound) then
                    !This point is in the positive k half of the BZ. Ignore it for the fitting.
                    exit
                endif
            enddo
            ineghalf = i-1  !Negative k-values are found in the values 1:ineghalf
            write(6,"(A,I6)") "Number of kpoints read in from negative half of BZ: ",ineghalf

            !Set the *first* derivative of the spline to be zero at the BZ boundary and gamma point
            yp1 = zero
            ypn = zero
            !Set the *second* derivative of the spline to be zero at the BZ boundary and gamma point
            !yp1 = 1.0e30
            !ypn = 1.0e30
            allocate(y2_spline(ineghalf))
            call spline(kmesh_read(1:ineghalf),evals_read(1:ineghalf),ineghalf,yp1,ypn,y2_spline)

            !Now, interpolate for the eigenvalues on the new k-point mesh
            do i = 1,nKPnts
                if(KPnts(1,i).lt.kmesh_read(1)) then
                    !We have a point outside the spline. Set it to be the same
                    Couplings(i,1) = evals_read(1)
                elseif((KPnts(1,i).gt.kmesh_read(1)).and.(KPnts(1,i).lt.kmesh_read(ineghalf))) then
                    !Interpolate
                    call splint(kmesh_read(1:ineghalf),evals_read(1:ineghalf),y2_spline,ineghalf,KPnts(1,i),Couplings(i,1))
                elseif((KPnts(1,i).gt.kmesh_read(ineghalf)).and.(KPnts(1,i).lt.GammaBound)) then
                    !Outside point on the spline (but still in the negative k half of the BZ. Set to be same as final point
                    Couplings(i,1) = evals_read(ineghalf)
                else
                    !In positive half of BZ. Ignore - we will make them k-symmetric later on if necessary
                    exit
                endif
            enddo

            !Finally, if necessary, k-symmetrize
            if(tKPntSymFit.and.(iLatParams.ne.i-1)) then
                call stop_all(t_r,'Wrong number of kpoints interpolated?')
            endif
            if(.not.tKPntSymFit) then
                !We need to mirror the kpoints to the other half of the plane
                allocate(neghalf_evals(ineghalf))
                neghalf_evals(:) = Couplings(1:ineghalf,1)
                Couplings(:,:) = zero
                call FindFullLatEvals(neghalf_evals,ineghalf,Couplings(:,1))
                deallocate(neghalf_evals)
            endif
            deallocate(y2_spline)
        endif

        close(iunit)

        !potentially, impose ph symmetry in the read in eigenvalues (the interpolation may have broken it)
        call Imposephsym(iLatParams,Couplings(:,1))

        do i = 2,nImp
            !Copy them to the other impurities
            Couplings(:,i) = Couplings(:,1)
        enddo
        write(6,*) "Lattice eigenvalues read in from disk, and initialised to: ",Couplings(:,1)

    end subroutine ReadLatticeEigenvalues

    !Write out the bandstructure from the fit one-electron hamiltonian
    subroutine WriteBandstructure(evals,nvals)
        implicit none
        integer, intent(in) :: nvals
        real(dp), intent(in) :: evals(nvals)
        real(dp), allocatable :: evals_full(:)
        integer :: iunit
        real(dp) :: yp1,ypn,kval,band
        real(dp), allocatable :: y2_spline(:)
        character(len=*), parameter :: t_r='WriteBandstructure'

        if(nImp.gt.1) call stop_all(t_r,'Cannot currently deal with > 1 impurity')

        write(6,"(A)") "Writing to disk the bandstructure from the fit lattice hamiltonian"

        iunit = get_free_unit()
        open(unit=iunit,file='Bandstructure',status='unknown')

        allocate(evals_full(nKPnts))
        call FindFullLatEvals(evals,nvals,evals_full)

        !Set the *first* derivative of the spline to be zero at the BZ boundarys
        yp1 = zero
        ypn = zero
        !Set the *second* derivative of the spline to be zero at the BZ boundarys
        !yp1 = 1.0e30
        !ypn = 1.0e30
        allocate(y2_spline(nKPnts))
        call spline(KPnts(1,:),evals_full,nKPnts,yp1,ypn,y2_spline)

        kval = KPnts(1,1)
        do while(.true.)
            if(kval.gt.KPnts(1,nKPnts)) exit
                    
            call splint(KPnts(1,:),evals_full,y2_spline,nKPnts,kval,band)

            write(iunit,"(2G20.10)") kval,band
            kval = kval + 0.001
        enddo

        deallocate(y2_spline)
        close(iunit)

    end subroutine WriteBandstructure



    !From num rec
    subroutine splint(xa,ya,y2a,n,x,y)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: x, xa(n),y2a(n),ya(n)
        real(dp), intent(out) :: y
        !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa_i's in order),
        !and given the array y2a(1:n), which is the output from spline, and given a value of x, this routine
        !returns a cubic-spline interpolated value y.
        integer :: k,khi,klo
        real(dp) :: a, b, h
        character(len=*), parameter :: t_r='splint'

        klo = 1
        khi = n
1       if (khi-klo.gt.1) then
            !Bisection to find right place in table
            k = (khi+klo)/2
            if(xa(k).gt.x) then
                khi = k
            else
                klo = k
            endif
            goto 1
        endif
        !klo and khi now bracket the value of x
        h = xa(khi)-xa(klo)
        if(h.eq.zero) call stop_all(t_r,'bad xa input in splint. kpoints should be distinct')
        !Evaluate spline
        a = (xa(khi)-x)/h
        b = (x-xa(klo))/h
        y = a*ya(klo)+b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
    end subroutine splint


    !From numerical recipies
    subroutine spline(x,y,n,yp1,ypn,y2)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: yp1,ypn,x(n),y(n)
        real(dp), intent(out) :: y2(n)
        !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e. y_i = f(x_i), with
        !x_1 < x_2 < x_3 < ... < x_n, and given values yp1 and ypn for the first derivative of the
        !interpolating function at points 1 and n, respectively, this routine returns an array y2(1:n)
        !of length n which contains the second derivatives of the interpolating function at the tabulated
        !points x_i. If yp1 and/or ypn are > 1e30, the routine is signaled to set the corresponding boundary
        !condition for a natural spline, with zero second derivative on that boundary.
        integer :: i,k
        real(dp) :: p,qn,sig,un
        real(dp), allocatable :: u(:)

        allocate(u(n))
        if(yp1.gt.0.99e30_dp) then
            !The lower boundary condition is set to be 'natural'
            y2(1) = zero
            u(1) = zero
        else
            !A specified first derivative
            y2(1)=-0.5_dp
            u(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        do i=2,n-1
            !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temp storage of the decomposed factors
            sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
            p=sig*y2(i-1)+2.0_dp
            y2(i)=(sig-one)/p
            u(i) = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        enddo
        if(ypn.gt.0.99e30_dp) then
            !Upper boundary condition is 'natural'
            qn = zero
            un = zero
        else
            !Specified first derivative
            qn = 0.5_dp
            un = (3.0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)
        do k=n-1,1,-1
            !This is the backsubstitution loop of the tridiagonal algorithm
            y2(k) = y2(k)*y2(k+1)+u(k)
        enddo
        deallocate(u)
    end subroutine spline
        

    subroutine Imposephsym(nCoups,vars)
        use sort_mod, only: sort_real
        implicit none
        integer, intent(in) :: nCoups
        real(dp), intent(inout) :: vars(nCoups)
        integer :: i,j,k,EOrder(nCoups),nSpec
        logical :: tNotTaken(nCoups)
        real(dp) :: diff,mu,dCurrentE,av
        character(len=*), parameter :: t_r='Imposephsym'

        if((.not.tKPntSymFit).and.tDiag_kSpace) then
            !After fitting, add back in momentum inversion symmetry (if we are working in k-space)
            !This means that e(k) = e(-k)
            !For shifted meshes, this is easy.
            !For Gamma-centered meshes, two k-points are only singly degenerate
            if(nCoups.ne.nSites) call stop_all(t_r,'Wrong number of variables')

            if(tShift_Mesh) then
                do i = 1,nSites/2
                    av = (vars(i) + vars(nSites-i+1))/2.0_dp
                    vars(i) = av
                    vars(nSites-i+1) = av
                enddo
            else
                !Do not average the gamma point, or the first kpoint (BZ boundary) - these are singly degenerate
                do i = 2,nSites/2
                    av = (vars(i) + vars(nSites-i+2))/2.0_dp
                    vars(i) = av
                    vars(nSites-i+2) = av
                enddo
            endif
        endif

        if(tphsym) then

            if(.not.tDiag_kSpace) then
                !Impose ph symmetry on the system
                !In this way, we only consider the lowest n/2 eigenvalues, and we reflect them around the chemical potential
                !We cannot do this if we are doing k-space diagonalizations, since the vectors are ordered by k-label, rather than energetically
                mu = U/2.0_dp
                call sort_real(vars,nSites)
                j = 1
                do i = nSites/2+1,nSites
                    vars(i) = 2*mu - vars(i-j)
                    j = j + 2
                enddo

            else
                !Impose ph symmetry on the system
                !In this way of doing it, we pick pairs which are closest, and average them
                mu = U/2.0_dp

                if(tShift_Mesh) then
                    nSpec = nSites/2
                else
                    nSpec = nSites/2 + 1
                endif
                if(nSpec.gt.nCoups) call stop_all(t_r,'Error here')
                if(mod(nSpec,2).ne.0) call stop_all(t_r,'nSpec should be even?')
                !Since we have already imposed k-point sym, we can just optimize within the restricted set we already have here,
                !and ignore the fact that we have already got doubly degenerate sets everywhere.
                !Do this *incredibly* crudely. Search for nearest 'pairs', and symmetrize them

    !            write(6,*) "nSpec: ",nSpec
    !            write(6,*) "nCoups: ",nCoups
    !            write(6,*) "tShift_Mesh: ",tShift_Mesh
    !            write(6,*) "vars: ",vars(:)

                !Find the energetic ordering of the eigenvalues
                tNotTaken(:) = .true.
                do i = 1,nSpec
                    !Each loop, find the position of the next lowest eigenvalue
                    dCurrentE = huge(0.0_dp)
                    do j = 1,nSpec
                        !Find lowest energy eigenvalue not already taken
                        if((vars(j).lt.dCurrentE).and.(tNotTaken(j))) then
                            EOrder(i) = j
                            dCurrentE = vars(j)
                        endif
                    enddo
                    tNotTaken(EOrder(i)) = .false.
                enddo
                !We now have the energetic ordering. EOrder(i) gives the ith lowest energy eigenvalue
    !            write(6,*) "EOrder: ",EOrder(:)
    !            write(6,*) "tNotTaken: ",tNotTaken(:)

                do i = 1,nSpec
                    if(tNotTaken(i)) then
                        call writevector(vars,'vars')
                        call stop_all(t_r,'Not all evals accounted for')
                    endif
                enddo

                do i = 2,nSpec
                    if(vars(EOrder(i)).lt.vars(EOrder(i-1))) then
                        call writevector(vars,'vars')
                        call stop_all(t_r,'Energy order not correct')
                    endif
                enddo

                !Now, symmetrize them
                do i = 1,nSpec/2
                    !Run through all occupied orbitals by going up in energy
                    do j = 1,nSpec
                        if(EOrder(j).eq.i) exit
                    enddo
                    !i'th lowest energy eigenvalue is in slot j
                    !Find the corresponding 'virtual'
                    do k = 1,nSpec
                        if(EOrder(k).eq.(nSpec-i+1)) then
                            !This is the corresponding 'virtual'
                            exit
                        endif
                    enddo
                    !j and k are pairs. Average the differences from the chemical potential
                    diff = 0.5_dp*( (mu - vars(j) ) + (vars(k) - mu))
                    vars(j) = mu - diff
                    vars(k) = mu + diff
                enddo

                if((.not.tKPntSymFit).and.tDiag_kSpace) then
                    !We might have just broken the momentum symmetry that we just imposed!
                    !Reimpose it - not by averaging, but by simply using the negative k-point values.
                    if(tShift_Mesh) then
                        do i = 1,nSites/2
                            vars(nSites-i+1) = vars(i)
                        enddo
                    else
                        !Do not average the gamma point, or the first kpoint (BZ boundary) - these are singly degenerate
                        do i = 2,nSites/2
                            vars(nSites-i+2) = vars(i)
                        enddo
                    endif

                endif
            endif
        endif

    end subroutine Imposephsym

    !Attempt for self-consistency on the matsubara axis via simple application of dysons equation
    subroutine SC_Imaginary_Dyson()
        implicit none
        real(dp) :: GFChemPot,Omega
        integer :: nESteps_Im,nESteps_Re,OmegaVal,iter
        complex(dp), allocatable :: SE_Im(:,:,:),SE_Old_Im(:,:,:),G_Mat_Im(:,:,:)
        complex(dp), allocatable :: G_Mat_Re(:,:,:),SE_Re_Err(:,:,:),SE_Re(:,:,:)
        complex(dp), allocatable :: LatticeGF_Im(:,:,:)
        real(dp) :: MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
        real(dp) :: dChangeSE_Tol
        
        dChangeSE_Tol = 1.0e-7_dp

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp

        !How many frequency points are there exactly?
        nESteps_Im = 0
        nESteps_Re = 0
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.true.)
            !write(6,*) "Counting: ",OmegaVal,Omega
            if(OmegaVal.lt.0) exit
            nESteps_Im = nESteps_Im + 1
        enddo
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.false.)
            if(OmegaVal.lt.0) exit
            nESteps_Re = nESteps_Re + 1
        enddo
        write(6,*) "Number of real frequency points: ",nESteps_Re
        write(6,*) "Number of Matsubara frequency points: ",nESteps_Im
        allocate(SE_Im(nImp,nImp,nESteps_Im))
        allocate(SE_Old_Im(nImp,nImp,nESteps_Im))
        allocate(G_Mat_Im(nImp,nImp,nESteps_Im))
        call Init_SE(nESteps_Im,SE_Im)
        call writedynamicFunction(nESteps_Im,SE_Im,'SE_Im',tag=0,tMatbrAxis=.true.)

        allocate(LatticeGF_Im(nImp,nImp,nESteps_Im))

        iter = 0
        !Very simple attempted convergence of the self energy, via iterative applications of Dyson's equation on the imaginary axis
        do while(.true.)

            iter = iter + 1

            !High level calculation on matsubara axis
            call SchmidtGF_wSE(G_Mat_Im,GFChemPot,SE_Im,nESteps_Im,tMatbrAxis=.true.)
            call writedynamicfunction(nESteps_Im,G_Mat_Im,'G_Imp_Im',tag=iter,tMatbrAxis=.true.)

            !Calculate lattice greens function with imaginary self-energy on matsubara axis
            call FindLocalMomGF(nESteps_Im,SE_Im,LatticeGF_Im,tMatbrAxis=.true.)
            call writedynamicfunction(nESteps_Im,LatticeGF_Im,'G_0_Im',tag=iter,tMatbrAxis=.true.)

            !Now update the self-energy via Dysons equation
            SE_Old_Im(:,:,:) = SE_Im(:,:,:)
            call Calc_SE_Update(nESteps_Im,SE_Im,G_Mat_Im,LatticeGF_Im)
            call writedynamicFunction(nESteps_Im,SE_Im,'SE_Im',tag=iter,tMatbrAxis=.true.)

            call FindSEDiffs(SE_Im,SE_Old_Im,G_Mat_Im,LatticeGF_Im,nESteps_Im,MaxDiffSE,MinDiffSE,MeanDiffSE,   &
                MaxDiffGF,MeanDiffGF)

            !Not strictly converged, but if we want to change parameters...
            call WriteConvSE(nESteps_Im,SE_Im,tMatbrAxis=.true.)

            write(6,'(A,I6,5F14.8)') "Iter finished: ",iter,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
            if((MaxDiffSE/damping_SE).lt.dChangeSE_Tol) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",dChangeSE_Tol
                exit
            endif
        enddo

        deallocate(SE_Old_Im)

        !Try analytically continuing the self energy
        allocate(SE_Re(nImp,nImp,nESteps_Re))
        allocate(SE_Re_Err(nImp,nImp,nESteps_Re))
        call FromMatsubaraToRealFreq(nESteps_Im,nESteps_Re,SE_Im,SE_Re,SE_Re_Err)
        call writedynamicfunction(nESteps_Re,SE_Re,'SE_Re_cont',tMatbrAxis=.false.,ErrMat=SE_Re_Err,    &
            tCheckCausal=.true.,tCheckOffDiagHerm=.true.,tWarn=.true.)

        !Now, use this analytically continued self energy to compute the real frequency spectral function 
        allocate(G_Mat_Re(nImp,nImp,nESteps_Re))
        call SchmidtGF_wSE(G_Mat_Re,GFChemPot,SE_Re,nESteps_Re,tMatbrAxis=.false.)
        call writedynamicfunction(nESteps_Re,G_Mat_Re,'G_Imp_Final',tMatbrAxis=.false.,tCheckCausal=.true., &
            tCheckOffDiagHerm=.true.,tWarn=.true.)

        deallocate(SE_Re,SE_Re_Err,SE_Im,G_Mat_Re,G_Mat_Im,LatticeGF_Im)

    end subroutine SC_Imaginary_Dyson

    subroutine Calc_SE_Update(n,SE,G_Imp,G_Lat)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(inout) :: SE(nImp,nImp,n)
        complex(dp), intent(in) :: G_Imp(nImp,nImp,n)
        complex(dp), intent(in) :: G_Lat(nImp,nImp,n)

        complex(dp), allocatable :: InvG_Imp(:,:,:),InvG_Lat(:,:,:)

        !First, invert both greens functions
        allocate(InvG_Imp(nImp,nImp,n))
        InvG_Imp(:,:,:) = G_Imp(:,:,:)
        call InvertLocalNonHermGF(n,InvG_Imp)

        allocate(InvG_Lat(nImp,nImp,n))
        InvG_Lat(:,:,:) = G_Lat(:,:,:)
        call InvertLocalNonHermGF(n,InvG_Lat)

        !Now, apply dysons equation
        SE(:,:,:) = SE(:,:,:) + (Damping_SE * ( InvG_Lat(:,:,:) - InvG_Imp(:,:,:) ))

        deallocate(InvG_Imp,InvG_Lat)

    end subroutine Calc_SE_Update


    !Real space self-consistency without dyson equation
    subroutine SC_Mom_LR_NoDyson()
        implicit none
        real(dp) :: Omega,GFChemPot,OmegaVal
        integer :: nESteps,iter,i,k,l,j
        complex(dp), allocatable :: SE(:,:,:),G_Mat(:,:,:),SE_Old(:,:,:),LocalMomGF(:,:,:)
        real(dp) :: dChangeSE_Tol,dFuncTol,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
        logical, allocatable :: tOmegaConv(:)
        logical :: tSuccess
        character(len=*), parameter :: t_r='SC_Mom_LR_NoDyson'

        dChangeSE_Tol = 1.0e-9_dp
        dFuncTol = 1.0e-8_dp
        
        !How many frequency points are there exactly?
        Omega = Start_Omega
        nESteps = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            nESteps = nESteps + 1
            Omega = Omega + Omega_Step
        enddo
        
        !Find the k-independent self-consistent self-energy
        allocate(SE(nImp,nImp,nESteps))
        allocate(SE_Old(nImp,nImp,nESteps))
        allocate(G_Mat(nImp,nImp,nESteps))
        allocate(LocalMomGF(nImp,nImp,nESteps))

        !Array of all frequency points, and whether they have correctly converged or not
        allocate(tOmegaConv(nESteps))

        call Init_SE(nESteps,SE)

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp
        iter = 0

        do while(.true.)

            iter = iter + 1

            write(6,*) "Starting iteration: ",iter

            !First, do high-level calculation, with SE for bath construction
            call SchmidtGF_wSE(G_Mat,GFChemPot,SE,nESteps)
            write(6,*) "1"
            call flush(6)
            call writedynamicfunction(nESteps,G_Mat,'G_Imp',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)
            write(6,"(A)") "Impurity Greens function causal"

            if(tSE_Scan) then
                call PlotSigmaSurface(nESteps,G_Mat,GFChemPot)
                call stop_all(t_r,'Finished Scanning Self-energy')
            endif

            !Now, numerically calculate the self energy, st. the lattice and impurity greens function match
            SE_Old(:,:,:) = SE(:,:,:)
            call ConvergeGlobalCoupling(nESteps,G_Mat,SE,GFChemPot,dFuncTol,dChangeSE_Tol,tOmegaConv,tSuccess,150)
            write(6,*) "2"
            call flush(6)
            !Now, simply set the values of the frequency outside the causal range to simply be equal to the largest frequency causal self-energy.
            call fit_noncausal_SE(nESteps,SE,tOmegaConv,tSuccess)
            !if(.not.tSuccess) then
                !Turn to a simplex algorithm to converge the points that failed with newton-raphson.
                !call ConvergeGlobalCoupling_Direct(nESteps,G_Mat,SE,GFChemPot,dFuncTol,dChangeSE_Tol,tOmegaConv,tSuccess)
            !endif

            !Now check whether we are actually converged, i.e. does the lattice greens function match the impurity greens function
            !Construct FT of k-space GFs and take zeroth part.
            call FindLocalMomGF(nESteps,SE,LocalMomGF)
            !call FindRealSpaceLocalMomGF(nESteps,SE,LocalMomGF)
            write(6,*) "3"
            call flush(6)
            call writedynamicfunction(nESteps,LocalMomGF,'LocalMomGF',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)
            write(6,"(A)") "Lattice greens function causal"
            !call CheckGFsSame(nESteps,G_Mat,LocalMomGF,1.0e-8_dp)
            !write(6,"(A)") "Lattice GF = Impurity GF"

            call writedynamicfunction(nESteps,SE,'SelfEnergy',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)
            write(6,"(A)") "Self energy causal"

            !Find the change in self-energy (and greens functions)
            call FindSEDiffs(SE,SE_Old,G_Mat,LocalMomGF,nESteps,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF)

            write(6,"(I7,4G15.7)") iter,MaxDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF

            if(MaxDiffSE.lt.dChangeSE_Tol) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",dChangeSE_Tol
                exit
            endif

            !if(iter.eq.1) call stop_all(t_r,'SELF ENERGY CAUSAL (in first iteration)! YAY!')
        enddo


        deallocate(SE,SE_Old,G_Mat,LocalMomGF,tOmegaConv)

    end subroutine SC_Mom_LR_NoDyson

    !Write out a converged self-energy for potential re-reading in later
    subroutine WriteConvSE(n,SE,tMatbrAxis)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        logical, intent(in), optional :: tMatbrAxis
        logical :: tMatbrAxis_
        integer :: iunit,i,k,l
        real(dp) :: Omega
        character(64) :: filename,filename2

        if(present(tMatbrAxis)) then
            tMatbrAxis_ = tMatbrAxis
        else
            tMatbrAxis_ = .false.
        endif
        
        write(6,"(A)") "Writing out converged self-energy"
        iunit = get_free_unit()
        if(tMatbrAxis_) then
            call append_ext_real('Converged_SE_Im',U,filename)
        else
            call append_ext_real('Converged_SE',U,filename)
        endif
        if(.not.tHalfFill) then
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')

        i = 0
        do while(.true.)
            call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_)
            if(i.lt.0) exit
            !Write out *lower* triangle
            write(iunit,"(G25.10)",advance='no') Omega
            do k = 1,nImp
                do l = k,nImp
                    if((k.eq.nImp).and.(l.eq.nImp)) then
                        exit
                    endif
                    write(iunit,"(2G25.10)",advance='no') real(SE(l,k,i),dp),aimag(SE(l,k,i))
                enddo
            enddo
            write(iunit,"(2G25.10)") real(SE(nImp,nImp,i),dp),aimag(SE(nImp,nImp,i))
        enddo
        close(iunit)

    end subroutine WriteConvSE
    
    !Attempt to get k-space spectral functions by self-consistenly calculating k-independent hybridization and self-energy contributions
    !As opposed to the routine below, this one attempts to use a global hybridization in the construction of the bath, rather than the 
    !self-energy
    subroutine SC_Mom_LR_Z()
        implicit none
        real(dp) :: Omega,GFChemPot,OmegaVal,MaxDiffSE,DiffSE,reSE,imSE,r(3)
        complex(dp) :: DiffMatSE(nImp,nImp)
        integer :: nESteps,iter,i,k,l,iunit,j
        complex(dp), allocatable :: SE(:,:,:),G_Mat(:,:,:),SE_Old(:,:,:)
        complex(dp), allocatable :: Hybrid(:,:,:),InvLocalMomGF(:,:,:),LocalMomGF(:,:,:)
        complex(dp), allocatable :: LocalCoupFn(:,:,:),GlobalCoup(:,:,:)
        logical, allocatable :: tOmegaConv(:)
        logical :: tSuccess
        character(64) :: filename,filename2
        logical :: exists
        character(len=*), parameter :: t_r='SC_Mom_LR_Z'
        
        !How many frequency points are there exactly?
        Omega = Start_Omega
        nESteps = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            nESteps = nESteps + 1
            Omega = Omega + Omega_Step
        enddo
        
        !Find the k-independent self-consistent self-energy
        allocate(SE(nImp,nImp,nESteps))
        allocate(SE_Old(nImp,nImp,nESteps))
        allocate(G_Mat(nImp,nImp,nESteps))
        allocate(Hybrid(nImp,nImp,nESteps))
        allocate(InvLocalMomGF(nImp,nImp,nESteps))
        allocate(LocalMomGF(nImp,nImp,nESteps))
        allocate(LocalCoupFn(nImp,nImp,nESteps))    !X'
        allocate(GlobalCoup(nImp,nImp,nESteps))     !Z

        Hybrid(:,:,:) = zzero
        InvLocalMomGF(:,:,:) = zzero
        LocalMomGF(:,:,:) = zzero
        LocalCoupFn(:,:,:) = zzero
        GlobalCoup(:,:,:) = zzero

        allocate(tOmegaConv(nESteps))
        
        call Init_SE(nESteps,SE)

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp

        iter = 0

        do while(.true.)

            iter = iter + 1

            write(6,*) "Starting iteration: ",iter
            !Construct FT of k-space GFs and take zeroth part.
            !write(6,*) "FT'ing mean-field greens functions for local part"
            call FindLocalMomGF(nESteps,SE,LocalMomGF)
            !call FindRealSpaceLocalMomGF(nESteps,SE,LocalMomGF)
            write(6,*) "1"
            call flush(6)
            call writedynamicfunction(nESteps,LocalMomGF,'LocalMomGF',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)

            !Invert the matrix of non-interacting local greens functions.
            !write(6,*) "Inverting Local greens function"
            InvLocalMomGF(:,:,:) = LocalMomGF(:,:,:)
            write(6,*) "2"
            call flush(6)
            call InvertLocalNonHermGF(nESteps,InvLocalMomGF)

            !Now find hybridization
            !This is given by (omega + mu + idelta - e_0 - SE - InvGF)
            !write(6,*) "Calculating Hybridization function to local greens function"
            write(6,*) "3"
            call flush(6)
            call FindHybrid(nESteps,InvLocalMomGF,SE,Hybrid)
            call writedynamicfunction(nESteps,Hybrid,'Hybrid',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)

            !Check that G^0 = G^0' now.
            write(6,*) "4"
            call flush(6)
            call CheckNIGFsSame(nESteps,LocalMomGF,SE,Hybrid,GFChemPot)

            !Now calculate X', the local coupling function (LocalCoupFn), as
            ![omega + mu + i delta - h00 - Delta]^-1
            write(6,*) "5"
            call flush(6)
            call CalcLocalCoupling(nESteps,Hybrid,LocalCoupFn,GFChemPot)
            call writedynamicfunction(nESteps,LocalCoupFn,'LocalCouplingFn',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)

            !Iteratively converge the global, k-independent coupling quantity 'GlobalCoup' (Z), 
            !which mimics the effect of the hybridization on the whole lattice.
            write(6,*) "6"
            call flush(6)
            tOmegaConv(:) = .false.
            call ConvergeGlobalCoupling(nESteps,LocalCoupFn,GlobalCoup,GFChemPot,1.0e-9_dp,1.0e-10_dp,tOmegaConv,tSuccess,150)

            !Check here that the two coupling functions are identical
            write(6,*) "7"
            call flush(6)
            !call CheckNICoupFnsSame(nESteps,LocalCoupFn,GlobalCoup,GFChemPot)
            call writedynamicfunction(nESteps,GlobalCoup,'GlobalCoupling_Z',tag=iter,   &
                tCheckCausal=.true.,tCheckOffDiagHerm=.true.,tWarn=.false.)

            !Construct interacting greens function from the global coupling
            !calculate G_00
            write(6,"(A)") "Calculating high-level correlation function..."
            write(6,*) "8"
            call flush(6)
            !TEST!
            !GlobalCoup(:,:,:) = zzero
            call SchmidtGF_wSE(G_Mat,GFChemPot,GlobalCoup,nESteps)
            write(6,"(A)") "High-level correlation function obtained."

            !Now calculate the change in the self energy from a damped application
            !of the dyson equation. This is not iterative, and the non-interating
            !and iteractive greens function do not necessarily agree from the outset.
            !This takes the old SE, and outputs the new one.
            SE_Old(:,:,:) = SE(:,:,:)
            call Calc_SE(nESteps,Hybrid,G_Mat,GFChemPot,SE,MaxDiffSE)
            call writedynamicfunction(nESteps,SE,'SelfEnergy',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.,tWarn=.false.)

            !Finally, should we do this all in a larger self-consistency, 
            !such that the self energy is used for the frequency dependent bath?
            if(MaxDiffSE.lt.1.0e-4_dp) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",1.0e-4_dp
                exit
            endif

            !if(iter.eq.1) call stop_all(t_r,'SELF ENERGY CAUSAL (in first iteration)! YAY!')

        enddo

        call WriteConvSE(nESteps,SE,tMatbrAxis=.false.)

        deallocate(LocalCoupFn,SE,SE_Old,G_Mat,Hybrid,InvLocalMomGF,LocalMomGF,GlobalCoup,tOmegaConv)

    end subroutine SC_Mom_LR_Z

    !Attempt to get k-space spectral functions by self-consistenly calculating k-independent hybridization and self-energy contributions
    subroutine SC_Mom_LR()
        implicit none
        real(dp) :: Omega,GFChemPot,OmegaVal,MaxDiffSE,DiffSE,reSE,imSE,r(3)
        complex(dp) :: DiffMatSE(nImp,nImp)
        integer :: nESteps,iter,i,k,l,iunit,j
        complex(dp), allocatable :: SE(:,:,:),G_Mat(:,:,:),SE_Old(:,:,:)
        character(64) :: filename,filename2
        logical :: exists
        character(len=*), parameter :: t_r='SC_Mom_LR'
        
        !How many frequency points are there exactly?
        Omega = Start_Omega
        nESteps = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            nESteps = nESteps + 1
            Omega = Omega + Omega_Step
        enddo
        
        !Find the k-independent self-consistent self-energy
        allocate(SE(nImp,nImp,nESteps))
        allocate(SE_Old(nImp,nImp,nESteps))
        allocate(G_Mat(nImp,nImp,nESteps))

        call Init_SE(nESteps,SE,tMatbrAxis=.false.)

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp

        iter = 0

        do while(.true.)

            iter = iter + 1

            !calculate G_00
!            call NonIntExCont_TDA_MCLR_Charged_Cmprs()
!            call NonIntExCont_TDA_MCLR_Charged()
            write(6,"(A)") "Calculating high-level correlation function..."
            call SchmidtGF_wSE(G_Mat,GFChemPot,SE,nESteps)
            write(6,"(A)") "High-level correlation function obtained."
            write(6,"(A,I8,A)") "Now attempting self-consistent determination of self-energy function from ", &
                nESteps," frequency points."
            call flush(6)

            if((iter.eq.1).and.tRandom_Init_SE) then
                !Randomize the initial SE for the microiterations on the first iteration
                CALL random_seed()
                do i = 1,nESteps
                    do j = 1,nImp
                        do k = 1,nImp
                            CALL Random_Number(r(:))
                            if(r(1).gt.0.5_dp) then
                                SE(k,j,i) = dcmplx(r(2),-r(3))
                            else
                                SE(k,j,i) = dcmplx(-r(2),-r(3))
                            endif
                        enddo
                    enddo
                enddo
            endif
            !tester
            !call testNIGFs(SE,nESteps)

            !Now calculate the (hybridization and) self-energy self-consistently
            !This will read back in the greens function
            !The returned self-energy is k-independent, but will reproduce the correlated local greens function
            SE_Old(:,:,:) = SE(:,:,:)
            call Converge_SE(SE,nESteps)
            !call Converge_SE_NoHybrid(G_Mat,SE,nESteps,iter)

            !What is the overall change in self-energy
            MaxDiffSE = zero
            do i = 1,nESteps
                DiffMatSE(:,:) = SE(:,:,i) - SE_Old(:,:,i)
                DiffSE = zero
                do j = 1,nImp
                    do k = 1,nImp
                        DiffSE = DiffSE + real(dconjg(DiffMatSE(k,j))*DiffMatSE(k,j),dp)
                    enddo
                enddo
                if(DiffSE.gt.MaxDiffSE) MaxDiffSE = DiffSE
            enddo

            !Finally, should we do this all in a larger self-consistency, 
            !such that the self energy is used for the frequency dependent bath?
            if(MaxDiffSE.lt.1.0e-4_dp) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",1.0e-4_dp
                exit
            endif

        enddo

        call WriteConvSE(nESteps,SE,tMatbrAxis=.false.)

    end subroutine SC_Mom_LR
    
    !Initialise the self energy (to zero generally, but can read in one)
    subroutine Init_SE(n,SE,tMatbrAxis)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(out) :: SE(nImp,nImp,n)
        logical, intent(in), optional :: tMatbrAxis
        real(dp) :: OmegaVal,Omega,reSE,imSE
        integer :: i,iunit,k,l
        logical :: exists,tMatbrAxis_
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='Init_SE'

        if(present(tMatbrAxis)) then
            tMatbrAxis_ = tMatbrAxis
        else
            tMatbrAxis_ = .false.
        endif

        SE(:,:,:) = zzero
    
        if(tRead_SelfEnergy) then
            write(6,"(A)") "Reading in self-energy..."
            if(tMatbrAxis_) then
                call append_ext_real('Converged_SE_Im',U,filename)
            else
                call append_ext_real('Converged_SE',U,filename)
            endif
            if(.not.tHalfFill) then
                call append_ext(filename,nOcc,filename2)
            else
                filename2 = filename
            endif
            inquire(file=filename2,exist=exists)
            if(.not.exists) then
                write(6,*) "Converged self-energy does not exist, filename: ",filename2
                call stop_all(t_r,'Cannot find file with converged selfenergy')
            endif
            
            iunit = get_free_unit()
            open(unit=iunit,file=filename2,status='unknown')
            i = 0
            do while(.true.)
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_)
                if(i.lt.0) exit
                !Read in *lower* triangle
                read(iunit,"(G25.10)",advance='no') OmegaVal
                if(abs(OmegaVal-Omega).gt.1.0e-8_dp) then
                    write(6,*) OmegaVal,Omega
                    call stop_all(t_r,'Omega values do not match up for read-in self-energy')
                endif
                do k = 1,nImp
                    do l = k,nImp
                        if((k.eq.nImp).and.(l.eq.nImp)) then
                            exit
                        endif
                        read(iunit,"(2G25.10)",advance='no') reSE,imSE
                        SE(l,k,i) = dcmplx(reSE,imSE)
                    enddo
                enddo
                read(iunit,"(2G25.10)") reSE,imSE
                SE(nImp,nImp,i) = dcmplx(reSE,imSE)
            enddo
            close(iunit)
        endif

    end subroutine Init_SE

    !Routine to read in a correlated greens function, and self-consistently converge
    !a self energy from the 1-electron hemailtonian to get it to match this.
    !MacroIter is the iteration of the main convergence
    subroutine Converge_SE_NoHybrid(G00,SE,nESteps,MacroIter)
        implicit none
        integer, intent(in) :: nESteps,MacroIter
        complex(dp), intent(in) :: G00(nImp,nImp,nESteps)
        complex(dp), intent(inout) :: SE(nImp,nImp,nESteps)
        integer :: iunit,i,iter,iunit2,iunit3,CurrSlot,MaxSlot,j,k
        logical :: exists
        complex(dp), allocatable :: OldSE(:,:,:)
        complex(dp), allocatable :: InvLocalMomGF(:,:,:),LocalMomGF(:,:,:)
        complex(dp), allocatable :: InvG00(:,:,:),RealSpaceLocGF(:,:,:)
        real(dp) :: Omega,Re_LR,Im_LR,Omega_Val,SE_Thresh,MaxDiffSE,MinDiffSE,MeanDiffSE
        real(dp) :: MaxDiffGF,MeanDiffGF,SpecWeight(nImp,nImp),IsoSpecWeight,MeanMax,MeanMean
        real(dp), allocatable :: MaxSEDiffs(:),MeanSEDiffs(:)
        character(64) :: filename,filename2,filename3,filename4,filename5,filename6
        character(128) :: header
        character(len=*), parameter :: t_r='Converge_SE_NoHybrid'

        write(6,*) "Entering Converge_SE_NoHybrid"

        MaxSlot = 10
        SE_Thresh = 1.0e-4_dp

!        if(nImp.ne.1) call stop_all(t_r,'Can currently only do self-consistency with one impurity site')

        allocate(OldSE(nImp,nImp,nESteps))
        allocate(LocalMomGF(nImp,nImp,nESteps)) !The local 1-electron GF from the non-interacting H
        allocate(InvLocalMomGF(nImp,nImp,nESteps))
        allocate(InvG00(nImp,nImp,nESteps))

        call CheckCausality(G00,nImp,nESteps,1)
        call CalcSumRule(G00,nImp,nESteps,SpecWeight,IsoSpecWeight)
        call writematrix(SpecWeight,'Integrated spectral weights for all HL GFs',.true.)
        write(6,*) "Integrated weight for isotropic high-level spectral function: ",IsoSpecWeight
        call CalcSumRule(SE,nImp,nESteps,SpecWeight,IsoSpecWeight)
        call writematrix(SpecWeight,'Integrated spectral weights for all initial Self-energies',.true.)
        write(6,*) "Integrated weight for isotropic initial self-energy: ",IsoSpecWeight
        call flush(6)

        LocalMomGF(:,:,:) = zzero
        OldSE(:,:,:) = SE(:,:,:)
            
        !write(6,*) "Inverting high-level greens function"
        InvG00(:,:,:) = G00(:,:,:)
        call InvertLocalNonHermGF(nESteps,InvG00)

        write(6,"(A)") "      Iter   Max[delta(SE)]         Min[delta(SE)]      Mean[delta(SE)]    Max[delta(GF)]   Mean[delta(GF)]"
        call flush(6)
        iter = 0
        CurrSlot = 0
        allocate(MaxSEDiffs(MaxSlot))
        allocate(MeanSEDiffs(MaxSlot))
        MaxSEDiffs(:) = zero
        MeanSEDiffs(:) = zero
        do while(.true.)
            iter = iter + 1
            write(6,*) "Beginning iteration: ",iter

            !Construct FT of k-space GFs and take zeroth part.
            !write(6,*) "FT'ing mean-field greens functions for local part"
            call FindLocalMomGF(nESteps,SE,LocalMomGF)
            call CheckCausality(LocalMomGF,nImp,nESteps,2)

            if(NIGF_WeightFac.ne.zero) then
                write(6,*) "Increasing weight of NI GFs by: ",NIGF_WeightFac
                LocalMomGF(:,:,:) = LocalMomGF(:,:,:)*NIGF_WeightFac
            endif

            !call FindRealSpaceLocalMomGF(nESteps,SE,LocalMomGF)
            !We now have the updated local greens function from the 1-electron hamiltonian
            !call writevectorcomp(LocalMomGF(1,1,:),'local FT of NI GF')

            !Invert the matrix of non-interacting local greens functions.
            !write(6,*) "Inverting Local greens function"
            InvLocalMomGF(:,:,:) = LocalMomGF(:,:,:)
            call InvertLocalNonHermGF(nESteps,InvLocalMomGF)
            !call writevectorcomp(InvLocalMomGF(1,1,:),'inverse of local FT of NI GF')
            
            !Save previous Self-energy
            OldSE(:,:,:) = SE(:,:,:)

            SE(:,:,:) = SE(:,:,:) + Damping_SE*(InvLocalMomGF(:,:,:) - InvG00(:,:,:))

!            do i = 1,nESteps
!                if((aimag(SE(1,1,i)).gt.zero)) then
!                    SE(1,1,i) = dcmplx(real(SE(1,1,i),dp),zero)
!                endif
!            enddo

            !DEBUG
            open(unit=69,file='SE',status='unknown')
            Omega = Start_Omega
            i = 0
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                if(i.gt.nESteps) call stop_all(t_r,'Too many frequency points')
                write(69,"(5F18.10)") Omega,real(SE(1,1,i),dp),aimag(SE(1,1,i)),real(LocalMomGF(1,1,i),dp),aimag(LocalMomGF(1,1,i))
                Omega = Omega + Omega_Step
            enddo
            close(69)
            call CheckCausality(SE,nImp,nESteps,3)

            !Find maximum difference between self-energies
            call FindSEDiffs(SE,OldSE,G00,LocalMomGF,nESteps,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF)
            write(6,"(I8,5F20.10)") Iter,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
            call flush(6)

            !Calculate mean MAD over last 20 iterations
            if(CurrSlot.ge.MaxSlot) then
                !We have filled all slots
                !Move all points down one
                MeanMax = 0.0_dp
                MeanMean = 0.0_dp
                do i = 2,MaxSlot
                    MaxSEDiffs(i-1) = MaxSEDiffs(i)
                    MeanSEDiffs(i-1) = MeanSEDiffs(i)
                    MeanMax = MeanMax + MaxSEDiffs(i-1)
                    MeanMean = MeanMean + MeanSEDiffs(i-1)
                enddo
                MaxSEDiffs(MaxSlot) = MaxDiffSE
                MeanSEDiffs(MaxSlot) = MeanDiffSE
                MeanMax = (MeanMax + MaxDiffSE) / real(MaxSlot,dp)
                MeanMean = (MeanMean + MeanDiffSE) / real(MaxSlot,dp)
                write(6,"(A,I7,A,2F20.10)") "Running average of maximum/mean SE change over last ", &
                    MaxSlot," iterations: ",MeanMax,MeanMean
            else
                !We haven't done MaxSlot iterations yet
                !Just put into next slot
                MaxSEDiffs(CurrSlot+1) = MaxDiffSE
                MeanSEDiffs(CurrSlot+1) = MeanDiffSE
                CurrSlot = CurrSlot + 1
            endif

            if(MaxDiffSE.lt.(Damping_SE*SE_Thresh)) then
                write(6,"(A,G12.5)") "Success! Convergence to maximum self-energy difference of ",SE_Thresh*Damping_SE
                write(6,"(A,G12.5)") "Mean squared difference between greens functions over all frequencies: ",MeanDiffGF
                exit
            endif

            if((iter.ge.max_SE_iter).and.(max_SE_iter.ne.0)) then
                write(6,*) "Hit max iters for self-energy fit. Exiting..."
                exit
            endif

        enddo
        
!        SE(:,:,:) = SE(:,:,:) * Damping_SE
        call FindLocalMomGF(nESteps,SE,LocalMomGF)
        call CalcSumRule(LocalMomGF,nImp,nESteps,SpecWeight,IsoSpecWeight)
        call writematrix(SpecWeight,'Integrated spectral weights for all converged mean-field GFs',.true.)
        write(6,*) "Integrated weight for isotropic converged mean-field spectral function: ",IsoSpecWeight
        call CalcSumRule(SE,nImp,nESteps,SpecWeight,IsoSpecWeight)
        call writematrix(SpecWeight,'Integrated spectral weights for all converged Self-energies',.true.)
        write(6,*) "Integrated weight for isotropic final self-energy: ",IsoSpecWeight
        call flush(6)

        write(6,"(A)") "Writing out local mean-field greens function with converged self-energy, and Self Energy"
        iunit = get_free_unit()
        call append_ext_real('MF_GF_wSE',U,filename)
        !Also append iteration number 
        call append_ext(filename,MacroIter,filename2)
        open(unit=iunit,file=filename2,status='unknown')

        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            write(iunit,"(4G25.10)") Omega,aimag(LocalMomGF(1,1,i)),aimag(G00(1,1,i)),aimag(SE(1,1,i))
            Omega = Omega + Omega_Step
        enddo
        close(iunit)
        
        if(tCheck) then
            !Use the self-energy striped through the space to find the NI greens function for each frequency.
            !This *I presume* is EXACTLY the same as what I am doing when I calculate the local NI GF in k-space and FT for the local part.
            !This is just more expensive!
            allocate(RealSpaceLocGF(nImp,nImp,nESteps))
            RealSpaceLocGF(:,:,:) = zzero
            call FindRealSpaceLocalMomGF(nESteps,SE,RealSpaceLocGF)

            write(6,"(A)") "Calculating and Writing out converged real-space mean-field greens function"
            call append_ext_real('RealSpace_NIGF',U,filename)
            if(.not.tHalfFill) then
                !Also append occupation of lattice to the filename
                call append_ext(filename,nOcc,filename2)
            else
                filename2 = filename
            endif
            open(unit=iunit,file=filename2,status='unknown')
            Omega = Start_Omega
            i = 0
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                write(iunit,"(3G25.10)") Omega,real(RealSpaceLocGF(1,1,i),dp),aimag(RealSpaceLocGF(1,1,i))
                Omega = Omega + Omega_Step
            enddo
            close(iunit)
            deallocate(RealSpaceLocGF)

        endif   !End check
        
        deallocate(OldSE,InvG00,InvLocalMomGF,LocalMomGF)

    end subroutine Converge_SE_NoHybrid
        
    subroutine CalcSumRule(G,nGF,n,SpecWeight,IsoSpecWeight)
        implicit none
        integer, intent(in) :: n,nGF
        complex(dp), intent(in) :: G(nGF,nGF,n)
        real(dp), intent(out) :: SpecWeight(nGF,nGF),IsoSpecWeight
        real(dp) :: Prev_Spec(nGF,nGF),Prev_IsoGF,IsoGF
        integer :: i,j,k
    
        SpecWeight(:,:) = zero
        IsoSpecWeight = zero

        Prev_Spec(:,:) = zero
        Prev_IsoGF = zero
        do i=1,n
            IsoGF = zero
            do j = 1,nGF
                IsoGF = IsoGF + aimag(G(j,j,i))
                do k = 1,nGF
                    if(aimag(G(k,j,i)).gt.zero) then
                        write(6,*) "WARNING!! Local greens function not causal for freq: ",k,j,i
                    endif
                    SpecWeight(k,j) = SpecWeight(k,j) + Omega_Step*(Prev_Spec(k,j) - aimag(G(k,j,i)) )/(2.0_dp*pi)
                    Prev_Spec(k,j) = -aimag(G(k,j,i))
                enddo
            enddo
            IsoGF = IsoGF / real(nGF,dp)
            IsoSpecWeight = IsoSpecWeight + Omega_Step*(Prev_IsoGF-IsoGF)/(2.0_dp*pi)
            Prev_IsoGF = -IsoGF
        enddo

    end subroutine CalcSumRule

    subroutine CheckCausality(G,nGF,n,flag)
        implicit none
        integer, intent(in) :: nGF,n,flag
        complex(dp), intent(in) :: G(nGF,nGF,n)
        integer :: i,j,k
        character(len=*), parameter :: t_r='CheckCausality'

        do i = 1,n
            do j = 1,nGF
                do k = 1,nGF
                    if(aimag(G(k,j,i)).gt.zero) then
                        write(6,*) "CAUSALITY BROKEN FOR FUNCTION: ",flag
                        write(6,*) "Function component: ",j,k
                        write(6,*) "Frequency point: ",i
                        write(6,*) "Function value: ",G(k,j,i)
                        call warning(t_r,'CAUSALITY BROKEN')
                        return
                    endif
                enddo
            enddo
        enddo

    end subroutine CheckCausality

    !Routine to read in a correlated greens function, and self-consistently converge
    !a self energy from the 1-electron hemailtonian to get it to match this.
    subroutine Converge_SE(SE,nESteps)
        implicit none
        integer, intent(in) :: nESteps
        complex(dp), intent(inout) :: SE(nImp,nImp,nESteps)
        integer :: iunit,i,iter,iunit2,iunit3
        logical :: exists
        complex(dp), allocatable :: G00(:,:,:),Hybrid(:,:,:),OldSE(:,:,:)
        complex(dp), allocatable :: InvLocalMomGF(:,:,:),LocalMomGF(:,:,:)
        complex(dp), allocatable :: InvG00(:,:,:),RealSpaceLocGF(:,:,:)
        real(dp) :: Omega,Re_LR,Im_LR,Omega_Val,SE_Thresh,MaxDiffSE,MinDiffSE,MeanDiffSE
        real(dp) :: MaxDiffGF,MeanDiffGF
        character(64) :: filename,filename2,filename3,filename4,filename5,filename6
        character(128) :: header
        character(len=*), parameter :: t_r='Converge_SE'
        
        write(6,*) "Entering Converge_SE"

        SE_Thresh = 1.0e-4_dp

        if(nImp.ne.1) call stop_all(t_r,'Can currently only do self-consistency with one impurity site')

        !First, read back in the G_00 (real and complex)
        iunit = get_free_unit()
        call append_ext_real('EC-TDA_GFResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        inquire(file=filename2,exist=exists)
        if(.not.exists) then
            call stop_all(t_r,'Cannot find local greens function file')
        endif
        open(unit=iunit,file=filename2,status='old',action='read')

        read(iunit,*) header
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            Omega = Omega + Omega_Step
        enddo
        if(i.ne.nESteps) call stop_all(t_r,'Wrong number of frequency points read in')

        allocate(G00(nImp,nImp,nESteps))    !The correlated greens function
        allocate(Hybrid(nImp,nImp,nESteps)) !The hybridization
        allocate(OldSE(nImp,nImp,nESteps))
        allocate(LocalMomGF(nImp,nImp,nESteps)) !The local 1-electron GF from the non-interacting H
        allocate(InvLocalMomGF(nImp,nImp,nESteps))
        allocate(InvG00(nImp,nImp,nESteps))
        G00(:,:,:) = zzero
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.nESteps) call stop_all(t_r,'Too many frequency points')
            read(iunit,*) Omega_val,Re_LR,Im_LR
            !write(6,*) Omega_val,Re_LR,Im_LR
            !HACK - can only read in 1 x 1 impurity matrix
            G00(1,1,i) = dcmplx(Re_LR,-Im_LR)
            Omega = Omega + Omega_Step
        enddo
        close(iunit)

        write(6,*) "High-level greens function sucessfully read back in..."
        call flush(6)

        Hybrid(:,:,:) = zzero
        LocalMomGF(:,:,:) = zzero
        OldSE(:,:,:) = SE(:,:,:)
            
        !write(6,*) "Inverting high-level greens function"
        InvG00(:,:,:) = G00(:,:,:)
        call InvertLocalNonHermGF(nESteps,InvG00)

        write(6,"(A)") "      Iter   Max[delta(SE)]         Min[delta(SE)]      Mean[delta(SE)]    Max[delta(GF)]   Mean[delta(GF)]"
        iter = 0
        do while(.true.)
            iter = iter + 1
            !write(6,*) "Beginning iteration: ",iter

            !Construct FT of k-space GFs and take zeroth part.
            !write(6,*) "FT'ing mean-field greens functions for local part"
            call FindLocalMomGF(nESteps,SE,LocalMomGF)
            !call FindRealSpaceLocalMomGF(nESteps,SE,LocalMomGF)
            !We now have the updated local greens function from the 1-electron hamiltonian
            !call writevectorcomp(LocalMomGF(1,1,:),'local FT of NI GF')

            !Invert the matrix of non-interacting local greens functions.
            !write(6,*) "Inverting Local greens function"
            InvLocalMomGF(:,:,:) = LocalMomGF(:,:,:)
            call InvertLocalNonHermGF(nESteps,InvLocalMomGF)
            !call writevectorcomp(InvLocalMomGF(1,1,:),'inverse of local FT of NI GF')

            !Now find hybridization
            !This is given by (omega + mu + idelta - e_0 - SE - InvGF)
            !write(6,*) "Calculating Hybridization function to local greens function"
            call FindHybrid(nESteps,InvLocalMomGF,SE,Hybrid)

            !Now to converge the k-independent Self-energy
            !Now find Self-energy
            !This is given by (omega + mu + idelta)I - e_0 - Hybrid - InvFullGF)
            !write(6,*) "Obtaining Self-energy by equating local mean-field GF to high level one"
            call FindSE(nESteps,InvG00,Hybrid,SE)

            !Find maximum difference between self-energies
            call FindSEDiffs(SE,OldSE,G00,LocalMomGF,nESteps,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF)
            OldSE(:,:,:) = SE(:,:,:)
            write(6,"(I8,5F20.10)") Iter,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
            call flush(6)
            if(MaxDiffSE.lt.SE_Thresh) then
                write(6,"(A,G12.5)") "Success! Convergence to maximum self-energy difference of ",SE_Thresh
                write(6,"(A,G12.5)") "Mean squared difference between greens functions over all frequencies: ",MeanDiffGF
                exit
            endif

        enddo

        write(6,"(A)") "Writing out local mean-field greens function with converged self-energy, Self Energy and Hybridization"
        iunit = get_free_unit()
        call append_ext_real('MF_GF_wSE',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        iunit2 = get_free_unit()
        call append_ext_real('SelfEnergy',U,filename3)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename3,nOcc,filename4)
        else
            filename4 = filename3
        endif
        open(unit=iunit2,file=filename4,status='unknown')
        iunit3 = get_free_unit()
        call append_ext_real('Hybridization',U,filename5)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename5,nOcc,filename6)
        else
            filename6 = filename5
        endif
        open(unit=iunit3,file=filename6,status='unknown')

        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            write(iunit,"(3G25.10)") Omega,real(LocalMomGF(1,1,i),dp),aimag(LocalMomGF(1,1,i))
            write(iunit2,"(3G25.10)") Omega,real(SE(1,1,i),dp),aimag(SE(1,1,i))
            write(iunit3,"(3G25.10)") Omega,real(Hybrid(1,1,i),dp),aimag(Hybrid(1,1,i))
            Omega = Omega + Omega_Step
        enddo
        close(iunit)
        close(iunit2)
        close(iunit3)
        
        if(tCheck) then
            !Use the self-energy striped through the space to find the NI greens function for each frequency.
            !This *I presume* is EXACTLY the same as what I am doing when I calculate the local NI GF in k-space and FT for the local part.
            !This is just more expensive!
            allocate(RealSpaceLocGF(nImp,nImp,nESteps))
            RealSpaceLocGF(:,:,:) = zzero
            call FindRealSpaceLocalMomGF(nESteps,SE,RealSpaceLocGF)

            write(6,"(A)") "Calculating and Writing out converged real-space mean-field greens function"
            call append_ext_real('RealSpace_NIGF',U,filename)
            if(.not.tHalfFill) then
                !Also append occupation of lattice to the filename
                call append_ext(filename,nOcc,filename2)
            else
                filename2 = filename
            endif
            open(unit=iunit,file=filename2,status='unknown')
            Omega = Start_Omega
            i = 0
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                write(iunit,"(3G25.10)") Omega,real(RealSpaceLocGF(1,1,i),dp),aimag(RealSpaceLocGF(1,1,i))
                Omega = Omega + Omega_Step
            enddo
            close(iunit)
            deallocate(RealSpaceLocGF)

        endif   !End check
        
        deallocate(G00,OldSE,InvG00,Hybrid,InvLocalMomGF,LocalMomGF)

    end subroutine Converge_SE 
            
    !Check that the local and global coupling functions are the same after the NR optimization.
    subroutine CheckNICoupFnsSame(n,LocalCoupFn,Z,mu)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: LocalCoupFn(nImp,nImp,n)
        complex(dp), intent(in) :: Z(nImp,nImp,n)
        real(dp), intent(in) :: mu

        integer :: SS_Period,i,j,k,lwork,info,kPnt,ind_1,ind_2
        real(dp) :: MaxDiff,Diff,Omega
        complex(dp), allocatable :: Fn(:,:,:),k_Ham(:,:),CompHam(:,:)
        complex(dp), allocatable :: RVec(:,:),LVec(:,:),W_Vals(:),Work(:),ztemp(:,:)
        complex(dp), allocatable :: InvMat(:,:),ztemp2(:,:),cWork(:)
        character(len=*), parameter :: t_r='CheckNICoupFnsSame'

        if(.not.tDiag_kspace) call stop_all(t_r,'k-space diagonalizations must be used here')

        allocate(Fn(nImp,nImp,n))
        Fn(:,:,:) = zzero
        SS_Period = nImp

        allocate(k_Ham(SS_Period,SS_Period))
        allocate(CompHam(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                CompHam(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo

        !Data for the diagonalization of the one-electron k-dependent Greens functions
        allocate(RVec(SS_Period,SS_Period))
        allocate(LVec(SS_Period,SS_Period))
        allocate(W_Vals(SS_Period))
        allocate(Work(max(1,2*SS_Period)))
        allocate(ztemp(nSites,SS_Period))
        allocate(InvMat(nImp,nImp))
        allocate(ztemp2(nImp,nImp))

        do kPnt = 1,nKPnts
            !Run through all k-points
            ind_1 = ((kPnt-1)*SS_Period) + 1
            ind_2 = SS_Period*kPnt
                
            !Transform the one-electron hamiltonian into this k-basis
            call ZGEMM('N','N',nSites,SS_Period,nSites,zone,CompHam,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ztemp,nSites)
            call ZGEMM('C','N',SS_Period,SS_Period,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ztemp,nSites,zzero,k_Ham,SS_Period)
            !k_Ham is the complex, one-electron hamiltonian (with correlation potential) for this k-point

            i = 0
            Omega = Start_Omega
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                if(i.gt.n) call stop_all(t_r,'Too many freq points')

                InvMat(:,:) = - k_Ham(:,:) - Z(:,:,i)
                do j = 1,SS_Period
                    InvMat(j,j) = InvMat(j,j) + dcmplx(Omega + mu,dDelta)
                enddo

                !Now, diagonalize this
                !Is *not* hermitian
                allocate(cWork(1))
                lwork = -1
                info = 0
                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(abs(cWork(1)))+1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
                deallocate(cWork)

                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
                call Order_zgeev_vecs(W_Vals,LVec,RVec)
                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
                call Orthonorm_zgeev_vecs(SS_Period,W_Vals,LVec,RVec)
                !Calculate greens function for this k-vector
                InvMat(:,:) = zzero
                do j = 1,SS_Period
                    InvMat(j,j) = zone / W_Vals(j)
                enddo
                !Now rotate this back into the original basis
                call zGEMM('N','N',SS_Period,SS_Period,SS_Period,zone,RVec,SS_Period,InvMat,SS_Period,zzero,ztemp2,SS_Period)
                call zGEMM('N','C',SS_Period,SS_Period,SS_Period,zone,ztemp2,SS_Period,LVec,SS_Period,zzero,InvMat,SS_Period)
                !InvMat is now the non-interacting greens function for this k: TEST THIS!
                !Sum this into the *local* greens function (i.e. a fourier transform of r=0 component)
                Fn(:,:,i) = Fn(:,:,i) + InvMat(:,:)

                Omega = Omega + Omega_Step
            enddo   !Enddo frequency point i
        enddo   !enddo k-point

        !Divide the entire local greens function by the 'volume' of the brillouin zone
        Fn(:,:,:) = Fn(:,:,:) / real(nSites/nImp,dp)

        !Now, subtract the Local coupling function from all frequencies.
        Fn(:,:,:) = Fn(:,:,:) - LocalCoupFn(:,:,:)

        !Test
        MaxDiff = zero
        do i = 1,n
            Diff = zero
            do j = 1,nImp
                do k = 1,nImp
                    Diff = Diff + abs(Fn(k,j,i))
                enddo
            enddo
            Diff = Diff / real(nImp**2,dp)
            if(Diff.gt.MaxDiff) MaxDiff = Diff
        enddo
        if(MaxDiff.gt.1.0e-6_dp) then
            call stop_all(t_r,'Local/Global Coupling functions not the same after NR optimization of Z...')
        else
            write(6,"(A)") "Local (w Hybrid) / Global (w Z) coupling functions the same after NR optimization of Z..."
        endif

        deallocate(Fn,k_Ham,CompHam,RVec,LVec,W_Vals,Work,ztemp,InvMat,ztemp2)

    end subroutine CheckNICoupFnsSame

    !Given the local part of a global non-interacting greens function, a self-energy and hybridization function,
    !check that the local greens function is the same as the local part of the global greens function.
    subroutine CheckNIGFsSame(n,LocalMomGF,SE,Hybrid,mu)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: LocalMomGF(nImp,nImp,n)
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(in) :: Hybrid(nImp,nImp,n)
        real(dp), intent(in) :: mu
        complex(dp), allocatable :: LocalH(:,:),TrueLocGF(:,:,:)
        real(dp) :: Omega
        integer :: i,j,k
        character(len=*), parameter :: t_r='CheckNIGFsSame'

        allocate(LocalH(nImp,nImp))
        allocate(TrueLocGF(nImp,nImp,n))
        TrueLocGF(:,:,:) = zzero
        do j = 1,nImp
            do k = 1,nImp
                LocalH(k,j) = dcmplx(h0v(k,j),zero)
            enddo
        enddo
        
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Too many freq points')

            !Calculate the local greens function for this frequency
            TrueLocGF(:,:,i) = - LocalH(:,:) - SE(:,:,i) - Hybrid(:,:,i)
            do j = 1,nImp
                TrueLocGF(j,j,i) = TrueLocGF(j,j,i) + dcmplx(Omega + mu,dDelta)
            enddo
            Omega = Omega + Omega_Step
        enddo

        call InvertLocalNonHermGF(n,TrueLocGF)

        do i = 1,n
            !Now, check they're the same
            do j = 1,nImp
                do k = 1,nImp
                    if(abs(LocalMomGF(k,j,i)-TrueLocGF(k,j,i)).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Local and Global Non-interacting greens functions not equal')
                    endif
                enddo
            enddo
        enddo

        write(6,"(A)") "Local (w Hybrid) and Global non-interacting (w SE) greens functions identical..."
        deallocate(LocalH,TrueLocGF)

    end subroutine CheckNIGFsSame
            
    subroutine CheckGFsSame(n,G,G_,dTol)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: G(nImp,nImp,n)
        complex(dp), intent(in) :: G_(nImp,nImp,n)
        real(dp) , intent(in) :: dTol
        character(len=*), parameter :: t_r='CheckGFsSame'
        integer :: i,j,k

        do i = 1,n
            do j = 1,nImp
                do k = 1,nImp
                    if(abs(G(k,j,i)-G_(k,j,i)).gt.dTol) then
                        write(6,*) "Frequency point: ",i
                        write(6,*) "Element: ",j,k
                        write(6,*) "Values: ",G(k,j,i),G_(k,j,i)
                        call stop_all(t_r,'Local functions not the same')
                    endif
                enddo
            enddo
        enddo

    end subroutine CheckGFsSame

    !Now calculate X', the local coupling function, as
    ![omega + mu + i delta - h00 - Delta]^-1
    subroutine CalcLocalCoupling(n,Hybrid,LocalCoup,mu)
        integer, intent(in) :: n
        complex(dp), intent(in) :: Hybrid(nImp,nImp,n)
        complex(dp), intent(out) :: LocalCoup(nImp,nImp,n)
        real(dp), intent(in) :: mu

        complex(dp), allocatable :: LocalH(:,:)
        real(dp) :: Omega
        integer :: i,j,k
        character(len=*), parameter :: t_r='CalcLocalCoupling'

        allocate(LocalH(nImp,nImp))
        do j = 1,nImp
            do k = 1,nImp
                LocalH(k,j) = dcmplx(h0v(k,j),zero)
            enddo
        enddo
        
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Too many freq points')
            
            !Calculate the local coupling function for this frequency
            LocalCoup(:,:,i) = - LocalH(:,:) - Hybrid(:,:,i)
            do j = 1,nImp
                LocalCoup(j,j,i) = LocalCoup(j,j,i) + dcmplx(Omega + mu,dDelta)
            enddo
            Omega = Omega + Omega_Step
        enddo

        call InvertLocalNonHermGF(n,LocalCoup)

        deallocate(LocalH)
    end subroutine CalcLocalCoupling

    !Write out the isotropic average of a dynamic function in the impurity space.
    subroutine writedynamicfunction(n,Func,FileRoot,tag,tCheckCausal,tCheckOffDiagHerm,tWarn,tMatbrAxis,ErrMat,FreqPoints)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: Func(nImp,nImp,n)
        character(len=*), intent(in) :: FileRoot
        integer, intent(in), optional :: tag
        logical, intent(in), optional :: tCheckCausal
        logical, intent(in), optional :: tCheckOffDiagHerm
        logical, intent(in), optional :: tWarn  !If true, don't die if non-causal
        logical, intent(in), optional :: tMatbrAxis !If true, then this correlation function is defined on the Im Axis
        complex(dp), intent(in), optional :: ErrMat(nImp,nImp,n)    !Optional errors on the function
        real(dp), intent(in), optional :: FreqPoints(n)

        character(64) :: filename
        logical :: tCheckOffDiagHerm_,tCheckCausal_,tWarn_,tMatbrAxis_
        integer :: iunit,i,j,k
        real(dp) :: Omega,Prev_Spec,SpectralWeight
        complex(dp) :: IsoAv,IsoErr
        character(len=*), parameter :: t_r='writedynamicfunction'

        if(present(tCheckCausal)) then
            tCheckCausal_ = tCheckCausal
        else
            tCheckCausal_ = .false.
        endif

        if(present(tCheckOffDiagHerm)) then
            tCheckOffDiagHerm_ = tCheckOffDiagHerm
        else
            tCheckOffDiagHerm_ = .false.
        endif

        if(present(tWarn)) then
            tWarn_ = tWarn
        else
            tWarn_ = .false.    !i.e. die by default
        endif

        if(present(tMatbrAxis)) then
            tMatbrAxis_ = tMatbrAxis
        else
            tMatbrAxis_ = .false.
        endif

        if(present(tag)) then
            call append_ext(FileRoot,tag,filename)
        else
            filename = FileRoot
        endif
        iunit = get_free_unit()
        open(unit=iunit,file=filename,status='unknown')

        SpectralWeight = 0.0_dp
        Prev_Spec = 0.0_dp
        i=0
        do while(.true.)
            if(present(FreqPoints)) then
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_,nFreqPoints=n,FreqPoints=FreqPoints)
            else
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_)
            endif
            if(i.lt.0) exit
            if(i.gt.n) call stop_all(t_r,'Wrong number of frequency points used')
            !Find isotropic part
            IsoAv = zzero
            IsoErr = zzero
            do j = 1,nImp
                IsoAv = IsoAv + Func(j,j,i)
                if(present(ErrMat)) IsoErr = IsoErr + ErrMat(j,j,i)
                if(tCheckOffDiagHerm_) then
                    do k = 1,nImp
                        if(k.ne.j) then
                            if(abs(Func(k,j,i)-dconjg(Func(j,k,i))).gt.1.0e-8) then
                                write(6,*) "While writing file: ",filename
                                if(present(tag)) write(6,*) "Filename extension: ",tag
                                write(6,*) "Element of function: ",k,j
                                write(6,*) "Frequency point: ",Omega
                                write(6,*) "Values: ",Func(k,j,i),Func(j,k,i)
                                if(tWarn) then
                                    call warning(t_r,'Function no longer off-diagonal hermitian')
                                else
                                    call stop_all(t_r,'Function no longer off-diagonal hermitian')
                                endif
                            endif
                        endif
                    enddo
                endif
            enddo
            IsoAv = IsoAv / real(nImp,dp)
            IsoErr = IsoErr / real(nImp,dp)
            if(tCheckCausal_.and.(aimag(IsoAv).gt.1.0e-8_dp)) then
                write(6,*) "While writing file: ",filename
                if(present(tag)) write(6,*) "Filename extension: ",tag
                write(6,*) "Frequency point: ",Omega
                write(6,*) "Isotropic value: ",IsoAv,Func(1,1,i)
                if(tWarn_) then
                    call warning(t_r,'Function not causal')
                else
                    call stop_all(t_r,'Function not causal')
                endif
            endif

            if(present(ErrMat)) then
                write(iunit,"(5G25.10)") Omega,real(IsoAv,dp),aimag(IsoAv),real(IsoErr,dp),aimag(IsoErr)
            else
                if(i.ne.1) then
                    if(tMatbrAxis_) then
                        SpectralWeight = SpectralWeight + Omega_Step_Im*(Prev_Spec-aimag(IsoAv))/(2.0_dp*pi)
                    else
                        SpectralWeight = SpectralWeight + Omega_Step*(Prev_Spec-aimag(IsoAv))/(2.0_dp*pi)
                    endif
                    Prev_Spec = -aimag(IsoAv)
                endif
                write(iunit,"(3G25.10)") Omega,real(IsoAv,dp),aimag(IsoAv)
            endif
        enddo
        if(tMatbrAxis_) then
            write(6,"(A,F17.10)") "Total approximate spectral weight on Matsubara axis: ",SpectralWeight
        else
            write(6,"(A,F17.10)") "Total approximate spectral weight on real axis: ",SpectralWeight
        endif
        close(iunit)

    end subroutine writedynamicfunction

    !A simplex algorithm for a direct minimization ofr the cost function
    subroutine ConvergeGlobalCoupling_Direct(n,G,SE,mu,dFuncTol,dChangeSE_Tol,tOmegaConv,tSuccess)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: G(nImp,nImp,n)
        complex(dp), intent(inout) :: SE(nImp,nImp,n)
        real(dp), intent(in) :: mu,dFuncTol,dChangeSE_Tol
        logical, intent(inout) :: tOmegaConv(n)
        logical, intent(out) :: tSuccess

        integer :: i,nFreqVal
        character(len=*), parameter :: t_r='ConvergeGlobalCoupling_Direct'

        tSuccess = .false.

        nFreqVal = 0
        do i = 1,n
            if(.not.tOmegaConv(i)) nFreqVal = nFreqVal + 1
        enddo

        if(nFreqVal.eq.0) call stop_all(t_r,'All frequency values seem converged. Should not be in simplex routine?')
        write(6,"(I7,A,I7,A)") nFreqVal," frequency points out of ",n," failed to converge with NR. "
        write(6,"(A)") "Attempting to converge them with simplex minimization of cost function."

    end subroutine ConvergeGlobalCoupling_Direct

    !Via a NR algorithm, iteratively converge the global coupling (Z) of the cluster to the lattice
    !LocalCoupFunc is X'(omega)
    subroutine ConvergeGlobalCoupling(n,LocalCoupFunc,GlobalCoup,mu,dFuncTol,dConvTol,tOmegaConv,tSuccess,iMaxIter)
        use matrixops, only: mat_inv
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: LocalCoupFunc(nImp,nImp,n)
        complex(dp), intent(inout) :: GlobalCoup(nImp,nImp,n)
        real(dp), intent(in) :: mu,dConvTol,dFuncTol
        logical, intent(out) :: tSuccess,tOmegaConv(n)
        integer, intent(in) :: iMaxIter

        real(dp) :: Omega
        real(dp) :: AbsChange,MaxAbsChange
        complex(dp), allocatable :: FuncVal(:,:,:),Grad(:,:,:)
        complex(dp), allocatable :: ChangeZ(:,:),k_Hams(:,:,:)
        complex(dp), allocatable :: Fn(:,:,:),InvMat(:,:)
        integer :: i,iter,j,k,kPnt
        character(len=*), parameter :: t_r='ConvergeGlobalCoupling'
        
!        dDampedNR = 0.2_dp
!        dConvTol = 1.0e-10_dp
        write(6,*) "Converging global coupling matrix..."
        tOmegaConv(:) = .false.

        allocate(FuncVal(nImp,nImp,n))
        allocate(Grad(nImp,nImp,n))
        allocate(ChangeZ(nImp,nImp))

        allocate(k_Hams(nImp,nImp,nKPnts))
        call CreateKHamBlocks(k_Hams)

        iter = 0
        do while(.true.)
            iter = iter + 1
            write(6,*) "Starting NR iteration: ",iter

            !Get function values and Jacobian (diagonal) for all omega values
            call GetFuncValsandGrads(n,mu,LocalCoupFunc,GlobalCoup,k_Hams,FuncVal,Grad)

            MaxAbsChange = zero
            Omega = Start_Omega
            tSuccess = .true.
            do i = 1,n
                !Run over omega values

                call FindChangeZ(Omega,mu,LocalCoupFunc(:,:,i),GlobalCoup(:,:,i),k_Hams,Grad(:,:,i),    &
                    FuncVal(:,:,i),ChangeZ,AbsChange,tOmegaConv(i),dFuncTol,dConvTol)

                if(.not.tOmegaConv(i)) then
                    tSuccess = .false.  !We have not converged all freq
                    if(AbsChange.gt.MaxAbsChange) MaxAbsChange = AbsChange

                    !Update the global coupling parameter
                    !GlobalCoup(:,:,i) = GlobalCoup(:,:,i) + (dDampedNR*ChangeZ(:,:))
                    GlobalCoup(:,:,i) = GlobalCoup(:,:,i) + ChangeZ(:,:)
                endif

                Omega = Omega + Omega_Step
            enddo

            write(6,*) "MaxAbsChange: ",MaxAbsChange
            call writedynamicfunction(n,GlobalCoup,'Conv_SE')
        
            !TEST! calulate and write out G0, to compare to G_imp
            allocate(InvMat(nImp,nImp))
            allocate(Fn(nImp,nImp,n))
            do kPnt = 1,nKPnts
                i = 0
                Omega = Start_Omega
                do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                    i = i + 1

                    ChangeZ(:,:) = - k_Hams(:,:,kPnt) - GlobalCoup(:,:,i)
                    do j = 1,nImp
                        ChangeZ(j,j) = ChangeZ(j,j) + dcmplx(Omega + mu,dDelta)
                    enddo
                    call mat_inv(ChangeZ,InvMat)
                    !Sum this into the *local* greens function (i.e. a fourier transform of r=0 component)
                    Fn(:,:,i) = Fn(:,:,i) + InvMat(:,:)
                    Omega = Omega + Omega_Step
                enddo   !Enddo frequency point i
            enddo   !enddo k-point
            !Divide the entire local greens function by the 'volume' of the brillouin zone
            Fn(:,:,:) = Fn(:,:,:) / real(nSites/nImp,dp)
            call writedynamicfunction(n,Fn,'G0_NR')
            call writedynamicfunction(n,FuncVal,'FuncVal')
            call writedynamicfunction(n,Grad,'Grad')
            deallocate(InvMat)
            deallocate(Fn)
            !END TEST

!            if(MaxAbsChange.lt.dConvTol) then
            if(tSuccess) then
                write(6,"(A,G14.7)") "Convergence of global coupling function achieved to accuracy of: ",dConvTol
                write(6,"(A,I9)") "Number of microiterations required: ",iter
                exit
            elseif(iter.gt.iMaxIter) then
                write(6,"(A,I8,A)") "Newton-Raphson root finding failed after ",iMaxIter," iterations." 
                write(6,"(A)") "Turning to simplex minimization for remaining frequency points"
                exit
            endif
        enddo
        deallocate(FuncVal,Grad,ChangeZ,k_Hams)

    end subroutine ConvergeGlobalCoupling
                
    subroutine FindChangeZ(Omega,mu,X_prime,Z,k_Hams,Grad,Func,ChangeZ,AbsChange,tConv,dFuncTol,dChangeTol)
        use matrixops, only: mat_inv
        implicit none
        real(dp), intent(in) :: Omega,mu
        complex(dp), intent(in) :: X_prime(nImp,nImp)
        complex(dp), intent(in) :: Z(nImp,nImp)
        complex(dp), intent(in) :: k_Hams(nImp,nImp,nKPnts)
        complex(dp), intent(in) :: Grad(nImp,nImp),Func(nImp,nImp)
        complex(dp), intent(out) :: ChangeZ(nImp,nImp)
        real(dp), intent(out) :: AbsChange
        logical, intent(inout) :: tConv !Is this frequency point converged?
        real(dp), intent(in) :: dFuncTol,dChangeTol !Tolerance for function magnitude, and change in function

        integer :: OptType,i,j,kPnt
        real(dp) :: OptDamp,CurrOptRes,r(3)
        complex(dp), allocatable :: NewFn(:,:),Fn(:,:),NewZ(:,:),InvMat(:,:)

        if(tConv) then 
            AbsChange = zero
            return    !We don't have anything to do here
        endif

        ChangeZ(:,:) = zzero
        AbsChange = zero
        allocate(NewFn(nImp,nImp))
        allocate(Fn(nImp,nImp))
        allocate(NewZ(nImp,nImp))
        allocate(InvMat(nImp,nImp))

        OptType = 0
        do i = 1,nImp
            do j = 1,nImp
                if(abs(Grad(j,i)).lt.(4.0_dp*eps)) then
                    !Gradient numerically zero.
                    !Just move a little bit in an arbitrary direction.
                    OptType = 1
                elseif(abs(Grad(j,i)).lt.1.0e-6_dp) then
                    OptType = 2
                endif
            enddo
        enddo
        
        if(OptType.eq.0) then
            !If the gradient is not small, then simply use an undamped newton raphson step
            ChangeZ(:,:) = - Func(:,:) / Grad(:,:)
            AbsChange = sum(abs(ChangeZ(:,:))) 
            if((sum(abs(Func)).lt.dFuncTol).and.(AbsChange.lt.dChangeTol)) then
                !Changes small, and function small.

                !However, also check that the value of Z (self-energy) is actually causal
                do i = 1,nImp
                    if(aimag(Z(i,i)).gt.zero) then
                        !not causal
                        tConv = .false.
                        exit
                    else
                        tConv = .true.
                    endif
                enddo
            endif
        elseif(OptType.eq.2) then
            !Do a linesearch along the gradient direction.
            CurrOptRes = sum(abs(Func(:,:)))
            OptDamp = 0.0_dp
            do i = 1,12

                NewZ(:,:) = Z(:,:) - dcmplx(0.1_dp*real(i*dp),zero) * Func(:,:) / Grad(:,:)

                !Calculate new function
                NewFn(:,:) = zzero
                do kPnt = 1,nKPnts
                    Fn(:,:) = - k_Hams(:,:,kPnt) - NewZ(:,:)
                    do j = 1,nImp
                        Fn(j,j) = Fn(j,j) + dcmplx(Omega + mu,dDelta)
                    enddo
                    call mat_inv(Fn,InvMat)
                    NewFn(:,:) = NewFn(:,:) + InvMat(:,:)
                enddo
                NewFn(:,:) = NewFn(:,:) / real(nKPnts,dp)
                NewFn(:,:) = NewFn(:,:) - X_prime(:,:)

                if((sum(abs(NewFn(:,:)))).lt.CurrOptRes) then
                    !This is a better value of Z
                    OptDamp = 0.1_dp*real(i*dp)
                    CurrOptRes = sum(abs(NewFn(:,:)))
                endif
            enddo

            if(OptDamp.gt.0.05_dp) then
                !We have moved, and the result is better! Use this
                ChangeZ(:,:) = - dcmplx(OptDamp,zero) * Func(:,:) / Grad(:,:)
                AbsChange = sum(abs(ChangeZ(:,:)))
            else
                !We have not moved.
                !Is where we are actually a root?
                if(CurrOptRes.lt.1.0e-7_dp) then
                    !Fine. Don't move then. This just unhelpfully happens to also be a stationary point
                    ChangeZ(:,:) = zzero
                    AbsChange = zero
                else
                    !We have a problem. Gradient is small, and every direction along the gradient direction seems to increase the value of the function.
                    !Change to a random value of Z, and try again next iteration
                    !I guess we could choose a finer mesh?
                    write(6,*) "Choosing random new Z as at stationary point...",Z(:,:)
                    do i = 1,nImp
                        do j = 1,nImp
                            call random_number(r(:))
                            if(r(1).gt.0.5) then
                                NewZ(j,i) = dcmplx(r(2),-r(3))
                            else
                                NewZ(j,i) = dcmplx(-r(2),-r(3))
                            endif
                        enddo
                    enddo
                    ChangeZ(:,:) = NewZ(:,:) - Z(:,:)
                    AbsChange = sum(abs(ChangeZ(:,:)))
                endif
            endif
        else
            !Gradient is numerically zero. Bad times. Choose random new Z value.
            write(6,*) "Choosing random new Z as at stationary point..."
            do i = 1,nImp
                do j = 1,nImp
                    call random_number(r(:))
                    if(r(1).gt.0.5) then
                        NewZ(j,i) = dcmplx(r(2),-r(3))
                    else
                        NewZ(j,i) = dcmplx(-r(2),-r(3))
                    endif
                enddo
            enddo
            ChangeZ(:,:) = NewZ(:,:) - Z(:,:)
            AbsChange = sum(abs(ChangeZ(:,:)))
        endif

        deallocate(NewFn,Fn,NewZ,InvMat)
    end subroutine FindChangeZ

    !Write out the surface of |F(Sigma)|^2 as a function of sigma.
    !If this doesn't fit zero, there is no solution.
    !This is only done for the frequency Start_Omega
    subroutine PlotSigmaSurface(n,G,mu)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: G(nImp,nImp,n)
        real(dp), intent(in) :: mu
        complex(dp) :: A,Factor,MinSE,MaxSE
        integer :: iunit,kPnt
        real(dp) :: Omega,SE_Real,SE_Im,MaxValue,MinValue
        complex(dp), allocatable :: k_Hams(:,:,:)
        character(len=*), parameter :: t_r='PlotSigmaSurface'

        if(nImp.gt.1) call stop_all(t_r,'Routine only designed for one impurity')
        
        !Open file
        iunit = get_free_unit()
        open(iunit,file='SigmaScan',status='unknown')

        Omega = Start_Omega
        Factor = dcmplx(real(nSites/nImp,dp),zero)
        
        allocate(k_Hams(nImp,nImp,nKPnts))
        call CreateKHamBlocks(k_Hams)

        MaxValue = zero
        MinValue = huge(1.0_dp)
        MinSE = zzero
        MaxSE = zzero

        SE_Real = Start_SE_Real
        do while((SE_Real.lt.max(Start_SE_Real,End_SE_Real)+1.0e-6_dp).and.(SE_Real.gt.min(Start_SE_Real,End_SE_Real)-1.0e-5_dp))
            SE_Im = Start_SE_Im
            do while((SE_Im.lt.max(Start_SE_Im,End_SE_Im)+1.0e-6_dp).and.(SE_Im.gt.min(Start_SE_Im,End_SE_Im)-1.0e-5_dp))
                A = zzero
                do kPnt = 1,nKPnts
                    A = A + zone/(dcmplx(Omega+mu-SE_Real,dDelta-SE_Im) - k_Hams(1,1,kPnt))
                enddo
                A = A / Factor

                A = A - G(1,1,1)
                A = A * conjg(A)

                write(iunit,"(3G20.10)") SE_Real,SE_Im,real(A,dp)
                if(real(A,dp).gt.MaxValue) then
                    MaxValue = real(A,dp)
                    MaxSE = dcmplx(SE_Real,SE_Im)
                endif
                if(real(A,dp).lt.MinValue) then
                    MinValue = real(A,dp)
                    MinSE = dcmplx(SE_Real,SE_Im)
                endif
                SE_Im = SE_Im + SE_Step_Im
            enddo
            write(iunit,"(A)") ""
            SE_Real = SE_Real + SE_Step_Real
        enddo

        write(6,*) "Maximum value of function was: ",MaxValue
        write(6,*) "This was found at a SE of: ",MaxSE
        write(6,*) 
        write(6,*) "Minimum value of function was: ",MinValue
        write(6,*) "This was found at a SE of: ",MinSE

        close(iunit)
        deallocate(k_Hams)

    end subroutine PlotSigmaSurface

!For frequencies that do not converge, or are not causal, simply set their value to be equal to the 
!last causal self-energy
    subroutine fit_noncausal_SE(n,SE,tOmegaConv,tSuccess)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(inout) :: SE(nImp,nImp,n)
        logical, intent(inout) :: tOmegaConv(n)
        logical, intent(in) :: tSuccess

        integer :: i,iZeroFreq,j
        real(dp) :: Omega
        complex(dp) :: dFreqFinalCausal(nImp),SE_LastConv(nImp)
        integer :: iFreqCausalFail(nImp)
        character(len=*), parameter :: t_r='fit_noncausal_SE'


        if(.false.) then
            !This option makes the self-energy the same from the first non-causal/unconverged value.
            if(Omega_Step.lt.zero) call stop_all(t_r,'Currently only working ramping up in frequency')

            iZeroFreq = -1
            iFreqCausalFail(:) = -1
            dFreqFinalCausal(:) = zzero

            Omega = Start_Omega
            i = 0
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1

                if(Omega.ge.zero) then
                    if(iZeroFreq.eq.-1) iZeroFreq = i
                    !Are we non-causal in any of the self-energy components?
                    do j = 1,nImp
                        if(((aimag(SE(j,j,i)).gt.zero).or.(.not.tOmegaConv(i))).and.(iFreqCausalFail(j).eq.-1)) then
                            !This is the first non-causal/converged frequency in the positive region for GF j
                            iFreqCausalFail(j) = i
                            dFreqFinalCausal(j) = SE(j,j,i-1)
                            write(6,"(A,I6,A,F13.7)") "First non-causal positive frequency point for self-energy component ", & 
                                j," at: ",Omega
                        endif
                    enddo
                endif
                Omega = Omega + Omega_Step
            enddo

            if(sum(iFreqCausalFail(:)).eq.-nImp) then
                write(6,*) "Self-energy appears causal for all frequencies!"
                if(.not.tSuccess) then
                    write(6,*) "However, some frequency points are unconverged (though causal)." 
                    write(6,*) "These will be set to be equal to their neighbouring self-energy."
                endif
            endif

            do i = 1,nImp
                if(iFreqCausalFail(i).ne.-1) then
                    write(6,*) "For component ",i," non-causal/converged energy at datapoint ",iFreqCausalFail(i)
                    write(6,*) "Fixing subsequent datapoints to equal the last causal and converged self-energy of: ", &
                        dFreqFinalCausal(i)
                endif
            enddo

            !Now, in the simplest approximation, just set the non-causal frequency points to equal the final causal frequency point.
            do i = 1,nImp
                do j = iFreqCausalFail(i),n
                    SE(i,i,j) = dFreqFinalCausal(i)
                    tOmegaConv(j) = .true.
                enddo
            enddo

            !What about unconverged, but causal frequency points. Set them to be the last converged datapoint.
            if(.not.tOmegaConv(iZeroFreq)) call stop_all(t_r,'Zero frequency value did not converge')
            do i = iZeroFreq,n
                do j = 1,nImp
                    if(.not.tOmegaConv(i)) then 
                        call stop_all(t_r,'We should not get here any more')
                        SE(j,j,i) = SE_LastConv(j)
                    else
                        SE_LastConv(j) = SE(j,j,i)
                    endif
                enddo
            enddo

            iFreqCausalFail(:) = -1
            dFreqFinalCausal(:) = zzero
            !Now, repeat for the negative frequencies
            do i = iZeroFreq,1,-1
                do j = 1,nImp
                    if(((aimag(SE(j,j,i)).gt.zero).or.(.not.tOmegaConv(i))).and.(iFreqCausalFail(j).eq.-1)) then
                        !This is the first non-causal frequency in the negative region for GF j
                        iFreqCausalFail(j) = i
                        dFreqFinalCausal(j) = SE(j,j,i+1)
                        write(6,"(A,I6,A,F13.7)") "First non-causal/converged negative frequency point for self-energy component ",j," at: ",i
                    endif
                enddo
            enddo

            if(sum(iFreqCausalFail(:)).eq.-nImp) then
                write(6,*) "Self-energy appears causal for all negative frequencies!"
                if(.not.tSuccess) then
                    write(6,*) "However, some frequency points are unconverged (though causal)." 
                    write(6,*) "These will be set to their nearest neighbouring value"
                endif
            endif

            do i = 1,nImp
                if(iFreqCausalFail(i).ne.-1) then
                    write(6,*) "For component ",i," non-causal energy at datapoint ", iFreqCausalFail(i)
                    write(6,*) "Fixing subsequent datapoints to equal the last causal self-energy of: ",dFreqFinalCausal(i)
                endif
            enddo

            do i = 1, nImp
                do j = iFreqCausalFail(i),1,-1
                    SE(i,i,j) = dFreqFinalCausal(i)
                    tOmegaConv(j) = .true.
                enddo
            enddo

            !What about unconverged, but causal frequency points. Set them to be the last converged datapoint.
            do i = iZeroFreq,1,-1
                do j = 1,nImp
                    if(.not.tOmegaConv(i)) then
                        call stop_all(t_r,'Should not get here any more')
                        SE(j,j,i) = SE_LastConv(j)
                    else
                        SE_LastConv(j) = SE(j,j,i)
                    endif
                enddo
            enddo
        else
            !Here, we instead only make the unconverged/non-causal self-energy points equal to their previous value.

            !Find the midpoint (only for half-filling!)
            Omega = Start_Omega
            i = 0
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1

                if(Omega.ge.zero) then
                    iZeroFreq = i
                    exit
                endif
                Omega = Omega + Omega_Step
            enddo

            do i = 1,nImp
                SE_LastConv(i) = SE(i,i,iZeroFreq)
            enddo

            do i = iZeroFreq,n
                do j = 1,nImp
                    if((.not.tOmegaConv(i)).or.(aimag(SE(j,j,i)).gt.zero)) then
                        SE(j,j,i) = SE_LastConv(j)
                    else
                        SE_LastConv(j) = SE(j,j,i)
                    endif
                enddo
            enddo
            
            do i = iZeroFreq-1,1,-1
                do j = 1,nImp
                    if((.not.tOmegaConv(i)).or.(aimag(SE(j,j,i)).gt.zero)) then
                        SE(j,j,i) = SE_LastConv(j)
                    else
                        SE_LastConv(j) = SE(j,j,i)
                    endif
                enddo
            enddo

        endif

    end subroutine fit_noncausal_SE

    !Get function values and gradients for all values of omega
    !LocalCoupFunc is the X'(omega) function with the local GF - local hybridization
    !GlobalCoup is the striped coupling function through k-space
    subroutine GetFuncValsandGrads(n,mu,LocalCoupFunc,GlobalCoup,k_Hams,FuncVal,Grad)
        use matrixops, only: mat_inv
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: mu
        complex(dp), intent(in) :: LocalCoupFunc(nImp,nImp,n)
        complex(dp), intent(in) :: GlobalCoup(nImp,nImp,n)
        complex(dp), intent(in) :: k_Hams(nImp,nImp,nkPnts)
        complex(dp), intent(out) :: FuncVal(nImp,nImp,n)
        complex(dp), intent(out) :: Grad(nImp,nImp,n)

        real(dp) :: Omega
        complex(dp) :: Factor
        complex(dp), allocatable :: A(:,:),A_Inv(:,:)
        complex(dp), allocatable :: GradTmp(:,:)
        integer :: i,j,ind_1,ind_2,kPnt
        character(len=*), parameter :: t_r='GetFuncValsandGrads'

        if(.not.tDiag_kspace) call stop_all(t_r,'Real space diagonalizations not available here')

        FuncVal(:,:,:) = zzero
        Grad(:,:,:) = zzero
        Factor = dcmplx(real(nSites/nImp,dp),zero)

        allocate(A(nImp,nImp))
        allocate(A_Inv(nImp,nImp))
        allocate(GradTmp(nImp,nImp))

        do kPnt = 1,nKPnts
            !Run through all k-points
            !k_Hams contain the complex, one-electron hamiltonian (with correlation potential) for this k-point
            i = 0
            Omega = Start_Omega
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                if(i.gt.n) call stop_all(t_r,'Too many freq points')

                !2) Construct A_k(omega)
                A(:,:) = zzero
                do j = 1,nImp
                    A(j,j) = dcmplx(Omega+mu,dDelta)
                enddo
                A(:,:) = A(:,:) - k_Hams(:,:,kPnt) - GlobalCoup(:,:,i)
                A(:,:) = A(:,:) * Factor
        
                !3) Invert A_k(omega)
                A_Inv(:,:) = zzero
                call mat_inv(A,A_Inv)

                !4) Sum into function evaluation for this frequency
                FuncVal(:,:,i) = FuncVal(:,:,i) + A_Inv(:,:)

                !5) Sum the square of the inverse into the jacobian construction
                call ZGEMM('N','N',nImp,nImp,nImp,zone,A_Inv,nImp,A_Inv,nImp,zzero,GradTmp,nImp)
                Grad(:,:,i) = Grad(:,:,i) + GradTmp(:,:)

                Omega = Omega + Omega_Step
            enddo
        enddo

        !The Function evaluation needs to have the X' value removed from it at all points
        FuncVal(:,:,:) = FuncVal(:,:,:) - LocalCoupFunc(:,:,:)
        Grad(:,:,:) = Grad(:,:,:)*Factor
        
        deallocate(A,A_Inv,GradTmp)
    end subroutine GetFuncValsandGrads

    !Damped update of the self-energy via the dyson equation
    subroutine Calc_SE(n,Hybrid,G,mu,SE,MaxDiffSE)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: Hybrid(nImp,nImp,n)
        complex(dp), intent(in) :: G(nImp,nImp,n)
        real(dp), intent(in) :: mu
        complex(dp), intent(inout) :: SE(nImp,nImp,n)
        real(dp), intent(out) :: MaxDiffSE

        complex(dp), allocatable :: InvG00(:,:,:),DeltaSE(:,:),LocalH(:,:)
        real(dp) :: Omega,DiffSE
        integer :: i,j,k
        character(len=*), parameter :: t_r='Calc_SE'

        !First, invert the interacting greens function
        allocate(InvG00(nImp,nImp,n))
        InvG00(:,:,:) = G(:,:,:)
        call InvertLocalNonHermGF(n,InvG00)
        
        allocate(LocalH(nImp,nImp))
        do j = 1,nImp
            do k = 1,nImp
                LocalH(k,j) = dcmplx(h0v(k,j),zero)
            enddo
        enddo

        allocate(DeltaSE(nImp,nImp))

        Omega = Start_Omega
        i=0
        MaxDiffSE = zero
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Wrong number of frequency points used')

            DeltaSE(:,:) = zzero
            do j = 1,nImp
                DeltaSE(j,j) = dcmplx(Omega+mu,dDelta)
            enddo
            DeltaSE(:,:) = DeltaSE(:,:) - SE(:,:,i) - LocalH(:,:) - Hybrid(:,:,i) - InvG00(:,:,i)

            SE(:,:,i) = SE(:,:,i) + Damping_SE*DeltaSE(:,:)

            DiffSE = zero
            do j = 1,nImp
                do k = 1,nImp
                    DiffSE = DiffSE + real(dconjg(DeltaSE(k,j))*DeltaSE(k,j),dp)
                enddo
            enddo
            if(DiffSE.gt.MaxDiffSE) MaxDiffSE = DiffSE

            Omega = Omega + Omega_Step
        enddo

        deallocate(DeltaSE,LocalH,InvG00)
    end subroutine Calc_SE

!    !Random test for my own edification. At Sigma=0, can we get the real-space and k-space NI GFs to agree?
!    subroutine TestNIGFs(SE,nESteps)
!        use utils, only: get_free_unit
!        implicit none
!        integer, intent(in) :: nESteps
!        complex(dp), intent(in) :: SE(nImp,nImp,nESteps)
!        complex(dp), allocatable :: LocalMomGF(:,:,:),RealSpaceLocGF(:,:,:)
!        integer :: iunit,i
!        real(dp) :: Omega,intweight,Prev_Spec
!
!!        allocate(LocalMomGF(nImp,nImp,nESteps))
!!        allocate(RealSpaceLocGF(nImp,nImp,nESteps))
!!
!!        call FindLocalMomGF(nESteps,SE,LocalMomGF)
!!        call FindRealSpaceLocalMomGF(nESteps,SE,RealSpaceLocGF)
!!        intweight = zero
!!        Prev_Spec = zero
!!        do i=1,nESteps
!!            intweight = intweight + Omega_Step*(Prev_Spec-aimag(RealSpaceLocGF(1,1,i)))/(2.0_dp*pi)
!!            Prev_Spec = -aimag(RealSpaceLocGF(1,1,i))
!!        enddo
!!        write(6,*) "Integrated weight for r-space NI spectral function: ",intweight
!!
!!        iunit = get_free_unit()
!!        open(unit=iunit,file='temp',status='unknown')
!!
!!        Omega = Start_Omega
!!        i = 0
!!        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
!!            i = i + 1
!!            write(iunit,"(7G25.10)") Omega,real(LocalMomGF(1,1,i),dp),aimag(LocalMomGF(1,1,i)),real(RealSpaceLocGF(1,1,i),dp),  &
!!                    aimag(RealSpaceLocGF(1,1,i)),real(LocalMomGF(1,1,i),dp)/real(RealSpaceLocGF(1,1,i),dp), &
!!                    aimag(localMomGF(1,1,i))/aimag(RealSpaceLocGF(1,1,i))
!!            Omega = Omega + Omega_Step
!!        enddo
!!        close(iunit)
!
!
!!Diagonalize a self-energy in r-space
!        allocate(CompHam(nSites,nSites))
!
!        call add_localpot_comp_inplace(CompHam,SE)
!        allocate(RVec(nSites,nSites))
!        allocate(LVec(nSites,nSites))
!        allocate(W_Vals(nSites))
!        allocate(Work(max(1,2*nSites)))
!        allocate(cWork(1))
!        lwork = -1
!        info = 0
!        call zgeev('V','V',nSites,CompHam,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
!        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
!        lwork = int(abs(cWork(1)))+1
!        deallocate(cWork)
!        allocate(cWork(lWork))
!        call zgeev('V','V',nSites,Mat,nSites,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
!        if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
!        deallocate(cWork)
!
!
!
!
!
!
!        !Data for the diagonalization of the one-electron k-dependent Greens functions
!        allocate(RVec(SS_Period,SS_Period))
!        allocate(LVec(SS_Period,SS_Period))
!        allocate(W_Vals(SS_Period))
!        allocate(Work(max(1,2*SS_Period)))
!        allocate(ztemp(nSites,SS_Period))
!
!        do kPnt = 1,nKPnts
!            !Run through all k-points
!            ind_1 = ((kPnt-1)*SS_Period) + 1
!            ind_2 = SS_Period*kPnt
!                
!            !Transform the one-electron hamiltonian into this k-basis
!            RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
!            call ZGEMM('N','N',nSites,SS_Period,nSites,zone,CompHam,nSites,RotMat,nSites,zzero,ztemp,nSites)
!            call ZGEMM('C','N',SS_Period,SS_Period,nSites,zone,RotMat,nSites,ztemp,nSites,zzero,k_Ham,SS_Period)
!            !k_Ham is the complex, one-electron hamiltonian (with correlation potential) for this k-point
!
!            !write(6,*) "For k-point: ",kPnt
!            !call writematrixcomp(k_Ham,'k-space ham is: ',.false.)
!
!            i = 0
!            Omega = Start_Omega
!            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
!                i = i + 1
!                if(i.gt.n) call stop_all(t_r,'Too many freq points')
!
!                !Find inverse greens function for this k-point
!                InvMat(:,:) = - k_Ham(:,:) - SE(:,:,i)
!                do j = 1,SS_Period
!                    InvMat(j,j) = InvMat(j,j) + dcmplx(Omega + mu,dDelta)
!                enddo
!
!                !Now, diagonalize this
!                !The self-energy is *not* hermitian
!                allocate(cWork(1))
!                lwork = -1
!                info = 0
!                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
!                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
!                lwork = int(abs(cWork(1)))+1
!                deallocate(cWork)
!                allocate(cWork(lWork))
!                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
!                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
!                deallocate(cWork)
!
!                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
!                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
!                call Order_zgeev_vecs(W_Vals,LVec,RVec)
!                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
!                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
!                call Orthonorm_zgeev_vecs(SS_Period,W_Vals,LVec,RVec)
!                !Calculate greens function for this k-vector
!                InvMat(:,:) = zzero
!                do j = 1,SS_Period
!                    InvMat(j,j) = zone / W_Vals(j)
!                enddo
!                !Now rotate this back into the original basis
!                call zGEMM('N','N',SS_Period,SS_Period,SS_Period,zone,RVec,SS_Period,InvMat,SS_Period,zzero,ztemp2,SS_Period)
!                call zGEMM('N','C',SS_Period,SS_Period,SS_Period,zone,ztemp2,SS_Period,LVec,SS_Period,zzero,InvMat,SS_Period)
!                !InvMat is now the non-interacting greens function for this k: TEST THIS!
!                !Sum this into the *local* greens function (i.e. a fourier transform of r=0 component)
!                LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + InvMat(:,:)
!
!                !write(6,*) "For frequency: ",Omega
!                !call writematrixcomp(InvMat,'k-space greens function is: ',.false.)
!
!                Omega = Omega + Omega_Step
!            enddo   !Enddo frequency point i
!
!        enddo   !enddo k-point
!
!        deallocate(LocalMomGF,RealSpaceLocGF)
!    end subroutine TestNIGFs


    subroutine FindSEDiffs(SE,OldSE,G00,LocalMomGF,n,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: SE(nImp,nImp,n),OldSE(nImp,nImp,n),G00(nImp,nImp,n),LocalMomGF(nImp,nImp,n)
        real(dp), intent(out) :: MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
        integer :: i,j,k
        complex(dp) :: DiffMatSE(nImp,nImp),DiffMatGF(nImp,nImp)
        real(dp) :: DiffSE,DiffGF

        MaxDiffSE = zero
        MinDiffSE = 1000000.0_dp
        MeanDiffSE = zero
        MaxDiffGF = zero
        MeanDiffGF = zero
        do i = 1,n
    
            DiffMatSE(:,:) = SE(:,:,i) - OldSE(:,:,i)
            DiffMatGF(:,:) = G00(:,:,i) - LocalMomGF(:,:,i)

            DiffSE = zero
            DiffGF = zero
            do j = 1,nImp
                do k = 1,nImp
                    DiffSE = DiffSE + real(dconjg(DiffMatSE(k,j))*DiffMatSE(k,j),dp)
                    DiffGF = DiffGF + real(dconjg(DiffMatGF(k,j))*DiffMatGF(k,j),dp)
                enddo
            enddo

            if(DiffSE.gt.MaxDiffSE) MaxDiffSE = DiffSE
            if(DiffSE.lt.MinDiffSE) MinDiffSE = DiffSE
            MeanDiffSE = MeanDiffSE + DiffSE

            if(DiffGF.gt.MaxDiffGF) MaxDiffGF = DiffGF
            MeanDiffGF = MeanDiffGF + DiffGF
        enddo

        MeanDiffSE = MeanDiffSE / real(n,dp)
        MeanDiffGF = MeanDiffGF / real(n,dp)

    end subroutine FindSEDiffs

    !This finds the self energy required to match the local greens function to the interacting local one.
    ! (omega + mu +- idelta)I - e_00 - hybrid - InvG00      Gives the local self-energy
    ! e_00 is taken from the local h0v part of tha hamiltonian
    subroutine FindSE(n,InvG00,Hybrid,SE)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: InvG00(nImp,nImp,n)
        complex(dp), intent(in) :: Hybrid(nImp,nImp,n)
        complex(dp), intent(out) :: SE(nImp,nImp,n)
        complex(dp) :: LocalH(nImp,nImp)
        real(dp) :: Omega,mu
        integer :: i,j
        character(len=*), parameter :: t_r='FindSE'

        !This is *NOT* the correct chemical potential for non-half filled systems, but nevermind...
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        SE(:,:,:) = zzero

        do i = 1,nImp
            do j = 1,nImp
                LocalH(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo

        i = 0
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Too many freq points')

            SE(:,:,i) = - LocalH(:,:) - Hybrid(:,:,i) - InvG00(:,:,i)
            do j = 1,nImp
                SE(j,j,i) = SE(j,j,i) + dcmplx(Omega + mu,dDelta)
            enddo

            Omega = Omega + Omega_Step
        enddo
    end subroutine FindSE

    !This hybridization is a nImp x nImp matrix, for all frequency points, calculated from:
    ! [omega + mu + i\eta]I - e_00 - SE - Inv_GF     where e_00 is the impurity block of the 'non-interacting' hamiltonian
    ! Reminder that here we are performing this self-consistency on the R£AL frequency axis, and the e_00 matrix simply comes from h0v
    ! Hybridization returned in Hybrid
    subroutine FindHybrid(n,InvLocalMomGF,SE,Hybrid)
        implicit none
        integer, intent(in) :: n    !number of frequency points
        complex(dp), intent(in) :: InvLocalMomGF(nImp,nImp,n)
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(out) :: Hybrid(nImp,nImp,n)
        complex(dp) :: LocalH(nImp,nImp)
        real(dp) :: mu,Omega
        integer :: i,j
        character(len=*), parameter :: t_r='FindHybrid'
        
        !This is *NOT* the correct chemical potential for non-half filled systems, but nevermind...
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        Hybrid(:,:,:) = zzero

        do i = 1,nImp
            do j = 1,nImp
                LocalH(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo

        i = 0
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Too many freq points')

            Hybrid(:,:,i) = - LocalH(:,:) - SE(:,:,i) - InvLocalMomGF(:,:,i)
            do j = 1,nImp
                Hybrid(j,j,i) = Hybrid(j,j,i) + dcmplx(Omega + mu,dDelta)
            enddo

            Omega = Omega + Omega_Step
        enddo

    end subroutine FindHybrid
            
    !This function inverts a local greens function matrix, allowing it to be non-hermitian
    !Function is sent in as the normal greens function, n, nImp x nImp matrices, and the inverse
    !for each frequency is returned
    subroutine InvertLocalNonHermGF(n,InvGF)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        use matrixops, only: mat_inv
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(inout) :: InvGF(nImp,nImp,n)
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: LVec(:,:),RVec(:,:),W_Vals(:),cWork(:)
        complex(dp) :: ztemp2(nImp,nImp),IGF(nImp,nImp)
        integer :: i,lwork,info,j
        character(len=*), parameter :: t_r='InvertLocalNonHermGF'

        if(.true.) then
            do i = 1,n
                call mat_inv(InvGF(:,:,i),IGF)
                InvGF(:,:,i) = IGF(:,:)
            enddo
        else
            !Invert via diagonalization
            allocate(LVec(nImp,nImp))
            allocate(RVec(nImp,nImp))
            allocate(W_Vals(nImp))
            allocate(Work(max(1,2*nImp)))

            do i = 1,n

                !Diagonalize
                LVec(:,:) = zzero
                RVec(:,:) = zzero
                W_Vals(:) = zzero
                IGF(:,:) = InvGF(:,:,i)

                allocate(cWork(1))
                lwork = -1
                info = 0
                call zgeev('V','V',nImp,IGF,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(abs(cWork(1)))+1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zgeev('V','V',nImp,IGF,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
                deallocate(cWork)

                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
                call Order_zgeev_vecs(W_Vals,LVec,RVec)
                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
                call Orthonorm_zgeev_vecs(nImp,W_Vals,LVec,RVec)
                !Calculate greens function for this k-vector
                IGF = zzero
                do j = 1,nImp
                    IGF(j,j) = zone / W_Vals(j)
                enddo
                !Now rotate this back into the original basis
                call zGEMM('N','N',nImp,nImp,nImp,zone,RVec,nImp,IGF,nImp,zzero,ztemp2,nImp)
                call zGEMM('N','C',nImp,nImp,nImp,zone,ztemp2,nImp,LVec,nImp,zzero,IGF,nImp)

                InvGF(:,:,i) = IGF(:,:)
            enddo

            deallocate(LVec,RVec,W_Vals,Work)
        endif

    end subroutine InvertLocalNonHermGF

    !The same as FindRealSpaceLocaMomGF below, but doing the diagonalizations in k-space rather than r-space.
    subroutine FindRealSpaceLocalMomGF_Kspace(n,SE,RealSpaceLocGF)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(out) :: RealSpaceLocGF(nImp,nImp,n)
        complex(dp), allocatable :: OneE_Ham(:,:),RotMat(:,:),cWork(:)
        complex(dp), allocatable :: k_Ham(:,:),LVec(:,:),RVec(:,:),W_Vals(:)
        complex(dp), allocatable :: r_vecs_R(:,:),r_vecs_L(:,:),ztemp(:,:)
        real(dp), allocatable :: Work(:)
        real(dp) :: Omega
        integer :: w,i,j,k,pertsite,pertBra,ind_1,ind_2,lwork,info
        character(len=*), parameter :: t_r='FindRealSpaceLocalMomGF_KSpace'

        allocate(OneE_Ham(nSites,nSites))
        allocate(RotMat(nSites,nImp))
        allocate(k_Ham(nImp,nImp))
        allocate(W_Vals(nImp))
        allocate(LVec(nImp,nImp))
        allocate(RVec(nImp,nImp))
        allocate(Work(2*nImp))
        allocate(r_vecs_R(nSites,nImp))
        allocate(r_vecs_L(nSites,nImp))
        allocate(ztemp(nSites,nImp))

        w = 0
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            w = w + 1

!            write(6,*) "*** kspace: ",w,n

            do i = 1,nSites
                do j = 1,nSites
                    OneE_Ham(j,i) = dcmplx(h0v(j,i),zero)
                enddo
            enddo
            !Stripe the complex (-)self-energy through the AO one-electron hamiltonian
            call add_localpot_comp_inplace(OneE_Ham,SE(:,:,w),tAdd=.true.)
            do k = 1,nKPnts
                ind_1 = ((k-1)*nImp) + 1
                ind_2 = nImp*k
                RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
                call ZGEMM('N','N',nSites,nImp,nSites,zone,OneE_Ham,nSites,RotMat,nSites,zzero,ztemp,nSites)
                call ZGEMM('C','N',nImp,nImp,nSites,zone,RotMat,nSites,ztemp,nSites,zzero,k_Ham,nImp)
            
                !Now, diagonalize the resultant non-hermitian one-electron hamiltonian
                RVec = zzero
                LVec = zzero
                W_Vals = zzero
                allocate(cWork(1))
                lwork = -1
                info = 0
                call zgeev('V','V',nImp,k_Ham,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(abs(cWork(1)))+1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zgeev('V','V',nImp,k_Ham,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diag of H - SE failed')
                deallocate(cWork)

                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
                call Order_zgeev_vecs(W_Vals,LVec,RVec)
                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
                call Orthonorm_zgeev_vecs(nImp,W_Vals,LVec,RVec)

                !write(6,*) W_Vals,LVec(:,:),RVec(:,:)

                !Rotate vectors back into r-space
                call ZGEMM('N','N',nSites,nImp,nImp,zone,RtoK_Rot(:,ind_1:ind_2),nSites,  &
                    RVec(:,:),nImp,zzero,r_vecs_R(:,:),nSites)
                call ZGEMM('N','N',nSites,nImp,nImp,zone,RtoK_Rot(:,ind_1:ind_2),nSites,  &
                    LVec(:,:),nImp,zzero,r_vecs_L(:,:),nSites)
            
                do pertsite = 1,nImp
                    do pertBra = 1,nImp
                        do i = 1,nImp
                            RealSpaceLocGF(pertsite,pertBra,w) = RealSpaceLocGF(pertsite,pertBra,w) +   &
                                r_vecs_R(pertBra,i)*dconjg(r_vecs_L(pertsite,i)) / (dcmplx(Omega,dDelta) - W_Vals(i))
                        enddo
                    enddo
                enddo

            enddo
            Omega = Omega + Omega_Step
        enddo
        deallocate(OneE_Ham,W_Vals,RVec,LVec,Work,r_vecs_L,r_vecs_R,k_Ham)
        deallocate(ztemp,RotMat)

    end subroutine FindRealSpaceLocalMomGF_KSpace

    subroutine FindRealSpaceLocalMomGF(n,SE,RealSpaceLocGF)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(out) :: RealSpaceLocGF(nImp,nImp,n)
        complex(dp), allocatable :: OneE_Ham(:,:),W_Vals(:),RVec(:,:),LVec(:,:),cWork(:)
        complex(dp), allocatable :: SelfE(:,:),ztemp(:,:),RotMat(:,:)
        real(dp), allocatable :: Work(:)
        real(dp) :: Omega,mu
        complex(dp) :: NI_LRMat_Ann(nImp,nImp),NI_LRMat_Cre(nImp,nImp)
        integer :: i,j,w,info,lwork,pertBra,pertsite,a,k,ind_1,ind_2
        character(len=*), parameter :: t_r='FindRealSpaceLocalMomGF'
        
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif
            
        allocate(OneE_Ham(nSites,nSites))
        allocate(W_Vals(nSites))
        allocate(RVec(nSites,nSites))
        allocate(LVec(nSites,nSites))
        allocate(Work(max(1,2*nSites)))

        allocate(SelfE(nSites,nSites))
        allocate(ztemp(nSites,nImp))
        allocate(RotMat(nSites,nImp))

        w = 0
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            w = w + 1
            write(6,*) "***",w,n

            do i = 1,nSites
                do j = 1,nSites
                    OneE_Ham(j,i) = dcmplx(h0v(j,i),zero)
                enddo
            enddo
            !Stripe the complex (-)self-energy through the AO one-electron hamiltonian
            call add_localpot_comp_inplace(OneE_Ham,SE(:,:,w),tAdd=.true.)
            !Instead, rotate each k-vector component back into r-space and add to the hamiltonain
            !do k = 1,nKPnts
            !    ind_1 = ((k-1)*nImp) + 1
            !    ind_2 = nImp*k
            !    RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
            !    call ZGEMM('N','N',nSites,nImp,nImp,zone,RotMat,nSites,SE(:,:,w),nImp,zzero,ztemp,nSites)
            !    call ZGEMM('N','C',nSites,nSites,nImp,zone,ztemp,nSites,RotMat,nSites,zzero,SelfE,nSites)
            !    !Add SelfE to the hamiltonian
            !    OneE_Ham(:,:) = OneE_Ham(:,:) + SelfE(:,:)
            !enddo


            !Now, diagonalize the resultant non-hermitian one-electron hamiltonian
            RVec = zzero
            LVec = zzero
            W_Vals = zzero
            allocate(cWork(1))
            lwork = -1
            info = 0
            call zgeev('V','V',nSites,OneE_Ham,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Workspace query failed')
            lwork = int(abs(cWork(1)))+1
            deallocate(cWork)
            allocate(cWork(lWork))
            call zgeev('V','V',nSites,OneE_Ham,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diag of H - SE failed')
            deallocate(cWork)

            !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
            !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
            call Order_zgeev_vecs(W_Vals,LVec,RVec)
            !call writevectorcomp(W_Vals,'Eigenvalues ordered')
            !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
            call Orthonorm_zgeev_vecs(nSites,W_Vals,LVec,RVec)

            NI_LRMat_Ann(:,:) = zzero
            NI_LRMat_Cre(:,:) = zzero
            do pertsite = 1,nImp
                do pertBra = 1,nImp
                    !Now perform the set of dot products of <0|V* with |1> for all combinations of sites
                    do i = 1,nOcc
                        NI_LRMat_Ann(pertsite,pertBra) = NI_LRMat_Ann(pertsite,pertBra) +   &
                            RVec(pertBra,i)*dconjg(LVec(pertsite,i)) / (dcmplx(Omega + mu,dDelta) - W_Vals(i))
                    enddo
                    do a = nOcc+1,nSites
                        NI_LRMat_Cre(pertsite,pertBra) = NI_LRMat_Cre(pertsite,pertBra) +   &
                            RVec(pertBra,a)*dconjg(LVec(pertsite,a)) / (dcmplx(Omega + mu,dDelta) - W_Vals(a))
                    enddo
                enddo
            enddo

            RealSpaceLocGF(:,:,w) = NI_LRMat_Ann(:,:) + NI_LRMat_Cre(:,:)
            Omega = Omega + Omega_Step
        enddo
        deallocate(OneE_Ham,W_Vals,RVec,LVec,Work)
        deallocate(SelfE,ztemp,RotMat)

    end subroutine FindRealSpaceLocalMomGF
        
    !Get the eigenvalues from the hamiltonian
    subroutine CouplingsToEVals(hLat,EVals)
        use mat_tools, only: DiagOneEOp
        implicit none
        real(dp), intent(in) :: hLat(nSites,nSites)
        real(dp), intent(out) :: EVals(nSites)
        integer :: SS_Period,kPnt,ind_1,ind_2,info,lWork,i,j
        complex(dp), allocatable :: zTemp(:,:),k_Ham(:,:),cWork(:),CompHam(:,:)
        real(dp), allocatable :: rWork(:),ham_temp(:,:)
        character(len=*), parameter :: t_r='CouplingsToEVals'

        if(.not.tDiag_KSpace) then
            allocate(ham_temp(nSites,nSites))
            ham_temp(:,:) = hLat(:,:)
            call DiagOneEOp(ham_temp,EVals,nImp,nSites,.false.)
            deallocate(ham_temp)    !ham_temp contains the real space eigenvectors
            return
        endif

        SS_Period = nImp

        allocate(zTemp(nSites,SS_Period))
        allocate(rWork(max(1,3*SS_Period-2)))
        allocate(k_Ham(SS_Period,SS_Period))
        allocate(CompHam(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                CompHam(j,i) = dcmplx(hLat(j,i),zero)
            enddo
        enddo

        do kPnt = 1,nKPnts
            ind_1 = ((kPnt-1)*SS_Period) + 1
            ind_2 = SS_Period*kPnt

            call ZGEMM('N','N',nSites,SS_Period,nSites,zone,CompHam,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ztemp,nSites)
            call ZGEMM('C','N',SS_Period,SS_Period,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ztemp,nSites,zzero,k_Ham,SS_Period)

            do i = 1,SS_Period
                if(abs(aimag(k_ham(i,i))).gt.1.0e-7_dp) then
                    call writematrixcomp(k_Ham,'k_ham',.false.)
                    call stop_all(t_r,"Real part of k-transformed hamiltonian diagonal too large")
                else
                    k_ham(i,i) = dcmplx(real(k_ham(i,i)),zero)
                endif
            enddo

            if(SS_Period.gt.1) then
                allocate(cWork(1))
                lWork = -1
                info = 0 
                !write(6,*) "KPnt: ",kPnt
                !call writematrixcomp(k_Ham,'k_ham',.false.)
                call zheev('V','U',SS_Period,k_Ham,SS_Period,EVals(ind_1:ind_2),cWork,lWork,rWork,info)
                if(info.ne.0) then
                    write(6,*) "INFO = ", info
                    write(6,*) "CWORK = ",cWork(:)
                    call stop_all(t_r,'Error in diag 1')
                endif
                lWork = int(real(cWork(1))) + 1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zheev('V','U',SS_Period,k_Ham,SS_Period,EVals(ind_1:ind_2),cWork,lWork,rWork,info)
                if(info.ne.0) call stop_all(t_r,'Error in diag 2')
                deallocate(cWork)
            else
                EVals(ind_1:ind_2) = real(k_Ham(1,1),dp)
            endif
        enddo

        deallocate(k_Ham,zTemp,CompHam,rWork)

    end subroutine CouplingsToEVals
            
    subroutine FindFullLatEvals(evals,nInd_evals,evals_full)
        implicit none
        integer, intent(in) :: nInd_evals
        real(dp), intent(in) :: evals(nInd_evals)
        real(dp), intent(out) :: evals_full(nSites)
        integer :: i
        character(len=*), parameter :: t_r='FindFullLatEvals'

        if(tKPntSymFit) then
            if(tShift_Mesh) then
                if(nInd_evals.ne.(nSites/2)) call stop_all(t_r,'Number of independent evals incorrect?')
            else
                if(nInd_evals.ne.((nSites/2)+1)) call stop_all(t_r,'# independent evals incorrect?')
            endif

            !write(6,*) "tShift_Mesh: ",tShift_Mesh,nInd_evals,nSites,size(evals_full),size(evals)
            evals_full(1:nSites/2) = evals(1:nSites/2)
            if(tShift_Mesh) then
                !No gamma point. All k-points symmetric
                do i = 1,nSites/2
                    evals_full(i+(nSites/2)) = evals((nSites/2)-i+1)
                enddo
            else
                !Put in the gamma point
                evals_full((nSites/2)+1) = evals((nSites/2)+1)
                do i = 2,nSites/2
                    evals_full(i+(nSites/2)) = evals((nSites/2)-i+2)
                enddo
            endif
        else
            if(nInd_evals.ne.nSites) call stop_all(t_r,'Something wrong here')
            evals_full(:) = evals(:)
        endif

    end subroutine FindFullLatEvals
        
    !This calculates the sum of all the k-space greens functions over the 'n' frequency points
    !The hamiltonian in this case is taken to be h0v (i.e. with the correlation potential)
    !1/V \sum_k [w + mu + i/eta - h_k - SE]^-1      where V is volume of the BZ
    !Check that this gives as expected
    !Can optionally be computed in real (default) or matsubara axis
    !Can optionally take a periodic lattice hamiltonian on which to calculate greens function (will then not use self-energy)
    !Can optionally take a coupling length for the non-local terms of the matrix, which should speed up the FT to k-space
    subroutine FindLocalMomGF(n,SE,LocalMomGF,tMatbrAxis,ham,CouplingLength,evals,FreqPoints)
        use matrixops, only: mat_inv
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(out) :: LocalMomGF(nImp,nImp,n)
        logical, intent(in), optional :: tMatbrAxis
        real(dp), intent(in), optional :: ham(nSites,nSites)
        integer, intent(in), optional :: CouplingLength
        real(dp), intent(in), optional :: evals(:)
        real(dp), intent(in), optional :: FreqPoints(n)
        complex(dp), allocatable :: k_Ham(:,:),CompHam(:,:),cWork(:)
        complex(dp), allocatable :: RVec(:,:),LVec(:,:),W_Vals(:),ztemp(:,:)
        complex(dp) :: InvMat(nImp,nImp),ztemp2(nImp,nImp)
        real(dp), allocatable :: Work(:), evals_full(:)
        integer :: kPnt,ind_1,ind_2,i,j,SS_Period,lwork,info,k
        real(dp) :: Omega,mu
        logical :: tMatbrAxis_,tCoupledHam_,tSparseCoupling,tUseEVals
        character(len=*), parameter :: t_r='FindLocalMomGF'

        if(present(tMatbrAxis)) then
            tMatbrAxis_ = tMatbrAxis
        else
            tMatbrAxis_ = .false.
        endif

        if(present(ham)) then
            tCoupledHam_ = .true.
        else
            tCoupledHam_ = .false.
        endif

        if(present(CouplingLength)) then
            tSparseCoupling = .true.
        else
            tSparseCoupling = .false.
        endif

        if(present(evals).and.tOptGF_EVals) then
            tUseEVals = .true.
        else
            tUseEVals = .false.
        endif

!        if(.not.tDiag_kspace) then
!            if(tMatbrAxis_.or.tCoupledHam_) call stop_all(t_r,'Fixme')
!            call FindRealSpaceLocalMomGF(n,SE,LocalMomGF)
!            return
!        endif

        LocalMomGF(:,:,:) = zzero
        SS_Period = nImp
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            !This should really be passed in
            mu = U/2.0_dp
        else
            mu = zero 
        endif

        if(tUseEVals) then
            if(nImp.gt.1) call stop_all(t_r,'Not working for > 1 imp')

            allocate(evals_full(nSites))
            if(.not.present(CouplingLength)) call stop_all(t_r,'Error')
            call FindFullLatEvals(evals,CouplingLength,evals_full)
            
            i = 0
            do while(.true.)
                if(present(FreqPoints)) then
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_,nFreqPoints=n,FreqPoints=FreqPoints)
                else
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_)
                endif
                if(i.lt.0) exit
                if(i.gt.n) then
                    write(6,*) "i: ",i
                    write(6,*) "n: ",n
                    call stop_all(t_r,'Too many freq points')
                endif

                if(tDiag_kspace) then
                    if(tMatbrAxis_) then
                        do k = 1,nSites
                            LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + dconjg(RtoK_Rot(1,k))*RtoK_Rot(1,k) /   &
                                dcmplx(mu - evals_full(k),Omega)
                        enddo
                    else
                        do k = 1,nSites
                            LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + dconjg(RtoK_Rot(1,k))*RtoK_Rot(1,k) /   &
                                dcmplx(Omega + mu - evals_full(k),dDelta)
                        enddo
                    endif
                else
                    if(tMatbrAxis_) then
                        do k = 1,nSites
                            LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + dcmplx(HFOrbs(1,k)**2,zero) /   &
                                dcmplx(mu - evals_full(k),Omega)
                        enddo
                    else
                        do k = 1,nSites
                            LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + dcmplx(HFOrbs(1,k)**2,zero) /   &
                                dcmplx(Omega + mu - evals_full(k),dDelta)
                        enddo
                    endif
                endif
            enddo

            deallocate(evals_full)
            return
        endif

        !allocate(RotMat(nSites,SS_Period))
        allocate(k_Ham(SS_Period,SS_Period))
        allocate(CompHam(nSites,nSites))
        if(tCoupledHam_) then
            !Use passed in lattice hamiltonian
            do i = 1,nSites
                do j = 1,nSites
                    CompHam(j,i) = dcmplx(ham(j,i),zero)
                enddo
            enddo
        else
            do i = 1,nSites
                do j = 1,nSites
                    CompHam(j,i) = dcmplx(h0v(j,i),zero)
                enddo
            enddo
        endif

        !Data for the diagonalization of the one-electron k-dependent Greens functions
!        allocate(RVec(SS_Period,SS_Period))
!        allocate(LVec(SS_Period,SS_Period))
!        allocate(W_Vals(SS_Period))
!        allocate(Work(max(1,2*SS_Period)))
        allocate(ztemp(nSites,SS_Period))

        do kPnt = 1,nKPnts
            !Run through all k-points
            ind_1 = ((kPnt-1)*SS_Period) + 1
            ind_2 = SS_Period*kPnt
                
            !Transform the one-electron hamiltonian into this k-basis
            !RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
            if(tSparseCoupling) then
                !Exploit the fact that the real-space lattice is sparse, to only FT the non-zero elements in real space
                ztemp(:,:) = zzero
                do i = 1,nSites
                    do j = 1,SS_Period

                        !Sparse summation over non-zero elements of the real space lattice only, rather than full coupling
                        do k = max(i-CouplingLength,1),min(i+CouplingLength,nSites)
                            ztemp(i,j) = ztemp(i,j) + CompHam(i,k)*RtoK_Rot(k,ind_1+j-1)
                            !write(6,*) "Main contrib: ",k,CompHam(i,k),RtoK_Rot(k,ind_1+j-1),CompHam(i,k)*RtoK_Rot(k,ind_1+j-1)
                        enddo

                        if(i-CouplingLength.lt.1) then
                            !We are wrapping around the lattice. Include the 'top right' triangle of couplings
                            do k = nSites-(CouplingLength-i),nSites
                                ztemp(i,j) = ztemp(i,j) + CompHam(i,k)*RtoK_Rot(k,ind_1+j-1)
                                !write(6,*) "TR: ",k,CompHam(i,k),RtoK_Rot(k,ind_1+j-1),CompHam(i,k)*RtoK_Rot(k,ind_1+j-1)
                            enddo

                        elseif(i+CouplingLength.gt.nSites) then
                            !We are wrapping around the lattice. Include the 'bottom left' triangle of couplings
                            do k = 1,(i-nSites)+CouplingLength
                                ztemp(i,j) = ztemp(i,j) + CompHam(i,k)*RtoK_Rot(k,ind_1+j-1)
                                !write(6,*) "BL: ",k,CompHam(i,k),RtoK_Rot(k,ind_1+j-1),CompHam(i,k)*RtoK_Rot(k,ind_1+j-1)
                            enddo

                        endif
                    enddo
                enddo
            else
                call ZGEMM('N','N',nSites,SS_Period,nSites,zone,CompHam,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ztemp,nSites)
            endif
!            call writematrixcomp(RtoK_Rot(:,ind_1:ind_2),'Rotation',.true.)
!            call writematrixcomp(ztemp,'Intermediate',.true.)
            call ZGEMM('C','N',SS_Period,SS_Period,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ztemp,nSites,zzero,k_Ham,SS_Period)
            !k_Ham is the complex, one-electron hamiltonian (with correlation potential) for this k-point

            !write(6,*) "For k-point: ",kPnt
            !call writematrixcomp(k_Ham,'k-space ham is: ',.false.)

            i = 0
            do while(.true.)
                if(present(FreqPoints)) then
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_,nFreqPoints=n,FreqPoints=FreqPoints)
                else
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_)
                endif
                if(i.lt.0) exit
                if(i.gt.n) call stop_all(t_r,'Too many freq points')

                !Find inverse greens function for this k-point
                if(tCoupledHam_) then
                    !No self-energy
                    InvMat(:,:) = - k_Ham(:,:) 
                else
                    InvMat(:,:) = - k_Ham(:,:) - SE(:,:,i)
                endif
                do j = 1,SS_Period
                    if(tMatbrAxis_) then
                        InvMat(j,j) = InvMat(j,j) + dcmplx(mu,Omega)
                    else
                        InvMat(j,j) = InvMat(j,j) + dcmplx(Omega + mu,dDelta)
                    endif
                enddo

!                !Now, diagonalize this
!                !The self-energy is *not* hermitian
!                allocate(cWork(1))
!                lwork = -1
!                info = 0
!                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
!                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
!                lwork = int(abs(cWork(1)))+1
!                deallocate(cWork)
!                allocate(cWork(lWork))
!                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
!                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
!                deallocate(cWork)
!
!                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
!                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
!                call Order_zgeev_vecs(W_Vals,LVec,RVec)
!                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
!                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
!                call Orthonorm_zgeev_vecs(SS_Period,W_Vals,LVec,RVec)
!                !Calculate greens function for this k-vector
!                InvMat(:,:) = zzero
!                do j = 1,SS_Period
!                    InvMat(j,j) = zone / W_Vals(j)
!                enddo
!                !Now rotate this back into the original basis
!                call zGEMM('N','N',SS_Period,SS_Period,SS_Period,zone,RVec,SS_Period,InvMat,SS_Period,zzero,ztemp2,SS_Period)
!                call zGEMM('N','C',SS_Period,SS_Period,SS_Period,zone,ztemp2,SS_Period,LVec,SS_Period,zzero,InvMat,SS_Period)

                !Instead of explicit diagonalization, we can just invert
                call mat_inv(InvMat,ztemp2)

                !InvMat is now the non-interacting greens function for this k: TEST THIS!
                !Sum this into the *local* greens function (i.e. a fourier transform of r=0 component)
                LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + ztemp2(:,:)

                !write(6,*) "For frequency: ",Omega
                !call writematrixcomp(InvMat,'k-space greens function is: ',.false.)
            enddo   !Enddo frequency point i

        enddo   !enddo k-point

        !write(6,*) "BZ volume: ",BZVol

        !Divide the entire local greens function by the 'volume' of the brillouin zone
!        LocalMomGF(:,:,:) = LocalMomGF(:,:,:) / BZVol
        LocalMomGF(:,:,:) = LocalMomGF(:,:,:) / real(nSites/nImp,dp)
        
        deallocate(k_Ham,CompHam,ztemp)

    end subroutine FindLocalMomGF
    
    !Find the jacobian matrix for the optimization (in a few cases)
    subroutine CalcJacobian2(n,CouplingLength,DiffMat,Jacobian,evals,tMatbrAxis,FreqPoints,Weights)
        implicit none
        integer, intent(in) :: n    !Number of frequency points
        integer, intent(in) :: CouplingLength   !The number of variables
        complex(dp), intent(in) :: DiffMat(nImp,nImp,n) !The (complex) difference between the greens functions
        real(dp), intent(out) :: Jacobian(CouplingLength) !The Jacobian
        real(dp), intent(in) :: evals(CouplingLength)
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in), optional :: FreqPoints(n),Weights(n)
        
        real(dp) :: mu,Omega,denom1,denom2,num,termr,termi
        complex(dp) :: compval
        real(dp), allocatable :: evals_full(:),FullJac(:)
        integer :: i,k
        logical :: tNonStandardGrid
        character(len=*), parameter :: t_r='CalcJacobian2'

        if(nImp.gt.1) call stop_all(t_r,'Cannot do for more than 1 impurity yet')
        if(.not.tMatbrAxis) call stop_all(t_r,'Cannot do gradients on real axis yet')
        if(iLatticeFitType.ne.1) call stop_all(t_r,'Cannot do gradients with non-linear objective functions yet')
        if(.not.tOptGF_EVals) call stop_all(t_r,'Cannot do gradients when optimizing lattice couplings rather than eigenvalues')
        if(present(FreqPoints)) then
            tNonStandardGrid = .true.
            if(.not.(present(Weights))) then
                call stop_all(t_r,"No Weights present")
            endif
        else
            tNonStandardGrid = .false.
        endif
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            !This should really be passed in
            mu = U/2.0_dp
        else
            mu = zero 
        endif

        allocate(evals_full(nSites))
        call FindFullLatEvals(evals,CouplingLength,evals_full)
        !Calculate the *full* hamiltonian
        allocate(FullJac(nSites))
        FullJac(:) = zero

        do k = 1,nSites
            i=0
            do while(.true.)
                if(tNonStandardGrid) then
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis,nFreqPoints=n,FreqPoints=FreqPoints)
                else
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis)
                endif
                if(i.lt.0) exit
                
                if(tDiag_kSpace) then
                    num = abs(RtoK_Rot(1,k))**2
                else
                    num = HFOrbs(1,k)**2
                endif

                compval = num*dconjg(DiffMat(1,1,i)) / (cmplx((mu - evals_full(k)),Omega,dp)**2) 

                !denom1 = Omega**2 + (mu - evals_full(k))**2
                !denom2 = denom1**2

                !termr = 2.0_dp*num*((mu - evals_full(k))**2)/denom2 - num/denom1
                !termi = -2.0_dp*Omega*(mu - evals_full(k))*num/denom2
                !compval = cmplx(termr,termi,dp)*dconjg(DiffMat(1,1,i))
                
                if(tNonStandardGrid) then
                    !Add the weighting factor
                    if(iFitGFWeighting.eq.0) then
                        FullJac(k) = FullJac(k) + real(compval + dconjg(compval),dp) * Weights(i)
                    elseif(iFitGFWeighting.eq.1) then
                        FullJac(k) = FullJac(k) + real(compval + dconjg(compval),dp) * Weights(i) / Omega
                    elseif(iFitGFWeighting.eq.2) then
                        FullJac(k) = FullJac(k) + real(compval + dconjg(compval),dp) * Weights(i) / (Omega**2)
                    endif
                else
                    if(iFitGFWeighting.eq.0) then
                        FullJac(k) = FullJac(k) + real(compval + dconjg(compval),dp)
                    elseif(iFitGFWeighting.eq.1) then
                        FullJac(k) = FullJac(k) + real(compval + dconjg(compval),dp) / Omega
                    elseif(iFitGFWeighting.eq.2) then
                        FullJac(k) = FullJac(k) + real(compval + dconjg(compval),dp) / (Omega**2)
                    endif
                endif
            enddo
        enddo

        !Now, compress the jacobian if necessary
        if(tKPntSymFit) then
            !However, we need now to package it up. Since we are only optimizing certain eigenvalue, we need
            !to add the two gradient contributions together
            Jacobian(:) = zero
            !Include the negative k terms
            Jacobian(1:nSites/2) = FullJac(1:nSites/2)
            if(tShift_Mesh) then
                !No gamma point. All k-points symmetric
                do i = 1,nSites/2
                    Jacobian((nSites/2)-i+1) = Jacobian((nSites/2)-i+1) + FullJac(i+(nSites/2))
                enddo
            else
                !The gamma point is unchanged (and therefore the gradient will on average be half the other values)
                Jacobian((nSites/2)+1) = FullJac((nSites/2)+1)
                do i = 2,nSites/2
                    Jacobian(((nSites/2)-i+2)) = Jacobian(((nSites/2)-i+2)) + FullJac(i+(nSites/2))
                enddo
            endif
        else
            if(CouplingLength.ne.nSites) call stop_all(t_r,'Error here')
            Jacobian(:) = FullJac(:)
            !write(6,*) "Jacobian: ",Jacobian(:)
        endif
        deallocate(FullJac)

    end subroutine CalcJacobian2
                
    !Find the jacobian matrix for the optimization (in a few cases)
    subroutine CalcJacobian(n,CouplingLength,DiffMatr,DiffMat,Jacobian,evals,tMatbrAxis,FreqPoints,Weights)
        implicit none
        integer, intent(in) :: n    !Number of frequency points
        integer, intent(in) :: CouplingLength   !The number of variables
        complex(dp), intent(in) :: DiffMat(nImp,nImp,n) !The (complex) difference between the greens functions
        real(dp), intent(in) :: DiffMatr(nImp,nImp,n)   !The abs difference between the greens functions
        real(dp), intent(out) :: Jacobian(n,CouplingLength) !The compressed Jacobian
        real(dp), intent(in) :: evals(CouplingLength)
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in), optional :: FreqPoints(n),Weights(n)

        real(dp) :: mu,Omega,denom1,denom2,num,termr,termi
        real(dp), allocatable :: evals_full(:),FullJac(:,:)
        integer :: i,k
        logical :: tNonStandardGrid
        character(len=*), parameter :: t_r='CalcJacobian'

        if(nImp.gt.1) call stop_all(t_r,'Cannot do for more than 1 impurity yet')
        if(.not.tMatbrAxis) call stop_all(t_r,'Cannot do gradients on real axis yet')
        if(iLatticeFitType.ne.1) call stop_all(t_r,'Cannot do gradients with non-linear objective functions yet')
        if(.not.tOptGF_EVals) call stop_all(t_r,'Cannot do gradients when optimizing lattice couplings rather than eigenvalues')
        if(present(FreqPoints)) then
            tNonStandardGrid = .true.
            if(.not.(present(Weights))) then
                call stop_all(t_r,"No Weights present")
            endif
        else
            tNonStandardGrid = .false.
        endif
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            !This should really be passed in
            mu = U/2.0_dp
        else
            mu = zero 
        endif

        allocate(evals_full(nSites))
        call FindFullLatEvals(evals,CouplingLength,evals_full)
        !Calculate the *full* hamiltonian
        allocate(FullJac(n,nSites))
        FullJac(:,:) = zero

        do k = 1,nSites
            i=0
            do while(.true.)
                if(tNonStandardGrid) then
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis,nFreqPoints=n,FreqPoints=FreqPoints)
                else
                    call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis)
                endif
                if(i.lt.0) exit

                denom1 = Omega**2 + (mu - evals_full(k))**2
                denom2 = denom1**2
                if(tDiag_kSpace) then
                    num = abs(RtoK_Rot(1,k))**2
                else
                    num = HFOrbs(1,k)**2
                endif

                termr = 2.0_dp*num*((mu - evals_full(k))**2)/denom2 - num/denom1
                termi = 2.0_dp*Omega*(mu - evals_full(k))*num/denom2

                FullJac(i,k) = (one/DiffMatr(1,1,i))*(real(DiffMat(1,1,i),dp)*termr - aimag(DiffMat(1,1,i))*termi)
                !FullJac(i,k) = 2.0_dp*(real(DiffMat(1,1,i),dp)*termr - aimag(DiffMat(1,1,i))*termi)

                !if(abs(Omega).lt.1.0e-9_dp) then
                !    write(6,*) "mod: ",DiffMatr(1,1,i)
                !    write(6,*) "Real: ",real(DiffMat(1,1,i),dp)
                !    write(6,*) "termr: ",termr
                !endif

                if(tNonStandardGrid) then
                    !Add the weighting factor
                    FullJac(i,k) = FullJac(i,k) * sqrt(Weights(i))
                    !FullJac(i,k) = FullJac(i,k) * Weights(i)
                endif
                if(iFitGFWeighting.eq.1) then
                    !FullJac(i,k) = FullJac(i,k)/abs(Omega)
                    FullJac(i,k) = FullJac(i,k)/sqrt(abs(Omega))
                elseif(iFitGFWeighting.eq.2) then
                    !FullJac(i,k) = FullJac(i,k)/(Omega**2)
                    FullJac(i,k) = FullJac(i,k)/abs(Omega)
                endif
            enddo
        enddo

        if(tKPntSymFit) then
            !However, we need now to package it up. Since we are only optimizing certain eigenvalue, we need
            !to add the two gradient contributions together
            Jacobian(:,:) = zero
            !Include the negative k terms
            Jacobian(:,1:nSites/2) = FullJac(:,1:nSites/2)
            if(tShift_Mesh) then
                !No gamma point. All k-points symmetric
                do i = 1,nSites/2
                    Jacobian(:,(nSites/2)-i+1) = Jacobian(:,(nSites/2)-i+1) + FullJac(:,i+(nSites/2))
                enddo
            else
                !The gamma point is unchanged (and therefore the gradient will on average be half the other values)
                Jacobian(:,(nSites/2)+1) = FullJac(:,(nSites/2)+1)
                do i = 2,nSites/2
                    Jacobian(:,((nSites/2)-i+2)) = Jacobian(:,((nSites/2)-i+2)) + FullJac(:,i+(nSites/2))
                enddo
            endif
        else
            if(CouplingLength.ne.nSites) call stop_all(t_r,'Error here')
            Jacobian(:,:) = FullJac(:,:)
        endif
        deallocate(FullJac)

    end subroutine CalcJacobian

    !From numerical recepies - Find Legendre points and weights between interval
    subroutine gauleg(x1,x2,x,w,n)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: x1,x2
        real(dp), intent(out) :: x(n),w(n)
        real(dp), parameter :: eps=3.0e-14_dp   !Relative precision
        !Given the lower and upper limits of integration x1 and x2,
        !and given n, this routine returns arrays x(1:n) and w(1:n) of length
        !n, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula
        integer :: i,j,m
        real(dp) :: p1,p2,p3,pp,xl,xm,z,z1
        
        m = (n+1)/2     !Roots are symmetric in the interval, so we only have to find half of them
        xm = 0.5_dp*(x2+x1)
        xl = 0.5_dp*(x2-x1)
        do i = 1,m      !Loop over desired roots
            z = cos(pi*(i-0.25_dp)/(n + 0.5_dp))
            !Starting with the above approximation to the ith root, we enter the
            !main loop of refinement by Newtons method
1           continue
            p1 = one
            p2 = zero
            do j = 1,n
                p3 = p2
                p2 = p1
                p1 = ((2.0_dp*j-1.0_dp)*z*p2-(j-one)*p3)/j
            enddo
            !p1 is now the desired Legendre polynomial. We next computer pp,
            !its derivative, by a standard relation involving also p2, the
            !polynomial of one lower order
            pp = n*(z*p1-p2)/(z*z-one)
            z1 = z
            z = z1-p1/pp                    !Newtons method
            if(abs(z-z1).gt.eps) goto 1
            x(i) = xm-xl*z                  !Scale the root by the desired inteval
            x(n+1-i) = xm+xl*z              !and put in its symmetric counterpart
            w(i) = 2.0_dp*xl/((one-z*z)*pp*pp)  !Compute weight
            w(n+1-i) = w(i)                     !and symmetric counterpart
        enddo

    end subroutine gauleg 

    !Scale the legendre abscissas and weights to the range \pm infty
    subroutine ScaleFreqPointsWeights(n,f,w)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(inout) :: w(n),f(n)
        integer :: i

        do i = 1,n
            w(i) = w(i)*((one + f(i)**2)/(one - f(i)**2)**2)
            f(i) = -f(i)/(one-f(i)**2)
        enddo
    end subroutine ScaleFreqPointsWeights

end module SelfConsistentLR
