module SelfConsistentLR2
    use const
    use errors
    use LinearResponse
    use globals
    use utils, only: get_free_unit,append_ext_real,append_ext
    use SC_Data                  
    use SelfConsistentUtils
    use matrixops, only: mat_inv
    implicit none

    contains

    !Self-consistently fit a frequency *independent* lattice coupling to the impurity in order to
    !match the matsubara greens functions
    subroutine SC_Spectrum_Opt()            
        implicit none
        real(dp) :: GFChemPot,Omega
        integer :: nFreq_Re,nFreq_Im,nFitPoints
        complex(dp), allocatable :: h_lat_fit(:,:)
        real(dp), allocatable :: LatParams(:)
        complex(dp), allocatable :: CorrFn_Fit(:,:,:),CorrFn_HL(:,:,:),CorrFn_HL_Old(:,:,:)
        complex(dp), allocatable :: CorrFn_Fit_Old(:,:,:),DiffImpCorrFn(:,:,:)
        complex(dp), allocatable :: dummy_Im(:,:,:),dummy_Re(:,:,:),CorrFn_HL_Re(:,:,:),CorrFn_Re(:,:,:)
        complex(dp), allocatable :: Debug_Lat_CorrFn_Re(:,:,:),Debug_Lat_CorrFn_Fit(:,:,:)
        real(dp), allocatable :: AllDiffs(:,:)

        real(dp) :: FinalDist,LowFreq,HighFreq
        integer :: i,iter,iLatParams,j,iCorrFnTag
        real(dp), parameter :: dDeltaImpThresh = 1.0e-4_dp 
        character(len=*), parameter :: t_r='SC_Spectrum_Opt'

        write(6,"(A)") ""
        write(6,"(A)") " *** Entering code of self-consistent optimizating of spectral functions ... ***"
        write(6,"(A)") ""

        iCorrFnTag = 1    !1 for greens function optimization
        
        if(iCorrFnTag.eq.1) then
            write(6,"(A)") "Single particle greens functions to be optimized."
            if(tDiagonalSC) then
                write(6,"(A)") "Fitting only the diagonal impurity greens function contributions."
            else
                write(6,"(A)") "Fitting the entire impurity local greens function matrix."
            endif
        endif

        if(.not.tDiag_kSpace) then
            call stop_all(t_r,"Unfortunately, this is only for systems with a kspace representation currently")
        endif

        !This will set FreqPoints and Weights. nFitPoints is the size of FreqPoints/Weights, 
        !and may refer to real or im axis based on the tFitRealFreq flag.
        call SetFreqPoints(nFreq_Re,nFreq_Im,nFitPoints)
        call SetChemPot(GFChemPot)

        !How many lattice parameters are there to fit?
        call SetLatticeParams(iLatParams)

        write(6,"(A,I7)") "Number of lattice hamiltonian parameters to optimize: ",iLatParams
        allocate(LatParams(iLatParams))
        allocate(h_lat_fit(nSites,nSites))  !The real-space (complex) lattice hamiltonian

        !Initialize the lattice parameters to be those from the ground-state optimization,
        !and return the initial lattice hamiltonian
        call InitLatticeParams(iLatParams,LatParams,GFChemPot,h_lat_fit)

        write(6,"(A)") "Starting *frequency independent* lattice parameters: "
        do i = 1,iLatParams
            write(6,"(I6,G20.13)") i,LatParams(i)
        enddo
        write(6,"(A)") "" 

        allocate(CorrFn_Fit(nImp,nImp,nFitPoints))
        allocate(CorrFn_Fit_Old(nImp,nImp,nFitPoints))
        allocate(CorrFn_HL(nImp,nImp,nFitPoints))
        allocate(CorrFn_HL_Old(nImp,nImp,nFitPoints))
        allocate(DiffImpCorrFn(nImp,nImp,nFitPoints))
        allocate(dummy_Im(nImp,nImp,nFitPoints))
        allocate(Debug_Lat_CorrFn_Fit(nImp,nImp,nFitPoints))
        dummy_Im(:,:,:) = zzero
        CorrFn_HL(:,:,:) = zzero
        if(tCalcRealSpectrum) then
            allocate(dummy_Re(nImp,nImp,nFreq_Re))
            allocate(CorrFn_HL_Re(nImp,nImp,nFreq_Re))
            allocate(CorrFn_Re(nImp,nImp,nFreq_Re))
            allocate(Debug_Lat_CorrFn_Re(nImp,nImp,nFreq_Re))
            dummy_Re(:,:,:) = zzero
        endif

!        call writematrixcomp(h_lat_fit,'Initial real space matrix',.true.)
        call CalcLatticeSpectrum(iCorrFnTag,nFitPoints,CorrFn_Fit,GFChemPot,tMatbrAxis=tFitMatAxis, &
            iLatParams=iLatParams,LatParams=LatParams,FreqPoints=FreqPoints)

        allocate(AllDiffs(3,0:iMaxIter_MacroFit+1))
        AllDiffs(:,:) = zero

        iter = 0
        do while(.not.tSkip_Lattice_Fit)
            iter = iter + 1

            call writedynamicfunction(nFitPoints,CorrFn_Fit,'G_Lat_Fit',tag=iter,tCheckCausal=.true.,   &
                tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)

!            write(6,*) "For iteration ",iter
!            write(6,*) "First nImp rows of the lattice hamiltonian: "
!            do i = 1,nImp
!                write(6,*) h_lat_fit(i,:)
!                write(6,*) "***"
!            enddo

            CorrFn_HL_Old(:,:,:) = CorrFn_HL(:,:,:)
            if(iCorrFnTag.eq.1) then
!                call writematrixcomp(h_lat_fit,'real space matrix sent to LR',.true.)
                call SchmidtGF_wSE(CorrFn_HL,GFChemPot,dummy_Im,nFitPoints,tMatbrAxis=tFitMatAxis,  &
                    cham=h_lat_fit,FreqPoints=FreqPoints,Lat_G_Mat=Debug_Lat_CorrFn_Fit)
                call writedynamicfunction(nFitPoints,Debug_Lat_CorrFn_Fit,'G_LatSchmidt_Fit',tag=iter,    &
                    tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
                !call CheckGFsSame(nFitPoints,Debug_Lat_CorrFn_Fit,CorrFn_Fit,1.0e-7_dp)
                call writedynamicfunction(nFitPoints,CorrFn_HL,'G_Imp_Fit',tag=iter,    &
                    tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
            else
                call stop_all(t_r,'Non GF correlation functions not yet coded up')
            endif

            if(tCalcRealSpectrum) then
                call SchmidtGF_wSE(CorrFn_HL_Re,GFChemPot,dummy_Re,nFreq_Re,tMatbrAxis=.false., &
                    cham=h_lat_fit,Lat_G_Mat=Debug_Lat_CorrFn_Re)
                call writedynamicfunction(nFreq_Re,CorrFn_HL_Re,'G_Imp_Re',tag=iter,    &
                    tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.)
                call writedynamicfunction(nFreq_Re,Debug_Lat_CorrFn_Re,'G_LatSchmidt_Re',tag=iter,    &
                    tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.)
                call CalcLatticeSpectrum(iCorrFnTag,nFreq_Re,CorrFn_Re,GFChemPot,tMatbrAxis=.false.,    &
                    iLatParams=iLatParams,LatParams=LatParams)
                call writedynamicfunction(nFreq_Re,CorrFn_Re,'G_Lat_Re',tag=iter,   &
                    tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.)
                !call CheckGFsSame(nFreq_Re,Debug_Lat_CorrFn_Re,CorrFn_Re,1.0e-7_dp)
            endif

            if(iLatticeFitType.eq.2) then
                !Invert HL Corr Fn, since we are fitting the residual of the inverses
                call InvertLocalNonHermFunc(nFitPoints,CorrFn_HL)
            endif

            if((abs(Damping_SE-one).gt.1.0e-8_dp).and.(iter.gt.1)) then
                !Take an admixture of the previous two high-level calculations to damp the fit
                CorrFn_HL(:,:,:) = (Damping_SE*CorrFn_HL(:,:,:)) + ((one-Damping_SE)*CorrFn_HL_Old(:,:,:))
            endif

            if(iter.eq.1) then
                !calculate the initial residual
                call CalcLatticeFitResidual_2(iCorrFnTag,CorrFn_HL,nFitPoints,GFChemPot,iLatParams,LatParams,   &
                    AllDiffs(1,0),tFitMatAxis,FreqPoints=FreqPoints,Weights=Weights)
                write(6,"(A,F20.10)") "Initial spectrum residual: ",AllDiffs(1,0)
                AllDiffs(2,0) = zero
            endif

            call FitLatParams(iCorrFnTag,CorrFn_HL,nFitPoints,GFChemPot,iLatParams, &
                LatParams,FinalDist,tFitMatAxis,FreqPoints,Weights,iter)

            !What is the lattice hamiltonian
            call LatParams_to_ham(iLatParams,LatParams,GFChemPot,h_lat_fit)

            !Now, calculate the local lattice greens function
            CorrFn_Fit_Old(:,:,:) = CorrFn_Fit(:,:,:)
            call CalcLatticeSpectrum(iCorrFnTag,nFitPoints,CorrFn_Fit,GFChemPot,tMatbrAxis=tFitMatAxis, &
                iLatParams=iLatParams,LatParams=LatParams,FreqPoints=FreqPoints)

            !Write out the lattice couplings
            call WriteLatticeParams(iLatParams,LatParams)

            !Calculate & write out stats (all of them)
            AllDiffs(1,iter) = FinalDist
            DiffImpCorrFn(:,:,:) = CorrFn_HL_Old(:,:,:) - CorrFn_HL(:,:,:)
            !write(6,*) DiffImpCorrFn(:,:,:)
            DiffImpCorrFn(:,:,:) = DiffImpCorrFn(:,:,:) * dconjg(DiffImpCorrFn(:,:,:))
            AllDiffs(2,iter) = sum(real(DiffImpCorrFn(:,:,:),dp))
            !And diff for LatticeGF
            DiffImpCorrFn(:,:,:) = CorrFn_Fit_Old(:,:,:) - CorrFn_Fit(:,:,:)
            DiffImpCorrFn(:,:,:) = DiffImpCorrFn(:,:,:) * dconjg(DiffImpCorrFn(:,:,:))
            AllDiffs(3,iter) = sum(real(DiffImpCorrFn(:,:,:),dp))

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
            do i = 1,iLatParams
                if(abs(LatParams(i)).gt.1000.0_dp) then
                    write(6,*) LatParams(:)
                    call stop_all(t_r,'Lattice eigenvalues diverging. Convergence failed.')
                endif
            enddo
        enddo

        !We should write these out to a file so that we can restart from them at a later date (possibly with more variables).
        !TODO: Actually write out k-blocks
        write(6,"(A)") "Final *frequency independent* lattice parameters: "
        do i = 1,iLatParams
            write(6,"(I6,G20.13)") i,LatParams(i)
        enddo
        write(6,"(A)") "" 

        !TODO
        !if(tOptGF_EVals.and.tDiag_kSpace) then
        !    !Write out the converged one-electron dispersion / bandstructure
        !    call WriteBandstructure(Couplings,iLatParams)
        !endif
        call writedynamicfunction(nFitPoints,CorrFn_Fit,'G_Lat_Fit_Final',      &
            tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)

        deallocate(DiffImpCorrFn,AllDiffs,CorrFn_Fit_Old,CorrFn_Fit,CorrFn_HL_Old,dummy_Im,Debug_Lat_CorrFn_Fit)
        if(tCalcRealSpectrum) deallocate(dummy_Re,CorrFn_HL_Re,CorrFn_Re,Debug_Lat_CorrFn_Re)

        if(tFitMatAxis) then
            allocate(CorrFn_Fit(nImp,nImp,nFreq_Re))
            allocate(dummy_Re(nImp,nImp,nFreq_Re))
            allocate(CorrFn_HL_Re(nImp,nImp,nFreq_Re))
            dummy_Re = zzero
                
            write(6,"(A)") "Now calculating spectral functions on the REAL axis"
            call flush(6)
            !What does the final lattice greens function look like on the real axis?
            call CalcLatticeSpectrum(iCorrFnTag,nFreq_Re,CorrFn_Fit,GFChemPot,tMatbrAxis=.false.,iLatParams=iLatParams, &
                LatParams=LatParams)
            !Write out lattice greens function
            !Is it the same as the impurity greens function on the real axis. This would be nice.
            call writedynamicfunction(nFreq_Re,CorrFn_Fit,'G_Lat_Re_Final',     &
                tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.)

            !Finally, calculate the greens function in real-frequency space.
            call SchmidtGF_wSE(CorrFn_HL_Re,GFChemPot,dummy_Re,nFreq_Re,tMatbrAxis=.false.,cham=h_lat_fit)
            call writedynamicfunction(nFreq_Re,CorrFn_HL_Re,'G_Imp_Re_Final',   &
                tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.)
            deallocate(CorrFn_HL_Re,CorrFn_Fit,dummy_Re)
        else
            call writedynamicfunction(nFreq_Re,CorrFn_HL,'G_Imp_Fit_Final',     &
                tCheckCausal=.true.,tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=tFitMatAxis,FreqPoints=FreqPoints)
        endif

        deallocate(CorrFn_HL)

        deallocate(h_lat_fit)
        deallocate(FreqPoints,Weights)

    end subroutine SC_Spectrum_Opt
        
    subroutine CalcLatticeSpectrum(iCorrFn,n,CorrFn,mu,tMatbrAxis,iLatParams,LatParams,FreqPoints,ham,SE)
        implicit none
        integer, intent(in) :: iCorrFn,n
        complex(dp), intent(out) :: CorrFn(nImp,nImp,n)
        real(dp), intent(in) :: mu
        logical, intent(in), optional :: tMatbrAxis
        integer, intent(in), optional :: iLatParams
        real(dp), intent(in), optional :: LatParams(:)
        real(dp), intent(in), optional :: FreqPoints(n)
        complex(dp), intent(in), optional :: ham(nSites,nSites)
        complex(dp), intent(in), optional :: SE(nImp,nImp,n)

        integer :: i,j,k,ind_1,ind_2!,ii,jj
        real(dp) :: Omega
!        complex(dp) :: compval
        complex(dp), allocatable :: KBlocks(:,:,:)!,EVecs(:,:,:),EVals(:)
        complex(dp) :: InvMat(nImp,nImp),InvMat2(nImp,nImp),num(nImp,nImp),GFContrib(nImp,nImp)!,hamtmp(nSites,nSites)
        logical :: tMatbrAxis_
        character(len=*), parameter :: t_r='CalcLatticeSpectrum'

        if((.not.present(LatParams)).and.(.not.present(ham))) then
            call stop_all(t_r,'No hamiltonian in real or lattice parameters passed in')
        endif
        if(present(LatParams).and.(present(ham))) then
            call stop_all(t_r,'Lattice parameters and hamiltonian present - which one to use?!')
        endif
        if(present(LatParams).and.(.not.present(iLatParams))) then
            call stop_all(t_r,'Lattice parameters sent in, but not suze')
        endif
        if(present(SE)) then
            call stop_all(t_r,'Sorry - cannot currently deal with self-energy')
        endif
        if(present(tMatbrAxis)) then
            tMatbrAxis_=tMatbrAxis
        else
            tMatbrAxis_=.false.
        endif

        allocate(KBlocks(nImp,nImp,nKPnts))

        CorrFn = zzero
        KBlocks = zzero

        if(present(LatParams)) then
            call LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)
        elseif(present(ham)) then
            call ham_to_KBlocks(ham,KBlocks)
        else
            call stop_all(t_r,'No hamiltonian (r or k space) read in')
        endif

!        call LatParams_to_ham(iLatParams,LatParams,mu,hamtmp)
!        allocate(EVecs(nImp,nImp,nKPnts))
!        allocate(EVals(nSites))
!        call KBlocks_to_diag(KBlocks,EVecs,EVals)
!        call writematrixcomp(hamtmp,'real space ham',.true.)
!        call writevectorcomp(EVals,'EVals')
!        deallocate(EVecs,EVals)
!        hamtmp = zzero
!        do i = 1,nKPnts
!            ind_1 = ((i-1)*nImp) + 1
!            ind_2 = nImp*i
!            hamtmp(ind_1:ind_2,ind_1:ind_2) = KBlocks(:,:,i)
!        enddo
!        compval = zzero
!        do i = 1,nSites
!            do j = 1,nSites
!                compval = compval + RtoK_Rot(1,i)*hamtmp(i,j)*dconjg(RtoK_Rot(2,j))
!            enddo
!        enddo
!        write(6,*) "Element 1,2 of ham: ",compval

        i = 0
        do while(.true.)
            if(present(FreqPoints)) then
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_,nFreqPoints=n,FreqPoints=FreqPoints)
            else
                call GetNextOmega(Omega,i,tMatbrAxis=tMatbrAxis_)
            endif
            if(i.lt.0) exit
            if(i.gt.n) call stop_all(t_r,'Too many freq points')

            do k = 1,nKPnts 

                InvMat(:,:) = - KBlocks(:,:,k)

!                write(6,"(A,I6,A,2G25.10)") "For kpoint: ",k," h Off diagonal matrix element: ",InvMat(1,2)
!                write(6,"(A,I6,A,G25.10)") "For kpoint: ",k," h Off diagonal hermiticity: ",abs(InvMat(1,2)-dconjg(InvMat(2,1)))
                if(iCorrFn.eq.1) then
                    !Greens function
                    do j = 1,nImp
                        if(tMatbrAxis_) then
                            InvMat(j,j) = InvMat(j,j) + cmplx(mu,Omega,dp)
                        else
                            !To get the off-diagonals to be hermitian, we have to be careful with the sign of the broadening.
                            !Do not worry about this for the moment, because for the diagonals it should be fine, and we are not fitting the real spectrum atm.
!                            if(Omega.gt.zero) then
                                InvMat(j,j) = InvMat(j,j) + cmplx(Omega + mu,dDelta,dp)
!                            else
!                                InvMat(j,j) = InvMat(j,j) + cmplx(Omega + mu,-dDelta,dp)
!                            endif
                        endif
                    enddo
                    call mat_inv(InvMat,InvMat2)
                
!                    call writematrixcomp(InvMat2,'Inverse matrix from mat_inv',.true.)
!                    write(6,"(A,I6,A,G25.10)") "For kpoint: ",k," (w-h)^-1 Off diagonal hermiticity: ",abs(InvMat2(1,2)-dconjg(InvMat2(2,1)))
!                    do ii = 1,nImp
!                        do jj = 1,nImp
!                            if(ii.eq.jj) then
!                                cycle
!                            else
!                                if(abs(aimag(InvMat2(ii,jj))+aimag(InvMat2(jj,ii))).gt.1.0e-6_dp) then
!                                    call writematrixcomp(InvMat,'Mat',.true.)
!                                    call writematrixcomp(InvMat2,'InvMat',.true.)
!                                    write(6,*) "Omega: ",Omega
!                                    write(6,*) "k: ",k
!                                    write(6,*) "Error is: ",abs(aimag(InvMat2(ii,jj))+aimag(InvMat2(jj,ii)))
!                                endif
!                            endif
!                        enddo
!                    enddo
                    ind_1 = ((k-1)*nImp) + 1
                    ind_2 = nImp*k

        !            C * C^+ here?  For 1 imp, should be |RtoK_Rot(1,k)|^2
!                    call ZGEMM('N','C',nImp,nImp,nImp,zone,RtoK_Rot(1:nImp,ind_1:ind_2),nImp,   &
!                        RtoK_Rot(1:nImp,ind_1:ind_2),nImp,zzero,num,nImp)
!                    CorrFn(:,:,i) = CorrFn(:,:,i) + ( InvMat2(:,:) * num(:,:) )
                    !Use InvMat as a temporary storage
                    InvMat(:,:) = RtoK_Rot(1:nImp,ind_1:ind_2)
                    call ZGEMM('N','C',nImp,nImp,nImp,zone,InvMat2,nImp,InvMat,nImp,    &
                        zzero,num,nImp)
                    call ZGEMM('N','N',nImp,nImp,nImp,zone,InvMat,nImp,num,nImp,  &
                        zzero,GFContrib,nImp)
                    
!                    write(6,"(A,I6,A,2G25.10)") "For kpoint: ",k," G real space Off diagonal matrix element: ",GFContrib(1,2)
!                    write(6,"(A,I6,A,G25.10)") "For kpoint: ",k," G real space Off diagonal hermiticity: ",abs(GFContrib(1,2)-dconjg(GFContrib(2,1)))
                    CorrFn(:,:,i) = CorrFn(:,:,i) + GFContrib(:,:)
                        
                else
                    call stop_all(t_r,'Cannot deal with non-greens functions right now')
                endif
            enddo
        enddo

!        if(iCorrFn.eq.1) then
            !num term Equivalent to below for one impurity
            !CorrFn(:,:,:) = CorrFn(:,:,:) / real(nKPnts,dp)
!        endif

        deallocate(KBlocks)

    end subroutine CalcLatticeSpectrum
        
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
    subroutine CalcLatticeFitResidual_2(iCorrFnTag,G_Imp,nESteps,mu,iLatParams,LatParams,dist,tMatbrAxis,     &
            FreqPoints,Weights,dJacobian)
        implicit none
        integer, intent(in) :: nESteps,iLatParams,iCorrFnTag
        real(dp), intent(in) :: mu
        complex(dp), intent(in) :: G_Imp(nImp,nImp,nESteps)
        real(dp), intent(in) :: LatParams(iLatParams)
        real(dp), intent(out) :: dist
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in), optional :: FreqPoints(nESteps)
        real(dp), intent(in), optional :: Weights(nESteps)
        real(dp), intent(out), optional :: dJacobian(iLatParams)
        complex(dp), allocatable :: Lattice_Func(:,:,:),DiffMat(:,:,:)
        real(dp), allocatable :: DiffMatr(:,:,:)
        real(dp) :: Omega,LattWeight
        logical :: tNonStandardGrid
        integer :: i,j,k
        integer, parameter :: iNormPower = 2    !The power of the matrix norm for the residual
        logical, parameter :: tTestDerivs = .false.
        character(len=*), parameter :: t_r='CalcLatticeFitResidual_2'

        !TODO: Fix this, so that the self-energy is an optional argument
        allocate(Lattice_Func(nImp,nImp,nESteps))

        if(present(FreqPoints)) then
            tNonStandardGrid = .true.
            call CalcLatticeSpectrum(iCorrFnTag,nESteps,Lattice_Func,mu,tMatbrAxis=tMatbrAxis,iLatParams=iLatParams, &
                LatParams=LatParams,FreqPoints=FreqPoints)
        else
            tNonStandardGrid = .false.
            call CalcLatticeSpectrum(iCorrFnTag,nESteps,Lattice_Func,mu,tMatbrAxis=tMatbrAxis,iLatParams=iLatParams, &
                LatParams=LatParams,FreqPoints=FreqPoints)
        endif

        if(iLatticeFitType.eq.2) then
            !If required (i.e. we are fitting the inverses of the functions), invert the lattice greens function
            call InvertLocalNonHermFunc(nESteps,Lattice_Func)
        endif
        
        allocate(DiffMat(nImp,nImp,nESteps))
        allocate(DiffMatr(nImp,nImp,nESteps))

        !Now, take the difference between the functions
        DiffMat(:,:,:) = Lattice_Func(:,:,:) - G_Imp(:,:,:)

        !Take abs values 
        DiffMatr(:,:,:) = abs(DiffMat(:,:,:))

        if(present(dJacobian)) then
            !Calculate the derivative of dist wrt each eigenvalue
            if(tNonStandardGrid) then
                call CalcJacobian_3(iCorrFnTag,nESteps,iLatParams,DiffMat,dJacobian,LatParams,mu,tMatbrAxis,     &
                    FreqPoints=FreqPoints,Weights=Weights)
            else
                call CalcJacobian_3(iCorrFnTag,nESteps,iLatParams,DiffMat,dJacobian,LatParams,mu,tMatbrAxis)
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
                    if(tDiagonalSC.and.(j.ne.k)) cycle  !In the diagonal approximation, the off diagonal GFs are not included
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
                        LattWeight = LattWeight + (abs(Lattice_Func(k,j,i)))**iNormPower
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
            
        if(present(dJacobian).and.tTestDerivs) then
            if(present(FreqPoints)) then
                call TestGradients(dist,iCorrFnTag,nESteps,iLatParams,LatParams,mu,G_Imp,tMatbrAxis,    &
                    dJacobian,FreqPoints=FreqPoints,Weights=Weights)
            else
                call TestGradients(dist,iCorrFnTag,nESteps,iLatParams,LatParams,mu,G_Imp,tMatbrAxis,dJacobian)
            endif
        endif

        deallocate(DiffMat,Lattice_Func,DiffMatr)

    end subroutine CalcLatticeFitResidual_2
        
    subroutine TestGradients(dist,iCorrFnTag,nESteps,iLatParams,LatParams,mu,G_Imp,tMatbrAxis,dJac,FreqPoints,Weights)
        implicit none
        real(dp), intent(in) :: dist
        integer, intent(in) :: iCorrFnTag,nESteps,iLatParams
        real(dp), intent(in) :: LatParams(iLatParams),mu
        complex(dp), intent(in) :: G_Imp(nImp,nImp,nESteps)
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in) :: dJac(iLatParams)
        real(dp), intent(in), optional :: FreqPoints(nESteps),Weights(nESteps)

        integer :: iunit,ivar
        real(dp) :: diff,NumDiff,Dist2
        real(dp), allocatable :: varsTemp(:)
        character(len=*), parameter :: t_r='TestGradients'

        iunit = 66
                
        write(6,*) "Testing derivatives of parameters: ",LatParams(:)
        
        allocate(varsTemp(iLatParams))
        do ivar = 1,iLatParams
            diff = 0.0001_dp
            do while(diff.gt.1.0e-16_dp)

                varsTemp(:) = LatParams(:)
                varsTemp(ivar) = varsTemp(ivar) + diff
                if(present(FreqPoints)) then
                    call CalcLatticeFitResidual_2(iCorrFnTag,G_Imp,nESteps,mu,iLatParams,varsTemp,dist2,    &
                        tMatbrAxis,FreqPoints=FreqPoints,Weights=Weights)
                else
                    call CalcLatticeFitResidual_2(iCorrFnTag,G_Imp,nESteps,mu,iLatParams,varsTemp,dist2,tMatbrAxis)
                endif
                NumDiff = (dist2 - dist) / diff
                write(iunit,"(I8,6G25.14)") ivar,diff,NumDiff,dJac(ivar),abs((NumDiff-dJac(ivar))/NumDiff), &
                    abs(NumDiff-dJac(ivar)),dist2
                if(abs((NumDiff-dJac(ivar))/NumDiff).gt.one) then
                    write(6,"(A,I8,6G25.14)") "***",ivar,diff,NumDiff,dJac(ivar),abs((NumDiff-dJac(ivar))/NumDiff), &
                        abs(NumDiff-dJac(ivar)),dist2
                endif
                diff = diff/2.0_dp
            enddo
            write(iunit,"(A)") ""
        enddo
        deallocate(varsTemp)
        call stop_all(t_r,'End of Gradients Test')

    end subroutine TestGradients

    !Use a minimization routine to fit the greens functions by adjusting the lattice coupling
    subroutine FitLatParams(iCorrFnTag,CorrFn_HL,n,mu,iLatParams,LatParams,FinalErr,tMatbrAxis,FreqPoints,Weights,iter)
        use MinAlgos
        use Levenberg_Marquardt
        use lbfgs
        use sort_mod, only: sort_real
        implicit none
        integer, intent(in) :: iCorrFnTag
        integer, intent(in) :: n    !Number of frequency points
        integer, intent(in) :: iLatParams    !Number of independent coupling parameters
        complex(dp), intent(in) :: CorrFn_HL(nImp,nImp,n)   !This may be the inverse
        real(dp), intent(inout) :: LatParams(iLatParams)    !The lattice couplings
        real(dp), intent(out) :: FinalErr
        real(dp), intent(in) :: mu
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in) :: FreqPoints(n),Weights(n)
        integer, intent(in) :: iter

        real(dp), allocatable :: step(:),vars(:),var(:),FinalVec(:),Jac(:,:),vars_temp(:)
        real(dp), allocatable :: Freqs_dum(:),Weights_dum(:),low(:),upp(:),grad(:),wa(:)
        integer, allocatable :: iwork(:),nbd(:),iwa(:)
        integer :: nop,i,maxf,iprint,nloop,iquad,ierr,nFuncs,j,ierr_tmp,corrs
        real(dp) :: stopcr,simp,rhobeg,rhoend,InitErr,dsave(29),pgtol,FinalErr2
        logical :: tfirst,tOptEVals_,tNonStandardGrid
        complex(dp), allocatable :: PreSymCorr(:,:,:)
        character(len=60) :: task,csave
        integer :: isave(44)
        logical :: lsave(4)
        character(len=*), parameter :: t_r='FitLatParams'

        !Initialize parameters
        allocate(vars(iLatParams))
        vars(:) = LatParams(:)

        if(iFitAlgo.eq.4) then
            !Use L-BFGS(-constrained) to fit
            if(.not.tAnalyticDerivs) call stop_all(t_r,'Need analytic derivatives if optimizing with BFGS')

            !We can specify upper and lower bounds on the solutions. Set these to be /pm1000 
            allocate(low(iLatParams))   !Lower limits
            allocate(upp(iLatParams))   !Upper limits
            allocate(nbd(iLatParams))   !Type of limit on each variable: 0 - unbounded, 1 - lower only, 2 - both bounds, 3 - upper only
            nbd(:) = 0
            low(:) = -1000.0_dp
            upp(:) = 1000.0_dp

            !rhoend = 1.0e12_dp  !Low accuracy
            !rhoend = 1.0e7_dp   !Mid accuracy
            rhoend = 10.0_dp    !High accuracy
            pgtol = dFitTol_SC   !Maximum gradient exit criterion. Default=1E-5

            allocate(grad(iLatParams))
            
            !corrs is the maximum number of variable metric corrections
            !used to define the limited memory matrix.
            corrs = 6

            !Working arrays
            allocate(wa(2*iLatParams*corrs + 5*iLatParams + 11*corrs*corrs + 8*corrs))
            allocate(iwa(3*iLatParams))

            iprint = 25

            task='START'
            do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START')

                call setulb(iLatParams,corrs,vars,low,upp,nbd,FinalErr,grad,rhoend,   &
                    pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)

                if(task(1:2).eq.'FG') then
                    !Request the function (FinalErr), and the gradient (grad) at the current parameters (vars)

                    !Compute function and gradient
                    call CalcLatticeFitResidual_2(iCorrFnTag,CorrFn_HL,n,mu,iLatParams,vars,FinalErr,  &
                        tMatbrAxis,FreqPoints=FreqPoints,Weights=Weights,dJacobian=grad)
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
            deallocate(low,upp,nbd,grad,wa,iwa) 
                
        elseif(iFitAlgo.eq.3) then
            !Use a Levenberg-Marquardt algorithm for non-linear optimization with either finite-difference
            !or analytic derivatives
            call stop_all(t_r,'Levenberg-Marquardt not working. Not coded up for nImp > 1')
        elseif(iFitAlgo.eq.2) then
            !Use modified Powell algorithm for optimization (based on fitting to quadratics)
            write(6,"(A)") "Optimizing lattice hamiltonian with modified Powell algorithm"
            rhobeg = 0.05_dp
            rhoend = 1.0e-6_dp
            iprint = 2
            if(iMaxFitMicroIter.eq.0) then
                maxf = 20*iLatParams
            else
                maxf = iMaxFitMicroIter
            endif

            !Set some globals, so we don't need to transfer them through
            SCF_mu = mu
            SCF_iCorrFn = iCorrFnTag

            !Vars is updated to be the best value to minimize the residual
            call uobyqa(iLatParams,vars,rhobeg,rhoend,iprint,maxf,FinalErr,Min_Interface,n,CorrFn_HL,   &
                tMatbrAxis,.true., FreqPoints, Weights)
            write(6,"(A,G20.10)") "Residual after fit: ",FinalErr

        elseif(iFitAlgo.eq.1) then
            !Use simplex method for optimization without derivatives
            write(6,"(A)") "Optimizing lattice hamiltonian with Simplex algorithm"

            !Starting values are assumed to be read in
            allocate(step(iLatParams))
            allocate(var(iLatParams))
            step(:) = 0.05_dp
            nop = iLatParams

            !Set max no of function evaluations. Default = 20*variables, print every 25
            if(iMaxFitMicroIter.eq.0) then
                maxf = 20*iLatParams
            else
                maxf = iMaxFitMicroIter
            endif
            iprint = 25

            !Set value for stopping criterion. Stopping occurs when the 
            !standard deviation of the values of the objective function at
            !the points of the current simplex < stopcr
            stopcr = 1.0e-4_dp
            nloop = 2*iLatParams !This is how often the stopping criterion is checked

            !Fit a quadratic surface to be sure a minimum has been found
            iquad = 1

            !As function value is being evaluated in real(dp), it should
            !be accurate to about 15 dp. If we set simp = 1.d-6, we should
            !get about 9dp accuracy in fitting the surface
            simp = 1.0e-6_dp
            
            !Set some globals, so we don't need to transfer them through
            SCF_mu = mu
            SCF_iCorrFn = iCorrFnTag

            !Now call minim to do the work
            tfirst = .true.
            do while(.true.)
                call minim(vars, step, iLatParams, FinalErr, maxf, iprint, stopcr, nloop, &
                    iquad, simp, var, Min_Interface, ierr, n, CorrFn_HL, tMatbrAxis, .true., FreqPoints,Weights)

                if(ierr.eq.0) exit
                if(.not.tFirst) exit
                tFirst = .false.    !We have found a minimum, but run again with a small number of iterations to check it is stable
                maxf = 3*iLatParams
            enddo

            !On output, Err is the minimized objective function, var contains the diagonal of the inverse of the information matrix, whatever that is
            if(ierr.eq.0) then
                write(6,"(A)") "Simplex optimization successful."
                write(6,"(A,F14.7)") "Minimized residual: ",FinalErr
            elseif(ierr.eq.4) then
                call stop_all(t_r,'nloop < 1')
            elseif(ierr.eq.3) then
                call stop_all(t_r,'iLatParams < 1')
            elseif(ierr.eq.2) then
                call warning(t_r,'Information matrix is not +ve semi-definite')
                write(6,"(A,F14.7)") "Final residual: ",FinalErr
            elseif(ierr.eq.1) then
                call warning(t_r,'Max number of Simplex function evaluations reached.')
                write(6,"(A,F14.7)") "Final residual: ",FinalErr
            endif

            deallocate(step,var)
        endif
            
        if(tImposephsym.or.tImposeksym.or.tConstrainphsym) then
!            allocate(PreSymCorr(nImp,nImp,n))
!            call CalcLatticeSpectrum(iCorrFnTag,n,PreSymCorr,mu,tMatbrAxis=tMatbrAxis,    &
!                iLatParams=iLatParams,LatParams=vars,FreqPoints=FreqPoints)
!            call writedynamicfunction(n,PreSymCorr,'G_Lat_Presym',tag=iter+1,tMatbrAxis=tMatbrAxis,FreqPoints=FreqPoints)
!            deallocate(PreSymCorr)
            if(tImposephsym.or.tImposeksym) then
                !Impose momentum inversion, and potentially ph symmetry
                !Symmetrize the resulting eigenvalues
                call Imposesym_2(iLatParams,vars,mu)
            elseif(tConstrainphsym) then
                call stop_all(t_r,'Constraining ph symmetry not currently working. Fix me?')
                !Eigenvalues may have drifted to be occupied orbitals. 
                !Flip them. This should not change the *local* lattice greens function at all, but will change the bandstructure
                !and the coupling to the impurity site. However, since this is not what we are fitting, we should not worry if we
                !change this?
                do i=1,iLatParams
                    if(vars(i).lt.mu) then
                        !It is below the chemical potential. Flip it to get it to become a virtual again
                        !This should not change anything, since a corresponding virtual will also be flipped
                        !write(6,*) "old e: ",vars(i)
                        vars(i) = 2.0_dp*mu - vars(i)
                        !write(6,*) "new e: ",vars(i)
                    endif
                enddo
                call sort_real(vars,iLatParams)
                !Actually, this has sorted things the wrong way around. Swap them
                allocate(vars_temp(iLatParams))
                vars_temp(:) = vars(:)
                do i = 1,iLatParams
                    vars(i) = vars_temp(iLatParams-i+1)
                enddo
                deallocate(vars_temp)
            endif

            !Calculate the new final residual
            call CalcLatticeFitResidual_2(iCorrFnTag,CorrFn_HL,n,mu,iLatParams,vars,FinalErr2,.true., &
                    FreqPoints=FreqPoints, Weights=Weights)
            write(6,"(A,G20.10)") "Final residual after symmetrization of lattice h: ",FinalErr2
            write(6,"(A,G20.10)") "Change in residual from symmetrizing variables: ",FinalErr2-FinalErr
            FinalErr = FinalErr2
            write(6,"(A)") "Lattice after symmetrization of parameters: "
            write(6,*) vars(:)
        endif
            
        !Update couplings
        LatParams(:) = vars(:)
        deallocate(vars)

    end subroutine FitLatParams

    !A wrapper for the simplex & powell routines
    !Two variables are taken from globals
    subroutine Min_Interface(LatParams,dist,nFreq,iLatParams,CorrFn_HL,tMatAxis,tGrid,Freqs,Weights,dJac)
        implicit none
        integer, intent(in) :: nFreq,iLatParams
        real(dp), intent(in) :: LatParams(iLatParams)
        real(dp), intent(out) :: dist
        complex(dp), intent(in) :: CorrFn_HL(nImp,nImp,nFreq)
        logical, intent(in) :: tMatAxis,tGrid
        real(dp), intent(in), optional :: Freqs(nFreq),Weights(nFreq)
        real(dp), intent(out), optional :: dJac(iLatParams)

        character(len=*), parameter :: t_r='Min_Interface'

        if(present(dJac)) then
            dJac(:) = 100.0_dp
            call stop_all(t_r,'Should not be here')
        endif
        if((.not.present(Freqs)).or.(.not.present(Weights))) call stop_all(t_r,'Should pass in freqs and weights')
        !write(6,*) "Input Vars: ",LatParams(:)
        call CalcLatticeFitResidual_2(SCF_iCorrFn,CorrFn_HL,nFreq,SCF_mu,iLatParams,LatParams,dist,  &
            tMatAxis,FreqPoints=FreqPoints,Weights=Weights)
        !write(6,*) "Output resid: ",dist

    end subroutine Min_Interface

    !Find the jacobian matrix for the optimization (in a few cases)
    subroutine CalcJacobian_3(iCorrFnTag,n,iLatParams,DiffMat,Jacobian,LatParams,mu,tMatbrAxis,FreqPoints,Weights)
        implicit none
        integer, intent(in) :: iCorrFnTag   !What spectrum are we considering
        integer, intent(in) :: n    !Number of frequency points
        integer, intent(in) :: iLatParams   !The number of variables
        complex(dp), intent(in) :: DiffMat(nImp,nImp,n) !The (complex) difference between the greens functions
        real(dp), intent(out) :: Jacobian(iLatParams) !The Jacobian
        real(dp), intent(in) :: LatParams(iLatParams)
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in) :: mu
        real(dp), intent(in), optional :: FreqPoints(n),Weights(n)
        
        complex(dp), allocatable :: KBlocks(:,:,:),FullJac_Re(:,:,:),FullJac_Im(:,:,:)
        complex(dp) :: ExtractMat(nImp,nImp),Mat(nImp,nImp),InvMat(nImp,nImp),ztmp(nImp,nImp),ztmp2(nImp,nImp)
        complex(dp) :: compval,realval,num(nImp,nImp)
        real(dp) :: Omega
        integer :: i,j,jj,k,w,ind,ks,ii,ind_1,ind_2
        logical :: tNonStandardGrid
        character(len=*), parameter :: t_r='CalcJacobian_3'

        if(iLatticeFitType.ne.1) call stop_all(t_r,'Cannot do gradients with non-linear objective functions yet')
        if(iCorrFnTag.ne.1) call stop_all(t_r,'Can currently only do greens functions')
        if(present(FreqPoints)) then
            tNonStandardGrid = .true.
            if(.not.(present(Weights))) then
                call stop_all(t_r,"No Weights present")
            endif
        else
            tNonStandardGrid = .false.
        endif

        !Expand the k-space hamiltonians
        allocate(KBlocks(nImp,nImp,nKPnts))
        call LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)

        allocate(FullJac_Re(nImp,nImp,nKPnts)) !Derivative of the cost function wrt. real parts of h
        FullJac_Re(:,:,:) = zzero  !Complex?
        allocate(FullJac_Im(nImp,nImp,nKPnts)) !Derivative of the cost function wrt. imag parts of h
        FullJac_Im(:,:,:) = zzero  !Complex?

        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k

!            C * C^+ here?  For 1 imp, should be |RtoK_Rot(1,k)|^2
!            call ZGEMM('N','C',nImp,nImp,nImp,zone,RtoK_Rot(1:nImp,ind_1:ind_2),nImp,   &
!                RtoK_Rot(1:nImp,ind_1:ind_2),nImp,zzero,num,nImp)
            num = RtoK_Rot(1:nImp,ind_1:ind_2)

            do i = 1,nImp
                do j = 1,nImp

                    !Pick real part
                    ExtractMat(:,:) = zzero
                    ExtractMat(j,i) = cmplx(one,zero,dp)

                    w = 0
                    do while(.true.)
                        if(tNonStandardGrid) then
                            call GetNextOmega(Omega,w,tMatbrAxis=tMatbrAxis,nFreqPoints=n,FreqPoints=FreqPoints)
                        else
                            call GetNextOmega(Omega,w,tMatbrAxis=tMatbrAxis)
                        endif
                        if(w.lt.0) exit

                        Mat(:,:) = - KBlocks(:,:,k)
                        do jj = 1,nImp
                            if(tMatbrAxis) then
                                Mat(jj,jj) = Mat(jj,jj) + cmplx(mu,Omega,dp)
                            else
                                Mat(jj,jj) = Mat(jj,jj) + cmplx(Omega + mu,dDelta,dp)
                            endif
                        enddo
                        InvMat = zzero
                        call mat_inv(Mat,InvMat)

                        call ZGEMM('N','N',nImp,nImp,nImp,zone,ExtractMat,nImp,InvMat,nImp,zzero,ztmp,nImp)
                        call ZGEMM('N','N',nImp,nImp,nImp,zone,InvMat,nImp,ztmp,nImp,zzero,ztmp2,nImp)

                        call ZGEMM('N','N',nImp,nImp,nImp,zone,num,nImp,ztmp2,nImp,zzero,ztmp,nImp)
                        call ZGEMM('N','C',nImp,nImp,nImp,zone,ztmp,nImp,num,nImp,zzero,ztmp2,nImp)

                        !Either divide by nKPnts, *or* multiply by num - they should be the same (only for 1 impurity?)
!                        ztmp2(:,:) = ztmp2(:,:) / real(nKPnts,dp)
!                        ztmp(:,:) = ztmp2(:,:) * dconjg(DiffMat(:,:,w)) * num(:,:)
                        ztmp(:,:) = ztmp2(:,:) * dconjg(DiffMat(:,:,w))

                        ztmp2(:,:) = ztmp(:,:) + dconjg(ztmp(:,:))

                        compval = zzero
                        do jj = 1,nImp
                            if(tDiagonalSC) then
                                compval = compval + ztmp2(jj,jj)
                            else
                                do ii = 1,nImp
                                    compval = compval + ztmp2(ii,jj)
                                enddo
                            endif
                        enddo
                        if(tNonStandardGrid) then
                            FullJac_Re(j,i,k) = FullJac_Re(j,i,k) + (compval * Weights(w))
                        else
                            FullJac_Re(j,i,k) = FullJac_Re(j,i,k) + compval
                        endif
                    enddo

                    if(i.ne.j) then
                        !Pick imaginary part
                        ExtractMat(:,:) = zzero
                        ExtractMat(j,i) = cmplx(zero,one,dp)

                        w = 0
                        do while(.true.)
                            if(tNonStandardGrid) then
                                call GetNextOmega(Omega,w,tMatbrAxis=tMatbrAxis,nFreqPoints=n,FreqPoints=FreqPoints)
                            else
                                call GetNextOmega(Omega,w,tMatbrAxis=tMatbrAxis)
                            endif
                            if(w.lt.0) exit

                            Mat(:,:) = - KBlocks(:,:,k)
                            do jj = 1,nImp
                                if(tMatbrAxis) then
                                    Mat(jj,jj) = Mat(jj,jj) + cmplx(mu,Omega,dp)
                                else
                                    Mat(jj,jj) = Mat(jj,jj) + cmplx(Omega + mu,dDelta,dp)
                                endif
                            enddo
                            InvMat = zzero
                            call mat_inv(Mat,InvMat)

                            call ZGEMM('N','N',nImp,nImp,nImp,zone,ExtractMat,nImp,InvMat,nImp,zzero,ztmp,nImp)
                            call ZGEMM('N','N',nImp,nImp,nImp,zone,InvMat,nImp,ztmp,nImp,zzero,ztmp2,nImp)
                        
                            call ZGEMM('N','N',nImp,nImp,nImp,zone,num,nImp,ztmp2,nImp,zzero,ztmp,nImp)
                            call ZGEMM('N','C',nImp,nImp,nImp,zone,ztmp,nImp,num,nImp,zzero,ztmp2,nImp)

    !                        ztmp2(:,:) = ztmp2(:,:) / real(nKPnts,dp)
                            !ztmp(:,:) = ztmp2(:,:) * dconjg(DiffMat(:,:,w)) * num(:,:)
                            ztmp(:,:) = ztmp2(:,:) * dconjg(DiffMat(:,:,w))

                            ztmp2(:,:) = ztmp(:,:) + dconjg(ztmp(:,:))

                            compval = zzero
                            do jj = 1,nImp
                                if(tDiagonalSC) then
                                    compval = compval + ztmp2(jj,jj)
                                else
                                    do ii = 1,nImp
                                        compval = compval + ztmp2(ii,jj)
                                    enddo
                                endif
                            enddo
                            if(tNonStandardGrid) then
                                FullJac_Im(j,i,k) = FullJac_Im(j,i,k) + (compval * Weights(w))
                            else
                                FullJac_Im(j,i,k) = FullJac_Im(j,i,k) + compval
                            endif
                        enddo
                    endif

                enddo
            enddo
        enddo
        deallocate(KBlocks)

        !Now, package up the jacobian back into the individual lattice parameters
        Jacobian(:) = zero  !This should be strictly real

        if(tConstrainKSym) then
            !e(k) = e(-k), and we are always using a uniform mesh
            !If gamma-centered mesh, then have nSites/2 + 1 independent parameters (we are sampling k=0 and BZ boundary which dont pair)
            !If Shifted mesh, then we have nSites/2 independent parameters
            if(tShift_Mesh) then
                ks = nKPnts/2
            else
                ks = (nKPnts/2) + 1
            endif
        else
            ks = nKPnts
        endif

        ind = 1
        do k = 1,ks
            if(tConstrainphsym) then
                do i = 1,nImp
                    do j = i,nImp-i+1
                        call stop_all(t_r,'Constraining ph sym not yet available for gradients')
                        if(i.eq.j) then
                        else
                        endif
                    enddo
                enddo
            else
                do i = 1,nImp
                    do j = i,nImp
                        if(i.eq.j) then
                            if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            realval = FullJac_Re(j,i,k)     !No transpose element to add
                            if(tConstrainKSym) then
                                if(tShift_Mesh) then
                                    realval = realval + FullJac_Re(j,i,nKPnts-k+1)
                                else
                                    if((k.ne.1).and.(k.ne.ks)) then
                                        realval = realval + FullJac_Re(j,i,nKPnts-k+2)
                                    endif
                                endif
                            endif
                            if(abs(aimag(realval)).gt.1.0e-7_dp) call stop_all(t_r,'This should be real')
                            Jacobian(ind) = real(realval,dp)
                            ind = ind + 1
                            !Completely ignore derivate wrt imaginary component of diagonal ham matrix element
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            realval = FullJac_Re(j,i,k) + FullJac_Re(i,j,k) !Add the transpose element
                            if(tConstrainKSym) then
                                !Add the other kpoint which matches
                                if(tShift_Mesh) then
                                    realval = realval + FullJac_Re(j,i,nKPnts-k+1) + FullJac_Re(i,j,nKPnts-k+1)
                                else
                                    if((k.ne.1).and.(k.ne.ks)) then
                                        realval = realval + FullJac_Re(j,i,nKPnts-k+2) + FullJac_Re(i,j,nKPnts-k+2)
                                    endif
                                endif
                            endif
                            if(abs(aimag(realval)).gt.1.0e-7_dp) call stop_all(t_r,'This should be real')
                            Jacobian(ind) = real(realval,dp)

                            !Remember that the transpose element will have opposite imaginary sign, therefore will move in the other direction?
                            realval = FullJac_Im(j,i,k) - FullJac_Im(i,j,k)
                            if(tConstrainKSym) then
                                !Add the other kpoint which matches
                                if(tShift_Mesh) then
                                    realval = realval - FullJac_Im(j,i,nKPnts-k+1) + FullJac_Im(i,j,nKPnts-k+1)
                                else
                                    if((k.ne.1).and.(k.ne.ks)) then
                                        realval = realval - FullJac_Im(j,i,nKPnts-k+2) + FullJac_Im(i,j,nKPnts-k+2)
                                    endif
                                endif
                            endif
                            if(abs(aimag(realval)).gt.1.0e-7_dp) call stop_all(t_r,'This should be real')
                            Jacobian(ind+1) = real(realval,dp)
                            ind = ind + 2
                        endif
                    enddo
                enddo
            endif
        enddo
        deallocate(FullJac_Re,FullJac_Im)

    end subroutine CalcJacobian_3
                
end module SelfConsistentLR2
