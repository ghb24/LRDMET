module SelfConsistentLR3
    use const
    use timing
    use errors
    use LinearResponse
    use globals
    use utils, only: get_free_unit,append_ext_real,append_ext
    use SC_Data                  
    use SelfConsistentUtils
    use matrixops, only: mat_inv
    use writedata
    implicit none

    contains

    subroutine SC_Spectrum_Static()
        implicit none
        real(dp) :: GFChemPot
        integer :: i,j,nFreq,iter
        complex(dp), allocatable :: h_lat_fit(:,:)
        complex(dp), allocatable :: LatVecs(:,:),Lat_CorrFn(:,:,:),CorrFn_HL(:,:,:)
        complex(dp), allocatable :: Prev_CorrFn_HL(:,:,:),Prev_Lat_CorrFn(:,:,:)
        complex(dp), allocatable :: SE_Update(:,:,:),UpdatePotential(:,:),TotalPotential(:,:)
        complex(dp), allocatable :: ctemp(:,:),LatSelfEnergy_i(:,:),LatSelfEnergy_j(:,:)
        real(dp), allocatable :: AllDiffs(:,:),LatVals(:) 
        integer, allocatable :: LatFreqs(:)
        real(dp), parameter :: dDeltaImpThresh = 1.0e-4_dp
        character(len=*), parameter :: t_r='SC_Spectrum_Static'

        call set_timer(SelfCon_LR)

        !TODO:  Non-contracted bath space
        !       Fit results to Pade to remove broadening from self-energy, and
        !           ensure we don't need to include exactly the frequency points of
        !           the one-electron eigenvalues.
        !       Are there advantages to imaginary-frequency self-consistency?
        !       Thermal quantities
        !       Reintroduce frequency dependence into self-consistency (perhaps only into potential?)
        !       Get energy from greens function! (imp or lattice?)
        !       Care with the sign of the broadening
        !       G(iw) is antisymmetric - efficiency gain
        !       Self-consistent optimization of chemical potential?
        !       Include self-energy at the end?

        tFitPoints_Legendre = .false.
        tFitRealFreq = .true.
        tFitMatAxis = .false.
        tLR_ReoptGS = .true.
        tRemakeStaticBath = .true.
        tFullReoptGS = .false.
        tCompressedMats = .false.
        tMinRes_NonDir = .true.
        tPrecond_MinRes = .false.
        tReuse_LS = .false.

        call SetChemPot(GFChemPot)

        !Let h be complex to allow for easier integration with previous code
        allocate(h_lat_fit(nSites,nSites))
        h_lat_fit(:,:) = zzero
        !Do we want to start from a prior GS DMET calculation, or 
        do i = 1,nSites
            do j = 1,nSites
                if(tSC_StartwGSCorrPot) then
                    h_lat_fit(j,i) = cmplx(h0v(j,i),zero,dp)
                else
                    h_lat_fit(j,i) = cmplx(h0(j,i),zero,dp)
                endif 
            enddo
            if(.not.tSC_StartwGSCorrPot) h_lat_fit(i,i) = h_lat_fit(i,i) + GFChemPot    !Add chemical potential
        enddo

        allocate(LatVals(nSites))
        allocate(LatVecs(nSites,nSites))
        LatVecs(:,:) = h_lat_fit(:,:)
        LatVals(:) = zero
        call DiagOneEOp(LatVecs,LatVals,nImp,nSites,tDiag_kspace,.false.)

        !LatFreqs tells you which frequency point corresponds to which
        !eigenvalue of the lattice
        allocate(LatFreqs(nSites))
        !nFreq should not change in this version of SetReFreqPoints
        call SetReFreqPoints(LatVals,nFreq,LatFreqs)

        allocate(Lat_CorrFn(nImp,nImp,nFreq))
        allocate(CorrFn_HL(nImp,nImp,nFreq))

        allocate(Prev_CorrFn_HL(nImp,nImp,nFreq))
        allocate(Prev_Lat_CorrFn(nImp,nImp,nFreq))

        allocate(SE_Update(nImp,nImp,nFreq))
        Lat_CorrFn(:,:,:) = zzero ; CorrFn_HL(:,:,:) = zzero
        Prev_CorrFn_HL(:,:,:) = zzero ; Prev_Lat_CorrFn(:,:,:) = zzero
        SE_Update(:,:,:) = zzero

        allocate(UpdatePotential(nSites,nSites))
        allocate(TotalPotential(nSites,nSites))
        TotalPotential(:,:) = zzero

        allocate(AllDiffs(3,0:iMaxIter_MacroFit+1))
        AllDiffs(:,:) = zero

        iter = 0
        do while(.not.tSkip_Lattice_Fit)
            iter = iter + 1

            if(iter.ne.1) call SetReFreqPoints(LatVals,nFreq,LatFreqs)
        
            call CalcLatticeSpectrum(1,nFreq,Lat_CorrFn,GFChemPot,tMatbrAxis=.false.,    &
                Freqpoints=FreqPoints,ham=h_lat_fit)
            
            call writedynamicfunction(nFreq,Lat_CorrFn,'Initial_G_Lat',tCheckCausal=.true.,   &
                tCheckOffDiagHerm=.true.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)

            call SchmidtGF_FromLat(CorrFn_HL,GFChemPot,nFreq,tFitMatAxis,h_lat_fit,FreqPoints)

            call writedynamicfunction(nFreq,CorrFn_HL,'G_Imp',tag=iter,tCheckCausal=.true.,  &
                tCheckOffDiagHerm=.true.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)
            
            !Now use Dysons equation: Sigma = G_lat^-1 - G_HL^-1
            SE_Update(:,:,:) = Lat_CorrFn(:,:,:)
            call InvertLocalNonHermFunc(nFreq,SE_Update)
            call InvertLocalNonHermFunc(nFreq,CorrFn_HL)
            SE_Update(:,:,:) = SE_Update(:,:,:) - CorrFn_HL(:,:,:)
            
            allocate(LatSelfEnergy_i(nSites,nSites))
            allocate(LatSelfEnergy_j(nSites,nSites))
            allocate(ctemp(nSites,nSites))
            UpdatePotential(:,:) = zzero
            !Use Quasi-particle self-consistency to find best static, hermitian
            !potential approximation to the value: IS THIS LOCAL?
            do i = 1,nSites
                
                !Find the energy according to this orbital
                !This is the frequency of the FreqPoints(LatFreqs(i))
                if(abs(FreqPoints(LatFreqs(i))-LatVals(i)).gt.1.0e-7_dp) then
                    call stop_all(t_r,'Error in assigning frequencies')
                endif

                !Stripe the self energy from the LatFreqs(i) local self energy
                !across the lattice
                LatSelfEnergy_i(:,:) = zzero
                call add_localpot_comp_inplace(LatSelfEnergy_i,SE_Update(:,:,LatFreqs(i)),.true.)

                !Rotate LatSelfEnergy to the MO basis
                call ZGEMM('C','N',nSites,nSites,nSites,zone,LatVecs,nSites,LatSelfEnergy_i,nSites,zzero,ctemp,nSites)
                call ZGEMM('N','N',nSites,nSites,nSites,zone,ctemp,nSites,LatVecs,nSites,zzero,LatSelfEnergy_i,nSites)

                do j = 1,nSites
                    !Find the energy according to this orbital
                    !This is the frequency of the FreqPoints(LatFreqs(i))
                    if(abs(FreqPoints(LatFreqs(j))-LatVals(j)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Error in assigning frequencies')
                    endif

                    !Stripe the self energy from the LatFreqs(i) local self energy
                    !across the lattice
                    LatSelfEnergy_j(:,:) = zzero
                    call add_localpot_comp_inplace(LatSelfEnergy_j,SE_Update(:,:,LatFreqs(j)),.true.)

                    !Rotate LatSelfEnergy to the MO basis
                    call ZGEMM('C','N',nSites,nSites,nSites,zone,LatVecs,nSites,LatSelfEnergy_j,nSites,zzero,ctemp,nSites)
                    call ZGEMM('N','N',nSites,nSites,nSites,zone,ctemp,nSites,LatVecs,nSites,zzero,LatSelfEnergy_j,nSites)

                    !Add to the update potential
                    UpdatePotential(i,j) = 0.25_dp * (LatSelfEnergy_i(i,j) + conjg(LatSelfEnergy_i(j,i))    &
                        + LatSelfEnergy_j(i,j) + conjg(LatSelfEnergy_j(j,i)) )
                enddo
            enddo

            !Rotate the new static potential (Update Potential) to the AO basis
            call ZGEMM('N','N',nSites,nSites,nSites,zone,LatVecs,nSites,UpdatePotential,nSites,zzero,ctemp,nSites)
            call ZGEMM('N','C',nSites,nSites,nSites,zone,ctemp,nSites,LatVecs,nSites,zzero,UpdatePotential,nSites)
            deallocate(ctemp)
            deallocate(LatSelfEnergy_i,LatSelfEnergy_j)

            do i = 1,nSites
                do j = 1,nSites
                    if(abs(UpdatePotential(j,i)-conjg(UpdatePotential(i,j))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Update potential not hermitian')
                    endif
                enddo
            enddo

            !Is Update potential local?!
            call writematrix(UpdatePotential,'Update potential in the AO basis',.true.)

            !Add new potential to lattice hamiltonian
            h_lat_fit(:,:) = h_lat_fit(:,:) + UpdatePotential(:,:)
            TotalPotential(:,:) = TotalPotential(:,:) + UpdatePotential(:,:)

            !Rediagonalize
            LatVecs(:,:) = h_lat_fit(:,:)
            LatVals(:) = zero
            call DiagOneEOp(LatVecs,LatVals,nImp,nSites,tDiag_kspace,.false.)

            !What is change in update potential and G (Use SE_Update to store
            !differences)
            AllDiffs(1,iter) = sum(real(UpdatePotential(:,:)*dconjg(UpdatePotential(:,:)))) / real(nSites,dp)
            SE_Update(:,:,:) = Prev_CorrFn_HL(:,:,:) - CorrFn_HL(:,:,:)
            SE_Update(:,:,:) = SE_Update(:,:,:) * dconjg(SE_Update(:,:,:))
            AllDiffs(2,iter) = sum(real(SE_Update(:,:,:),dp))
            SE_Update(:,:,:) = Prev_Lat_CorrFn(:,:,:) - Lat_CorrFn(:,:,:)
            SE_Update(:,:,:) = SE_Update(:,:,:) * dconjg(SE_Update(:,:,:))
            AllDiffs(3,iter) = sum(real(SE_Update(:,:,:),dp))

            !Update previous correlation functions
            Prev_CorrFn_HL(:,:,:) = CorrFn_HL(:,:,:)
            Prev_Lat_CorrFn(:,:,:) = Lat_CorrFn(:,:,:)

            write(6,"(A)") ""
            write(6,"(A,I7,A)") "***   COMPLETED MACROITERATION ",iter," ***"
            write(6,"(A)") "     Iter.  PotentialChange   Delta_GF_Imp(iw)    Delta_GF_Lat(iw)"
            do i = 0,iter
                write(6,"(I7,3G20.13)") i,AllDiffs(1,i),AllDiffs(2,i),AllDiffs(3,i)
            enddo
            write(6,"(A)") ""
            call flush(6)

            if(iter.ge.iMaxIter_MacroFit) then
                write(6,"(A,I9)") "Exiting. Max iters hit of: ",iMaxIter_MacroFit
                exit
            endif

            if(AllDiffs(1,iter).lt.dDeltaImpThresh) then
                write(6,"(A)") "Success! Static potential converged"
                write(6,"(A,G20.13)") "Correlation potential changing by less than: ",dDeltaImpThresh
                exit
            endif

        enddo

        call writedynamicfunction(nFreq,CorrFn_HL,'G_HL_Final',tCheckCausal=.true., &
            tCheckOffDiagHerm=.true.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)
        
        call writedynamicfunction(nFreq,Lat_CorrFn,'G_Lat_Final',tCheckCausal=.true., &
            tCheckOffDiagHerm=.true.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)

        !Calculate final self-energy

        !Calculate HL with self-energy as environment potential?

        !Write out G + self-energy for final bit (to match greens functions?)

        !Calculate the bandstructure

        deallocate(h_lat_fit,LatVals,LatVecs,Lat_CorrFn,CorrFn_HL,Prev_CorrFn_HL,Prev_Lat_CorrFn)
        deallocate(UpdatePotential,TotalPotential,AllDiffs,SE_Update,LatFreqs)

        call halt_timer(SelfCon_LR)

    end subroutine SC_Spectrum_Static
                
end module SelfConsistentLR3
