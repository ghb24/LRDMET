module SelfConsistentLR3
    use const
    use timing
    use errors
    use LinearResponse
    use globals
    use utils, only: get_free_unit,append_ext_real,append_ext
    use lattices, only: zero_localpot_comp
    use SC_Data                  
    use SelfConsistentUtils
    use matrixops, only: mat_inv
    use writedata
    implicit none

    contains

    subroutine SC_Spectrum_Static()
        implicit none
        real(dp) :: GFChemPot,NI_ChemPot
        integer :: i,j,nFreq,iter
        complex(dp), allocatable :: h_lat_fit(:,:)
        complex(dp), allocatable :: LatVecs(:,:),Lat_CorrFn(:,:,:),CorrFn_HL(:,:,:)
        complex(dp), allocatable :: Prev_CorrFn_HL(:,:,:),Prev_Lat_CorrFn(:,:,:)
        complex(dp), allocatable :: SE_Update(:,:,:),UpdatePotential(:,:),TotalPotential(:,:)
        complex(dp), allocatable :: ctemp(:,:),LatSelfEnergy_i(:,:),LatSelfEnergy_j(:,:)
        complex(dp), allocatable :: CorrFn_HL_Inv(:,:,:),HalfContract_i(:),HalfContract_j(:)
        complex(dp) :: Ei_ij_val,Ei_ji_val,Ej_ij_val,Ej_ji_val,ZDOTC,MeanImpDiag
        real(dp) :: MaxOffLocalEl
        real(dp), allocatable :: AllDiffs(:,:),LatVals(:) 
        integer, allocatable :: LatFreqs(:)
        real(dp), parameter :: dDeltaImpThresh = 1.0e-7_dp
        character(len=*), parameter :: t_r='SC_Spectrum_Static'

        call set_timer(SelfCon_LR)

        write(6,"(A)") "Entering quasiparticle self-consistent DMET..."

        !TODO:  Is potential local?
        !       How tight should convergence be?
        !       Damp potential, instead of/as well as self energy?
        !       Consistent gaps?
        !       Optimize for larger lattices
        !       Ground state energy from potential (via lat GF, imp GF and ground-state DMET)
        !       Get retarded greens function for the spectrum
        !       Update chemical potential correctly
        !       Constraints on correlation potential to aid convergence
        !           (real diags, hermitian, diagonals average to zero? Translationally invariant)

        !TODO:  Non-contracted bath space
        !       Fit results to Pade to remove broadening from self-energy, and
        !           ensure we don't need to include exactly the frequency points of
        !           the one-electron eigenvalues.
        !       Are there advantages to imaginary-frequency self-consistency?
        !       Thermal quantities with contracted ground state space
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
        tMinRes_NonDir = .false.
        tPrecond_MinRes = .false.
        tReuse_LS = .false.

        call SetChemPot(GFChemPot)
        write(6,"(A,G20.15)") "Chemical potential set to: ",GFChemPot
        if(tSC_StartwGSCorrPot) then
            write(6,"(A)") "Starting from correlation potential from ground-state calculation..."
        else
            write(6,"(A)") "Starting from bare lattice hamiltonian..."
        endif

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
!            if(.not.tSC_StartwGSCorrPot) h_lat_fit(i,i) = h_lat_fit(i,i) + GFChemPot    !Add chemical potential
        enddo

        !NI_ChemPot is fixed for the calculation as a shift for the grid
        NI_ChemPot = GFChemPot

        allocate(LatVals(nSites))
        allocate(LatVecs(nSites,nSites))
        LatVecs(:,:) = h_lat_fit(:,:)
        LatVals(:) = zero
        call DiagOneEOp(LatVecs,LatVals,nImp,nSites,tDiag_kspace,.false.)
        !Note that the chemical potential is *not* included in the definition of h 

        call writevector(LatVals,'LatVals')

        !LatFreqs tells you which frequency point corresponds to which
        !eigenvalue of the lattice
        allocate(LatFreqs(nSites))
        !nFreq should not change in this version of SetReFreqPoints
        call SetReFreqPoints(LatVals,NI_ChemPot,nFreq,LatFreqs)

        allocate(Lat_CorrFn(nImp,nImp,nFreq))
        allocate(CorrFn_HL(nImp,nImp,nFreq))
        allocate(CorrFn_HL_Inv(nImp,nImp,nFreq))

        allocate(Prev_CorrFn_HL(nImp,nImp,nFreq))
        allocate(Prev_Lat_CorrFn(nImp,nImp,nFreq))

        allocate(SE_Update(nImp,nImp,nFreq))
        Lat_CorrFn(:,:,:) = zzero ; CorrFn_HL(:,:,:) = zzero
        Prev_CorrFn_HL(:,:,:) = zzero ; Prev_Lat_CorrFn(:,:,:) = zzero
        SE_Update(:,:,:) = zzero

        allocate(UpdatePotential(nSites,nSites))
        allocate(TotalPotential(nSites,nSites))
        TotalPotential(:,:) = zzero

        !1: Change in correlation potential
        !2: Change in impurity greens function
        !3: Change in lattice greens function
        !4: Largest off-local part of correlation potential
        !5: Mean diagonal of local part of correlation potential
        allocate(AllDiffs(5,0:iMaxIter_MacroFit+1))
        AllDiffs(:,:) = zero
           
        if(.false.) then
            call CalcLatticeSpectrum(1,nFreq,Lat_CorrFn,GFChemPot,tMatbrAxis=.false.,    &
                Freqpoints=FreqPoints,ham=h_lat_fit)
                
            call writedynamicfunction(nFreq,Lat_CorrFn,'G_Lat',tag=iter,tCheckCausal=.false.,   &
                tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)
        endif

        iter = 0
        do while(.not.tSkip_Lattice_Fit)
            iter = iter + 1

            if(iter.ne.1) call SetReFreqPoints(LatVals,NI_ChemPot,nFreq,LatFreqs)
            call writevector(LatFreqs,'LatFreqs')
        
            call CalcLatticeSpectrum(1,nFreq,Lat_CorrFn,GFChemPot,tMatbrAxis=.false.,    &
                Freqpoints=FreqPoints,ham=h_lat_fit)
            
            call writedynamicfunction(nFreq,Lat_CorrFn,'G_Lat',tag=iter,tCheckCausal=.false.,   &
                tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)

            call SchmidtGF_FromLat(CorrFn_HL,GFChemPot,nFreq,tFitMatAxis,h_lat_fit,FreqPoints)

            call writedynamicfunction(nFreq,CorrFn_HL,'G_Imp',tag=iter,tCheckCausal=.false.,  &
                tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)
            
            !Now use Dysons equation: Sigma = G_lat^-1 - G_HL^-1
            SE_Update(:,:,:) = Lat_CorrFn(:,:,:)
            call InvertLocalNonHermFunc(nFreq,SE_Update)
            CorrFn_HL_Inv(:,:,:) = CorrFn_HL(:,:,:)
            call InvertLocalNonHermFunc(nFreq,CorrFn_HL_Inv)
            SE_Update(:,:,:) = Damping_SE*(SE_Update(:,:,:) - CorrFn_HL_Inv(:,:,:))

            if(.false.) then
                !Test - is the lattice greens function with the self energy
                !*exactly* the same as the HL greens function?
                call writedynamicfunction(nFreq,SE_Update,'SE_Imp',tag=iter,tCheckCausal=.false., &
                    tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)

                call CalcLatticeSpectrum(1,nFreq,CorrFn_HL_Inv,GFChemPot,tMatbrAxis=.false., &
                    FreqPoints=FreqPoints,ham=h_lat_fit,SE=SE_Update)

                call writedynamicfunction(nFreq,CorrFn_HL_Inv,'G_Lat_wSE',tag=iter,tCheckCausal=.false., &
                    tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)
            endif
            
            allocate(LatSelfEnergy_i(nSites,nSites))
            allocate(LatSelfEnergy_j(nSites,nSites))
            allocate(ctemp(nSites,nSites))
            allocate(HalfContract_i(nSites))
            allocate(HalfContract_j(nSites))
            UpdatePotential(:,:) = zzero
            !Use Quasi-particle self-consistency to find best static, hermitian
            !potential approximation to the value: IS THIS LOCAL?
            do i = 1,nSites
                
!                write(6,"(A,I7)") "Site: ",i
!                call flush(6)
                !Find the energy according to this orbital
                !This is the frequency of the FreqPoints(LatFreqs(i))
                if(abs(FreqPoints(LatFreqs(i))-(LatVals(i)-GFChemPot)).gt.1.0e-7_dp) then
                    write(6,*) i,LatFreqs(i),FreqPoints(LatFreqs(i)),LatVals(i)
                    call stop_all(t_r,'Error in assigning frequencies')
                endif

                !Stripe the self energy from the LatFreqs(i) local self energy
                !across the lattice
                LatSelfEnergy_i(:,:) = zzero
                call add_localpot_comp_inplace(LatSelfEnergy_i,SE_Update(:,:,LatFreqs(i)),.true.)

                !Rotate the ij and ji elements of LatSelfEnergy to the MO basis
!                call ZGEMV('C',nSites,nSites,zone,LatSelfEnergy_i,nSites,LatVecs(:,i),1,zzero,HalfContract_i,1)
                call ZGEMM('C','N',nSites,nSites,nSites,zone,LatVecs,nSites,LatSelfEnergy_i,nSites,zzero,ctemp,nSites)
                call ZGEMM('N','N',nSites,nSites,nSites,zone,ctemp,nSites,LatVecs,nSites,zzero,LatSelfEnergy_i,nSites)

                do j = 1,nSites
                    !Find the energy according to this orbital
                    !This is the frequency of the FreqPoints(LatFreqs(i))
                    if(abs(FreqPoints(LatFreqs(j))-(LatVals(j)-GFChemPot)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Error in assigning frequencies')
                    endif
                
                    !call ZGEMV('C',nSites,nSites,zone,LatSelfEnergy_i,nSites,LatVecs(:,j),1,zzero,HalfContract_j,1)
                    !Ei_ij_val = zdotc(nSites,HalfContract_i,1,LatVecs(:,j),1)
                    !Ei_ji_val = zdotc(nSites,HalfContract_j,1,LatVecs(:,i),1)

                    !Stripe the self energy from the LatFreqs(i) local self energy
                    !across the lattice
                    LatSelfEnergy_j(:,:) = zzero
                    call add_localpot_comp_inplace(LatSelfEnergy_j,SE_Update(:,:,LatFreqs(j)),.true.)

                    !Rotate LatSelfEnergy to the MO basis
                    call ZGEMM('C','N',nSites,nSites,nSites,zone,LatVecs,nSites,LatSelfEnergy_j,nSites,zzero,ctemp,nSites)
                    call ZGEMM('N','N',nSites,nSites,nSites,zone,ctemp,nSites,LatVecs,nSites,zzero,LatSelfEnergy_j,nSites)
!                    call ZGEMV('C',nSites,nSites,zone,LatSelfEnergy_j,nSites,LatVecs(:,i),1,zzero,HalfContract_i,1)
!                    call ZGEMV('C',nSites,nSites,zone,LatSelfEnergy_j,nSites,LatVecs(:,j),1,zzero,HalfContract_j,1)
!                    Ej_ij_val = zdotc(nSites,HalfContract_i,1,LatVecs(:,j),1)
!                    Ej_ji_val = zdotc(nSites,HalfContract_j,1,LatVecs(:,i),1)
                    
                    !Add to the update potential
                    UpdatePotential(i,j) = 0.25_dp * (LatSelfEnergy_i(i,j) + conjg(LatSelfEnergy_i(j,i))    &
                        + LatSelfEnergy_j(i,j) + conjg(LatSelfEnergy_j(j,i)) )
!                    UpdatePotential(i,j) = 0.25_dp * (Ei_ij_val + conjg(Ei_ji_val)    &
!                        + Ej_ij_val + conjg(Ej_ji_val) )
                enddo
            enddo

            !Rotate the new static potential (Update Potential) to the AO basis
            call ZGEMM('N','N',nSites,nSites,nSites,zone,LatVecs,nSites,UpdatePotential,nSites,zzero,ctemp,nSites)
            call ZGEMM('N','C',nSites,nSites,nSites,zone,ctemp,nSites,LatVecs,nSites,zzero,UpdatePotential,nSites)
            deallocate(LatSelfEnergy_i,LatSelfEnergy_j)
            deallocate(HalfContract_i,HalfContract_j)

            do i = 1,nSites
                do j = 1,nSites
                    if(abs(UpdatePotential(j,i)-conjg(UpdatePotential(i,j))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Update potential not hermitian')
                    endif
                enddo
            enddo

            call writematrix(UpdatePotential,'Update potential in the AO basis',.true.)

            !Add new potential to lattice hamiltonian
            h_lat_fit(:,:) = h_lat_fit(:,:) + UpdatePotential(:,:)
            TotalPotential(:,:) = TotalPotential(:,:) + UpdatePotential(:,:)
            
            !Is Update potential local?!
            ctemp(:,:) = TotalPotential(:,:)
            call zero_localpot_comp(ctemp)
            MaxOffLocalEl = zero
            do i = 1,nSites
                do j = 1,nSites
                    if(abs(ctemp(j,i)).gt.MaxOffLocalEl) MaxOffLocalEl = abs(ctemp(j,i))
                enddo
            enddo
            AllDiffs(4,iter) = MaxOffLocalEl 
            MeanImpDiag = TotalPotential(1,1)
            do i = 2,nImp
                MeanImpDiag = MeanImpDiag + TotalPotential(i,i)
            enddo
            AllDiffs(5,iter) = abs(MeanImpDiag / nImp )
            !write(6,*) "Largest off-local part of the correlation potential: ",MaxOffLocalEl
            deallocate(ctemp)

            !Rediagonalize
            LatVecs(:,:) = h_lat_fit(:,:)
            LatVals(:) = zero
            call DiagOneEOp(LatVecs,LatVals,nImp,nSites,tDiag_kspace,.false.)

            !What is change in update potential and G (Use SE_Update to store
            !differences)
            AllDiffs(1,iter) = sum(real(UpdatePotential(:,:)*dconjg(UpdatePotential(:,:)))) / real(nSites**2,dp)
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
            write(6,"(A)") "     Iter.  PotentialChange   Delta_GF_Imp(iw)    Delta_GF_Lat(iw)   MaxOff-local_v     MeanDiag_v"
            do i = 0,iter
                write(6,"(I7,5G20.13)") i,AllDiffs(1,i),AllDiffs(2,i),AllDiffs(3,i),AllDiffs(4,i),AllDiffs(5,i)
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

        call writedynamicfunction(nFreq,CorrFn_HL,'G_HL_Final',tCheckCausal=.false., &
            tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)
        
        call writedynamicfunction(nFreq,Lat_CorrFn,'G_Lat_Final',tCheckCausal=.false., &
            tCheckOffDiagHerm=.false.,tWarn=.true.,tMatbrAxis=.false.,FreqPoints=FreqPoints)

        !Calculate final self-energy

        !Calculate HL with self-energy as environment potential?

        !Write out G + self-energy for final bit (to match greens functions?)

        !Calculate the bandstructure

        deallocate(h_lat_fit,LatVals,LatVecs,Lat_CorrFn,CorrFn_HL,Prev_CorrFn_HL,Prev_Lat_CorrFn)
        deallocate(UpdatePotential,TotalPotential,AllDiffs,SE_Update,LatFreqs,CorrFn_HL_Inv)

        call halt_timer(SelfCon_LR)

    end subroutine SC_Spectrum_Static
                
end module SelfConsistentLR3
