Program RealHub
    use const
    use timing
    use errors, only: stop_all
    use globals
    use Solvers
    use SchmidtDecomp
    use LRDriver 
    use fitting
    use mat_tools
    use readinput
    use Lattices
    use report  
    implicit none

    call init_calc()
    call set_timer(full_timer)
    call run_DMETcalc()
    call halt_timer(full_timer)
    call end_calc()

    contains
    
    subroutine init_calc()
        implicit none

        write(6,"(A)") "***  Starting real-space hubbard/anderson calculation  ***"
        call init_timing()
        call name_timers()
        call environment_report()
        call set_defaults()
        call read_input()
        call check_input()
        call check_openmp()
        call Setup_Lattice()
        call write_output_header()

    end subroutine init_calc

    subroutine run_DMETcalc()
        implicit none
        real(dp) :: ElecPerSite,FillingFrac,VarVloc,ErrRDM,mean_vloc
        integer :: i,it,Occ,CurrU,DMETfile
        logical :: tFinishedU,t2RDM
        character(len=*), parameter :: t_r="run_DMETcalc"

        !Set up initial conditions, i.e. starting potential
        allocate(v_loc(nImp,nImp))
        allocate(h0(nSites,nSites))     !The core hamiltonian
        allocate(h0v(nSites,nSites))    !The core hamiltonian with the local potential
        allocate(HFEnergies(nSites))    !The fock energies
        v_loc(:,:) = zero 
        h0(:,:) = zero 
        h0v(:,:) = zero 
        HFEnergies(:) = zero 
        if(tUHF) then
            allocate(v_loc_b(nImp,nImp))
            allocate(h0_b(nSites,nSites))
            allocate(h0v_b(nSites,nSites))
            allocate(HFEnergies_b(nSites))
            v_loc_b(:,:) = zero 
            h0_b(:,:) = zero 
            h0v_b(:,:) = zero 
            HFEnergies_b(:) = zero 
        endif

        CurrU = 0
        do while(.true.)

            !Find the next U value
            call GetNextUVal(CurrU,tFinishedU)
            if(tFinishedU) exit
            if(.not.tSingFiss) write(6,*) "Running DMET calculation with U = ",U
        
            !Calculate the core hamiltonian based on the hopping matrix of the hubbard model in real space
            !If reading in the hopping matrix, it is done here and stored in h0
            call make_hop_mat()
        
            if((tDiag_KSpace.or.tProjectHFKPnts.or.tKSpaceOrbs).and.(.not.allocated(KPnts))) then
                call setup_kspace()
            endif

            !Diagonalize the mean-field hamiltonian
            !Get occupations with unique GS
            call find_occs()

            !Loop over occupation numbers 
            do Occ=1,N_Occs

                allocate(MeanFieldDM(nSites,nSites))    !DM from mean-field
                MeanFieldDM(:,:) = zero
                if(tUHF) then
                    allocate(MeanFieldDM_b(nSites,nSites))
                    MeanFieldDM_b(:,:) = zero
                endif

                call OpenDMETFile(DMETfile)
                write(DMETfile,"(A)") " #Iteration  E_DMET/Imp   E_HL   d[V]   Initial_Err[1RDM]   "    &
                    //"Filling   Filling_Err   mean_diag_correlation"

                !These occupations refer to number of closed shell orbitals, so total electrons is 2 x nOcc
                nOcc = allowed_occs(Occ)    !Number of occupied orbitals in this iteration
                NEl = 2*nOcc    !Total electrons in system
                ElecPerSite = (2.0_dp*nOcc)/real(nSites,dp) !Average number of electrons per site
                FillingFrac = ElecPerSite/2.0_dp    !A filling fraction of 1/2 is half filling (duh)

                !Write out some stats
                write(6,*) 
                write(6,"(A,F8.3,A,I5,A,I5,A)") "Electrons per site:   ",ElecPerSite," (in ", &
                    nOcc," doubly occupied orbitals on ",nSites," sites)"
                write(6,"(A,F10.5)")          "Filling Fraction:     ",FillingFrac
                write(6,"(A,F8.3)")           "Hubbard U:            ",U
                write(6,"(A,I5,A)")           "Embedded system size: ",nImp," sites"
                if(tAnderson) then
                    write(6,"(A,I5,A)")           "Anderson lattice of ",nSites," sites"
                else
                    write(6,"(A,I5,A)")           "1D Hubbard lattice of ",nSites," sites"
                endif
                write(6,*) 
                
                if(tDiagFullSystem) then
                    !Find all eigenvalues of system
                    call DiagFullSystem()
                endif

                !Calculate full hf, including mean-field on-site repulsion (which is included in correlation potential in DMET
                call set_timer(FullSCF)
                if(tSCFHF) then
                    call run_true_hf()
                else
                    call run_hf(0)
                endif
                call halt_timer(FullSCF)

                if(tDumpFCIDUMP) then
                    call set_timer(FCIDUMP)
                    call DumpFCIDUMP()
                    call halt_timer(FCIDUMP)
                endif

                !Calculate single reference linear response - non-interacting, TDA and RPA
                if(tMFResponse) then
                    call SR_LinearResponse()
                endif

                if(tReadInCorrPot.or.tReadSystem) then
                    !Read in the correlation potential from another source
                    call read_in_corrpot()
                elseif(tExactCorrPot) then
                    !We know the exact correlation potential for the hubbard model: U/2
                    v_loc(:,:) = U/2.0_dp
                endif

                !At this point, we have h0, U and a set of system sites (the first nImp indices), as well as a local potential
                do it=1,iMaxIterDMET

                    !Write out header file for the convergence
                    write(6,"(A,I6)") "Iteration: ",it

                    !Do iMaxIterDMET microiterations to converge the DMET for this occupation number
                    call add_localpot(h0,h0v,v_loc,core_b=h0_b,core_v_b=h0v_b,CorrPot_b=v_loc_b)

                    !Now run a HF calculation by constructing and diagonalizing the fock matrix
                    !This will also return the RDM in the AO basis
                    call set_timer(DiagT)
                    if(tReadSystem) then
                        call read_orbitals()
                    else
                        if(tT1SCF) then
                            if(it.eq.1) then
                                call run_hf(it)
                            else
                                call RotDet_T1(it)
                            endif
                        else
                            call run_hf(it)
                        endif
                    endif
                    call halt_timer(DiagT)

                    !Construct the embedded basis
                    call set_timer(ConstEmb)
                    if(tThermal) then
                        call ConstructThermalEmbedding()
                    else
                        if(tConstructFullSchmidtBasis) then
                            call ConstructFullSchmidtBasis()
                        else
                            call CalcEmbedding()
                        endif
                    endif
                    if(tWriteOut) then
                        call writematrix(EmbeddedBasis,'EmbeddedBasis',.true.)
                        if(tUHF) call writematrix(EmbeddedBasis_b,'BetaEmbeddedBasis',.true.)
                    endif
                    call halt_timer(ConstEmb)
                    
                    !Now transform the 1 electron quantities into the embedded basis
                    !This should be exactly correct, i.e. we can now diagonalize the fock matrix in this basis
                    !to get the same result. We could also check that the number of electrons adds to the correct thing
                    call set_timer(Trans1e)
                    call Transform1e()
                    call halt_timer(Trans1e)
                    
                    call set_timer(HL_Time)
                    !Construct the two electron integrals in the system, and solve embedded system with high-level method
                    t2RDM = .false.
                    if(tFCIQMC) t2RDM = .false.
                    call SolveSystem(t2RDM)
                    call halt_timer(HL_Time)

                    if((.not.tAnderson).and.tContinueConvergence) then
                        !Fit new potential
                        !vloc_change (global) is updated in here to reflect the optimal change
                        !VarVloc is a meansure of the change in the potential
                        !ErrRDM is a measure of the initial difference in the RDMs
                        call set_timer(Fit_v_time)
                        call Fit_vloc(VarVloc,ErrRDM)

                        if(tDebug) call writematrix(vloc_change,'vloc_change',.true.)

                        !Mean vloc is actually for the old vloc for consistency with Geralds code
                        mean_vloc = 0.0_dp
                        do i=1,nImp
                            mean_vloc = mean_vloc + v_loc(i,i)
                        enddo
                        mean_vloc = mean_vloc/real(nImp)
                        call halt_timer(Fit_v_time)

                        !Write out stats:
                        !   Iter    E/Site  d[V]    Initial_ERR[RDM]    ERR[Filling]    mean[corr_pot]      Some RDM stuff...?
                        write(6,"(I7,5G22.10)") it,TotalE_Imp,VarVloc,ErrRDM,FillingError,mean_vloc
                        write(DMETfile,"(I7,7G22.10)") it,TotalE_Imp,HL_Energy,VarVloc,ErrRDM,  &
                            Actualfilling_Imp,FillingError,mean_vloc

                        !Update vloc
                        v_loc(:,:) = v_loc(:,:) + vloc_change(:,:)

                        if(VarVloc.lt.dTolDMET) then
                            write(6,"(A)") "...correlation potential converged" 
                            exit
                        endif
                    elseif(tT1SCF) then
                        !A new self-consistency. Find the best single determinant analytically which matches the HL vector
                        call T1SCF()
                    else
                        !Write out stats:
                        !   Iter    E/Site  d[V]    ERR[RDM]    ERR[Filling]    mean[corr_pot]      Some RDM stuff...?
                        VarVloc = 0.0_dp
                        ErrRDM = 0.0_dp
                        mean_vloc = 0.0_dp

                        write(6,"(I7,5G22.10)") it,TotalE_Imp,VarVloc,ErrRDM,FillingError,mean_vloc
                        write(DMETfile,"(I7,7G22.10)") it,TotalE_Imp,HL_Energy,VarVloc,ErrRDM,  &
                            Actualfilling_Imp,FillingError,mean_vloc
                        call flush(6)

                        exit    !Anderson model, so we do not want to iterate
                    endif

                enddo   !DMET convergence

                if(it.gt.iMaxIterDMET) call warning(t_r,'DMET Convergence failed - try increasing MAXITER_DMET ?')
                    
                !Potentially run FCI again now to get correlation functions from 2RDMs?
                write(6,"(A,F10.4,A,G20.10)") "FINAL energy per site for U=",U,' is: ',TotalE_Imp
                call flush(6)
                
                if(.not.tAnderson) then
                    close(DMETfile)
                    if(tHalfFill) then
                        call WriteCorrPot()
                    endif
                    call writematrix(v_loc,'Converged Correlation Potential',.true.)
                endif
        
                deallocate(MeanFieldDM)
                if(tUHF) deallocate(MeanFieldDM_b)

                if(tProjectHFKPnts) then
                    call ProjectHFontoK()
                endif

                if(tKSpaceOrbs) then
                    call GetKSpaceOrbs()
                endif

                if(tCorrNI_Spectra) then
                    !Calculates single reference spectral functions using the correlated 1-electron hamiltonian
                    call Correlated_SR_LR()
                endif

                if(tLR_DMET) then
                    !Perform linear response on the resulting DMET state
                    call MR_LinearResponse()
                endif
                deallocate(HFOrbs)
                if(tUHF) deallocate(HFOrbs_b)

                !Set potential for the next occupation number, or wipe it?
                if(.not.tSaveCorrPot) then
                    v_loc(:,:) = zero
                    if(tUHF) v_loc_b(:,:) = zero
                endif

            enddo   !Loop over occupations

            if(.not.tHalfFill) then
                !Wipe correlation potential if we have ramped through occupations
                !We can potentially keep it though if we are just doing half-filling
                v_loc(:,:) = zero
                if(tUHF) v_loc_b(:,:) = zero
            endif

        enddo   !Loop over U values

    end subroutine run_DMETcalc

    !CurrU is zero on entry for the first time, and is returned as non-zero
    !tFinished = true when we have run through all U values.
    subroutine GetNextUVal(CurrU,tFinished)
        implicit none
        logical, intent(out) :: tFinished
        integer, intent(inout) :: CurrU
        character(len=*), parameter :: t_r='GetNextUVal'

        tFinished = .false.

        if((CurrU.ne.0).and.(tSingFiss)) then
            tFinished = .true.
            return !We don't want to loop through U!
        endif

        if(nU_Vals.eq.0) then
            !We are sweeping through U values, rather than specifying them individually
            if(CurrU.eq.0) then
                !First value
                U = StartU
                CurrU = 1   !So we don't go into this block again
            else
                !Carry on sweeping through.
                !Increment from the last value
                U = U + UStep
                !Find out if we are still in range
                if((U.gt.max(StartU,EndU)+1.0e-5_dp).or.(U.lt.min(StartU,EndU)-1.0e-5_dp)) then
                    !We have finished
                    tFinished = .true.
                endif
            endif
        else
            !We are running through specified U values
            if(.not.allocated(U_Vals)) call stop_all(t_r,'U_Vals array not allocated')
            CurrU = CurrU + 1
            if(CurrU.gt.nU_Vals) then
                !We have run through all U values we want
                tFinished = .true.
            else
                U = U_Vals(CurrU)
            endif
        endif

    end subroutine GetNextUVal

    subroutine WriteCorrPot()
        use utils, only: get_free_unit,append_ext_real
        implicit none
        integer :: iunit
        character(64) :: filename
!        character(len=*), parameter :: t_r='WriteCorrPot'

        write(6,*) "Writing out converged correlation potential..."

        call append_ext_real("CORRPOTS",U,filename)
        iunit = get_free_unit()
        open(iunit,file=filename,status='unknown')

        write(iunit,*) U,nOcc,v_loc(:,:)
        close(iunit)

    end subroutine WriteCorrPot

    subroutine read_in_corrpot()
        use utils, only: get_free_unit
        use DetTools, only: tospat
        implicit none
        integer :: iunit,Occ_val,ios,i,j,k
        logical :: texist,tFoundCorrPot
        real(dp) :: U_val,CorrPot_tmp(nImp*nImp)
        real(dp), allocatable :: temp(:)
        character(len=*), parameter :: t_r='read_in_corrpot'

        write(6,*) "Reading in correlation potential..."

        if(tReadSystem) then

            inquire(file='FinalVCorr.dat',exist=texist)
            if(.not.texist) call stop_all(t_r,'Correlation potential file cannot be found') 
            iunit = get_free_unit()
            open(iunit,file='FinalVCorr.dat',status='old',action='read')
            if(tUHF) then
                allocate(temp(nImp*2))
                do i = 1,nImp*2
                    read(iunit,*) temp(:)
                    do j = 1,nImp*2
                        if(mod(i,2).eq.1) then
                            !i is alpha spin-orbital
                            if((mod(j,2).eq.0).and.(abs(temp(j)).gt.1.0e-8_dp)) then
                                call stop_all(t_r,'Coupling between different spin types in correlation potential?!')
                            elseif(mod(j,2).eq.1) then
                                !j is also alpha
                                v_loc(tospat(i),tospat(j)) = temp(j)
                            endif
                        else
                            !i is beta spin-orbital
                            if((mod(j,2).eq.1).and.(abs(temp(j)).gt.1.0e-8)) then
                                call stop_all(t_r,'Coupling between different spin types in correlation potential?!')
                            elseif(mod(j,2).eq.0) then
                                !j is also beta
                                v_loc_b(tospat(i),tospat(j)) = temp(j)
                            endif
                        endif
                    enddo
                enddo
                deallocate(temp)
            else
                do i = 1,nImp
                    read(iunit,*) v_loc(i,:)
                enddo
            endif
            close(iunit)
        else
            iunit = get_free_unit()
            inquire(file=CorrPot_file,exist=texist)
            if(.not.texist) then
                write(6,*) "correlation potential filename: ",CorrPot_file
                call stop_all(t_r,'Expecting to read in a file with a converged '    &
     &              //'correlation potential, but unable to find appropriate file')
            endif

            open(iunit,file=CorrPot_file,status='old')

            tFoundCorrPot = .false.
            do while(.true.)
                read(iunit,*,iostat=ios) U_val,Occ_val,CorrPot_tmp(1:nImp*nImp)
                if(ios.gt.0) call stop_all(t_r,'Error reading in correlation potential')
                if(ios.lt.0) exit   !EOF
                if((abs(U_val-U).lt.1.0e-7_dp).and.(Occ_val.eq.nOcc)) then
                    tFoundCorrPot = .true.
                    exit
                endif
            enddo

            if(.not.tFoundCorrPot.and.(.not.tContinueConvergence)) then
                call stop_all(t_r,'Did not read in correlation potential corresponding to this run')
            else
                k=1
                do i=1,nImp
                    do j=1,nImp
                        v_loc(j,i) = CorrPot_tmp(k)
                        k = k+1
                    enddo
                enddo
            endif

            close(iunit)
        endif

        write(6,"(A)") "Read in correlation potential: "
        call writematrix(v_loc,"v_loc",.true.)
        do i=1,nImp
            do j=1,nImp
                if(abs(v_loc(i,j)-v_loc(j,i)).gt.1.0e-6_dp) then
                    call stop_all(t_r,'correlation potential not symmetric')
                endif
            enddo
        enddo
        if(tUHF) then
            !Test beta component too
            call writematrix(v_loc_b,"v_loc_b",.true.)
            do i=1,nImp
                do j=1,nImp
                    if(abs(v_loc_b(i,j)-v_loc_b(j,i)).gt.1.0e-6_dp) then
                        call stop_all(t_r,'beta correlation potential not symmetric')
                    endif
                enddo
            enddo
        endif

    end subroutine read_in_corrpot

    !Open an output file for the DMET convergence
    subroutine OpenDMETFile(iunit)
        use utils, only: get_free_unit,append_ext,append_ext_real
        implicit none
        integer, intent(out) :: iunit
        character(64) :: filename,filename2

        iunit = get_free_unit()
        call append_ext_real("DMET",U,filename)
        if(.not.tHalfFill) then
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')

    end subroutine OpenDMETFile
    
    subroutine DumpFCIDUMP()
        use utils, only: get_free_unit,append_ext_real
        use DetTools, only: GetHFInt_spinorb
        implicit none
        integer :: iunit,i,j,k,l,A,B,ex(2,2)
        real(dp) :: hel
        real(dp), allocatable :: temp(:,:),h0HF(:,:)
        character(len=64) :: filename
        
        iunit = get_free_unit()
        call append_ext_real('FCIDUMP',U,filename)
        open(unit=iunit,file=filename,status='unknown')
        write(iunit,'(2A6,I3,A9,I3,A6,I2,A)') '&FCI ','NORB=',nSites,', NELEC=',NEl,', MS2=',0,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,nSites
            write(iunit,'(I1,A1)',advance='no') 1,','
        enddo
        write(iunit,*) ""
        WRITE(iunit,'(A7,I1)') 'ISYM=',1
        WRITE(iunit,'(A5)') '&END'

        !Calculate hopping matrix in MO basis
        if(allocated(FullHFOrbs)) then
            write(6,'(A)') "Writing out FCIDUMP file corresponding to hartree--fock orbitals"

            allocate(temp(nSites,nSites))
            allocate(h0HF(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,h0HF,nSites)
            deallocate(temp)

            do i=1,nSites
                do j=1,nSites
                    A=(i*(i-1))/2+j
                    DO k=1,nSites
                        DO l=1,nSites
                            B=(k*(k-1))/2+l

                            !IF(B.lt.A) CYCLE
                            !IF((i.lt.j).and.(k.lt.l)) CYCLE
                            !IF((i.gt.j).and.(k.lt.l)) CYCLE

                            ex(1,1) = 2*i
                            ex(1,2) = 2*k
                            ex(2,1) = 2*j
                            ex(2,2) = 2*l
                            hel = GetHFInt_spinorb(ex,FullHFOrbs)
                            if(abs(hel).gt.1.0e-8_dp) then
                                WRITE(iunit,'(1X,G20.14,4I3)') hel,i,j,k,l
                            endif
                        enddo
                    enddo
                enddo
            enddo
            do i=1,nSites
                do j=1,i
                    if(abs(h0HF(i,j)).gt.1.0e-8_dp) then
                        WRITE(iunit,'(1X,G20.14,4I3)') h0HF(i,j),i,j,0,0
                    endif
                enddo
            enddo
            do i=1,nSites
                WRITE(iunit,'(1X,G20.14,4I3)') FullHFEnergies(i),i,0,0,0
            enddo
            WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
            close(iunit)
            deallocate(h0HF)
        else
            write(6,'(A)') "Writing out FCIDUMP file corresponding to non-interacting orbitals"
            allocate(temp(nSites,nSites))
            allocate(h0HF(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,HFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,HFOrbs,nSites,0.0_dp,h0HF,nSites)
            deallocate(temp)

            do i=1,nSites
                do j=1,nSites
                    A=(i*(i-1))/2+j
                    DO k=1,nSites
                        DO l=1,nSites
                            B=(k*(k-1))/2+l

                            !IF(B.lt.A) CYCLE
                            !IF((i.lt.j).and.(k.lt.l)) CYCLE
                            !IF((i.gt.j).and.(k.lt.l)) CYCLE

                            ex(1,1) = 2*i
                            ex(1,2) = 2*k
                            ex(2,1) = 2*j
                            ex(2,2) = 2*l
                            hel = GetHFInt_spinorb(ex,HFOrbs)
                            if(abs(hel).gt.1.0e-8_dp) then
                                WRITE(iunit,'(1X,G20.14,4I3)') hel,i,j,k,l
                            endif
                        enddo
                    enddo
                enddo
            enddo
            do i=1,nSites
                do j=1,i
                    if(abs(h0HF(i,j)).gt.1.0e-8_dp) then
                        WRITE(iunit,'(1X,G20.14,4I3)') h0HF(i,j),i,j,0,0
                    endif
                enddo
            enddo
            do i=1,nSites
                WRITE(iunit,'(1X,G20.14,4I3)') HFEnergies(i),i,0,0,0
            enddo
            WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
            close(iunit)
            deallocate(h0HF)

        endif

        !Now dump in the AO basis
        call append_ext_real('AO_FCIDUMP',U,filename)
        open(unit=iunit,file=filename,status='unknown')
        write(iunit,'(2A6,I3,A9,I3,A6,I2,A)') '&FCI ','NORB=',nSites,', NELEC=',NEl,', MS2=',0,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,nSites
            write(iunit,'(I1,A1)',advance='no') 1,','
        enddo
        write(iunit,*) ""
        WRITE(iunit,'(A7,I1)') 'ISYM=',1
        WRITE(iunit,'(A5)') '&END'
        if(tAnderson) then
            !Only first site correlated
            write(iunit,'(1X,G20.14,4I3)') U,1,1,1,1
        else
            do i=1,nSites
                write(iunit,'(1X,G20.14,4I3)') U,i,i,i,i
            enddo
        endif
        do i=1,nSites
            do j=1,i
!                if(tChemPot.and.(i.eq.j).and.(i.eq.1)) then
!                    write(iunit,'(1X,G20.14,4I3)') h0(i,j)-(U/2.0_dp),i,j,0,0
!                elseif(abs(h0(i,j)).gt.1.0e-8_dp) then
                if(abs(h0(i,j)).gt.1.0e-8_dp) then
                    WRITE(iunit,'(1X,G20.12,4I3)') h0(i,j),i,j,0,0
                endif
            enddo
        enddo
        WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine DumpFCIDUMP

    subroutine check_openmp
!$      use omp_lib
        implicit none
!$      integer(kind=OMP_integer_kind) :: maxthreads

        tOpenMP = .false.
!$      tOpenMP = .true.

        if(tOpenMP) then
            write(6,"(A)") " *** OpenMP compile detected *** "
!$          maxthreads = OMP_get_max_threads()
!$          max_omp_threads = maxthreads
!$          write(6,"(A,I7)") "Maxiumum number of threads to be used: ",max_omp_threads
!!$          OpenMP_wallstart = OMP_get_wtime()
!           Do not allow nested threads to start
!$          call OMP_set_nested(.false.)
        else
            write(6,"(A)") " No OpenMP optimizations"
            max_omp_threads = 1
        endif

        !Lets test it out
!$OMP PARALLEL
!$      write(6,*) "My thread number is: ",OMP_get_thread_num()
!$OMP END PARALLEL

        write(6,*) ""

    end subroutine check_openmp

    subroutine write_output_header()
        use utils, only: get_free_unit
        implicit none
        integer :: i,minres_unit
        real(dp) :: U_tmp
        logical :: exists
        character(len=*), parameter :: t_r='write_output_header'

        inquire(file='zMinResQLP.txt',exist=exists)
        if(exists) then
            minres_unit = get_free_unit()
            open(minres_unit,file='zMinResQLP.txt',status='old')
            close(minres_unit,status='delete')
        endif

        if(tAnderson) then
            write(6,"(A)") "Running:    o Anderson Model (single site)"
        elseif(tReadSystem) then
            write(6,"(A)") "Reading unknown system from files..."
        elseif(tSingFiss) then
            write(6,"(A)") "Reading Singlet Fission parameters..."
        else
            write(6,"(A)") "Running:    o Hubbard Model"
        endif
        if(tChemPot) then
            write(6,"(A)") "            o Chemical potential of -U/2 at impurity site for interacting system"
        else
            write(6,"(A)") "            o No chemical potential applied at impurity site"
        endif
        if(.not.tReadSystem) then
            if(tPeriodic) then
                write(6,"(A)") "            o PBCs employed"
            elseif(tAntiPeriodic) then
                write(6,"(A)") "            o APBCs employed"
            elseif(.not.(tPeriodic.or.tAntiPeriodic)) then
                write(6,"(A)") "            o Open boundary conditions employed"
            endif
            if(LatticeDim.eq.2) then
                write(6,"(A)") "            o 2-dimensional model" 
                if(CellShape.eq.1) then
                    write(6,"(A,I7,A,I7,A)") "            o Size of lattice: ",TDLat_Ni,' x ',TDLat_Nj,' at 45 degrees' 
                elseif(CellShape.eq.2) then
                    write(6,"(A,I7,A,I7,A)") "            o Size of lattice: ",nSites_x,' x ',nSites_y,' without tilt' 
                elseif(CellShape.eq.3) then
                    write(6,"(A,I7,A,I7,A)") "            o Size of lattice: ",nSites_x,' x ',nSites_y,' at 45 degrees tilt' 
                endif
                write(6,"(A,I7)") "            o Total lattice sites: ",nSites
            elseif(LatticeDim.eq.1) then
                write(6,"(A)") "            o 1-dimensional model" 
                write(6,"(A,I7)") "            o Size of lattice: ",nSites 
            else
                call stop_all(t_r,"Cannot determine dimensionality of system")
            endif
        else
            write(6,"(A,I9)") "            o Number of lattice sites: ",nSites
        endif
        if(tUHF) then
            write(6,"(A)") "            o *Unrestricted* bath construction: (Anti-)Ferromagnetic phase"
        else
            write(6,"(A)") "            o *Restricted* bath construction: Paramagnetic phase"
        endif
        write(6,"(A)") "            o Range of U values to consider: " 
        if(nU_Vals.eq.0) then
            !Sweeping through
            U_tmp=StartU
            do while((U_tmp.lt.max(StartU,EndU)+1.0e-5_dp).and.(U_tmp.gt.min(StartU,EndU)-1.0e-5_dp))
                write(6,"(A,F10.5)") "            o U = ",U_tmp 
                U_tmp=U_tmp+UStep
            enddo
        else
            do i = 1,nU_Vals
                write(6,"(A,F10.5)") "            o U = ",U_Vals(i) 
            enddo
        endif
        if(.not.tReadSystem) then
            if(tReadInCorrPot) then
                write(6,"(A,A)") "            o Correlation potentials for system will be read from file: ", trim(CorrPot_file)
                write(6,"(A)") "            o No DMET self-consistency of correlation potential"
            else
                if(.not.tAnderson) then
                    write(6,"(A,I8)") "            o Maximum iterations for DMET self-consistency: ", iMaxIterDMET
                    if(tSaveCorrPot) then
                        if(tHalfFill) then
                            write(6,"(A,I8)") "            o The correlation potential from the previous U value will " &
                                //"be used as a starting point for self-consistency" 
                        else
                            write(6,"(A,I8)") "            o The correlation potential from the previous electron number will " &
                                //"be used as a starting point for self-consistency"
                        endif
                    endif
                endif
            endif
        endif
        if(tHalfFill) then
            write(6,"(A)") "            o Only half filling to be considered"
        else
            if(tRampDownOcc) then
                write(6,"(A)") "            o Ramping down from half filling occupation of lattice"
            else
                write(6,"(A)") "            o Ramping up to half filling occupation of lattice"
            endif
        endif
        if(tSCFHF.and..not.tReadSystem) then
            write(6,"(A)") "            o Full self-consistent Hartree--Fock orbitals will be calculated"
        endif
        write(6,"(A,I7)") "            o Number of impurity sites: ",nImp 
        if(tDumpFCIDump) then
            write(6,"(A)") "            o Creating FCIDUMPs for system" 
        endif
        if(tFCIQMC) then
            if(nNECICores.eq.0) then
                write(6,"(A)") "            o Impurity solver: FCIQMC via call to serial NECI code"
            else
                write(6,"(A,I6)") "            o Impurity solver: FCIQMC via call to parallel NECI code. nCores: ",nNECICores
            endif
        elseif(tNonDirDavidson) then
            write(6,"(A)") "            o Impurity solver: Non-direct iterative Davidson diagonalizer" 
        elseif(tCompleteDiag) then
            write(6,"(A)") "            o Impurity solver: Complete diagonalization"
        else
            write(6,"(A)") "            o Impurity solver: Direct Davidson diagonalizer via call to FCI code"
        endif
        if(tNIResponse) then
            write(6,"(A)") "            o Calculating non-interacting linear response function" 
        endif
        if(tTDAResponse) then
            write(6,"(A)") "            o Calculating TDA linear response function" 
        endif
        if(tRPAResponse) then
            write(6,"(A)") "            o Calculating RPA linear response function" 
        endif
        if(tEC_TDA_Response) then
            write(6,"(A)") "            o Calculating externally-contracted MC-TDA DMET linear response function" 
            if(tDDResponse) then
                write(6,"(A)") "                o Local density response calculated" 
            endif
            if(tChargedResponse) then
                write(6,"(A)") "                o Local Greens function calculated" 
            endif
            if(tCharged_MomResponse) then
                write(6,"(A)") "                o Momentum-resolved Greens functions calculated" 
            endif

        endif
        if(tIC_TDA_Response) then
            write(6,"(A)") "            o Calculating internally-contracted MC-TDA DMET linear response function" 
        endif
        if(tMFResponse.or.tLR_DMET) then
            write(6,"(A,F13.8)") "            o Spectral broadening for linear response functions: ",dDelta
        endif
        write(6,"(A)") ""
        write(6,"(A)") ""

    end subroutine write_output_header
    
    subroutine end_calc()
!!$      use omp_lib
        implicit none
        real(dp) :: OpenMP_wallend 
            
        call deallocate_mem()
        call end_timing()
        call print_timing_report()
!        if(tOpenMP) then
!!$          OpenMP_wallend = OMP_get_wtime()
!            write(6,"(A)") ""
!            write(6,"(A,F14.2)") "Total calculational walltime per OpenMP thread: ",  &
!                OpenMP_wallend - OpenMP_wallstart
!        endif

    end subroutine end_calc

    !Deallocate memory which is still allocated
    subroutine deallocate_mem()
        implicit none

        if(allocated(U_Vals)) deallocate(U_Vals)
        if(allocated(TD_Imp_Lat)) deallocate(TD_Imp_Lat,TD_Imp_Phase)

    end subroutine deallocate_mem


End Program RealHub
