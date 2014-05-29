module SelfConsistentUtils
    use const
    use errors
    use Globals
    use utils, only: get_free_unit,append_ext
    use SC_Data
    use mat_tools, only: MakeBlockHermitian,writematrixcomp,writevector
    use mat_tools, only: Add_Nonlocal_comp_inplace,var_to_couplingind
    use matrixops, only: mat_inv
    use LRSolvers, only: GetNextkval
    implicit none

    contains

    !Calculate and write out the bandstructure
    !If SelfEnergy specified, include the self-energy contribution (therefore a correlated bandstructure)
    !Else, it is just a single-particle.
    !If h_lat is not specified, then it will just use 
    !TODO: Include spline interpolation for a finer resolution of the bandstructure?
    !Note for gnuplot: for [IDX=0:623:50] 'Bandstructure_SE' i IDX u 1:($3 + $2 +0.5*IDX) w lp
    subroutine CalcBandstructure(n,mu,fileroot,tag,h_lat,SelfEnergy,FreqPoints)
        integer, intent(in) :: n    !Number of frequency points
        real(dp), intent(in) :: mu  !Chemical potential
        character(len=*), intent(in) :: fileroot
        integer, intent(in), optional :: tag
        complex(dp), intent(in), optional :: h_lat(nSites,nSites)
        complex(dp), intent(in), optional :: SelfEnergy(nImp,nImp,n)
        real(dp), intent(in), optional :: FreqPoints(n)
        
        integer :: iunit,i,j,k,lWork,info
        real(dp) :: Omega,EVals(nImp)
        real(dp), allocatable :: rWork(:)
        complex(dp) :: Mat(nImp,nImp),MatInv(nImp,nImp)
        complex(dp), allocatable :: KBlocks(:,:,:),cWork(:)
        complex(dp), allocatable :: h0c(:,:)
        character(64) :: filename
        logical :: tFinishedk
        character(len=*), parameter :: t_r='CalcBandstructure'

        write(6,*) "Calculating bandstructure..."

        if(present(SelfEnergy).and.present(h_lat)) then
            call stop_all(t_r,'Both a lattice hamiltonian and a self-energy specified. Doesnt seem right')
        endif

        !Open file
        if(present(tag)) then
            call append_ext(fileroot,tag,filename)
        else
            filename = fileroot
        endif
        iunit = get_free_unit()
        open(unit=iunit,file=filename,status='unknown')

        allocate(KBlocks(nImp,nImp,nKPnts))
        if(.not.present(h_lat)) then
            allocate(h0c(nSites,nSites))
            do i = 1,nSites
                do j = 1,nSites
                    h0c(j,i) = cmplx(h0(j,i),zero,dp)
                enddo
                h0c(i,i) = h0c(i,i) + cmplx(mu,zero,dp)
            enddo
            call ham_to_KBlocks(h0c,KBlocks)
            deallocate(h0c)
        else
            call ham_to_KBlocks(h_lat,KBlocks)
        endif

        if(present(SelfEnergy)) then
            !This can actually provide a correlated bandstructure with non-unit spectral weight.
            !We have to write out the whole range of frequencies for each kpoint
            !By default, we only write out 40 kpoints, unless specified otherwise
            !This can be controlled with the KPNT_CALCS option in the LINEAR_RESPONSE block
            k = 0
            do while(.true.)
                call GetNextkVal(k,tFinishedk)
                if(tFinishedk) exit 
                
                !Write out the header
                write(iunit,"(A,F8.4)",advance='no') '"k = ',KPnts(1,k)
                do i = 2,LatticeDim
                    write(iunit,"(A,F8.4)",advance='no') ', ',KPnts(i,k)
                enddo
                write(iunit,"(A)") ' "'

                !Now loop over frequencies.
                do i = 1,n
                    if(present(FreqPoints)) then
                        call GetOmega(Omega,i,.false.,FreqPoints=FreqPoints)
                    else
                        call GetOmega(Omega,i,.false.)
                    endif

                    !Construct the matrix
                    Mat(:,:) = - KBlocks(:,:,k) - SelfEnergy(:,:,i)
                    do j = 1,nImp
                        Mat(j,j) = Mat(j,j) + cmplx(Omega + mu,dDelta)
                    enddo
                    call mat_inv(Mat,MatInv,nImp)

                    write(iunit,"(G22.10)",advance='no') Omega
                    do j = 1,nImp-1 
                        write(iunit,"(F20.12)",advance='no') -aimag(MatInv(j,j))
                    enddo
                    write(iunit,"(F20.12)") -aimag(MatInv(nImp,nImp))
                enddo
                write(iunit,"(A)") ""
                write(iunit,"(A)") ""
            enddo
        else
            !Easy peasy. Just write out the lattice eigenvalues with their klabel
            !TODO: Splines
            allocate(rWork(max(1,3*nImp-2)))
            do k = 1,nKPnts
                !Diagonalize
                allocate(cWork(1))
                lWork = -1
                info = 0
                call zheev('V','U',nImp,KBlocks(:,:,k),nImp,EVals,cWork,lWork,rWork,info)
                if(info.ne.0) call stop_all(t_r,'Error in diag 1')
                lWork = int(real(cWork(1))) + 1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zheev('V','U',nImp,KBlocks(:,:,k),nImp,EVals,cWork,lWork,rWork,info)
                if(info.ne.0) call stop_all(t_r,'Error in diag 2')
                deallocate(cWork)

                !Write out the eigenvalues.
                do i = 1,LatticeDim
                    write(iunit,"(F20.12)",advance='no') KPnts(i,k)
                enddo
                do i = 1,nImp-1
                    write(iunit,"(F20.12)",advance='no') EVals(i)
                enddo
                write(iunit,"(F20.12)") EVals(nImp)
            enddo
            deallocate(rWork)
        endif
        deallocate(KBlocks)
        close(iunit)

    end subroutine CalcBandstructure
    
    subroutine SetChemPot(GFChemPot)
        implicit none
        real(dp), intent(out) :: GFChemPot
        character(len=*), parameter :: t_r='SetChemPot'

        if(.not.tHalfFill) call stop_all(t_r,'Chemical potential currently only set for half-filling')
        !Set chemical potential (This will only be right for half-filling)
        if(.not.tAnderson) then
            !Hubbard chemical potential
            GFChemPot = U/2.0_dp
        else
            !For anderson model, we want h0v to be h0 with -U/2 on impurity
            !Then, we want no chemical potential for the linear response.
            GFChemPot = zero
        endif

    end subroutine SetChemPot

    !Check options, and write out what it going on
    subroutine CheckSCOptions(iCorrFnTag)
        implicit none
        integer, intent(in) :: iCorrFnTag   !In the end, this will be a global input variable
        character(len=*), parameter :: t_r='CheckSCOptions'
    
        write(6,"(A)") ""
        write(6,"(A)") " *** Entering code of self-consistent optimizating of spectral functions ***"
        write(6,"(A)") ""

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
            
        if((.not.tSC_StartwGSCorrPot).and.tRealSpaceSC) then
            write(6,"(A)") "With real-space couplings being optimized and starting from the " &
                //"uncorrelated lattice h, no local couplings will be updated."
            write(6,"(A)") "Are you sure this is what you want to do? No local interactions in lattice at all."
            call warning(t_r,'No local lattice couplings over impurity')
        endif

    end subroutine CheckSCOptions
        
    subroutine SetLatticeParams(iLatParams)
        implicit none
        integer, intent(out) :: iLatParams
        integer :: ks,params_per_k,i,MaxBlocks
        character(len=*), parameter :: t_r='SetLatticeParams'

        if(tRealSpaceSC) then

            !Now, add non-local terms, which retain periodicity
            if(mod(nKPnts,2).eq.0) then
                !Even number of kpoints overall (Odd number of nonlocal coupling blocks)
                MaxBlocks = (nSites/nImp)/2
            else
                !Odd number of kpoints overall (Even number of nonlocal coupling blocks)
                MaxBlocks = (nSites/nImp - 1)/2
            endif
            write(6,"(A,I8)") "Maximum possible number of non-local real-space coupling unit cells: ",MaxBlocks
            if((iNonLocBlocks.eq.0).or.(iNonLocBlocks.ge.MaxBlocks)) then
                !We want all the possible non-redundant couplings blocks
                iNonLocBlocks = 0
                if(mod((nSites-nImp)/nImp,2).ne.0) then
                    !Even number of kpoints overall (Odd number of nonlocal coupling blocks)
                    tOddFullNonlocCoups = .true.    !We have one coupling block which has only one copy. This has special properties (only triangular number of independent variables)
                    iLatParams = nImp*((nSites-(2*nImp))/2) + (nImp*(nImp-1)/2)
                else
                    !Odd number of kpoints overall (Even number of nonlocal coupling blocks)
                    tOddFullNonlocCoups = .false.
                    iLatParams = nImp*((nSites-nImp)/2)
                endif
                write(6,"(A,I8)") "All non-local couplings will be optimized."
            else
                write(6,"(A,I8)") "Truncated real-space coupling. Number of unit cells: ",iNonLocBlocks
                tOddFullNonlocCoups = .false.
                iLatParams = iNonLocBlocks*nImp*nImp
            endif

            write(6,"(A,I8)") "Total number of (real) adjustable parameters in non-local couplings: ",iLatParams

        else
            
            if(tConstrainKSym) then
                !e(k) = e(-k), and we are always using a uniform mesh
                if(mod(nKPnts,2).eq.0) then
                    !Even number of kpoints 
                    !If gamma-centered mesh, then have nSites/2 + 1 independent parameters 
                    !(we are sampling k=0 and BZ boundary which dont pair)
                    !If Shifted mesh, then we have nSites/2 independent parameters
                    if(tShift_Mesh) then
                        nIndKPnts = nKPnts/2
                    else
                        nIndKPnts = (nKPnts/2) + 1
                    endif
                else
                    !This is independent of whether we have a shifted mesh or not
                    nIndKPnts = (nKPnts+1)/2    
                endif

            else
                nIndKPnts = nKPnts
            endif

            if(tConstrainphsym) then
                !First nImp is the first column
                !Second nImp-2 is the second column (taking into account herm)
                !etc...
                params_per_k = 0
                do i = nImp,0,-2
                    params_per_k = params_per_k + i
                enddo
            else
                !Still constrain hermiticity. Allow for complex off-diagonal matrix elements
                params_per_k = (nImp*(nImp+1)/2) + (nImp*(nImp-1)/2)
            endif

            iLatParams = params_per_k * nIndKPnts

            write(6,"(A,I8)") "Total number of independent k-points to optimize: ",nIndKPnts
            write(6,"(A,I8)") "Total number of independent (real) parameters per kpoint: ",params_per_k
            write(6,"(A,I8)") "Total number of (real) adjustable parameters in non-local couplings: ",iLatParams

        endif

    end subroutine SetLatticeParams

    subroutine WriteLatticeParams(iLatParams,LatParams)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: LatParams(iLatParams)
        integer :: iunit,i,j,k,ks,kspace_params
        complex(dp), allocatable :: KBlocks(:,:,:)
            
        write(6,"(A)") "Writing lattice k-space hamiltonian..."
        iunit = get_free_unit()

        open(unit=iunit,file='LatkBlocks',status='unknown')
        allocate(KBlocks(nImp,nImp,nKPnts))
        call LatParams_to_KBlocks(iLatParams,LatParams,U/2.0_dp,KBlocks)
!        write(iunit,*) iLatParams
        write(iunit,"(3I9)") nImp,nKPnts,LatticeDim
        do k = 1,nKPnts
            !First, write out k-vector
            do i = 1,LatticeDim-1
                write(iunit,"(F20.12)",advance='no') KPnts(i,k)
            enddo
            write(iunit,"(I9,F20.12)") k, KPnts(LatticeDim,k)

            !Then write out the block
            do i = 1,nImp
                do j = 1,nImp-1
                    write(iunit,"(2G20.12)",advance='no') real(KBlocks(i,j,k),dp),aimag(KBlocks(i,j,k))
                enddo
                write(iunit,"(2G20.12)") real(KBlocks(i,nImp,k),dp),aimag(KBlocks(i,nImp,k))
            enddo
        enddo
        deallocate(KBlocks)
        close(iunit)

!        if(tRealSpaceSC) then
!            !Convert to k-space... HACK
!            !TODO: In the future, just write out real space coupling
!            if(tShift_Mesh) then
!                ks = nKPnts/2
!            else
!                ks = (nKPnts/2) + 1
!            endif
!            !because they're complex
!            kspace_params = ((nImp*(nImp+1)/2) + (nImp*(nImp-1)/2)) * ks
!            write(iunit,*) kspace_params
!            do k = 1,ks
!                do i = 1,nImp
!                    !Run through columns
!                    do j = i,nImp
!                        !Run through rows
!                        if(i.eq.j) then
!                            write(iunit,*) real(KBlocks(j,i,k),dp)
!                        else
!                            write(iunit,*) real(KBlocks(j,i,k),dp)
!                            write(iunit,*) aimag(KBlocks(j,i,k))
!                        endif
!                    enddo
!                enddo
!            enddo
!            deallocate(KBlocks)
!        else
!            !Write out k-space hamiltonian
!            allocate(KBlocks(nImp,nImp,nKPnts))
!            call LatParams_to_KBlocks(iLatParams,LatParams,U/2.0_dp,KBlocks)
!        endif

    end subroutine WriteLatticeParams

    subroutine InitLatticeParams(iLatParams,LatParams,mu,ham)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(out) :: ham(nSites,nSites)  !Return the lattice hamiltonian too

        complex(dp), allocatable :: KBlocks(:,:,:)
        real(dp) :: KPnt_read(LatticeDim),kBlk_read(2*nImp)
        integer :: i,j,iloclatcoups,iunit,k_read,nImp_read,nKPnts_read,LatticeDim_read,k
        logical :: exists
        character(len=*), parameter :: t_r='InitLatticeParams'

        if(tReadCouplings) then

            if(tRealSpaceSC) then
                call stop_all(t_r,'Reading parameters not available for real-space opt yet')
            endif

            write(6,"(A)") "Reading lattice k-space hamiltonian..."
            LatParams(:) = zero

            inquire(file='LatkBlocks_Read',exist=exists)
            if(.not.exists) call stop_all(t_r,'LatkBlocks_Read file does not exist')
            iunit = get_free_unit()
            open(unit=iunit,file='LatkBlocks_Read',status='old')

            read(iunit,*) nImp_read,nKPnts_read,LatticeDim_read
            if(nImp_read.ne.nImp) then
                call stop_all(t_r,'Number of impurities has changed from LatkBlock_read. Cannot read in lattice h')
            elseif(LatticeDim_read.ne.LatticeDim) then
                call stop_all(t_r,'Dimensionality of lattice has changed from LatkBlock_read. Can not read in lattice h')
            endif
            if(nKPnts_read.ne.nKPnts) then
                call stop_all(t_r,'Reading in from smaller lattice? This should be possible. Bug ghb24 to code up splines')
            else
                write(6,"(A,I9)") "Reading in a commensurate lattice with number of k-points: ",nKPnts
                allocate(KBlocks(nImp,nImp,nKPnts))
                KBlocks(:,:,:) = zzero

                do k = 1,nKPnts
                    read(iunit,*) k_read,KPnt_read(1:LatticeDim)
                    if(k_read.ne.k) call stop_all(t_r,'Not reading kpoints sequentially')
                    do i = 1,LatticeDim
                        if(abs(KPnt_read(i)-KPnts(i,k)).gt.1.0e-7_dp) then
                            call stop_all(t_r,'I thought the kpoint meshes were commensurate, but they are not')
                        endif
                    enddo
                    !Now read this kBlock
                    do i = 1,nImp
                        read(iunit,*) kBlk_read(1:2*nImp)   !real and imaginary
                        do j = 1,nImp
                            KBlocks(i,j,k) = cmplx(kBlk_read(j*2-1),kBlk_read(j*2),dp)
                        enddo
                    enddo
                enddo

                call KBlocks_to_LatParams(iLatParams,LatParams,KBlocks)
                deallocate(KBlocks)
                if(tImposephsym.or.tImposeKSym) then
                    write(6,"(A)") "Imposing symmetry constrains on read in lattice hamiltonian"
                    call ImposeSym_2(iLatParams,LatParams,mu)
                endif
                call LatParams_to_ham(iLatParams,LatParams,mu,ham)
                !call WriteLatticeParams(iLatParams,LatParams)
            endif

!            read(iunit,*) iloclatcoups
!            if(iloclatcoups.gt.iLatParams) then
!                call stop_all(t_r,'Cannot read in from larger lattice / check inputs compatible')
!            elseif(iloclatcoups.eq.iLatParams) then
!                write(6,"(A)") "In lattice read, number of free parameters the same. Assuming that the "    &
!                    //"k-point mesh is commensurate, and just using these values."
!                do i = 1,iloclatcoups
!                    read(iunit,*) LatParams(i)
!                enddo
!                if(abs(dShiftLatticeEvals).gt.1.0e-8_dp) then
!                    call ShiftLatParams(iLatParams,LatParams,dShiftLatticeEvals)
!                endif
!            else
!            endif
!
!            if(tImposephsym.or.tImposeKSym) then
!                call ImposeSym_2(iLatParams,LatParams,mu)
!            endif
!            call LatParams_to_ham(iLatParams,LatParams,mu,ham)
        else
            !Not reading in
            if(tSC_StartwGSCorrPot) then
                !Just set to h0v
                ham(:,:) = zzero
                do i = 1,nSites
                    do j = 1,nSites
                        ham(j,i) = cmplx(h0v(j,i),zero,dp)
                    enddo
                enddo
            else
                !Just set to h0
                !Note that if we are doing real-space couplings, this will *never* allow for updates of the local hamiltonian
                ham(:,:) = zzero
                do i = 1,nSites
                    do j = 1,nSites
                        ham(j,i) = cmplx(h0(j,i),zero,dp)
                    enddo
                    ham(i,i) = ham(i,i) + mu    !Add the correlation potential at least
                enddo
            endif
            call ham_to_LatParams(iLatParams,LatParams,ham)
        endif

    end subroutine InitLatticeParams

    subroutine ShiftLatParams(iLatParams,LatParams,dShift)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(inout) :: LatParams(iLatParams)
        real(dp), intent(in) :: dShift

    end subroutine ShiftLatParams

    subroutine ham_to_KBlocks(ham,KBlocks)
        complex(dp), intent(in) :: ham(nSites,nSites) 
        complex(dp), intent(out) :: KBlocks(nImp,nImp,nKPnts)
        complex(dp), allocatable :: ctemp(:,:)
        integer :: k,ind_1,ind_2
        
        KBlocks(:,:,:) = zzero
!$OMP PARALLEL DO PRIVATE(ind_1,ind_2,ctemp)
        do k = 1,nKPnts
            allocate(ctemp(nSites,nImp))
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            call ZGEMM('N','N',nSites,nImp,nSites,zone,ham,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ctemp,nSites)
            call ZGEMM('C','N',nImp,nImp,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ctemp,nSites,zzero,KBlocks(:,:,k),nImp)
            deallocate(ctemp)
        enddo
!$OMP END PARALLEL DO

    end subroutine ham_to_KBlocks

    !Real space hamiltonian to non-local lattice parameters
    subroutine ham_to_LatParams(iLatParams,LatParams,ham,tUpdateLocalHam)
        integer, intent(in) :: iLatParams
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(in) :: ham(nSites,nSites)
        logical, intent(in), optional :: tUpdateLocalHam
        complex(dp), allocatable :: KBlocks(:,:,:)
        integer :: i,j,ind,ind_r,ind_c,ind_block
        logical :: tUpdateLocalHam_
        character(len=*), parameter :: t_r='ham_to_LatParams'

        if(tRealSpaceSC) then

            if(present(tUpdateLocalHam)) then
                tUpdateLocalHam_ = tUpdateLocalHam
            else
                tUpdateLocalHam_ = .true.
            endif

            do i = 1,iLatParams
                call var_to_couplingind(i,nImp,ind_r,ind_c,ind_block)
                LatParams(i) = real(ham(ind_r,ind_c+(ind_block*nImp)),dp)
            enddo

            if(tUpdateLocalHam_) then
                !Also, update the global array which gives the block diagonals of the lattice hamiltonian in real space
                if(.not.allocated(LatDiagBlocks)) then
                    allocate(LatDiagBlocks(nImp,nImp))
                endif
                LatDiagBlocks(:,:) = ham(1:nImp,1:nImp)
            endif

        else

            !First, ham to k_blocks
            allocate(KBlocks(nImp,nImp,nKPnts))
            call ham_to_KBlocks(ham,KBlocks)

            call KBlocks_to_LatParams(iLatParams,LatParams,KBlocks)
            deallocate(KBlocks)
        endif

    end subroutine ham_to_LatParams

    !Convert a set of variational parameters back into a real space one-electron hamiltonian
    subroutine LatParams_to_ham(iLatParams,LatParams,mu,ham)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: LatParams(iLatParams)
        complex(dp), intent(out) :: ham(nSites,nSites)
        complex(dp), allocatable :: KBlocks(:,:,:)
        complex(dp), allocatable :: ctemp(:,:),vars(:)
        integer :: k,ind_1,ind_2,kPnt,ind,i,j
        character(len=*), parameter :: t_r='LatParams_to_ham'
        
        ham(:,:) = zzero

        if(tRealSpaceSC) then
            !From real space lattice couplings to a hamiltonian

            !Put in diagonals
            do kPnt = 1,nKPnts
                ind_1 = ((kPnt-1)*nImp) + 1
                ind_2 = nImp*kPnt

                ham(ind_1:ind_2,ind_1:ind_2) = LatDiagBlocks(:,:)
            enddo
            !ham is complex, so unfortunately, we have to transfer variables to complex array
            allocate(vars(iLatParams))
            vars(:) = zzero
            do i = 1,iLatParams
                vars(i) = cmplx(LatParams(i),zero,dp)
            enddo
            call Add_Nonlocal_comp_inplace(ham,vars,iLatParams,tAdd=.true.)
            deallocate(vars)

        else

            allocate(KBlocks(nImp,nImp,nKPnts))
            call LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)

            !Now, rotate back into real space from kspace
            allocate(ctemp(nSites,nImp))
            do k = 1,nKPnts
                ind_1 = ((k-1)*nImp) + 1
                ind_2 = nImp*k
                call ZGEMM('N','N',nSites,nImp,nImp,zone,RtoK_Rot(:,ind_1:ind_2),nSites,KBlocks(:,:,k),nImp,zzero,ctemp,nSites)
                call ZGEMM('N','C',nSites,nSites,nImp,zone,ctemp,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zone,ham,nSites)
            enddo
            deallocate(ctemp)
            deallocate(KBlocks)

    !        do i = 1,nKPnts
    !            ind_1 = (i-1)*nImp + 1
    !            ind_2 = i*nImp
    !            ham(ind_1:ind_2,ind_1:ind_2) = KBlocks(:,:,i)
    !        enddo
    !        allocate(ctemp(nSites,nSites))
    !        call ZGEMM('N','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham,nSites,zzero,ctemp,nSites)
    !        call ZGEMM('N','C',nSites,nSites,nSites,zone,ctemp,nSites,RtoK_Rot,nSites,zzero,ham,nSites)
    !        deallocate(ctemp)

        endif

    end subroutine LatParams_to_ham

    subroutine KBlocks_to_LatParams(iLatParams,LatParams,KBlocks)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(in) :: KBlocks(nImp,nImp,nKPnts)

        integer :: i,j,k,ind,ks,ind_1,ind_2,ind_r,ind_c,ind_block
        complex(dp), allocatable :: LatticeCoups_(:,:),LattCoupTemp(:,:)
        complex(dp), allocatable :: BlockTemp(:,:),BlockTemp_2(:,:),BlockTemp_3(:,:)
        character(len=*), parameter :: t_r='KBlocks_to_LatParams'

        LatParams(:) = zero

        if(tRealSpaceSC) then
            
            allocate(LatticeCoups_(nImp,nSites))
            LatticeCoups_(:,:) = zzero

            !For testing
!            allocate(BlockTemp(nSites,nSites))
!            allocate(BlockTemp_2(nSites,nSites))
!            BlockTemp(:,:) = zzero
!            do k = 1,nKPnts
!                ind_1 = ((k-1)*nImp) + 1
!                ind_2 = nImp*k
!
!                BlockTemp(ind_1:ind_2,ind_1:ind_2) = KBlocks(:,:,k)
!            enddo
!            call ZGEMM('N','C',nSites,nSites,nSites,zone,BlockTemp,nSites,RtoK_Rot,nSites,zzero,BlockTemp_2,nSites)
!            call ZGEMM('N','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,BlockTemp_2,nSites,zzero,BlockTemp,nSites)
!
!            LatticeCoups_(:,:) = BlockTemp(1:nImp,nImp+1:nImp+iLatticeCoups)
!            deallocate(BlockTemp,BlockTemp_2)

            allocate(LattCoupTemp(nImp,nSites))
            LattCoupTemp(:,:) = zzero
            allocate(BlockTemp(nImp,nImp))
            allocate(BlockTemp_2(nImp,nImp))
            allocate(BlockTemp_3(nImp,nImp))

            do k = 1,nKPnts
                ind_1 = ((k-1)*nImp) + 1
                ind_2 = nImp*k

                BlockTemp_2(:,:) = RtoK_Rot(1:nImp,ind_1:ind_2)
                BlockTemp_3(:,:) = KBlocks(:,:,k)

                call ZGEMM('N','N',nImp,nImp,nImp,zone,BlockTemp_2,nImp,BlockTemp_3,nImp,zzero,BlockTemp,nImp)

                LattCoupTemp(:,ind_1:ind_2) = BlockTemp(:,:)
            enddo
            call ZGEMM('N','C',nImp,nSites,nSites,zone,LattCoupTemp,nImp,RtoK_Rot,nSites,zzero,LatticeCoups_,nImp)

            deallocate(BlockTemp,BlockTemp_2,BlockTemp_3,LattCoupTemp)
            
            do i = 1,iLatParams
                call var_to_couplingind(i,nImp,ind_r,ind_c,ind_block)
                if(abs(aimag(LatticeCoups_(ind_r,ind_c+(ind_block*nImp)))).gt.1.0e-7_dp) then
                    write(6,"(A,I7,G17.10)") "Removing complex component of coupling: ",    &
                        i,abs(aimag(LatticeCoups_(ind_r,ind_c+(ind_block*nImp))))
                endif
                LatParams(i) = real(LatticeCoups_(ind_r,ind_c+(ind_block*nImp)),dp)
            enddo
            deallocate(LatticeCoups_)
        else
            
            ind = 1
            do k = 1,nIndKPnts
                if(tConstrainphsym) then
                    do i = 1,nImp
                        do j = i,nImp-i+1
                            if(i.eq.j) then
                                if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                                LatParams(ind) = real(KBlocks(j,i,k),dp)
                                ind = ind + 1
                            else
                                if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                                LatParams(ind) = real(KBlocks(j,i,k),dp)
                                LatParams(ind+1) = aimag(KBlocks(j,i,k))
                                ind = ind + 2
                            endif
                        enddo
                    enddo
                else
                    do i = 1,nImp
                        !Run through columns
                        do j = i,nImp
                            !Run through rows
                            if(i.eq.j) then
                                if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                                LatParams(ind) = real(KBlocks(j,i,k),dp)
                                ind = ind + 1
                            else
                                if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                                LatParams(ind) = real(KBlocks(j,i,k),dp)
                                LatParams(ind+1) = aimag(KBlocks(j,i,k))
                                ind = ind + 2
                            endif
                        enddo
                    enddo
                endif
            enddo
        endif

    end subroutine KBlocks_to_LatParams

!    subroutine RealSpace_EValsVecs(EVecs_k,EVals_k,EVecs_r,EVals_r)
!        complex(dp), intent(in) :: EVecs_k(nImp,nImp,nKPnts)
!        complex(dp), intent(in) :: EVals_k(nSites)
!        complex(dp), intent(out) :: EVecs_r(nSites,nSites)
!        complex(dp), intent(out) :: EVals_r(nSites)
!
!        do k = 1,nKPnts
!
!            call ZGEMM('N','N',
!
!    end subroutine RealSpace_EValsVecs

!Function to return the index in the variables array from a k-space optimization from the kblocks
    elemental subroutine kspace_var_to_coupind(k_ind,ind_r,ind_c,blocksize,ind)
        implicit none
        integer, intent(in) :: k_ind,ind_r,ind_c,blocksize
        integer, intent(out) :: ind

        ind = k_ind*blocksize**2 - (blocksize-ind_c+1)**2 + ind_r - ind_c + 1 + max(ind_r-ind_c-1,0)
!        ind = (k_ind-1)*blocksize**2 + blocksize**2 - (blocksize-ind_c+1)**2 + ind_r - ind_c + 1 + max(ind_r-ind_c-1,0)

    end subroutine kspace_var_to_coupind


    subroutine KBlocks_to_diag(KBlocks,EVecs,EVals)
        complex(dp), intent(in) :: KBlocks(nImp,nImp,nKPnts)
        complex(dp), intent(out) :: EVecs(nImp,nImp,nKPnts)
        complex(dp), intent(out) :: EVals(nSites)
        complex(dp) :: KBlock(nImp,nImp)
        complex(dp), allocatable :: cWork(:)
        real(dp) :: Bands_k(nImp)
        real(dp), allocatable :: rWork(:)
        integer :: k,lWork,ierr,i
        character(len=*), parameter :: t_r='KBlocks_to_diag'

        EVals = zzero
        EVecs = zzero

        allocate(rWork(max(1,3*nImp-2)))
        do k = 1,nKPnts

            EVecs(:,:,k) = KBlocks(:,:,k)
            !Diagonalize each kBlock
            allocate(cWork(1))
            lWork = -1
            ierr = 0
            call zheev('V','U',nImp,EVecs(:,:,k),nImp,Bands_k,cWork,lWork,rWork,ierr)
            if(ierr.ne.0) call stop_all(t_r,'Error in diag')
            lWork = int(real(cWork(1))) + 1
            deallocate(cWork)
            allocate(cWork(lWork))
            call zheev('V','U',nImp,EVecs(:,:,k),nImp,Bands_k,cWork,lWork,rWork,ierr)
            if(ierr.ne.0) call stop_all(t_r,'Error in diag')
            deallocate(cWork)

            do i = 1,nImp
                EVals(((k-1)*nImp)+i) = cmplx(Bands_k(i),zero,dp)
            enddo

        enddo
        deallocate(rWork)

    end subroutine KBlocks_to_diag

    subroutine LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: LatParams(iLatParams)
        complex(dp), intent(out) :: KBlocks(nImp,nImp,nKPnts)
        complex(dp), allocatable :: ham(:,:),ztemp(:,:),ztemp_2(:,:),temp_block(:,:)
        integer :: ind,i,j,ks,k,ind_1,ind_2
        character(len=*), parameter :: t_r='LatParams_to_KBlocks'

        if((nImp.ne.1).and.(mod(nImp,2).ne.0)) then
            call stop_all(t_r,'Cannot deal with non-1 odd numbers of impurity sites')
        endif

        KBlocks(:,:,:) = zzero

        if(tRealSpaceSC) then
            !Real-space optimization
            allocate(ham(nSites,nSites))    !Temp storage of complex hamiltonian
            call LatParams_to_ham(iLatParams,LatParams,mu,ham)
            allocate(ztemp(nSites,nImp))
            allocate(ztemp_2(nSites,nImp))
            allocate(temp_block(nImp,nImp))
            do k = 1,nKPnts
                ind_1 = ((k-1)*nImp) + 1
                ind_2 = nImp*k
                ztemp_2(:,:) = RtoK_Rot(:,ind_1:ind_2)
                call ZGEMM('N','N',nSites,nImp,nSites,zone,ham,nSites,ztemp_2,nSites,zzero,ztemp,nSites)
                call ZGEMM('C','N',nImp,nImp,nSites,zone,ztemp_2,nSites,ztemp,nSites,zzero,temp_block,nImp)
                KBlocks(:,:,k) = temp_block(:,:)
            enddo
            deallocate(ham,ztemp,ztemp_2,temp_block)
            !TODO: ph symmetry imposed?
        else
            !k-space hamiltonian optimization
            
            ind = 1
            do k = 1,nIndKPnts
                if(tConstrainphsym) then
                    do i = 1,nImp
                        do j = i,nImp-i+1
                            if(i.eq.j) then
                                if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                                KBlocks(j,i,k) = cmplx(LatParams(ind),zero,dp)
                                ind = ind + 1
                                !For the diagonals, we want to flip the energy about the chemical potential
                                KBlocks((nImp-i)+1,(nImp-j)+1,k) = cmplx(2.0_dp*mu,zero,dp) - KBlocks(j,i,k)
                            else
                                if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                                KBlocks(j,i,k) = cmplx(LatParams(ind),LatParams(ind+1),dp)
                                ind = ind + 2
                                !Put in the same coupling
                                !Swap indices & reverse order
                                KBlocks((nImp-i)+1,(nImp-j)+1,k) = KBlocks(j,i,k)
                            endif
                        enddo
                    enddo
                else
                    !The variables are the lower triangle of the matrix for this kpoint
                    do i = 1,nImp
                        !Run through columns
                        do j = i,nImp
                            !Run through rows
                            if(i.eq.j) then
                                if(ind.gt.iLatParams) then
                                    write(6,*) "i,j,k: ",i,j,k
                                    write(6,*) "ind: ",ind
                                    write(6,*) "iLatParams: ",iLatParams
                                    call stop_all(t_r,'Incorrect indexing 1')
                                endif
                                KBlocks(j,i,k) = cmplx(LatParams(ind),zero,dp)
                                ind = ind + 1
                            else
                                if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing 2')
                                KBlocks(j,i,k) = cmplx(LatParams(ind),LatParams(ind+1),dp)
                                ind = ind + 2
                            endif
                        enddo
                    enddo
                endif
                call MakeBlockHermitian(KBlocks(:,:,k),nImp)
            enddo

            if(tConstrainKSym) then
                !We have only filled up half the states. Fill the others.
                !Update this so that it uses a function to map to equivalent kpoints
                !Currently this is only going to work for 1D

                if(mod(nKPnts,2).eq.0) then
                    !Even number of kpoints
                    if(tShift_Mesh) then
                        !No gamma point sampled. All k-points symmetric.
                        !Mirror the k-space hamiltonian
                        do i = 1,nIndKPnts
                            KBlocks(:,:,i+nIndKPnts) = dconjg(KBlocks(:,:,nIndKPnts-i+1))
                        enddo
                    else
                        !Mirror the kpoints, but ignore the gamma point and BZ boundary
                        do i = 2,nIndKPnts-1
                            KBlocks(:,:,i+nIndKPnts-1) = dconjg(KBlocks(:,:,nIndKPnts-i+1))
                        enddo
                    endif
                else
                    !Odd number of kpoints
                    !This means that we can never actually sample the chemical potential
                    !tShift=T will mean that we sample the Gamma point (and it has no equivalent point)
                    !tShift=F will mean that we sample the BZ boundary (and it has no equivalent point)
                    do i = 1,nIndKPnts-1
                        if(tShift_Mesh) then
                            KBlocks(:,:,i+nIndKPnts) = dconjg(KBlocks(:,:,nIndKPnts-i))
                        else
                            KBlocks(:,:,i+nIndKPnts) = dconjg(KBlocks(:,:,nIndKPnts-i+1))
                        endif
                    enddo
                endif
            endif

        endif

    end subroutine LatParams_to_KBlocks

    !This will set FreqPoints and Weights. nFitPoints is the size of FreqPoints/Weights, 
    !and may refer to real or im axis based on the tFitRealFreq flag.
    subroutine SetFreqPoints(nFreq_Re,nFreq_Im,nFitPoints)
        implicit none
        integer, intent(out) :: nFreq_Re,nFreq_Im
        integer, intent(out) :: nFitPoints
        logical :: tUsePoswFreqPoints
        integer :: OmegaVal,i,j
        real(dp) :: Omega,LowFreq,HighFreq
        real(dp), allocatable :: FreqPoints_temp(:),Weights_temp(:)
        character(len=*), parameter :: t_r='SetFreqPoints'

        if(tFitRealFreq.and.tFitPoints_Legendre) then
            call stop_all(t_r,'Trying to fit the real frequency spectrum to Legendre points?')
        endif
        if(.not.tHalfFill.and.tUsePoswFreqPoints) then
            call stop_all(t_r,'Trying to just optimize the positive frequency points, but not half filled system')
        endif
!        if(tHalfFill.and.(nImp.eq.1)) then
!            tUsePoswFreqPoints = .true. !For ph symmetric systems, this seems to make no difference (but is faster!)
!        else
            tUsePoswFreqPoints = .false.
!        endif

        nFreq_Re = 0
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,.false.)
            if(OmegaVal.lt.0) exit
            nFreq_Re = nFreq_Re + 1
        enddo
        
        if(tFitPoints_Legendre) then
            !Calculate the Legendre polynomials for Legendre integration points and weights
            allocate(FreqPoints(nFreqPoints))
            allocate(Weights(nFreqPoints))
            write(6,"(A)") "Calculating Legendre frequency points..."
            LowFreq = -1.0_dp
            HighFreq = 1.0_dp
            call gauleg(LowFreq,HighFreq,FreqPoints,Weights,nFreqPoints)
            call ScaleFreqPointsWeights(nFreqPoints,FreqPoints,Weights)
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
                        Weights_temp(j) = Weights(i)
                    endif
                enddo
                deallocate(FreqPoints,Weights)
                nFreqPoints = j
                allocate(FreqPoints(nFreqPoints))
                allocate(Weights(nFreqPoints))
                FreqPoints(:) = FreqPoints_temp(:)
                Weights(:) = Weights_temp(:)
                deallocate(FreqPoints_temp,Weights_temp)
                write(6,"(A,I7)") "New number of frequency points: ",nFreqPoints
            endif
            nFreq_Im = nFreqPoints
            write(6,"(A)") "Frequency integration points and weights: "
            do i = 1,nFreqPoints
                write(6,"(I9,2G20.13)") i,FreqPoints(i),Weights(i)
            enddo
            !Which axis do we want to do the fitting on?
            if(tFitRealFreq) then
                nFitPoints = nFreq_Re   
                tCalcRealSpectrum = .false.
                tFitMatAxis = .false.
            else
                nFitPoints = nFreq_Im
                tFitMatAxis = .true.
            endif
        
        else
            nFreq_Im = 0
            OmegaVal = 0
            do while(.true.)
                call GetNextOmega(Omega,OmegaVal,.true.)
                if(OmegaVal.lt.0) exit
                nFreq_Im = nFreq_Im + 1
            enddo
            !Which axis do we want to do the fitting on?
            if(tFitRealFreq) then
                nFitPoints = nFreq_Re   
                tCalcRealSpectrum = .false.
                tFitMatAxis = .false.
            else
                nFitPoints = nFreq_Im
                tFitMatAxis = .true.
            endif
            allocate(FreqPoints(nFitPoints))
            allocate(Weights(nFitPoints))
            OmegaVal = 0
            do while(.true.)
                call GetNextOmega(Omega,OmegaVal,tFitMatAxis)
                if(OmegaVal.lt.0) exit
                FreqPoints(OmegaVal) = Omega
                Weights(OmegaVal) = one
            enddo
        endif

        !Initially, just see if we can fit the two different Matsubara spectral functions
        write(6,*) "Number of frequency points to fit across: ",nFitPoints
        call flush(6)

    end subroutine SetFreqPoints

    !Write out the isotropic average, and individual elements of a dynamic function in the impurity space.
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
        integer :: iunit,i,j,k,imp1,imp2
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
        if(tMatbrAxis_) tCheckCausal_ = .false. !The imaginary part won't be negative on the matsubara axis (fn is antisym)

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
                call GetNextOmega(Omega,i,tMatbrAxis_,nFreqPoints=n,FreqPoints=FreqPoints)
            else
                call GetNextOmega(Omega,i,tMatbrAxis_)
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
                            if(abs(Func(k,j,i)-dconjg(Func(j,k,i))).gt.1.0e-5) then
                                write(6,*) "While writing file: ",filename
                                if(present(tag)) write(6,*) "Filename extension: ",tag
                                write(6,*) "Element of function: ",k,j
                                write(6,*) "Frequency point: ",Omega
                                write(6,*) "Values: ",Func(k,j,i),Func(j,k,i),abs(Func(k,j,i)-dconjg(Func(j,k,i)))
                                if(tWarn_) then
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
                write(iunit,"(3G25.10)",advance='no') Omega,real(IsoAv,dp),aimag(IsoAv)
                do imp1=1,nImp
                    do imp2=1,nImp
                        write(iunit,"(2G25.10)",advance='no') real(Func(imp2,imp1,i),dp),aimag(Func(imp2,imp1,i))
                    enddo
                enddo
                write(iunit,*) 
            endif
        enddo
!        write(6,"(A,A)") "For function written to file: ",filename
        if(tMatbrAxis_) then
            write(6,"(A,A,F17.10)") trim(filename)," Total approximate spectral weight on Matsubara axis: ",SpectralWeight
        else
            write(6,"(A,A,F17.10)") trim(filename)," Total approximate spectral weight on real axis: ",SpectralWeight
        endif
        close(iunit)

    end subroutine writedynamicfunction
            
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
    
    !This function inverts a local greens function matrix, allowing it to be non-hermitian
    !Function is sent in as the normal greens function, n, nImp x nImp matrices, and the inverse
    !for each frequency is returned
    subroutine InvertLocalNonHermFunc(n,InvGF)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(inout) :: InvGF(nImp,nImp,n)
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: LVec(:,:),RVec(:,:),W_Vals(:),cWork(:)
        complex(dp) :: ztemp2(nImp,nImp),IGF(nImp,nImp)
        integer :: i,lwork,info,j
        character(len=*), parameter :: t_r='InvertLocalNonHermFunc'

        if(.true.) then
            do i = 1,n
                call mat_inv(InvGF(:,:,i),IGF,nImp)
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

    end subroutine InvertLocalNonHermFunc
    
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

    !Rather than getting the next Omega value, get the current one for the 'i'th
    !index through the list
    pure subroutine GetOmega(Omega,i,tMatbrAxis,FreqPoints)
        implicit none
        real(dp), intent(out) :: Omega
        integer, intent(in) :: i
        logical, intent(in) :: tMatbrAxis
        real(dp), intent(in), optional :: FreqPoints(:)
        
        if(present(FreqPoints)) then
            Omega = FreqPoints(i)
        else
            if(tMatbrAxis) then
                Omega = Start_Omega_Im + (i-1)*Omega_Step_Im
            else
                Omega = Start_Omega + (i-1)*Omega_Step
            endif
        endif

    end subroutine GetOmega
    
    !Get next omega value for optionally real/matsubara axis
    !OmegaVal = 0 on input initializes
    !OmegaVal = -1 on exit means we have finished
    subroutine GetNextOmega(Omega,OmegaVal,tMatbrAxis,nFreqPoints,FreqPoints)
        implicit none
        real(dp), intent(inout) :: Omega
        integer, intent(inout) :: OmegaVal
        logical, intent(in) :: tMatbrAxis
        integer, intent(in), optional :: nFreqPoints
        real(dp), intent(in), optional :: FreqPoints(:)
        logical :: tUseFreqArray
        character(len=*), parameter :: t_r='GetNextOmega'

!        write(6,*) "In Get Next Omega: ",OmegaVal
        if(present(FreqPoints)) then
            tUseFreqArray = .true.
        else
            tUseFreqArray = .false.
        endif

        if(OmegaVal.eq.0) then
            if(tUseFreqArray) then
                Omega = FreqPoints(1)
            else
                if(tMatbrAxis) then
                    Omega = Start_Omega_Im
                else
                    Omega = Start_Omega
                endif
            endif
            OmegaVal = OmegaVal + 1
!            write(6,*) "Returned Omega: ",Omega
            return
        endif

        if(tUseFreqArray) then
            if(.not.present(nFreqPoints)) call stop_all(t_r,'Error here')
            if(OmegaVal.eq.nFreqPoints) then
                OmegaVal = -1
                return
            else
                Omega = FreqPoints(OmegaVal+1)
            endif
        elseif(tMatbrAxis) then
            !Omega = Omega + Omega_Step_Im
            Omega = Start_Omega_Im + OmegaVal*Omega_Step_Im
            if((Omega.gt.(max(Start_Omega_Im,End_Omega_Im)+1.0e-6_dp))  &
                    .or.(Omega.lt.(min(Start_Omega_Im,End_Omega_Im)-1.0e-6_dp))) then
                !Hit exit criterion
                OmegaVal = -1
                return
            endif
        else
            !Omega = Omega + Omega_Step
            Omega = Start_Omega + OmegaVal*Omega_Step
            if((Omega.gt.(max(Start_Omega,End_Omega)+1.0e-6_dp))  &
                    .or.(Omega.lt.(min(Start_Omega,End_Omega)-1.0e-6_dp))) then
                !Hit exit criterion
                OmegaVal = -1
                return
            endif
        endif
        OmegaVal = OmegaVal + 1
!        write(6,*) "Returned Omega: ",Omega

    end subroutine GetNextOmega
    
    subroutine Imposesym_2(iLatParams,vars,mu)
        use sort_mod, only: sort_real
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(inout) :: vars(iLatParams)
        real(dp), intent(in) :: mu
        integer :: i,k,lWork,iErr
        real(dp) :: DistFromMu,Bands_k(nImp)
        integer :: ks
        real(dp), allocatable :: rWork(:)
        complex(dp) :: KBlock(nImp,nImp),KBlock2(nImp,nImp),cTemp(nImp,nImp)
        complex(dp), allocatable :: KBlocks(:,:,:),cWork(:)
        character(len=*), parameter :: t_r='Imposesym_2'

        if(tRealSpaceSC) then
            !The code will work if it goes through the other bit, but it will
            !have to throw away the complex coefficients at the end...
            call Imposesym_RS(iLatParams,vars,mu)
            return
        endif
            
        allocate(KBlocks(nImp,nImp,nKPnts))
        call LatParams_to_KBlocks(iLatParams,vars,mu,KBlocks)
            
        if(tImposeKSym) then
            !After fitting, add back in momentum inversion symmetry (if we are working in k-space)
            !This means that h(k) = h*(-k)
            !For shifted meshes, this is easy.
            !For Gamma-centered meshes, two k-points are only singly degenerate
            !We should not be able to be constraining any syms
            if(tRealSpaceSC) call stop_all(t_r,'K symmetry should be automatically conserved rather '   &
                &   //'than imposed with real-space lattice opt')
            if(mod(nKPnts,2).eq.0) then
                !Even number of kpoints
                if(tShift_Mesh) then
                    do i = 1,nIndKPnts  !nKPnts/2
                        KBlock(:,:) = (KBlocks(:,:,i) + dconjg(KBlocks(:,:,nKPnts-i+1))) / 2.0_dp
                        KBlocks(:,:,i) = KBlock(:,:)
                        KBlocks(:,:,nKPnts-i+1) = dconjg(KBlock(:,:))
                    enddo
                else
                    do i = 2,nIndKPnts-1
                        KBlock(:,:) = (KBlocks(:,:,i) + dconjg(KBlocks(:,:,nKPnts-i+2))) / 2.0_dp
                        KBlocks(:,:,i) = KBlock(:,:)
                        KBlocks(:,:,nKPnts-i+2) = dconjg(KBlock(:,:))
                    enddo
                endif
            else
                if(tShift_Mesh) then
                    do i = 1,nIndKPnts-1
                        KBlock(:,:) = (KBlocks(:,:,i) + dconjg(KBlocks(:,:,nKPnts-i+1))) / 2.0_dp
                        KBlocks(:,:,i) = KBlock(:,:)
                        KBlocks(:,:,nKPnts-i+1) = dconjg(KBlock(:,:))
                    enddo
                else
                    do i = 2,nIndKPnts
                        KBlock(:,:) = (KBlocks(:,:,i) + dconjg(KBlocks(:,:,nKPnts-i+2))) / 2.0_dp
                        KBlocks(:,:,i) = KBlock(:,:)
                        KBlocks(:,:,nKPnts-i+2) = dconjg(KBlock(:,:))
                    enddo
                endif
            endif
        endif

        if(tImposephsym) then
        
!            write(6,*) "*** IMPOSING PH SYMMETRY ***",mu
            if(mod(nImp,2).eq.0) then
                !Multiple of two bands per kpoint. Constrain them so that they are in pairs
!$OMP PARALLEL DO PRIVATE(KBlock,cWork,rWork,lWork,ierr,Bands_k,DistFromMu,KBlock2,cTemp)
                do k = 1,nKPnts
                    allocate(rWork(max(1,3*nImp-2)))
                    KBlock(:,:) = KBlocks(:,:,k)
                    !Diagonalize each kBlock
                    allocate(cWork(1))
                    lWork = -1
                    ierr = 0
                    Bands_k(:) = zero
                    call zheev('V','U',nImp,KBlock,nImp,Bands_k,cWork,lWork,rWork,ierr)
                    if(ierr.ne.0) call stop_all(t_r,'Error in diag')
                    lWork = int(real(cWork(1))) + 1
                    deallocate(cWork)
                    allocate(cWork(lWork))
                    call zheev('V','U',nImp,KBlock,nImp,Bands_k,cWork,lWork,rWork,ierr)
                    if(ierr.ne.0) call stop_all(t_r,'Error in diag')
                    deallocate(cWork)

                    !write(6,*) "For k point: ",k
                    !call writematrixcomp(KBlocks(:,:,k),'k-space ham',.true.)
                    !call writevector(Bands_k,'eigenvalues at kpoint')

                    !Now pair up the bands in this kblock
                    do i = 1,nImp/2
                        DistFromMu = 0.5_dp * (- Bands_k(i) + Bands_k(nImp-i+1) )
                        !write(6,*) "Dist from mu: ",DistFromMu 
                        Bands_k(i) = mu - DistFromMu
                        Bands_k(nImp-i+1) = mu + DistFromMu
                    enddo
                    !call writevector(Bands_k,'ph symmetrized eigenvalues at kpoint')

                    !Now, rotate the block back into k-space using the same eigenvectors
                    KBlock2(:,:) = zzero
                    do i = 1,nImp
                        KBlock2(i,i) = cmplx(Bands_k(i),zero,dp)
                    enddo
                    call ZGEMM('N','N',nImp,nImp,nImp,zone,KBlock,nImp,KBlock2,nImp,zzero,cTemp,nImp)
                    call ZGEMM('N','C',nImp,nImp,nImp,zone,cTemp,nImp,KBlock,nImp,zzero,KBlocks(:,:,k),nImp)
                    !call writematrixcomp(KBlocks(:,:,k),'new k-space ham',.true.)
                    deallocate(rWork)
                enddo
!$OMP END PARALLEL DO
            elseif(nImp.eq.1) then
                !If there is only one band per kpoint, then the pairs correspond to different kpoints.
                if(mod(nKPnts,2).eq.1) call stop_all(t_r,'Hardcoded for even number of kpoints for 1-band problem? FIXME')
                do i = 1,nIndKPnts/2
                    DistFromMu = 0.5_dp*( (mu - real(KBlocks(1,1,nIndKPnts-i+1)) ) + (real(KBlocks(1,1,i)) - mu))
                    KBlocks(1,1,nIndKPnts-i+1) = cmplx(mu - DistFromMu,zero,dp)
                    KBlocks(1,1,i) = cmplx(mu + DistFromMu,zero,dp)
                enddo

                !Since we are mixing kpoint, this might break k-symmetry. Reimpose it if it is being constrained
                if(tImposeKSym) then
                    !After fitting, add back in momentum inversion symmetry (if we are working in k-space)
                    !This means that e(k) = e(-k)
                    !For shifted meshes, this is easy.
                    !For Gamma-centered meshes, two k-points are only singly degenerate
                    !We should not be able to be constraining any syms
                    if(tShift_Mesh) then
                        do i = 1,nKPnts/2
                            KBlock(:,:) = (KBlocks(:,:,i) + KBlocks(:,:,nKPnts-i+1)) / 2.0_dp
                            KBlocks(:,:,i) = KBlock(:,:)
                            KBlocks(:,:,nKPnts-i+1) = KBlock(:,:)
                        enddo
                    else
                        do i = 2,nKPnts/2
                            KBlock(:,:) = (KBlocks(:,:,i) + KBlocks(:,:,nKPnts-i+2)) / 2.0_dp
                            KBlocks(:,:,i) = KBlock(:,:)
                            KBlocks(:,:,nKPnts-i+2) = KBlock(:,:)
                        enddo
                    endif
                endif
            else
                call stop_all(t_r,'Cannot impose ph symmetry if nImp != 1 or multiple of 2')
            endif
        endif
            
        call KBlocks_to_LatParams(iLatParams,vars,KBlocks)
        deallocate(KBlocks)

    end subroutine Imposesym_2
            
    !Do diag in real space to avoid getting a complex hamiltonian on back rotation
    !Diagonalize real-space hamiltonian
    subroutine Imposesym_RS(iLatParams,vars,mu)
        integer, intent(in) :: iLatParams
        real(dp), intent(inout) :: vars(iLatParams)
        real(dp), intent(in) :: mu
        complex(dp), allocatable :: ham_comp(:,:)
        real(dp), allocatable :: ham_real(:,:),vals(:),work(:),ham_real_2(:,:),temp(:,:)
        integer :: i,j,info,lWork
        real(dp) :: DistFromMu
        character(len=*), parameter :: t_r='Imposesym_RS'

        allocate(ham_comp(nSites,nSites))

        call LatParams_to_ham(iLatParams,vars,mu,ham_comp)

        allocate(ham_real(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                ham_real(j,i) = real(ham_comp(j,i),dp)
            enddo
        enddo

        !Diagonalize
        allocate(vals(nSites))
        vals(:) = zero
        allocate(work(1))
        lWork = -1
        info = 0
        call dsyev('V','L',nSites,ham_real,nSites,vals,work,lwork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','L',nSites,ham_real,nSites,vals,work,lwork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        !Organise pairs of orbitals
        do i = 1,nSites/2
            DistFromMu = 0.5_dp * (- vals(i) + vals(nSites-i+1) )
            !write(6,*) "Dist from mu: ",DistFromMu 
            vals(i) = mu - DistFromMu
            vals(nSites-i+1) = mu + DistFromMu
        enddo

        !Now, rotate the block back into real space using the same eigenvectors
        allocate(ham_real_2(nSites,nSites))
        ham_real_2(:,:) = zzero
        do i = 1,nSites
            ham_real_2(i,i) = vals(i)
        enddo
        deallocate(vals)
        allocate(temp(nSites,nSites))
        call DGEMM('N','N',nSites,nSites,nSites,one,ham_real,nSites,ham_real_2,nSites,zero,temp,nSites)
        call DGEMM('N','T',nSites,nSites,nSites,one,temp,nSites,ham_real,nSites,zero,ham_real_2,nSites)
        deallocate(temp,ham_real)
        !ham_real_2 now contains the new hamiltonian

        !Do not update the local part of the hamiltonians if it has changed
        !Wants to take a complex hamiltonian
        ham_comp(:,:) = zzero
        do i = 1,nSites
            do j = 1,nSites
                ham_comp(j,i) = cmplx(ham_real_2(j,i),zero,dp)
            enddo
        enddo
        call ham_to_LatParams(iLatParams,vars,ham_comp,tUpdateLocalHam=.false.)
        deallocate(ham_real_2,ham_comp)

    end subroutine Imposesym_RS
    
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

end module SelfConsistentUtils
