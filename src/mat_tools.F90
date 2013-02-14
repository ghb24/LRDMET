!tools for matrix construction, manipulation and HF solver for DMET
module mat_tools
    use const
    use errors, only: stop_all, warning
    use globals
    implicit none

    contains
    
    !Make mean-field real-space hubbard matrix
    subroutine make_hop_mat()
        implicit none
        integer :: i
        character(len=*), parameter :: t_r='make_hop_mat'

        h0(:,:) = 0.0_dp
        if(LatticeDim.eq.1) then
            !Tridiagonal matrix
            do i=1,nSites-1
                h0(i,i+1) = -1.0_dp
                h0(i+1,i) = -1.0_dp
            enddo

            if(tPeriodic) then
                h0(1,nSites) = -1.0_dp
                h0(nSites,1) = -1.0_dp
            elseif(tAntiPeriodic) then
                h0(1,nSites) = 1.0_dp
                h0(nSites,1) = 1.0_dp
            endif
        else
            call stop_all(t_r,'Higher dimensional models not coded up')
        endif

    end subroutine make_hop_mat
    
    !Diagonalizes the core hamiltonian, and stores the allowed CS occupations in terms of
    !number of fully occupied sites.
    subroutine find_occs()
        implicit none
        real(dp), allocatable :: Work(:),h0Eigenvecs(:,:),h0Eigenvals(:)
        integer :: lWork,info,i,j,k
        integer :: iImp,nHoppingsImp,nHoppingsEnv
        character(len=*), parameter :: t_r="find_occs"

!        call writematrix(h0,'Core hamil',.false.)

        allocate(h0Eigenvecs(nSites,nSites))
        allocate(h0Eigenvals(nSites))
        h0Eigenvecs = 0.0_dp
        !First, diagonalize one-body hamiltonian
        h0Eigenvecs(:,:) = h0(:,:)
        h0Eigenvals(:) = 0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('N','U',nSites,h0Eigenvecs,nSites,h0Eigenvals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('N','U',nSites,h0Eigenvecs,nSites,h0Eigenvals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

!        call writevector(h0Eigenvals,'Core eigenvals')

        if(allocated(allowed_occs)) deallocate(allowed_occs)
        if(tHalfFill) then
            N_occs = 1
            allocate(allowed_occs(N_occs))
            allowed_occs(1) = nSites/2 
        else
            !Now count the number of allowed occupations corresponding to degenerate CS mean-field solutions
            !This counts number of totally occupied sites rather than electrons
            N_occs = 0 
            i = 1
            do while(i.le.nSites)   !i wants to always be pointing at the first element in a degenerate set
                j = i + 1
    !            write(6,*) i
                do while((j.le.nSites).and.(abs(h0Eigenvals(min(nSites,j))-h0Eigenvals(i)).lt.1.0e-10_dp))
                    !This loop will count the degenerate space
                    j = j + 1   
                enddo
                !j should not point at the end of the degenerate set which starts at i and ends at j-1
                !Now we have a CS number of sites
                if((j-1.lt.((nSites/2)+1)).and.(j-1.gt.int(0.2_dp*(nSites/2)))) then
                    !Constraint so that we only look at occupations from 0.2 x half filling to half filling
                    N_occs = N_occs + 1
                endif
                i = j    !Skip to end of degenerate set
            enddo

            write(6,"(A,I6)") "Number of allowed unique closed shell orbital configurations: ",N_occs

            allocate(allowed_occs(N_occs))
            allowed_occs(:) = 0
            !Now actually fill the array with these allowed occupations
            k = N_occs   !The slot for the occupation number
            i = 1
            do while(i.le.nSites)
                j = i + 1
                do while((j.le.nSites).and.(abs(h0Eigenvals(min(nSites,j))-h0Eigenvals(i)).lt.1.0e-10_dp))
                    !This loop will count the degenerate space
                    j = j + 1   
                enddo
                !j should now point at the end of the degenerate set which starts at i and ends at j-1
                !Now we have a CS number of sites
                if((j-1.lt.((nSites/2)+1)).and.(j-1.gt.int(0.2_dp*(nSites/2)))) then
                    !Constraint so that we only look at occupations from 0.2 x half filling to half filling
                    !Fill things up from low occupation to large, so low occupations go into final slot
                    if(tRampDownOcc) then
                        allowed_occs(k) = j - 1
                    else
                        allowed_occs(N_occs-k+1) = j - 1
                    endif
                    k = k - 1
                endif
                i = j    !Skip to end of degenerate set
            enddo

            write(6,*) "AllowedOccs: ",allowed_occs(:)

            if(k.ne.0) call stop_all(t_r,"Error in calculating occupations")
        endif
        deallocate(h0Eigenvecs,h0Eigenvals)
        if(mod(nSites,nImp).ne.0) call stop_all(t_r,'Number of sites should be factor of impurity size')

        !Count the hopping parameters available
        nHoppingsImp = 0
        nHoppingsEnv = 0
        do iImp=1,nImp
            do i=1,nImp
                !This counts hopping elements within the impurity cluster, which is defined as the first nImp sites
                if(abs(h0(i,iImp)).gt.1.0e-8_dp) nHoppingsImp = nHoppingsImp + 1
            enddo
            do i=nImp+1,nSites
                !This counts hopping elements from impurity to the environment, which is defined as the indices after the nImp sites
                if(abs(h0(i,iImp)).gt.1.0e-8_dp) nHoppingsEnv = nHoppingsEnv + 1
            enddo
        enddo
        write(6,"(A,I6,A,I6,A)") "Connections of impurity sites via t: ",nHoppingsImp," within emp, ",nHoppingsEnv," to env"

    end subroutine find_occs

    !Add the local potential striped across the core hamiltonian 
    !(only possible with translational invariance)
    !If tAdd, then the correlation potential is added to the core potential, otherwise it is subtracted
    subroutine add_localpot(core,core_v,CorrPot,tAdd)
        implicit none
        real(dp) , intent(in) :: core(nSites,nSites)
        real(dp) , intent(out) :: core_v(nSites,nSites)
        real(dp) , intent(in) :: CorrPot(nImp,nImp)
        logical , intent(in), optional :: tAdd
        real(dp) :: CorrPot_Flip(nImp,nImp)
        integer :: i,j,k,i_ind,j_ind
        logical :: tAdd_

        if(present(tAdd)) then
            tAdd_ = tAdd
        else
            tAdd_ = .true.
        endif

        if(tFlipUTiling) then
            write(6,*) "Flipping correlation potential when tiling mean-field hamiltonian"
            CorrPot_Flip(:,:) = 0.0_dp
            do i=1,nImp
                do j=1,nImp
                    j_ind = nImp - j + 1
                    i_ind = nImp - i + 1
                    CorrPot_Flip(j_ind,i_ind) = CorrPot(i,j)
                enddo
            enddo
        endif

        if(LatticeDim.eq.1) then
            !Construct new hamiltonian which is block diagonal in the local potential
            core_v = 0.0_dp
            do k=0,(nSites/nImp)-1
                if(tFlipUTiling.and.(mod(k,2).eq.1)) then
                    do i=1,nImp
                        do j=1,nImp
                            if(tAdd_) then
                                core_v((k*nImp)+i,(k*nImp)+j)=CorrPot_Flip(i,j)
                            else
                                core_v((k*nImp)+i,(k*nImp)+j)=-CorrPot_Flip(i,j)
                            endif
                        enddo
                    enddo
                else
                    do i=1,nImp
                        do j=1,nImp
                            if(tAdd_) then
                                core_v((k*nImp)+i,(k*nImp)+j)=CorrPot(i,j)
                            else
                                core_v((k*nImp)+i,(k*nImp)+j)=-CorrPot(i,j)
                            endif
                        enddo
                    enddo
                endif
            enddo

            !Add this to the original mean-field hamiltonian
            do i=1,nSites
                do j=1,nSites
                    core_v(i,j) = core_v(i,j) + core(i,j)
                enddo
            enddo
        endif

    end subroutine add_localpot

    !Run a full HF, including mean-field on-site repulsion term in the fock matrix
    !These are stored in FullHFOrbs and FullHFEnergies
    subroutine run_true_hf()
        implicit none
        real(dp) :: HEl,GetHFAntisymInt_spinorb,PDiff,fockel
        real(dp), allocatable :: Work(:),OccOrbs_HF(:,:),PMatrix_old(:,:),PMatrix(:,:)
        real(dp), allocatable :: fock(:,:),temp(:,:),h0HF(:,:)
        integer :: i,lWork,info,ex(2,2),j,nIter
        logical :: tFailedSCF
        character(len=*), parameter :: t_r='run_true_hf'

        write(6,"(A)") "Constructing full HF solution. DMET will start from core hamiltonian solution."

        !Construct fock matrix
        !The fock matrix is just the core hamiltonian (without the fitted potential) + diag(1/2 U * rdm(i,i)) on the diagonals
        allocate(fock(nSites,nSites))
        fock(:,:) = h0(:,:) !Core hamiltonian
        if(.not.tAnderson) then
            do i=1,nSites
                !Include the on-site repulsion
                fock(i,i) = fock(i,i) + U * 0.5_dp * (NEl/real(nSites))
            enddo
        elseif(tChemPot) then
            fock(1,1) = fock(1,1) - (U/2.0_dp)
        endif
        
        if(allocated(FullHFOrbs)) then
            deallocate(FullHFOrbs,FullHFEnergies)
        endif
        allocate(FullHFOrbs(nSites,nSites)) !The orbitals from the diagonalization of the fock matrix
        allocate(FullHFEnergies(nSites))
        !call writematrix(fock,'fock',.true.)

        !Now just diagonalise this fock matrix, rather than use diis
        !First, diagonalize one-body hamiltonian
        FullHFOrbs(:,:) = Fock(:,:)
        FullHFEnergies(:) = 0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        PDiff = 1.0_dp
        allocate(OccOrbs_HF(nSites,nOcc))
        OccOrbs_HF(:,:) = FullHFOrbs(:,1:nOcc)

        allocate(PMatrix_old(nSites,nSites))
        allocate(PMatrix(nSites,nSites))

        !Calculate initial trial P matrix:
        call dgemm('N','T',nSites,nSites,nOcc,1.0_dp,OccOrbs_HF,nSites,OccOrbs_HF,nSites,0.0_dp,PMatrix_old,nSites)
        !call writevector(FullHFEnergies(1:10),'Initial HF eigenvalues')

        tFailedSCF = .false.
        nIter = 0
        do while(PDiff.gt.1.0e-8_dp)
            nIter = nIter + 1
            if(nIter.gt.1000) then
                call warning(t_r,'Failed to converge SCF. Exiting calculation of true HF orbitals...')
                tFailedSCF = .true.
                exit
            endif
            FullHFOrbs(:,:) = h0(:,:)
            if(.not.tAnderson) then
                do i=1,nSites
                    FullHFOrbs(i,i) = FullHFOrbs(i,i) + PMatrix_old(i,i)*U
                enddo
            else
                !In the anderson model, we only have an on-site repulsion on the first site
                FullHFOrbs(1,1) = FullHFOrbs(1,1) + PMatrix_old(1,1)*U
            endif
            FullHFEnergies(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
        
            OccOrbs_HF(:,:) = FullHFOrbs(:,1:nOcc)

            !Create new PMatrix
            call dgemm('N','T',nSites,nSites,nOcc,1.0_dp,OccOrbs_HF,nSites,OccOrbs_HF,nSites,0.0_dp,PMatrix,nSites)

            PDiff = 0.0_dp
            do i=1,nSites
                do j=1,nSites
                    PDiff = PDiff + abs(PMatrix(i,j)-PMatrix_old(i,j))
                enddo
            enddo
            PMatrix_old(:,:) = PMatrix(:,:)
            !write(6,*) "PDiff: ",PDiff
            !call writevector(FullHFEnergies(1:10),'HF eigenvalues')
        enddo
        deallocate(PMatrix,PMatrix_old,OccOrbs_HF)

        !write(6,*) "Full HF Orbs: "
        !call writematrix(FullHFOrbs,'FullHFOrbs',.true.)
            
        if(.not.tFailedSCF) then
            write(6,*) "nOCC", nOcc
            write(6,*) "*True* Fock eigenvalues around fermi level: "
            do i=max(1,nOcc-3),nOcc
                write(6,*) FullHFEnergies(i),"*"
            enddo
            do i=nOcc+1,min(nSites,nOcc+3)
                write(6,*) FullHFEnergies(i)
            enddo

            !Convert core hamiltonian into HF basis
            allocate(temp(nSites,nSites))
            allocate(h0HF(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,h0HF,nSites)
            deallocate(temp)

            if(.true.) then
                !Generate fock eigenvalues and see if they are the same
                do i=1,nSites
                    fockel = h0HF(i,i)
                    do j=1,nel
                        ex(1,1) = j
                        ex(1,2) = i*2
                        ex(2,1) = j
                        ex(2,2) = i*2
                        HEl = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                        fockel = fockel + HEl
                    enddo
                    !write(6,*) "Fock eigenvalue calculated: ",i,fockel
                    if(abs(fockel-FullHFEnergies(i)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'HF solution not correct - fock eigenvalues do not agree')
                    endif
                enddo
            endif

            !Now calculate HF energy:
            HFEnergy = 0.0_dp
            do i=1,nOcc
                HFEnergy = HFEnergy + h0HF(i,i)*2.0_dp
            enddo
            do i=1,nel
                do j=1,nel
                    ex(1,1) = i
                    ex(1,2) = j
                    ex(2,1) = i
                    ex(2,2) = j
                    HEl = GetHFAntisymInt_spinorb(ex,FullHFOrbs)

                    HFEnergy = HFEnergy + 0.5_dp*HEl 
                enddo
            enddo
            write(6,*) "HF energy from core hamiltonian: ",HFEnergy

            HFEnergy = 0.0_dp
            do i=1,nOcc
                HFEnergy = HFEnergy + 2.0_dp*FullHFEnergies(i)
            enddo
            do i=1,nel
                do j=1,nel
                    ex(1,1) = i
                    ex(1,2) = j
                    ex(2,1) = i
                    ex(2,2) = j
                    HEl = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                    HFEnergy = HFEnergy - HEl*0.5_dp
                enddo
            enddo
            write(6,*) "HF energy from fock eigenvalues: ",HFEnergy

            deallocate(h0HF)
        endif

        deallocate(fock)

    end subroutine run_true_hf

    !Run a HF calculation on the entire system. In this case, it just consists of just diagonalizing the system rather than iterative DIIS (add later)
    subroutine run_hf(it)
        implicit none
        integer, intent(in) :: it
        real(dp), allocatable :: Work(:)
        real(dp), allocatable :: fock(:,:),OccOrbs(:,:)
        logical :: tRotateOrbs
        integer :: lWork,info,iFirst,iLast,i
        real(dp) :: ThrDeg
        !complex(dp), allocatable :: k_ham(:,:)
        character(len=*), parameter :: t_r="run_hf"

        !Construct fock matrix
        !The fock matrix is just the core hamiltonian (with the fitted potential) + diag(1/2 U * rdm(i,i)) on the diagonals
        allocate(fock(nSites,nSites))
        fock(:,:) = h0v(:,:)
                    
!        call writematrix(h0v,'h0v',.true.)

        if(allocated(HFOrbs)) then
            !If HFOrbs is already allocated, then we want to maximize the overlap between orbitals of different iterations
            tRotateOrbs = .true.
        else
            tRotateOrbs = .false.
            allocate(HFOrbs(nSites,nSites)) !The orbitals from the diagonalization of the fock matrix
        endif
        !Don't include the 2 electron terms in the mean field. Therefore fock matrix is core hamiltonian
        !The diagonal on-site repulsion will just be mopped up by the correlation potential
!        do i=1,nSites
!            !Include the on-site repulsion
!            fock(i,i) = fock(i,i) + U * 0.5_dp * (NEl/real(nSites))
!        enddo

        !Now just diagonalise this fock matrix, rather than use diis
        !First, diagonalize one-body hamiltonian
        HFOrbs(:,:) = Fock(:,:)
        HFEnergies(:) = 0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,HFOrbs,nSites,HFEnergies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,HFOrbs,nSites,HFEnergies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
            
        write(6,*) "nOCC", nOcc
        write(6,*) "Fock eigenvalues around fermi level: "
        do i=max(1,nOcc-7),nOcc
            write(6,*) HFEnergies(i),"*"
        enddo
        do i=nOcc+1,min(nSites,nOcc+7)
            write(6,*) HFEnergies(i)
        enddo
        
!        allocate(k_ham(nSites,nSites))
!        if(it.eq.1) then
!            !No correlation potential applied - periodicity is just 1
!            call Convert1DtoKSpace(fock,nSites,1,k_ham) 
!        else
!            call Convert1DtoKSpace(fock,nSites,nImp,k_ham) 
!        endif
!        deallocate(k_ham)


        if(tRotateOrbs.and.(it.ge.4)) then
            !Rotate the orbitals so that we can maximise overlap between these orbitals and the orbitals of the previous iteration
            !at the Fermi level if not a unique CS solution there.

            iFirst = max(1,nOcc-40)
            iLast = min(nSites,nOcc+40)
            ThrDeg = 1.0e-4
            do while(abs(HFEnergies(iFirst) - HFEnergies(nOcc)) .gt. ThrDeg)
                iFirst = iFirst + 1
            enddo
            do while(abs(HFEnergies(iLast) - HFEnergies(nOcc)) .gt.ThrDeg)
                iLast = iLast - 1
            enddo
            if(iLast .lt. iFirst) call stop_all(t_r,'Error in rotating orbitals')
            if((iLast - iFirst).ge.1) then
                write(6,*) "Warning: May potentially need to rotate orbitals at the fermi level for maximum overlap"
!                call stop_all(t_r,'Need to rotate orbitals at Fermi level')
            endif
        endif

        !Now calculate the density matrix from the calculation based on double occupancy of the lowest lying nOcc orbitals
        !First, extract the occupied orbitals. Since eigenvalues are ordered in increasing order, these will be the first nOcc
        allocate(OccOrbs(nSites,nOcc))
        OccOrbs(:,:) = HFOrbs(:,1:nOcc)
        !Now construct the density matrix in the original AO basis. The eigenvectors are given as AO x MO, so we want to contract out the
        !MO contributions in order to get the 1DM in the AO basis.
        call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,OccOrbs,nSites,OccOrbs,nSites,0.0_dp,MeanFieldDM,nSites)

        ChemPot = (HFEnergies(nOcc) + HFEnergies(nOcc+1))/2.0_dp  !Chemical potential is half way between HOMO and LUMO
        HLGap = HFEnergies(nOcc+1)-HFEnergies(nOcc)   !one body HOMO-LUMO Gap

        if(HLGap.lt.1.0e-6_dp) then
            write(6,"(A,G15.5,A)") "Warning. HL gap is: ",HLGap," Possible failure in assigning orbitals in degenerate set."
        endif

        deallocate(Fock,OccOrbs)

    end subroutine run_hf

    !Convert a real-space symmetric operator into a k-space operator.
    !In: The operator in real space. Size = nSuperCell x nSupercell
    !   However, it assumes periodicity, such that only the first unit cell, and its connections to the other unit cells
    !   is referenced. I.e. R_Ham(1:nUnitCell,:). In the future, this should be changed so that only this is called.
    !Out: The operator in k-space. In this, the operator is block-diagonal, with each block being of size nUnitCell x nUnitCell.
    subroutine Convert1DtoKSpace(R_Ham,nSuperCell,nUnitCell,k_Ham)
        implicit none
        integer, intent(in) :: nSuperCell,nUnitCell
        real(dp), intent(in) :: R_Ham(nSuperCell,nSuperCell)
        complex(dp), intent(out) :: k_Ham(nSuperCell,nSuperCell)
        complex(dp) :: KPntHam(nUnitCell,nUnitCell)
        integer :: k,i,j,a,b,nKpnts,cell_start,cell_end,k_start,k_end
        complex(dp) :: phase
        real(dp) :: KPnt_val
        real(dp), allocatable :: K_Vals(:),Orbs(:,:),W(:),Work(:),Vals(:),KPnts(:)
        complex(dp), allocatable :: cWork(:)
        integer :: lWork,info
        logical, parameter :: tCheck = .false.
        character(len=*), parameter :: t_r='Convert1DtoKSpace'

        if(tWriteOut) then
            write(6,*) "Converting real space operator to k-space: "
        endif
        if(tAntiPeriodic) call stop_all(t_r,'Cannot convert to k-space with anti-periodic boundary conditions')
        if(.not.tPeriodic) call stop_all(t_r,'Need periodic boundary conditions to convert to k-space')

        k_Ham(:,:) = dcmplx(0.0_dp,0.0_dp)

        !Number of kpoints = nSupercell/nUnitCell
        if(mod(nSupercell,nUnitCell).ne.0) call stop_all(t_r,'Not integer number of unit cells in supercell')
        nKpnts = nSupercell / nUnitCell
        if(mod(nKpnts,2).eq.1) then
            call stop_all(t_r,'For some reason, I am not getting hermitian operators with odd numbers of kpoints. " &
     &           //"Debug this routine.')
        endif

        allocate(KPnts(nKpnts))
        !Allocate values for k-vectors -> just 1D to start with here
!        write(6,*) "Number of k-points: ",nKpnts

!        if(mod(nKpnts,2).eq.0) then
            !Number of k-points even. Just use equally spaced mesh starting at -pi/a, and working our way across
            do k = 1,nKPnts
                KPnts(k) = -pi/real(nUnitCell,dp) + (k-1)*(2.0_dp*pi/nKpnts)/real(nUnitCell,dp)
            enddo
!        else
!            !Ensure that there is a kpoint at the Gamma point, and then equally spaced (don't sample BZ boundary)
!            do k = 1,nKPnts
!                KPnts(k) = -pi/real(nUnitCell,dp) + k*(2.0_dp*pi/(nKpnts+1))/real(nUnitCell,dp)
!            enddo
!        endif
!        call writevector(KPnts,'KPoint values')

        do k = 1,nKpnts
            !Create each block
            KPntHam(:,:) = R_Ham(1:nUnitCell,1:nUnitCell)   !The operator of the unit cell
            !write(6,*) "KPnt ",k,0,R_Ham(1:nUnitCell,1:nUnitCell)
            KPnt_val = KPnts(k)
            !write(6,*) "KPnt_val: ",KPnt_val

            do i = 1,nKpnts-1   !Real space translation lattice vectors to the (i+1)th cell
                cell_start = i*nUnitCell + 1
                cell_end = (i+1)*nUnitCell
                phase = exp(dcmplx(0.0_dp,KPnt_val*real(nUnitCell*i,dp))) !Phase between cell 1 and cell i+1
                !Add to the current kpoint, the opertor over the translated sites x by the phase change to them
                KPntHam(:,:) = KPntHam(:,:) + phase*R_Ham(1:nUnitCell,cell_start:cell_end)
                !write(6,*) "KPnt ",k,i,phase,phase*R_Ham(1:nUnitCell,cell_start:cell_end)
            enddo
                
            !TEST: Check hermiticity at all points
            do a = 1,nUnitCell
                do b = a,nUnitCell
                    if(abs(KPntHam(b,a)-conjg(KPntHam(a,b))).gt.1.0e-9_dp) then
                        write(6,*) a,b,KPntHam(a,b),KPntHam(b,a)
                        call writematrix(R_Ham(1:nUnitCell,:),'Coupling in real space',.true.)
                        call writematrix(R_Ham(:,:),'R_Ham',.true.)
                        call writematrixcomp(KPntHam,'Hamiltonian at kpoint',.true.)
                        call stop_all(t_r,'k-space operator not hermitian. Does input operator have periodicity')
                    endif
                enddo
            enddo

            !Now add this k-point to the full k-space operator
            k_start = (k-1)*nUnitCell + 1
            k_end = k*nUnitCell

            k_Ham(k_start:k_end,k_start:k_end) = KPntHam(:,:)
        enddo

        if(tWriteOut) then
            write(6,*) "Diagonal part of k-space operator: "
            do i = 1,nSuperCell
                write(6,*) i,k_Ham(i,i)
            enddo
        endif

        if(tCheck) then
            !TEST: Diagonalize real-space hamiltonian and check that eigenvalues are the same
            !We should now have the eigenvalues of the hamiltonian
            !Diagonalize the real-space hamiltonian and check that we have got this
            allocate(Orbs(nSuperCell,nSuperCell))
            Orbs(:,:) = R_Ham(:,:)
            allocate(W(nSuperCell))
            W(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','L',nSuperCell,Orbs,nSuperCell,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','L',nSuperCell,Orbs,nSuperCell,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            allocate(K_Vals(nSuperCell))
            if(nUnitCell.eq.1) then
                do i = 1,nSuperCell
                    K_Vals(i) = real(k_Ham(i,i),dp)
                enddo
            else
                lWork = max(1,2*nUnitCell-1)
                allocate(cWork(lWork))
                allocate(Work(max(1,3*nUnitCell-2)))
                allocate(Vals(nUnitCell))
                do k = 1,nKPnts
                    !Diagonalize block
                    k_start = (k-1)*nUnitCell + 1
                    k_end = k*nUnitCell
                    KPntHam(:,:) = k_Ham(k_start:k_end,k_start:k_end)
                    !Hermitian matrix diagonalization
                    call writematrixcomp(KPntHam,'KPntHam',.true.)
                    call ZHEEV('N','U',nUnitCell,KPntHam,nUnitCell,Vals,cWork,lWork,Work,info)
                    if(info.ne.0) call stop_all(t_r,'Diag Failed')
                    do j = 1,nUnitCell
                        K_Vals(k_start+j-1) = Vals(j)
                    enddo
                enddo
                deallocate(cWork,Work,Vals)
            endif

            !Sort values
            call sort_real(K_Vals,nSuperCell)

            do i = 1,nSuperCell
                !write(6,*) i,K_Vals(i),W(i)
                if(abs(K_Vals(i)-W(i)).gt.1.0e-9_dp) then
                    write(6,*) i,K_Vals(i),W(i)
                    call stop_all(t_r,'Conversion to k-space failed')
                endif
            enddo
            deallocate(K_Vals,W,Orbs)
        endif

        deallocate(KPnts)

    end subroutine Convert1DtoKSpace
    
    !The error metric used for the fitting of the vloc in order to match the RDMs
    !The error metric is not actually calculated, but can be considered as the squared sum of the elements in the
    !returned matrix. The matrix is then the gradients in each direction, and the jacobian is made numerically.
    !IN: vloc over impurity sites in triangular packed form
    !OUT: Error matrix between the systems (just difference over all embedded sys) (triangular packed)
    subroutine RDMErr(v,ErrMat_packed)
        implicit none
        real(dp), intent(in) :: v(nImpCombs)
        real(dp), intent(out) :: ErrMat_packed(EmbCombs)
        real(dp) :: ErrMat_unpacked(EmbSize,EmbSize)
        real(dp) :: MFRdm(EmbSize,EmbSize)

        call mkrdm(v,MFRdm)
        ErrMat_unpacked(:,:) = MFRdm(:,:) - HL_1RDM(:,:)
!        call writematrix(MFRdm,'MFRdm',.true.)
        call ToTriangularPacked(EmbSize,ErrMat_unpacked,ErrMat_packed) 
    end subroutine RDMErr

    !Construct rdm over the embedded system from the fock + potential v (over impurity sites)
    !This will add the v to the impurity sites on the fock matrix, before diagonalizing, and rotating back to the RDM in the embedded system
    !with the original meanfield+vloc occupation numbers
    subroutine mkrdm(v,rdm)
        implicit none
        real(dp), intent(in) :: v(nImpCombs)
        real(dp), intent(out) :: rdm(EmbSize,EmbSize)
        real(dp) :: EValues(EmbSize),EVectors(EmbSize,EmbSize)
        real(dp) :: temp(EmbSize,EmbSize),temp2(EmbSize,EmbSize)
        integer :: i

        !Take the diagonalized system from mkorb, and construct the RDM in the basis of meanfield solution,
        call mkorb(v,EValues,EVectors)

!        call writevector(Evalues,'Evalues')
!        call writematrix(EVectors,'EVectors',.true.)

        !Now transform the occupation numbers from the embedded system into the mean-field+new_vloc natural orbital embedded basis
        ! U,diag(),U^T. This is what we want to match.
        !Create diagonal matrix
        temp(:,:) = 0.0_dp
        do i=1,EmbSize
            temp(i,i) = MFEmbOccs(i)
        enddo
!        call writematrix(temp,'MFEmbOccs',.true.)
        call DGEMM('N','N',EmbSize,EmbSize,EmbSize,1.0_dp,EVectors,EmbSize,temp,EmbSize,0.0_dp,temp2,EmbSize)
        call DGEMM('N','T',EmbSize,EmbSize,EmbSize,1.0_dp,temp2,EmbSize,EVectors,EmbSize,0.0_dp,rdm,EmbSize)

    end subroutine mkrdm

    !Take a potential over the impurity system (in triangular form), and add it to the fock matrix in the embedded system (only over the impurity sites) 
    !Diagonalize this embedded system and return the eigenvalues and vectors
    subroutine mkorb(v,EValues,EVectors)
        implicit none
        real(dp), intent(in) :: v(nImpCombs)
        real(dp), intent(out) :: EValues(EmbSize),EVectors(EmbSize,EmbSize)
        real(dp) :: vloc_unpacked(nImp,nImp)
        real(dp), allocatable :: work(:)
        integer :: lWork,info
        character(len=*), parameter :: t_r='mkorb'

        !Unpack potential
        call FromTriangularPacked(nImp,v,vloc_unpacked)
        !Add the potential over the impurity sites to the fock matrix over the entire embedded system
        EVectors(:,:) = Emb_Fock(:,:)
        EVectors(1:nImp,1:nImp) = EVectors(1:nImp,1:nImp) + vloc_unpacked(:,:)

        !EVectors now contains Fock + the vlocal over the impurity sites
        !Now diagonalize
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',EmbSize,EVectors,EmbSize,EValues,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',EmbSize,EVectors,EmbSize,EValues,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

    end subroutine mkorb

    !Routine to triangular pack a matrix and return it in 'Packed'.
    pure subroutine ToTriangularPacked(Length,Unpacked,Packed)
        implicit none
        integer, intent(in) :: Length
        real(dp) , intent(in) :: Unpacked(Length,Length)
        real(dp) , intent(out) :: Packed((Length*(Length+1))/2)
        integer :: i,j,k

        Packed(:) = 0.0_dp
        k=1
        do i=1,Length
            do j=1,i
                Packed(k) = Unpacked(i,j)
                k=k+1
            enddo
        enddo

    end subroutine ToTriangularPacked

    pure subroutine FromTriangularPacked(Length,Packed,Unpacked)
        implicit none
        integer, intent(in) :: Length
        real(dp) , intent(out) :: Unpacked(Length,Length)
        real(dp) , intent(in) :: Packed((Length*(Length+1))/2)
        integer :: i,j,k

        k=1
        do i=1,Length
            do j=1,i
                Unpacked(i,j) = Packed(k)
                Unpacked(j,i) = Packed(k)
                k=k+1
            enddo
        enddo

    end subroutine FromTriangularPacked

    subroutine WriteMatrixcomp(mat,matname,tOneLine)
        implicit none
        complex(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        integer :: i,j
        logical :: tOneLine

        write(6,*) "Writing out matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(2G25.10)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,2G25.10)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrixcomp

    subroutine WriteMatrix(mat,matname,tOneLine)
        implicit none
        real(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        integer :: i,j
        logical :: tOneLine

        write(6,*) "Writing out matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(G25.10)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,G25.10)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrix

    subroutine WriteVector(vec,vecname)
        implicit none
        real(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVector
    
    subroutine WriteVectorcomp(vec,vecname)
        implicit none
        complex(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(2G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVectorcomp
    
    subroutine WriteVectorInt(vec,vecname)
        implicit none
        integer, intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(I12)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVectorInt

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*****************************************************************************
    function znrm2 ( n, x, incx )
    !*****************************************************************************
    !
    !! SCNRM2 returns the euclidean norm of a complex(kind=dp) vector.
    !
    !
    !  Discussion:
    !
    !    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
    !            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vector.
    !
    !    Input, complex(kind=dp) X(*), the vector.
    !
    !    Input, integer INCX, the increment between successive entries of X.
    !
    !    Output, real(kind=dp) SCNRM2, the norm of the vector.
    !
      implicit none
    !
      integer(ip), intent(in)      :: incx
      integer(ip)                  :: ix
      integer(ip), intent(in)      :: n
      real(kind=dp)                :: norm
     !real(kind=dp), parameter     :: one = 1.0_dp
      real(kind=dp)                :: scale
      real(kind=dp)                :: znrm2
      real(kind=dp)                :: ssq
      real(kind=dp)                :: temp
      complex(kind=dp), intent(in) :: x(*)
     
    !
      if ( n < 1 .or. incx < 1 ) then

        norm  = zero

      else

        scale = zero
        ssq = one

        do ix = 1, 1 + ( n - 1 ) * incx, incx
          if ( real(x(ix), dp) /= zero ) then  
            temp = abs ( real(x(ix), dp) )   
            if ( scale < temp ) then
              ssq = one + ssq * ( scale / temp )**2
              scale = temp
            else
              ssq = ssq + ( temp / scale )**2
            end if
          end if

          if ( aimag ( x(ix) ) /= zero ) then
            temp = abs ( aimag ( x(ix) ) )
            if ( scale < temp ) then
              ssq = one + ssq * ( scale / temp )**2
              scale = temp
            else
              ssq = ssq + ( temp / scale )**2
            end if

          end if

        end do

        norm  = scale * sqrt ( ssq )

      end if

      znrm2 = norm

      return
    end function znrm2

end module mat_tools
