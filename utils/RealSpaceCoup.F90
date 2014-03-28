!Code to test how to include real-space one-electron couplings to the lattice hamiltonian while preserving k-point symmetry.
!Also, test for derivatives of the greens function wrt these derivatives.
program RealSpaceCoup
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)   !For double precision real numbers

    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
    complex(dp), parameter :: zzero = cmplx(0.0_dp,0.0_dp,dp)
    complex(dp), parameter :: zone = cmplx(1.0_dp,0.0_dp,dp)
    real(dp), parameter :: pi = 3.1415926535897931_dp

    integer :: nSites

    integer :: nImp,nKPnts
    logical :: tAntiPeriodic,tPeriodic,tShift_Mesh
    real(dp) :: BZVol,r2
    real(dp), allocatable :: h0(:,:),KPnts(:),RecipLattVecs(:)
    real(dp) :: PrimLattVec(1)
    complex(dp), allocatable :: RtoK_Rot(:,:)
    save

    call main()

    contains

    subroutine main()
        implicit none
        integer :: i,j,k,ind_1,ind_2
        real(dp) :: phase,r
        complex(dp), allocatable :: temp(:,:),ham_temp(:,:),Random_Corrpot(:,:)
        integer :: NonLocCoupLength
        complex(dp), allocatable :: NonLocCoup(:,:)

        write(6,*) "Enter total number of sites: "
        read(*,*) nSites
        write(6,*) "Total number of sites: ",nSites 

        write(6,*) "Enter number of impurity sites: "
        read(*,*) nImp
        write(6,*) "Number of impurity sites: ",nImp

        if(mod(nSites,2).ne.0) stop 'must have even number of sites'
        if(mod(nSites,nImp).ne.0) stop 'must have an integer number of impurity replicas through space'

        if(mod(nSites/2,2).eq.0) then
            tAntiPeriodic = .true.
            tPeriodic = .false.
            write(6,*) "Antiperiodic boundary conditions included"
        else
            tAntiPeriodic = .false.
            tPeriodic = .true.
            write(6,*) "Periodic boundary conditions included"
        endif

        !First, set up the non-interacting one-electron hamiltonian
        allocate(h0(nSites,nSites))
        h0 = 0.0_dp
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

        !Now, construct k-space
        nKPnts = nSites/nImp
        write(6,*) "Number of kpoints: ",nKPnts
        allocate(KPnts(nKPnts))
        allocate(RecipLattVecs(1))
        KPnts(:) = zero
        RecipLattVecs(:) = zero
        if(mod(nKPnts,2).eq.0) then
            if(tPeriodic) then
                !We want to use a Gamma centered mesh
                tShift_Mesh = .false.
            else
                !We want a Monkhort-Pack mesh
                tShift_Mesh = .true.
            endif
        else
            !Odd number of k-points
            if(tPeriodic) then
                !We want to use a Gamma centered mesh
                tShift_Mesh = .true. 
            else
                !Monkhorst-Pack mesh
                tShift_Mesh = .false.
            endif
        endif
        if(tShift_Mesh) then
            write(6,"(A)") "Using a Monkhorst-pack kpoint mesh of 1st BZ - no gamma point"
        else
            write(6,"(A)") "Using a gamma-centered kpoint mesh of 1st BZ - gamma point included"
        endif

        BZVol = 2.0_dp*pi/real(nImp,dp)
        !write(6,"(A,G21.14)") "Brillouin zone volume: ",BZVol

        RecipLattVecs(1) = 2.0_dp*pi/real(nImp,dp)
        !Just use equally spaced mesh starting at -pi/SS_Period, and working our way across
        do k = 1,nKPnts
            KPnts(k) = -RecipLattVecs(1)/2.0_dp + (k-1)*RecipLattVecs(1)/nKPnts
            if(tShift_Mesh) then
                !Shift kpoint mesh by half the kpoint spacing
                KPnts(k) = KPnts(k) + RecipLattVecs(1)/(2.0_dp*real(nKPnts,dp))
            endif
        enddo
            
        !Setup rotation matrix from site basis to k-space
        !First index r, second k
        allocate(RtoK_Rot(nSites,nSites))
        RtoK_Rot(:,:) = zzero
        !Run though all kpoints
        do k = 1,nKPnts
            !Construct rotation
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            do i = 1,nSites
                PrimLattVec(1) = real(i-1,dp)   !The real-space translation to this site
                phase = KPnts(k)*PrimLattVec(1)
                !phase = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec,1)
                RtoK_Rot(i,ind_1+mod(i,nImp)) = exp(dcmplx(zero,phase))/sqrt(real(nKPnts,dp))
            enddo
        enddo

        !Now, is RtoK_Rot unitary?
        allocate(temp(nSites,nSites))
        !Check unitarity of matrix
        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,RtoK_Rot,nSites,zzero,temp,nSites) 
        do i = 1,nSites
            do j = 1,nSites
                if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                    write(6,*) "i,j: ",i,j
                    call writematrixcomp(temp,'Identity?',.true.)
                    stop 'Rotation matrix not unitary'
                elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                    write(6,*) "i,j: ",i,j
                    call writematrixcomp(temp,'Identity?',.true.)
                    stop 'Rotation matrix not unitary 2'
                endif
            enddo
        enddo
        !Try other way...
        call ZGEMM('N','C',nSites,nSites,nSites,zone,RtoK_Rot,nSites,RtoK_Rot,nSites,zzero,temp,nSites) 
        do i = 1,nSites
            do j = 1,nSites
                if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                    write(6,*) "i,j: ",i,j
                    call writematrixcomp(temp,'Identity?',.true.)
                    stop 'Rotation matrix not unitary'
                elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                    write(6,*) "i,j: ",i,j
                    call writematrixcomp(temp,'Identity?',.true.)
                    stop 'Rotation matrix not unitary 2'
                endif
            enddo
        enddo
        write(6,*) "Rotation matrix unitary... :)"
        deallocate(temp)

        !Right, just as a sanity check, see if this rotation block diagonalizes the hopping matrix
        allocate(ham_temp(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                ham_temp(j,i) = cmplx(h0(j,i),zero,dp)
            enddo
        enddo
        !Add a random correlation potential
        allocate(Random_Corrpot(nImp,nImp))
        Random_CorrPot(:,:) = zzero
        do i = 1,nImp
            do j = i,nImp
                if(i.eq.j) then
                    call random_number(r)
                    Random_CorrPot(j,i) = cmplx(r,zero,dp)
                else
                    call random_number(r)
                    call random_number(r2)
                    Random_CorrPot(j,i) = cmplx(r,r2,dp)
                endif
            enddo
        enddo
        call MakeBlockHermitian(Random_CorrPot,nImp)
        call add_localpot_comp_inplace(ham_temp,Random_CorrPot,tAdd=.true.)
        allocate(temp(nSites,nSites))
        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham_temp,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RtoK_Rot,nSites,zzero,ham_temp,nSites)
        !ham_temp is now in kspace
        !Zero the diagonal blocks
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            ham_temp(ind_1:ind_2,ind_1:ind_2) = zzero
        enddo
        !Now check to see if there are any non-zero elements. If so, the unitary rotation is not correct
        do i = 1,nSites
            do j = 1,nSites
                if(abs(ham_temp(j,i)).gt.1.0e-7_dp) then
                    write(6,*) "i,j: ",j,i
                    write(6,*) "ham in kspace: ",ham_temp(j,i)
                    stop 'kspace rotations not correctly set up. ' &
                        //'One-electron hamiltonian is not block diagonal in kspace'
                endif
            enddo
        enddo
        deallocate(temp,Random_Corrpot,ham_temp)
        write(6,"(A)") "k-space rotations correctly block diagonalize general one-electron matrix"

        !Now, add arbitrary non-local couplings in real space
        allocate(ham_temp(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                ham_temp(j,i) = cmplx(h0(j,i),zero,dp)
            enddo
        enddo
        !Add a random, real correlation potential
        allocate(Random_Corrpot(nImp,nImp))
        Random_CorrPot(:,:) = zzero
        do i = 1,nImp
            do j = i,nImp
                if(i.eq.j) then
                    call random_number(r)
                    Random_CorrPot(j,i) = cmplx(r,zero,dp)
                else
                    call random_number(r)
!                    call random_number(r2)
                    Random_CorrPot(j,i) = cmplx(r,zero,dp)
                endif
            enddo
        enddo
        call MakeBlockHermitian(Random_CorrPot,nImp)
        call add_localpot_comp_inplace(ham_temp,Random_CorrPot,tAdd=.true.)
        
        !Now, add non-local terms, which retain periodicity
        if(mod(nSites-nImp,2).ne.0) then
            NonLocCoupLength = (nSites-nImp-1)/2
        else
            NonLocCoupLength = (nSites-nImp)/2
        endif

        write(6,*) "Non local hamiltonian coupling length: ",NonLocCoupLength

        allocate(NonLocCoup(nImp,NonLocCoupLength))
        NonLocCoup(:,:) = zzero
        do j = 1,NonLocCoupLength
            do i = 1,nImp
                call random_number(r)
                NonLocCoup(i,j) = cmplx(r,zero,dp)
            enddo
        enddo

        call WriteMatrixcomptoreal(ham_temp,'Before non-loc coups',.true.)

        call Add_Nonlocal_comp_inplace(ham_temp,NonLocCoup,NonLocCoupLength,tAdd=.true.)

        call WriteMatrixcomptoreal(ham_temp,'After non-loc coups',.true.)

        !Now, check whether it is still k-space kosher
        !Rotate ham_temp into k-space
        allocate(temp(nSites,nSites))
        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham_temp,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RtoK_Rot,nSites,zzero,ham_temp,nSites)
        !ham_temp is now in kspace
        !Zero the diagonal blocks
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            ham_temp(ind_1:ind_2,ind_1:ind_2) = zzero
        enddo
        !Now check to see if there are any non-zero elements. If so, the unitary rotation is not correct
        do i = 1,nSites
            do j = 1,nSites
                if(abs(ham_temp(j,i)).gt.1.0e-7_dp) then
                    write(6,*) "i,j: ",j,i
                    write(6,*) "ham in kspace: ",ham_temp(j,i)
                    stop 'kspace rotations not correctly set up. ' &
                        //'One-electron hamiltonian is not block diagonal in kspace'
                endif
            enddo
        enddo
        deallocate(temp,Random_Corrpot,ham_temp)
        write(6,"(A)") "k-space rotations correctly block diagonalize one-electron matrix with non-local interactions"

    end subroutine main

    !Given the lower triangle, make this matrix hermitian
    subroutine MakeBlockHermitian(Block,length)
        implicit none
        integer, intent(in) :: length
        complex(dp), intent(inout) :: Block(length,length)
        integer :: i,j

        do i = 1,nImp
            do j = i+1,nImp
                Block(i,j) = dconjg(Block(j,i))
            enddo
        enddo

    end subroutine MakeBlockHermitian

    subroutine add_Nonlocal_comp_inplace(ham,NonLocCoup,NonLocCoupLength,tAdd)
        implicit none
        integer, intent(in) :: NonLocCoupLength
        complex(dp), intent(inout) :: ham(nSites,nSites)
        complex(dp), intent(in) :: NonLocCoup(nImp,NonLocCoupLength)
        logical , intent(in), optional :: tAdd
        real(dp) :: phase
        integer :: k,i,j,ind_1,ind_2,ind
            
        if(tPeriodic) then
            phase = 1.0_dp
        elseif(tAntiPeriodic) then
            phase = -1.0_dp
        endif

        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k

!            write(6,*) "k: ",k
!            write(6,*) "ind_1: ",ind_1
!            write(6,*) "ind_2: ",ind_2
            do i = 1,NonLocCoupLength

                !Fill up coupling to the right
                ind = 0
                do j = ind_2+1,ind_2+NonLocCoupLength

                    ind = ind + 1
                    if(j.le.nSites) then
                        ham(ind_1:ind_2,j) = NonLocCoup(:,ind)
                    else
                        !Wrap around matrix with appropriate boundary conditions
                        ham(ind_1:ind_2,j-nSites) = NonLocCoup(:,ind)*phase
                    endif

                enddo
                if(ind.ne.NonLocCoupLength) then
                    write(6,*) "NonLocCoupLength: ",NonLocCoupLength
                    write(6,*) "ind: ",ind
                    stop 'error in indexing'
                endif

                !Fill up coupling to the left
                ind = 0
                do j = ind_1-1,ind_1-NonLocCoupLength,-1

                    ind = ind + 1
                    if(j.ge.1) then
                        ham(ind_1:ind_2,j) = NonLocCoup(:,ind)
                    else
                        !Wrap around matrix with appropriate boundary conditions
                        ham(ind_1:ind_2,j+nSites) = NonLocCoup(:,ind)*phase
                    endif
                enddo
                if(ind.ne.NonLocCoupLength) then
                    write(6,*) "NonLocCoupLength: ",NonLocCoupLength
                    write(6,*) "ind: ",ind
                    stop 'Error in indexing'
                endif
            enddo
        enddo

    end subroutine add_Nonlocal_comp_inplace
    
    !Add a *complex* local potential striped across a *real* hamiltonian 
    !(only possible with translational invariance)
    !If tAdd, then the correlation potential is added to the core potential, otherwise it is subtracted
    subroutine add_localpot_comp_inplace(core_v,CorrPot,tAdd)
        implicit none
        complex(dp) , intent(inout) :: core_v(nSites,nSites)
        complex(dp) , intent(in) :: CorrPot(nImp,nImp)
        logical , intent(in), optional :: tAdd
        integer :: i,j,k,a,b
        logical :: tAdd_

        if(present(tAdd)) then
            tAdd_ = tAdd
        else
            tAdd_ = .true.
        endif

        !Construct new hamiltonian which is block diagonal in the local potential
        do k=0,(nSites/nImp)-1
            do i=1,nImp
                do j=1,nImp
                    if(tAdd_) then
                        core_v((k*nImp)+i,(k*nImp)+j)=core_v((k*nImp)+i,(k*nImp)+j) + CorrPot(i,j)
                    else
                        core_v((k*nImp)+i,(k*nImp)+j)=core_v((k*nImp)+i,(k*nImp)+j) - CorrPot(i,j)
                    endif
                enddo
            enddo
        enddo
    end subroutine add_localpot_comp_inplace
    
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
                    write(6,"(2G18.7)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,2G18.7)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrixcomp

    subroutine WriteMatrixcomptoreal(mat,matname,tOneLine)
        implicit none
        complex(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        integer :: i,j
        logical :: tOneLine

        write(6,*) "Writing out real part of matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(G18.7)",advance='no') real(mat(i,j),dp)
                else
                    write(6,"(I6,2G18.7)") i,j,real(mat(i,j),dp)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrixcomptoreal

end program RealSpaceCoup
