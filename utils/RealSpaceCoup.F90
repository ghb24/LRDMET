!Code to test how to include real-space one-electron couplings to the lattice hamiltonian while preserving k-point symmetry.
!Also, test for derivatives of the greens function wrt these derivatives.
program RealSpaceCoup
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)   !For double precision real numbers

    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
    complex(dp), parameter :: zzero = cmplx(0.0_dp,0.0_dp,dp)
    complex(dp), parameter :: zone = cmplx(1.0_dp,0.0_dp,dp)
    real(dp), parameter :: pi = 3.1415926535897931_dp

    integer, parameter :: nFreqPoints = 100

    integer :: nSites

    integer :: nImp,nKPnts
    logical :: tAntiPeriodic,tPeriodic,tShift_Mesh
    real(dp) :: BZVol,r2
    real(dp), allocatable :: h0(:,:),KPnts(:),RecipLattVecs(:)
    real(dp) :: PrimLattVec(1)
    complex(dp), allocatable :: RtoK_Rot(:,:)
    complex(dp), allocatable :: FreqPoints(:)
    complex(dp), allocatable :: Derivatives(:,:,:,:,:)
    save

    call main()

    contains

    subroutine main()
        implicit none
        integer :: i,j,k,ind_1,ind_2,w,ii,jj
        real(dp) :: phase,r,diff,r2,err
        complex(dp), allocatable :: temp(:,:),ham_temp(:,:),Random_Corrpot(:,:)
        integer :: NonLocCoupLength
        complex(dp), allocatable :: NonLocCoup(:,:),ham_k(:,:)
        complex(dp), allocatable :: GF_k(:,:,:),GF_r(:,:,:),InvHam(:,:),InvHam_2(:,:)
        complex(dp), allocatable :: ham_temp_2(:,:),ham_temp_3(:,:),num(:,:)
        complex(dp), allocatable :: NonLocCoup_tmp(:,:)

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
                call random_number(r2)
                if(r2.gt.0.5_dp) then
                    NonLocCoup(i,j) = cmplx(-r,zero,dp)
                else
                    NonLocCoup(i,j) = cmplx(r,zero,dp)
                endif
            enddo
        enddo

        !call WriteMatrixcomptoreal(ham_temp,'Before non-loc coups',.true.)

        call Add_Nonlocal_comp_inplace(ham_temp,NonLocCoup,NonLocCoupLength,tAdd=.true.)

        !call WriteMatrixcomptoreal(ham_temp,'After non-loc coups',.true.)

        !Now, check whether it is still k-space kosher
        !Rotate ham_temp into k-space
        allocate(temp(nSites,nSites))
        allocate(ham_k(nSites,nSites))
        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham_temp,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RtoK_Rot,nSites,zzero,ham_k,nSites)
        !ham_temp is now in kspace
        !call WriteMatrixcomp(ham_temp,'k-space ham after non-loc coups',.true.)
        !Zero the diagonal blocks
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            ham_k(ind_1:ind_2,ind_1:ind_2) = zzero
        enddo
        !Now check to see if there are any non-zero elements. If so, the unitary rotation is not correct
        do i = 1,nSites
            do j = 1,nSites
                if(abs(ham_k(j,i)).gt.1.0e-7_dp) then
                    write(6,*) "i,j: ",j,i
                    write(6,*) "ham in kspace: ",ham_k(j,i)
                    stop 'kspace rotations not correctly set up. ' &
                        //'One-electron hamiltonian is not block diagonal in kspace'
                endif
            enddo
        enddo
        write(6,"(A)") "k-space rotations correctly block diagonalize one-electron matrix with non-local interactions"

        write(6,*) ""

        !Check real-space vs k-space definitions of the greens function
        !First, stupidly remake k-space ham
        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham_temp,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RtoK_Rot,nSites,zzero,ham_k,nSites)
        !Set frequency (imaginary) points
        allocate(FreqPoints(nFreqPoints))
        FreqPoints(:) = zzero
        FreqPoints(1) = cmplx(0.0_dp,-2.0_dp,dp)
        do i = 2,nFreqPoints
            FreqPoints(i) = FreqPoints(i-1) + cmplx(0.0_dp,4.0_dp/nFreqPoints,dp)
        enddo

        allocate(GF_k(nImp,nImp,nFreqPoints))
        allocate(GF_r(nImp,nImp,nFreqPoints))
        GF_k(:,:,:) = zzero
        GF_r(:,:,:) = zzero

        allocate(InvHam(nImp,nImp))
        allocate(InvHam_2(nImp,nImp))
        allocate(num(nImp,nImp))

        !k-space definition
        do i = 1,nFreqPoints
            do k = 1,nKPnts
                ind_1 = ((k-1)*nImp) + 1
                ind_2 = nImp*k

                InvHam(:,:) = - ham_k(ind_1:ind_2,ind_1:ind_2)
                do j = 1,nImp
                    InvHam(j,j) = InvHam(j,j) + FreqPoints(i)
                enddo
                call mat_inv(InvHam,InvHam_2)
                InvHam(:,:) = RtoK_Rot(1:nImp,ind_1:ind_2)
                call ZGEMM('N','C',nImp,nImp,nImp,zone,InvHam_2,nImp,InvHam,nImp,zzero,num,nImp)
                call ZGEMM('N','N',nImp,nImp,nImp,zone,InvHam,nImp,num,nImp,zone,GF_k(:,:,i),nImp)
            enddo
        enddo
        deallocate(InvHam,InvHam_2,num)

        !r-space definition
        allocate(ham_temp_2(nSites,nSites))
        allocate(ham_temp_3(nSites,nSites))
        do i = 1,nFreqPoints

            ham_temp_2(:,:) = - ham_temp(:,:)
            do j = 1,nSites
                ham_temp_2(j,j) = ham_temp_2(j,j) + FreqPoints(i)
            enddo
            call mat_inv(ham_temp_2,ham_temp_3)
            GF_r(:,:,i) = ham_temp_3(1:nImp,1:nImp)
        enddo
            
        !Are they the same?
        call writeGF(GF_r,'RealSpace_GF')
        call writeGF(GF_k,'KSpace_GF')

        do i = 1,nFreqPoints
            do j = 1,nImp
                do k = 1,nImp
                    if(abs(GF_k(k,j,i)-GF_r(k,j,i)).gt.1.0e-8_dp) then
                        write(6,*) "R and k space greens functions not the same"
                        write(6,*) "i: ",i,FreqPoints(i)
                        write(6,*) "k-space: ",GF_k(k,j,i)
                        write(6,*) "r-space: ",GF_r(k,j,i)
                        stop 'greens functions not the same'
                    endif
                enddo
            enddo
        enddo

        write(6,*) "r and k-space greens functions the same"

        !Now find the derivative of the real-space greens function wrt changing the non-local coupling matrix elements
        !Index 1 & 2 label the matrix
        !Index 3 & 4 label the differential that we are taking
        !Index 5 is the frequency
        allocate(Derivatives(nImp,nImp,nImp,NonLocCoupLength,nFreqPoints))
        Derivatives(:,:,:,:,:) = zzero

        if(tAntiPeriodic) then
            phase = -1.0_dp
        elseif(tPeriodic) then
            phase = 1.0_dp
        endif

        allocate(num(nSites,nSites))

        do w = 1,nFreqPoints

            !B is independent of the specific derivative. Therefore calculate it here.

            ham_temp_2(:,:) = - ham_temp(:,:)
            do j = 1,nSites
                ham_temp_2(j,j) = ham_temp_2(j,j) + FreqPoints(w)
            enddo
            call mat_inv(ham_temp_2,ham_temp_3)
            !ham_temp_3 is now the inverse of the greens function (we could do this with 2x k-space rotations instead?)
            
            do i = 1,NonLocCoupLength
                do j = 1,nImp
                    !Now, calculate the derivative of the inverse of the GF
                    ham_temp_2(:,:) = zzero

                    !Set up matrix with 1's or -1s for values of the coupling matrix element which we are differentiating wrt.
                    do k = 1,nKPnts
                        ind_1 = ((k-1)*nImp) + 1
                        ind_2 = nImp*k

                        !Coupling to the right
                        if((ind_2+i).le.nSites) then
                            ham_temp_2(ind_1+j-1,ind_2+i) = cmplx(-1.0_dp,0.0_dp,dp)
                        else
                            ham_temp_2(ind_1+j-1,ind_2+i-nSites) = cmplx(-phase,0.0_dp,dp)
                        endif

                        !Coupling to the left
                        if((ind_1-i).ge.1) then
                            ham_temp_2(ind_1+j-1,ind_1-i) = cmplx(-1.0_dp,0.0_dp,dp)
                        else
                            ham_temp_2(ind_1+j-1,ind_1-i+nSites) = cmplx(-phase,0.0_dp,dp)
                        endif

                    enddo
!                    write(6,*) "imp: ",j
!                    write(6,*) "CoupParam: ",i
!                    call writematrixcomptoreal(ham_temp_2,'Deriv of B',.true.)

                    !Right, now we can get the full derivative of the inverse of B (i.e. the greens function wrt changing the coefficients).
                    !This currently does not use/require any periodicity. However, the matrix inverse could be improved with periodicity.

                    !The operations below could *certainly* be sped up. Just for testing.
                    call ZGEMM('N','N',nSites,nSites,nSites,zone,ham_temp_2,nSites,ham_temp_3,nSites,zzero,num,nSites)
                    call ZGEMM('N','N',nSites,nSites,nSites,-zone,ham_temp_3,nSites,num,nSites,zzero,ham_temp_2,nSites)

                    !ham_temp_2 is now the derivative of the *full* greens function wrt changing the coupling.
                    !Seems a little profligate - I'm sure the local greens function could be done more efficiently.
                    Derivatives(:,:,j,i,w) = ham_temp_2(1:nImp,1:nImp)
                enddo
            enddo
        enddo
        deallocate(num)

        !Now we have all the derivatives, test if they are right or not!
        allocate(NonLocCoup_tmp(nImp,NonLocCoupLength))
        allocate(num(nImp,nImp))
        open(78,file='DerivsTest',status='unknown')
        do w = 1,nFreqPoints
            do i = 1,NonLocCoupLength
                do j = 1,nImp

                    diff = 0.0001_dp
                    do while(diff.gt.1.0e-12_dp)

                        NonLocCoup_tmp(:,:) = NonLocCoup(:,:)
                        NonLocCoup_tmp(j,i) = NonLocCoup_tmp(j,i) + diff

                        !Calculate greens function for this frequency with new non-local coupling
                        temp(:,:) = ham_temp(:,:)   !Real space hamiltonian
                        call Add_Nonlocal_comp_inplace(temp,NonLocCoup_tmp,NonLocCoupLength,tAdd=.true.)
                        ham_temp_2(:,:) = - temp(:,:)
                        do jj = 1,nSites
                            ham_temp_2(jj,jj) = ham_temp_2(jj,jj) + FreqPoints(w)
                        enddo
                        call mat_inv(ham_temp_2,ham_temp_3)

                        !New greens function is now the 1:nImp block of ham_temp_3
                        !Find the derivative of each element of the greens function
                        num(:,:) = (ham_temp_3(1:nImp,1:nImp) - GF_r(:,:,w)) / diff
                        !num is now the numerical gradient at each point
                        !Find the difference as a single value over all the greens functions
                        err = zero
                        do ii = 1,nImp
                            do jj = 1,nImp

                                err = err + real((num(jj,ii)-Derivatives(jj,ii,j,i,w))*dconjg(num(jj,ii)-Derivatives(jj,ii,j,i,w)),dp)

                            enddo
                        enddo

                        write(78,"(3I7,7G25.14)") w,i,j,diff,real(num(1,1)),aimag(num(1,1)),real(Derivatives(1,1,j,i,w)),aimag(Derivatives(1,1,j,i,w)),err,err/abs(Derivatives(1,1,j,i,w))

                        diff = diff/2.0_dp
                    enddo
                    write(78,"(A)") ""
                enddo
            enddo
        enddo
        close(78)

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

    subroutine WriteGF(gf,gfname)
        implicit none
        complex(dp), intent(in) :: gf(nImp,nImp,nFreqPoints)
        character(len=*), intent(in) :: gfname
        integer :: i,j,k

        open(77,file=gfname,status='unknown')
        do i = 1,nFreqPoints
            write(77,"(F12.6)",advance='no') aimag(FreqPoints(i))
            do j = 1,nImp
                do k = 1,nImp
                    if(k.eq.nImp.and.j.eq.nImp) cycle
                    write(77,"(2F12.6)",advance='no') real(gf(k,j,i),dp),aimag(gf(k,j,i))
                enddo
            enddo
            write(77,"(2F12.6)") real(gf(nImp,nImp,i),dp),aimag(gf(nImp,nImp,i))
        enddo
        close(77)

    end subroutine WriteGF

    
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
  
    SUBROUTINE mat_inv(mat,matinv)
    complex(dp), intent(in) :: mat(:,:)
    complex(dp), dimension(size(mat,1),size(mat,2)), intent(out) :: matinv
    integer, dimension(size(mat,1)) :: ipiv
    integer :: msize,nsize,lwork,info
    complex(dp), allocatable :: cWork(:)

    msize=size(mat,1)
    nsize=size(mat,2)
    if(msize.ne.nsize) stop 'error in z_inv2'
    matinv = mat

    info=0
    call ZGETRF(msize,nsize,matinv,nsize,ipiv,info)
    IF (INFO /= 0) then
        write(6,*) "info: ",info
        STOP 'Error with z_inv2 1'
    endif
    allocate(cWork(1))
    lwork = -1
    call ZGETRI(msize,matinv,msize,ipiv,cwork,lwork,info)
    if(info.ne.0) stop 'error with workspace query in z_inv2'
    lWork = int(real(cWork(1))) + 1
    deallocate(cWork)
    allocate(cWork(lWork))
    call ZGETRI(msize,matinv,msize,ipiv,cwork,lwork,info)
    if(info.ne.0) stop 'error with inversion in z_inv2'
    deallocate(cWork)

  end subroutine mat_inv

end program RealSpaceCoup
