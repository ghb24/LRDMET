module LinearResponse
    use const
    use errors, only: stop_all
    use mat_tools, only: WriteVector,WriteMatrix
    use globals
    implicit none
    contains

    !This is the high level routine to work out how we want to do the linear response
    !TODO:
    !   Code up all of these options!
    !   Really work out difference between non-interacting LR, TDA and RPA, and look at how quality changes in response functions as U increased
    !   Look at difference in quality between TDA-type and RPA-type MCLR methods 
    !   Look at difference in quality between full MCLR and the fully contracted type
    !   Work out self-consistency condition to optimise both the full MCLR and the fully contracted type - response of correlation potential?
    !   Consider using different methods to obtain contraction coefficients in the fully contracted methods - TDA/RPA. Does this improve things?
    !   In the fully contracted case, should we split between impurity and environment, rather than embedded system and core? More semi-internal yuckiness...
    !   Perhaps look at CC2 to get a deeper understanding
    subroutine MR_LinearResponse()
        implicit none
        !Create contracted single excitation space using the non-interacting reference for the contractions
        !The matrix is then created in a CI fashion
        !The difference between this and the one below is that it does not include any operators in the active space, and therefore relies on coupling between
        !the N and N+1 and N-1 active spaces.
        call NonIntContracted_TDA_MCLR_2()

        !Create contracted single excitation space using the non-interacting reference for the contractions
        !The matrix is then created in a CI fashion
        !call NonIntContracted_TDA_MCLR()

        !Create contracted single excitation space using the non-interacting reference for the contractions
        !The matrix is then created in an RPA fashion
!        call NonIntContracted_RPA_MCLR()

        !Full MCLR, creating excitations in a CI fashion, rather than with commutators. Should reduce to TDA in single reference limit
!        call TDA_MCLR()

        !Full MCLR, with excitations and deexcitations. Should reduce to RPA in single reference limit
!        call RPA_MCLR()

    end subroutine MR_LinearResponse

    !Run single reference linear response calculations, based on true HF calculation.
    subroutine SR_LinearResponse()
        implicit none
        
        !Non-interacting linear response
        call NonInteractingLR()
        !Single reference TDA
        call TDA_LR()
        !Single reference RPA
        call RPA_LR()

    end subroutine SR_LinearResponse

    !Create contracted single excitation space using the non-interacting reference for the contractions
    !The matrix is then created in a CI fashion
    !The difference between this and the one below is that it does not include any operators in the active space, and therefore relies on coupling between
    !the N and N+1 and N-1 active spaces.
    subroutine NonIntContracted_TDA_MCLR_2()
        use utils, only: get_free_unit
        use DetToolsData
        use DetBitOps 
        implicit none
        real(dp), allocatable :: LinearSystem(:,:),temp(:,:)
        integer :: i,ierr,iunit,nLinearSystem,a,i_orb,IC,j,k,orb,orbdiffs,spat_orb
        integer :: gtid
        real(dp) :: Omega
        character(len=*), parameter :: t_r='NonIntContracted_TDA_MCLR'
        
        write(6,*) "Calculating non-interacting IC MR-TDA LR version II..."
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(.not.tCompleteDiag) call stop_all(t_r,'To solve LR, must perform complete diag')
        
        iunit = get_free_unit()
        open(unit=iunit,file='IC-TDA_DDResponse',status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"
        
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            write(6,*) "Calculating linear response for frequency: ",Omega

            !First, find the non-interacting solution expressed in the schmidt basis
            call FindSchmidtPert(.false.,Omega)

            !Enumerate excitations for fully coupled space
            call GenDets(Elec,EmbSize,.true.,.true.)
            write(6,*) "Number of determinants in {N,N+1,N-1} FCI space: ",ECoupledSpace

            write(6,*) "N-electron space: "
            do i=1,nFCIDet
                write(6,*) FCIDetList(:,i),FCIBitList(i)
            enddo
            write(6,*) "N-1 electron space: "
            do i=1,nNm1FCIDet
                write(6,*) Nm1FCIDetList(:,i),Nm1BitList(i)
            enddo
            write(6,*) "N+1 electron space: "
            do i=1,nNp1FCIDet
                write(6,*) Np1FCIDetList(:,i),Np1BitList(i)
            enddo

            !Calculate the size of the hamiltonian matrix including the N-1 and N+1 electron active spaces
            !Also include 1 fully contracted core-virtual excitation,
            !plus 2*nImp fully contracted core-active excitation, and 2*nImp fully contracted active-virtual 
            !excitations
            nLinearSystem = ECoupledSpace+1+4*nImp 

            !Allocate memory for hmailtonian in this system:
            write(6,"(A,F14.6,A)") "Allocating memory for the LR hessian: ",real((nLinearSystem**2)*8,dp)/1048576.0_dp," Mb"
            allocate(LinearSystem(nLinearSystem,nLinearSystem),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Error allocating')
            LinearSystem = 0.0_dp

            write(6,"(A)") "Constructing hessian matrix..."

            !***************************************************************************************
            ! Uncontracted functions
            !***************************************************************************************
            
            !First, construct FCI space, in determinant basis
            !This is the first block
            do i=1,nFCIDet
                LinearSystem(i,i) = Spectrum(i) 
            enddo
            !Now transform this block back into the determinant basis
            allocate(temp(nFCIDet,nFCIDet))
            call DGEMM('N','N',nFCIDet,nFCIDet,nFCIDet,1.0_dp,FullHamil,nFCIDet,LinearSystem(1:nFCIDet,1:nFCIDet),  &
                nFCIDet,0.0_dp,temp,nFCIDet)
            call DGEMM('N','T',nFCIDet,nFCIDet,nFCIDet,1.0_dp,temp,nFCIDet,FullHamil,nFCIDet,0.0_dp,    &
                LinearSystem(1:nFCIDet,1:nFCIDet),nFCIDet)
            deallocate(temp)

            !Calculate the FCI space within the N+1 and N-1 space blocks

            !TODO!!!

            !Now, calculate the hamiltonian coupling between the N-1 and N+1 determinant spaces
            !First, coupling between the N-1 and N electron spaces
            do i=1,nFCIDet
                do j=1,nNm1FCIDet

                    !Calculate whether there is a single orbital between the spaces
                    orbdiffs = ieor(FCIBitList(i),Nm1BitList(j))
                    IC = CountBits(orbdiffs)
                    !write(6,"(B32.32,I6)") orbdiffs,IC
                    if(IC.eq.1) then
                        !There is just one orbital different between the two spaces.
                        !Find out which this orbital is!
                        orb = 0
                        do k=0,end_n_int
                            if(btest(orbdiffs,k)) then
                                !An electron in this orbital
                                if(orb.ne.0) call stop_all(t_r,'Should only be one non-zero bit set')
                                orb = k+1
                            endif
                        enddo
                        if(orb.eq.0) call stop_all(t_r,'Did not find orbital')

                        spat_orb = gtid(orb)

                        !Calculate the coupling between these determinants
                        do a = nOcc+nImp+1,nSites 
                            !write(6,*) "Summing into ",i,nFCIDet+j,FockSchmidt(a,spat_orb)
                            LinearSystem(i,nFCIDet+j) = LinearSystem(i,nFCIDet+j) + FockSchmidt(a,spat_orb)
                        enddo
                        LinearSystem(nFCIDet+j,i) = LinearSystem(i,nFCIDet+j)
                    endif

                enddo
            enddo

            !Now for coupling to N+1 space
            do i=1,nFCIDet
                do j=1,nNp1FCIDet

                    !Calculate whether there is a single orbital between the spaces
                    orbdiffs = ieor(FCIBitList(i),Np1BitList(j))
                    IC = CountBits(orbdiffs)
                    if(IC.eq.1) then
                        !There is just one orbital different between the two spaces.
                        !Find out which this orbital is!
                        orb = 0
                        do k=0,end_n_int
                            if(btest(orbdiffs,k)) then
                                !An electron in this orbital
                                if(orb.ne.0) call stop_all(t_r,'Should only be one non-zero bit set')
                                orb = k+1
                            endif
                        enddo
                        if(orb.eq.0) call stop_all(t_r,'Did not find orbital')

                        spat_orb = gtid(orb)

                        !Calculate the coupling between these determinants
                        do i_orb = 1,nOcc-nImp
                            LinearSystem(i,nFCIDet+nNm1FCIDet+j) = LinearSystem(i,nFCIDet+nNm1FCIDet+j) +   &
                                FockSchmidt(i_orb,spat_orb)
                        enddo
                        LinearSystem(nFCIDet+nNm1FCIDet+j,i) = LinearSystem(i,nFCIDet+nNm1FCIDet+j)
                    endif

                enddo
            enddo

            !********************************************************************************************
            !  CORE-VIRTUAL excitations
            !********************************************************************************************

            !Now put in coupling to fully contracted core-virtual excitation
            !Now, calculate the fully IC sum of all core-virtual excitations
            !The diagonal hamiltonian matrix element is given by: G_ai G_bi F_ab - G_ai G_aj F_ij
            !There is no coupling to the active space via the mean-field hamiltonian, since we cant have just a single active space index
            !The contracted excitation is orthogonal to the FCI space uncontracted determinants (since cores are always orthogonal), but 
            !the contracted function itself is not normalized. Find the normalization, and renomalize all its contributions.
            CoreVirtualNorm = 0.0_dp
            do i=1,nOcc-nImp
                do a=nOcc+nImp+1,nSites
                    CoreVirtualNorm = CoreVirtualNorm + SchmidtPert(i,a)*SchmidtPert(a,i)
                enddo
            enddo
            CoreVirtualNorm = CoreVirtualNorm * 2.0_dp  !Spin integration 
            write(6,*) "Fully contracted core-virtual excitations have a normalization of: ",CoreVirtualNorm
            !Strongly contracted function is defined as:
            !   1/sqrt(Norm) |0> \sum_{ai(spat)} G_ai(w) a_a^+ a_i |core>
            do i=1,nOcc-nImp
                do j=1,nOcc-nImp
                    do a=nOcc+nImp+1,nSites
                        LinearSystem(ECoupledSpace+1,ECoupledSpace+1) = LinearSystem(ECoupledSpace+1,ECoupledSpace+1) - &
                            SchmidtPert(a,i)*SchmidtPert(a,j)*FockSchmidt(i,j)
                    enddo
                enddo
            enddo
            do i=1,nOcc-nImp
                do a=nOcc+nImp+1,nSites
                    do b=nOcc+nImp+1,nSites
                        LinearSystem(ECoupledSpace+1,ECoupledSpace+1) = LinearSystem(ECoupledSpace+1,ECoupledSpace+1) +   &
                            SchmidtPert(a,i)*SchmidtPert(b,i)*FockSchmidt(a,b)
                    enddo
                enddo
            enddo
            LinearSystem(ECoupledSpace+1,ECoupledSpace+1) =     &
                LinearSystem(ECoupledSpace+1,ECoupledSpace+1)*2.0_dp/CoreVirtualNorm    !for the other spin type.
            write(6,*) "Diagonal hamiltonian contribution from fully contracted core-virtual function: ",   &
                LinearSystem(ECoupledSpace+1,ECoupledSpace+1)
            !We are not going to add on the active space FCI energy, since we assume that we have offset the hamiltonian by the 
            !zeroth order energy.
            
            !Now, we need to define the coupling to the uncontracted *N-electron* determinant space.
            !We first want to calculate <core|F|\sum_{ia}G_{ai}a_a^+a_i core>
            !The F_ai component must match the G_ai one to be non-zero
            !Therefore we know that this is \sum_ia F_ia G_ai
            CoreCoupling = 0.0_dp
            do i=1,nOcc-nImp
                do a=nOcc+nImp+1,nSites
                    CoreCoupling = CoreCoupling + FockSchmidt(i,a)*SchmidtPert(a,i)
                enddo
            enddo
            CoreCoupling = CoreCoupling * 2.0_dp    !For the other spin type.
            !Now we need to include the contribution from \sum_k C_k <D_j|H|0>
            do i=1,nFCIDet
                do j=1,nFCIDet
                    LinearSystem(i,nFCIDet+1) = LinearSystem(i,ECoupledSpace+1) + LinearSystem(j,i)*FullHamil(j,1)
                enddo
                LinearSystem(i,ECoupledSpace+1) = LinearSystem(i,ECoupledSpace+1) + CoreCoupling    !Add core coupling contribution
                LinearSystem(i,ECoupledSpace+1) = LinearSystem(i,ECoupledSpace+1)/sqrt(CoreVirtualNorm) !Normalise IC single excits
                LinearSystem(ECoupledSpace+1,i) = LinearSystem(i,ECoupledSpace+1)   !Hermiticity
            enddo

            !Now, what about coupling from the core-virtual excitations to the N-1 space...
            !Assume for now that there is no coupling, but potentially there is, via annihilation of the extra electron in the core in the uncontracted function
            !TODO: Work out if this is actually justified
        
            !********************************************************************************************
            !  ACTIVE-VIRTUAL semi-internal excitations
            !********************************************************************************************
            !Precompute the normalization constants for each semi-internal excitation
            allocate(ActiveVirtualNorm(nOcc-nImp+1:nOcc+nImp))
            ActiveVirtualNorm(:) = 0.0_dp
            do alpha = nOcc-nImp+1,nOcc+nImp
                do a=nOcc+nImp+1,nSites
                    ActiveVirtualNorm(alpha) = ActiveVirtualNorm(alpha) +   &
                        SchmidtPert(a,alpha)*SchmidtPert(alpha,a)*2.0_dp
                enddo
            enddo

            !Now calculate the diagonal hamiltonian matrix elements
            ExcitInd = ECoupledSpace + 1
            do alpha=nOcc-nImp+1,nOcc+nImp
                ExcitInd = ExcitInd + 1

                !Now the diagonal hamiltonian matrix element, having subtracted the FCI energy away from it
                do i = 1,nOcc-nImp
                    LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd) + 2.0_dp*FockSchmidt(i,i)
                enddo
                do a = nOcc+nImp+1,nSites
                    do b = nOcc+nImp+1,nSites
                        LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd) + 2.0_dp*SchmidtPert(a,alpha)*    &
                            SchmidtPert(b,alpha)*FockSchmidt(b,a)/ActiveVirtualNorm(alpha)
                    enddo
                enddo

                !Now for coupling to other semi-internal excitations
                ExcitInd2 = ECoupledSpace + 1
                do alphap = nOcc-nImp+1,nOcc+nImp
                    ExcitInd2 = ExcitInd2 + 1
                    if(alphap.eq.alpha) cycle   !Already done diagonals

                    tmp = 0.0_dp
                    do a = nOcc+nImp+1,nSites
                        tmp = tmp + SchmidtPert(a,alpha)*SchmidtPert(a,alphap)
                    enddo
                    LinearSystem(ExcitInd,ExcitInd2) = 2.0_dp*OneElecE*tmp + 2.0_dp*TwoElecE*tmp
                    do a = nOcc+nImp+1,nSites
                        do b = nOcc+nImp+1,nSites
                            LinearSystem(ExcitInd,ExcitInd2) + LinearSystem(ExcitInd,ExcitInd2) + &
                                SchmidtPert(a,alpha)*SchmidtPert(b,alphap)*FockSchmidt(b,a)*2.0_dp
                        enddo
                    enddo
                    LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) /   &
                        sqrt(ActiveVirtualNorm(alpha)*ActiveVirtualNorm(alphap))
                    !Should be able to use hermiticity here, but can check it for a test
                enddo

                !Now for coupling to the N+1 uncontracted determinant block (this is the only block
                !it connects to)
                do dj = 1,nNp1FCIDet

                    !Calculate coupling
                    do di = 1,nFCIDet

                        orbdiffs = ieor(FCIBitList(di),Np1BitList(dj))
                        IC = CountBits(orbdiffs)
                        if(IC.eq.1) then
                            do j=0,end_n_int
                                if(btest(orbdiffs,k)) then
                                    orb = k+1
                                    exit
                                endif
                            enddo
                            spat_orb = gtid(orb)

                            do a = nOcc+nImp+1,nSites
                                LinearSystem(ExcitInd,nFCIDet+nNm1FCIDet+dj) = LinearSystem(ExcitInd,nFCIDet+nNm1FCIDet+dj) + &
                                    SchmidtPert(a,alpha)*FockSchmidt(spat_orb,a)*FullHamil(di,1)
                            enddo

                        endif
                    enddo
                    LinearSystem(ExcitInd,nFCIDet+nNm1FCIDet+dj) = LinearSystem(ExcitInd,nFCIDet+nNm1FCIDet+dj) / sqrt(ActiveVirtualNorm(alpha)
                    LinearSystem(nFCIDet+nNm1FCIDet+dj,ExcitInd) = LinearSystem(ExcitInd,nFCIDet+nNm1FCIDet+dj)
                enddo

            enddo

            !********************************************************************************************
            !  CORE-ACTIVE semi-internal excitations
            !********************************************************************************************

            call stop_all(t_r,'End of test')

        enddo


    end subroutine NonIntContracted_TDA_MCLR_2

    !Contract the basis of single excitations, by summing together all the uncontracted parts with the non-interacting LR coefficients
    !This results in requiring the solution of a system which *only* scales with the size of the active hilbert space, not the lattice
    subroutine NonIntContracted_TDA_MCLR()
        use utils, only: get_free_unit
        use DetToolsData, only: nFCIDet
        implicit none
        integer :: nLinearSystem,ierr,i,j,a,highbound,b,info,lwork,n,iunit
        integer :: i_spat,a_spat,gtid,alpha,beta,gam,dj,ExcitInd
        real(dp) :: CoreCoupling,CoreVirtualNorm,Omega,Res1,Res2,EDiff,ResponseFn 
        real(dp) :: testham,testnorm
        real(dp), allocatable :: LinearSystem(:,:),temp(:,:),Work(:),W(:),Residues(:)
        real(dp), allocatable :: RDM1(:,:),RDM2(:,:),ActiveVirtualNorm(:)
        real(dp), allocatable :: Trans1RDM_bra(:,:),Trans1RDM_ket(:,:)
        character(len=*), parameter :: t_r='NonIntContracted_TDA_MCLR'
        logical, parameter :: tNonIntTest = .false.

        write(6,*) "Calculating non-interacting IC MR-TDA LR system..."
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(.not.tCompleteDiag) call stop_all(t_r,'To solve LR, must perform complete diag')
        
        iunit = get_free_unit()
        open(unit=iunit,file='IC-TDA_DDResponse',status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"

        allocate(Trans1RDM_bra(EmbSize,EmbSize))
        allocate(Trans1RDM_ket(EmbSize,EmbSize))
        
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            write(6,*) "Calculating linear response for frequency: ",Omega

            !First, find the non-interacting solution expressed in the schmidt basis
            call FindSchmidtPert(tNonIntTest,Omega)

            !Calculate the size of the hamiltonian matrix
            !This is simply the normal active space size, plus 1 fully contracted core-virtual excitation,
            !plus 2*nImp fully contracted core-active excitation, and 2*nImp fully contracted active-virtual 
            !excitations
            nLinearSystem = nFCIDet+1+4*nImp !just initially, do without the semi-internal excitations

            !Test that we reduce to the non-interacting limit
            if(tNonIntTest) then
                nLinearSystem = 1
                nFCIDet = 0
                nImp = 0
            endif

            !Allocate memory for hmailtonian in this system:
            write(6,"(A,F14.6,A)") "Allocating memory for the LR hessian: ",real((nLinearSystem**2)*8,dp)/1048576.0_dp," Mb"
            allocate(LinearSystem(nLinearSystem,nLinearSystem),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Error allocating')
            LinearSystem = 0.0_dp

            write(6,"(A)",advance='no') "Constructing hessian matrix..."
            
            !First, construct FCI space, in determinant basis
            !This is the first block
            do i=1,nFCIDet
                LinearSystem(i,i) = Spectrum(i) 
            enddo
            !Now transform this block back into the determinant basis
            allocate(temp(nFCIDet,nFCIDet))
            if(.not.tNonIntTest) then
                call DGEMM('N','N',nFCIDet,nFCIDet,nFCIDet,1.0_dp,FullHamil,nFCIDet,LinearSystem(1:nFCIDet,1:nFCIDet),  &
                    nFCIDet,0.0_dp,temp,nFCIDet)
                call DGEMM('N','T',nFCIDet,nFCIDet,nFCIDet,1.0_dp,temp,nFCIDet,FullHamil,nFCIDet,0.0_dp,    &
                    LinearSystem(1:nFCIDet,1:nFCIDet),nFCIDet)
            endif
            deallocate(temp)

            !Now, calculate the fully IC sum of all core-virtual excitations
            !The diagonal hamiltonian matrix element is given by: G_ai G_bi F_ab - G_ai G_aj F_ij
            !There is no coupling to the active space via the mean-field hamiltonian, since we cant have just a single active space index
            !The contracted excitation in orthogonal to the FCI space uncontracted determinants (since cores are always orthogonal), but 
            !the contracted function itself is not normalized. Find the normalization, and renomalize all its contributions.
            CoreVirtualNorm = 0.0_dp
            do i=1,nOcc-nImp
                do a=nOcc+nImp+1,nSites
                    CoreVirtualNorm = CoreVirtualNorm + SchmidtPert(i,a)*SchmidtPert(a,i)
                enddo
            enddo
            CoreVirtualNorm = CoreVirtualNorm * 2.0_dp  !Spin integration 
            write(6,*) "Fully contracted core-virtual excitations have a normalization of: ",CoreVirtualNorm
            !Strongly contracted function is defined as:
            !   1/sqrt(Norm) |0> \sum_{ai(spat)} G_ai(w) a_a^+ a_i |core>
            !Anything involving this function should come with a 1/sqrt(CoreVirtualNorm) in front of it
            do i=1,nOcc-nImp
                do j=1,nOcc-nImp
                    do a=nOcc+nImp+1,nSites
                        LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1) -     &
                            SchmidtPert(a,i)*SchmidtPert(a,j)*FockSchmidt(i,j)
                    enddo
                enddo
            enddo
            do i=1,nOcc-nImp
                do a=nOcc+nImp+1,nSites
                    do b=nOcc+nImp+1,nSites
                        LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1) +     &
                            SchmidtPert(a,i)*SchmidtPert(b,i)*FockSchmidt(a,b)
                    enddo
                enddo
            enddo
            LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1)*2.0_dp/CoreVirtualNorm    !for the other spin type.
            write(6,*) "Diagonal hamiltonian contribution from fully contracted core-virtual function: ",   &
                LinearSystem(nFCIDet+1,nFCIDet+1)
            if(tNonIntTest) then
                !Check that it gives what we expect
                testham = 0.0_dp
                testnorm = 0.0_dp
                do i=1,nel
                    do a=nel+1,nSites*2
                        if(mod(i,2).ne.mod(a,2)) cycle
                        i_spat = gtid(i)
                        a_spat = gtid(a)
                        testham = testham + SchmidtPert(a_spat,i_spat)*SchmidtPert(i_spat,a_spat)   &
                            *(FullHFEnergies(a_spat)-FullHFEnergies(i_spat))
                        testnorm = testnorm + SchmidtPert(a_spat,i_spat)*SchmidtPert(i_spat,a_spat)
                    enddo
                enddo
                if(abs(testnorm-CoreVirtualNorm).gt.1.0e-7_dp) then
                    write(6,*) "testnorm: ",testnorm
                    write(6,*) "CoreVirtualNorm: ",CoreVirtualNorm
                    !call stop_all(t_r,'non interacting test norm fail')
                endif
                testham = testham/testnorm
                if(abs(testham-LinearSystem(nFCIDet+1,nFCIDet+1)).gt.1.0e-7_dp) then
                    write(6,*) "testham ",testham
                    write(6,*) "Ham matrix el: ",LinearSystem(nFCIDet+1,nFCIDet+1)
                    call stop_all(t_r,'non interacting test hamiltonian fail')
                endif
            endif
            !We are not going to add on the active space energy, since we assume that we have offset the hamiltonian by the 
            !zeroth order energy.

            !Now, we need to define the coupling to the uncontracted determinant space.
            !We first want to calculate <core|F|\sum_{ia}G_{ai}a_a^+a_i core>
            !The F_ai component must match the G_ai one to be non-zero
            !Therefore we know that this is \sum_ia F_ia G_ai
            CoreCoupling = 0.0_dp
            do i=1,nOcc-nImp
                do a=nOcc+nImp+1,nSites
                    CoreCoupling = CoreCoupling + FockSchmidt(i,a)*SchmidtPert(a,i)
                enddo
            enddo
            CoreCoupling = CoreCoupling * 2.0_dp    !For the other spin type.
            !Now we need to include the contribution from \sum_k C_k <D_j|H|0>
            do i=1,nFCIDet
                do j=1,nFCIDet
                    LinearSystem(i,nFCIDet+1) = LinearSystem(i,nFCIDet+1) + LinearSystem(j,i)*FullHamil(j,1)
                enddo
                LinearSystem(i,nFCIDet+1) = LinearSystem(i,nFCIDet+1) + CoreCoupling    !Add core coupling contribution
                LinearSystem(i,nFCIDet+1) = LinearSystem(i,nFCIDet+1)/sqrt(CoreVirtualNorm) !Normalise IC single excits
                LinearSystem(nFCIDet+1,i) = LinearSystem(i,nFCIDet+1)   !Hermiticity
            enddo

            !*****************************************************************************************
            !Now for the active-virtual semi-internal excitations
            !*****************************************************************************************
            !Precompute the normalization constants for each semi-internal excitation
            allocate(ActiveVirtualNorm(nOcc-nImp+1:nOcc+nImp))
            ActiveVirtualNorm(:) = 0.0_dp
            do alpha = nOcc-nImp+1,nOcc+nImp
                do a=nOcc+nImp+1,nSites
                    ActiveVirtualNorm(alpha) = ActiveVirtualNorm(alpha) +   &
                        SchmidtPert(a,alpha)*SchmidtPert(alpha,a)*HL_1RDM(alpha-nOcc+nImp,alpha-nOcc+nImp)
                enddo
            enddo
            !First, creating a particle in the virtual manifold, for each annihilation in the active space
            ExcitInd = nFCIDet + 1
            do alpha=nOcc-nImp+1,nOcc+nImp
                ExcitInd = ExcitInd + 1

                !Now the diagonal hamiltonian matrix element
                do a = nOcc+nImp+1,nSites
                    do b = nOcc+nImp+1,nSites
                        LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd) +     &
                          SchmidtPert(a,alpha)*SchmidtPert(b,alpha)*FockSchmidt(b,a)*HL_1RDM(alpha-nOcc+nImp,alpha-nOcc+nImp)
                    enddo
                enddo
                do i = 1,nOcc-nImp
                    do a = nOcc+nImp+1,nSites
                        LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd) +     &
                          2.0_dp*SchmidtPert(a,alpha)*SchmidtPert(a,alpha)*FockSchmidt(i,i)*    &
                          HL_1RDM(alpha-nOcc+nImp,alpha-nOcc+nImp)
                    enddo
                enddo
                do a = nOcc+nImp+1,nSites
                    do beta = nOcc-nImp+1,nOcc+nImp
                        do gam = nOcc-nImp+1,nOcc+nImp
                            LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd) +     &
                            SchmidtPert(a,alpha)*SchmidtPert(a,alpha)*FockSchmidt(beta,gam)*        &
                            HL_2RDM(alpha-nOcc+nImp,alpha-nOcc+nImp,beta-nOcc+nImp,gam-nOcc+nImp)
                        enddo
                    enddo
                enddo
                !Now add 2electron energy TODO: **** THIS IS NOT CORRECT **** Should be using 3RDMs BOOO...
                LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd) + Two_ElecE 
                !Normalize
                LinearSystem(ExcitInd,ExcitInd) = LinearSystem(ExcitInd,ExcitInd)/ActiveVirtualNorm(alpha)

                !Now for the coupling to the determinant space (dj)
                !Do very inefficiently - can optimize this a lot!
                do dj = 1,nFCIDet
                    !Calculate a transition matrix
                    !Find <Dj|p^+q|0> if tGSKet=T, or the other way around if false.
                    !det is the label to the single determinant function Dj
                    call FindDeterminantTransRDM(dj,Trans1RDM_ket,.true.)
                    call FindDeterminantTransRDM(dj,Trans1RDM_bra,.false.)

                    do a = nOcc+nImp+1,nSites
                        do beta = nOcc-nImp+1,nOcc+nImp
                            LinearSystem(ExcitInd,dj) = LinearSystem(ExcitInd,dj) + SchmidtPert(a,alpha)* &
                                FockSchmidt(beta,a)*Trans1RDM_ket(beta-nOcc+nImp,alpha-nOcc+nImp)
                            LinearSystem(dj,ExcitInd) = LinearSystem(dj,ExcitInd) + SchmidtPert(a,alpha)* &
                                FockSchmidt(beta,a)*Trans1RDM_bra(alpha-nOcc+nImp,beta-nOcc+nImp)
                        enddo
                    enddo
                    LinearSystem(ExcitInd,dj) = LinearSystem(ExcitInd,dj)/sqrt(ActiveVirtualNorm(alpha))
                    LinearSystem(dj,ExcitInd) = LinearSystem(dj,ExcitInd)/sqrt(ActiveVirtualNorm(alpha))
                enddo

                !Now, for the coupling to the strongly contracted excitation
                do i = 1,nOcc-nImp
                    do a = nOcc+nImp+1,nSites
                        do beta = nOcc-nImp+1,nOcc+nImp
                            LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) + &
                                2.0_dp*SchmidtPert(a,i)*SchmidtPert(a,alpha)*FockSchmidt(i,beta)* &
                                HL_1RDM(alpha-nOcc+nImp,beta-nOcc+nImp)
                         enddo
                     enddo
                 enddo
                 LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) /  &
                    sqrt(CoreVirtualNorm*ActiveVirtualNorm(alpha))
                 LinearSystem(nFCIDet+1,ExcitInd) = LinearSystem(ExcitInd,nFCIDet+1)

                 !Now for the coupling to the other semi-internal excitations

             enddo  !Finish looping over active-virtual semi-internal excitations
             deallocate(ActiveVirtualNorm)



                    

            !Is there an overlap matrix we need?
            !call writematrix(SchmidtPert,'SchmidtPert',.true.)
            
            !Only now do we remove the ground state from the FCI space
            do i=1,nFCIDet
                do j=1,nFCIDet
                    LinearSystem(j,i) = LinearSystem(j,i) - Spectrum(1)*FullHamil(i,1)*FullHamil(j,1)
                enddo
                !Finally, subtract the ground state energy from the diagonals, since we want to offset it.
                LinearSystem(i,i) = LinearSystem(i,i) - Spectrum(1)
            enddo

            !Now we have the full hamiltonian. Diagonalize this fully
            allocate(Work(1))
            if(ierr.ne.0) call stop_all(t_r,"alloc err")
            allocate(W(nLinearSystem))
            W(:)=0.0_dp
            lWork=-1
            info=0
            call dsyev('V','U',nLinearSystem,LinearSystem,nLinearSystem,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nLinearSystem,LinearSystem,nLinearSystem,W,Work,lWork,info)
            if (info.ne.0) call stop_all(t_r,"Diag failed")
            deallocate(work)
            
            write(6,*) "First 10 MR-TDA-LR transition frequencies: "
            highbound = min(nLinearSystem,100)
            call writevector(W(1:highbound),'transition frequencies')

            write(6,*) "Calculating residues: "

            !Now find the transition moments
            allocate(RDM1(nSites,nSites))
            allocate(RDM2(nSites,nSites))
            allocate(Residues(nLinearSystem))
            Residues(:) = 0.0_dp

            do n=1,nLinearSystem    !loop over states
                if(.not.tNonIntTest) then
                    !If we have no FCI space at all, we should not even have core contributions
                    !First, find the transition RDM between the ground FCI state, and state n
                    call Calc1RDM(FullHamil(:,1),LinearSystem(1:nFCIDet,n),RDM1)
                    !...and then vice versa
                    call Calc1RDM(LinearSystem(1:nFCIDet,n),FullHamil(:,1),RDM2)
                else
                    RDM1(:,:) = 0.0_dp
                    RDM2(:,:) = 0.0_dp
                endif

                !Now add to the 1RDM the conributions from the fully IC core-active excitations
                do i=1,nOcc-nImp
                    do a=nOcc+nImp+1,nSites

                        RDM1(i,a) = RDM1(i,a) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(i,a)/sqrt(CoreVirtualNorm)
                        !RDM1(a,i) = RDM1(a,i) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(a,i)/sqrt(CoreVirtualNorm)

                        !RDM2(i,a) = RDM2(i,a) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(a,i)/sqrt(CoreVirtualNorm)
                        RDM2(a,i) = RDM2(a,i) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(a,i)/sqrt(CoreVirtualNorm)

                    enddo
                enddo

                !Now, calculate the (pertsite,pertsite) component of this in the AO basis
                Res1 = 0.0_dp
                Res2 = 0.0_dp
                do i = 1,nSites
                    do j = 1,nSites
                        Res1 = Res1 + FullSchmidtBasis(pertsite,i)*RDM1(i,j)*FullSchmidtBasis(pertsite,j)
                        Res2 = Res2 + FullSchmidtBasis(pertsite,i)*RDM2(i,j)*FullSchmidtBasis(pertsite,j)
                    enddo
                enddo

                write(6,*) "Residues for state: ",n,Res1,Res2
                Residues(n) = Res1*Res2*Lambda
            enddo

            ResponseFn = 0.0_dp
            do i=1,nLinearSystem
                if(tNonIntTest) then
                    EDiff = W(i)
                else
                    EDiff = W(i) !-Spectrum(1)
                endif
                ResponseFn = ResponseFn + Residues(i)/(Omega-EDiff)
                ResponseFn = ResponseFn - Residues(i)/(Omega+EDiff)
            enddo
            write(iunit,*) Omega,ResponseFn

            if(tNonIntTest) then
                !Debug comparison info
                if(abs(W(1)-testham).gt.1.0e-7_dp) then
                    write(6,*) "testham: ",testham
                    write(6,*) "Eigenvalue: ",W(1)
                    call stop_all(t_r,'eigenvalue not as expected')
                endif
                testham = 0.0_dp
                do i=1,nel
                    do a=nel+1,nSites*2
                        if(mod(i,2).ne.mod(a,2)) cycle
                        i_spat = gtid(i)
                        a_spat = gtid(a)

                        testham = testham + ((FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2) / &
                            ( (Omega**2/(FullHFEnergies(a_spat)-FullHFEnergies(i_spat))) - 2*Omega + &
                                (FullHFEnergies(a_spat)-FullHFEnergies(i_spat)))

                    enddo
                enddo
                testham = testham / testnorm
                testham = Omega - testham
                !write(6,*) "****",abs(testham-(Omega-EDiff)),abs(testham),abs(testham-(Omega-EDiff))/abs(testham)
                if(abs(testham-(Omega-EDiff)).gt.1.0e-8_dp) then
                    !write(6,*) "test denominator: ",testham
                    !write(6,*) "Calculated Denominator: ",Omega-EDiff
                    !call stop_all(t_r,'non interacting test denominator fail')
                endif

                testham = 0.0_dp
                do i=1,nel
                    do a=nel+1,nSites*2
                        if(mod(i,2).ne.mod(a,2)) cycle
                        i_spat = gtid(i)
                        a_spat = gtid(a)

                        testham = testham + ((FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2) / &
                            (Omega - (FullHFEnergies(a_spat)-FullHFEnergies(i_spat)))

                    enddo
                enddo
                testham = testham / sqrt(testnorm)
                testham = testham*testham
                if(abs(testham-Residues(1)).gt.1.0e-7_dp) then
                    write(6,*) "test residue: ",testham
                    write(6,*) "Calculated Residue: ",Residues(1)
                    call stop_all(t_r,'non interacting test residue fail')
                endif
            endif

            Omega = Omega + Omega_Step
        
            deallocate(Residues,W,LinearSystem,RDM1,RDM2)

        enddo   !Enddo loop over omega
        close(iunit)
        if(tNonIntTest) call stop_all(t_r,'End of NonInt test')

    end subroutine NonIntContracted_TDA_MCLR
            
    !Find <Dj|p^+q|0> if tGSKet=T, or the other way around if false.
    !det is the label to the single determinant function Dj
    subroutine FindDeterminantTransRDM(det,Trans1RDM,tGSKet)
        use DetToolsData, only: FCIDetList,nFCIDet
        implicit none
        integer, intent(in) :: det
        logical, intent(in) :: tGSKet
        real(dp), intent(out) :: Trans1RDM(EmbSize,EmbSize)
        logical :: tSign
        integer :: IC,Ex(2),k,i,iGetExcitLevel,gtid
        character(len=*), parameter :: t_r='FindDeterminantTransRDM'

        Trans1RDM(:,:) = 0.0_dp

        do i=1,nFCIDet
            IC = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,det),Elec)
            if(IC.eq.1) then
                Ex(1) = 1
                if(tGSKet) then
                    call GetExcitation(FCIDetList(:,i),FCIDetList(:,det),Elec,Ex,tSign)
                else
                    call GetExcitation(FCIDetList(:,det),FCIDetList(:,i),Elec,Ex,tSign)
                endif
                if(tSign) then
                    Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) = Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) - FullHamil(i,1)
                else
                    Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) = Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) + FullHamil(i,1)
                endif
            elseif(IC.eq.0) then
                !Same det
                if(i.ne.det) call stop_all(t_r,'Error here')
                do k=1,Elec
                    Trans1RDM(gtid(FCIDetList(k,i)),gtid(FCIDetList(k,i))) =        &
                        Trans1RDM(gtid(FCIDetList(k,i)),gtid(FCIDetList(k,i))) + FullHamil(i,1)
                enddo
            endif
        enddo

    end subroutine FindDeterminantTransRDM

    !Calculate the (transition) 1RDM.
    !Do this very naively to start with.
    subroutine Calc1RDM(Bra,Ket,RDM)
        use DetToolsData, only: nFCIDet,FCIDetList
        implicit none
        real(dp), intent(in) :: Bra(nFCIDet),Ket(nFCIDet)
        real(dp), intent(out) :: RDM(nSites,nSites)
        !real(dp) :: trace
        integer :: nCore,i,j,ex(2),k,i_spat,a_spat,k_spat,IC,igetexcitlevel
        integer :: gtid
        logical :: tSign
        !character(len=*), parameter :: t_r='Calc1RDM'

        nCore = nOcc-nImp

        RDM(:,:) = 0.0_dp
        do i=1,nFCIDet
            do j=1,nFCIDet
                IC = igetexcitlevel(FCIDetList(:,i),FCIDetList(:,j),elec)
                if(IC.eq.1) then
                    !Calculate parity
                    ex(1) = 1
                    call getexcitation(FCIDetList(:,i),FCIDetList(:,j),elec,ex,tSign)
                    i_spat = gtid(ex(1)) + nCore
                    a_spat = gtid(ex(2)) + nCore
                    if(tSign) then
                        RDM(i_spat,a_spat) = RDM(i_spat,a_spat) - Bra(i)*Ket(j)
                    else
                        RDM(i_spat,a_spat) = RDM(i_spat,a_spat) + Bra(i)*Ket(j)
                    endif
                elseif(IC.eq.0) then
                    do k=1,elec
                        k_spat = gtid(FCIDetList(k,i)) + nCore
                        RDM(k_spat,k_spat) = RDM(k_spat,k_spat) + Bra(i)*Ket(i)
                    enddo
                endif
            enddo
        enddo
        do i=1,nCore
            RDM(i,i) = 2.0_dp
        enddo
!
!        trace = 0.0_dp
!        do i=1,nOcc
!            trace = trace + RDM(i,i)
!        enddo
!        if(abs(trace-nel).gt.1.0e-7_dp) then
!            write(6,*) "trace: ",trace
!            write(6,*) "nel: ",nel
!            call stop_all(t_r,'Trace of 1RDM incorrect')
!        endif

    end subroutine Calc1RDM

    !Solve the response equations in the basis of Psi^(0) + all single excits.
    !This is the basis of Psi^(0), the basis of internally contracted single excitations of it into
    !the virtual space, the basis of single excitations of the core into virtual space, and the basis
    !of single excitations of core into active space. This will all be constructed explicitly initially.
    subroutine TDA_MCLR()
        use DetToolsData, only: nFCIDet,FCIDetList
        implicit none
        integer :: nCoreVirt,nCoreActive,nActiveVirt,nLinearSystem,ierr,info
        integer :: i,j,CoreNEl,ind2,ind1,a,x,Ex(2),b
        logical :: tSign_a,tSign_b
        integer, allocatable :: Pivots(:),RefCore(:),Excit1_a(:),Excit1_b(:)
        integer, allocatable :: Excit2_a(:),Excit2_b(:)
        real(dp), allocatable :: LinearSystem(:,:),temp(:,:),Overlap(:,:),Response(:)
        real(dp), allocatable :: ResponseSaved(:), LinearSystemSaved(:,:)
        real(dp) :: Omega,ResponseVal
        character(len=*), parameter :: t_r='SolveDMETResponse'

        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(.not.tCompleteDiag) call stop_all(t_r,'To solve LR, must perform complete diag')

        !Sizes of subblocks
        !FCI space is of size nFCIDet
        nCoreVirt = 2*((nOcc-nImp)**2)
        write(6,*) "Number of core-virtual single excitations: ",nCoreVirt

        !Assume that for the non core-particle conserving excitations, that all have at least some weight for all possible excitations
        nCoreActive = 2*(nOcc-nImp)*(2*nImp)
        write(6,*) "Number of core-active IC-single excitations: ",nCoreActive
        nActiveVirt = 4*nImp*(nSites-nImp-nOcc)
        write(6,*) "Number of active-virtual IC-single excitations: ",nActiveVirt

        nLinearSystem = nFCIDet+nCoreVirt!+nCoreActive+nActiveVirt
        write(6,*) "Size of linear response system: ",nLinearSystem

        !Allocate memory for hmailtonian in this system:
        write(6,"(A,F14.6,A)") "Allocating memory for the LR hessian: ",real((nLinearSystem**2)*8,dp)/1048576.0_dp," Mb"
        allocate(LinearSystem(nLinearSystem,nLinearSystem),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        LinearSystem = 0.0_dp

        write(6,"(A)",advance='no') "Constructing hessian matrix..."

        !First, construct FCI space, in determinant basis
        !This is the first block
        do i=1,nFCIDet
            LinearSystem(i,i) = Spectrum(i) 
        enddo
        !Now transform this block back into the determinant basis
        allocate(temp(nFCIDet,nFCIDet))
        call DGEMM('N','N',nFCIDet,nFCIDet,nFCIDet,1.0_dp,FullHamil,nFCIDet,LinearSystem(1:nFCIDet,1:nFCIDet),  &
            nFCIDet,0.0_dp,temp,nFCIDet)
        call DGEMM('N','T',nFCIDet,nFCIDet,nFCIDet,1.0_dp,temp,nFCIDet,FullHamil,nFCIDet,0.0_dp,    &
            LinearSystem(1:nFCIDet,1:nFCIDet),nFCIDet)
        deallocate(temp)
        !Finally, subtract the ground state energy from the diagonals, since we want to offset it.
        do i=1,nFCIDet
            LinearSystem(i,i) = LinearSystem(i,i) - Spectrum(1)
        enddo

        !TODO: Check here that this is the same as the original determinant basis

        !Allocate memory for a core reference function from which to calculate parities
        CoreNEl = 2*(nOcc-nImp)
        allocate(RefCore(CoreNEl))
        do i=1,CoreNEl
            RefCore(i) = i
        enddo
        allocate(Excit1_a(CoreNEl))
        allocate(Excit1_b(CoreNEl))
        allocate(Excit2_a(CoreNEl))
        allocate(Excit2_b(CoreNEl))

        !Now for the uncontracted core-virtual single excitations
        ind2 = 0    !index for second index excitations
        do i=1,nOcc-nImp    !Run over core orbitals
            do a=nOcc+nImp+1,nSites
                ind2 = ind2 + 1
                if(ind2.gt.nCoreVirt) call stop_all(t_r,'ind2 indexing error')
                !Excitation i -> j

                Excit1_a(:) = RefCore(:)
                Excit1_b(:) = RefCore(:)
                do x=1,CoreNEl
                    if(Excit1_a(x).eq.((2*i)-1)) then
                        Excit1_a(x) = (2*a)-1   !alpha-alpha excitation
                    elseif(Excit1_b(x).eq.(2*i)) then
                        Excit1_b(x) = 2*a       !beta-beta excitation
                    endif
                enddo
                call sort_int(Excit1_a,CoreNEl)
                call sort_int(Excit1_b,CoreNEl)
                !Find parity of single excitations
                Ex(1)=1
                call GetExcitation(RefCore,Excit1_a,CoreNEl,Ex,tSign_a)
                Ex(1)=1
                call GetExcitation(RefCore,Excit1_b,CoreNEl,Ex,tSign_b)

                !First calculate the block connecting it to the FCI space
                !By othogonality, we will pick out the determinant coefficient from the corresponding FCI vector,
                !multiplied by the fock matrix element of the excitation.
                !We also need to parity for the sign of the matrix element
                do j=1,nFCIDet
                    !alpha-alpha is 2*ind2-1
                    if(tSign_a) then
                        LinearSystem(j,nFCIDet+(2*ind2)-1) = -FullHamil(j,1)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2)-1,j) = -FullHamil(j,1)*FockSchmidt(i,a)
                    else
                        LinearSystem(j,nFCIDet+(2*ind2)-1) = FullHamil(j,1)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2)-1,j) = FullHamil(j,1)*FockSchmidt(i,a)
                    endif
                    !beta-beta excitation is 2*ind2
                    if(tSign_b) then
                        LinearSystem(j,nFCIDet+(2*ind2)) = -FullHamil(j,1)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2),j) = -FullHamil(j,1)*FockSchmidt(i,a)
                    else
                        LinearSystem(j,nFCIDet+(2*ind2)) = FullHamil(j,1)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2),j) = FullHamil(j,1)*FockSchmidt(i,a)
                    endif
                enddo

                !Now for diagonal block with other uncontracted excitations
                ind1 = 0    !index for first index excitations
                do j=1,nOcc-nImp
                    do b=nOcc+nImp+1,nSites
                        !Excitation j -> b
                        ind1 = ind1 + 1
                        if(ind1.gt.nCoreVirt) then
                            write(6,*) j,b,ind1
                            call stop_all(t_r,'ind1 indexing error')
                        endif

                        !Find the determinant
                        Excit2_a(:) = RefCore(:)
                        Excit2_b(:) = RefCore(:)
                        do x=1,CoreNEl
                            if(Excit2_a(x).eq.((2*j)-1)) then
                                Excit2_a(x) = (2*b)-1   !alpha-alpha excitation
                            elseif(Excit2_b(x).eq.(2*j)) then
                                Excit2_b(x) = 2*b       !beta-beta excitation
                            endif
                        enddo
                        call sort_int(Excit2_a,CoreNEl)
                        call sort_int(Excit2_b,CoreNEl)

                        if((i.eq.j).and.(a.ne.b)) then
                            ! i = j, a /= b

                            ! alpha/alpha block
                            Ex(1)=1
                            call GetExcitation(Excit1_a,Excit2_a,CoreNEl,Ex,tSign_a)
                            if(tSign_a) then
                                LinearSystem(nFCIDet+(2*ind2)-1,nFCIDet+(2*ind1)-1) = -FockSchmidt(a,b)
                                LinearSystem(nFCIDet+(2*ind1)-1,nFCIDet+(2*ind2)-1) = -FockSchmidt(a,b)
                            else
                                LinearSystem(nFCIDet+(2*ind2)-1,nFCIDet+(2*ind1)-1) = FockSchmidt(a,b)
                                LinearSystem(nFCIDet+(2*ind1)-1,nFCIDet+(2*ind2)-1) = FockSchmidt(a,b)
                            endif

                            !beta/beta excitations
                            Ex(1)=1
                            call GetExcitation(Excit1_b,Excit2_b,CoreNEl,Ex,tSign_b)
                            if(tSign_b) then
                                LinearSystem(nFCIDet+(2*ind2),nFCIDet+(2*ind1)) = -FockSchmidt(a,b)
                                LinearSystem(nFCIDet+(2*ind1),nFCIDet+(2*ind2)) = -FockSchmidt(a,b)
                            else
                                LinearSystem(nFCIDet+(2*ind2),nFCIDet+(2*ind1)) = FockSchmidt(a,b)
                                LinearSystem(nFCIDet+(2*ind1),nFCIDet+(2*ind2)) = FockSchmidt(a,b)
                            endif

                        elseif((a.eq.b).and.(i.ne.j)) then
                            ! a = b, i /= j
                            
                            ! alpha/alpha block
                            Ex(1)=1
                            call GetExcitation(Excit1_a,Excit2_a,CoreNEl,Ex,tSign_a)
                            if(tSign_a) then
                                LinearSystem(nFCIDet+(2*ind2)-1,nFCIDet+(2*ind1)-1) = -FockSchmidt(i,j)
                                LinearSystem(nFCIDet+(2*ind1)-1,nFCIDet+(2*ind2)-1) = -FockSchmidt(i,j)
                            else
                                LinearSystem(nFCIDet+(2*ind2)-1,nFCIDet+(2*ind1)-1) = FockSchmidt(i,j)
                                LinearSystem(nFCIDet+(2*ind1)-1,nFCIDet+(2*ind2)-1) = FockSchmidt(i,j)
                            endif

                            !beta/beta excitations
                            Ex(1)=1
                            call GetExcitation(Excit1_b,Excit2_b,CoreNEl,Ex,tSign_b)
                            if(tSign_b) then
                                LinearSystem(nFCIDet+(2*ind2),nFCIDet+(2*ind1)) = -FockSchmidt(i,j)
                                LinearSystem(nFCIDet+(2*ind1),nFCIDet+(2*ind2)) = -FockSchmidt(i,j)
                            else
                                LinearSystem(nFCIDet+(2*ind2),nFCIDet+(2*ind1)) = FockSchmidt(i,j)
                                LinearSystem(nFCIDet+(2*ind1),nFCIDet+(2*ind2)) = FockSchmidt(i,j)
                            endif

                        elseif((a.eq.b).and.(i.eq.j)) then
                            ! a = b and i = j
                            if(ind1.ne.ind2) call stop_all(t_r,'core-virtual indexing error')

                            !alpha/alpha blocks are diagonal, so it is just faa-fii, since we subtract the core and active diagonal energy terms.
                            LinearSystem(nFCIDet+(2*ind2)-1,nFCIDet+(2*ind1)-1) = FockSchmidt(a,a)-FockSchmidt(i,i)
                            !similarly for the beta beta diagonal term
                            LinearSystem(nFCIDet+(2*ind2),nFCIDet+(2*ind1)) = FockSchmidt(a,a)-FockSchmidt(i,i)

                        endif
                    enddo
                enddo
            enddo
        enddo

!        !Now consider the core-active, and active-virtual IC excitations
!        do j=1,nFCIDet
!            do i=1,nFCIDet
!                ic = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,j),elec)
!                if(ic.eq.1) then
!                    Ex(1)=1
!                    call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
!                    beta = Ex(1)
!                    alpha = Ex(2)
!
!                    !sum over excitations from the core
!                    do core_i = 1,nOcc-nImp
!                        
!                        if(mod(beta,2).eq.1) then
!                            !beta, and hence core_i is an alpha orbital (!)
!                            ExcitIndex = nFCIDet+nCoreVirt+(((alpha+1)/2)*core_i*2)-1
!                            beta_spat = (beta + 1)/2
!                        else
!                            !beta is an beta orbital (!)
!                            ExcitIndex = nFCIDet+nCoreVirt+((alpha/2)*core_i*2)
!                            beta_spat = beta/2
!                        endif
!                        if(ExcitIndex.gt.(nFCIDet+nCoreVirt+nCoreActive)) call stop_all(t_r,'Indexing error')
!
!                        if(tSign) then
!                            LinearSystem(j,ExcitIndex) = LinearSystem(j,ExcitIndex) +   &
!                                FullHamil(i,1)*FockSchmidt(core_i,beta_spat+(nOcc-nImp))
!                            LinearSystem(ExcitIndex,j) = LinearSystem(ExcitIndex,j) +   &
!                                FullHamil(i,1)*FockSchmidt(core_i,beta_spat+(nOcc-nImp))
!                        else
!                            LinearSystem(j,ExcitIndex) = LinearSystem(j,ExcitIndex) -   &
!                                FullHamil(i,1)*FockSchmidt(core_i,beta_spat+(nOcc-nImp))
!                            LinearSystem(ExcitIndex,j) = LinearSystem(ExcitIndex,j) -   &
!                                FullHamil(i,1)*FockSchmidt(core_i,beta_spat+(nOcc-nImp))
!                        endif
!
!                    enddo
!                    
!                    !sum over excitations into the virtual space between the core parts of the wavefunction
!                    do core_a = nOcc+nImp+1,nSites
!                        
!                        if(mod(beta,2).eq.1) then
!                            !beta, and hence core_a is an alpha orbital (!)
!                            ExcitIndex = nFCIDet+nCoreVirt+nCoreActive+(((alpha+1)/2)*(core_a-(nOcc+nImp))*2)-1
!                            beta_spat = (beta + 1)/2
!                        else
!                            !beta is an beta orbital (!)
!                            ExcitIndex = nFCIDet+nCoreVirt+nCoreActive+((alpha/2)*(core_a-(nOcc+nImp))*2)
!                            beta_spat = beta/2
!                        endif
!                        if(ExcitIndex.gt.nLinearSystem) call stop_all(t_r,'Indexing error')
!
!                        if(tSign) then
!                            LinearSystem(j,ExcitIndex) = LinearSystem(j,ExcitIndex) +   &
!                                FullHamil(i,1)*FockSchmidt(core_a,beta_spat+(nOcc-nImp))
!                            LinearSystem(ExcitIndex,j) = LinearSystem(ExcitIndex,j) +   &
!                                FullHamil(i,1)*FockSchmidt(core_a,beta_spat+(nOcc-nImp))
!                        else
!                            LinearSystem(j,ExcitIndex) = LinearSystem(j,ExcitIndex) -   &
!                                FullHamil(i,1)*FockSchmidt(core_a,beta_spat+(nOcc-nImp))
!                            LinearSystem(ExcitIndex,j) = LinearSystem(ExcitIndex,j) -   &
!                                FullHamil(i,1)*FockSchmidt(core_a,beta_spat+(nOcc-nImp))
!                        endif
!
!                    enddo
!
!                elseif(ic.eq.0) then
!
!                    !Need to run through all occupied spin orbtials = alpha
!                    do core_i = 1,nOcc-nImp
!                        do iel = 1,elec
!                            alpha = FCIDetList(iel,i)
!
!                            if(mod(alpha,2).eq.1) then
!                                !alpha is alpha-spin orbital
!                                ExcitIndex = nFCIDet+nCoreVirt+(((alpha+1)/2)*core_i*2)-1
!                                beta_spat = (alpha + 1)/2
!                            else
!                                !alpha is beta-spin orbital
!                                ExcitIndex = nFCIDet+nCoreVirt+((alpha/2)*core_i*2)
!                                beta_spat = alpha/2
!                            endif
!                            if(ExcitIndex.gt.(nFCIDet+nCoreVirt+nCoreActive)) call stop_all(t_r,'Indexing error')
!                            
!                            LinearSystem(j,ExcitIndex) = LinearSystem(j,ExcitIndex) -   &
!                                FullHamil(i,1)*FockSchmidt(core_i,beta_spat+(nOcc-nImp))
!                            LinearSystem(ExcitIndex,j) = LinearSystem(ExcitIndex,j) -   &
!                                FullHamil(i,1)*FockSchmidt(core_i,beta_spat+(nOcc-nImp))
!
!                        enddo
!                    enddo
!                    
!                    !Now do it for active-virtual excitations
!                    do core_a = nOcc+nImp+1,nSites
!                        do iel = 1,elec
!                            alpha = FCIDetList(iel,i)
!
!                            if(mod(alpha,2).eq.1) then
!                                !alpha is alpha-spin orbital
!                                ExcitIndex = nFCIDet+nCoreVirt+nCoreActive+(((alpha+1)/2)*(core_a-(nOcc+nImp))*2)-1
!                                beta_spat = (alpha + 1)/2
!                            else
!                                !alpha is beta-spin orbital
!                                ExcitIndex = nFCIDet+nCoreVirt+nCoreActive+((alpha/2)*(core_a-(nOcc+nImp))*2)
!                                beta_spat = alpha/2
!                            endif
!                            if(ExcitIndex.gt.nLinearSystem) call stop_all(t_r,'Indexing error')
!                            
!                            LinearSystem(j,ExcitIndex) = LinearSystem(j,ExcitIndex) -   &
!                                FullHamil(i,1)*FockSchmidt(core_a,beta_spat+(nOcc-nImp))
!                            LinearSystem(ExcitIndex,j) = LinearSystem(ExcitIndex,j) -   &
!                                FullHamil(i,1)*FockSchmidt(core_a,beta_spat+(nOcc-nImp))
!
!                        enddo
!                    enddo
!
!                endif
!            enddo
!        enddo

        !The FCI and core-virtual excitation blocks are orthogonal, therefore let the overlap initially be unit in these blocks
        allocate(Overlap(nLinearSystem,nLinearSystem))
        Overlap = 0.0_dp
        do i=1,nFCIDet+nCoreVirt
            Overlap(i,i) = 1.0_dp
        enddo

        !The V|Psi^0> only has weight over the active space, since it spans the same FCI space, and doesn't excite into the rest of the space
        allocate(Response(nLinearSystem))
        Response = 0.0_dp
        do i=1,nFCIDet
            do j=1,elec
                if(FCIDetList(j,i).eq.(pertsite*2)-1) then
                    !alpha spin of the perturbation
                    Response(i) = Response(i) + Lambda*FullHamil(i,1)
                elseif(FCIDetList(j,i).eq.(pertsite*2)) then
                    !beta spin of the perturbation
                    Response(i) = Response(i) + Lambda*FullHamil(i,1)
                endif
            enddo
        enddo
        write(6,"(A)") "done."
        
        !Save linear system for use with multiple omegas
        allocate(LinearSystemSaved(nLinearSystem,nLinearSystem))
        LinearSystemSaved(:,:) = LinearSystem(:,:)
        allocate(ResponseSaved(nLinearSystem))
        ResponseSaved(:) = Response(:)
        allocate(Pivots(nLinearSystem))

        !Now solve the equations...
        write(6,"(A)") "Solving linear response equations..."
        
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            write(6,*) "Omega = ",Omega

            LinearSystem(:,:) = LinearSystemSaved(:,:) - Omega*Overlap(:,:)
            Response(:) = ResponseSaved(:)

            !Now, solve the linear equation Ax = b, where b is Response, A is the hessian and x will be the 1st order wavefunction
            call DGESV(nLinearSystem,1,LinearSystem,nLinearSystem,Pivots,Response,nLinearSystem,info)
            if(info.ne.0) call stop_all(t_r,'Error with solving linear system')
                        
            !Response is now |Psi^1>
            !Now, apply the perturbation again, and project back onto the zeroth order wavefunction.
            !Since the zeroth order wavefunction only spans the active space, we only need to consider how V|Psi^1> changes the active space
            ResponseVal = 0.0_dp
            do i=1,nFCIDet
                do j=1,elec
                    if(FCIDetList(j,i).eq.(pertsite*2)-1) then
                        ResponseVal = ResponseVal + FullHamil(i,1)*Response(i)*Lambda
                    elseif(FCIDetList(j,i).eq.(pertsite*2)) then
                        ResponseVal = ResponseVal + FullHamil(i,1)*Response(i)*Lambda
                    endif
                enddo
            enddo
                        
            write(6,*) "Response function: ",Omega,ResponseVal

            Omega = Omega + Omega_Step

        enddo


    end subroutine TDA_MCLR

    !Find the non-interacting perturbation, and project this operator into the schmidt basis of phi^0 + its virtual space
    subroutine FindSchmidtPert(tNonIntTest,Omega)
        implicit none
        logical, intent(in) :: tNonIntTest  !Test to return just the perturbation in the normal HF basis
        real(dp), allocatable :: HFPertBasis(:,:),temp(:,:)
        real(dp) :: EDiff
        real(dp), intent(in) :: Omega
        integer :: a,i
        character(len=*), parameter :: t_r='FindSchmidtBasis'

        if(tNonIntTest) then
            !Screw everything up by ensuring that HFOrbs = FullHFOrbs (likewise for energy eigenvalues)
            HFOrbs(:,:) = FullHFOrbs(:,:)
            HFEnergies(:) = FullHFEnergies(:)
        endif

        allocate(HFPertBasis(nSites,nSites))
        HFPertBasis(:,:) = 0.0_dp

        !Assume perturbation is local to the first impurity site (pertsite = 1).
        if(pertsite.ne.1) call stop_all(t_r,'Perturbation is not local to the impurity site')
        !Assuming that MF DD operator is V_ia/e_a-e_i-w + V_ia/e_a-e_i+w
        do i=1,nOcc
            do a=nOcc+1,nSites
                EDiff = HFEnergies(a)-HFEnergies(i)
                HFPertBasis(i,a) = HFOrbs(pertsite,i)*HFOrbs(pertsite,a)*Lambda*(1.0_dp/(Omega-EDiff))
                HFPertBasis(a,i) = HFPertBasis(i,a)
            enddo
        enddo

        !call writematrix(HFPertBasis,'Perturbation in HF basis',.true.)
        !write(6,*) "Transforming non-interacting response operator into full schmidt basis..."

        allocate(temp(nSites,nSites))
        
        !call DGEMM('N','T',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,HFtoSchmidtTransform,nSites,0.0_dp,temp,nSites)
        !call writematrix(temp,'Test of unitarity of HF to Schmidt Transform',.true.)
        !do i=1,nSites
        !    do j=1,nSites
        !        if((i.ne.j).and.(abs(temp(i,j)).gt.1.0e-7_dp)) then
        !            call stop_all(t_r,'Transformation matrix not unitary')
        !        elseif((i.eq.j).and.(abs(temp(i,j)-1.0_dp).gt.1.0e-7)) then
        !            call stop_all(t_r,'Transformation matrix not unitary')
        !        endif
        !    enddo
        !enddo
        
        if(allocated(SchmidtPert)) deallocate(SchmidtPert)
        allocate(SchmidtPert(nSites,nSites))
        if(tNonIntTest) then
            SchmidtPert(:,:) = HFPertBasis(:,:)
            FockSchmidt(:,:) = 0.0_dp
            do i=1,nSites
                FockSchmidt(i,i) = HFEnergies(i)
            enddo
            FullSchmidtBasis(:,:) = HFOrbs(:,:) 
        else
            call DGEMM('T','N',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,HFPertBasis,nSites,0.0_dp,temp,nSites)
            call DGEMM('N','N',nSites,nSites,nSites,1.0_dp,temp,nSites,HFtoSchmidtTransform,nSites,0.0_dp,SchmidtPert,nSites)
        endif
        deallocate(temp,HFPertBasis)

        !SchmidtPert is now the perturbation in the schmidt basis
        !call writematrix(SchmidtPert,'Perturbation in schmidt basis',.true.)

    end subroutine FindSchmidtPert

    subroutine NonInteractingLR()
        use utils, only: get_free_unit
        implicit none
        integer :: ov_space,virt_start,i,a,a_spat,i_spat,ai_ind,gtid,iunit
        integer :: highbound
        real(dp) :: Omega,EDiff,ResponseFn
        real(dp), allocatable :: transitions(:,:)   !(ov_space,2)   !1 = transition frequencies, 2 = moments
        !character(len=*), parameter :: t_r='NonInteractingLR'

        write(6,*) "Calculating the non-interacting linear response function"

        !First, just enumerate transitions
        ov_space =2*nOcc*(nSites-nOcc)
        virt_start = (2*nOcc)+1
        allocate(transitions(ov_space,2))
        transitions(:,:) = 0.0_dp
        do i=1,nel
            do a=virt_start,2*nSites
                if(mod(i,2).ne.mod(a,2)) cycle      !Only want same spin excitations 
                ai_ind = ov_space_spinind(a,i)
                i_spat = gtid(i)
                a_spat = gtid(a)

                transitions(ai_ind,1) = FullHFEnergies(a_spat)-FullHFEnergies(i_spat)
                !Now calculate the moment
                transitions(ai_ind,2) = (FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2
!                write(6,*) "i, a: ",i,a
!                write(6,*) "moment: ",(FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2
                !write(6,*) Transitions(ai_ind,1),Transitions(ai_ind,2)
            enddo
        enddo
        call sort_real2(transitions,ov_space,2)

        iunit = get_free_unit()
        open(unit=iunit,file='NonInt_Transitions',status='unknown')
        write(iunit,"(A)") "#Excitation     Transition_Frequency       Transition_Moment"
        do i=1,ov_space
            write(iunit,"(I8,2G22.12)") i,transitions(i,1),transitions(i,2)
        enddo
        close(iunit)
        write(6,*) "First 10 non-interacting transition frequencies: "
        highbound = min(ov_space,10)
        call writevector(transitions(1:highbound,1),'transition frequencies')
        deallocate(transitions)

        write(6,*) "Writing non-interacting linear response function to disk..."

        open(unit=iunit,file='NonInt_DDResponse',status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ResponseFn = 0.0_dp
            do i=1,nel
                do a=virt_start,2*nSites
                    if(mod(i,2).ne.mod(a,2)) cycle      !Only want same spin excitations 
                    i_spat = gtid(i)
                    a_spat = gtid(a)

                    EDiff = FullHFEnergies(a_spat)-FullHFEnergies(i_spat)
                    ResponseFn = ResponseFn + ((FullHFOrbs(pertsite,a_spat)*FullHFOrbs(pertsite,i_spat))**2)/(Omega-EDiff)
                    ResponseFn = ResponseFn - ((FullHFOrbs(pertsite,a_spat)*FullHFOrbs(pertsite,i_spat))**2)/(Omega+EDiff)
                enddo
            enddo
            ResponseFn = ResponseFn*Lambda
            write(iunit,*) Omega,ResponseFn

            Omega = Omega + Omega_Step

        enddo
        close(iunit)

    end subroutine NonInteractingLR

    subroutine TDA_LR()
        use utils, only: get_free_unit
        use DetToolsData, only: tmat,umat
        implicit none
        integer :: ov_space,virt_start,ierr,i,j,n,m,nj_ind,mi_ind,ex(2,2),gtid
        integer :: m_spat,i_spat,lwork,info,k,umatind,l,orbpairs,umatsize,ai_ind,a
        integer :: state,iunit,a_spat,highbound
        logical :: tSign
        integer, allocatable :: detHF(:),detR(:),detL(:)
        real(dp) :: HEl1,GetHFAntisymInt_spinorb,GetHFInt_spinorb,Omega,ResponseFn
        real(dp), allocatable :: A_mat(:,:),W(:),Work(:),temp(:,:),Residues(:)
        real(dp), allocatable :: DM(:,:),DM_conj(:,:),DM_AO(:,:),DM_AO_conj(:,:)
        character(len=*), parameter :: t_r='TDA_LR'

        write(6,*) "Calculating the linear response function via the Tamm-Dancoff approximation"

        ov_space =2*nOcc*(nSites-nOcc)
        virt_start = (2*nOcc)+1

        if(.false.) then
            !Temporarily create & store TMAT & UMAT
            if(allocated(tmat)) deallocate(tmat)
            if(allocated(umat)) deallocate(umat)
            allocate(tmat(nSites,nSites))
            allocate(temp(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,tmat,nSites)
            deallocate(temp)
            OrbPairs = (nSites*(nSites+1))/2
            umatsize = (OrbPairs*(OrbPairs+1))/2 
            allocate(umat(umatsize))
            do i=1,nSites
                do j=1,nSites
                    do k=1,nSites
                        do l=1,nSites
                            ex(1,1) = i*2
                            ex(1,2) = j*2
                            ex(2,1) = k*2
                            ex(2,2) = l*2
                            umat(umatind(i,j,k,l)) = GetHFInt_spinorb(ex,FullHFOrbs)
                        enddo
                    enddo
                enddo
            enddo

            !calculate matrix brute force to check it is right
            allocate(A_mat(ov_space,ov_space),stat=ierr)    !This is the singles hamiltonian
            if(ierr.ne.0) call stop_all(t_r,'alloc error')
            A_mat(:,:) = 0.0_dp

            allocate(detL(nel))
            allocate(detR(nel))
            allocate(detHF(nel))
            do k=1,nel
                detHF(k) = k
            enddo

            do j=1,nel
                do n=virt_start,2*nSites
                    if(mod(j,2).ne.mod(n,2)) cycle      !Only want same spin excitations for j and n
                    detL(:) = detHF(:)
                    do k=1,nel
                        if(detL(k).eq.j) then
                            detL(k) = n
                            exit
                        endif
                    enddo
                    call sort_int(detL,nel)
                    nj_ind = ov_space_spinind(n,j)

                    do i=1,nel
                        do m=virt_start,2*nSites
                            if(mod(i,2).ne.mod(m,2)) cycle      !Only want same spin excitations for i and m
                            detR(:) = detHF(:)
                            do k=1,nel
                                if(detR(k).eq.i) then
                                    detR(k) = m
                                    exit
                                endif
                            enddo
                            call sort_int(detR,nel)
                            mi_ind = ov_space_spinind(m,i)

                            call GetHElement(detL,detR,nel,HEl1)
                            A_mat(nj_ind,mi_ind) = HEl1
                        enddo
                    enddo
                enddo
            enddo

            !Remove the HF energy from the diagonals
            do i=1,ov_space
                A_mat(i,i) = A_mat(i,i) - HFEnergy
            enddo

            !Check that A is hermition 
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7) then
                        call stop_all(t_r,'A not hermitian')
                    endif
                enddo
            enddo

            !Diagonalize A
            allocate(Work(1))
            if(ierr.ne.0) call stop_all(t_r,"alloc err")
            allocate(W(ov_space))
            W(:)=0.0_dp
            lWork=-1
            info=0
            call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lwork,info)
            if (info.ne.0) call stop_all(t_r,"Diag failed")
            deallocate(work)

            write(6,*) "The first 10 transition frequencies are: "
            highbound = min(ov_space,10)
            call writevector(W(1:highbound),'Transition frequencies')

            deallocate(A_Mat,W,umat,tmat,detL,detR,detHF)

        endif
        
        allocate(A_mat(ov_space,ov_space),stat=ierr)    !This is the singles hamiltonian
        if(ierr.ne.0) call stop_all(t_r,'alloc error')
        A_mat(:,:) = 0.0_dp

        !First, construct A and B, which correspond to:
        !  <HF| [a*_i a_m [H, a*_n a_j]] |HF> = A
        do j=1,nel
            ex(1,2) = j     !second index in integral
            do n=virt_start,2*nSites
                if(mod(j,2).ne.mod(n,2)) cycle      !Only want same spin excitations for j and n
                nj_ind = ov_space_spinind(n,j)
                ex(2,2) = n !4th index in integral
                do i=1,nel
                    ex(2,1) = i !3rd index in integral
                    do m=virt_start,2*nSites
                        if(mod(i,2).ne.mod(m,2)) cycle      !Only want same spin excitations for i and m
                        mi_ind = ov_space_spinind(m,i)
                        ex(1,1) = m !First index in integral

                        !Calculate the antisymmetrized integral, < m j || i n > for A_mat and < m n || i j > for B_mat
                        HEl1 = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                        A_mat(mi_ind,nj_ind) = HEl1
                    enddo
                enddo
            enddo
        enddo

        !Now add the diagonal part to A
        do i=1,nel
            do m=virt_start,2*nSites
                if(mod(i,2).ne.mod(m,2)) cycle  !Only want same spin excitations
                mi_ind = ov_space_spinind(m,i)

                m_spat = gtid(m)
                i_spat = gtid(i)

                A_mat(mi_ind,mi_ind) = A_mat(mi_ind,mi_ind) + (FullHFEnergies(m_spat)-FullHFEnergies(i_spat))
            enddo
        enddo

        !Check that A is hermition 
        do i=1,ov_space
            do j=1,ov_space
                if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7) then
                    call stop_all(t_r,'A not hermitian')
                endif
            enddo
        enddo

        !Diagonalize A
        allocate(Work(1))
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        allocate(W(ov_space))
        W(:)=0.0_dp
        lWork=-1
        info=0
        call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lwork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)

        write(6,*) "The first 10 transition frequencies are: "
        highbound = min(ov_space,10)
        call writevector(W(1:highbound),'Transition frequencies')
        write(6,*) "Calculating TDA Transition moments..."
        call flush(6)

        !call writematrix(A_mat,'evecs',.true.)

        !We now have our eigenvalues and eigenvectors. Calculate the response functions
        !First, calculate the residues
        allocate(Residues(ov_space))
        Residues(:) = 0.0_dp
        allocate(DM(nSites,nSites))
        allocate(DM_conj(nSites,nSites))
        allocate(DM_AO(nSites,nSites))
        allocate(DM_AO_conj(nSites,nSites))
        allocate(detL(nel))
        allocate(detHF(nel))
        do k=1,nel
            detHF(k) = k
        enddo
        allocate(temp(nSites,nSites))

        do i=1,nel
            do a=virt_start,2*nSites
                if(mod(i,2).ne.mod(a,2)) cycle  !Only want same spin excitations
                ai_ind = ov_space_spinind(a,i)
                detL(:) = detHF(:)
                do k=1,nel
                    if(detL(k).eq.i) then
                        detL(k) = a
                    endif
                enddo
                call sort_int(detL(:),nel)
                i_spat = gtid(i)
                a_spat = gtid(a)

                !Calculate permutation
                ex(1,1) = 1
                call GetExcitation(detL,detHF,nel,ex,tSign)

                DM(:,:) = 0.0_dp
                DM_conj(:,:) = 0.0_dp
                if(tSign) then
                    DM(a_spat,i_spat) = -1.0_dp
                    DM_conj(i_spat,a_spat) = -1.0_dp
                else
                    DM(a_spat,i_spat) = 1.0_dp
                    DM_conj(i_spat,a_spat) = 1.0_dp
                endif

                !Transfer to AO basis
                call dgemm('n','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,DM,nSites,0.0_dp,temp,nSites)
                call dgemm('n','t',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,DM_AO,nSites)
                !Do it also for the conjugate density matrix
                call dgemm('n','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,DM_conj,nSites,0.0_dp,temp,nSites)
                call dgemm('n','t',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,DM_AO_conj,nSites)

                !Extract pertsite,pertsite component
                !Now run over all states
                do state=1,ov_space
                    Residues(state) = Residues(state) + DM_AO(pertsite,pertsite)*DM_AO_conj(pertsite,pertsite)* &
                        A_mat(ai_ind,state)*A_mat(ai_ind,state)
!                    if(abs(A_mat(ai_ind,state)).gt.0.5_dp) then
!                        write(6,*) "i, a: ",i,a
!                        write(6,*) "moment: ",(FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2
!                        write(6,*) "Actually: ",DM_AO(pertsite,pertsite)*DM_AO_conj(pertsite,pertsite)
!                        write(6,*) "coefficient: ",A_mat(ai_ind,state)
!                        write(6,*) "ai_ind: ",ai_ind
!                        write(6,*) DM_AO(pertsite,pertsite),DM_AO_conj(pertsite,pertsite)
!                        write(6,*) FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat)
!                    endif
                enddo
            enddo
        enddo
        deallocate(detHF,detL,temp,DM,DM_conj,DM_AO,DM_AO_conj)
        Residues(:) = Residues(:)*Lambda
        
        iunit = get_free_unit()
        open(unit=iunit,file='TDA_Transitions',status='unknown')
        write(iunit,"(A)") "#Excitation     Transition_Frequency       Transition_Moment"
        do i=1,ov_space
            write(iunit,"(I8,2G22.12)") i,W(i),Residues(i)
        enddo
        close(iunit)

        write(6,*) "Writing TDA linear response function to disk..."

        open(unit=iunit,file='TDA_DDResponse',status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ResponseFn = 0.0_dp
            do i=1,ov_space
                ResponseFn = ResponseFn + ((Residues(i)/(Omega-W(i))) - (Residues(i)/(Omega+W(i))))
            enddo
            write(iunit,*) Omega,ResponseFn

            Omega = Omega + Omega_Step

        enddo
        close(iunit)

        deallocate(A_Mat,W,Residues)

    end subroutine TDA_LR
    
    !Set up the RPA equations and solve them.
    !Finally, create the density-density linear response function from the resulting excitations/deexcitations
    !This is done in the spin-orbital space
    subroutine RPA_LR()
        use utils, only: get_free_unit
        use matrixops, only: d_inv
        implicit none
        integer :: ov_space,virt_start,ierr,j,ex(2,2),ex2(2,2),n,i,m,nj_ind,mi_ind,info,lwork
        integer :: m_spat,i_spat,StabilitySize,mu,gtid,j_spat,ai_ind,iunit,a,excit,highbound
        real(dp) :: HEl1,HEl2,X_norm,Y_norm,norm,Energy_stab,DMEl1,DMEl2,Omega,ResponseFn
        real(dp) :: GetHFAntisymInt_spinorb
        real(dp), allocatable :: A_mat(:,:),B_mat(:,:),Stability(:,:),StabilityCopy(:,:),W(:),Work(:)
        real(dp), allocatable :: S_half(:,:),temp(:,:),temp2(:,:),W2(:),X_stab(:,:),Y_stab(:,:)
        real(dp), allocatable :: trans_moment(:),AOMO_Spin(:,:),DM(:,:)
        character(len=*), parameter :: t_r='RPA_LR'

        ov_space = 2*nOcc*(nSites-nOcc)
        virt_start = (2*nOcc)+1
        allocate(A_mat(ov_space,ov_space),stat=ierr)
        allocate(B_mat(ov_space,ov_space),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'alloc error')
        A_mat(:,:) = 0.0_dp
        B_mat(:,:) = 0.0_dp

        !First, construct A and B, which correspond to:
        !  <HF| [a*_i a_m [H, a*_n a_j]] |HF> = A
        ! -<HF| [a*_i a_m [H, a*_j a_n]] |HF> = B
        do j=1,nel
            ex(1,2) = j     !second index in integral
            ex2(2,2) = j
            do n=virt_start,2*nSites
                if(mod(j,2).ne.mod(n,2)) cycle      !Only want same spin excitations for j and n
                nj_ind = ov_space_spinind(n,j)
                ex(2,2) = n !4th index in integral
                ex2(1,2) = n
                do i=1,nel
                    ex(2,1) = i !3rd index in integral
                    ex2(2,1) = i
                    do m=virt_start,2*nSites
                        if(mod(i,2).ne.mod(m,2)) cycle      !Only want same spin excitations for i and m
                        mi_ind = ov_space_spinind(m,i)
                        ex(1,1) = m !First index in integral
                        ex2(1,1) = m

                        !Calculate the antisymmetrized integral, < m j || i n > for A_mat and < m n || i j > for B_mat
                        HEl1 = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                        HEl2 = GetHFAntisymInt_spinorb(ex2,FullHFOrbs)
                        A_mat(mi_ind,nj_ind) = HEl1
                        B_mat(mi_ind,nj_ind) = HEl2
                    enddo
                enddo
            enddo
        enddo

        !Now add the diagonal part to A
        do i=1,nel
            do m=virt_start,2*nSites
                if(mod(i,2).ne.mod(m,2)) cycle  !Only want same spin excitations
                mi_ind = ov_space_spinind(m,i)

                m_spat = gtid(m)
                i_spat = gtid(i)

                A_mat(mi_ind,mi_ind) = A_mat(mi_ind,mi_ind) + (FullHFEnergies(m_spat)-FullHFEnergies(i_spat))
            enddo
        enddo

        !Check that A is hermition and B is symmetric
        do i=1,ov_space
            do j=1,ov_space
                if(abs(B_mat(i,j)-B_mat(j,i)).gt.1.0e-7_dp) then
                    call stop_all(t_r,'B not symmetric')
                endif
                if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7) then
                    call stop_all(t_r,'A not hermitian')
                endif
            enddo
        enddo

        !Calculate here via direct diagonalization of the stability matrix
        write(6,*) "Calculating RPA from stability matrix"
        call flush(6)

        !Stability = ( A  B  )
        !            ( B* A* )
        !Assume all integrals real to start with
        StabilitySize=2*ov_space
        allocate(Stability(StabilitySize,StabilitySize),stat=ierr)
        Stability(:,:)=0.0_dp
        Stability(1:ov_space,1:ov_space)=A_mat(1:ov_space,1:ov_space)
        Stability(ov_space+1:StabilitySize,1:ov_space)=B_mat(1:ov_space,1:ov_space)
        Stability(1:ov_space,ov_space+1:StabilitySize)=B_mat(1:ov_space,1:ov_space)
        Stability(ov_space+1:StabilitySize,ov_space+1:StabilitySize)=A_mat(1:ov_space,1:ov_space)

        !Now diagonalize
        !Find optimal space
        allocate(StabilityCopy(StabilitySize,StabilitySize))
        StabilityCopy(:,:)=Stability(:,:)
        allocate(W(StabilitySize),stat=ierr)    !Eigenvalues of stability matrix
        allocate(Work(1))
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        W(:)=0.0_dp
        lWork=-1
        info=0
        call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lwork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)

        do i=1,StabilitySize
            if(W(i).lt.0.0_dp) then
                write(6,*) i,W(i)
                call stop_all(t_r,"HF solution not stable. Not local minimum. Recompute HF.")
            endif
        enddo
        write(6,"(A)") "Stability matrix positive definite. HF solution is minimum. RPA stable"

        !Now compute S^(1/2), and transform into original basis
        allocate(S_half(StabilitySize,StabilitySize))
        S_half(:,:)=0.0_dp
        do i=1,StabilitySize
            S_half(i,i)=sqrt(W(i))
        enddo
        allocate(temp(StabilitySize,StabilitySize))
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,    &
            S_half,StabilitySize,0.0_dp,temp,StabilitySize)
        call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp,StabilitySize,Stability,   &
            StabilitySize,0.0_dp,S_half,StabilitySize)
        !S_half is now S^1/2 in the original basis

        !Check this by squaring it.
        call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,S_half,    &
            StabilitySize,0.0_dp,temp,StabilitySize)
        do i=1,StabilitySize
            do j=1,StabilitySize
                if(abs(StabilityCopy(i,j)-temp(i,j)).gt.1.0e-7) then
                    call stop_all(t_r,'S^1/2 not calculated correctly in original basis')
                endif
            enddo
        enddo

        temp(:,:)=0.0_dp
        do i=1,ov_space
            temp(i,i)=1.0_dp
        enddo
        do i=ov_space+1,StabilitySize
            temp(i,i)=-1.0_dp
        enddo
        allocate(temp2(StabilitySize,StabilitySize))
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,temp,  &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,S_half, &
            StabilitySize,0.0_dp,temp,StabilitySize)
        !Now diagonalize temp = S^(1/2) (1 0 \\ 0 -1 ) S^(1/2)

        lWork=-1
        allocate(W2(StabilitySize))
        allocate(Work(1))
        W2(:)=0.0_dp
        call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lwork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)
!            call writevector(W2,'Excitation energies')
        ! temp now holds the eigenvectors X~ Y~
        ! W2 runs over StabilitySize eigenvalues (ov_space*2). Therefor we expect redundant pairs of +-W2, corresponding
        ! to pairs of eigenvectors (X^v Y^v) and (X^v* Y^v*) (Same in real spaces).
        do i=1,ov_space
            !This they are listed in order of increasing eigenvalue, we should be able to easily check that they pair up
            if(abs(W2(i)+W2(StabilitySize-i+1)).gt.1.0e-7_dp) then
                write(6,*) i,StabilitySize-i+1, W2(i), W2(StabilitySize-i+1), abs(W2(i)-W2(StabilitySize-i+1))
                call stop_all(t_r,"Excitation energy eigenvalues do not pair")
            endif
        enddo

        !We actually have everything we need for the energy already now. However, calculate X and Y too.
        !Now construct (X Y) = S^(-1/2) (X~ Y~)
        !First get S^(-1/2) in the original basis
        S_half(:,:)=0.0_dp
        do i=1,StabilitySize
            S_half(i,i)=-sqrt(W(i))
        enddo
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,S_half, &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,Stability,  &
            StabilitySize,0.0_dp,S_half,StabilitySize)
        !S_half is now S^(-1/2) in the original basis

        !Now multiply S^(-1/2) (X~ y~)
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,temp,  &
            StabilitySize,0.0_dp,temp2,StabilitySize)

        !Check that eigenvectors are also paired.
        !Rotations among degenerate sets will screw this up though
!            do i=1,ov_space
!                write(6,*) "Eigenvectors: ",i,StabilitySize-i+1,W2(i),W2(StabilitySize-i+1)
!                do j=1,StabilitySize
!                    write(6,*) j,temp2(j,i),temp2(j,StabilitySize-i+1)
!                enddo
!            enddo
!            call writematrix(temp2,'X Y // Y X',.true.)
        !temp2 should now be a matrix of (Y X)
!                                            (X Y)
!           This is the other way round to normal, but due to the fact that our eigenvalues are ordered -ve -> +ve
!           TODO: Are the signs of this matrix correct?
        allocate(X_stab(ov_space,ov_space)) !First index is (m,i) compound index. Second is the eigenvector index.
        allocate(Y_stab(ov_space,ov_space))
        X_stab(:,:)=0.0_dp
        Y_stab(:,:)=0.0_dp
        !Put the eigenvectors corresponding to *positive* eigenvalues into the X_stab and Y_stab arrays.
        X_stab(1:ov_space,1:ov_space)=temp2(1:ov_space,ov_space+1:StabilitySize)
        Y_stab(1:ov_space,1:ov_space)=-temp2(ov_space+1:StabilitySize,ov_space+1:StabilitySize)
        deallocate(temp2)

        !Normalize the eigenvectors appropriately
        do mu=1,ov_space
            norm=0.0_dp
            Y_norm = 0.0_dp
            X_norm = 0.0_dp
            do i=1,ov_space
                norm = norm + X_stab(i,mu)*X_stab(i,mu) - Y_stab(i,mu)*Y_stab(i,mu)
                Y_norm = Y_norm + Y_stab(i,mu)*Y_stab(i,mu)
                X_norm = X_norm + X_stab(i,mu)*X_stab(i,mu)
            enddo
            if(norm.le.0.0_dp) then
                write(6,*) "Norm^2 for vector ",mu," is: ",norm
                call stop_all(t_r,'norm undefined')
            endif
            norm = sqrt(norm)
            do i=1,ov_space
                X_stab(i,mu) = X_stab(i,mu)/norm
                Y_stab(i,mu) = Y_stab(i,mu)/norm
            enddo
            if(Y_norm.gt.X_norm/2.0_dp) then
                write(6,*) "Warning: hole amplitudes large for excitation: ",mu,    &
                    " Quasi-boson approximation breaking down."
                write(6,*) "Norm of X component: ",X_norm
                write(6,*) "Norm of Y component: ",Y_norm
            endif
        enddo
!            call writematrix(X_stab,'X',.true.)

        !Now check orthogonality 
        !call Check_XY_orthogonality(X_stab,Y_stab)

!            call writevector(W2,'Stab_eigenvalues')

        !Now check that we satisfy the original RPA equations
        !For the *positive* eigenvalue space (since we have extracted eigenvectors corresponding to this), check that:
!           ( A B ) (X) = E_v(X )
!           ( B A ) (Y)      (-Y)
        deallocate(temp)
        allocate(temp(StabilitySize,ov_space))
        temp(1:ov_space,1:ov_space) = X_stab(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = Y_stab(:,:)
        allocate(temp2(StabilitySize,ov_space))
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-6_dp) then
                    write(6,*) i,j,temp2(j,i),(W2(i+ov_space)*X_stab(j,i)),W2(i+ov_space)
                    call stop_all(t_r,"RPA equations not satisfied for positive frequencies in X matrix")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j+ov_space,i)-(-W2(i+ov_space)*Y_stab(j,i))).gt.1.0e-6_dp) then
                    write(6,*) i,j,temp2(j+ov_space,i),(-W2(i+ov_space)*Y_stab(j,i)),-W2(i+ov_space)
                    call stop_all(t_r,"RPA equations not satisfied for positive frequencies in Y matrix")
                endif
            enddo
        enddo
        deallocate(temp,temp2)

        !Is is also satisfied the other way around?
        !Check that we also satisfy (still for the *positive* eigenvalues):
!           ( A B ) (Y) = -E_v(Y )
!           ( B A ) (X)       (-X)
        allocate(temp(StabilitySize,ov_space))
        temp(1:ov_space,1:ov_space) = Y_stab(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = X_stab(:,:)
        allocate(temp2(StabilitySize,ov_space))
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(-W2(i+ov_space)*Y_stab(j,i))).gt.1.0e-6_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for negative frequencies in X matrix")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j+ov_space,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-6_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for negative frequencies in Y matrix")
                endif
            enddo
        enddo
        deallocate(temp,temp2)
        !TODO: Finally, check that we satisfy eq. 1 in the Scuseria paper for X and Y defined for positive eigenvalues...
        do i=ov_space+1,StabilitySize
            do j=1,StabilitySize
                StabilityCopy(i,j)=-StabilityCopy(i,j)
            enddo
        enddo
        !Stability copy is now (A B // -B -A)
        allocate(temp(StabilitySize,ov_space))
        allocate(temp2(StabilitySize,ov_space))
        temp=0.0_dp
        temp(1:ov_space,1:ov_space) = X_stab(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = Y_stab(:,:)
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-7) then
                    call stop_all(t_r,"RPA equations not satisfied for X")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=ov_space+1,StabilitySize
                if(abs(temp2(j,i)-(W2(i+ov_space)*Y_stab(j-ov_space,i))).gt.1.0e-7) then
                    call stop_all(t_r,"RPA equations not satisfied for Y")
                endif
            enddo
        enddo
        deallocate(temp,temp2)

        !Now calculate energy, in two different ways:
        !1. -1/2 Tr[A] + 1/2 sum_v E_v(positive)
        Energy_stab=0.0_dp
        do i=1,ov_space
            Energy_stab = Energy_stab + W2(ov_space+i) - A_mat(i,i)
        enddo
        Energy_stab = Energy_stab/2.0_dp

        write(6,"(A,G25.10)") "Full RPA energy from stability analysis (plasmonic RPA-TDA excitation energies): ",  &
            Energy_stab

        Energy_stab = 0.0_dp
        !E = 0.25 * Tr[BZ] where Z = Y X^-1

        allocate(temp2(ov_space,ov_space))
        temp2(:,:) = 0.0_dp
        !Find X^-1 
        call d_inv(X_stab,temp2)
!            call writematrix(temp2,'X^-1',.true.)
        allocate(temp(ov_space,ov_space))
        !Find Z (temp)
        call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,Y_stab,ov_space,temp2,ov_space,0.0_dp,temp,ov_space)
        !Find BZ (temp2)
        call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,B_mat,ov_space,temp,ov_space,0.0_dp,temp2,ov_space)
        !Take trace of BZ
        do i=1,ov_space
            Energy_stab = Energy_stab + temp2(i,i)
        enddo
        Energy_stab = Energy_stab/2.0_dp
        write(6,"(A,G25.10)") "Full RPA energy from stability analysis (Ring-CCD: 1/2 Tr[BZ]): ",Energy_stab

        Energy_stab = 0.0_dp
        do i=1,ov_space
            Y_norm = 0.0_dp
            do j=1,ov_space
                Y_norm = Y_norm + Y_stab(j,i)**2
            enddo
            Energy_stab = Energy_stab - W2(i+ov_space)*Y_norm
        enddo
        write(6,"(A,G25.10)") "Full RPA energy from stability analysis (Y-matrix): ",Energy_stab

        !Now, calculate the response functions for expectation value A and perturbation V
        !This is (for positive frequencies):
        !\sum_nu (<0|[A,Q_nu^+]|0><0|[Q_nu,V]|0> / (omega - W_nu)) - (<0|[V,Q_nu^+]|0><0|[Q_nu,A]|0> / (omega + W_nu))
        !Calculate the transition moments first, for a density-density response at site pertsite 
        allocate(trans_moment(ov_space))    
        trans_moment(:) = 0.0_dp

        !Construct an MO-AO orbital rotation matrix for spin-orbitals
        allocate(AOMO_Spin(nSites*2,nSites*2))
        AOMO_Spin(:,:) = 0.0_dp
        do i=1,nSites*2
            do j=1,nSites*2
                i_spat = gtid(i)
                j_spat = gtid(j)
                AOMO_Spin(j,i) = FullHFOrbs(gtid(j),gtid(i))
            enddo
        enddo
        allocate(DM(nSites*2,nSites*2))
        deallocate(temp)
        allocate(temp(nSites*2,nSites*2))

        write(6,*) "Calculating RPA Transition moments..."
        call flush(6)

        !Calculate <|[V,Q_nu^+]|0><0|[Q_nu,V]|0> and store for each nu
        do i=1,nel
            do a=virt_start,2*nSites
                if(mod(i,2).ne.mod(a,2)) cycle      !Only want same spin excitations 
                ai_ind = ov_space_spinind(a,i)      !This is the index in the array

                DM(:,:) = 0.0_dp
                if(mod(i,2).eq.1) then
                    !Alpha -> alpha transition.
                    !Parity is -1
                    DM(i,a) = -1.0_dp
                else
                    !Beta -> beta transition
                    !Parity is 1
                    DM(i,a) = 1.0_dp
                endif
                !Now, rotate this density matrix into the AO basis
                call dgemm('n','n',nSites*2,nSites*2,nSites*2,1.0_dp,AOMO_Spin,nSites*2,DM,nSites*2,0.0_dp,temp,nSites*2)
                call dgemm('n','t',nSites*2,nSites*2,nSites*2,1.0_dp,temp,nSites*2,AOMO_Spin,nSites*2,0.0_dp,DM,nSites*2)

                !Extract the pertsite,pertsite component, and sum the two spins
!                pertsite_alpha = pertsite*2 - 1
!                pertsite_beta = pertsite*2 
                DMEl1 = DM(pertsite*2-1,pertsite*2-1) + DM(pertsite*2,pertsite*2)

                !Now do <D_i^a|a_a^+ a_i|D_0> element, which is the one on the other side
                DM(:,:) = 0.0_dp
                if(mod(i,2).eq.1) then
                    DM(a,i) = -1.0_dp
                else
                    DM(i,a) = 1.0_dp
                endif
                !Now, rotate this density matrix into the AO basis
                call dgemm('n','n',nSites*2,nSites*2,nSites*2,1.0_dp,AOMO_Spin,nSites*2,DM,nSites*2,0.0_dp,temp,nSites*2)
                call dgemm('n','t',nSites*2,nSites*2,nSites*2,1.0_dp,temp,nSites*2,AOMO_Spin,nSites*2,0.0_dp,DM,nSites*2)

                !Extract the pertsite,pertsite component, and sum the two spins
                DMEl2 = DM(pertsite*2-1,pertsite*2-1) + DM(pertsite*2,pertsite*2)

                do excit=1,ov_space
                    trans_moment(excit) = trans_moment(excit) + ((X_stab(ai_ind,excit)*DMEl1 - Y_stab(ai_ind,excit)*DMEl2)*  &
                        (X_stab(ai_ind,excit)*DMEl2 - Y_stab(ai_ind,excit)*DMEl1))/4.0_dp   !Divide by 4 since each commutator is /2
                enddo
            enddo
        enddo
        trans_moment(:) = trans_moment(:)*Lambda

        iunit = get_free_unit()
        open(unit=iunit,file='RPA_Transitions',status='unknown')
        write(iunit,"(A)") "#Excitation     Transition_Frequency       Transition_Moment"
        do i=1,ov_space
            write(iunit,"(I8,2G22.12)") i,W2(ov_space+i),trans_moment(i)
        enddo
        close(iunit)
        write(6,*) "First 10 RPA transition frequencies: "
        highbound = min(ov_space,10)
        call writevector(W2(ov_space+1:ov_space+highbound),'transition frequencies')

        write(6,*) "Writing RPA linear response function to disk..."

        open(unit=iunit,file='RPA_DDResponse',status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ResponseFn = 0.0_dp
            do i=1,ov_space
                ResponseFn = ResponseFn + ((trans_moment(i)/(Omega-W2(ov_space+i))) - (trans_moment(i)/(Omega+W2(ov_space+i))))
            enddo
            write(iunit,*) Omega,ResponseFn

            Omega = Omega + Omega_Step

        enddo
        close(iunit)

        deallocate(W2,W,temp,temp2,StabilityCopy,Stability,A_mat,B_Mat,trans_moment,S_half)
        deallocate(X_stab,Y_stab,AOMO_Spin,DM)

    end subroutine RPA_LR

    !Only want to consider single excitation space, consisting of i -> a
    !First list all alpha excitations, then beta excitations
    !Within each spin-type, it is virtual fast
    integer function ov_space_spinind(a,i)
        implicit none
        integer, intent(in) :: i,a
        integer :: a_spat,i_spat,nVirt_spat,gtid

        if(mod(i,2).ne.mod(a,2)) ov_space_spinind = -1  !*Should* be easy to see where this goes wrong

        !Convert to spatial. Index the virtual excitations starting at 1
        a_spat = gtid(a-NEl)    !Runs from 1 -> number of spatial virtual orbitals
        i_spat = gtid(i)        !Runs from 1 -> nOcc
        nVirt_spat = nSites - nOcc

        if(mod(i,2).eq.1) then
            !It is an alpha -> alpha transition
            !These are indexed first
            ov_space_spinind = (i_spat-1)*nVirt_spat + a_spat
        else
            !It is a beta -> beta transition
            !Add on the entire set of alpha -> alpha transitions which come first
            ov_space_spinind = (i_spat-1)*nVirt_spat + a_spat + (nVirt_spat*nOcc)
        endif
    end function ov_space_spinind

    !dPsi/dLambda for Static DD response (i.e. omega -> 0)
    subroutine StaticMF_DD()
        real(dp) , allocatable :: Orbs(:,:),Energies(:),Work(:),PertHamil(:,:),PertOrbs(:,:)
        real(dp) , allocatable :: TempRDM(:,:),PertBath(:),GSBath(:),PertDM(:,:)
        real(dp) :: StaticResponse,Overlap,dStep,DDOT,PertNorm
        character(len=*), parameter :: t_r='StaticMF_DD'
        integer :: lWork, info,i,j
        
        allocate(Orbs(nSites,nSites))
        allocate(Energies(nSites))
        allocate(TempRDM(nSites,nSites))
        allocate(GSBath(nSites))
        Orbs(:,:) = h0(:,:)
        Energies(:) = 0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,Orbs,nSites,Energies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,Orbs,nSites,Energies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        
        !Determine GS bath orbital
        call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,Orbs(:,1:nOcc),nSites,Orbs(:,1:nOcc),nSites,0.0_dp,TempRDM,nSites)
        GSBath(:) = 0.0_dp
        GSBath(nImp+1:nSites) = TempRDM(nImp+1:nSites,1)
        PertNorm = DDOT(nSites,GSBath(:),1,GSBath(:),1)
        GSBath(:) = GSBath(:) / sqrt(PertNorm)
        deallocate(TempRDM)
        
        dStep = 0.01

        do while(.true.)

            dStep = dStep/2.0_dp

            if(dStep.lt.1.0e-8_dp) exit

            allocate(PertBath(nSites))
            allocate(PertHamil(nSites,nSites))
            allocate(PertDM(nSites,nSites))
            PertHamil(:,:) = h0(:,:)
            PertHamil(pertsite,pertsite) = PertHamil(pertsite,pertsite) + dStep
            Energies(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nSites,PertHamil,nSites,Energies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nSites,PertHamil,nSites,Energies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            !Determine perturbed bath orbital
            !Calc RDM
            call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,PertHamil(:,1:nOcc),nSites,    &
                PertHamil(:,1:nOcc),nSites,0.0_dp,PertDM,nSites)
            PertBath(:) = 0.0_dp
            PertBath(nImp+1:nSites) = PertDM(nImp+1:nSites,1)
            PertNorm = DDOT(nSites,PertBath(:),1,PertBath(:),1)
            PertBath(:) = PertBath(:) / sqrt(PertNorm)

            !Calculate derivative
            PertBath(:) = PertBath(:) - GSBath(:)
            PertBath(:) = PertBath(:) / dStep

            call writevector(PertBath,'Static Pert Bath')

            allocate(PertOrbs(nSites,nOcc))
            do i=1,nOcc
                PertOrbs(:,i) = (PertHamil(:,i) - Orbs(:,i))/dStep
            enddo

            !Matrix elements are now PertHamil(pertsite,i)
            StaticResponse = 0.0_dp
            do i=1,nOcc
                Overlap = 0.0_dp
                do j=1,nSites
                    Overlap = Overlap + PertOrbs(j,i)*Orbs(j,i)
                enddo
                write(6,*) "Overlap: ",i,Overlap
                StaticResponse = StaticResponse + Overlap*Lambda*PertOrbs(pertsite,i)*Orbs(pertsite,i)
            enddo

            write(6,*) "Mean field static resonse = ",StaticResponse,dStep

            deallocate(PertOrbs,PertHamil,PertBath,PertDM)

        enddo

        deallocate(Orbs,Energies,GSBath)

                    
    end subroutine StaticMF_DD


    !Just calculate the response at the mean-field level to compare to as U -> 0
    subroutine non_interactingLR()
        implicit none
        real(dp) :: MFDD_Response,EDiff,Omega
        integer :: n,a

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            MFDD_Response = 0.0_dp

            do n=1,nOcc
                do a=nOcc+1,nSites
                    EDiff = HFEnergies(a)-HFEnergies(n)
                    MFDD_Response = MFDD_Response + ((HFOrbs(pertsite,n)*HFOrbs(pertsite,a))**2)*Lambda*  &
                        ((1.0_dp/(EDiff-Omega))+(1.0_dp/(EDiff+Omega)))
                enddo
            enddo

            write(6,*) "Mean-field DD response: ",Omega,MFDD_Response

            Omega = Omega + Omega_Step

        enddo

    end subroutine non_interactingLR

                    
!    !Calculate density density response to perturbation of frequency omega at site pertsite 
!    subroutine calc_mf_dd_response()
!        implicit none
!        integer :: i,x,a,j
!        real(dp) :: CheckOrthog,DDOT    !,StepSize
!        real(dp) , allocatable :: temp(:,:),Pert(:,:),NormB0(:)
!        character(len=*) , parameter :: t_r='calc_mf_dd_response'
!
!        if(allocated(ResponseBasis)) deallocate(ResponseBasis)
!        allocate(ResponseBasis(nSites,2))
!
!        write(6,*) "Perturbation response for orbital: ",pertsite
!        write(6,*) "Frequency of perturbation: ",Omega
!        write(6,*) "Strength of perturbation: ",Lambda
!
!        if(nImp.gt.1) call stop_all(t_r,"Response not yet coded up for > 1 impurity site")
!
!        !The response operator is (\sum_{ia} V_ai |phi_a><psi_i| + V_ia|phi_i><phi_a|) / omega - (e_a-e_i)
!        !where V_ai = <phi_a|a_pertsite^+ a_pertsite|phi_i>
!
!        !The vector corresponding to this perturbation is calculated from the impurity to the environment sites
!        !Therefore, it is (c_imp,env_a)^(1) = <orb_imp| response operator | orb_env_a>
!
!
!        ResponseBasis(:,:) = 0.0_dp !Response over impurity sites = 0
!        do x=nImp+1,nSites
!            do i=1,nOcc
!                do a=nOcc+1,nSites
!                    !This is <phi_a| a_pertsite^+ a_pertsite |phi_i> * <imp|phi_a><phi_i|x>/omega-(e_a-e_i)
!                    ResponseBasis(x,2) = ResponseBasis(x,2) + (HFOrbs(pertsite,a)*HFOrbs(pertsite,i)*HFOrbs(1,a)*HFOrbs(x,i)/ &
!                        ((HFEnergies(a)-HFEnergies(i)) - omega))
!                    ResponseBasis(x,2) = ResponseBasis(x,2) + (HFOrbs(pertsite,i)*HFOrbs(pertsite,a)*HFOrbs(1,a)*HFOrbs(x,i)/ &
!                        ((HFEnergies(a)-HFEnergies(i)) + omega))
!                enddo
!            enddo
!        enddo
!
!        !Analytically calculate new bath orbital
!        !Renormalize the change in the first order bath orbital, so that it overall noramlized (to 1st order)
!        ResponseBasis(:,2) = ResponseBasis(:,2) / sqrt(ZerothBathNorm)
!
!!        !Add the newly normalized zeroth order orbital - do we need to do that if we just want the first-order change?
!!
!!!        ResponseBasis(:,2) = ResponseBasis(:,2) + EmbeddedBasis(:,2)*   &
!!!            (1.0_dp - DDOT(nSites,EmbeddedBasis(:,2),1,ResponseBasis(:,2),1)/ZerothBathNorm)
!!
!!
!!
!!
!!       !Numerically differentiate
!!        StepSize = 0.0001
!!
!!        ResponseBasis(:,2) = ResponseBasis(:,2) * StepSize 
!!
!!        ResponseBasis(1:nSites,2) = ResponseBasis(1:nSites,2) + EmbeddedBasis(:,2)  !Add original bath orbital
!!
!!        ResponseBasis(1:nSites,2) = ResponseBasis(1:nSites,2) - EmbeddedBasis(:,2)
!!        ResponseBasis(:,2) = ResponseBasis(:,2) / StepSize
!
!        call writevector(ResponseBasis(:,2),'ResponseBasis')
!
!        CheckOrthog = DDOT(nSites,ResponseBasis(:,2),1,ResponseBasis(:,2),1)
!        write(6,*) "norm: ",CheckOrthog
!
!        !ResponseBasis is now the bath orbital for first order change in the MF solution
!        !It should be orthogonal to the original bath orbital 
!        !However, since we have got a misture of the first and second order orbital in the solution, we have to project out the first
!        !order bath orbital from the original bath
!        !B^(0)/norm[B^(0)] * (1 - <B^(0)|B^(1)>/<B^(0)|B^(0)>
!        !We have to 'unnormalize' the states
!        CheckOrthog = DDOT(nSites,EmbeddedBasis(:,2)*sqrt(ZerothBathNorm),1,ResponseBasis(:,2)*sqrt(ZerothBathNorm),1)
!        CheckOrthog = 1.0_dp - CheckOrthog/ZerothBathNorm
!        allocate(NormB0(nSites))
!        NormB0(:) = EmbeddedBasis(:,2)*CheckOrthog
!        !only *now* can we correctly check for orthogonality
!        CheckOrthog = DDOT(nSites,NormB0(:),1,ResponseBasis,1)
!        write(6,*) "Projection against other bath: ",CheckOrthog
!        deallocate(NormB0)
!
!        !Add the impurity orbital to zero. We don't want to include impurity -> impurity or impurity -> bath^(0) coupling 
!        ResponseBasis(:,1) = 0.0_dp
!!        ResponseBasis(1,1) = 1.0_dp
!
!        !Now calculate the one-electron perturbations
!        !The standard 1-electron perturbation is 1/2 Lambda a_pertsite^+ a_pertsite.
!        !We calculate this first in the HF basis, and then transform into the zeroth-order embedding basis
!        if(allocated(Pert)) deallocate(Pert)
!        allocate(Pert(nSites,nSites))
!        allocate(temp(EmbSize,nSites))
!
!        if(allocated(Emb_Pert)) deallocate(Emb_Pert)
!        allocate(Emb_Pert(EmbSize,EmbSize))
!        if(allocated(Emb_h1)) deallocate(Emb_h1)
!        allocate(Emb_h1(EmbSize,EmbSize))
!
!        Pert(:,:) = 0.0_dp
!        do i=1,nSites
!            do j=1,nSites
!                Pert(i,j) = HFOrbs(pertsite,i)*HFOrbs(pertsite,j)
!            enddo
!        enddo
!        !Transform into embedding basis
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,Pert,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_Pert,EmbSize)
!
!        !Now we need to calculate H^(1)
!        !Transform h0 into the embedding basis
!        !We want C^(1)T h0 C^(0) + C^(0)T h0 C^(1)
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,ResponseBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_h1,EmbSize)
!
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,ResponseBasis,nSites,1.0_dp,Emb_h1,EmbSize)
!
!        call writematrix(Emb_h1,'Emb_H1',.true.)
!
!        !We now have the perturbations delta H and V in the embedding basis.
!        deallocate(temp,Pert)
!
!    end subroutine calc_mf_dd_response

!test schmidt decomposition
    subroutine test()
        implicit none
        real(dp), allocatable :: Simp_Eigenvals(:),ImpOverlap(:,:)
        real(dp), allocatable :: Senv_Eigenvals(:),EnvOverlap(:,:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: RotEnv(:,:),RotImp(:,:)
        integer :: lWork,info,i
        character(len=*), parameter :: t_r='test'


        allocate(ImpOverlap(nOcc,nOcc))

        !Extract the components of the MOs on the system sites, and construct the overlap matrix from them
!        call DGEMM('N','T',nImp,nImp,nOcc,1.0_dp,HFOrbs(1:nImp,1:nOcc),nImp,HFOrbs(1:nImp,1:nOcc),nImp,0.0_dp,ImpOverlap,nImp)
        call DGEMM('T','N',nOcc,nOcc,nImp,1.0_dp,HFOrbs(1:nImp,1:nOcc),nImp,HFOrbs(1:nImp,1:nOcc),nImp,0.0_dp,ImpOverlap,nOcc)

        allocate(Simp_Eigenvals(nOcc))
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nOcc,ImpOverlap,nOcc,Simp_Eigenvals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nOcc,ImpOverlap,nOcc,Simp_Eigenvals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        
        do i=1,nOcc
            if(Simp_Eigenvals(i).lt.1.0e-10_dp) then
                Simp_Eigenvals(i) = 0.0_dp
                ImpOverlap(:,i) = 0.0_dp
            else
                ImpOverlap(:,i) = ImpOverlap(:,i) / sqrt(Simp_Eigenvals(i))
            endif
        enddo

        call writevector(Simp_Eigenvals,'Simp_Eigenvals')
        call writematrix(ImpOverlap,'S_imp Eigenvectors',.true.)

        !Now do the same for the enivronment sites
        allocate(EnvOverlap(nOcc,nOcc))

!        call DGEMM('N','T',nSites-nImp,nSites-nImp,nOcc,1.0_dp,HFOrbs(nImp+1:nSites,1:nOcc),nSites-nImp,    &
!            HFOrbs(nImp+1:nSites,1:nOcc),nSites-nImp,0.0_dp,EnvOverlap,nSites-nImp)
        call DGEMM('T','N',nOcc,nOcc,nSites-nImp,1.0_dp,HFOrbs(nImp+1:nSites,1:nOcc),nSites-nImp,    &
            HFOrbs(nImp+1:nSites,1:nOcc),nSites-nImp,0.0_dp,EnvOverlap,nOcc)

        !Diagonalize S_env
        allocate(Senv_Eigenvals(nOcc))
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nOcc,EnvOverlap,nOcc,Senv_Eigenvals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nOcc,EnvOverlap,nOcc,Senv_Eigenvals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        do i=1,nOcc
            if(Senv_Eigenvals(i).lt.1.0e-10_dp) then
                Senv_Eigenvals(i) = 0.0_dp
                EnvOverlap(:,i) = 0.0_dp
            else
                EnvOverlap(:,i) = EnvOverlap(:,i) / sqrt(Senv_Eigenvals(i))
            endif
        enddo

        !These represent the vectors in the MO (i.e. occupied) space
        call writevector(Senv_Eigenvals,'Senv_Eigenvals')
        call writematrix(EnvOverlap,'S_env Eigenvectors',.true.)

        !Now, rotate the determinant, such that there is only one orbital spanning both sites
        !This will give the final vectors in the AO space
        !First, rotate the impurity sites:
        allocate(RotImp(nImp,nOcc))
        call DGEMM('N','N',nImp,nOcc,nOcc,1.0_dp,HFOrbs(1:nImp,1:nOcc),nImp,    &
            ImpOverlap,nOcc,0.0_dp,RotImp,nImp)
        call writematrix(RotImp,'Rotated Impurity sites',.true.)

        allocate(RotEnv(nSites-nImp,nOcc))
        call DGEMM('N','N',nSites-nImp,nOcc,nOcc,1.0_dp,HFOrbs(nImp+1:nSites,1:nOcc),   &
            nSites-nImp,EnvOverlap,nOcc,0.0_dp,RotEnv,nSites-nImp)
        call writematrix(RotEnv,'Rotated Environment sites',.true.)

!        allocate(OverallR(nOcc,nOcc))
!        OverallR(:,:) = 0.0_dp
!        RotMat(:,:) = 0.0_dp
!        RotMat(1:nImp,1:nImp) = ImpOverlap(:,:)
!        RotMat(nImp+1:nSites,nImp+1:nSites) = EnvOverlap(:,:)

        !Now rotate the original determinant


    end subroutine test
end module LinearResponse
