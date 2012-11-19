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
    !   Work out self-consistency condition to optimise both the full MCLR and the fully contracted type
    !   Consider using different methods to obtain contraction coefficients in the fully contracted methods - TDA/RPA. Does this improve things?
    !   Perhaps look at CC2 to get a deeper understanding
    subroutine MR_LinearResponse()
        implicit none

        !Create contracted single excitation space using the non-interacting reference for the contractions
        !The matrix is then created in a CI fashion
!        call NonIntContracted_TDA_MCLR()

        !Create contracted single excitation space using the non-interacting reference for the contractions
        !The matrix is then created in an RPA fashion
!        call NonIntContracted_RPA_MCLR()

        !Full MCLR, creating excitations in a CI fashion, rather than with commutators. Should reduce to TDA in single reference limit
        call TDA_MCLR()

        !Full MCLR, with excitations and deexcitations. Should reduce to RPA in single reference limit
!        call RPA_MCLR()

    end subroutine MR_LinearResponse

    !Run single reference linear response calculations, based on true HF calculation.
    subroutine SR_LinearResponse()
        implicit none
        
        !Single reference RPA
        call RPA_LR()
        !Single reference TDA
!        call TDA_LR()
        !Non-interacting linear response
!        call NonInteractingLR()

    end subroutine SR_LinearResponse


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
    subroutine FindSchmidtPert()
        implicit none
        real(dp), allocatable :: HFPertBasis(:,:),temp(:,:)
        real(dp) :: EDiff
        integer :: a,i
        character(len=*), parameter :: t_r='FindSchmidtBasis'

        if(allocated(HFPertBasis)) deallocate(HFPertBasis)
        allocate(HFPertBasis(nSites,nSites))
        HFPertBasis(:,:) = 0.0_dp

        !Assume perturbation is local to the first impurity site (pertsite = 1).
        if(pertsite.ne.1) call stop_all(t_r,'Perturbation is not local to the impurity site')
        !Assuming that MF DD operator is V_ia/e_a-e_i-w + V_ia/e_a-e_i+w
        do i=1,nOcc
            do a=nOcc+1,nSites
                EDiff = HFEnergies(a)-HFEnergies(i)
                HFPertBasis(i,a) = HFOrbs(pertsite,i)*HFOrbs(pertsite,a)*Lambda*((1.0_dp/(EDiff-Omega))+(1.0_dp/(EDiff+Omega)))
                HFPertBasis(a,i) = HFPertBasis(i,a)
            enddo
        enddo

        !write(6,*) "Mean-field DD response: ",Omega,MFDD_Response

        !write(6,*) "Transforming non-interacting response operator into full schmidt basis..."

        allocate(temp(nSites,nSites))
        if(allocated(SchmidtPert)) deallocate(SchmidtPert)
        allocate(SchmidtPert(nSites,nSites))
        call DGEMM('T','N',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,HFPertBasis,nSites,0.0_dp,temp,nSites)
        call DGEMM('N','N',nSites,nSites,nSites,1.0_dp,temp,nSites,HFPertBasis,nSites,0.0_dp,SchmidtPert,nSites)
        deallocate(temp,HFPertBasis)

        !SchmidtPert is now the perturbation in the schmidt basis
        !call writematrix(SchmidtPert,'Perturbation in schmidt basis',.true.)

    end subroutine FindSchmidtPert
    
    !Set up the RPA equations and solve them.
    !Finally, create the density-density linear response function from the resulting excitations/deexcitations
    !This is done in the spin-orbital space
    subroutine RPA_LR()
        use utils, only: get_free_unit
        use matrixops, only: d_inv
        implicit none
        integer :: ov_space,virt_start,ierr,j,ex(2,2),ex2(2,2),n,i,m,nj_ind,mi_ind,info,lwork
        integer :: m_spat,i_spat,StabilitySize,mu,gtid,j_spat,ai_ind,iunit,a,excit
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
                        (X_stab(ai_ind,excit)*DMEl2 - Y_stab(ai_ind,excit)*DMEl1))
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

                    
    !Calculate density density response to perturbation of frequency omega at site pertsite 
    subroutine calc_mf_dd_response()
        implicit none
        integer :: i,x,a,j
        real(dp) :: CheckOrthog,DDOT    !,StepSize
        real(dp) , allocatable :: temp(:,:),Pert(:,:),NormB0(:)
        character(len=*) , parameter :: t_r='calc_mf_dd_response'

        if(allocated(ResponseBasis)) deallocate(ResponseBasis)
        allocate(ResponseBasis(nSites,2))

        write(6,*) "Perturbation response for orbital: ",pertsite
        write(6,*) "Frequency of perturbation: ",Omega
        write(6,*) "Strength of perturbation: ",Lambda

        if(nImp.gt.1) call stop_all(t_r,"Response not yet coded up for > 1 impurity site")

        !The response operator is (\sum_{ia} V_ai |phi_a><psi_i| + V_ia|phi_i><phi_a|) / omega - (e_a-e_i)
        !where V_ai = <phi_a|a_pertsite^+ a_pertsite|phi_i>

        !The vector corresponding to this perturbation is calculated from the impurity to the environment sites
        !Therefore, it is (c_imp,env_a)^(1) = <orb_imp| response operator | orb_env_a>


        ResponseBasis(:,:) = 0.0_dp !Response over impurity sites = 0
        do x=nImp+1,nSites
            do i=1,nOcc
                do a=nOcc+1,nSites
                    !This is <phi_a| a_pertsite^+ a_pertsite |phi_i> * <imp|phi_a><phi_i|x>/omega-(e_a-e_i)
                    ResponseBasis(x,2) = ResponseBasis(x,2) + (HFOrbs(pertsite,a)*HFOrbs(pertsite,i)*HFOrbs(1,a)*HFOrbs(x,i)/ &
                        ((HFEnergies(a)-HFEnergies(i)) - omega))
                    ResponseBasis(x,2) = ResponseBasis(x,2) + (HFOrbs(pertsite,i)*HFOrbs(pertsite,a)*HFOrbs(1,a)*HFOrbs(x,i)/ &
                        ((HFEnergies(a)-HFEnergies(i)) + omega))
                enddo
            enddo
        enddo

        !Analytically calculate new bath orbital
        !Renormalize the change in the first order bath orbital, so that it overall noramlized (to 1st order)
        ResponseBasis(:,2) = ResponseBasis(:,2) / sqrt(ZerothBathNorm)

!        !Add the newly normalized zeroth order orbital - do we need to do that if we just want the first-order change?
!
!!        ResponseBasis(:,2) = ResponseBasis(:,2) + EmbeddedBasis(:,2)*   &
!!            (1.0_dp - DDOT(nSites,EmbeddedBasis(:,2),1,ResponseBasis(:,2),1)/ZerothBathNorm)
!
!
!
!
!       !Numerically differentiate
!        StepSize = 0.0001
!
!        ResponseBasis(:,2) = ResponseBasis(:,2) * StepSize 
!
!        ResponseBasis(1:nSites,2) = ResponseBasis(1:nSites,2) + EmbeddedBasis(:,2)  !Add original bath orbital
!
!        ResponseBasis(1:nSites,2) = ResponseBasis(1:nSites,2) - EmbeddedBasis(:,2)
!        ResponseBasis(:,2) = ResponseBasis(:,2) / StepSize

        call writevector(ResponseBasis(:,2),'ResponseBasis')

        CheckOrthog = DDOT(nSites,ResponseBasis(:,2),1,ResponseBasis(:,2),1)
        write(6,*) "norm: ",CheckOrthog

        !ResponseBasis is now the bath orbital for first order change in the MF solution
        !It should be orthogonal to the original bath orbital 
        !However, since we have got a misture of the first and second order orbital in the solution, we have to project out the first
        !order bath orbital from the original bath
        !B^(0)/norm[B^(0)] * (1 - <B^(0)|B^(1)>/<B^(0)|B^(0)>
        !We have to 'unnormalize' the states
        CheckOrthog = DDOT(nSites,EmbeddedBasis(:,2)*sqrt(ZerothBathNorm),1,ResponseBasis(:,2)*sqrt(ZerothBathNorm),1)
        CheckOrthog = 1.0_dp - CheckOrthog/ZerothBathNorm
        allocate(NormB0(nSites))
        NormB0(:) = EmbeddedBasis(:,2)*CheckOrthog
        !only *now* can we correctly check for orthogonality
        CheckOrthog = DDOT(nSites,NormB0(:),1,ResponseBasis,1)
        write(6,*) "Projection against other bath: ",CheckOrthog
        deallocate(NormB0)

        !Add the impurity orbital to zero. We don't want to include impurity -> impurity or impurity -> bath^(0) coupling 
        ResponseBasis(:,1) = 0.0_dp
!        ResponseBasis(1,1) = 1.0_dp

        !Now calculate the one-electron perturbations
        !The standard 1-electron perturbation is 1/2 Lambda a_pertsite^+ a_pertsite.
        !We calculate this first in the HF basis, and then transform into the zeroth-order embedding basis
        if(allocated(Pert)) deallocate(Pert)
        allocate(Pert(nSites,nSites))
        allocate(temp(EmbSize,nSites))

        if(allocated(Emb_Pert)) deallocate(Emb_Pert)
        allocate(Emb_Pert(EmbSize,EmbSize))
        if(allocated(Emb_h1)) deallocate(Emb_h1)
        allocate(Emb_h1(EmbSize,EmbSize))

        Pert(:,:) = 0.0_dp
        do i=1,nSites
            do j=1,nSites
                Pert(i,j) = HFOrbs(pertsite,i)*HFOrbs(pertsite,j)
            enddo
        enddo
        !Transform into embedding basis
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,Pert,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_Pert,EmbSize)

        !Now we need to calculate H^(1)
        !Transform h0 into the embedding basis
        !We want C^(1)T h0 C^(0) + C^(0)T h0 C^(1)
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,ResponseBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_h1,EmbSize)

        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,ResponseBasis,nSites,1.0_dp,Emb_h1,EmbSize)

        call writematrix(Emb_h1,'Emb_H1',.true.)

        !We now have the perturbations delta H and V in the embedding basis.
        deallocate(temp,Pert)

    end subroutine calc_mf_dd_response

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
