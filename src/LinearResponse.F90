module LinearResponse
    use const
    use errors, only: stop_all
    use mat_tools, only: WriteVector,WriteMatrix
    use globals
    implicit none
    contains

    !Solve the response equations in the basis of Psi^(0) + all single excits.
    !This is the basis of Psi^(0), the basis of internally contracted single excitations of it into
    !the virtual space, the basis of single excitations of the core into virtual space, and the basis
    !of single excitations of core into active space. This will all be constructed explicitly initially.
    subroutine SolveDMETResponse()
        use DetToolsData, only: nFCIDet,FCIDetList
        implicit none
        integer :: nCoreVirt,nCoreActive,nActiveVirt,nLinearSystem,ierr,info
        integer :: i,j,CoreNEl,ind2,ind1,a,x,Ex(2),b
        logical :: tSign,tSign_a,tSign_b
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


    end subroutine SolveDMETResponse

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
