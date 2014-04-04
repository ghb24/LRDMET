module SchmidtDecomp
    use const
    use timing
    use errors, only: stop_all
    use globals
    use mat_tools, only: WriteVector,WriteMatrix,WriteMatrixComp
    implicit none

    contains

    !Schmidt decompose a determinant resulting from a complex one-electron hamiltonian passed in, including virtual and core spaces.
    subroutine SchmidtDecompose_C(matc)
        use mat_tools, only: DiagOneeOp
        implicit none
        complex(dp), intent(in) :: matc(nSites,nSites)

        complex(dp), allocatable :: orbs(:,:),ProjOverlap(:,:),cWork(:)
        complex(dp), allocatable :: RotOccOrbs(:,:),ImpurityOrbs(:,:),temp(:,:)
        complex(dp), allocatable :: ProjOverlapVirt(:,:),OverlapVirt(:,:)
        complex(dp), allocatable :: VirtSpace(:,:),HFtoSchmidtTransform_c(:,:)
        complex(dp), allocatable :: EmbeddedBasis_c(:,:),h0_c(:,:)
        complex(dp) :: normc,Overlap
        real(dp) :: norm
        real(dp), allocatable :: energies(:),ProjOverlapEVals(:),rWork(:)
        integer :: i,j,k,lWork,info,nbath,nVirt
        character(len=*), parameter :: t_r='SchmidtDecompose_C'

        write(6,"(A)") "Schmidt decomposing new complex matrix..."

        if(tUHF) call stop_all(t_r,'Cannot currently deal with UHF')

        !First, diagonalize the hamiltonian for a set of orbitals.
        if(tSCFHF) then
            call stop_all(t_r,'Not currently set up for HF-SCF')
        endif
        if(tChempot.or.tAnderson) call stop_all(t_r,'Cannot deal with anderson models')

        if(tCheck) then
            !Check matrix is hermitian
            do i = 1,nSites
                do j = i+1,nSites
                    if(abs(matc(i,j)-dconjg(matc(j,i))).gt.1.0e-8_dp) then
                        call writematrixcomp(matc,'matrix',.true.)
                        call stop_all(t_r,'matrix not hermitian')
                    endif
                enddo
                if(abs(aimag(matc(i,i))).gt.1.0e-8_dp) then
                    call stop_all(t_r,'matrix not diagonal real')
                endif
            enddo
        endif

        allocate(orbs(nSites,nSites))
        allocate(energies(nSites))
        orbs(:,:) = matc(:,:)
        energies(:) = zero
        if(tDiag_kspace) then
            call DiagOneEOp(Orbs,Energies,nImp,nSites,tDiag_kspace)
        else
            call DiagOneEOp(Orbs,Energies,nImp,nSites,.false.)
        endif

        write(6,*) "nOCC", nOcc
        write(6,*) "Fock eigenvalues around fermi level: "
        do i=max(1,nOcc-7),nOcc
            write(6,*) Energies(i),"*"
        enddo
        do i=nOcc+1,min(nSites,nOcc+7)
            write(6,*) Energies(i)
        enddo

        write(6,"(A,F20.10)") "Chemical potential is: ",(Energies(nOcc) + Energies(nOcc+1))/2.0_dp
        write(6,"(A,F20.10)") "HL Gap is: ",Energies(nOcc+1)-Energies(nOcc)

        !Now we have (complex) orbitals, lets calculate the overlap matrix
        allocate(ProjOverlap(nOcc,nOcc))
        call ZGEMM('C','N',nOcc,nOcc,nImp,zone,Orbs(1:nImp,1:nOcc),nImp,Orbs(1:nImp,1:nOcc),nImp,zzero,ProjOverlap,nOcc)
        
        if(tCheck) then
            !Check overlap matrix is hermitian
            do i = 1,nOcc
                do j = i+1,nOcc
                    if(abs(ProjOverlap(i,j)-dconjg(ProjOverlap(j,i))).gt.1.0e-8_dp) then
                        call writematrixcomp(ProjOverlap,'projected overlap matrix for occupied orbitals',.true.)
                        write(6,*) "i,j: ",i,j
                        write(6,*) "ProjOverlap(i,j): ",ProjOverlap(i,j)
                        write(6,*) "ProjOverlap(j,i): ",ProjOverlap(j,i)
                        call stop_all(t_r,'matrix not hermitian')
                    endif
                enddo
                if(abs(aimag(ProjOverlap(i,i))).gt.1.0e-8_dp) then
                    call stop_all(t_r,'matrix not diagonal real')
                endif
            enddo
        endif

        !Diagonalize this overlap matrix
        allocate(ProjOverlapEVals(nOcc))
        ProjOverlapEVals(:) = zero
        allocate(cWork(1))
        allocate(rWork(max(1, 3*nOcc-2)))
        lWork=-1
        info=0
        call zheev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,cWork,lWork,rWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
        lwork=int(abs(cwork(1)))+1
        deallocate(cwork)
        allocate(cwork(lwork))
        call zheev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,cWork,lWork,rWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(cWork,rWork)

        if(tWriteout) then
            call writevector(ProjOverlapEVals,'Projected overlap eigenvalues')
        endif

        !We should only have nImp non-zero evals
        nbath = 0
        do i = 1,nOcc
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.ne.nImp) call stop_all(t_r,'error here')

        !Now rotate original occ space into new entangled basis
        allocate(RotOccOrbs(nSites,nOcc))
        call ZGEMM('N','N',nSites,nOcc,nOcc,zone,Orbs(1:nSites,1:nOcc),nSites,ProjOverlap,nOcc,zzero,RotOccOrbs,nSites)

        !Construct bath by projecting out impurity and renormalizing
        do i = nOcc,nOcc-nImp+1,-1
            RotOccOrbs(1:nImp,i) = zzero
            normc = zzero
            do j = 1,nSites
                normc = normc + dconjg(RotOccOrbs(j,i))*RotOccOrbs(j,i)
            enddo
            if(abs(aimag(normc)).gt.1.0e-8_dp) then
                write(6,*) "norm for bath orbitals: ",normc
                call stop_all(t_r,'complex norm?')
            endif
            norm = sqrt(real(normc,dp))
            RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
        enddo

        !These states are now the bath states
!        call writematrixcomp(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)
        
        allocate(ImpurityOrbs(nSites,nImp))
        ImpurityOrbs(:,:) = zzero
        do i = 1,nImp
            ImpurityOrbs(i,i) = zone
        enddo

        !We now have all the orbitals. Which are orthogonal to which?
        do i = nOcc,nOcc-nImp+1,-1
            do j = nOcc,nOcc-nImp+1,-1
                Overlap = zzero
                do k = 1,nSites
                    Overlap = Overlap + dconjg(RotOccOrbs(k,i))*RotOccOrbs(k,j)
                enddo
!                Overlap = ZDOTC(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(i.eq.j) then
                    if(abs(Overlap-zone).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Bath orbitals not normalized set')
                    endif
                else
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Bath orbitals not orthogonal set')
                    endif
                endif
            enddo
        enddo

        !Now onto 'core' orbitals
        do i = 1,nOcc-nImp
            do j = nOcc,nOcc-nImp+1,-1
                !Overlap of core (i) with bath (j)
                Overlap = zzero
                do k = 1,nSites
                    Overlap = Overlap + dconjg(RotOccOrbs(k,i))*RotOccOrbs(k,j)
                enddo
!                Overlap = ZDOTC(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'bath orbitals with core not orthogonal set')
                endif
            enddo
        enddo

        !Need to sort virtual space, and remove redundancy
        nVirt = nSites-nOcc
        allocate(ProjOverlapVirt(nVirt,nVirt))
        
        !This array is used to calculate the overlap of each virtual space function with each impurity
        !function
        allocate(OverlapVirt(2*nImp,nVirt))
        !The first nImp correspond to the impurity orbital, and the next two correspond to the bath orbitals
        OverlapVirt(:,:) = zzero 
        do i = nOcc+1,nSites  !run through virtual space
            do j = 1,nImp     !run through impurity orbitals
                OverlapVirt(j,i-nOcc) = Orbs(j,i)
            enddo
        enddo
        
        call ZGEMM('C','N',nImp,nVirt,nSites,zone,RotOccOrbs(:,nOcc-nImp+1:nOcc),nSites,Orbs(:,nOcc+1:nSites), &
            nSites,zzero,OverlapVirt(nImp+1:2*nImp,1:nVirt),2*nImp-nImp)

        !Combine overlaps to get full projected overlap matrix
        call ZGEMM('C','N',nVirt,nVirt,2*nImp,zone,OverlapVirt,2*nImp,OverlapVirt,2*nImp,zzero,ProjOverlapVirt,nVirt)

        !Diagonalize virtual projected overlap matrix
        deallocate(ProjOverlapEVals)
        allocate(ProjOverlapEVals(nVirt))
        ProjOverlapEVals(:) = zero
        allocate(cWork(1))
        allocate(rWork(max(1, 3*nOcc-2)))
        lWork=-1
        info=0
        call zheev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,cWork,lWork,rWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
        lwork=int(abs(cwork(1)))+1
        deallocate(cwork)
        allocate(cwork(lwork))
        call zheev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,cWork,lWork,rWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(cWork,rWork)

        nbath = 0   !Count the number of virtual functions which span the space of the embedded system
        do i=1,nVirt
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.ne.nImp) then
            call stop_all(t_r,'Virtual space redundancy not as expected')
        endif
        
        !Now rotate orbitals such that they are orthogonal, while deleting the redundant ones
        !Assume the last nImp are redundant
        allocate(VirtSpace(nSites,nVirt-nImp))
        call ZGEMM('N','N',nSites,nVirt-nImp,nVirt,zone,Orbs(:,nOcc+1:nSites),nSites,   &
            ProjOverlapVirt(1:nVirt,1:nVirt-nImp),nVirt,zzero,VirtSpace,nSites)

        !Check now that the virtual space is now orthogonal to all occupied HF orbitals, as well as the impurity and bath sites

        !First against the impurity sites
        do i=1,nImp
            do j=1,nVirt-nImp
                Overlap = zzero
                do k = 1,nSites
                    Overlap = Overlap + dconjg(ImpurityOrbs(k,i))*VirtSpace(k,j)
                enddo
!                Overlap = ZDOTC(nSites,ImpurityOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'virtual orbitals not orthogonal to impurity orbitals')
                endif
            enddo
        enddo

        do i=1,nOcc
            do j=1,nVirt-nImp
                Overlap = zzero
                do k = 1,nSites
                    Overlap = Overlap + dconjg(RotOccOrbs(k,i))*VirtSpace(k,j)
                enddo
!                Overlap = ZDOTC(nSites,RotOccOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    if(i.gt.(nOcc-nImp)) then
                        call stop_all(t_r,'virtual orbitals not orthogonal to bath orbitals')
                    else
                        call stop_all(t_r,'virtual orbitals not orthogonal to core orbitals')
                    endif
                endif
            enddo
        enddo

        if(allocated(FullSchmidtBasis_c)) deallocate(FullSchmidtBasis_c)
        allocate(FullSchmidtBasis_c(nSites,nSites))   ! (Atomicbasis,SchmidtBasis)
        FullSchmidtBasis_c(:,:) = zzero
        FullSchmidtBasis_c(:,1:nOcc-nImp) = RotOccOrbs(:,1:nOcc-nImp)    !The first orbitals consist of core  
        FullSchmidtBasis_c(:,nOcc-nImp+1:nOcc) = ImpurityOrbs(:,:)    !Impurity Orbs
        FullSchmidtBasis_c(:,nOcc+1:nOcc+nImp) = RotOccOrbs(:,nOcc-nImp+1:nOcc)   !Bath
        FullSchmidtBasis_c(:,nOcc+nImp+1:nSites) = VirtSpace(:,:)     !Virtual space
        
        allocate(EmbeddedBasis_c(nSites,2*nImp))
        EmbeddedBasis_c(:,:) = zzero  
        EmbeddedBasis_c(:,1:nImp) = ImpurityOrbs(:,:)
        EmbeddedBasis_c(:,nImp+1:2*nImp) = FullSchmidtBasis_c(:,nOcc+1:nOcc+nImp)

        !call writematrix(FullSchmidtBasis,'FullSchmidtBasis',.true.)

        !Construct unitary basis transformation matrix from HF to embedded basis
        allocate(HFtoSchmidtTransform_c(nSites,nSites))
        call ZGEMM('C','N',nSites,nSites,nSites,zone,Orbs,nSites,FullSchmidtBasis_c,nSites,zzero,HFtoSchmidtTransform_c,nSites)

        allocate(temp(nSites,nSites))
        !Check that this operator is unitary
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

        !Use this to rotate the fock operator into this new basis
        if(allocated(FockSchmidt_c)) deallocate(FockSchmidt_c)
        allocate(FockSchmidt_c(nSites,nSites))
        !Set up FockSchmidt to temporarily to be the HF basis fock operator (i.e. diagonal)
        FockSchmidt_c(:,:) = zero
        do i=1,nSites
            FockSchmidt_c(i,i) = cmplx(Energies(i),zero,dp)
        enddo
        call ZGEMM('C','N',nSites,nSites,nSites,zone,HFtoSchmidtTransform_c,nSites,FockSchmidt_c,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,HFtoSchmidtTransform_c,nSites,zzero,FockSchmidt_c,nSites)
            
        deallocate(temp)

        !do i=1,nSites
        !    write(6,*) "FOCKSCHMIDT: ",i,FockSchmidt(i,i)
        !enddo
        !FockSchmidt has been overwritten with the fock matrix in the schmidt basis
!        call writematrix(FockSchmidt,'Fock in schmidt basis',.true.)

!       Calculate the non-interacting core energy of the DMET wavefunction
        CoreEnergy = zero
        do i = 1,nOcc-nImp
            if(abs(aimag(FockSchmidt_c(i,i))).gt.1.0e-7_dp) then
                write(6,*) "i: ",i
                write(6,*) "FockSchmidt_c(i,i): ",FockSchmidt_c(i,i)
                call writematrixcomp(FockSchmidt_c,'FockSchmidt',.true.)
                call stop_all(t_r,'complex energies')
            endif
            CoreEnergy = CoreEnergy + real(FockSchmidt(i,i),dp)
        enddo
        write(6,*) "Non-interacting core energy for DMET ground state wavefunction is: ",CoreEnergy


        !Now to also calculate Emb_h0_c and Emb_h0v_c
        if(allocated(Emb_h0_c)) deallocate(Emb_h0_c)
        if(allocated(Emb_h0v_c)) deallocate(Emb_h0v_c)
        allocate(Emb_h0_c(EmbSize,EmbSize))
        allocate(Emb_h0v_c(EmbSize,EmbSize))

        allocate(h0_c(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                h0_c(j,i) = cmplx(h0(j,i),zero,dp)
            enddo
        enddo
        
        allocate(temp(EmbSize,nSites))
        !Core hamiltonian
        call ZGEMM('C','N',EmbSize,nSites,nSites,zone,EmbeddedBasis_c,nSites,h0_c,nSites,zzero,temp,EmbSize)
        call ZGEMM('N','N',EmbSize,EmbSize,nSites,zone,temp,EmbSize,EmbeddedBasis_c,nSites,zzero,Emb_h0_c,EmbSize)

        !Embedded hamiltonian. Just h0 over impurity and its coupling to bath.
        call ZGEMM('C','N',EmbSize,nSites,nSites,zone,EmbeddedBasis_c,nSites,matc,nSites,zzero,temp,EmbSize)
        call ZGEMM('N','N',EmbSize,EmbSize,nSites,zone,temp,EmbSize,EmbeddedBasis_c,nSites,zzero,Emb_h0v_c,EmbSize)

        !Remove fit hamiltonian from impurity space and couplings
        Emb_h0v_c(1:nImp,1:EmbSize) = Emb_h0_c(1:nImp,1:EmbSize)
        Emb_h0v_c(nImp+1:EmbSize,1:nImp) = Emb_h0_c(nImp+1:EmbSize,1:nImp)
        
        deallocate(orbs,energies,ProjOverlap,ProjOverlapEVals,RotOccOrbs,ImpurityOrbs)
        deallocate(ProjOverlapVirt,OverlapVirt,VirtSpace,HFtoSchmidtTransform_c)
        deallocate(h0_c,temp,EmbeddedBasis_c)

        if(tCheck) then
            !Is the full schmidt basis a unitary transformation?
            !Is the FockSchmidt_c a hermitian operator
            allocate(temp(nSites,nSites))
            call ZGEMM('C','N',nSites,nSites,nSites,zone,FullSchmidtBasis_c,nSites,FullSchmidtBasis_c,nSites,zzero,temp,nSites)

            do i = 1,nSites
                do j = i+1,nSites
                    if(abs(temp(i,j)).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Schmidt transform not unitary')
                    elseif(abs(temp(j,i)).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Schmidt transform not unitary 2')
                    endif
                    if(abs(FockSchmidt_c(i,j)-dconjg(FockSchmidt_c(j,i))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'FockSchmidt_c not hermitian')
                    endif
                enddo
                if(abs(aimag(FockSchmidt_c(i,i))).gt.1.0e-8_dp) then
                    call stop_all(t_r,'FockSchmidt_c not diagonal real')
                endif
                if(abs(temp(i,i)-zone).gt.1.0e-8_dp) then
                    call stop_all(t_r,'Schmidt transform not unitary')
                endif
            enddo
            deallocate(temp)
        endif

    end subroutine SchmidtDecompose_C

    !Construct full embedding basis, along with orthogonal core and virtual set, and check orthonormality of orbital space
    subroutine ConstructFullSchmidtBasis()
        implicit none
        real(dp), allocatable :: ProjOverlap(:,:),ProjOverlapEVals(:),Work(:)
        real(dp), allocatable :: RotOccOrbs(:,:),ImpurityOrbs(:,:),ProjOverlapVirt(:,:)
        real(dp), allocatable :: OverlapVirt(:,:),VirtSpace(:,:),temp(:,:) 
        real(dp) :: norm,DDOT,Overlap
        integer :: lwork,info,i,j,nVirt
        character(len=*), parameter :: t_r='ConstructFullSchmidtBasis'
        integer :: nbath

        write(6,*) "Constructing full schmidt basis of HF det, including core + virt spaces"

        !HFOrbs defines HF basis.
        !Construct full projected overlap basis and diagonalize to see redundancies
        !Just do this for the occupied orbitals, and check orthogonality to the core HF and virtual orbitals
        allocate(ProjOverlap(nOcc,nOcc))
        call DGEMM('T','N',nOcc,nOcc,nImp,1.0_dp,HFOrbs(1:nImp,1:nOcc),nImp,HFOrbs(1:nImp,1:nOcc),nImp,0.0_dp,ProjOverlap,nOcc)

        !Diagonalize this
        allocate(ProjOverlapEVals(nOcc))
        ProjOverlapEVals(:)=0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        if(tWriteOut) then
            call writevector(ProjOverlapEVals,'Projected overlap eigenvalues alpha')
        endif

        !We should only have nImp non-zero eigenvalues
        nbath = 0
        do i=1,nOcc
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.gt.nImp) call stop_all(t_r,'error here')

        !Now rotate original occupied orbitals into this new orthogonal basis
        allocate(RotOccOrbs(nSites,nOcc))
        call DGEMM('N','N',nSites,nOcc,nOcc,1.0_dp,HFOrbs(1:nSites,1:nOcc),nSites,ProjOverlap,nOcc,0.0_dp,RotOccOrbs,nSites)

        !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
        !Only the last nImp orbitals will have any coupling to the impurity on them.
        !These RotOccOrbs constitute a legitamate HF wavefunction, are orthonormal to all other orbitals. Just simple rotation.
!        call writematrix(RotOccOrbs,'Occupied Orbs in schmidt basis',.true.)

        !Construct bath states by projecting out component on impurity and renormalizing
        !Assume last nImp states are the states with overlap with impurity only
        !Also normalize these orbitals
        do i=nOcc,nOcc-nImp+1,-1
            RotOccOrbs(1:nImp,i) = 0.0_dp
            norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
            norm = sqrt(norm)
            RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
        enddo

        !These states are now the bath states.
!        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

        allocate(ImpurityOrbs(nSites,nImp))
        ImpurityOrbs(:,:) = 0.0_dp
        do i=1,nImp
            ImpurityOrbs(i,i) = 1.0_dp
        enddo

        !We now have all the orbitals. Which are orthogonal to which?
!        write(6,*) "All alpha impurity/bath orbitals orthogonal by construction"
        do i=nOcc,nOcc-nImp+1,-1
            do j=nOcc,nOcc-nImp+1,-1
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(i.eq.j) then
                    if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not normalized set')
                    endif
                else
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not orthogonal set')
                    endif
                endif
            enddo
        enddo
!        write(6,*) "All alpha bath/bath orbitals orthonormal"

        !Now, consider the rotated core orbitals
        !Are they orthonormal to the bath orbitals (they are othogonal to impurity, since they have no component on them)
        !They must also be orthogonal to the bath orbitals, since they never had any component on the impurity, 
        !which is the only bit which has changed in constructing the bath (norm doesn't affect)
        do i=1,nOcc-nImp
            do j=nOcc,nOcc-nImp+1,-1
                !Overlap of core (i) with bath (j)
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'bath orbitals with core not orthogonal set')
                endif
            enddo
        enddo

        !However, the virtual space is *not* orthogonal to the embedded system (though it is wrt the core).
        !Now create orthogonal set of orbitals from the virtual space. There will now be a redundancy.
        !Calculate the overlap of the virtual space with a projection onto the embedded system.
        !From the diagonalization of this, we expect exactly *nImp* non-zero eigenvalues, which are the redundant orbitals.
        !Remove these, and the rest are the now non-canonical virtual orbital space.
        nVirt = nSites - nOcc
        allocate(ProjOverlapVirt(nVirt,nVirt))

        !This array is used to calculate the overlap of each virtual space function with each impurity
        !function
        allocate(OverlapVirt(2*nImp,nVirt))
        !The first nImp correspond to the impurity orbital, and the next two correspond to the bath orbitals
        OverlapVirt(:,:) = 0.0_dp
        do i=nOcc+1,nSites  !run through virtual space
            do j=1,nImp     !run through impurity orbitals
                OverlapVirt(j,i-nOcc) = HFOrbs(j,i)
            enddo
        enddo

        !Now calculate overlap with bath orbitals
        call DGEMM('T','N',nImp,nVirt,nSites,1.0_dp,RotOccOrbs(:,nOcc-nImp+1:nOcc),nSites,HFOrbs(:,nOcc+1:nSites), &
            nSites,0.0_dp,OverlapVirt(nImp+1:2*nImp,1:nVirt),2*nImp-nImp)

        !Combine overlaps to get full projected overlap matrix
        call DGEMM('T','N',nVirt,nVirt,2*nImp,1.0_dp,OverlapVirt,2*nImp,OverlapVirt,2*nImp,0.0_dp,ProjOverlapVirt,nVirt)

        !Diagonalize this
        deallocate(ProjOverlapEVals)
        allocate(ProjOverlapEVals(nVirt))
        ProjOverlapEVals(:)=0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

!        call writevector(ProjOverlapEVals,'virtual overlap eigenvalues')

        nbath = 0   !Count the number of virtual functions which span the space of the embedded system
        do i=1,nVirt
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.ne.nImp) then
            call stop_all(t_r,'Virtual space redundancy not as expected')
        endif

        !Now rotate orbitals such that they are orthogonal, while deleting the redundant ones
        !Assume the last nImp are redundant
        allocate(VirtSpace(nSites,nVirt-nImp))
        call DGEMM('N','N',nSites,nVirt-nImp,nVirt,1.0_dp,HFOrbs(:,nOcc+1:nSites),nSites,   &
            ProjOverlapVirt(1:nVirt,1:nVirt-nImp),nVirt,0.0_dp,VirtSpace,nSites)

        !Check now that the virtual space is now orthogonal to all occupied HF orbitals, as well as the impurity and bath sites

        !First against the impurity sites
        do i=1,nImp
            do j=1,nVirt-nImp
                Overlap = DDOT(nSites,ImpurityOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'virtual orbitals not orthogonal to impurity orbitals')
                endif
            enddo
        enddo

        do i=1,nOcc
            do j=1,nVirt-nImp
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    if(i.gt.(nOcc-nImp)) then
                        call stop_all(t_r,'virtual orbitals not orthogonal to bath orbitals')
                    else
                        call stop_all(t_r,'virtual orbitals not orthogonal to core orbitals')
                    endif
                endif
            enddo
        enddo

        if(allocated(FullSchmidtBasis)) deallocate(FullSchmidtBasis)
        allocate(FullSchmidtBasis(nSites,nSites))   ! (Atomicbasis,SchmidtBasis)
        FullSchmidtBasis(:,:) = 0.0_dp
        FullSchmidtBasis(:,1:nOcc-nImp) = RotOccOrbs(:,1:nOcc-nImp)    !The first orbitals consist of core  
        FullSchmidtBasis(:,nOcc-nImp+1:nOcc) = ImpurityOrbs(:,:)    !Impurity Orbs
        FullSchmidtBasis(:,nOcc+1:nOcc+nImp) = RotOccOrbs(:,nOcc-nImp+1:nOcc)   !Bath
        FullSchmidtBasis(:,nOcc+nImp+1:nSites) = VirtSpace(:,:)     !Virtual space

        !call writematrix(FullSchmidtBasis,'FullSchmidtBasis',.true.)

        !Construct unitary basis transformation matrix from HF to embedded basis
        if(allocated(HFtoSchmidtTransform)) deallocate(HFtoSchmidtTransform)
        allocate(HFtoSchmidtTransform(nSites,nSites))
        call DGEMM('T','N',nSites,nSites,nSites,1.0_dp,HFOrbs,nSites,FullSchmidtBasis,nSites,0.0_dp,HFtoSchmidtTransform,nSites)

        allocate(temp(nSites,nSites))
        !Check that this operator is unitary
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

        !Use this to rotate the fock operator into this new basis
        if(allocated(FockSchmidt)) deallocate(FockSchmidt)
        allocate(FockSchmidt(nSites,nSites))
        !Set up FockSchmidt to temperarily to be the HF basis fock operator (i.e. diagonal)
        FockSchmidt(:,:) = zero
        do i=1,nSites
            FockSchmidt(i,i) = HFEnergies(i)
        enddo
        call DGEMM('T','N',nSites,nSites,nSites,one,HFtoSchmidtTransform,nSites,FockSchmidt,nSites,zero,temp,nSites)
        call DGEMM('N','N',nSites,nSites,nSites,one,temp,nSites,HFtoSchmidtTransform,nSites,zero,FockSchmidt,nSites)
            
        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,2*nImp))
        EmbeddedBasis(:,:) = 0.0_dp
        EmbeddedBasis(:,1:nImp) = ImpurityOrbs(:,:)
        EmbeddedBasis(:,nImp+1:2*nImp) = FullSchmidtBasis(:,nOcc+1:nOcc+nImp)

        if(tUHF) then
            !Do all the same for the beta orbital space
            write(6,"(A)") "Constructing beta bath space"
            !HFOrbs defines HF basis.
            !Construct full projected overlap basis and diagonalize to see redundancies
            !Just do this for the occupied orbitals, and check orthogonality to the core HF and virtual orbitals
            call DGEMM('T','N',nOcc,nOcc,nImp,one,HFOrbs_b(1:nImp,1:nOcc),nImp,HFOrbs_b(1:nImp,1:nOcc),nImp,zero,ProjOverlap,nOcc)

            !Diagonalize this
            deallocate(ProjOverlapEVals)
            allocate(ProjOverlapEVals(nOcc))
            ProjOverlapEVals(:) = zero
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            if(tWriteOut) then
                call writevector(ProjOverlapEVals,'Projected overlap eigenvalues beta')
            endif

            !We should only have nImp non-zero eigenvalues
            nbath = 0
            do i=1,nOcc
                if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                    nbath = nbath + 1
                endif
            enddo
            if(nbath.gt.nImp) call stop_all(t_r,'error here')

            !Now rotate original beta occupied orbitals into this new orthogonal basis
            call DGEMM('N','N',nSites,nOcc,nOcc,one,HFOrbs_b(1:nSites,1:nOcc),nSites,ProjOverlap,nOcc,zero,RotOccOrbs,nSites)

            !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
            !Only the last nImp orbitals will have any coupling to the impurity on them.
            !These RotOccOrbs constitute a legitamate HF wavefunction, are orthonormal to all other orbitals. Just simple rotation.
    !        call writematrix(RotOccOrbs,'Occupied Orbs in schmidt basis',.true.)

            !Construct bath states by projecting out component on impurity and renormalizing
            !Assume last nImp states are the states with overlap with impurity only
            !Also normalize these orbitals
            do i=nOcc,nOcc-nImp+1,-1
                RotOccOrbs(1:nImp,i) = zero
                norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
                norm = sqrt(norm)
                RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
            enddo

            !These states are now the bath states.
    !        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

            ImpurityOrbs(:,:) = zero 
            do i=1,nImp
                ImpurityOrbs(i,i) = one 
            enddo

            !We now have all the orbitals. Which are orthogonal to which?
            write(6,*) "All beta impurity/bath orbitals orthogonal by construction"
            do i=nOcc,nOcc-nImp+1,-1
                do j=nOcc,nOcc-nImp+1,-1
                    Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                    if(i.eq.j) then
                        if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                            call stop_all(t_r,'beta bath orbitals not normalized set')
                        endif
                    else
                        if(abs(Overlap).gt.1.0e-7_dp) then
                            call stop_all(t_r,'beta bath orbitals not orthogonal set')
                        endif
                    endif
                enddo
            enddo
            write(6,*) "All beta bath/bath orbitals orthonormal"

            !Now, consider the rotated core orbitals
            !Are they orthonormal to the bath orbitals (they are othogonal to impurity, since they have no component on them)
            !They must also be orthogonal to the bath orbitals, since they never had any component on the impurity, 
            !which is the only bit which has changed in constructing the bath (norm doesn't affect)
            do i=1,nOcc-nImp
                do j=nOcc,nOcc-nImp+1,-1
                    !Overlap of core (i) with bath (j)
                    Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'beta bath orbitals with core not orthogonal set')
                    endif
                enddo
            enddo

            !However, the virtual space is *not* orthogonal to the embedded system (though it is wrt the core).
            !Now create orthogonal set of orbitals from the virtual space. There will now be a redundancy.
            !Calculate the overlap of the virtual space with a projection onto the embedded system.
            !From the diagonalization of this, we expect exactly *nImp* non-zero eigenvalues, which are the redundant orbitals.
            !Remove these, and the rest are the now non-canonical virtual orbital space.
            nVirt = nSites - nOcc
            !allocate(ProjOverlapVirt(nVirt,nVirt))

            !This array is used to calculate the overlap of each virtual space function with each impurity
            !function
            !allocate(OverlapVirt(2*nImp,nVirt))
            !The first nImp correspond to the impurity orbital, and the next two correspond to the bath orbitals
            OverlapVirt(:,:) = zero
            do i=nOcc+1,nSites  !run through virtual space
                do j=1,nImp     !run through impurity orbitals
                    OverlapVirt(j,i-nOcc) = HFOrbs_b(j,i)
                enddo
            enddo

            !Now calculate overlap with bath orbitals
            call DGEMM('T','N',nImp,nVirt,nSites,one,RotOccOrbs(:,nOcc-nImp+1:nOcc),nSites,HFOrbs_b(:,nOcc+1:nSites), &
                nSites,zero,OverlapVirt(nImp+1:2*nImp,1:nVirt),2*nImp-nImp)

            !Combine overlaps to get full projected overlap matrix
            call DGEMM('T','N',nVirt,nVirt,2*nImp,one,OverlapVirt,2*nImp,OverlapVirt,2*nImp,zero,ProjOverlapVirt,nVirt)

            !Diagonalize this
            deallocate(ProjOverlapEVals)
            allocate(ProjOverlapEVals(nVirt))
            ProjOverlapEVals(:) = zero
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

    !        call writevector(ProjOverlapEVals,'virtual overlap eigenvalues')

            nbath = 0   !Count the number of virtual functions which span the space of the embedded system
            do i=1,nVirt
                if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                    nbath = nbath + 1
                endif
            enddo
            if(nbath.ne.nImp) then
                call stop_all(t_r,'Virtual beta space redundancy not as expected')
            endif

            !Now rotate orbitals such that they are orthogonal, while deleting the redundant ones
            !Assume the last nImp are redundant
            !allocate(VirtSpace(nSites,nVirt-nImp))
            call DGEMM('N','N',nSites,nVirt-nImp,nVirt,one,HFOrbs_b(:,nOcc+1:nSites),nSites,   &
                ProjOverlapVirt(1:nVirt,1:nVirt-nImp),nVirt,zero,VirtSpace,nSites)

            !Check now that the virtual space is now orthogonal to all occupied HF orbitals, as well as the impurity and bath sites

            !First against the impurity sites
            do i=1,nImp
                do j=1,nVirt-nImp
                    Overlap = DDOT(nSites,ImpurityOrbs(:,i),1,VirtSpace(:,j),1)
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'virtual beta orbitals not orthogonal to impurity orbitals')
                    endif
                enddo
            enddo

            do i=1,nOcc
                do j=1,nVirt-nImp
                    Overlap = DDOT(nSites,RotOccOrbs(:,i),1,VirtSpace(:,j),1)
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        if(i.gt.(nOcc-nImp)) then
                            call stop_all(t_r,'virtual beta orbitals not orthogonal to bath orbitals')
                        else
                            call stop_all(t_r,'virtual beta orbitals not orthogonal to core orbitals')
                        endif
                    endif
                enddo
            enddo

            if(allocated(FullSchmidtBasis_b)) deallocate(FullSchmidtBasis_b)
            allocate(FullSchmidtBasis_b(nSites,nSites))   ! (Atomicbasis,SchmidtBasis)
            FullSchmidtBasis_b(:,:) = zero 
            FullSchmidtBasis_b(:,1:nOcc-nImp) = RotOccOrbs(:,1:nOcc-nImp)    !The first orbitals consist of core  
            FullSchmidtBasis_b(:,nOcc-nImp+1:nOcc) = ImpurityOrbs(:,:)    !Impurity Orbs
            FullSchmidtBasis_b(:,nOcc+1:nOcc+nImp) = RotOccOrbs(:,nOcc-nImp+1:nOcc)   !Bath
            FullSchmidtBasis_b(:,nOcc+nImp+1:nSites) = VirtSpace(:,:)     !Virtual space

            !call writematrix(FullSchmidtBasis,'FullSchmidtBasis',.true.)

            !Construct unitary basis transformation matrix from HF to embedded basis
            if(allocated(HFtoSchmidtTransform_b)) deallocate(HFtoSchmidtTransform_b)
            allocate(HFtoSchmidtTransform_b(nSites,nSites))
            call DGEMM('T','N',nSites,nSites,nSites,one,HFOrbs_b,nSites,FullSchmidtBasis_b, &
                nSites,zero,HFtoSchmidtTransform_b,nSites)

            !Check that this operator is unitary
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

            !Use this to rotate the fock operator into this new basis
            if(allocated(FockSchmidt_b)) deallocate(FockSchmidt_b)
            allocate(FockSchmidt_b(nSites,nSites))
            !Set up FockSchmidt to temperarily to be the HF basis fock operator (i.e. diagonal)
            FockSchmidt_b(:,:) = zero 
            do i=1,nSites
                FockSchmidt_b(i,i) = HFEnergies_b(i)
            enddo
            call DGEMM('T','N',nSites,nSites,nSites,one,HFtoSchmidtTransform_b,nSites,FockSchmidt_b,nSites,zero,temp,nSites)
            call DGEMM('N','N',nSites,nSites,nSites,one,temp,nSites,HFtoSchmidtTransform_b,nSites,zero,FockSchmidt_b,nSites)
        
            if(allocated(EmbeddedBasis_b)) deallocate(EmbeddedBasis_b)
            allocate(EmbeddedBasis_b(nSites,2*nImp))
            EmbeddedBasis_b(:,:) = zero
            EmbeddedBasis_b(:,1:nImp) = ImpurityOrbs(:,:)
            EmbeddedBasis_b(:,nImp+1:2*nImp) = FullSchmidtBasis_b(:,nOcc+1:nOcc+nImp)
        endif
        deallocate(temp)

        !do i=1,nSites
        !    write(6,*) "FOCKSCHMIDT: ",i,FockSchmidt(i,i)
        !enddo
        !FockSchmidt has been overwritten with the fock matrix in the schmidt basis
!        call writematrix(FockSchmidt,'Fock in schmidt basis',.true.)

!       Calculate the non-interacting core energy of the DMET wavefunction
        CoreEnergy = zero
        do i = 1,nOcc-nImp
            CoreEnergy = CoreEnergy + FockSchmidt(i,i)
            if(tUHF) CoreEnergy = CoreEnergy + FockSchmidt_b(i,i)
        enddo
        if(.not.tUHF) CoreEnergy = CoreEnergy * 2.0_dp
        write(6,*) "Non-interacting core energy for DMET wavefunction is: ",CoreEnergy

        EmbSize = 2*nImp      !This is the total size of the embedded system with which to do the high-level calculation on 
        if(tUHF) then
            EmbSizeSpin = 2*EmbSize
        else
            EmbSizeSpin = EmbSize
        endif
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(ProjOverlap,ProjOverlapEVals,RotOccOrbs,ImpurityOrbs,ProjOverlapVirt,OverlapVirt,VirtSpace)
    end subroutine ConstructFullSchmidtBasis


    !Transform the one-electron integrals into the embedded basis by unitary transformation
    subroutine Transform1e()
        implicit none
        real(dp), allocatable :: CorrPotential(:,:),temp(:,:),CorrPotential_b(:,:)
!        real(dp), allocatable :: FockPotential(:,:)
        integer :: i,j
!        character(len=*), parameter :: t_r="Transform1e"
        
        if(allocated(Emb_Fock)) deallocate(Emb_Fock)
        if(allocated(Emb_h0)) deallocate(Emb_h0)                 !Core hamiltonian
        if(allocated(Emb_MF_DM)) deallocate(Emb_MF_DM)        !Mean field RDM
        if(allocated(Emb_FockPot)) deallocate(Emb_FockPot)      !The fock potential transforms h0v into the fock matrix in the AO basis
        if(allocated(Emb_CorrPot)) deallocate(Emb_CorrPot)      !The correlation potential is the block diagonal v_pot in the AO basis
        !Emb_h0v is the core hamiltonian in the embedded basis, with correlation potential (not over the impurity sites)
        if(allocated(Emb_h0v)) deallocate(Emb_h0v)

        !Transform all of these quantities into the embedded basis
        allocate(Emb_h0(EmbSize,EmbSize))
        allocate(Emb_MF_DM(EmbSize,EmbSize))
        allocate(Emb_FockPot(EmbSize,EmbSize))
        allocate(Emb_CorrPot(EmbSize,EmbSize))
        allocate(Emb_Fock(EmbSize,EmbSize))     !For hub, this is just the core + corrPot
        allocate(Emb_h0v(EmbSize,EmbSize))

        !Rotate them
        allocate(temp(EmbSize,nSites))
        !Core hamiltonian
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_h0,EmbSize)

        if(tDebug) call writematrix(Emb_h0,'Embedded Core hamil',.true.)

        !Mean field RDM
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,MeanFieldDM,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_MF_DM,EmbSize)
        
        if(tDebug) call writematrix(Emb_MF_DM,'Embedded RDM',.true.)

        !Since we don't actually store the fock potential (the diagonal repulsion/2 electron terms), construct it here in the AO basis
        !This is the bit of the fock matrix not in the core hamiltonian
!        allocate(FockPotential(nSites,nSites))
!        FockPotential(:,:) = 0.0_dp
!        do i=1,nSites
!            !Include the on-site repulsion
!            FockPotential(i,i) = FockPotential(i,i) + U * 0.5_dp * (NEl/real(nSites))
!        enddo
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,FockPotential,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_FockPot,EmbSize)
!        deallocate(FockPotential)
        Emb_FockPot(:,:)=0.0_dp     !We don't include the 2e terms in the fock matrix here, instead, leaving them to be captured by the correlation potential
!        call writematrix(Emb_FockPot,"Embedded Fock potential",.true.)

        !We also do not store the "Correlation potential", which is the potential which is added to the fock matrix to make the DMs match
        allocate(CorrPotential(nSites,nSites))
        CorrPotential(:,:) = 0.0_dp
        !It is just h0v - h0
        do i=1,nSites
            do j=1,nSites
                CorrPotential(i,j) = h0v(i,j) - h0(i,j)
            enddo
        enddo
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,CorrPotential,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_CorrPot,EmbSize)
        deallocate(CorrPotential)
        deallocate(temp)

        if(tDebug) call writematrix(Emb_CorrPot,"Embedded correlation potential",.true.)

        !The 2e terms left out of fock matrix everywhere, although included basically in CorrPot
        Emb_Fock(:,:) = Emb_h0(:,:) + Emb_CorrPot(:,:) ! + Emb_FockPot  
        
        !Set the correlation potential to zero over the impurity, since we do not need to include the fitted potential over the 
        !orbitals that we are going to do the high-level calculation on.
        Emb_h0v(:,:) = Emb_CorrPot(:,:)
        Emb_h0v(1:nImp,1:nImp) = 0.0_dp     !Set correlation potential over the impurity sites to zero
        Emb_h0v(:,:) = Emb_h0v(:,:) + Emb_h0(:,:)    !Add the embedded core hamiltonian over all sites

        if(tUHF) then
            !Now for the beta orbitals!!
            if(allocated(Emb_Fock_b)) deallocate(Emb_Fock_b)
            if(allocated(Emb_h0_b)) deallocate(Emb_h0_b)                 !Core hamiltonian
            if(allocated(Emb_MF_DM_b)) deallocate(Emb_MF_DM_b)        !Mean field RDM
            if(allocated(Emb_FockPot_b)) deallocate(Emb_FockPot_b)      !The fock potential transforms h0v into the fock matrix in the AO basis
            if(allocated(Emb_CorrPot_b)) deallocate(Emb_CorrPot_b)      !The correlation potential is the block diagonal v_pot in the AO basis
            !Emb_h0v is the core hamiltonian in the embedded basis, with correlation potential (not over the impurity sites)
            if(allocated(Emb_h0v_b)) deallocate(Emb_h0v_b)

            !Transform all of these quantities into the embedded basis
            allocate(Emb_h0_b(EmbSize,EmbSize))
            allocate(Emb_MF_DM_b(EmbSize,EmbSize))
            allocate(Emb_FockPot_b(EmbSize,EmbSize))
            allocate(Emb_CorrPot_b(EmbSize,EmbSize))
            allocate(Emb_Fock_b(EmbSize,EmbSize))     !For hub, this is just the core + corrPot
            allocate(Emb_h0v_b(EmbSize,EmbSize))

            !Rotate them
            allocate(temp(EmbSize,nSites))
            !Core hamiltonian
            call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis_b,nSites,h0_b,nSites,0.0_dp,temp,EmbSize)
            call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis_b,nSites,0.0_dp,Emb_h0_b,EmbSize)

            if(tDebug) call writematrix(Emb_h0_b,'Embedded Core beta hamil',.true.)

            !Mean field RDM
            call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis_b,nSites,MeanFieldDM_b,nSites,0.0_dp,temp,EmbSize)
            call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis_b,nSites,0.0_dp,Emb_MF_DM_b,EmbSize)
            
            if(tDebug) call writematrix(Emb_MF_DM_b,'Embedded beta RDM',.true.)

            !Since we don't actually store the fock potential (the diagonal repulsion/2 electron terms), construct it here in the AO basis
            !This is the bit of the fock matrix not in the core hamiltonian
    !        allocate(FockPotential(nSites,nSites))
    !        FockPotential(:,:) = 0.0_dp
    !        do i=1,nSites
    !            !Include the on-site repulsion
    !            FockPotential(i,i) = FockPotential(i,i) + U * 0.5_dp * (NEl/real(nSites))
    !        enddo
    !        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,FockPotential,nSites,0.0_dp,temp,EmbSize)
    !        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_FockPot,EmbSize)
    !        deallocate(FockPotential)
            Emb_FockPot_b(:,:)=0.0_dp     !We don't include the 2e terms in the fock matrix here, instead, leaving them to be captured by the correlation potential
    !        call writematrix(Emb_FockPot,"Embedded Fock potential",.true.)

            !We also do not store the "Correlation potential", which is the potential which is added to the fock matrix to make the DMs match
            allocate(CorrPotential_b(nSites,nSites))
            CorrPotential_b(:,:) = 0.0_dp
            !It is just h0v - h0
            do i=1,nSites
                do j=1,nSites
                    CorrPotential_b(i,j) = h0v_b(i,j) - h0_b(i,j)
                enddo
            enddo
            call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis_b,nSites,CorrPotential_b,nSites,0.0_dp,temp,EmbSize)
            call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis_b,nSites,0.0_dp,Emb_CorrPot_b,EmbSize)
            deallocate(CorrPotential_b)
            deallocate(temp)

            if(tDebug) call writematrix(Emb_CorrPot_b,"Embedded beta correlation potential",.true.)

            !The 2e terms left out of fock matrix everywhere, although included basically in CorrPot
            Emb_Fock_b(:,:) = Emb_h0_b(:,:) + Emb_CorrPot_b(:,:) ! + Emb_FockPot  
            
            !Set the correlation potential to zero over the impurity, since we do not need to include the fitted potential over the 
            !orbitals that we are going to do the high-level calculation on.
            Emb_h0v_b(:,:) = Emb_CorrPot_b(:,:)
            Emb_h0v_b(1:nImp,1:nImp) = 0.0_dp     !Set correlation potential over the impurity sites to zero
            Emb_h0v_b(:,:) = Emb_h0v_b(:,:) + Emb_h0_b(:,:)    !Add the embedded core hamiltonian over all sites

        endif
        
    end subroutine Transform1e

    !Calculate the embedded basis through decomposition of the 1RDM, and transform the 1e operators into this new basis.
    !Core and virtual basis not included.
    subroutine CalcEmbedding()
        implicit none
        real(dp) :: ImpurityOverlap(nImp,nImp),OverlapEVs(nImp)
        real(dp), allocatable :: RDMonImp(:,:),Work(:),SminHalf(:,:),temp(:,:)
        integer :: info,lWork,i,j,nDelete,nSys
        character(len=*), parameter :: t_r="CalcEmbedding"
            
        !Take density matrix over impurity: nSites-nImp x nImp
        allocate(RDMonImp(nSites-nImp,nImp))
        do i=nImp+1,nSites
            do j=1,nImp
                RDMonImp(i-nImp,j) = MeanFieldDM(i,j)
            enddo
        enddo

!        call writematrix(RDMonImp,'RDMonImp',.true.)

        !Now, we want to Lowdin orthogonalize these orbitals, i.e. |psi_alpha> = |psi_beta> S_beta,alpha^(-1/2)

        !Find overlap of this density over the impurity (contract out nSites-nImp)
        call DGEMM('T','N',nImp,nImp,nSites-nImp,1.0_dp,RDMonImp,nSites-nImp,RDMonImp,nSites-nImp,0.0_dp,ImpurityOverlap,nImp)
!        call writematrix(ImpurityOverlap,'s',.true.)
        !Now find the eigenbasis of this overlap in order to raise it to the power of -1/2
        !Diagonalize the system over the impurity sites
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        if(nImp.eq.1) then
            if(abs(ImpurityOverlap(1,1)-1.0_dp).gt.1.0e-7_dp) call stop_all(t_r,'Diag error')
        endif
        !call writevector(OverlapEVs,'S_EVs')

        nSys = 0    !nSys is the number of orbitals in the bath having removed linear dependencies
        nDelete = 0
        do i=1,nImp
            if(OverlapEVs(i).lt.1.0e-10_dp) then
                write(6,*) "Warning: Bath basis linearly dependent"
                write(6,*) "Overlap eigenvalue: ",OverlapEVs(i)
                nDelete = nDelete + 1
            else
                nSys = nSys + 1
            endif
        enddo

        if(nSys.ne.nImp) write(6,*) "Bath orbitals removed due to linear dependency: ",nDelete
        write(6,*) "Total bath orbitals: ",nSys
        if(nDelete.ne.(nImp-nSys)) call stop_all(t_r,'nDelete.ne.(nImp-nSys)')

        !Now, we need to construct the full embedded basis including orbitals on the impurity sites
        !The orbitals to remove if necessary will be first.
        !Construct S^(-1/2), and then rotate back into the original basis
        allocate(SminHalf(nSys,nSys))
        SminHalf = 0.0_dp
        do i=1,nSys
            SminHalf(i,i) = OverlapEVs(i+nDelete)**(-0.5_dp)
        enddo
!        call writematrix(ImpurityOverlap(1:nImp,nDelete+1:nImp),'ImpOverlap',.true.)
        !Now rotate back into original basis as U D U^T (Having removed the first nDelete orbitals from the Impurity overlap eigenvectors
        allocate(temp(nImp,nSys))
        call DGEMM('N','N',nImp,nSys,nSys,1.0_dp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,SminHalf,nSys,0.0_dp,temp,nImp)
        call DGEMM('N','T',nImp,nImp,nSys,1.0_dp,temp,nImp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,0.0_dp,SminHalf,nImp)
        deallocate(temp)
        !SminHalf in now in the original basis
!        call writematrix(SminHalf,'SminHalf',.true.)

        !Now, we need to multiply the original bath orbitals extracted from the RDM on the impurity by this matrix to complete the Lowdin orthogonalization
        allocate(temp(nSites-nImp,nImp))
        call DGEMM('N','N',nSites-nImp,nImp,nImp,1.0_dp,RDMonImp,nSites-nImp,SminHalf,nImp,0.0_dp,temp,nSites-nImp)
        !temp now holds the orthogonalized bath orbitals
!        call writematrix(temp,'Orthog bath',.true.)

        !Now include the impurity sites along with the bath orbitals to make the final embedded basis
        !The system orbitals are represented by the first nImp orbitals, and the bath orbitals are the last nSys
        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,nImp+nSys))
        EmbeddedBasis(:,:) = 0.0_dp
        do i=1,nImp
            EmbeddedBasis(i,i) = 1.0_dp
        enddo

        if(nImp.ne.nSys) call stop_all(t_r,'You need to check in the code here to make sure we have the right bath orbitals')
        !Here, we only want to add orbitals on the bath which haven't been removed due to being linearly dependent
        !Which orbitals in temp are the removed linear dependent ones? TODO: Check this is being done correctly
        !Currently, we are just taking the first nSys of them
        EmbeddedBasis(nImp+1:nSites,nImp+1:2*nImp) = temp(:,1:nSys)

        if(tUHF) then
            !Now for beta orbitals
            !Take density matrix over impurity: nSites-nImp x nImp
            RDMonImp = zero
            do i=nImp+1,nSites
                do j=1,nImp
                    RDMonImp(i-nImp,j) = MeanFieldDM_b(i,j)
                enddo
            enddo

    !        call writematrix(RDMonImp,'RDMonImp',.true.)

            !Now, we want to Lowdin orthogonalize these orbitals, i.e. |psi_alpha> = |psi_beta> S_beta,alpha^(-1/2)

            !Find overlap of this density over the impurity (contract out nSites-nImp)
            call DGEMM('T','N',nImp,nImp,nSites-nImp,one,RDMonImp,nSites-nImp,RDMonImp,nSites-nImp,zero,ImpurityOverlap,nImp)
    !        call writematrix(ImpurityOverlap,'s',.true.)
            !Now find the eigenbasis of this overlap in order to raise it to the power of -1/2
            !Diagonalize the system over the impurity sites
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
            if(nImp.eq.1) then
                if(abs(ImpurityOverlap(1,1)-1.0_dp).gt.1.0e-7_dp) call stop_all(t_r,'Diag error')
            endif
            !call writevector(OverlapEVs,'S_EVs')

            nSys = 0    !nSys is the number of orbitals in the bath having removed linear dependencies
            nDelete = 0
            do i=1,nImp
                if(OverlapEVs(i).lt.1.0e-10_dp) then
                    write(6,*) "Warning: Bath basis linearly dependent"
                    write(6,*) "Overlap eigenvalue: ",OverlapEVs(i)
                    nDelete = nDelete + 1
                else
                    nSys = nSys + 1
                endif
            enddo

            if(nSys.ne.nImp) write(6,*) "Bath orbitals removed due to linear dependency: ",nDelete
            write(6,*) "Total beta bath orbitals: ",nSys
            if(nDelete.ne.(nImp-nSys)) call stop_all(t_r,'nDelete.ne.(nImp-nSys)')

            !Now, we need to construct the full embedded basis including orbitals on the impurity sites
            !The orbitals to remove if necessary will be first.
            !Construct S^(-1/2), and then rotate back into the original basis
            SminHalf = zero 
            do i=1,nSys
                SminHalf(i,i) = OverlapEVs(i+nDelete)**(-0.5_dp)
            enddo
    !        call writematrix(ImpurityOverlap(1:nImp,nDelete+1:nImp),'ImpOverlap',.true.)
            !Now rotate back into original basis as U D U^T (Having removed the first nDelete orbitals from the Impurity overlap eigenvectors
            deallocate(temp)
            allocate(temp(nImp,nSys))
            call DGEMM('N','N',nImp,nSys,nSys,1.0_dp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,SminHalf,nSys,0.0_dp,temp,nImp)
            call DGEMM('N','T',nImp,nImp,nSys,1.0_dp,temp,nImp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,0.0_dp,SminHalf,nImp)
            deallocate(temp)
            !SminHalf in now in the original basis
    !        call writematrix(SminHalf,'SminHalf',.true.)

            !Now, we need to multiply the original bath orbitals extracted from the RDM on the impurity by this matrix to complete the Lowdin orthogonalization
            allocate(temp(nSites-nImp,nImp))
            call DGEMM('N','N',nSites-nImp,nImp,nImp,1.0_dp,RDMonImp,nSites-nImp,SminHalf,nImp,0.0_dp,temp,nSites-nImp)
            !temp now holds the orthogonalized bath orbitals
    !        call writematrix(temp,'Orthog bath',.true.)

            !Now include the impurity sites along with the bath orbitals to make the final embedded basis
            !The system orbitals are represented by the first nImp orbitals, and the bath orbitals are the last nSys
            if(allocated(EmbeddedBasis_b)) deallocate(EmbeddedBasis_b)
            allocate(EmbeddedBasis_b(nSites,nImp+nSys))
            EmbeddedBasis_b(:,:) = zero 
            do i=1,nImp
                EmbeddedBasis_b(i,i) = one 
            enddo

            if(nImp.ne.nSys) call stop_all(t_r,'You need to check in the code here to make sure we have the right bath orbitals')
            !Here, we only want to add orbitals on the bath which haven't been removed due to being linearly dependent
            !Which orbitals in temp are the removed linear dependent ones? TODO: Check this is being done correctly
            !Currently, we are just taking the first nSys of them
            EmbeddedBasis_b(nImp+1:nSites,nImp+1:2*nImp) = temp(:,1:nSys)

        endif
        
        EmbSize = 2*nImp      !This is the total size of the embedded system with which to do the high-level calculation on 
        if(tUHF) then
            EmbSizeSpin = 2*EmbSize
        else
            EmbSizeSpin = EmbSize
        endif
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(temp,SminHalf)
        
    end subroutine CalcEmbedding

    !Construct thermal embedding basis. Core and virtual set not included initially
    subroutine ConstructThermalEmbedding()
        implicit none
        real(dp), allocatable :: ProjOverlap(:,:),ProjOverlapEVals(:),Work(:)
        real(dp), allocatable :: RotOccOrbs(:,:),ImpurityOrbs(:,:)
        real(dp) :: norm,DDOT,Overlap,ThermoPot
        integer :: lwork,info,i,j,imp,ispin,jspin,ni,nj
        character(len=*), parameter :: t_r='ConstructThermalEmbedding'
        integer :: nbath

        write(6,*) "Constructing schmidt basis of thermal zeroth order wavefunction"

        !First, calculate the thermodynamic potential
        !Assume that the chemical potential is zero
        ThermoPot = zero
        do i = 1,nSites
            ThermoPot = ThermoPot + log(one + exp(-HFEnergies(i)/Temperature))
        enddo
        ThermoPot = -ThermoPot * Temperature

        write(6,*) "Thermodynamic potential (assuming mu=0) : ",ThermoPot

        allocate(ProjOverlap(nSites,nSites))
        do i = 1,nSites
            do ispin = 1,2  !Run over spins of i
                do j = 1,nSites
                    do jspin = 1,2  !Run over spins of j    (Shouldn't actually make any difference
                        if(mod(ispin,2).ne.mod(jspin,2)) cycle
                        do ni = 0,1
                            do nj = 0,1
                                do imp = 1,nImp
                                    ProjOverlap(i,j) = ProjOverlap(i,j) + exp(-HFEnergies(i)*ni)*HFOrbs(imp,i)  &
                                        *HFOrbs(imp,j)*exp(-HFEnergies(j)*nj)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        ProjOverlap(:,:) = ProjOverlap(:,:) * exp(ThermoPot/Temperature)

        !Diagonalize this
        allocate(ProjOverlapEVals(nSites))
        ProjOverlapEVals(:)=zero
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,ProjOverlap,nSites,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,ProjOverlap,nSites,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        write(6,*) "Largest 2 x nImp Projected overlap eigenvalues:"
        do i=nSites,nSites-2*nImp,-1
            write(6,*) ProjOverlapEVals(i)
        enddo

        if(tWriteOut) then
            call writevector(ProjOverlapEVals,'Projected overlap eigenvalues alpha')
        endif

        !We should only have nImp non-zero eigenvalues
        nbath = 0
        do i=1,nSites
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.gt.nImp) call stop_all(t_r,'error here')

        !Now rotate original occupied orbitals into this new orthogonal basis
        allocate(RotOccOrbs(nSites,nSites))
        call DGEMM('N','N',nSites,nSites,nSites,one,HFOrbs,nSites,ProjOverlap,nSites,zero,RotOccOrbs,nSites)

        !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
        !Only the last nImp orbitals will have any coupling to the impurity on them.
        !These RotOccOrbs constitute a legitamate wavefunction, are orthonormal to all other orbitals. Just simple rotation.

        !Construct bath states by projecting out component on impurity and renormalizing
        !Assume last nImp states are the states with overlap with impurity only
        !Also normalize these orbitals
        do i=nSites,nSites-nImp+1,-1
            RotOccOrbs(1:nImp,i) = zero
            norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
            norm = sqrt(norm)
            RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
        enddo

        !These states are now the bath states.
!        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

        allocate(ImpurityOrbs(nSites,nImp))
        ImpurityOrbs(:,:) = zero 
        do i=1,nImp
            ImpurityOrbs(i,i) = one 
        enddo

        !We now have all the orbitals. Which are orthogonal to which?
        write(6,*) "All alpha impurity/bath orbitals orthogonal by construction"
        do i=nOcc,nOcc-nImp+1,-1
            do j=nOcc,nOcc-nImp+1,-1
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(i.eq.j) then
                    if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not normalized set')
                    endif
                else
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not orthogonal set')
                    endif
                endif
            enddo
        enddo
        write(6,*) "All alpha bath/bath orbitals orthonormal"

        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,2*nImp))
        EmbeddedBasis(:,:) = zero
        EmbeddedBasis(:,1:nImp) = ImpurityOrbs(:,:)
        EmbeddedBasis(:,nImp+1:2*nImp) = RotOccOrbs(:,nOcc+1:nOcc+nImp)

        EmbSize = 2*nImp      !This is the total size of the embedded system with which to do the high-level calculation on 
        EmbSizeSpin = EmbSize
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(ProjOverlap,ProjOverlapEVals,RotOccOrbs,ImpurityOrbs)
    end subroutine ConstructThermalEmbedding
            
end module SchmidtDecomp
