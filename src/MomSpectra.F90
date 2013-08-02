module MomSpectra
    use const
    use timing
    use globals
    use errors, only: stop_all,warning
    use mat_tools, only: WriteVector,WriteMatrix,WriteVectorInt,WriteMatrixComp,WriteVectorComp,znrm2
    implicit none

    contains

    subroutine MomGF_Ex()
        use utils, only: get_free_unit,append_ext_real,append_ext
        use DetBitOps, only: DecodeBitDet,SQOperator,CountBits
        use DetToolsData
        use DetTools, only : GetHElement_comp
        use DetTools, only: tospat,umatind,gendets
        use solvers, only: CreateIntMats
        implicit none
        complex(dp), allocatable :: NFCIHam(:,:),Np1FCIHam_alpha(:,:),Nm1FCIHam_beta(:,:)
        complex(dp), allocatable :: LinearSystem_p(:,:),LinearSystem_h(:,:),Psi1_p(:),Psi1_h(:)
        complex(dp), allocatable :: Cre_0(:),Ann_0(:),FockSchmidt_SE(:,:),FockSchmidt_SE_VV(:,:)
        complex(dp), allocatable :: FockSchmidt_SE_CC(:,:),FockSchmidt_SE_VX(:,:),FockSchmidt_SE_XV(:,:)
        complex(dp), allocatable :: FockSchmidt_SE_CX(:,:),FockSchmidt_SE_XC(:,:)
        complex(dp), allocatable :: Gc_a_F_ax_Bra(:),Gc_a_F_ax_Ket(:),Gc_b_F_ab(:),Gz_i_F_xi_Bra(:)
        complex(dp), allocatable :: Ga_i_F_xi_Bra(:),Ga_i_F_xi_Ket(:),Ga_i_F_ij(:)
        complex(dp), allocatable :: SchmidtkGF_Cre_Ket(:),SchmidtkGF_Ann_Ket(:)
        complex(dp), allocatable :: SchmidtkGF_Cre_Bra(:),SchmidtkGF_Ann_Bra(:)
        complex(dp), allocatable :: StatickGF_Cre_Ket(:),StatickGF_Ann_Ket(:)
        complex(dp), allocatable :: StatickGF_Cre_Bra(:),StatickGF_Ann_Bra(:)
        integer, allocatable :: Coup_Create_alpha(:,:,:),Coup_Ann_alpha(:,:,:)
        integer :: nLinearSystem,nLinearSystem_h,CoreEnd,VirtStart,VirtEnd,ActiveStart,ActiveEnd,nCore,nVirt
        integer :: VIndex,iunit,i,j,k,DiffOrb,nOrbs,orbdum(1),gam,tempK,KPnt,ierr
        complex(dp) :: ni_lr_ann,ni_lr_cre
        real(dp) :: SpectralWeight,Prev_Spec,mu,Omega
        logical :: tFirst,tParity,tFinishedK
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='MomGF_Ex'
        
        call set_timer(LR_EC_GF_Precom)
        
        SpectralWeight = zero 
        Prev_Spec = zero
        tFirst = .true.

        write(6,"(A)") "Calculating non-interacting EC MR-TDA LR system for *hole* & *particle* "   &
            &   //"alpha spin-orbital k-space perturbations..."
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(iSolveLR.eq.1) then
            write(6,"(A)") "Solving linear system with standard ZGESV linear solver"
        elseif(iSolveLR.eq.2) then
            write(6,"(A)") "Solving linear system with advanced ZGELS linear solver"
        elseif(iSolveLR.eq.3) then
            write(6,"(A)") "Solving linear system with direct inversion of hamiltonian"
        elseif(iSolveLR.eq.4) then
            write(6,"(A)") "Solving linear system via complete diagonalization of hamiltonian"
        else
            call stop_all(t_r,"Linear equation solver unknown")
        endif
        if(tLR_ReoptGS) call stop_all(t_r,'Cannot reopt GS for k-space GFs')
        if(tMinRes_NonDir.or.tGMRes_NonDir) call stop_all(t_r,'Cannot use iterative solvers for k-space GFs')

        !umat and tmat for the active space
        call CreateIntMats(tComp=.true.)
        
        !Enumerate excitations for fully coupled space
        !Seperate the lists into different Ms sectors in the N+- lists
        call GenDets(Elec,EmbSize,.true.,.true.,.true.) 
        write(6,*) "Number of determinants in {N,N+1,N-1} FCI space: ",ECoupledSpace
        write(6,"(A,F12.5,A)") "Memory required for det list storage: ",DetListStorage*RealtoMb, " Mb"

        !We are going to choose the uncontracted particle space to have an Ms of 1 (i.e. alpha spin)
        !Construct FCI hamiltonians for the N, N+1_alpha and N-1_beta spaces
        if(nNp1FCIDet.ne.nNm1FCIDet) call stop_all(t_r,'Active space not half-filled')
        if(nNp1FCIDet.ne.nNp1bFCIDet) call stop_all(t_r,'Cannot deal with open shell systems')
        write(6,"(A,F12.5,A)") "Memory required for N-electron hamil: ",(real(nFCIDet,dp)**2)*ComptoMb, " Mb"
        write(6,"(A,F12.5,A)") "Memory required for N+1-electron hamil: ",(real(nNp1FCIDet,dp)**2)*ComptoMb, " Mb"
        write(6,"(A,F12.5,A)") "Memory required for N-1-electron hamil: ",(real(nNm1bFCIDet,dp)**2)*ComptoMb, " Mb"
        allocate(NFCIHam(nFCIDet,nFCIDet))
        allocate(Np1FCIHam_alpha(nNp1FCIDet,nNp1FCIDet))
        allocate(Nm1FCIHam_beta(nNm1bFCIDet,nNm1bFCIDet))
        !call Fill_N_Np1_Nm1b_FCIHam(Elec,1,1,1,NHam=NFCIHam,Np1Ham=Np1FCIHam_alpha,Nm1Ham=Nm1FCIHam_beta)
        NFCIHam = zzero
        do i = 1,nFCIDet
            do j = 1,nFCIDet
                call GetHElement_comp(FCIDetList(:,i),FCIDetList(:,j),Elec,NFCIHam(i,j),  &
                    ilutnI=FCIBitList(i),ilutnJ=FCIBitList(j))
            enddo
        enddo
        Np1FCIHam_alpha = zzero
        do i = 1,nNp1FCIDet
            do j = 1,nNp1FCIDet
                call GetHElement_comp(Np1FCIDetList(:,i),Np1FCIDetList(:,j),Elec+1,Np1FCIHam_alpha(i,j),    &
                    ilutnI=Np1BitList(i),ilutnJ=Np1BitList(j))
            enddo
        enddo
        Nm1FCIHam_beta = zzero
        do i = 1,nNm1bFCIDet
            do j = 1,nNm1bFCIDet
                call GetHElement_comp(Nm1bFCIDetList(:,i),Nm1bFCIDetList(:,j),Elec-1,Nm1FCIHam_beta(i,j), &
                    ilutnI=Nm1bBitList(i),ilutnJ=Nm1bBitList(j))
            enddo
        enddo



        nLinearSystem = nNp1FCIDet + nFCIDet*2
        nLinearSystem_h = nNm1bFCIDet + nFCIDet*2

        if(nLinearSystem.ne.nLinearSystem_h) then
            call stop_all(t_r,'Should change restriction on the Linear system size being the same for particle/hole addition')
        endif
        
        write(6,"(A,F14.6,A)") "Memory required for the LR system: ",2*(real(nLinearSystem,dp)**2)*16/1048576.0_dp," Mb"
        
        !Set up orbital indices for SCHMIDT basis
        CoreEnd = nOcc-nImp
        VirtStart = nOcc+nImp+1
        VirtEnd = nSites
        ActiveStart = nOcc-nImp+1
        ActiveEnd = nOcc+nImp
        nCore = nOcc-nImp
        nVirt = nSites-nOcc-nImp   
        if(tHalfFill.and.(nCore.ne.nVirt)) then
            call stop_all(t_r,'Error in setting up half filled lattice')
        endif

        !Set up indices for the block of the linear system
        VIndex = nNp1FCIDet + 1   !Beginning of EC virtual excitations
        if(VIndex+nFCIDet-1.ne.nLinearSystem) call stop_all(t_r,'Indexing error')

        write(6,*) "External indices start from: ",VIndex
        write(6,*) "Total size of linear sys: ",nLinearSystem
        
        iunit = get_free_unit()
        call append_ext_real('EC-TDA_GFResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# 1.Frequency     2.GF_LinearResponse(Re)    3.GF_LinearResponse(Im)    " &
            & //"4.ParticleGF(Re)   5.ParticleGF(Im)   6.HoleGF(Re)   7.HoleGF(Im)    8.Old_GS    9.New_GS   " &
            & //"10.Particle_Norm  11.Hole_Norm   12.NI_GF(Re)   13.NI_GF(Im)  14.NI_GF_Part(Re)   15.NI_GF_Part(Im)   " &
            & //"16.NI_GF_Hole(Re)   17.NI_GF_Hole(Im)  18.SpectralWeight  19.Iters_p   20.Iters_h    " 
                
        allocate(SchmidtkGF_Cre_Ket(VirtStart:nSites))
        allocate(SchmidtkGF_Cre_Bra(VirtStart:nSites))
        allocate(SchmidtkGF_Ann_Ket(1:CoreEnd))
        allocate(SchmidtkGF_Ann_Bra(1:CoreEnd))
        
        !Allocate memory for hamiltonian in this system:
        allocate(LinearSystem_p(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(LinearSystem_h(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(Psi1_p(nLinearSystem),stat=ierr)
        allocate(Psi1_h(nLinearSystem),stat=ierr)
        Psi1_p = zzero
        Psi1_h = zzero
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        
        !Since this does not depend at all on the new basis, since the ground state has a completely
        !disjoint basis to |1>, the RHS of the equations can be precomputed
        allocate(Cre_0(nLinearSystem))
        allocate(Ann_0(nLinearSystem))
        
        !Store the fock matrix in complex form, so that we can ZGEMM easily
        !Also allocate the subblocks seperately. More memory, but will ensure we don't get stack overflow problems
        allocate(FockSchmidt_SE(nSites,nSites),stat=ierr)
        allocate(FockSchmidt_SE_VV(VirtStart:VirtEnd,VirtStart:VirtEnd),stat=ierr)    !Just defined over core-core blocks
        allocate(FockSchmidt_SE_CC(1:CoreEnd,1:CoreEnd),stat=ierr)    !Just defined over virtual-virtual blocks
        allocate(FockSchmidt_SE_VX(VirtStart:VirtEnd,ActiveStart:ActiveEnd),stat=ierr)
        allocate(FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:VirtEnd),stat=ierr)
        allocate(FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd),stat=ierr)
        allocate(FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd),stat=ierr)
        write(6,"(A,F12.5,A)") "Memory required for fock matrix: ", &
            (((2.0_dp*real(nSites,dp))**2)+(real(nOcc,dp)**2)+(real(nSites-nOcc,dp)**2))*ComptoMb, " Mb"
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        
        !Set the impurity parts to the correct values (even though we don't access them from FockSchmidt)
        !They are different since the correlation potential is not defined over the impurity sites.
        FockSchmidt(nOcc-nImp+1:nOcc+nImp,nOcc-nImp+1:nOcc+nImp) = Emb_h0v(:,:)
        !If we are doing self-consistent linear response, with a self-energy in the HL part, 
        !then this is calculated later, since it is now omega and iteration dependent
        do i = 1,nSites
            do j = 1,nSites
                FockSchmidt_SE(j,i) = dcmplx(FockSchmidt(j,i),0.0_dp)
            enddo
        enddo
        FockSchmidt_SE_CC(1:CoreEnd,1:CoreEnd) = FockSchmidt_SE(1:CoreEnd,1:CoreEnd)
        FockSchmidt_SE_VV(VirtStart:nSites,VirtStart:nSites) = FockSchmidt_SE(VirtStart:nSites,VirtStart:nSites)
        FockSchmidt_SE_VX(VirtStart:VirtEnd,ActiveStart:ActiveEnd) = FockSchmidt_SE(VirtStart:VirtEnd,ActiveStart:ActiveEnd)
        FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd) = FockSchmidt_SE(1:CoreEnd,ActiveStart:ActiveEnd)
        FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:VirtEnd) = FockSchmidt_SE(ActiveStart:ActiveEnd,VirtStart:VirtEnd)
        FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd) = FockSchmidt_SE(ActiveStart:ActiveEnd,1:CoreEnd)

        !Space for useful intermediates
        !For particle addition
        allocate(Gc_a_F_ax_Bra(ActiveStart:ActiveEnd),stat=ierr)
        allocate(Gc_a_F_ax_Ket(ActiveStart:ActiveEnd),stat=ierr)
        allocate(Gc_b_F_ab(VirtStart:VirtEnd),stat=ierr)
        !For hole addition
        allocate(Ga_i_F_xi_Bra(ActiveStart:ActiveEnd),stat=ierr)
        allocate(Ga_i_F_xi_Ket(ActiveStart:ActiveEnd),stat=ierr)
        allocate(Ga_i_F_ij(1:nCore),stat=ierr)

        !Allocate and precompute 1-operator coupling coefficients between the different sized spaces.
        !First number is the index of operator, and the second is the parity change when applying the operator
        allocate(Coup_Create_alpha(2,nFCIDet,nNm1bFCIDet))
        allocate(Coup_Ann_alpha(2,nFCIDet,nNp1FCIDet))
        Coup_Create_alpha(:,:,:) = 0
        Coup_Ann_alpha(:,:,:) = 0
        !TODO: This can be largly sped up by considering single excitations in isolation, rather than running through all elements
        !   This will result in an N_FCI * n_imp scaling, rather than N_FCI^2
        do J = 1,nFCIDet
            !Start with the creation of an alpha orbital
            do K = 1,nNm1bFCIDet
                DiffOrb = ieor(Nm1bBitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(mod(nOrbs,2).ne.1) then 
                    !There must be an odd number of orbital differences between the determinants
                    !since they are of different electron number by one
                    call stop_all(t_r,'Not odd number of electrons')
                endif
                if(nOrbs.ne.1) then
                    !We only want one orbital different
                    cycle
                endif
                !Now, find out what SPINorbital this one is from the bit number set
                call DecodeBitDet(orbdum,1,DiffOrb)
                gam = orbdum(1)
                if(mod(gam,2).ne.1) then
                    call stop_all(t_r,'differing orbital should be an alpha spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = Nm1bBitList(K)
                call SQOperator(tempK,gam,tParity,.false.)
                Coup_Create_alpha(1,J,K) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_alpha(2,J,K) = -1
                else
                    Coup_Create_alpha(2,J,K) = 1
                endif
            enddo
            !Now, annihilation of an alpha orbital
            do K = 1,nNp1FCIDet
                DiffOrb = ieor(Np1BitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(mod(nOrbs,2).ne.1) then 
                    !There must be an odd number of orbital differences between the determinants
                    !since they are of different electron number by one
                    call stop_all(t_r,'Not odd number of electrons 3')
                endif
                if(nOrbs.ne.1) then
                    !We only want one orbital different
                    cycle
                endif
                !Now, find out what SPINorbital this one is from the bit number set
                call DecodeBitDet(orbdum,1,DiffOrb)
                gam = orbdum(1)
                if(mod(gam,2).ne.1) then
                    call stop_all(t_r,'differing orbital should be an alpha spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = Np1BitList(K)
                call SQOperator(tempK,gam,tParity,.true.)
                Coup_Ann_alpha(1,J,K) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_alpha(2,J,K) = -1
                else
                    Coup_Ann_alpha(2,J,K) = 1
                endif
            enddo
        enddo
        i = 2*nFCIDet*nNm1bFCIDet + 2*nFCIDet*nNp1FCIDet
        write(6,"(A,F12.5,A)") "Memory required for coupling coefficient matrices: ",real(i,dp)*RealtoMb, " Mb"

        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        call halt_timer(LR_EC_GF_Precom)

        do while(.true.)
            !call GetNextkVal(kPnt,tFinishedk)
            if(tFinishedk) exit
            write(6,"(A)") "Calculating spectra with k = ",KPnts(:,kPnt)
        
            call FindStaticMomSchmidtPert(kPnt,StatickGF_Cre_Ket,StatickGF_Ann_Ket)

            Omega = Start_Omega
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                write(6,"(A)") "Calculating linear response for frequency: ",Omega
                call flush(6)

                call FindMomGFSchmidtPert(Omega,ni_lr_cre,ni_lr_ann,kPnt,SchmidtkGF_Cre_Ket,SchmidtkGF_Ann_Ket)

                Omega = Omega + Omega_Step
            enddo

        enddo


    end subroutine MomGF_Ex

    !Construct the contraction coefficients for the action of the static perturbation 
    subroutine FindStaticMomSchmidtPert(kPnt,StatickGF_Cre_Ket,StatickGF_Ann_Ket)
        implicit none
        integer, intent(in) :: kPnt
        complex(dp), intent(out) :: StatickGF_Cre_Ket(nOcc+nImp+1:nSites),StatickGF_Ann_Ket(1:nOcc-nImp)
        complex(dp), allocatable :: HFPertBasis_Ann(:),HFPertBasis_Cre(:)
        integer :: SS_Period,ind_1,ind_2,i,x,n
        
        SS_Period = nImp
        allocate(HFPertBasis_Ann(SS_Period))
        allocate(HFPertBasis_Cre(SS_Period))
        HFPertBasis_Ann(:) = zzero
        HFPertBasis_Cre(:) = zzero
            
        ind_1 = ((kPnt-1)*SS_Period) + 1
        ind_2 = SS_Period * kPnt

        do i = 1,SS_Period  !Run over vectors which span this kpoint
            if(KVec_EMapping(ind_1+i-1).le.nOcc) then
                !i is an occupied orbital
                do x = 1,SS_Period  !Run over components of this vector
                    HFPertBasis_Ann(i) = HFPertBasis_Ann(i) + dconjg(k_vecs(x,ind_1 + i - 1))
                enddo
            else
                !i is a virtual orbital
                do x = 1,SS_Period  !Run over components of this vector
                    HFPertBasis_Cre(i) = HFPertBasis_Cre(i) + dconjg(k_vecs(x,ind_1 + i - 1))
                enddo
            endif 
        enddo

        !Now rotate from the HF vectors, to the Schmidt basis
        StatickGF_Ann_Ket(:) = zzero
        do n = 1,nOcc-nImp  !Loop over schmidt basis of core orbitals
            do i = 1,SS_Period  !Contract over eigenvectors corresponding to this kpoint
                StatickGF_Ann_Ket(n) = StatickGF_Ann_Ket(n) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Ann(i)
            enddo
        enddo
        StatickGF_Cre_Ket(:) = zzero
        do n = nOcc+nImp+1,nSites
            do i = 1,SS_Period
                StatickGF_Cre_Ket(n) = StatickGF_Cre_Ket(n) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Cre(i)
            enddo
        enddo

        deallocate(HFPertBasis_Ann,HFPertBasis_Cre)

    end subroutine FindStaticMomSchmidtPert

    subroutine FindMomGFSchmidtPert(Omega,ni_lr_cre,ni_lr_ann,kPnt,SchmidtkGF_Cre_Ket,SchmidtkGF_Ann_Ket)
        implicit none
        integer, intent(in) :: kPnt !This is the index of the kpoint perturbation
        real(dp), intent(in) :: Omega
        complex(dp), intent(out) :: ni_lr_cre,ni_lr_ann
        complex(dp), intent(out) :: SchmidtkGF_Cre_Ket(nOcc+nImp+1:nSites),SchmidtkGF_Ann_Ket(1:nOcc-nImp)
        integer :: ind_1,ind_2,SS_Period,i,x,n
        complex(dp), allocatable :: HFPertBasis_Ann(:),HFPertBasis_Cre(:)

        ni_lr_cre = zzero
        ni_lr_ann = zzero

        SS_Period = nImp
        allocate(HFPertBasis_Ann(SS_Period))
        allocate(HFPertBasis_Cre(SS_Period))
        HFPertBasis_Ann(:) = zzero
        HFPertBasis_Cre(:) = zzero
            
        ind_1 = ((kPnt-1)*SS_Period) + 1
        ind_2 = SS_Period * kPnt

        do i = 1,SS_Period  !Run over vectors which span this kpoint
            if(KVec_EMapping(ind_1+i-1).le.nOcc) then
                !i is an occupied orbital
                do x = 1,SS_Period  !Run over components of this vector
                    HFPertBasis_Ann(i) = HFPertBasis_Ann(i) + dconjg(k_vecs(x,ind_1 + i - 1)) /     &
                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                    ni_lr_ann = ni_lr_ann + (k_vecs(x,ind_1 + i - 1)*dconjg(k_vecs(x,ind_1 + i - 1))) / &
                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                enddo
            else
                !i is a virtual orbital
                do x = 1,SS_Period  !Run over components of this vector
                    HFPertBasis_Cre(i) = HFPertBasis_Cre(i) + dconjg(k_vecs(x,ind_1 + i - 1)) /     &
                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                    ni_lr_cre = ni_lr_cre + (k_vecs(x,ind_1 + i - 1)*dconjg(k_vecs(x,ind_1 + i - 1))) / &
                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                enddo
            endif
        enddo

        !Now rotate from the HF vectors, to the Schmidt basis
        SchmidtkGF_Ann_Ket(:) = zzero
        do n = 1,nOcc-nImp  !Loop over schmidt basis of core orbitals
            do i = 1,SS_Period  !Contract over eigenvectors corresponding to this kpoint
                SchmidtkGF_Ann_Ket(n) = SchmidtkGF_Ann_Ket(n) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Ann(i)
            enddo
        enddo
        SchmidtkGF_Cre_Ket(:) = zzero
        do n = nOcc+nImp+1,nSites
            do i = 1,SS_Period
                SchmidtkGF_Cre_Ket(n) = SchmidtkGF_Cre_Ket(n) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Cre(i)
            enddo
        enddo

        deallocate(HFPertBasis_Ann,HFPertBasis_Cre)
        
    end subroutine FindMomGFSchmidtPert

end module MomSpectra
