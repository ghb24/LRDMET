module MomSpectra
    use const
    use timing
    use globals
    use errors, only: stop_all,warning
    use mat_tools, only: WriteVector,WriteMatrix,WriteVectorInt,WriteMatrixComp,WriteVectorComp,znrm2
    use LRSolvers
    implicit none

    contains

    !Routine to read in a correlated greens function, and self-consistently converge
    !a self energy from the 1-electron hemailtonian to get it to match this.
    subroutine SC_Mom_LR()
        use utils, only: get_free_unit,append_ext_real,append_ext
        implicit none
        integer :: iunit,i,nESteps
        logical :: exists
        complex(dp), allocatable :: G00(:),SE(:),Hybrid(:),LocalMomGF(:)
        real(dp) :: Omega,Re_LR,Im_LR,Omega_Val
        character(64) :: filename,filename2
        character(128) :: header
        character(len=*), parameter :: t_r='SC_Mom_LR'

        !First, read back in the G_00 (real and complex)
        iunit = get_free_unit()
        call append_ext_real('EC-TDA_GFResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        inquire(file=filename2,exist=exists)
        if(.not.exists) then
            call stop_all(t_r,'Cannot find local greens function file')
        endif
        open(unit=iunit,file=filename2,status='old',action='read')

        read(iunit,*) header
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
        enddo
        i = nESteps
        allocate(G00(nImp,nImp,nESteps))    !The correlated greens function
        allocate(SE(nImp,nImp,nESteps))     !The self energy
        allocate(Hybrid(nImp,nImp,nESteps)) !The hybridization
        allocate(LocalMomGF(nImp,nImp,nESteps)) !The local 1-electron GF from the non-interacting H
        G00(:,:,:) = zzero
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.nESteps) call stop_all(t_r,'Too many frequency points')
            read(iunit,*) Omega_val,Re_LR,Im_LR
            G00(1,1,i) = dcmplx(Re_LR,Im_LR)
        enddo
        close(iunit)

        SE(:,:,:) = zzero
        Hybrid(:,:,:) = zzero
        LocalMomGF(:,:,:) = zzero

        do while(.true.)

            !Construct FT of k-space GFs and take zeroth part.
            call FindLocalMomGF(nESteps,SE,LocalMomGF)

            !We now have the updated local greens function from the 1-electron hamiltonian
            !Exit criterion here once the local greens function is converged.

            !Invert the matrix of non-interacting local greens functions.
            InvLocalMomGF(:,:,:) = LocalMomGF(:,:,:)
            call InvertLocalNonHermGF(nESteps,InvLocalMomGF)

            !Now find hybridization
            !This is given by (omega + mu + idelta - e_0 - SE - InvGF)
            call FindHybrid(nESteps,InvLocalMomGF,SE,Hybrid)

            !Now to converge the k-independent Self-energy
            InvG00(:,:,:) = G00(:,:,:)
            call InvertLocalNonHermGF(nESteps,InvG00)
            !Now find Self-energy
            !This is given by (omega + mu + idelta)I - e_0 - Hybrid - InvFullGF)
            call FindSE(nESteps,InvG00,Hybrid,SE)

        enddo

    end subroutine SC_Mom_LR

    !This finds the self energy required to match the local greens function to the interacting local one.
    ! (omega + mu +- idelta)I - e_00 - hybrid - InvG00      Gives the local self-energy
    ! e_00 is taken from the local h0v part of tha hamiltonian
    subroutine FindSE(n,InvG00,Hybrid,SE)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: InvG00(nImp,nImp,n)
        complex(dp), intent(in) :: Hybrid(nImp,nImp,n)
        complex(dp), intent(out) :: SE(nImp,nImp,n)
        real(dp) :: Omega,mu
        character(len=*), parameter :: t_r='FindSE'

        !This is *NOT* the correct chemical potential for non-half filled systems, but nevermind...
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        SE(:,:,:) = zzero

        do i = 1,nImp
            do j = 1,nImp
                LocalH(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo

        i = 0
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Too many freq points')

            SE(:,:,i) = - LocalH(:,:) - Hybrid(:,:,i) - InvG00(:,:,i)
            do j = 1,nImp
                SE(j,j,i) = SE(j,j,i) + dcmplx(Omega + mu,dDelta)
            enddo

            Omega = Omega + Omega_Step
        enddo
    end subroutine FindSE

    !This hybridization is a nImp x nImp matrix, for all frequency points, calculated from:
    ! [omega + mu + i\eta]I - e_00 - SE - Inv_GF     where e_00 is the impurity block of the 'non-interacting' hamiltonian
    ! Reminder that here we are performing this self-consistency on the R£AL frequency axis, and the e_00 matrix simply comes from h0v
    ! Hybridization returned in Hybrid
    subroutine FindHybrid(n,InvLocalMomGF,SE,Hybrid)
        implicit none
        integer, intent(in) :: n    !number of frequency points
        complex(dp), intent(in) :: InvLocalMomGF(nImp,nImp,n)
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(out) :: Hybrid(nImp,nImp,n)
        complex(dp) :: LocalH(nImp,nImp)
        real(dp) :: mu,Omega
        integer :: i
        character(len=*), parameter :: t_r='FindHybrid'
        
        !This is *NOT* the correct chemical potential for non-half filled systems, but nevermind...
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        Hybrid(:,:,:) = zzero

        do i = 1,nImp
            do j = 1,nImp
                LocalH(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo

        i = 0
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            if(i.gt.n) call stop_all(t_r,'Too many freq points')

            Hybrid(:,:,i) = - LocalH(:,:) - SE(:,:,i) - InvLocalMomGF(:,:,i)
            do j = 1,nImp
                Hybrid(j,j,i) = Hybrid(j,j,i) + dcmplx(Omega + mu,dDelta)
            enddo

            Omega = Omega + Omega_Step
        enddo

    end subroutine FindHybrid
            
    !This function inverts a local greens function matrix, allowing it to be non-hermitian
    !Function is sent in as the normal greens function, n, nImp x nImp matrices, and the inverse
    !for each frequency is returned
    subroutine InvertLocalNonHermGF(n,InvGF)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(inout) :: InvGF(nImp,nImp,n)
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: LVec(:,:),RVec(:,:),W_Vals(:),cWork(:)
        integer :: i,lwork,info,j

        allocate(LVec(nImp,nImp))
        allocate(RVec(nImp,nImp))
        allocate(W_Vals(nImp))
        allocate(Work(max(1,2*nImp)))

        do i = 1,n

            !Diagonalize
            LVec(:,:) = zzero
            RVec(:,:) = zzero
            W_Vals(:) = zzero

            allocate(cWork(1))
            lwork = -1
            info = 0
            call zgeev('V','V',nImp,InvGF(:,:,i),nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Workspace query failed')
            lwork = int(abs(cWork(1)))+1
            deallocate(cWork)
            allocate(cWork(lWork))
            call zgeev('V','V',nImp,InvGF(:,:,i),nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
            deallocate(cWork)

            !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
            !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
            call Order_zgeev_vecs(W_Vals,LVec,RVec)
            !call writevectorcomp(W_Vals,'Eigenvalues ordered')
            !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
            call Orthonorm_zgeev_vecs(nImp,W_Vals,LVec,RVec)
            !Calculate greens function for this k-vector
            InvMat(:,:) = zzero
            do j = 1,nImp
                InvGF(j,j,i) = zone / W_Vals(j)
            enddo
            !Now rotate this back into the original basis
            call zGEMM('N','N',nImp,nImp,nImp,zone,RVec,nImp,InvGF(:,:,i),nImp,zzero,ztemp2,nImp)
            call zGEMM('N','C',nImp,nImp,nImp,zone,ztemp2,nImp,LVec,nImp,zzero,InvGF(:,:,i),nImp)

        enddo

        deallocate(LVec,RVec,W_Vals,Work)

    end subroutine InvertLocalNonHermGF
        
    !This calculates the sum of all the k-space greens functions over the 'n' frequency points
    !The hamiltonian in this case is taken to be h0v (i.e. with the correlation potential)
    !1/V \sum_k [w + mu + i/eta - h_k - SE]^-1      where V is volume of the BZ
    !TODO: Care needed with sign of broadening?
    !Check that this gives as expected
    subroutine FindLocalMomGF(n,SE,LocalMomGF)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: SE(nImp,nImp,n)
        complex(dp), intent(out) :: LocalMomGF(nImp,nImp,n)
        integer :: kPnt,ind_1,ind_2,i,j,SS_Period
        real(dp) :: Omega,mu
        complex(dp) :: ni_lr_ann(n),ni_lr_cre(n)
        character(len=*), parameter :: t_r='FindLocalMomGF'

        LocalMomGF(:,:,:) = zzero
        SS_Period = nImp
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        allocate(RotMat(nLat,SS_Period))
        allocate(k_Ham(SS_Period,SS_Period))
        allocate(CompHam(nLat,nLat))
        do i = 1,nLat
            do j = 1,nLat
                CompHam(j,i) = dcmplx(h0v,zero)
            enddo
        enddo

        !Data for the diagonalization of the one-electron k-dependent Greens functions
        allocate(RVec(SS_Period,SS_Period))
        allocate(LVec(SS_Period,SS_Period))
        allocate(W_Vals(SS_Period))
        allocate(Work(max(1,2*SS_Period)))

        do kPnt = 1,nKPnts
            !Run through all k-points
            ind_1 = ((kPnt-1)*SS_Period) + 1
            ind_2 = SS_Period*kPnt
                
            !Zero the individual greens functions
            ni_lr_ann(:) = zzero
            ni_lr_cre(:) = zzero
                
            !Transform the one-electron hamiltonian into this k-basis
            RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
            call ZGEMM('N','N',nLat,SS_Period,nLat,zone,CompHam,nLat,RotMat,nLat,zzero,ztemp,nLat)
            call ZGEMM('C','N',SS_Period,SS_Period,nLat,zone,RotMat,nLat,ztemp,nLat,zzero,k_Ham,SS_Period)
            !k_Ham is the complex, one-electron hamiltonian (with correlation potential) for this k-point

            i = 0
            Omega = Start_Omega
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                if(i.gt.n) call stop_all(t_r,'Too many freq points')

                !Find inverse greens function for this k-point
                InvMat(:,:) = - k_Ham(:,:) - SE(:,:,i)
                do i = 1,SS_Period
                    InvMat(i,i) = InvMat(i,i) + dcmplx(Omega + mu,dDelta)
                enddo

                !Now, diagonalize this
                !The self-energy is *not* hermitian
                allocate(cWork(1))
                lwork = -1
                info = 0
                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(abs(cWork(1)))+1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zgeev('V','V',SS_Period,InvMat,SS_Period,W_Vals,LVec,SS_Period,RVec,SS_Period,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
                deallocate(cWork)

                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
                call Order_zgeev_vecs(W_Vals,LVec,RVec)
                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
                call Orthonorm_zgeev_vecs(SS_Period,W_Vals,LVec,RVec)
                !Calculate greens function for this k-vector
                InvMat(:,:) = zzero
                do i = 1,SS_Period
                    InvMat(i,i) = zone / W_Vals(i)
                enddo
                !Now rotate this back into the original basis
                call zGEMM('N','N',SS_Period,SS_Period,SS_Period,zone,RVec,SS_Period,InvMat,SS_Period,zzero,ztemp2,SS_Period)
                call zGEMM('N','C',SS_Period,SS_Period,SS_Period,zone,ztemp2,SS_Period,LVec,SS_Period,zzero,InvMat,SS_Period)
                !InvMat is now the non-interacting greens function for this k: TEST THIS!
                !Sum this into the *local* greens function (i.e. a fourier transform of r=0 component)
                LocalMomGF(:,:,i) = LocalMomGF(:,:,i) + InvMat(:,:)

                Omega = Omega + Omega_Step
            enddo

        enddo

        !Divide the entire local greens function by the 'volume' of the brillouin zone
        LocalMomGF(:,:,:) = LocalMomGF(:,:,:) / BZVol

    end subroutine FindLocalMomGF

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
        complex(dp), allocatable :: Overlap_p(:,:),Overlap_h(:,:)
        complex(dp), allocatable :: FockSchmidt_SE(:,:),FockSchmidt_SE_VV(:,:)
        complex(dp), allocatable :: FockSchmidt_SE_CC(:,:),FockSchmidt_SE_VX(:,:),FockSchmidt_SE_XV(:,:)
        complex(dp), allocatable :: FockSchmidt_SE_CX(:,:),FockSchmidt_SE_XC(:,:)
        complex(dp), allocatable :: Gc_a_F_ax_Bra(:,:),Gc_a_F_ax_Ket(:,:),Gc_b_F_ab(:,:)
        complex(dp), allocatable :: Ga_i_F_xi_Bra(:,:),Ga_i_F_xi_Ket(:,:),Ga_i_F_ij(:,:)
        complex(dp), allocatable :: Xc_a_F_ax_Bra(:,:),Xc_a_F_ax_Ket(:,:),Xc_b_F_ab(:,:)
        complex(dp), allocatable :: Xa_i_F_xi_Bra(:,:),Xa_i_F_xi_Ket(:,:),Xa_i_F_ij(:,:)
        complex(dp), allocatable :: SchmidtkGF_Cre_Ket(:,:),SchmidtkGF_Ann_Ket(:,:)
        complex(dp), allocatable :: SchmidtkGF_Cre_Bra(:,:),SchmidtkGF_Ann_Bra(:,:)
        complex(dp), allocatable :: StatickGF_Cre_Ket(:,:),StatickGF_Ann_Ket(:,:)
        complex(dp), allocatable :: StatickGF_Cre_Bra(:,:),StatickGF_Ann_Bra(:,:)
        complex(dp), allocatable :: StatickGF_Ann_Emb_Ket(:,:),StatickGF_Ann_Emb_Bra(:,:)
        complex(dp), allocatable :: StatickGF_Cre_Emb_Ket(:,:),StatickGF_Cre_Emb_Bra(:,:)
        complex(dp), allocatable :: Psi_0(:),V0_Ann(:),V0_Cre(:)
        complex(dp), allocatable :: S_EigVec(:,:),tempc(:),DynamicLS(:,:)
        real(dp), allocatable :: S_EigVal(:),Work(:)
        integer, allocatable :: Coup_Create_alpha(:,:,:),Coup_Ann_alpha(:,:,:)
        integer :: nLinearSystem,nLinearSystem_h,CoreEnd,VirtStart,VirtEnd,ActiveStart,ActiveEnd,nCore,nVirt
        integer :: VIndex,iunit,i,j,k,DiffOrb,nOrbs,orbdum(1),gam,tempK,KPnt,ierr,nImp_GF,VStatIndex,info,lWork
        integer :: nSpan,a,LinSize
        complex(dp) :: ni_lr_h,ni_lr_p,ni_lr,CNorm,CStatNorm,VNorm,VstatNorm,StatCoup_Ann,StatCoup_Cre
        complex(dp) :: tempel,tempel2,tempel4,dNorm_p,dNorm_h,Response_h,Response_p,ResponseFn
        real(dp) :: SpectralWeight,Prev_Spec,mu,Omega,GFChemPot
        logical :: tFirst,tParity,tFinishedK,tNoCre,tNoAnn
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

        nImp_GF = 1 !Only want one perturbation at a time, but keep as matrices

        !Set up indices for the block of the linear system
        VIndex = nNp1FCIDet + 1   !Beginning of EC virtual excitations
        VStatIndex = VIndex + nFCIDet
        if(VStatIndex+nFCIDet-1.ne.nLinearSystem) call stop_all(t_r,'Indexing error')

        write(6,*) "Total size of linear sys: ",nLinearSystem
        
        iunit = get_free_unit()
        call append_ext_real('Mom_GF',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# 1.Frequency     2.GF(Re)    3.GF(Im)    " &
            & //"4.ParticleGF(Re)   5.ParticleGF(Im)   6.HoleGF(Re)   7.HoleGF(Im)    8.Old_GS    9.New_GS   " &
            & //"10.Particle_Norm  11.Hole_Norm   12.NI_GF(Re)   13.NI_GF(Im)  14.NI_GF_Part(Re)   15.NI_GF_Part(Im)   " &
            & //"16.NI_GF_Hole(Re)   17.NI_GF_Hole(Im)  18.SpectralWeight  19.Iters_p   20.Iters_h    " 
                
        allocate(SchmidtkGF_Cre_Ket(VirtStart:nSites,nImp_GF))
        allocate(SchmidtkGF_Cre_Bra(VirtStart:nSites,nImp_GF))
        allocate(SchmidtkGF_Ann_Ket(1:CoreEnd,nImp_GF))
        allocate(SchmidtkGF_Ann_Bra(1:CoreEnd,nImp_GF))

        allocate(StatickGF_Cre_Ket(VirtStart:nSites,nImp_GF))
        allocate(StatickGF_Cre_Bra(VirtStart:nSites,nImp_GF))
        allocate(StatickGF_Ann_Ket(1:CoreEnd,nImp_GF))
        allocate(StatickGF_Ann_Bra(1:CoreEnd,nImp_GF))
        !Also for the Emb space for the static perturbation
        allocate(StatickGF_Cre_Emb_Ket(ActiveStart:ActiveEnd,nImp_GF))
        allocate(StatickGF_Cre_Emb_Bra(ActiveStart:ActiveEnd,nImp_GF))
        allocate(StatickGF_Ann_Emb_Ket(ActiveStart:ActiveEnd,nImp_GF))
        allocate(StatickGF_Ann_Emb_Bra(ActiveStart:ActiveEnd,nImp_GF))
        
        !Allocate memory for hamiltonian in this system:
        allocate(LinearSystem_p(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(LinearSystem_h(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(Overlap_p(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(Overlap_h(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(Psi1_p(nLinearSystem),stat=ierr)
        allocate(Psi1_h(nLinearSystem),stat=ierr)
        Psi1_p = zzero
        Psi1_h = zzero
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        
        !Since this does not depend at all on the new basis, since the ground state has a completely
        !disjoint basis to |1>, the RHS of the equations can be precomputed
        allocate(Psi_0(nFCIDet))
        GFChemPot = HL_Energy
        do i = 1,nFCIDet
            Psi_0(i) = dcmplx(HL_Vec(i),zero)
        enddo
        allocate(V0_Cre(nLinearSystem))
        allocate(V0_Ann(nLinearSystem))
        
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

        !Space for useful intermediates for dynamic bath states
        !For particle addition
        allocate(Gc_a_F_ax_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Gc_a_F_ax_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Gc_b_F_ab(VirtStart:VirtEnd,nImp_GF),stat=ierr)
        !For hole addition
        allocate(Ga_i_F_xi_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Ga_i_F_xi_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Ga_i_F_ij(1:nCore,nImp_GF),stat=ierr)
        !Space for useful intermediates for static bath states
        !For particle addition
        allocate(Xc_a_F_ax_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Xc_a_F_ax_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Xc_b_F_ab(VirtStart:VirtEnd,nImp_GF),stat=ierr)
        !For hole addition
        allocate(Xa_i_F_xi_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Xa_i_F_xi_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Xa_i_F_ij(1:nCore,nImp_GF),stat=ierr)

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

        kPnt = 0
        do while(.true.)
            call GetNextkVal(kPnt,tFinishedk)
            if(tFinishedk) exit
            write(6,*) "Calculating spectra with k = ",KPnts(:,KPnt),KPnt
            call WriteKVecHeader(iunit,KPnts(:,KPnt))

            call FindStaticMomSchmidtPert(kPnt,nImp_GF,StatickGF_Cre_Ket,StatickGF_Cre_Emb_Ket,StatickGF_Ann_Ket,   &
                StatickGF_Ann_Emb_Ket,StatickGF_Cre_Bra,StatickGF_Ann_Bra,StatickGF_Cre_Emb_Bra,StatickGF_Ann_Emb_Bra,tNoCre,tNoAnn)

            !Construct useful intermediates, concerning the *Static* bath orbitals
            !sum_a Xc_a^* F_ax (Creation)
            call ZGEMM('T','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_VX(VirtStart:VirtEnd,ActiveStart:ActiveEnd),  &
                nVirt,StatickGF_Cre_Bra(VirtStart:VirtEnd,:),nVirt,zzero,Xc_a_F_ax_Bra,EmbSize)
            call ZGEMM('N','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:VirtEnd),  &
                EmbSize,StatickGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Xc_a_F_ax_Ket,EmbSize)
            !sum_b Gc_b F_ab  (Creation)
            call ZGEMM('N','N',nVirt,nImp_GF,nVirt,zone,FockSchmidt_SE_VV(VirtStart:VirtEnd,VirtStart:VirtEnd),    &
                nVirt,StatickGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Xc_b_F_ab,nVirt)

            !sum_i Ga_i^bra F_xi (Annihilation)
            call ZGEMM('N','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd),EmbSize,    &
                StatickGF_Ann_Bra(1:nCore,:),nCore,zzero,Xa_i_F_xi_Bra,EmbSize)
            call ZGEMM('T','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd),nCore,    &
                StatickGF_Ann_Ket(1:nCore,:),nCore,zzero,Xa_i_F_xi_Ket,EmbSize)
            !sum_i Ga_i F_ij (Annihilation)
            call ZGEMM('T','N',nCore,nImp_GF,nCore,zone,FockSchmidt_SE_CC(:,:),nCore,StatickGF_Ann_Ket(:,:),nCore,zzero, &
                Xa_i_F_ij,nCore)

            Omega = Start_Omega
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                write(6,"(A,F12.5,I6)") "Calculating linear response for frequency: ",Omega,KPnt
                call flush(6)

                call FindMomGFSchmidtPert(Omega,ni_lr_p,ni_lr_h,kPnt,nImp_GF,SchmidtkGF_Cre_Ket,SchmidtkGF_Ann_Ket, &
                    SchmidtkGF_Cre_Bra,SchmidtkGF_Ann_Bra)
                
                !Construct useful intermediates, concerning the *Dynamic* bath orbitals
                !sum_a Gc_a^* F_ax (Creation)
                call ZGEMM('T','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_VX(VirtStart:VirtEnd,ActiveStart:ActiveEnd),  &
                    nVirt,SchmidtkGF_Cre_Bra(VirtStart:VirtEnd,:),nVirt,zzero,Gc_a_F_ax_Bra,EmbSize)
                call ZGEMM('N','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:VirtEnd),  &
                    EmbSize,SchmidtkGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Gc_a_F_ax_Ket,EmbSize)
                !sum_b Gc_b F_ab  (Creation)
                call ZGEMM('N','N',nVirt,nImp_GF,nVirt,zone,FockSchmidt_SE_VV(VirtStart:VirtEnd,VirtStart:VirtEnd),    &
                    nVirt,SchmidtkGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Gc_b_F_ab,nVirt)

                !sum_i Ga_i^bra F_xi (Annihilation)
                call ZGEMM('N','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd),EmbSize,    &
                    SchmidtkGF_Ann_Bra(1:nCore,:),nCore,zzero,Ga_i_F_xi_Bra,EmbSize)
                call ZGEMM('T','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd),nCore,    &
                    SchmidtkGF_Ann_Ket(1:nCore,:),nCore,zzero,Ga_i_F_xi_Ket,EmbSize)
                !sum_i Ga_i F_ij (Annihilation)
                call ZGEMM('T','N',nCore,nImp_GF,nCore,zone,FockSchmidt_SE_CC(:,:),nCore,SchmidtkGF_Ann_Ket(:,:),nCore,zzero, &
                    Ga_i_F_ij,nCore)
                StatCoup_Ann = zzero
                StatCoup_Cre = zzero
                do i = 1,CoreEnd
                    !sum_i X*_i G_i
                    StatCoup_Ann = StatCoup_Ann + StatickGF_Ann_Bra(i,1)*SchmidtkGF_Ann_Ket(i,1)
                enddo
                do a = VirtStart,VirtEnd
                    StatCoup_Cre = StatCoup_Cre + StatickGF_Cre_Bra(a,1)*SchmidtkGF_Cre_Ket(a,1)
                enddo

                !Now, construct the hessians
                LinearSystem_h = zzero
                LinearSystem_p = zzero
                !Overlap matrix
                Overlap_h(:,:) = zzero
                Overlap_p(:,:) = zzero
                do i = 1,nLinearSystem
                    Overlap_h(i,i) = zone
                    Overlap_p(i,i) = zone
                enddo

                !Block 1
                !First, construct n + 1 (alpha) FCI space, in determinant basis
                LinearSystem_p(1:nNp1FCIDet,1:nNp1FCIDet) = Np1FCIHam_alpha(:,:)
                !Block 1 for hole hamiltonian
                LinearSystem_h(1:nNm1bFCIDet,1:nNm1bFCIDet) = Nm1FCIHam_beta(:,:)

                !Block 2
                !Copy the N electron FCI hamiltonian to this diagonal blockk
                LinearSystem_h(VIndex:VStatIndex-1,VIndex:VStatIndex-1) = nFCIHam(:,:)
                LinearSystem_p(VIndex:VStatIndex-1,VIndex:VStatIndex-1) = nFCIHam(:,:)

                CNorm = zzero
                VNorm = zzero
                tempel = zzero 
                tempel2 = zzero
                do i = 1,CoreEnd
                    !Calc normalization
                    CNorm = CNorm + SchmidtkGF_Ann_Bra(i,1)*SchmidtkGF_Ann_Ket(i,1)
                    !Calc diagonal correction
                    tempel = tempel + SchmidtkGF_Ann_Bra(i,1)*Ga_i_F_ij(i,1)
                enddo
                tempel = -tempel / CNorm
                do a = VirtStart,VirtEnd
                    !Calc particle norm
                    VNorm = VNorm + SchmidtkGF_Cre_Bra(a,1)*SchmidtkGF_Cre_Ket(a,1)
                    !Diag particle correction
                    tempel2 = tempel2 + SchmidtkGF_Cre_Bra(a,1)*Gc_b_F_ab(a,1)
                enddo
                tempel2 = tempel2 / VNorm
                !Add diagonal correction
                do i = VIndex,VStatIndex-1
                    LinearSystem_h(i,i) = LinearSystem_h(i,i) + tempel
                    LinearSystem_p(i,i) = LinearSystem_p(i,i) + tempel2
                enddo

                !Block 4
                !This is essentially the same, but now for the static bath
                LinearSystem_h(VStatIndex:nLinearSystem,VStatIndex:nLinearSystem) = nFCIHam(:,:)
                LinearSystem_p(VStatIndex:nLinearSystem,VStatIndex:nLinearSystem) = nFCIHam(:,:)

                CStatNorm = zzero
                VstatNorm = zzero
                tempel = zzero 
                tempel2 = zzero
                do i = 1,CoreEnd
                    !Calc normalization
                    CStatNorm = CStatNorm + StatickGF_Ann_Bra(i,1)*StatickGF_Ann_Ket(i,1)
                    !Calc diagonal correction
                    tempel = tempel + StatickGF_Ann_Bra(i,1)*Xa_i_F_ij(i,1)
                enddo
                tempel = -tempel / CStatNorm
                do a = VirtStart,VirtEnd
                    !Calc particle static normalization
                    VstatNorm = VstatNorm + StatickGF_Cre_Bra(a,1)*StatickGF_Cre_Ket(a,1)
                    !Diag particle correction
                    tempel2 = tempel2 + StatickGF_Cre_Bra(a,1)*Xc_b_F_ab(a,1)
                enddo
                tempel2 = tempel2 / VstatNorm
                !Add diagonal correction
                do i = VStatIndex,nLinearSystem
                    LinearSystem_h(i,i) = LinearSystem_h(i,i) + tempel
                    LinearSystem_p(i,i) = LinearSystem_p(i,i) + tempel2
                enddo

                !Block 3
                do J = 1,nNm1bFCIDet
                    do K = 1,nFCIDet
                        if(Coup_Create_alpha(1,K,J).ne.0) then
                            !Determinants are connected via a single SQ operator
                            LinearSystem_h(K+VIndex-1,J) = - Ga_i_F_xi_Bra(Coup_Create_alpha(1,K,J),1)   &
                                *Coup_Create_alpha(2,K,J)/sqrt(CNorm)
                            LinearSystem_h(J,K+VIndex-1) = - Ga_i_F_xi_Ket(Coup_Create_alpha(1,K,J),1)   &
                                *Coup_Create_alpha(2,K,J)/sqrt(CNorm)
                            !Block 6 is essentially the same
                            LinearSystem_h(K+VStatIndex-1,J) = - Xa_i_F_xi_Bra(Coup_Create_alpha(1,K,J),1)   &
                                *Coup_Create_alpha(2,K,J)/sqrt(CStatNorm)
                            LinearSystem_h(J,K+VStatIndex-1) = - Xa_i_F_xi_Ket(Coup_Create_alpha(1,K,J),1)   &
                                *Coup_Create_alpha(2,K,J)/sqrt(CStatNorm)
                        endif
                    enddo
                enddo

                !Block 3 - particle
                do J = 1,nNp1FCIDet
                    do K = 1,nFCIDet
                        if(Coup_Ann_alpha(1,K,J).ne.0) then
                            !Determinants are connected via a single SQ operator
                            !First index is the spatial index, second is parity
                            LinearSystem_p(K+VIndex-1,J) = Gc_a_F_ax_Bra(Coup_Ann_alpha(1,K,J),1)    &
                                *Coup_Ann_alpha(2,K,J)/sqrt(VNorm)
                            !LinearSystem_p(J,K+VIndex-1) = conjg(LinearSystem_p(K+VIndex-1,J))
                            LinearSystem_p(J,K+VIndex-1) = Gc_a_F_ax_Ket(Coup_Ann_alpha(1,K,J),1)    &
                                *Coup_Ann_alpha(2,K,J)/sqrt(VNorm)
                            !Block 6 essentially the same
                            LinearSystem_p(K+VStatIndex-1,J) = Xc_a_F_ax_Bra(Coup_Ann_alpha(1,K,J),1)    &
                                *Coup_Ann_alpha(2,K,J)/sqrt(VstatNorm)
                            LinearSystem_p(J,K+VStatIndex-1) = Xc_a_F_ax_Ket(Coup_Ann_alpha(1,K,J),1)    &
                                *Coup_Ann_alpha(2,K,J)/sqrt(VstatNorm)
                        endif
                    enddo
                enddo

                !Now, for block 5
                LinearSystem_h(VIndex:VStatIndex-1,VStatIndex:nLinearSystem) = nFCIHam(:,:)
                LinearSystem_h(VStatIndex:nLinearSystem,VIndex:VStatIndex-1) = nFCIHam(:,:)

                tempel = zzero 
                tempel2 = zzero
                do i = 1,CoreEnd
                    !Calc diagonal correction
                    tempel = tempel + StatickGF_Ann_Bra(i,1)*Ga_i_F_ij(i,1)
                enddo
                tempel4 = -tempel/sqrt(CStatNorm*CNorm)
                do a = VirtStart,VirtEnd
                    tempel2 = tempel2 + StatickGF_Cre_Bra(a,1)*Gc_b_F_ab(a,1)
                enddo
                tempel2 = tempel2/sqrt(VstatNorm*VNorm)
                !Add diagonal correction
                do i = VIndex,VStatIndex-1
                    LinearSystem_h(i+nFCIDet,i) = LinearSystem_h(i+nFCIDet,i) + tempel4
                    LinearSystem_h(i,i+nFCIDet) = LinearSystem_h(i,i+nFCIDet) + dconjg(tempel4)
                    LinearSystem_p(i+nFCIDet,i) = LinearSystem_p(i+nFCIDet,i) + tempel2
                    LinearSystem_p(i,i+nFCIDet) = LinearSystem_p(i,i+nFCIDet) + dconjg(tempel2)
                    Overlap_h(i+nFCIDet,i) = StatCoup_Ann/sqrt(CStatNorm*CNorm)
                    Overlap_h(i,i+nFCIDet) = dconjg(StatCoup_Ann)/sqrt(CStatNorm*CNorm)
                    Overlap_p(i+nFCIDet,i) = StatCoup_Cre/sqrt(VstatNorm*VNorm)
                    Overlap_p(i,i+nFCIDet) = dconjg(StatCoup_Cre)/sqrt(VstatNorm*VNorm)
                enddo

                call ApplyMomPert_GS(Psi_0,V0_Cre,V0_Ann,nImp_GF,nLinearSystem,nFCIDet,Coup_Ann_alpha,Coup_Create_alpha,    &
                    StatCoup_Ann,StatCoup_Cre,CNorm,CStatNorm,VNorm,VstatNorm,StatickGF_Cre_Emb_Ket,StatickGF_Ann_Emb_Ket) 

                if(tWriteOut) then
                    if(tNoCre) then
                        write(6,"(A)") "No virtual orbitals found at this kpoint"
                    elseif(tNoAnn) then
                        write(6,"(A)") "No occupied orbitals found at this kpoint"
                    endif
                endif

                if(tCheck) then
                    !Do a quick diagonalization of the overlap matrix, just to see
                    !how linearly dependent the basis is...
                    allocate(S_EigVec(nLinearSystem,nLinearSystem))
                    allocate(S_EigVal(nLinearSystem))
                    if(.not.tNoAnn) then
                        write(6,*) "Using Overlap_h"
                        S_EigVec(:,:) = Overlap_h(:,:)
                    else
                        write(6,*) "Using Overlap_p"
                        S_EigVec(:,:) = Overlap_p(:,:)
                    endif
                    !call writematrixcomp(S_EigVec,'S_EigVec',.true.)
                    allocate(Work(max(1,3*nLinearSystem-2)))
                    allocate(tempc(1))
                    lWork = -1
                    info = 0
                    call zheev('V','U',nLinearSystem,S_EigVec,nLinearSystem,S_EigVal,tempc,lWork,Work,info)
                    if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                    lwork = int(tempc(1))+1
                    deallocate(tempc)
                    allocate(tempc(lwork))
                    call zheev('V','U',nLinearSystem,S_EigVec,nLinearSystem,S_EigVal,tempc,lWork,Work,info)
                    if(info.ne.0) call stop_all(t_r,'S diag failed')
                    deallocate(work,tempc)
                    !Count linear dependencies
                    nSpan = 0
                    do i = 1,nLinearSystem
                        if(S_EigVal(i).gt.MinS_Eigval) then
                            nSpan = nSpan + 1
                        elseif(S_EigVal(i).lt.(-1.0e-8_dp)) then
                            write(6,*) "S_EigVal: ",i,S_EigVal(i)
                            call stop_all(t_r,'Overlap eigenvalues less than 0')
                        endif
                    enddo
                    if(nSpan.eq.nLinearSystem) then
                        write(6,"(A)") "No linear dependencies in basis?!"
                    else
                        write(6,"(A,I7,A,I7,A)") "Linearly dependent vectors in space: ",nLinearSystem-nSpan, &
                            " out of ",nFCIDet," possible."
                    endif
                    !call writevector(S_EigVal,'S_EigVal')
                    !I guess we could potentially rotate into this basis, but we don't really want to do that.
                    deallocate(S_EigVec,S_EigVal)
                endif

                !Construct final linear equations
                do i = 1,nLinearSystem
                    do j = 1,nLinearSystem
                        LinearSystem_h(j,i) = LinearSystem_h(j,i) + Overlap_h(j,i)*dcmplx(Omega+mu-GFChemPot,dDelta)
                    enddo
                enddo
                do i = 1,nLinearSystem
                    do j = 1,nLinearSystem
                        LinearSystem_p(j,i) = Overlap_p(j,i)*dcmplx(Omega+mu+GFChemPot,dDelta) - LinearSystem_p(j,i)
                    enddo
                enddo

                allocate(tempc(nLinearSystem))
                if(.not.tNoAnn) then
                    Psi1_h(:) = V0_Ann(:)
                    if(tNoStatickBasis) then
                        LinSize = VStatIndex-1
                        allocate(DynamicLS(LinSize,LinSize))
                        DynamicLS(:,:) = LinearSystem_h(1:LinSize,1:LinSize)
                        !call writematrixcomp(DynamicLS,'Hole Hamil',.true.)
                        !call writevectorcomp(Psi1_h(1:LinSize),'V0_Ann')
                        call SolveCompLinearSystem(DynamicLS,Psi1_h(1:LinSize),LinSize,info)
                        !call writevectorcomp(Psi1_h(1:LinSize),'|1>')
                        deallocate(DynamicLS)
                    else
                        LinSize = nLinearSystem
                        call SolveCompLinearSystem(LinearSystem_h,Psi1_h,nLinearSystem,info)
                    endif
                    if(info.ne.0) then 
                        write(6,*) "INFO: ",info
                        call warning(t_r,'Solving linear system failed for hole hamiltonian - skipping this frequency')
                        Omega = Omega + Omega_Step
                        cycle
                    endif
                    if(.not.tNoStatickBasis) then
                        !Apply S to |1>
                        call ZGEMV('N',nLinearSystem,nLinearSystem,zone,Overlap_h,nLinearSystem,Psi1_h,1,zzero,tempc,1)
                    else
                        tempc(1:LinSize) = Psi1_h(1:LinSize)
                    endif
                    !Norm is now dot product
                    dNorm_h = zzero
                    do i = 1,LinSize
                        dNorm_h = dNorm_h + dconjg(Psi1_h(i))*tempc(i)
                    enddo
                    !Now calculate spectrum
                    Response_h = zzero
                    do i = 1, LinSize
                        Response_h = Response_h + dconjg(V0_Ann(i))*tempc(i)
                    enddo
                else
                    Response_h = zzero
                    dNorm_h = zzero
                endif


                if(.not.tNoCre) then
                    !Now for particle
                    Psi1_p(:) = V0_Cre(:)
                    if(tNoStatickBasis) then
                        LinSize = VStatIndex-1
                        allocate(DynamicLS(LinSize,LinSize))
                        DynamicLS(:,:) = LinearSystem_p(1:LinSize,1:LinSize)
                        call SolveCompLinearSystem(DynamicLS,Psi1_p(1:LinSize),LinSize,info)
                        deallocate(DynamicLS)
                    else
                        LinSize = nLinearSystem
                        call SolveCompLinearSystem(LinearSystem_p,Psi1_p,nLinearSystem,info)
                    endif
                    if(info.ne.0) then 
                        write(6,*) "INFO: ",info
                        call warning(t_r,'Solving linear system failed for particle hamiltonian - skipping this frequency')
                        Omega = Omega + Omega_Step
                        cycle
                    endif

                    if(.not.tNoStatickBasis) then
                        call ZGEMV('N',nLinearSystem,nLinearSystem,zone,Overlap_p,nLinearSystem,Psi1_p,1,zzero,tempc,1)
                    else
                        tempc(1:LinSize) = Psi1_p(1:LinSize)
                    endif
                    !Norm is now dot product
                    dNorm_p = zzero
                    do i = 1,LinSize
                        dNorm_p = dNorm_p + dconjg(Psi1_p(i))*tempc(i)
                    enddo
                    !Now calculate spectrum
                    Response_p = zzero
                    do i = 1,LinSize
                        Response_p = Response_p + dconjg(V0_Cre(i))*tempc(i)
                    enddo
                else
                    Response_p = zzero
                    dNorm_p = zzero
                endif

                ResponseFn = Response_h + Response_p
                ni_lr = ni_lr_p + ni_lr_h
                deallocate(tempc)

                write(iunit,"(18G22.10,2I6)") Omega,real(ResponseFn),-aimag(ResponseFn),real(Response_p),   &
                    -aimag(Response_p),real(Response_h),-aimag(Response_h),HL_Energy,HL_Energy,abs(dNorm_p),    &
                    abs(dNorm_h),real(ni_lr),-aimag(ni_lr),real(ni_lr_p),   &
                    -aimag(ni_lr_p),real(ni_lr_h),-aimag(ni_lr_h),SpectralWeight,1,1
                call flush(iunit)

                Omega = Omega + Omega_Step
            enddo
            write(iunit,"(A)") ""
            write(iunit,"(A)") ""
        enddo

        deallocate(LinearSystem_h,Overlap_h,Overlap_p,LinearSystem_p)
        deallocate(NFCIHam,Np1FCIHam_alpha,Nm1FCIHam_beta)
        deallocate(Psi1_p,Psi1_h,FockSchmidt_SE,FockSchmidt_SE_VV,FockSchmidt_SE_CC,FockSchmidt_SE_VX)
        deallocate(FockSchmidt_SE_XV,FockSchmidt_SE_CX,FockSchmidt_SE_XC)
        deallocate(Gc_a_F_ax_Bra,Gc_a_F_ax_Ket,Gc_b_F_ab,Ga_i_F_xi_Bra,Ga_i_F_xi_Ket,Ga_i_F_ij)
        deallocate(Xc_a_F_ax_Bra,Xc_a_F_ax_Ket,Xc_b_F_ab,Xa_i_F_xi_Bra,Xa_i_F_xi_Ket,Xa_i_F_ij)
        deallocate(SchmidtkGF_Cre_Ket,SchmidtkGF_Ann_Ket,SchmidtkGF_Cre_Bra,SchmidtkGF_Ann_Bra)
        deallocate(StatickGF_Cre_Ket,StatickGF_Ann_Ket,StatickGF_Cre_Bra,StatickGF_Ann_Bra)
        deallocate(StatickGF_Ann_Emb_Ket,StatickGF_Ann_Emb_Bra,StatickGF_Cre_Emb_Ket,StatickGF_Cre_Emb_Bra)
        deallocate(Psi_0,V0_Ann,V0_Cre)
        deallocate(Coup_Create_alpha,Coup_Ann_alpha)
        !Deallocate determinant lists
        if(allocated(FCIDetList)) deallocate(FCIDetList)
        if(allocated(FCIBitList)) deallocate(FCIBitList)
        if(allocated(UMat)) deallocate(UMat)
        if(allocated(TMat)) deallocate(TMat)
        if(allocated(TMat_Comp)) deallocate(TMat_Comp)
        if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
        if(allocated(Nm1BitList)) deallocate(Nm1BitList)
        if(allocated(Np1FCIDetList)) deallocate(Np1FCIDetList)
        if(allocated(Np1BitList)) deallocate(Np1BitList)
        if(allocated(Nm1bFCIDetList)) deallocate(Nm1bFCIDetList)
        if(allocated(Nm1bBitList)) deallocate(Nm1bBitList)
        if(allocated(Np1bFCIDetList)) deallocate(Np1bFCIDetList)
        if(allocated(Np1bBitList)) deallocate(Np1bBitList)
        close(iunit)

    end subroutine MomGF_Ex

    !Apply the schmidt decomposed momentum perturbation to the ground state wavefunction.
    !Since this perturbation is non-local, it will excite to all components of the linear
    !response basis
    !If tLR_ReoptGS=F then Psi_0 will be of size nFCIDet, otherwise it must be larger.
    subroutine ApplyMomPert_GS(Psi_0,V0_Cre,V0_Ann,nImp_GF,nSizeLR,nSizeGS,Coup_Ann_alpha,Coup_Create_alpha,    &
        StatCoup_Ann,StatCoup_Cre,CNorm,CStatNorm,VNorm,VstatNorm,StatickGF_Cre_Emb_Ket,StatickGF_Ann_Emb_Ket)
        use DetToolsData
        implicit none
        integer, intent(in) :: nSizeLR,nSizeGS,nImp_GF
        complex(dp), intent(in) :: Psi_0(nSizeGS)
        complex(dp), intent(in) :: StatCoup_Ann,StatCoup_Cre,CNorm,CStatNorm,VNorm,VstatNorm
        complex(dp), intent(out) :: V0_Cre(nSizeLR),V0_Ann(nSizeLR)
        integer, intent(in) :: Coup_Ann_alpha(2,nFCIDet,nNp1FCIDet),Coup_Create_alpha(2,nFCIDet,nNm1bFCIDet)
        complex(dp), intent(in) :: StatickGF_Ann_Emb_Ket(nOcc-nImp+1:nOcc+nImp,nImp_GF)
        complex(dp), intent(in) :: StatickGF_Cre_Emb_Ket(nOcc-nImp+1:nOcc+nImp,nImp_GF)
        integer :: j,k
        character(len=*), parameter :: t_r='ApplyMomPert_GS'

        if(.not.tLR_ReoptGS) then
            if(nSizeGS.ne.nFCIDet) call stop_all(t_r,'Error in GS size')
        else
            call stop_all(t_r,'Not ready for reoptimizing GS yet')
        endif

        V0_Cre(:) = zzero
        V0_Ann(:) = zzero

        !call writevectorcomp(Psi_0,'Psi_0')
        !call writevectorcomp(StatickGF_Ann_Emb_Ket(:,1),'Statick_Ann_Emb')

        !Projection onto Block 1 of V0
        do K = 1,nFCIDet
            do J = 1,nNm1bFCIDet
                if(Coup_Create_alpha(1,K,J).ne.0) then
                    !Determinants are connected by a single SQ operator
                    V0_Ann(J) = V0_Ann(J) + StatickGF_Ann_Emb_Ket(Coup_Create_alpha(1,K,J),1)*Psi_0(K)* &
                        Coup_Create_alpha(2,K,J)
                endif
            enddo
        enddo

        !Projection onto Block 2 of V0_Ann
        do K = 1,nFCIDet
            V0_Ann(K+nNm1bFCIDet) = (dconjg(StatCoup_Ann)*Psi_0(K))/sqrt(CNorm)
        enddo
        !Projection onto Block 3 of V0_Ann
        do K = 1,nFCIDet
            V0_Ann(K+nNm1bFCIDet+nFCIDet) = sqrt(CStatNorm)*Psi_0(K)
        enddo

        !call writevectorcomp(V0_Ann,'V0_Ann')

        !Creation operator
        do K = 1,nFCIDet
            do J = 1,nNp1FCIDet
                if(Coup_Ann_alpha(1,K,J).ne.0) then
                    !Connected
                    V0_Cre(J) = V0_Cre(J) + StatickGF_Cre_Emb_Ket(Coup_Ann_alpha(1,K,J),1)*Psi_0(K)* &
                        Coup_Ann_alpha(2,K,J)
                endif
            enddo
        enddo
        !Block 2
        do K = 1,nFCIDet
            V0_Cre(K+nNp1FCIDet) = (dconjg(StatCoup_Cre)*Psi_0(K))/sqrt(VNorm)
        enddo
        !Block 3
        do K = 1,nFCIDet
            V0_Cre(K+nNp1FCIDet+nFCIDet) = sqrt(VstatNorm)*Psi_0(K)
        enddo

    end subroutine ApplyMomPert_GS

    !Construct the contraction coefficients for the action of the static perturbation 
    subroutine FindStaticMomSchmidtPert(kPnt,nImp_GF,StatickGF_Cre_Ket,StatickGF_Cre_Emb_Ket,   &
        StatickGF_Ann_Ket,StatickGF_Ann_Emb_Ket,StatickGF_Cre_Bra,   &
        StatickGF_Ann_Bra,StatickGF_Cre_Emb_Bra,StatickGF_Ann_Emb_Bra,tNoCre,tNoAnn)
        implicit none
        integer, intent(in) :: kPnt,nImp_GF
        logical, intent(out) :: tNoCre,tNoAnn
        complex(dp), intent(out) :: StatickGF_Cre_Ket(nOcc+nImp+1:nSites,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Ann_Ket(1:nOcc-nImp,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Cre_Bra(nOcc+nImp+1:nSites,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Ann_Bra(1:nOcc-nImp,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Cre_Emb_Ket(nOcc-nImp+1:nOcc+nImp,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Cre_Emb_Bra(nOcc-nImp+1:nOcc+nImp,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Ann_Emb_Ket(nOcc-nImp+1:nOcc+nImp,nImp_GF)
        complex(dp), intent(out) :: StatickGF_Ann_Emb_Bra(nOcc-nImp+1:nOcc+nImp,nImp_GF)
        complex(dp), allocatable :: HFPertBasis_Ann(:),HFPertBasis_Cre(:)
        integer :: SS_Period,ind_1,ind_2,i,n
        logical :: tCreFound,tAnnFound
        character(len=*), parameter :: t_r='FindStaticMomSchmidtPert'

        tCreFound = .false.
        tAnnFound = .false.
        
        SS_Period = nImp
        allocate(HFPertBasis_Ann(SS_Period))
        allocate(HFPertBasis_Cre(SS_Period))
        HFPertBasis_Ann(:) = zzero
        HFPertBasis_Cre(:) = zzero
            
        ind_1 = ((kPnt-1)*SS_Period) + 1
        ind_2 = SS_Period * kPnt

        do i = 1,SS_Period  !Run over vectors which span this kpoint
            if(KVec_InvEMap(ind_1+i-1).le.nOcc) then
                write(6,*) "For this kpoint, occupied orbital for vector: ",i
                !i is an occupied orbital
                !do x = 1,SS_Period  !Run over components of this vector
                !    HFPertBasis_Ann(i) = HFPertBasis_Ann(i) + dconjg(k_vecs(x,ind_1 + i - 1))
                !enddo
                HFPertBasis_Ann(i) = zone
                tAnnFound = .true.
            else
                write(6,*) "For this kpoint, virtual orbital for vector: ",i
                !i is a virtual orbital
                !do x = 1,SS_Period  !Run over components of this vector
                !    HFPertBasis_Cre(i) = HFPertBasis_Cre(i) + dconjg(k_vecs(x,ind_1 + i - 1))
                !enddo
                HFPertBasis_Cre(i) = zone
                tCreFound = .true.
            endif 
        enddo

        !call writevectorcomp(HFPertBasis_Cre,'HFPertBasis_Cre')
        !call writevectorcomp(HFPertBasis_Ann,'HFPertBasis_Ann')
        
        !Now rotate from the HF vectors, to the Schmidt basis
        StatickGF_Ann_Ket(:,:) = zzero
        do n = 1,nOcc-nImp  !Loop over schmidt basis of core orbitals
            do i = 1,SS_Period  !Contract over eigenvectors corresponding to this kpoint
                StatickGF_Ann_Ket(n,1) = StatickGF_Ann_Ket(n,1) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Ann(i)
            enddo
        enddo
        StatickGF_Cre_Ket(:,:) = zzero
        do n = nOcc+nImp+1,nSites
            do i = 1,SS_Period
                StatickGF_Cre_Ket(n,1) = StatickGF_Cre_Ket(n,1) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Cre(i)
            enddo
        enddo
        StatickGF_Cre_Emb_Ket(:,:) = zzero
        StatickGF_Ann_Emb_Ket(:,:) = zzero
        do n = nOcc-nImp+1,nOcc+nImp
            do i = 1,SS_Period
                StatickGF_Cre_Emb_Ket(n,1) = StatickGF_Cre_Emb_Ket(n,1) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Cre(i)
                StatickGF_Ann_Emb_Ket(n,1) = StatickGF_Ann_Emb_Ket(n,1) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Ann(i)
            enddo
        enddo

        StatickGF_Ann_Bra(:,1) = dconjg(StatickGF_Ann_Ket(:,1))
        StatickGF_Cre_Bra(:,1) = dconjg(StatickGF_Cre_Ket(:,1))
        StatickGF_Cre_Emb_Bra(:,1) = dconjg(StatickGF_Cre_Emb_Ket(:,1))
        StatickGF_Ann_Emb_Bra(:,1) = dconjg(StatickGF_Ann_Emb_Ket(:,1))

        !call writevectorcomp(StatickGF_Ann_Ket(:,1),'StatickGF_Ann_Ket')
        !call writevectorcomp(StatickGF_Cre_Ket(:,1),'StatickGF_Cre_Ket')

        deallocate(HFPertBasis_Ann,HFPertBasis_Cre)

        tNoCre = .not.tCreFound
        tNoAnn = .not.tAnnFound

        if(tNoCre.and.tNoAnn) call stop_all(t_r,'Cannot have no virtual or occupied bands at a given kpoint')

    end subroutine FindStaticMomSchmidtPert

    subroutine FindMomGFSchmidtPert(Omega,ni_lr_cre,ni_lr_ann,kPnt,nImp_GF,SchmidtkGF_Cre_Ket,SchmidtkGF_Ann_Ket,   &
        SchmidtkGF_Cre_Bra,SchmidtkGF_Ann_Bra)
        implicit none
        integer, intent(in) :: kPnt !This is the index of the kpoint perturbation
        integer, intent(in) :: nImp_GF
        real(dp), intent(in) :: Omega
        complex(dp), intent(out) :: ni_lr_cre,ni_lr_ann
        complex(dp), intent(out) :: SchmidtkGF_Cre_Ket(nOcc+nImp+1:nSites,nImp_GF),SchmidtkGF_Ann_Ket(1:nOcc-nImp,nImp_GF)
        complex(dp), intent(out) :: SchmidtkGF_Cre_Bra(nOcc+nImp+1:nSites,nImp_GF),SchmidtkGF_Ann_Bra(1:nOcc-nImp,nImp_GF)
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
            if(KVec_InvEMap(ind_1+i-1).le.nOcc) then
                !i is an occupied orbital
!                write(6,*) "Orbital is occupied",i,ind_1+i-1
                do x = 1,SS_Period  !Run over components of this vector
!                    HFPertBasis_Ann(i) = HFPertBasis_Ann(i) + dconjg(k_vecs(x,ind_1 + i - 1)) /     &
!                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                    ni_lr_ann = ni_lr_ann + (k_vecs(x,ind_1 + i - 1)*dconjg(k_vecs(x,ind_1 + i - 1))) / &
                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                enddo
                HFPertBasis_Ann(i) = zone / (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
            else
                !i is a virtual orbital
!                write(6,*) "Orbital is virtual",i,ind_1+i-1
                do x = 1,SS_Period  !Run over components of this vector
!                    HFPertBasis_Cre(i) = HFPertBasis_Cre(i) + dconjg(k_vecs(x,ind_1 + i - 1)) /     &
!                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                    ni_lr_cre = ni_lr_cre + (k_vecs(x,ind_1 + i - 1)*dconjg(k_vecs(x,ind_1 + i - 1))) / &
                        (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
                enddo
                HFPertBasis_Cre(i) = zone / (dcmplx(Omega,dDelta)-k_HFEnergies(ind_1 + i - 1))
            endif
        enddo

        !call writevectorcomp(HFPertBasis_Cre,'HFPertBasis_Cre')
        !call writevectorcomp(HFPertBasis_Ann,'HFPertBasis_Ann')

        !Now rotate from the HF vectors, to the Schmidt basis
        SchmidtkGF_Ann_Ket(:,:) = zzero
        do n = 1,nOcc-nImp  !Loop over schmidt basis of core orbitals
            do i = 1,SS_Period  !Contract over eigenvectors corresponding to this kpoint
                SchmidtkGF_Ann_Ket(n,1) = SchmidtkGF_Ann_Ket(n,1) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Ann(i)
            enddo
        enddo
        SchmidtkGF_Cre_Ket(:,:) = zzero
        do n = nOcc+nImp+1,nSites
            do i = 1,SS_Period
                SchmidtkGF_Cre_Ket(n,1) = SchmidtkGF_Cre_Ket(n,1) + k_HFtoSchmidtTransform(n,ind_1+i-1)*HFPertBasis_Cre(i)
            enddo
        enddo

        SchmidtkGF_Ann_Bra(:,:) = dconjg(SchmidtkGF_Ann_Ket(:,:))
        SchmidtkGF_Cre_Bra(:,:) = dconjg(SchmidtkGF_Cre_Ket(:,:))
        
        !call writevectorcomp(SchmidtkGF_Ann_Ket(:,1),'SchmidtkGF_Ann_Ket')

        deallocate(HFPertBasis_Ann,HFPertBasis_Cre)
        
    end subroutine FindMomGFSchmidtPert

end module MomSpectra
