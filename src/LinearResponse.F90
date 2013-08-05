module LinearResponse
    use const
    use timing
    use errors, only: stop_all,warning
    use mat_tools, only: WriteVector,WriteMatrix,WriteVectorInt,WriteMatrixComp,WriteVectorComp,znrm2
    use globals
    use LRSolvers
    implicit none
    integer :: CVIndex,AVIndex,CAIndex
    
    contains
    
    subroutine NonIntExCont_TDA_MCLR_DD_Cmprs()
        use utils, only: get_free_unit,append_ext_real,append_ext
        use DetBitOps, only: DecodeBitDet,SQOperator,CountBits
        use zminresqlpModule, only: MinresQLP  
        use DetToolsData
        use DetTools, only: tospat,umatind,GenDets,GetHElement
        use solvers, only: CountSizeCompMat,StoreCompMat,CreateIntMats
        implicit none
        integer :: a,i,j,k,b,beta
        integer :: CoreEnd,VirtStart,VirtEnd,iunit
        integer :: DiffOrb,nOrbs,gam,gam1_ind,gam1_spat,gam2,gam2_ind,gam2_spat,ierr
        integer :: gam_spat,nLinearSystem,tempK,ActiveEnd,ActiveStart
        integer :: orbdum(1),info,nCore,nVirt,nFullNm1,nFullNp1,iters
        integer :: maxminres_iter,minres_unit,InterMem,HamMem,CoupMem,gam_ind,i_Ann,i_Create,i_Ann_b,i_Create_b
        integer :: ind,ind_S,Nmax_Coup,Nmax_Lin,Nmax_Nm1,Nmax_Nm1b,Nmax_Np1,Nmax_Np1b,Nmax_N,Nmax_S,pertsite
        integer :: pertsitealpha,pertsitebeta,row
        logical :: tParity
        real(dp) :: Omega,CVNorm
        complex(dp) :: tempel,ResponseFn,dOrthog,dNorm,testc,ni_lr,matel,matel_S
        complex(dp) , allocatable :: Overlap(:)
        complex(dp) , allocatable :: temp_vecc(:),RHS(:),RHS2(:),Psi1(:)
        complex(dp) , allocatable :: G_ai_G_aj(:,:),G_ai_G_bi(:,:),G_xa_G_ya(:,:)
        complex(dp) , allocatable :: G_xa_F_ab(:,:),G_xa_G_yb_F_ab(:,:),FockSchmidtComp(:,:),G_ia_G_xa(:,:)
        complex(dp) , allocatable :: F_xi_G_ia_G_ya(:,:),G_xa_F_ya(:,:),G_xi_G_yi(:,:),G_ix_F_ij(:,:)
        complex(dp) , allocatable :: G_ix_G_jy_F_ji(:,:),G_ia_G_ix(:,:),F_ax_G_ia_G_iy(:,:),G_ix_F_iy(:,:)
        complex(dp) , allocatable , target :: LinearSystem(:)
        real(dp), allocatable :: NFCIHam_cmps(:),Nm1FCIHam_alpha_cmps(:),Nm1FCIHam_beta_cmps(:)
        real(dp), allocatable :: Np1FCIHam_alpha_cmps(:),Np1FCIHam_beta_cmps(:)
        integer, allocatable :: NFCIHam_inds(:),Np1FCIHam_alpha_inds(:),Np1FCIHam_beta_inds(:)
        integer, allocatable :: Nm1FCIHam_alpha_inds(:),Nm1FCIHam_beta_inds(:)
        real(dp), allocatable :: AVNorm(:),CANorm(:)
        integer, allocatable :: Coup_Create_alpha(:,:),Coup_Create_alpha_T(:,:),Coup_Create_beta(:,:)
        integer, allocatable :: Coup_Create_beta_T(:,:),Coup_Ann_alpha(:,:),Coup_Ann_alpha_T(:,:)
        integer, allocatable :: Coup_Ann_beta(:,:),Coup_Ann_beta_T(:,:)
        integer, allocatable :: Coup_Create_alpha_inds(:),Coup_Create_alpha_inds_T(:),Coup_Create_beta_inds(:) 
        integer, allocatable :: Coup_Create_beta_inds_T(:),Coup_Ann_alpha_inds(:),Coup_Ann_alpha_inds_T(:)
        integer, allocatable :: Coup_Ann_beta_inds(:),Coup_Ann_beta_inds_T(:)
        integer, allocatable :: Coup_Create_alpha_cum(:),Coup_Create_alpha_cum_T(:),Coup_Create_beta_cum(:) 
        integer, allocatable :: Coup_Create_beta_cum_T(:),Coup_Ann_alpha_cum(:),Coup_Ann_alpha_cum_T(:)
        integer, allocatable :: Coup_Ann_beta_cum(:),Coup_Ann_beta_cum_T(:)
        integer, allocatable :: Overlap_inds(:)
        integer, allocatable, target :: LinearSystem_inds(:)
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='NonIntExCont_TDA_MCLR_DD_Cmprs'
        integer(ip) :: maxminres_iter_ip,minres_unit_ip,nLinearSystem_ip,info_ip

        if(tUHF) then
            !What needs to be done for this:
            !   Construct G(w) for aa and bb excitations, and use appropriately
            !   Different fock sectors for aa and bb
            call stop_all(t_r,'UHF not yet correctly implemented for DD response - bug ghb')
        endif

        call set_timer(LR_EC_TDA_Precom)
        
        maxminres_iter = iMinRes_MaxIter
        iters = 0
        zDirMV_Mat => null()

        write(6,*) "Calculating non-interacting EC MR-TDA LR system for DD excitations..."
        write(6,*) "Perturbation is number operator for spatial orbitals"
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(tMinRes_NonDir) then
            if(tPreCond_MinRes) then
                write(6,"(A)") "Solving linear system with iterative non-direct preconditioned MinRes-QLP algorithm"
            else
                write(6,"(A)") "Solving linear system with iterative non-direct MinRes-QLP algorithm"
            endif
            write(6,"(A,G22.10)") "Tolerance for solution of linear system: ",rtol_LR
            write(6,"(A,G22.10)") "Maximum iterations for each solution: ",maxminres_iter
        elseif(tGMRES_NonDir) then
            if(tPreCond_MinRes) then
                write(6,'(A)') "Solving linear system with preconditioned iterative non-direct GMRES algorithm"
            else
                write(6,'(A)') "Solving linear system with iterative non-direct GMRES algorithm"
            endif
            write(6,"(A,G22.10)") "Tolerance for solution of linear system: ",rtol_LR
            write(6,"(A,G22.10)") "Maximum iterations for each solution: ",maxminres_iter
            write(6,"(A,I8)") "Number of krylov subspace vectors to store: ",nKrylov
        endif
        !umat and tmat for the active space
        call CreateIntMats(tComp=.false.)
        if(tMinRes_NonDir) then
            minres_unit = get_free_unit()
            open(minres_unit,file='zMinResQLP.txt',status='unknown')
        else
            minres_unit = get_free_unit()
            open(minres_unit,file='zGMRES.txt',status='unknown')
        endif
        
        !Enumerate excitations for fully coupled space
        !Seperate the lists into different Ms sectors in the N+- lists
        call GenDets(Elec,EmbSize,.true.,.true.,.true.) 
        write(6,*) "Number of determinants in {N,N+1,N-1} FCI space: ",ECoupledSpace
        write(6,"(A,F12.5,A)") "Memory required for det list storage: ",DetListStorage*RealtoMb, " Mb"

        !Construct FCI hamiltonians for the N, N-1_alpha, N-1_beta, N+1_alpha and N+1_beta spaces
        !N electron
        if(nNp1FCIDet.ne.nNm1FCIDet) call stop_all(t_r,'Active space not half-filled')
        if(nNp1FCIDet.ne.nNp1bFCIDet) call stop_all(t_r,'Cannot deal with open shell systems')
        write(6,"(A)") "Computing size of compressed Hamiltonian matrices"
        call flush(6)
        if(iHamSize_N.eq.0) then
            call CountSizeCompMat(FCIDetList,Elec,nFCIDet,Nmax_N,FCIBitList)
            iHamSize_N = Nmax_N
        else
            write(6,"(A,I16)") "Assuming size of compressed N electron hamiltonian is: ",iHamSize_N
            Nmax_N = iHamSize_N
        endif
        if(iHamSize_Nm1.eq.0) then
            call CountSizeCompMat(Nm1bFCIDetList,Elec-1,nNm1bFCIDet,Nmax_Nm1b,Nm1bBitList)
            iHamSize_Nm1 = Nmax_Nm1b
        else
            write(6,"(A,I16)") "Assuming size of compressed N-1 electron hamiltonian is: ",iHamSize_Nm1
            Nmax_Nm1b = iHamSize_Nm1
        endif
        if(iHamSize_Np1.eq.0) then
            call CountSizeCompMat(Np1FCIDetList,Elec+1,nNp1FCIDet,Nmax_Np1,Np1BitList)
            iHamSize_Np1 = Nmax_Np1
        else
            write(6,"(A,I16)") "Assuming size of compressed N+1 electron hamiltonian is: ",iHamSize_Np1
            Nmax_Np1  = iHamSize_Np1
        endif

        !Assume flipping spins for the open shells doesn't change the number of off diagonal matrix elements
        Nmax_Nm1 = Nmax_Nm1b
        Nmax_Np1b = Nmax_Np1
        write(6,"(A,F12.5,A)") "Memory required for N-electron hamil: ",(real(Nmax_N,dp)*2)*RealtoMb, " Mb"
        write(6,"(A,F12.5,A)") "Memory required for N+1-electron hamil: ",(real(Nmax_Nm1b,dp)*4)*RealtoMb, " Mb"
        write(6,"(A,F12.5,A)") "Memory required for N-1-electron hamil: ",(real(Nmax_Np1,dp)*4)*RealtoMb, " Mb"
        HamMem = 2*nNp1FCIDet**2 + 2*nNm1FCIDet**2 + nFCIDet**2
        write(6,'(A,F9.3,A)') "IF uncompressed, memory required for hamiltonians: ",real(HamMem,dp)*RealtoMb,' Mb'
        call flush(6)
        allocate(NFCIHam_cmps(Nmax_N),stat=ierr)
        allocate(Np1FCIHam_alpha_cmps(Nmax_Np1),stat=ierr)
        allocate(Np1FCIHam_beta_cmps(Nmax_Np1b),stat=ierr)
        allocate(Nm1FCIHam_alpha_cmps(Nmax_Nm1),stat=ierr)
        allocate(Nm1FCIHam_beta_cmps(Nmax_Nm1b),stat=ierr)
        !And for indexing arrays
        allocate(NFCIHam_inds(Nmax_N),stat=ierr)
        allocate(Np1FCIHam_alpha_inds(Nmax_Np1),stat=ierr)
        allocate(Np1FCIHam_beta_inds(Nmax_Np1b),stat=ierr)
        allocate(Nm1FCIHam_alpha_inds(Nmax_Nm1),stat=ierr)
        allocate(Nm1FCIHam_beta_inds(Nmax_Nm1b),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Allocation fail')

        NFCIHam_cmps(:) = zero
        Np1FCIHam_alpha_cmps(:) = zero
        Np1FCIHam_beta_cmps(:) = zero
        Nm1FCIHam_alpha_cmps(:) = zero
        Nm1FCIHam_beta_cmps(:) = zero
        NFCIHam_inds(:) = 0
        Np1FCIHam_alpha_inds(:) = 0
        Np1FCIHam_beta_inds(:) = 0
        Nm1FCIHam_alpha_inds(:) = 0
        Nm1FCIHam_beta_inds(:) = 0

        call StoreCompMat(FCIDetList,Elec,nFCIDet,Nmax_N,NFCIHam_cmps,NFCIHam_inds,FCIBitList)
        call StoreCompMat(Np1FCIDetList,Elec+1,nNp1FCIDet,Nmax_Np1,Np1FCIHam_alpha_cmps,Np1FCIHam_alpha_inds,Np1BitList)
        call StoreCompMat(Np1bFCIDetList,Elec+1,nNp1bFCIDet,Nmax_Np1b,Np1FCIHam_beta_cmps,Np1FCIHam_beta_inds,Np1bBitList)
        call StoreCompMat(Nm1FCIDetList,Elec-1,nNm1FCIDet,Nmax_Nm1,Nm1FCIHam_alpha_cmps,Nm1FCIHam_alpha_inds,Nm1BitList)
        call StoreCompMat(Nm1bFCIDetList,Elec-1,nNm1bFCIDet,Nmax_Nm1b,Nm1FCIHam_beta_cmps,Nm1FCIHam_beta_inds,Nm1bBitList)
        
        nLinearSystem = (2*nFCIDet) + nImp*2*(nNm1FCIDet+nNm1bFCIDet) + nImp*2*(nNp1FCIDet+nNp1bFCIDet) 
        write(6,*) "Total size of linear sys: ",nLinearSystem
        
        !Set up orbital indices
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
        !And full N pm 1 space sizes
        nFullNm1 = nNm1FCIDet + nNm1bFCIDet
        nFullNp1 = nNp1FCIDet + nNp1bFCIDet

        !Set up indices for the block of the linear system
        CVIndex = nFCIDet + 1   !Beginning of EC Core virtual excitations
        AVIndex = nFCIDet + nFCIDet + 1 !Beginning of EC Active Virtual excitations
        CAIndex = AVIndex + (nImp*2)*nFullNm1 !Beginning of EC Core Active excitations

        write(6,*) "CV indices start from: ",CVIndex
        write(6,*) "AV indices start from: ",AVIndex
        write(6,*) "CA indices start from: ",CAIndex
            
        !Allocate memory for normalization constants
        allocate(AVNorm(1:nImp*2))
        allocate(CANorm(1:nImp*2))

        !Allocate and precompute 1-operator coupling coefficients between the different sized spaces.
        !First number is the index of operator, and the second is the parity change when applying the operator

        !First, determine the size of these matrices in compressed form
        !Assume that all *8* matrices are the same size (this should be true for hole-particle symmetric systems)
        Nmax_Coup = 0
        do J = 1,nFCIDet
            do K = 1,nNm1bFCIDet
                DiffOrb = ieor(Nm1bBitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(nOrbs.eq.1) then
                    Nmax_Coup = Nmax_Coup + 1
                endif
            enddo
        enddo

        write(6,'(A,I8)') "Number of elements in coupling matrices: ",Nmax_Coup
        !Need to store 24 of this number - 3 for each coupling matrix (2 integers per list + indexing) and 8 different types (including transposes)
        write(6,'(A,F9.3,A)') "Memory required to store compressed coupling matrices: ",24*real(Nmax_Coup,dp)*RealtoMb,' Mb'
        CoupMem = 2*nFCIDet*nNm1bFCIDet + 2*nFCIDet*nNm1FCIDet + 2*nFCIDet*nNp1FCIDet + 2*nFCIDet*nNp1bFCIDet
        write(6,'(A,F9.3,A)') "Memory required to store uncompressed coupling matrices: ",real(CoupMem,dp)*RealtoMb,' Mb'

        allocate(Coup_Create_alpha(2,Nmax_Coup))
        allocate(Coup_Create_alpha_inds(Nmax_Coup))
        allocate(Coup_Create_alpha_cum(nFCIDet+1))  !Indexing positions of new rows
        Coup_Create_alpha(:,:) = 0
        Coup_Create_alpha_inds(:) = 0
        Coup_Create_alpha_cum(:) = 0
        Coup_Create_alpha_cum(1) = 1

        allocate(Coup_Create_alpha_T(2,Nmax_Coup))
        allocate(Coup_Create_alpha_inds_T(Nmax_Coup))
        allocate(Coup_Create_alpha_cum_T(nNm1bFCIDet+1))  !Indexing positions of new rows
        Coup_Create_alpha_T(:,:) = 0
        Coup_Create_alpha_inds_T(:) = 0
        Coup_Create_alpha_cum_T(:) = 0
        Coup_Create_alpha_cum_T(1) = 1

        allocate(Coup_Create_beta(2,Nmax_Coup))
        allocate(Coup_Create_beta_inds(Nmax_Coup))
        allocate(Coup_Create_beta_cum(nFCIDet+1))  !Indexing positions of new rows
        Coup_Create_beta(:,:) = 0
        Coup_Create_beta_inds(:) = 0
        Coup_Create_beta_cum(:) = 0
        Coup_Create_beta_cum(1) = 1
        
        allocate(Coup_Create_beta_T(2,Nmax_Coup))
        allocate(Coup_Create_beta_inds_T(Nmax_Coup))
        allocate(Coup_Create_beta_cum_T(nNm1FCIDet+1))  !Indexing positions of new rows
        Coup_Create_beta_T(:,:) = 0
        Coup_Create_beta_inds_T(:) = 0
        Coup_Create_beta_cum_T(:) = 0
        Coup_Create_beta_cum_T(1) = 1
        
        allocate(Coup_Ann_alpha(2,Nmax_Coup))
        allocate(Coup_Ann_alpha_inds(Nmax_Coup))
        allocate(Coup_Ann_alpha_cum(nFCIDet+1))  !Indexing positions of new rows
        Coup_Ann_alpha(:,:) = 0
        Coup_Ann_alpha_inds(:) = 0
        Coup_Ann_alpha_cum(:) = 0
        Coup_Ann_alpha_cum(1) = 1
        
        allocate(Coup_Ann_alpha_T(2,Nmax_Coup))
        allocate(Coup_Ann_alpha_inds_T(Nmax_Coup))
        allocate(Coup_Ann_alpha_cum_T(nNp1FCIDet+1))  !Indexing positions of new rows
        Coup_Ann_alpha_T(:,:) = 0
        Coup_Ann_alpha_inds_T(:) = 0
        Coup_Ann_alpha_cum_T(:) = 0
        Coup_Ann_alpha_cum_T(1) = 1
        
        allocate(Coup_Ann_beta(2,Nmax_Coup))
        allocate(Coup_Ann_beta_inds(Nmax_Coup))
        allocate(Coup_Ann_beta_cum(nFCIDet+1))  !Indexing positions of new rows
        Coup_Ann_beta(:,:) = 0
        Coup_Ann_beta_inds(:) = 0
        Coup_Ann_beta_cum(:) = 0
        Coup_Ann_beta_cum(1) = 1
        
        allocate(Coup_Ann_beta_T(2,Nmax_Coup))
        allocate(Coup_Ann_beta_inds_T(Nmax_Coup))
        allocate(Coup_Ann_beta_cum_T(nNp1bFCIDet+1))  !Indexing positions of new rows
        Coup_Ann_beta_T(:,:) = 0
        Coup_Ann_beta_inds_T(:) = 0
        Coup_Ann_beta_cum_T(:) = 0
        Coup_Ann_beta_cum_T(1) = 1

        i_Create = 0
        i_Create_b = 0
        i_Ann = 0
        i_Ann_b = 0
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

                i_Create = i_Create + 1
                if(i_Create.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Create_alpha(1,i_Create) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_alpha(2,i_Create) = -1
                else
                    Coup_Create_alpha(2,i_Create) = 1
                endif
                Coup_Create_alpha_inds(i_Create) = K
            enddo
            !Now the creation of a beta orbital
            do K = 1,nNm1FCIDet
                DiffOrb = ieor(Nm1BitList(K),FCIBitList(J))
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
                if(mod(gam,2).ne.0) then
                    call stop_all(t_r,'differing orbital should be an beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = Nm1BitList(K)
                call SQOperator(tempK,gam,tParity,.false.)

                i_Create_b = i_Create_b + 1
                if(i_Create_b.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Create_beta(1,i_Create_b) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_beta(2,i_Create_b) = -1
                else
                    Coup_Create_beta(2,i_Create_b) = 1
                endif
                Coup_Create_beta_inds(i_Create_b) = K
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

                i_Ann = i_Ann + 1
                if(i_Ann.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Ann_alpha(1,i_Ann) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_alpha(2,i_Ann) = -1
                else
                    Coup_Ann_alpha(2,i_Ann) = 1
                endif
                Coup_Ann_alpha_inds(i_Ann) = K
            enddo
            !Now, annihilation of a beta orbital
            do K = 1,nNp1bFCIDet
                DiffOrb = ieor(Np1bBitList(K),FCIBitList(J))
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
                if(mod(gam,2).ne.0) then
                    call stop_all(t_r,'differing orbital should be an beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = Np1bBitList(K)
                call SQOperator(tempK,gam,tParity,.true.)

                i_Ann_b = i_Ann_b + 1
                if(i_Ann_b.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Ann_beta(1,i_Ann_b) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_beta(2,i_Ann_b) = -1
                else
                    Coup_Ann_beta(2,i_Ann_b) = 1
                endif
                Coup_Ann_beta_inds(i_Ann_b) = K
            enddo

            Coup_Create_alpha_cum(J+1) = i_Create + 1
            Coup_Create_beta_cum(J+1) = i_Create_b + 1
            Coup_Ann_alpha_cum(J+1) = i_Ann + 1
            Coup_Ann_beta_cum(J+1) = i_Ann_b + 1
            if(mod(J,25000).eq.0) then
                write(6,*) "Constructing coupling matrices 1: ",J,nFCIDet
                call flush(6)
            endif
        enddo

        !Now to explicitly construct the transpose of the matrices too!
        i = 0
        do J = 1,nNm1bFCIDet
            !Start with the creation of an alpha orbital
            do K = 1,nFCIDet
                DiffOrb = ieor(Nm1bBitList(J),FCIBitList(K))
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
                tempK = FCIBitList(K)
                call SQOperator(tempK,gam,tParity,.true.)

                i = i + 1
                if(i.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Create_alpha_T(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_alpha_T(2,i) = -1
                else
                    Coup_Create_alpha_T(2,i) = 1
                endif
                Coup_Create_alpha_inds_T(i) = K
            enddo
            Coup_Create_alpha_cum_T(J+1) = i+1
        enddo
        
        i = 0
        do J = 1,nNm1FCIDet
            !Now the creation of an beta orbital
            do K = 1,nFCIDet
                DiffOrb = ieor(Nm1BitList(J),FCIBitList(K))
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
                if(mod(gam,2).ne.0) then
                    call stop_all(t_r,'differing orbital should be an beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = FCIBitList(K)
                call SQOperator(tempK,gam,tParity,.true.)

                i = i + 1
                if(i.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Create_beta_T(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_beta_T(2,i) = -1
                else
                    Coup_Create_beta_T(2,i) = 1
                endif
                Coup_Create_beta_inds_T(i) = K
            enddo
            Coup_Create_beta_cum_T(J+1) = i+1
        enddo

        i = 0
        do J = 1,nNp1FCIDet
            !Annihilation of alpha
            do K = 1,nFCIDet
                DiffOrb = ieor(Np1BitList(J),FCIBitList(K))
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
                tempK = FCIBitList(K)
                call SQOperator(tempK,gam,tParity,.false.)

                i = i + 1
                if(i.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Ann_alpha_T(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_alpha_T(2,i) = -1
                else
                    Coup_Ann_alpha_T(2,i) = 1
                endif
                Coup_Ann_alpha_inds_T(i) = K
            enddo
            Coup_Ann_alpha_cum_T(J+1) = i+1
        enddo

        i = 0
        do J = 1,nNp1bFCIDet
            !Annihilation of beta
            do K = 1,nFCIDet
                DiffOrb = ieor(Np1bBitList(J),FCIBitList(K))
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
                if(mod(gam,2).ne.0) then
                    call stop_all(t_r,'differing orbital should be an beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = FCIBitList(K)
                call SQOperator(tempK,gam,tParity,.false.)

                i = i + 1
                if(i.gt.Nmax_Coup) call stop_all(t_r,'Compressed array size too small')
                Coup_Ann_beta_T(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_beta_T(2,i) = -1
                else
                    Coup_Ann_beta_T(2,i) = 1
                endif
                Coup_Ann_beta_inds_T(i) = K
            enddo
            Coup_Ann_beta_cum_T(J+1) = i+1
        enddo

        write(6,"(A)") "Coupling matrices calculated."
        call flush(6)

        !Now, compute the size of the full linear system matrix, and the overlap
        Nmax_Lin = Nmax_N                   !Block 1
        Nmax_Lin = Nmax_Lin + Nmax_N        !Block 2 (N-electron hamiltonian with modified diagonals)
        Nmax_Lin = Nmax_Lin + (Nmax_Nm1+Nmax_Nm1b)*((nImp*2)**2)    !Block 4
        Nmax_Lin = Nmax_Lin + (Nmax_Np1+Nmax_Np1b)*((nImp*2)**2)    !Block 7
        Nmax_Lin = Nmax_Lin + nFCIDet*2     !Block 3 and 3^T    (Diagonal N-electron hamiltonian)
        Nmax_Lin = Nmax_Lin + 2*2*nMax_Coup*nImp*2      !Block 5 and 5^T
        Nmax_Lin = Nmax_Lin + 2*2*nMax_Coup*nImp*2      !Block 6 and 6^T
        Nmax_Lin = Nmax_Lin + 2*2*nMax_Coup*nImp*2      !Block 9 and 9^T
        Nmax_Lin = Nmax_Lin + 2*2*nMax_Coup*nImp*2      !Block 10 and 10^T

        !Now for the overlap matrix size
        Nmax_S = Nmax_N*2       !Block 1 and 2 (diagonal)
        Nmax_S = Nmax_S + (nNm1bFCIDet + nNm1FCIDet)*((nImp*2)**2)  !Block 4 (each subblock is diagonal)
        Nmax_S = Nmax_S + (nNp1bFCIDet + nNp1FCIDet)*((nImp*2)**2)  !Block 7 (each subblock is diagonal)

        write(6,"(A,I12)") "Maximum number of elements in linear system: ",Nmax_Lin
        write(6,"(A,I12)") "Maximum number of elements in overlap matrix: ",Nmax_S

        !Allocate memory for hmailtonian in this system:
        write(6,'(A,F13.3,A)') "Memory required to store linear system: ",    &
            1.5_dp*real(Nmax_Lin,dp)*ComptoMb,' Mb'
        write(6,'(A,F13.3,A)') "Memory required to store overlap: ",    &
            1.5_dp*real(Nmax_S,dp)*ComptoMb,' Mb'
        allocate(LinearSystem(Nmax_Lin),stat=ierr)
        allocate(Overlap(Nmax_S),stat=ierr)
        allocate(LinearSystem_inds(Nmax_Lin),stat=ierr)
        allocate(Overlap_inds(Nmax_S),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        
        iunit = get_free_unit()
        call append_ext_real('EC-TDA_DDResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# 1.Frequency     2.DD_LinearResponse(Re)    3.DD_LinearResponse(Im)    " &
            & //"4.Orthog    5.Norm    6.NewGS   7.OldGS  8.OrthogRHS    9.NI_LR(Re)   10.NI_LR(Im)   11.Iters"
        
        !Store the fock matrix in complex form, so that we can ZGEMM easily
        write(6,'(A,F9.3,A)') "Memory required to store fock matrix: ",real(nSites**2,dp)*ComptoMb,' Mb'
        write(6,'(A,F9.3,A)') "Memory required to store bath states: ",4.0_dp*real(nSites**2,dp)*ComptoMb,' Mb'
        allocate(FockSchmidtComp(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                FockSchmidtComp(j,i) = dcmplx(FockSchmidt(j,i),0.0_dp)
            enddo
        enddo

        InterMem = nCore**2 + nVirt**2 + 8*(EmbSize**2) + 2*(EmbSize*nVirt) + 2*(nCore*EmbSize)
        write(6,'(A,F9.3,A)') "Memory required to store contracted intermediates: ",real(InterMem,dp)*ComptoMb,' Mb'
        !Space for useful intermediates
        allocate(G_ai_G_aj(nCore,nCore))    !Core,Core, contracted over virtual space
        allocate(G_ai_G_bi(VirtStart:VirtEnd,VirtStart:VirtEnd))    !Virt,Virt, contracted over core space
        allocate(G_xa_G_ya(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))    !Active,Active, contracted over virtuals
        allocate(G_xa_F_ab(ActiveStart:ActiveEnd,VirtStart:VirtEnd))
        allocate(G_xa_G_yb_F_ab(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ia_G_xa(1:CoreEnd,ActiveStart:ActiveEnd))
        allocate(F_xi_G_ia_G_ya(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_xa_F_ya(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_xi_G_yi(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ix_F_ij(ActiveStart:ActiveEnd,1:CoreEnd))
        allocate(G_ix_G_jy_F_ji(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ia_G_ix(VirtStart:VirtEnd,ActiveStart:ActiveEnd))
        allocate(F_ax_G_ia_G_iy(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ix_F_iy(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))

        !Since we are not reoptimizing the ground state, the ground state lives in a space which does *not* change with frequency.
        !Therefore, we can compute the RHS of the equations for all frequencies
        allocate(RHS(nLinearSystem))
        RHS(:) = zzero
        do i = 1,nFCIDet
            RHS(i) = dcmplx(HL_Vec(i),0.0_dp)
        enddo
        !Apply V to |0>
        pertsite = 1
        pertsitealpha = 2*pertsite-1
        pertsitebeta = 2*pertsite
        allocate(temp_vecc(nLinearSystem))
        temp_vecc(:) = zzero
        do i = 1,nFCIDet
            if(btest(FCIBitList(i),pertsitealpha-1)) then
                !The perturbation site is occupied
                temp_vecc(i) = temp_vecc(i) + RHS(i)
            endif
            if(btest(FCIBitList(i),pertsitebeta-1)) then
                temp_vecc(i) = temp_vecc(i) + RHS(i)
            endif
        enddo
        RHS(:) = temp_vecc(:)

        !Now, project out the ground state
        !Remove the ground state from the hamiltonian
        !The projector also only lives in the orthogonal space, so it is easy to construct the RHS (it is also real!)
        temp_vecc(:) = RHS(:)
        do i = 1,nFCIDet
            do j = 1,nFCIDet
                temp_vecc(i) = temp_vecc(i) - dcmplx(HL_Vec(i)*HL_Vec(j),0.0_dp)*RHS(j)
            enddo
        enddo
        RHS(:) = dcmplx(-1.0_dp,0.0_dp) * temp_vecc(:)
        !We now have the RHS

        !Check that RHS is orthogonal to the ground state
        testc = dcmplx(0.0_dp,0.0_dp)
        do j = 1,nFCIDet
            testc = testc + RHS(j)*dcmplx(HL_Vec(j),0.0_dp)  
        enddo

        write(6,"(A,2G20.10)") "Orthogonality condition between RHS and |0>: ",testc
            
        if(tMinRes_NonDir) allocate(RHS2(nLinearSystem))
        if(tPrecond_MinRes.and.tMinRes_NonDir) allocate(Precond_Diag(nLinearSystem))
        allocate(Psi1(nLinearSystem))
        Psi1(:) = zzero

        call halt_timer(LR_EC_TDA_Precom)

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            write(6,*) "Calculating linear response for frequency: ",Omega
            call flush(6)
            call set_timer(LR_EC_TDA_HBuild)
        
            LinearSystem(:) = zzero
            Overlap(:) = zzero
            LinearSystem_inds(:) = 0
            Overlap_inds(:) = 0
            ind = 0 
            write(minres_unit,*) "Iteratively solving for frequency: ",Omega

            !**************************    CONTRACTION COEFFICIENTS    **********************

            !First, find the non-interacting solution expressed in the schmidt basis
            call FindSchmidtPert(.false.,Omega,ni_lr)

            !**************************    INTERMEDIATES    *****************************

            !sum_a G_ai^* G_aj
            call ZGEMM('C','N',nCore,nCore,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(VirtStart:VirtEnd,1:CoreEnd),nVirt, &
                SchmidtPert(VirtStart:VirtEnd,1:CoreEnd),nVirt,dcmplx(0.0_dp,0.0_dp),G_ai_G_aj,nCore)
            !sum_i G_ai^* G_bi = sum_i G_ia^* G_bi
            call ZGEMM('C','N',nVirt,nVirt,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,VirtStart:VirtEnd),nCore, &
                SchmidtPert(1:CoreEnd,VirtStart:VirtEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ai_G_bi,nVirt)
            !sum_a G_xa^* G_ya  (where (x,y) in active space)
            call ZGEMM('C','N',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),    &
                SchmidtPert(VirtStart:VirtEnd,ActiveStart:ActiveEnd),nVirt, &
                SchmidtPert(VirtStart:VirtEnd,ActiveStart:ActiveEnd),nVirt,dcmplx(0.0_dp,0.0_dp),G_xa_G_ya,EmbSize)
            !sum_a G_xa F_ab
            call ZGEMM('N','N',EmbSize,nVirt,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(ActiveStart:ActiveEnd,VirtStart:VirtEnd), &
                EmbSize,FockSchmidtComp(VirtStart:VirtEnd,VirtStart:VirtEnd),nVirt,dcmplx(0.0_dp,0.0_dp),G_xa_F_ab,EmbSize)
            !sum_ab G_xa^* G_yb F_ab
            call ZGEMM('C','T',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(VirtStart:VirtEnd,ActiveStart:ActiveEnd), &
                nVirt,G_xa_F_ab(ActiveStart:ActiveEnd,VirtStart:VirtEnd),EmbSize,dcmplx(0.0_dp,0.0_dp),G_xa_G_yb_F_ab,EmbSize)
            !sum_a G_ia^* G_xa
            call ZGEMM('C','T',nCore,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(VirtStart:VirtEnd,1:CoreEnd),nVirt,   &
                SchmidtPert(ActiveStart:ActiveEnd,VirtStart:VirtEnd),EmbSize,dcmplx(0.0_dp,0.0_dp),G_ia_G_xa,nCore)
            !sum_ia F_xi G_ia^* G_ya
            call ZGEMM('N','N',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),FockSchmidtComp(ActiveStart:ActiveEnd,1:CoreEnd),   &
                EmbSize,G_ia_G_xa,nCore,dcmplx(0.0_dp,0.0_dp),F_xi_G_ia_G_ya,EmbSize)
            !sum_a G_xa F_ya
            call ZGEMM('N','T',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(ActiveStart:ActiveEnd,VirtStart:VirtEnd), &
                EmbSize,FockSchmidtComp(ActiveStart:ActiveEnd,VirtStart:VirtEnd),EmbSize,dcmplx(0.0_dp,0.0_dp),G_xa_F_ya,EmbSize)
            !sum_i G_ix^* G_iy
            call ZGEMM('C','N',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore, &
                SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_xi_G_yi,EmbSize)
            !sum_i G_ix F_ij
            call ZGEMM('T','N',EmbSize,nCore,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore,   &
                FockSchmidtComp(1:CoreEnd,1:CoreEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ix_F_ij,EmbSize)
            !sum_ij G_ix^* G_jy F_ji
            call ZGEMM('C','T',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore, &
                G_ix_F_ij,EmbSize,dcmplx(0.0_dp,0.0_dp),G_ix_G_jy_F_ji,EmbSize)
            !sum_i G_ia^* G_ix
            call ZGEMM('C','N',nVirt,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,VirtStart:VirtEnd),nCore,   &
                SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ia_G_ix,nVirt)
            !sum_ia F_ax G_ia^* G_iy
            call ZGEMM('T','N',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),    &
                FockSchmidtComp(VirtStart:VirtEnd,ActiveStart:ActiveEnd),nVirt,G_ia_G_ix,nVirt, &
                dcmplx(0.0_dp,0.0_dp),F_ax_G_ia_G_iy,EmbSize)
            !sum_i G^ix F_iy
            call ZGEMM('T','N',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore, &
                FockSchmidtComp(1:CoreEnd,ActiveStart:ActiveEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ix_F_iy,EmbSize)

            write(6,*) "Intermediates formed..."
            call flush(6)

            !****************************    NORMALIZATION CONSTANTS    ************
            !Calc normalization for the CV block
            CVNorm = 0.0_dp
            do i=1,CoreEnd
                CVNorm = CVNorm + real(G_ai_G_aj(i,i))
            enddo
            CVNorm = CVNorm * 2.0_dp
            ! Active-Virtual block normalizations
            AVNorm(:) = 0.0_dp
            do a = VirtStart,VirtEnd
                do gam = 1,EmbSize
                    gam_spat = gam+CoreEnd
                    !AVNorm(gam) = AVNorm(gam) + real(SchmidtPert(gam_spat,a))**2 + aimag(SchmidtPert(gam_spat,a))**2
                    AVNorm(gam) = AVNorm(gam) + real(conjg(SchmidtPert(gam_spat,a))*SchmidtPert(gam_spat,a),dp)
                enddo
            enddo
            ! Core-Active block normalizations
            CANorm(:) = 0.0_dp
            do gam = 1,nImp*2
                gam_spat = gam+(nOcc-nImp)
                do i = 1,CoreEnd 
                    CANorm(gam) = CANorm(gam) + real(SchmidtPert(i,gam_spat))**2 + aimag(SchmidtPert(i,gam_spat))**2
                enddo
            enddo

            !****************************    DIAGONALS    *************************
            !BLOCK 1
            do i = 1,nFCIDet
                LinearSystem(i) = dcmplx(NFCIHam_cmps(i),0.0_dp) - (dcmplx(HL_Energy + Omega,dDelta))
                Overlap(i) = zone
            enddo
            !BLOCK 2 
            !Find the augmenting factor for the diagonals of this block 
            tempel = zzero 
            do i=1,CoreEnd
                do j=1,CoreEnd
                    tempel = tempel - FockSchmidt(j,i)*G_ai_G_aj(j,i)
                enddo
            enddo
            do b = VirtStart,VirtEnd
                do a = VirtStart,VirtEnd
                    tempel = tempel + FockSchmidt(a,b)*G_ai_G_bi(a,b)
                enddo
            enddo
            tempel = tempel * (2.0_dp/CVNorm)
            j = 0
            do i = CVIndex,AVIndex-1
                j = j + 1
                LinearSystem(i) = (NFCIHam_cmps(j) + tempel) - (dcmplx(HL_Energy + Omega,dDelta))
                Overlap(i) = zone
            enddo
            if(j.ne.nFCIDet) call stop_all(t_r,'Indexing problem')
            !BLOCK 4
            do gam = 1,EmbSize
                !Run through active space orbitals
                gam_spat = gam + (nOcc-nImp)
                !gam_ind gives the starting index of the AV diagonal block for orbital gam
                gam_ind = AVIndex + (gam-1)*nFullNm1
                
                do i = gam_ind,gam_ind+nNm1bFCIDet-1
                    !Run through beta block
                    LinearSystem(i) = ((dcmplx(Nm1FCIHam_beta_cmps(i-gam_ind+1),0.0_dp)*G_xa_G_ya(gam_spat,gam_spat) +  &
                        G_xa_G_yb_F_ab(gam_spat,gam_spat)) / (sqrt(AVNorm(gam)*AVNorm(gam)))) - (dcmplx(HL_Energy + Omega,dDelta))
                    Overlap(i) = zone
                enddo
                do i = gam_ind+nNm1bFCIDet,gam_ind+nFullNm1-1
                    !And alpha block
                    LinearSystem(i) = ((dcmplx(Nm1FCIHam_alpha_cmps(i-(gam_ind+nNm1bFCIDet)+1),0.0_dp)* &
                        G_xa_G_ya(gam_spat,gam_spat) +   &
                        G_xa_G_yb_F_ab(gam_spat,gam_spat)) / (sqrt(AVNorm(gam)*AVNorm(gam)))) -     &
                        (dcmplx(HL_Energy + Omega,dDelta))
                    Overlap(i) = zone
                enddo
            enddo
            !BLOCK 7 (similar to block 4)
            do gam = 1,EmbSize
                gam_spat = gam+(nOcc-nImp)
                gam_ind = CAIndex + (gam-1)*nFullNp1
                do i = gam_ind,gam_ind+nNp1FCIDet-1
                    LinearSystem(i) = ((dcmplx(Np1FCIHam_alpha_cmps(i-gam_ind+1),0.0_dp)*G_xi_G_yi(gam_spat,gam_spat) -     &
                        G_ix_G_jy_F_ji(gam_spat,gam_spat)) / (sqrt(CANorm(gam)*CANorm(gam)))) - (dcmplx(HL_Energy + Omega,dDelta))
                    Overlap(i) = zone
                enddo
                do i = gam_ind+nNp1FCIDet,gam_ind+nFullNp1-1
                    LinearSystem(i) = ((dcmplx(Np1FCIHam_beta_cmps(i-(gam_ind+nNp1FCIDet)+1),0.0_dp)*   &
                        G_xi_G_yi(gam_spat,gam_spat) -     &
                        G_ix_G_jy_F_ji(gam_spat,gam_spat)) / (sqrt(CANorm(gam)*CANorm(gam)))) -     &
                        (dcmplx(HL_Energy + Omega,dDelta))
                    Overlap(i) = zone
                enddo
            enddo

            !*****************************    OFF-DIAGONAL MATRIX ELEMENT    ***********

            !Precompute the value of all elements in block 3
            tempel = dcmplx(0.0_dp,0.0_dp)
            do i = 1,CoreEnd
                do a = VirtStart,VirtEnd
                    tempel = tempel + SchmidtPert(a,i)*FockSchmidt(a,i)
                enddo
            enddo
            tempel = tempel * (2.0_dp/(sqrt(CVNorm)))

            !Set up indexing for the off-diagonal elements of the linear system
            ind = nLinearSystem + 1
            LinearSystem_inds(1) = nLinearSystem + 2
            ind_S = nLinearSystem + 1
            Overlap_inds(1) = nLinearSystem + 2
            
            !First, the top nFCIDet rows, running through blocks 1, 3, 6, 10
            do i = 1,nFCIDet
                !Block 1
                do k = NFCIHam_inds(i),NFCIHam_inds(i+1)-1
                    ind = ind + 1
                    if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                    LinearSystem(ind) = NFCIHam_cmps(k)
                    LinearSystem_inds(ind) = NFCIHam_inds(k)
                enddo
                !Block 3 is diagonal, and so there is only *1* non-zero element per row
                ind = ind + 1
                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                LinearSystem(ind) = tempel
                LinearSystem_inds(ind) = nFCIDet + i
                !Block 6
                !Run through orbitals in active space
                do beta = 1,EmbSize
                    !Run through non-zero connections in that block
                    do k = Coup_Create_alpha_cum(i),Coup_Create_alpha_cum(i+1)-1

                        matel = G_xa_F_ya(beta+CoreEnd,Coup_Create_alpha(1,k))*Coup_Create_alpha(2,k) / sqrt(AVNorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = AVIndex + (beta-1)*nFullNm1 + Coup_Create_alpha_inds(k) - 1
                        endif
                    enddo
                    !Now for other spin type
                    do k = Coup_Create_beta_cum(i),Coup_Create_beta_cum(i+1)-1
                        matel = G_xa_F_ya(beta+CoreEnd,Coup_Create_beta(1,k))*Coup_Create_beta(2,k) / sqrt(AVNorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = AVIndex + (beta-1)*nFullNm1 + nNm1bFCIDet + Coup_Create_beta_inds(k) - 1
                        endif
                    enddo
                enddo
                !Block 10 (Similar to block 6)
                do beta = 1,EmbSize
                    do k = Coup_Ann_alpha_cum(i),Coup_Ann_alpha_cum(i+1)-1
                        matel = -G_ix_F_iy(beta+CoreEnd,Coup_Ann_alpha(1,k))*Coup_Ann_alpha(2,k) / sqrt(CANorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = CAIndex + (beta-1)*nFullNp1 + Coup_Ann_alpha_inds(k) - 1
                        endif
                    enddo
                    !Other spin
                    do k = Coup_Ann_beta_cum(i),Coup_Ann_beta_cum(i+1)-1
                        matel = -G_ix_F_iy(beta+CoreEnd,Coup_Ann_beta(1,k))*Coup_Ann_beta(2,k) / sqrt(CANorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = CAIndex + (beta-1)*nFullNp1 + nNp1FCIDet + Coup_Ann_beta_inds(k) - 1
                        endif
                    enddo
                enddo
                LinearSystem_inds(i+1) = ind + 1
                Overlap_inds(i+1) = ind_S + 1
            enddo   !Finish row i (Block 1 column)
            
            !Now, the next nFCIDet rows, running through blocks 3^T, 2, 5, 9
            do i = 1,nFCIDet
                !Block 3^T
                ind = ind + 1
                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                LinearSystem(ind) = conjg(tempel)
                LinearSystem_inds(ind) = i
                !Block 2
                do k = NFCIHam_inds(i),NFCIHam_inds(i+1)-1
                    ind = ind + 1
                    if(ind.gt.NMax_Lin) call stop_all(t_r,'Compressed array too small')
                    LinearSystem(ind) = NFCIHam_cmps(k)
                    LinearSystem_inds(ind) = NFCIHam_inds(k) + nFCIDet
                enddo
                !Block 5
                do beta = 1,EmbSize
                    do k = Coup_Create_alpha_cum(i),Coup_Create_alpha_cum(i+1)-1
                        matel = -F_xi_G_ia_G_ya(Coup_Create_alpha(1,k),beta+CoreEnd)*Coup_Create_alpha(2,k)     &
                            / sqrt(CVNorm*AVNorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = AVIndex + (beta-1)*nFullNm1 + Coup_Create_alpha_inds(k) - 1
                        endif
                    enddo
                    !Now for other spin type
                    do k = Coup_Create_beta_cum(i),Coup_Create_beta_cum(i+1)-1
                        matel = -F_xi_G_ia_G_ya(Coup_Create_beta(1,k),beta+CoreEnd)*Coup_Create_beta(2,k)   &
                            / sqrt(CVNorm*AVNorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = AVIndex + (beta-1)*nFullNm1 + nNm1bFCIDet +    &
                                Coup_Create_beta_inds(k) - 1
                        endif
                    enddo
                enddo
                !Block 9
                do beta = 1,EmbSize
                    do k = Coup_Ann_alpha_cum(i),Coup_Ann_alpha_cum(i+1)-1
                        matel = -F_ax_G_ia_G_iy(Coup_Ann_alpha(1,k),beta+CoreEnd)*Coup_Ann_alpha(2,k) /     &
                            sqrt(CVNorm*CANorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = CAIndex + (beta-1)*nFullNp1 + Coup_Ann_alpha_inds(k) - 1
                        endif
                    enddo
                    !Other spin
                    do k = Coup_Ann_beta_cum(i),Coup_Ann_beta_cum(i+1)-1
                        matel = -F_ax_G_ia_G_iy(Coup_Ann_beta(1,k),beta+CoreEnd)*Coup_Ann_beta(2,k) /   &
                            sqrt(CVNorm*CANorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = matel
                            LinearSystem_inds(ind) = CAIndex + (beta-1)*nFullNp1 + nNp1FCIDet +     &
                                Coup_Ann_beta_inds(k) - 1
                        endif
                    enddo
                enddo
                LinearSystem_inds(nFCIDet+i+1) = ind + 1
                Overlap_inds(nFCIDet+i+1) = ind_S + 1
            enddo   !Finish row i (Block 2 columns)

            !Now, the next nFullNm1*nImp*2 rows, running through blocks 6^T, 5^T, 4, 8
            row = AVIndex
            do beta = 1,EmbSize
                do i = 1,nNm1bFCIDet
                    !Cols through 6^T
                    do k = Coup_Create_alpha_cum_T(i),Coup_Create_alpha_cum_T(i+1)-1
                        matel = G_xa_F_ya(beta+CoreEnd,Coup_Create_alpha_T(1,k))*   &
                            Coup_Create_alpha_T(2,k) / sqrt(AVNorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Create_alpha_inds_T(k)
                        endif
                    enddo
                    !Cols through 5^T
                    do k = Coup_Create_alpha_cum_T(i),Coup_Create_alpha_cum_T(i+1)-1
                        matel = -F_xi_G_ia_G_ya(Coup_Create_alpha_T(1,k),beta+CoreEnd)* &
                            Coup_Create_alpha_T(2,k) / sqrt(CVNorm*AVNorm(beta)) 
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Create_alpha_inds_T(k) + NFCIDet
                        endif
                    enddo
                    !Cols through 4. This will not necessarily be in the correct order any more
                    gam1_spat = beta+(nOcc-nImp)
                    gam1_ind = AVIndex + (beta-1)*nFullNm1
                    do gam2 = 1,EmbSize
                        gam2_spat = gam2+(nOcc-nImp)
                        gam2_ind = AVIndex + (gam2-1)*nFullNm1
                        do k = Nm1FCIHam_beta_inds(i),Nm1FCIHam_beta_inds(i+1)-1
                            !Off-diagonals
                            matel = dcmplx(Nm1FCIHam_beta_cmps(k),0.0_dp)*G_xa_G_ya(gam1_spat,gam2_spat) /  &
                                sqrt(AVNorm(beta)*AVNorm(gam2)) 
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = Nm1FCIHam_beta_inds(k) + gam2_ind - 1
                            endif
                        enddo
                        if(gam2.ne.beta) then
                            !diagonal element for this gamma
                            !Don't include when gam2 = beta, since this is now a true diagonal of the matrix and has already been considered
                            !We also need to include the off-diagonal element for the overlap matrix here
                            matel_S = G_xa_G_ya(gam1_spat,gam2_spat) / sqrt(AVNorm(beta)*AVNorm(gam2))
                            if(abs(matel_S).gt.CompressThresh) then
                                ind_S = ind_S + 1
                                if(ind_S.gt.Nmax_S) call stop_all(t_r,'Compressed S array too small')
                                Overlap(ind_S) = matel_S
                                Overlap_inds(ind_S) = gam2_ind + i - 1
                            endif
                            !Now for the Hessian
                            matel = ((dcmplx(Nm1FCIHam_beta_cmps(i),0.0_dp)*G_xa_G_ya(gam1_spat,gam2_spat) +    &
                                G_xa_G_yb_F_ab(gam1_spat,gam2_spat)) / sqrt(AVNorm(beta)*AVNorm(gam2))) -   &
                                matel_S*dcmplx(HL_Energy + Omega,dDelta)
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = gam2_ind + i - 1
                            endif
                        endif
                    enddo
                    !Block 8 is mercifully 0, so we can ignore
                    row = row + 1
                    LinearSystem_inds(row) = ind + 1
                    Overlap_inds(row) = ind_S + 1
                enddo
                !Now for other spin
                do i = 1,nNm1FCIDet
                    !Cols through 6^T
                    do k = Coup_Create_beta_cum_T(i),Coup_Create_beta_cum_T(i+1)-1
                        matel = G_xa_F_ya(beta+CoreEnd,Coup_Create_beta_T(1,k))*Coup_Create_beta_T(2,k)     &
                            / sqrt(AVNorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Create_beta_inds_T(k)
                        endif
                    enddo
                    !Cols through 5^T
                    do k = Coup_Create_beta_cum_T(i),Coup_Create_beta_cum_T(i+1)-1
                        matel = -F_xi_G_ia_G_ya(Coup_Create_beta_T(1,k),beta+CoreEnd)*Coup_Create_beta_T(2,k) / &
                            sqrt(CVNorm*AVNorm(beta)) 
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Create_beta_inds_T(k) + NFCIDet
                        endif
                    enddo
                    !Cols through 4. This will not necessarily be in the correct order any more
                    gam1_spat = beta+(nOcc-nImp)
                    gam1_ind = AVIndex + (beta-1)*nFullNm1
                    do gam2 = 1,EmbSize
                        gam2_spat = gam2+(nOcc-nImp)
                        gam2_ind = AVIndex + (gam2-1)*nFullNm1
                        do k = Nm1FCIHam_alpha_inds(i),Nm1FCIHam_alpha_inds(i+1)-1
                            !Off-diagonals
                            matel = dcmplx(Nm1FCIHam_alpha_cmps(k),0.0_dp)*G_xa_G_ya(gam1_spat,gam2_spat) / &
                                sqrt(AVNorm(beta)*AVNorm(gam2)) 
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = Nm1FCIHam_alpha_inds(k) + gam2_ind + nNm1bFCIDet - 1
                            endif
                        enddo
                        if(gam2.ne.beta) then
                            !diagonal element for this gamma
                            !Don't include when gam2 = beta, since this is now a true diagonal of the matrix and has already been considered
                            !We also need to include the off-diagonal element for the overlap matrix here
                            matel_S = G_xa_G_ya(gam1_spat,gam2_spat) / sqrt(AVNorm(beta)*AVNorm(gam2))
                            if(abs(matel_S).gt.CompressThresh) then
                                ind_S = ind_S + 1
                                if(ind_S.gt.Nmax_S) call stop_all(t_r,'Compressed S array too small')
                                Overlap(ind_S) = matel_S
                                Overlap_inds(ind_S) = gam2_ind + nNm1bFCIDet + i - 1
                            endif
                            !Now for the Hessian
                            matel = ((dcmplx(Nm1FCIHam_alpha_cmps(i),0.0_dp)*G_xa_G_ya(gam1_spat,gam2_spat) +   &
                                G_xa_G_yb_F_ab(gam1_spat,gam2_spat)) / sqrt(AVNorm(beta)*AVNorm(gam2))) -   &
                                matel_S*dcmplx(HL_Energy + Omega,dDelta)
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = gam2_ind + nNm1bFCIDet + i - 1
                            endif
                        endif
                    enddo
                    !Block 8 is mercifully 0, so we can ignore
                    row = row + 1
                    LinearSystem_inds(row) = ind + 1
                    Overlap_inds(row) = ind_S + 1
                enddo
            enddo   !Finish row beta (Block 4 columns)

            !Now for the final nFullNp1*nImp*2 rows, running through blocks 10^T, 9^T, 8^T, 7
            row = CAIndex
            do beta = 1,EmbSize
                do i = 1,nNp1FCIDet
                    !Cols through 10^T
                    do k = Coup_Ann_alpha_cum_T(i),Coup_Ann_alpha_cum_T(i+1)-1
                        matel = -G_ix_F_iy(beta+CoreEnd,Coup_Ann_alpha_T(1,k))*Coup_Ann_alpha_T(2,k) /  &
                            sqrt(CANorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Ann_alpha_inds_T(k)
                        endif
                    enddo
                    !Cols through 9^T
                    do k = Coup_Ann_alpha_cum_T(i),Coup_Ann_alpha_cum_T(i+1)-1
                        matel = -F_ax_G_ia_G_iy(Coup_Ann_alpha_T(1,k),beta+CoreEnd)*Coup_Ann_alpha_T(2,k) / &
                            sqrt(CVNorm*CANorm(beta)) 
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Ann_alpha_inds_T(k) + NFCIDet
                        endif
                    enddo
                    !Block 8^T is mercifully 0, so we can ignore
                    !Cols through 7. This will not necessarily be in the correct order any more
                    gam1_spat = beta+(nOcc-nImp)
                    gam1_ind = CAIndex + (beta-1)*nFullNp1
                    do gam2 = 1,EmbSize
                        gam2_spat = gam2+(nOcc-nImp)
                        gam2_ind = CAIndex + (gam2-1)*nFullNp1
                        do k = Np1FCIHam_alpha_inds(i),Np1FCIHam_alpha_inds(i+1)-1
                            !Off-diagonals
                            matel = dcmplx(Np1FCIHam_alpha_cmps(k),0.0_dp)*G_xi_G_yi(gam1_spat,gam2_spat) / &
                                sqrt(CANorm(beta)*CANorm(gam2)) 
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = Np1FCIHam_alpha_inds(k) + gam2_ind - 1
                            endif
                        enddo
                        if(gam2.ne.beta) then
                            !diagonal element for this gamma
                            !Don't include when gam2 = beta, since this is now a true diagonal of the matrix and has already been considered
                            !We also need to include the off-diagonal element for the overlap matrix here
                            matel_S = G_xi_G_yi(gam1_spat,gam2_spat) / sqrt(CANorm(beta)*CANorm(gam2))
                            if(abs(matel_S).gt.CompressThresh) then
                                ind_S = ind_S + 1
                                if(ind_S.gt.Nmax_S) call stop_all(t_r,'Compressed S array too small')
                                Overlap(ind_S) = matel_S
                                Overlap_inds(ind_S) = gam2_ind + i - 1
                            endif
                            !Now for Hessian contribution
                            matel = ((dcmplx(Np1FCIHam_alpha_cmps(i),0.0_dp)*G_xi_G_yi(gam1_spat,gam2_spat) -   &
                                G_ix_G_jy_F_ji(gam1_spat,gam2_spat)) / sqrt(CANorm(beta)*CANorm(gam2))) -   &
                                matel_S*dcmplx(HL_Energy + Omega,dDelta)
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = gam2_ind + i - 1
                            endif
                        endif
                    enddo
                    row = row + 1
                    LinearSystem_inds(row) = ind + 1
                    Overlap_inds(row) = ind_S + 1
                enddo
                !Now for other spin
                do i = 1,nNp1bFCIDet
                    !Cols through 10^T
                    do k = Coup_Ann_beta_cum_T(i),Coup_Ann_beta_cum_T(i+1)-1
                        matel = -G_ix_F_iy(beta+CoreEnd,Coup_Ann_beta_T(1,k))*Coup_Ann_beta_T(2,k) /    &
                            sqrt(CANorm(beta))
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Ann_beta_inds_T(k)
                        endif
                    enddo
                    !Cols through 9^T
                    do k = Coup_Ann_beta_cum_T(i),Coup_Ann_beta_cum_T(i+1)-1
                        matel = -F_ax_G_ia_G_iy(Coup_Ann_beta_T(1,k),beta+CoreEnd)*Coup_Ann_beta_T(2,k) /   &
                            sqrt(CVNorm*CANorm(beta)) 
                        if(abs(matel).gt.CompressThresh) then
                            ind = ind + 1
                            if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                            LinearSystem(ind) = conjg(matel)
                            LinearSystem_inds(ind) = Coup_Ann_beta_inds_T(k) + NFCIDet
                        endif
                    enddo
                    !Block 8^T is mercifully 0, so we can ignore
                    !Cols through 7. This will not necessarily be in the correct order any more
                    gam1_spat = beta+(nOcc-nImp)
                    gam1_ind = CAIndex + (beta-1)*nFullNp1
                    do gam2 = 1,EmbSize
                        gam2_spat = gam2+(nOcc-nImp)
                        gam2_ind = CAIndex + (gam2-1)*nFullNp1
                        do k = Np1FCIHam_beta_inds(i),Np1FCIHam_beta_inds(i+1)-1
                            !Off-diagonals
                            matel = dcmplx(Np1FCIHam_beta_cmps(k),0.0_dp)*G_xi_G_yi(gam1_spat,gam2_spat) /  &
                                sqrt(CANorm(beta)*CANorm(gam2)) 
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = Np1FCIHam_beta_inds(k) + nNp1FCIDet + gam2_ind - 1
                            endif
                        enddo
                        if(gam2.ne.beta) then
                            !diagonal element for this gamma
                            !Don't include when gam2 = beta, since this is now a true diagonal of the matrix and has already been considered
                            !We also need to include the off-diagonal element for the overlap matrix here
                            matel_S = G_xi_G_yi(gam1_spat,gam2_spat) / sqrt(CANorm(beta)*CANorm(gam2))
                            if(abs(matel_S).gt.CompressThresh) then
                                ind_S = ind_S + 1
                                if(ind_S.gt.Nmax_S) call stop_all(t_r,'Compressed S array too small')
                                Overlap(ind_S) = matel_S
                                Overlap_inds(ind_S) = gam2_ind + nNp1FCIDet + i - 1
                            endif
                            !Now for Hessian
                            matel = ((dcmplx(Np1FCIHam_beta_cmps(i),0.0_dp)*G_xi_G_yi(gam1_spat,gam2_spat) -    &
                                G_ix_G_jy_F_ji(gam1_spat,gam2_spat)) / sqrt(CANorm(beta)*CANorm(gam2))) -   &
                                matel_S*dcmplx(HL_Energy + Omega,dDelta)
                            if(abs(matel).gt.CompressThresh) then
                                ind = ind + 1
                                if(ind.gt.Nmax_Lin) call stop_all(t_r,'Compressed array too small')
                                LinearSystem(ind) = matel
                                LinearSystem_inds(ind) = gam2_ind + nNp1FCIDet + i - 1
                            endif
                        endif
                    enddo
                    row = row + 1
                    LinearSystem_inds(row) = ind + 1
                    Overlap_inds(row) = ind_S + 1
                enddo
            enddo   !Finish row beta (Block 4 columns)

            write(6,*) "Hessian constructed successfully..."
            call flush(6)
            call halt_timer(LR_EC_TDA_HBuild)

            !*********************   Hessian construction finished   **********************
            
            call set_timer(LR_EC_TDA_SolveLR)

            !Now solve these linear equations
            zDirMV_Mat_cmprs => LinearSystem
            zDirMV_Mat_cmprs_inds => LinearSystem_inds
            if(tMinRes_NonDir) then
                maxminres_iter_ip = int(maxminres_iter,ip)
                minres_unit_ip = int(minres_unit,ip)
                nLinearSystem_ip = int(nLinearSystem,ip)
                call setup_RHS(nLinearSystem,RHS,RHS2)
                if(tPrecond_MinRes) then
                    call FormPrecond(nLinearSystem)
                    call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS2,nout=minres_unit_ip,x=Psi1, &
                        itnlim=maxminres_iter_ip,Msolve=zPreCond,istop=info_ip,rtol=rtol_LR,itn=iters, &
                        startguess=tReuse_LS)
                else
                    call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS2,nout=minres_unit_ip,x=Psi1, &
                        itnlim=maxminres_iter_ip,istop=info_ip,rtol=rtol_LR,itn=iters,startguess=tReuse_LS)
                endif
                info = info_ip
                if(info.gt.7) write(6,*) "info: ",info
                if(info.eq.8) write(6,"(A,I9)") "Linear equation solver hit iteration limit: ",iters
                if((info.eq.9).or.(info.eq.10).or.(info.eq.11)) then
                    call stop_all(t_r,'Input matrices to linear solver incorrect')
                endif
                if(info.gt.11) call stop_all(t_r,'Linear equation solver failed')
                info = 0
            elseif(tGMRES_NonDir) then
                call GMRES_Solve(nLinearSystem,RHS,minres_unit,maxminres_iter,rtol_LR,tPrecond_MinRes,tReuse_LS,iters,Psi1,info)
            endif
            zDirMV_Mat_cmprs => null()
            zDirMV_Mat_cmprs_inds => null()
            if(info.ne.0) then 
                write(6,*) "INFO: ",info
                !call stop_all(t_r,'Solving Linear system failed') 
                call warning(t_r,'Solving linear system failed - skipping this frequency')
            else
                !Find normalization of first-order wavefunction
                !Test whether the perturbed wavefunction is orthogonal to the ground state wavefunction
                dNorm = dcmplx(0.0_dp,0.0_dp)
                dOrthog = dcmplx(0.0_dp,0.0_dp)

                !Apply the compressed overlap matrix to the ground state: S|1>
                do i = 1,nLinearSystem
                    temp_vecc(i) = Overlap(i)*Psi1(i)   !Diagonal terms
                    do j = Overlap_inds(i),Overlap_inds(i+1)-1
                        temp_vecc(i) = temp_vecc(i) + Overlap(j)*Psi1(Overlap_inds(j))
                    enddo
                enddo
                !Norm now just dot product
                do i = 1,nLinearSystem
                    dNorm = dNorm + conjg(Psi1(i))*temp_vecc(i)
                enddo

                !Since |0> just in uncontracted space, orthogonality easy to calculation
                do i = 1,nFCIDet
                    dOrthog = dOrthog + dcmplx(HL_Vec(i),0.0_dp)*Psi1(i)
                enddo
!                write(6,*) "Normalization of first-order wavefunction: ",dNorm
!                write(6,*) "Orthogonality of first-order wavefunction: ",dOrthog

                if(abs(dOrthog).gt.1.0e-4_dp) then
                    !call warning(t_r,'First order wavefunction not orthogonal - skipping this frequency')
                    write(6,"(A,G25.10)") 'First order solution not orthogonal',abs(dOrthog)
                elseif(abs(dOrthog).gt.1.0e-8_dp) then
                    write(6,"(A,G25.10)") 'First order solution not strictly orthogonal',abs(dOrthog)
                endif

                !Calculating the response is also trivial - <RHS|S|1>
                !RHS again only lives in space of |0>
                ResponseFn = dcmplx(0.0_dp,0.0_dp)
                do i = 1,nFCIDet
                    ResponseFn = ResponseFn - RHS(i)*temp_vecc(i)
                enddo
                
                write(iunit,"(10G22.10,I8)") Omega,real(ResponseFn),-aimag(ResponseFn), &
                    abs(dOrthog),abs(dNorm),HL_Energy,HL_Energy,abs(testc),real(ni_lr),-aimag(ni_lr),iters
                call flush(iunit)

            endif

            Omega = Omega + Omega_Step

            call halt_timer(LR_EC_TDA_SolveLR)
        enddo   !End loop over omega

        deallocate(RHS,Psi1,temp_vecc)
        if(tMinRes_NonDir) deallocate(RHS2)
        if(tPrecond_MinRes.and.tMinRes_NonDir) deallocate(Precond_Diag)
        deallocate(LinearSystem,Overlap,LinearSystem_inds,Overlap_inds)
        close(iunit)
        close(minres_unit)

        !Deallocate determinant lists
        if(allocated(FCIDetList)) deallocate(FCIDetList)
        if(allocated(FCIBitList)) deallocate(FCIBitList)
        if(allocated(UMat)) deallocate(UMat)
        if(allocated(TMat)) deallocate(TMat)
        if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
        if(allocated(Nm1BitList)) deallocate(Nm1BitList)
        if(allocated(Np1FCIDetList)) deallocate(Np1FCIDetList)
        if(allocated(Np1BitList)) deallocate(Np1BitList)
        if(allocated(Nm1bFCIDetList)) deallocate(Nm1bFCIDetList)
        if(allocated(Nm1bBitList)) deallocate(Nm1bBitList)
        if(allocated(Np1bFCIDetList)) deallocate(Np1bFCIDetList)
        if(allocated(Np1bBitList)) deallocate(Np1bBitList)

        !Stored intermediates
        deallocate(NFCIHam_cmps,Nm1FCIHam_alpha_cmps,Nm1FCIHam_beta_cmps,Np1FCIHam_alpha_cmps,Np1FCIHam_beta_cmps)
        deallocate(NFCIHam_inds,Nm1FCIHam_alpha_inds,Nm1FCIHam_beta_inds,Np1FCIHam_alpha_inds,Np1FCIHam_beta_inds)
        deallocate(AVNorm,CANorm)
        deallocate(Coup_Create_alpha,Coup_Create_alpha_T,Coup_Create_beta)
        deallocate(Coup_Create_beta_T,Coup_Ann_alpha,Coup_Ann_alpha_T)
        deallocate(Coup_Ann_beta,Coup_Ann_beta_T,Coup_Ann_beta_inds,Coup_Ann_beta_inds_T)
        deallocate(Coup_Create_alpha_inds,Coup_Create_alpha_inds_T,Coup_Create_beta_inds)
        deallocate(Coup_Create_beta_inds_T,Coup_Ann_alpha_inds,Coup_Ann_alpha_inds_T)
        deallocate(Coup_Create_alpha_cum,Coup_Create_alpha_cum_T,Coup_Create_beta_cum)
        deallocate(Coup_Create_beta_cum_T,Coup_Ann_alpha_cum,Coup_Ann_alpha_cum_T)
        deallocate(Coup_Ann_beta_cum,Coup_Ann_beta_cum_T)
        deallocate(FockSchmidtComp,G_ai_G_aj,G_ai_G_bi,G_xa_G_ya,G_xa_F_ab,G_xa_G_yb_F_ab)
        deallocate(G_ia_G_xa,F_xi_G_ia_G_ya,G_xa_F_ya,G_xi_G_yi,G_ix_F_ij,G_ix_G_jy_F_ji)
        deallocate(G_ia_G_ix,F_ax_G_ia_G_iy,G_ix_F_iy)


    end subroutine NonIntExCont_TDA_MCLR_DD_Cmprs
    
    !Calculate linear response for charged excitations - add the hole creation to particle creation
    !This routine uses a compressed matrix formalism
    subroutine NonIntExCont_TDA_MCLR_Charged_Cmprs()
        use mat_tools, only: add_localpot_comp_inplace
        use utils, only: get_free_unit,append_ext_real,append_ext
        use Davidson, only: Comp_NonDir_Davidson
        use fitting, only: Fit_SE
        use zminresqlpModule, only: MinresQLP  
        use DetToolsData
        use solvers, only: CountSizeCompMat,CreateIntMats
        use DetTools, only: tospat,umatind,gendets
        use DetBitOps, only: CountBits
        implicit none
        complex(dp), allocatable :: NFCIHam_cmps(:),Np1FCIHam_cmps(:),Nm1FCIHam_cmps(:)
        integer, allocatable :: NFCIHam_inds(:),Np1FCIHam_inds(:),Nm1FCIHam_inds(:)
        real(dp), allocatable :: dNorm_p(:),dNorm_h(:)
        complex(dp), allocatable , target :: LinearSystem_p(:),LinearSystem_h(:)
        integer, allocatable , target :: LinearSystem_p_inds(:),LinearSystem_h_inds(:)
        complex(dp), allocatable :: Cre_0(:,:),Ann_0(:,:),ResponseFn_p(:,:),ResponseFn_h(:,:)
        complex(dp), allocatable :: Gc_a_F_ax_Bra(:,:),Gc_a_F_ax_Ket(:,:),Gc_b_F_ab(:,:)
        complex(dp), allocatable :: ResponseFn_Mat(:,:),Ga_i_F_xi_Ket(:,:)
        complex(dp), allocatable :: Psi1_p(:),Psi1_h(:),Ga_i_F_xi_Bra(:,:),Ga_i_F_ij(:,:),ni_lr_Mat(:,:)
        complex(dp), allocatable :: Psi_0(:),RHS(:)
        complex(dp), allocatable :: NI_LRMat_Cre(:,:),NI_LRMat_Ann(:,:)
        integer, allocatable :: Coup_Ann(:,:),Coup_Create(:,:),Coup_Create_inds(:)
        integer, allocatable :: Coup_Ann_inds(:),Coup_Create_cum(:),Coup_Ann_cum(:)
        integer, allocatable :: Coup_Create_T(:,:),Coup_Ann_T(:,:),Coup_Create_inds_T(:)
        integer, allocatable :: Coup_Ann_inds_T(:),Coup_Create_cum_T(:),Coup_Ann_cum_T(:)
        integer :: i,a,j,k,ActiveEnd,ActiveStart,CoreEnd,ierr,info,iunit,nCore
        integer :: nLinearSystem,nOrbs,nVirt,VIndex,VirtStart,VirtEnd,ind_p,ind_h
        integer :: nLinearSystem_h,nGSSpace,Np1GSInd,Nm1GSInd,minres_unit
        integer :: maxminres_iter,nImp_GF,pertsite,isize,iunit_tmp,DiffOrb
        integer :: Nmax_N,Nmax_Nm1,Nmax_Np1,Nmax_Coup_Create,Nmax_Coup_Ann,Nmax_Lin_p,Nmax_Lin_h
        integer(ip) :: nLinearSystem_ip,minres_unit_ip,info_ip,maxminres_iter_ip,iters_p,iters_h
        real(dp) :: Omega,GFChemPot,mu,SpectralWeight,Prev_Spec,AvdNorm_p,AvdNorm_h
        real(dp) :: Diff_GF
        complex(dp) :: ResponseFn,tempel,ni_lr,ni_lr_p,ni_lr_h,AvResFn_p,AvResFn_h
        complex(dp) :: zdotc,VNorm,CNorm,matel
        logical :: tFirst,texist
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='NonIntExCont_TDA_MCLR_Charged_Cmprs'

        maxminres_iter = iMinRes_MaxIter
        iters_p = 0
        iters_h = 0

        SpectralWeight = 0.0_dp
        Prev_Spec = 0.0_dp
        tFirst = .true.
        zDirMV_Mat => null()

        call set_timer(LR_EC_GF_Precom)

        if(tBetaExcit) then
            write(6,"(A)") "Calculating non-interacting EC MR-TDA LR system for *hole* & *particle* alpha spin-orbital pert..."
        else
            write(6,"(A)") "Calculating non-interacting EC MR-TDA LR system for *hole* & *particle* beta spin-orbital pert..."
        endif
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(tMinRes_NonDir) then
            if(tPreCond_MinRes) then
                write(6,"(A)") "Solving linear system with iterative non-direct preconditioned MinRes-QLP algorithm"
            else
                write(6,"(A)") "Solving linear system with iterative non-direct MinRes-QLP algorithm"
            endif
            write(6,"(A,G22.10)") "Tolerance for solution of linear system: ",rtol_LR
            write(6,"(A,G22.10)") "Maximum iterations for each solution: ",maxminres_iter
        elseif(tGMRES_NonDir) then
            if(tPreCond_MinRes) then
                write(6,'(A)') "Solving linear system with preconditioned iterative non-direct GMRES algorithm"
            else
                write(6,'(A)') "Solving linear system with iterative non-direct GMRES algorithm"
            endif
            write(6,"(A,G22.10)") "Tolerance for solution of linear system: ",rtol_LR
            write(6,"(A,G22.10)") "Maximum iterations for each solution: ",maxminres_iter
            write(6,"(A,I8)") "Number of krylov subspace vectors to store: ",nKrylov
        else
            call stop_all(t_r,"Require iterative linear equation solver with compressed matrices")
        endif
        if(tLR_ReoptGS) call stop_all(t_r,'Cannot use compressed matrices with reoptimization of ground state currently')
        if(tSC_LR) call stop_all(t_r,'Cannot use compressed matrices with self-consistent greens functions currently')

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
        write(6,"(A)") "Computing size of compressed Hamiltonian matrices"
        call flush(6)
        if(tReadMats) then
            inquire(file='CompressHam_N',exist=texist)
            if(texist) then
                write(6,*) "N-electron hamiltonian found"
                iunit_tmp = get_free_unit()
                open(iunit_tmp,file='CompressHam_N',status='old')
                read(iunit_tmp,*) isize
                Nmax_N = isize
                rewind(iunit_tmp)
                close(iunit_tmp)
                write(6,*) "Read in NMAX_N = ",Nmax_N
            else
                if(iHamSize_N.eq.0) then
                    call CountSizeCompMat(FCIDetList,Elec,nFCIDet,Nmax_N,FCIBitList)
                    iHamSize_N = Nmax_N
                else
                    write(6,"(A,I16)") "Assuming size of compressed N electron hamiltonian is: ",iHamSize_N
                    Nmax_N = iHamSize_N
                endif
            endif
        else
            if(iHamSize_N.eq.0) then
                call CountSizeCompMat(FCIDetList,Elec,nFCIDet,Nmax_N,FCIBitList)
                iHamSize_N = Nmax_N
            else
                write(6,"(A,I16)") "Assuming size of compressed N electron hamiltonian is: ",iHamSize_N
                Nmax_N = iHamSize_N
            endif
        endif
        if(iHamSize_Nm1.eq.0) then
            if(tBetaExcit) then
                call CountSizeCompMat(Nm1FCIDetList,Elec-1,nNm1FCIDet,Nmax_Nm1,Nm1bBitList)
            else
                call CountSizeCompMat(Nm1bFCIDetList,Elec-1,nNm1bFCIDet,Nmax_Nm1,Nm1bBitList)
            endif
            iHamSize_Nm1 = Nmax_Nm1
        else
            write(6,"(A,I16)") "Assuming size of compressed N-1 electron hamiltonian is: ",iHamSize_Nm1
            Nmax_Nm1 = iHamSize_Nm1
        endif
        if(iHamSize_Np1.eq.0) then
            if(tBetaExcit) then
                call CountSizeCompMat(Np1bFCIDetList,Elec+1,nNp1bFCIDet,Nmax_Np1,Np1bBitList)
            else
                call CountSizeCompMat(Np1FCIDetList,Elec+1,nNp1FCIDet,Nmax_Np1,Np1BitList)
            endif
            iHamSize_Np1 = Nmax_Np1
        else
            write(6,"(A,I16)") "Assuming size of compressed N+1 electron hamiltonian is: ",iHamSize_Np1
            Nmax_Np1 = iHamSize_Np1
        endif

        write(6,"(A,F12.5,A)") "Memory required for N-electron hamil: ",(real(Nmax_N,dp)*2)*ComptoMb, " Mb"
        write(6,"(A,F12.5,A)") "Memory required for N+1-electron hamil: ",(real(Nmax_Np1,dp)*2)*ComptoMb, " Mb"
        write(6,"(A,F12.5,A)") "Memory required for N-1-electron hamil: ",(real(Nmax_Nm1,dp)*2)*ComptoMb, " Mb"
        call flush(6)
        allocate(NFCIHam_cmps(Nmax_N),stat=ierr)
        allocate(NFCIHam_inds(Nmax_N),stat=ierr)
        !Want Np1b for particle and Nm1 for hole if tBetaExcit, or vice versa otherwise
        !We have removed the distinction, and so will have to remember...
        allocate(Np1FCIHam_cmps(Nmax_Np1),stat=ierr)
        allocate(Np1FCIHam_inds(Nmax_Np1),stat=ierr)
        allocate(Nm1FCIHam_cmps(Nmax_Nm1),stat=ierr)
        allocate(Nm1FCIHam_inds(Nmax_Nm1),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Allocation fail')
        !Want Np1b for particle and Nm1 for hole if tBetaExcit, or vice versa otherwise
        !We have removed the distinction, and so will have to remember...
        call Fill_N_Np1_Nm1b_FCIHam(Elec,Nmax_N,Nmax_Np1,Nmax_Nm1,NHam_cmps=NFCIHam_cmps,Np1Ham_cmps=Np1FCIHam_cmps, &
            Nm1Ham_cmps=Nm1FCIHam_cmps,NHam_inds=NFCIHam_inds,Np1Ham_inds=Np1FCIHam_inds,   &
            Nm1Ham_inds=Nm1FCIHam_inds,tSwapExcits = tBetaExcit)

        if(tBetaExcit) then
            nLinearSystem = nNp1bFCIDet + nFCIDet
            nLinearSystem_h = nNm1FCIDet + nFCIDet
        else
            nLinearSystem = nNp1FCIDet + nFCIDet
            nLinearSystem_h = nNm1bFCIDet + nFCIDet
        endif

        if(nLinearSystem.ne.nLinearSystem_h) then
            call stop_all(t_r,'Should change restriction on the Linear system size being the same for particle/hole addition')
        endif

!        if(tWriteOut) then
!            call writevectorcomp(NFCIHam_cmps,'NFCIHam_cmps')
!            call writevectorcomp(Np1FCIHam_cmps,'Np1FCIHam_cmps')
!        endif
        
        !If doing full optimization of the GS problem
        nGSSpace = nFCIDet + nNp1FCIDet + nNm1bFCIDet
        Np1GSInd = nFCIDet + 1
        Nm1GSInd = np1GSInd + nNp1FCIDet
        
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
        if(tSC_LR) then
            write(iunit,"(A)") "# 1.Frequency     2.GF_LinearResponse(Re)    3.GF_LinearResponse(Im)    " &
                & //"4.ParticleGF(Re)   5.ParticleGF(Im)   6.HoleGF(Re)   7.HoleGF(Im)    8.Old_GS    9.New_GS   " &
                & //"10.Particle_Norm  11.Hole_Norm   12.NI_GF(Re)   13.NI_GF(Im)  14.NI_GF_Part(Re)   15.NI_GF_Part(Im)   " &
                & //"16.NI_GF_Hole(Re)   17.NI_GF_Hole(Im)  18.SpectralWeight  19.Iters_p   20.Iters_h    " &
                & //"21.Tr[SE](Re)  22.Tr[SE](Im) "
        else
            write(iunit,"(A)") "# 1.Frequency     2.GF_LinearResponse(Re)    3.GF_LinearResponse(Im)    " &
                & //"4.ParticleGF(Re)   5.ParticleGF(Im)   6.HoleGF(Re)   7.HoleGF(Im)    8.Old_GS    9.New_GS   " &
                & //"10.Particle_Norm  11.Hole_Norm   12.NI_GF(Re)   13.NI_GF(Im)  14.NI_GF_Part(Re)   15.NI_GF_Part(Im)   " &
                & //"16.NI_GF_Hole(Re)   17.NI_GF_Hole(Im)  18.SpectralWeight  19.Iters_p   20.Iters_h    " 
        endif
                
        if(tMinRes_NonDir) then
            minres_unit = get_free_unit()
            open(minres_unit,file='zMinResQLP.txt',status='unknown',position='append')
            allocate(RHS(nLinearSystem))
            if(tPrecond_MinRes) then
                allocate(Precond_Diag(nLinearSystem))
            endif
        else
            minres_unit = get_free_unit()
            open(minres_unit,file='zGMRES.txt',status='unknown',position='append')
        endif
        
        if(tAllImp_LR) then
            !We calculate all greens functions between all impurity sites - both in the non-interacting and interacting case
            nImp_GF = nImp 
        else
            nImp_GF = 1
        endif
        allocate(SchmidtPertGF_Cre_Ket(VirtStart:nSites,nImp_GF))
        allocate(SchmidtPertGF_Cre_Bra(VirtStart:nSites,nImp_GF))
        allocate(SchmidtPertGF_Ann_Ket(1:CoreEnd,nImp_GF))
        allocate(SchmidtPertGF_Ann_Bra(1:CoreEnd,nImp_GF))
        allocate(NI_LRMat_Cre(nImp_GF,nImp_GF))
        allocate(NI_LRMat_Ann(nImp_GF,nImp_GF))
        allocate(ni_lr_Mat(nImp_GF,nImp_GF))
        allocate(ResponseFn_p(nImp_GF,nImp_GF))
        allocate(ResponseFn_h(nImp_GF,nImp_GF))
        allocate(ResponseFn_Mat(nImp_GF,nImp_GF))
        allocate(dNorm_p(nImp_GF))
        allocate(dNorm_h(nImp_GF))
        
        !Allocate memory for hamiltonian in this system:
        allocate(Psi_0(nGSSpace),stat=ierr)   !If tLR_ReoptGS = .F. , then only the first nFCIDet elements will be used
        allocate(Psi1_p(nLinearSystem),stat=ierr)
        allocate(Psi1_h(nLinearSystem),stat=ierr)
        Psi1_p = zzero
        Psi1_h = zzero
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        i = nGSSpace + nLinearSystem*3
        write(6,"(A,F12.5,A)") "Memory required for wavefunctions: ",real(i,dp)*ComptoMb, " Mb"
        
        !Since this does not depend at all on the new basis, since the ground state has a completely
        !disjoint basis to |1>, the RHS of the equations can be precomputed
        allocate(Cre_0(nLinearSystem,nImp_GF))
        allocate(Ann_0(nLinearSystem,nImp_GF))
        if(.not.tLR_ReoptGS) then
            Psi_0(:) = zzero 
            Psi_0(1:nFCIDet) = HL_Vec(:)    
            GFChemPot = HL_Energy
            do i = 1,nImp_GF
                call ApplySP_PertGS_EC(Psi_0,nGSSpace,Cre_0(:,i),Ann_0(:,i),nLinearSystem,i,tSwapExcits=tBetaExcit)
            enddo
        endif
        
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
        if(tBetaExcit) then
            FockSchmidt(nOcc-nImp+1:nOcc+nImp,nOcc-nImp+1:nOcc+nImp) = Emb_h0v_b(:,:)
        else
            FockSchmidt(nOcc-nImp+1:nOcc+nImp,nOcc-nImp+1:nOcc+nImp) = Emb_h0v(:,:)
        endif
        !If we are doing self-consistent linear response, with a self-energy in the HL part, 
        !then this is calculated later, since it is now omega and iteration dependent
        do i = 1,nSites
            do j = 1,nSites
                if(tBetaExcit) then
                    FockSchmidt_SE(j,i) = dcmplx(FockSchmidt_b(j,i),0.0_dp)
                else
                    FockSchmidt_SE(j,i) = dcmplx(FockSchmidt(j,i),0.0_dp)
                endif
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
        allocate(Gc_a_F_ax_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Gc_a_F_ax_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Gc_b_F_ab(VirtStart:VirtEnd,nImp_GF),stat=ierr)
        !For hole addition
        allocate(Ga_i_F_xi_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Ga_i_F_xi_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Ga_i_F_ij(1:nCore,nImp_GF),stat=ierr)

        !Allocate and precompute 1-operator coupling coefficients between the different sized spaces.
        !First number is the index of operator, and the second is the parity change when applying the operator

        !First, determine the size of this matrix in compressed form
        Nmax_Coup_Create = 0
        Nmax_Coup_Ann = 0
        do J = 1,nFCIDet
            !Start with the creation of a spin orbital
            do K = 1,nNm1bFCIDet
                DiffOrb = ieor(Nm1bBitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(nOrbs.eq.1) then
                    Nmax_Coup_Create = Nmax_Coup_Create + 1
                endif
            enddo
            !Now, annihilation of a spin orbital
            do K = 1,nNp1FCIDet
                DiffOrb = ieor(Np1BitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(nOrbs.eq.1) then
                    Nmax_Coup_Ann = Nmax_Coup_Ann + 1
                endif
            enddo
        enddo
        
        allocate(Coup_Create(2,Nmax_Coup_Create))
        allocate(Coup_Ann(2,Nmax_Coup_Ann))
        !Also indexing arrays
        allocate(Coup_Create_inds(Nmax_Coup_Create))
        allocate(Coup_Ann_inds(Nmax_Coup_Ann))
        !_cum arrays index start of row block
        allocate(Coup_Create_cum(nFCIDet+1))
        allocate(Coup_Ann_cum(nFCIDet+1))

        !Unfortunately, it is not easy to calculate the transpose. Do it explicitly now (matrices are at least of the same size)
        allocate(Coup_Create_T(2,Nmax_Coup_Create))
        allocate(Coup_Ann_T(2,Nmax_Coup_Ann))
        !Also indexing arrays
        allocate(Coup_Create_inds_T(Nmax_Coup_Create))
        allocate(Coup_Ann_inds_T(Nmax_Coup_Ann))
        allocate(Coup_Create_cum_T(nNm1bFCIDet+1))
        allocate(Coup_Ann_cum_T(nNp1FCIDet+1))
        
        !Create coupling matrices of the appropriate spin
        call CreateCmprsCouplingMatrices(NMax_Coup_Create,NMax_Coup_Ann,Coup_Create,Coup_Ann,Coup_Create_inds,  &
            Coup_Ann_inds,Coup_Create_cum,Coup_Ann_cum,Coup_Create_T,Coup_Ann_T,Coup_Create_inds_T, &
            Coup_Ann_inds_T,Coup_Create_cum_T,Coup_Ann_cum_T,tSwapSpins=tBetaExcit)

        !Now allocate the compressed matrices for the linear systems.
        Nmax_Lin_p = Nmax_Np1 + Nmax_N + 2*Nmax_Coup_Ann
        Nmax_Lin_h = Nmax_Nm1 + Nmax_N + 2*Nmax_Coup_Create
        write(6,"(A,F14.6,A)") "Memory required for compressed LR systems: ",   &
            (3*(real(Nmax_Lin_p,dp)*8)+3*(real(Nmax_Lin_h,dp)*8))/1048576.0_dp," Mb"
        write(6,"(A,F14.6,A)") "Compare to memory required for the uncompressed system: ",  &
            2*(real(nLinearSystem,dp)**2)*16/1048576.0_dp," Mb"
        allocate(LinearSystem_p(Nmax_Lin_p))
        allocate(LinearSystem_h(Nmax_Lin_h))
        allocate(LinearSystem_p_inds(Nmax_Lin_p))
        allocate(LinearSystem_h_inds(Nmax_Lin_h))

        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        call halt_timer(LR_EC_GF_Precom)

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            call set_timer(LR_EC_GF_HBuild)
        
            write(6,*) "Calculating linear response for frequency: ",Omega
            write(minres_unit,*) "Iteratively solving for frequency: ",Omega

            !First, find the non-interacting solution expressed in the schmidt basis
            !This will only calculate it over the first impurity site
            call FindSchmidtPert_Charged(Omega+mu,NI_LRMat_Cre,NI_LRMat_Ann,tSwapSpins=tBetaExcit)

            !Construct useful intermediates, using zgemm
            !sum_a Gc_a^* F_ax (Creation)
            call ZGEMM('T','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_VX(VirtStart:VirtEnd,ActiveStart:ActiveEnd),  &
                nVirt,SchmidtPertGF_Cre_Bra(VirtStart:VirtEnd,:),nVirt,zzero,Gc_a_F_ax_Bra,EmbSize)
            call ZGEMM('N','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:VirtEnd),  &
                EmbSize,SchmidtPertGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Gc_a_F_ax_Ket,EmbSize)
            !sum_b Gc_b F_ab  (Creation)
            call ZGEMM('N','N',nVirt,nImp_GF,nVirt,zone,FockSchmidt_SE_VV(VirtStart:VirtEnd,VirtStart:VirtEnd),    &
                nVirt,SchmidtPertGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Gc_b_F_ab,nVirt)

            !sum_i Ga_i^bra F_xi (Annihilation)
            call ZGEMM('N','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd),EmbSize,    &
                SchmidtPertGF_Ann_Bra(1:nCore,:),nCore,zzero,Ga_i_F_xi_Bra,EmbSize)
!                call writematrixcomp(FockSchmidt_SE_XC,'FockSchmidt_SE_XC',.false.)
!                call writevectorcomp(SchmidtPertGF_Ann_Bra(:,1),'SchmidtPertGF_Ann_Bra')
!                call writevectorcomp(Ga_i_F_xi_Bra(:,1),'Ga_i_F_xi_Bra')
            call ZGEMM('T','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd),nCore,    &
                SchmidtPertGF_Ann_Ket(1:nCore,:),nCore,zzero,Ga_i_F_xi_Ket,EmbSize)
!                call writevectorcomp(Ga_i_F_xi_Bra(:,1),'Ga_i_F_xi_Ket')
            !sum_i Ga_i F_ij (Annihilation)
            call ZGEMM('T','N',nCore,nImp_GF,nCore,zone,FockSchmidt_SE_CC(:,:),nCore,SchmidtPertGF_Ann_Ket(:,:),nCore,zzero, &
                Ga_i_F_ij,nCore)

            !Find the first-order wavefunction in the interacting picture for all impurity sites 
            do pertsite = 1,nImp_GF

                LinearSystem_p(:) = zzero 
                LinearSystem_h(:) = zzero 
                LinearSystem_p_inds(:) = 0
                LinearSystem_h_inds(:) = 0

                VNorm = zzero
                tempel = zzero 
                do a = VirtStart,VirtEnd
                    !Calc normalization for the CV block
                    VNorm = VNorm + SchmidtPertGF_Cre_Bra(a,pertsite)*SchmidtPertGF_Cre_Ket(a,pertsite) 
                    !Calculate the diagonal correction
                    !tempel = tempel + conjg(SchmidtPertGF_Cre(a,pertsite))*Gc_b_F_ab(a,pertsite)
                    tempel = tempel + SchmidtPertGF_Cre_Bra(a,pertsite)*Gc_b_F_ab(a,pertsite)
                enddo
                tempel = tempel / VNorm

                !First, store diagonals
                LinearSystem_p(1:nNp1FCIDet) = Np1FCIHam_cmps(1:nNp1FCIDet)
                j = 0
                do i = nNp1FCIDet+1,nLinearSystem
                    j = j + 1
                    LinearSystem_p(i) = NFCIHam_cmps(j) + tempel
                enddo
                if(j.ne.nFCIDet) call stop_all(t_r,'Error with matrix indexing')

                !And for hole diagonals
                CNorm = zzero
                tempel = zzero 
                do i = 1,CoreEnd
                    !Calc normalization
                    CNorm = CNorm + SchmidtPertGF_Ann_Bra(i,pertsite)*SchmidtPertGF_Ann_Ket(i,pertsite)
                    !Calc diagonal correction
                    tempel = tempel + SchmidtPertGF_Ann_Bra(i,pertsite)*Ga_i_F_ij(i,pertsite)
                enddo
                tempel = -tempel / CNorm
                
                !Store diagonals
                LinearSystem_h(1:nNm1bFCIDet) = Nm1FCIHam_cmps(1:nNm1bFCIDet)
                j = 0
                do i = nNm1bFCIDet + 1,nLinearSystem
                    j = j + 1
                    LinearSystem_h(i) = NFCIHam_cmps(j) + tempel
                enddo
                if(j.ne.nFCIDet) call stop_all(t_r,'Error with matrix indexing')

                !Now for off-diagonal elements
                ind_h = nLinearSystem + 1
                LinearSystem_h_inds(1) = nLinearSystem + 2
                ind_p = nLinearSystem + 1
                LinearSystem_p_inds(1) = nLinearSystem + 2
                
                !Deal with particle matrix
                !First - the top rows (np1 rows)
                do i = 1,nNp1FCIDet
                    !Loop through rows - first through block 1
                    do k = Np1FCIHam_inds(i),Np1FCIHam_inds(i+1)-1
                        !k indexes the non-zero, off-diagonal matrix elements moving along the row of i
                        if(abs(Np1FCIHam_cmps(k)).ge.CompressThresh) then
                            ind_p = ind_p + 1
                            if(ind_p.gt.Nmax_Lin_p) call stop_all(t_r,'Compressed array p size too small')
                            LinearSystem_p(ind_p) = Np1FCIHam_cmps(k)
                            !But what determinant do these correspond to? This is given by Np1FCIHam_alpha_inds(k)
                            LinearSystem_p_inds(ind_p) = Np1FCIHam_inds(k)
                        endif
                    enddo
                    !Now loop through the same row of block 3
                    do k = Coup_Ann_cum_T(i),Coup_Ann_cum_T(i+1)-1
                        matel = Gc_a_F_ax_Ket(Coup_Ann_T(1,k),pertsite)*Coup_Ann_T(2,k)/sqrt(VNorm)
                        if(abs(matel).ge.CompressThresh) then
                            ind_p = ind_p + 1
                            if(ind_p.gt.Nmax_Lin_p) call stop_all(t_r,'Compressed array p size too small 2')
                            LinearSystem_p(ind_p) = matel 
                            LinearSystem_p_inds(ind_p) = Coup_Ann_inds_T(k) + nNp1FCIDet  !Add on cols from block 1
                        endif
                    enddo
                    LinearSystem_p_inds(i+1) = ind_p + 1
                enddo
                !Now, the rows from block 3_T and 2
                do i = 1,nFCIDet
                    do k = Coup_Ann_cum(i),Coup_Ann_cum(i+1)-1
                        matel = Gc_a_F_ax_Bra(Coup_Ann(1,k),pertsite)*Coup_Ann(2,k)/sqrt(VNorm)
                        if(abs(matel).gt.CompressThresh) then
                            ind_p = ind_p + 1
                            if(ind_p.gt.Nmax_Lin_p) call stop_all(t_r,'Compressed array p size too small 3')
                            LinearSystem_p(ind_p) = matel 
                            LinearSystem_p_inds(ind_p) = Coup_Ann_inds(k)
                        endif
                    enddo
                    !Now loop through the same row of block 2
                    do k = NFCIHam_inds(i),NFCIHam_inds(i+1)-1
                        if(abs(NFCIHam_cmps(k)).gt.CompressThresh) then
                            ind_p = ind_p + 1
                            if(ind_p.gt.Nmax_Lin_p) call stop_all(t_r,'Compressed array p size too small 4')
                            LinearSystem_p(ind_p) = NFCIHam_cmps(k)
                            LinearSystem_p_inds(ind_p) = NFCIHam_inds(k) + nNp1FCIDet   !Add on cols from block 3^T
                        endif
                    enddo
                    LinearSystem_p_inds(i+1+nNp1FCIDet) = ind_p + 1
                enddo

                !Now for hole hamiltonian
                do i = 1,nNm1bFCIDet
                    !Loop through rows
                    do k = Nm1FCIHam_inds(i),Nm1FCIHam_inds(i+1)-1
                        if(abs(Nm1FCIHam_cmps(k)).ge.CompressThresh) then
                            ind_h = ind_h + 1
                            if(ind_h.gt.Nmax_Lin_h) call stop_all(t_r,'Compressed array h size too small')
                            LinearSystem_h(ind_h) = Nm1FCIHam_cmps(k)
                            LinearSystem_h_inds(ind_h) = Nm1FCIHam_inds(k)
                        endif
                    enddo
                    do k = Coup_Create_cum_T(i),Coup_Create_cum_T(i+1)-1
                        matel = - Ga_i_F_xi_Ket(Coup_Create_T(1,k),pertsite)*Coup_Create_T(2,k)/sqrt(CNorm)
                        if(abs(matel).gt.CompressThresh) then
                            ind_h = ind_h + 1
                            if(ind_h.gt.Nmax_Lin_p) call stop_all(t_r,'Compressed array h size too small 2')
                            LinearSystem_h(ind_h) = matel 
                            LinearSystem_h_inds(ind_h) = Coup_Create_inds_T(k) + nNm1bFCIDet
                        endif
                    enddo
                    LinearSystem_h_inds(i+1) = ind_h + 1
                enddo
                !Now, the rows from block 3^T and 2
                do i = 1,nFCIDet
                    do k = Coup_Create_cum(i),Coup_Create_cum(i+1)-1
                        matel = - Ga_i_F_xi_Bra(Coup_Create(1,k),pertsite)*Coup_Create(2,k)/sqrt(CNorm)
                        if(abs(matel).gt.CompressThresh) then
                            ind_h = ind_h + 1
                            if(ind_h.gt.Nmax_Lin_h) call stop_all(t_r,'Compressed array h size too small 3')
                            LinearSystem_h(ind_h) = matel 
                            LinearSystem_h_inds(ind_h) = Coup_Create_inds(k)
                        endif
                    enddo
                    !Now loop through the same row of block 2
                    do k = NFCIHam_inds(i),NFCIHam_inds(i+1)-1
                        if(abs(NFCIHam_cmps(k)).gt.CompressThresh) then
                            ind_h = ind_h + 1
                            if(ind_h.gt.Nmax_Lin_h) call stop_all(t_r,'Compressed array h size too small 4')
                            LinearSystem_h(ind_h) = NFCIHam_cmps(k)
                            LinearSystem_h_inds(ind_h) = NFCIHam_inds(k) + nNm1bFCIDet
                        endif
                    enddo
                    LinearSystem_h_inds(i+1+nNm1bFCIDet) = ind_h + 1
                enddo

                !write(6,*) "Hessian constructed successfully...",Omega
                call halt_timer(LR_EC_GF_HBuild)
                call set_timer(LR_EC_GF_SolveLR)

                !Solve particle GF to start
                !The V|0> for particle and hole perturbations are held in Cre_0 and Ann_0
                !call writevectorcomp(Psi_0,'Psi_0')
                !call writevectorcomp(Cre_0(:,1),'Ann_0')
                !call writevectorcomp(Ann_0,'Ann_0')
                
                LinearSystem_p(:) = -LinearSystem_p(:)
                !Offset matrix
                do i = 1,nLinearSystem
                    LinearSystem_p(i) = LinearSystem_p(i) + dcmplx(Omega+mu+GFChemPot,dDelta)
                enddo

                !Now solve these linear equations
                !call writevectorcomp(Psi1_p,'Cre_0')
                zDirMV_Mat_cmprs => LinearSystem_p
                zDirMV_Mat_cmprs_inds => LinearSystem_p_inds
                if(tMinRes_NonDir) then
                    call setup_RHS(nLinearSystem,Cre_0(:,pertsite),RHS)
                    maxminres_iter_ip = int(maxminres_iter,ip)
                    minres_unit_ip = int(minres_unit,ip)
                    nLinearSystem_ip = int(nLinearSystem,ip)
                    if(tPrecond_MinRes) then
                        call FormPrecond(nLinearSystem)
                        call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_p, &
                            itnlim=maxminres_iter_ip,Msolve=zPreCond,istop=info_ip,rtol=rtol_LR,itn=iters_p, &
                            startguess=tReuse_LS)
                    else
                        call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_p, &
                            itnlim=maxminres_iter_ip,istop=info_ip,rtol=rtol_LR,itn=iters_p,startguess=tReuse_LS)
                    endif
                    info = info_ip
                    if(info.gt.7) write(6,*) "info: ",info
                    if(info.eq.8) write(6,"(A,I9)") "Linear equation solver hit iteration limit: ",iters_p
                    if((info.eq.9).or.(info.eq.10).or.(info.eq.11)) then
                        call stop_all(t_r,'Input matrices to linear solver incorrect')
                    endif
                    if(info.gt.11) call stop_all(t_r,'Linear equation solver failed')
                elseif(tGMRES_NonDir) then
                    call GMRES_Solve(nLinearSystem,Cre_0(:,pertsite),minres_unit,maxminres_iter,rtol_LR,tPrecond_MinRes,    &
                        tReuse_LS,iters_p,Psi1_p,info)
                endif
                zDirMV_Mat_cmprs => null() 
                zDirMV_Mat_cmprs_inds => null() 
                !call writevectorcomp(Psi1_p,'Psi1_p')

                !Now solve the LR for the hole addition
                do i = 1,nLinearSystem
                    LinearSystem_h(i) = dcmplx(Omega+mu,dDelta) + (LinearSystem_h(i) - dcmplx(GFChemPot,0.0_dp))
                enddo
                zDirMV_Mat_cmprs => LinearSystem_h
                zDirMV_Mat_cmprs_inds => LinearSystem_h_inds
                if(tMinRes_NonDir) then
                    call setup_RHS(nLinearSystem,Ann_0(:,pertsite),RHS)
                    maxminres_iter_ip = int(maxminres_iter,ip)
                    minres_unit_ip = int(minres_unit,ip)
                    nLinearSystem_ip = int(nLinearSystem,ip)
                    if(tPrecond_MinRes) then
                        call FormPrecond(nLinearSystem)
                        call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_h, &
                            itnlim=maxminres_iter_ip,Msolve=zPreCond,istop=info_ip,rtol=rtol_LR,itn=iters_h, &
                            startguess=tReuse_LS)
                    else
                        call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_h, &
                            itnlim=maxminres_iter_ip,istop=info_ip,rtol=rtol_LR,itn=iters_h,startguess=tReuse_LS)
                    endif
                    info = info_ip
                    if(info.gt.7) write(6,*) "info: ",info
                    if(info.eq.8) write(6,"(A,I9)") "Linear equation solver hit iteration limit: ",iters_h
                    if((info.eq.9).or.(info.eq.10).or.(info.eq.11)) then
                        call stop_all(t_r,'Input matrices to linear solver incorrect')
                    endif
                    if(info.gt.11) call stop_all(t_r,'Linear equation solver failed')
                elseif(tGMRES_NonDir) then
                    call GMRES_Solve(nLinearSystem,Ann_0(:,pertsite),minres_unit,maxminres_iter,rtol_LR,tPrecond_MinRes,    &
                        tReuse_LS,iters_h,Psi1_h,info)
                endif
        
                zDirMV_Mat_cmprs => null()
                zDirMV_Mat_cmprs_inds => null()
                !Find normalization of first-order wavefunctions
                dNorm_p(pertsite) = znrm2(nLinearSystem,Psi1_p,1)
                dNorm_h(pertsite) = znrm2(nLinearSystem,Psi1_h,1)
                !write(6,*) "Normalization of first-order wavefunction: ",dNorm_p,dNorm_h
                
                !Now, calculate all interacting greens funtions that we want, by dotting in the first-order space
                if((nImp_GF.gt.1).and.(tLR_ReoptGS)) then
                    call stop_all(t_r,'Cannot do dotting of vectors here, since the LHS vectors are still being computed')
                endif
                do j = 1,nImp_GF
                    ResponseFn_p(pertsite,j) = zdotc(nLinearSystem,Cre_0(:,j),1,Psi1_p,1)
                    ResponseFn_h(pertsite,j) = zdotc(nLinearSystem,Ann_0(:,j),1,Psi1_h,1)
                    ResponseFn_Mat(pertsite,j) = ResponseFn_p(pertsite,j) + ResponseFn_h(pertsite,j)
                    ni_lr_Mat(pertsite,j) = NI_LRMat_Cre(pertsite,j) + NI_LRMat_Ann(pertsite,j) 
                enddo

            enddo   !End do over pertsite. We now have all the greens funtions 

            !Now, calculate NI and interacting response function as trace over diagonal parts of the local greens functions
            ni_lr = zzero
            ResponseFn = zzero
            AvResFn_p = zzero
            AvResFn_h = zzero
            ni_lr_p = zzero
            ni_lr_h = zzero
            AvdNorm_p = zero
            AvdNorm_h = zero
            do i = 1,nImp_GF
                ni_lr = ni_lr + ni_lr_Mat(i,i)
                ResponseFn = ResponseFn + ResponseFn_Mat(i,i)
                AvResFn_p = AvResFn_p + ResponseFn_p(i,i)
                AvResFn_h = AvResFn_h + ResponseFn_h(i,i)
                ni_lr_p = ni_lr_p + NI_LRMat_Cre(i,i)
                ni_lr_h = ni_lr_h + NI_LRMat_Ann(i,i)
                AvdNorm_p = AvdNorm_p + dNorm_p(i)
                AvdNorm_h = AvdNorm_h + dNorm_h(i)
            enddo
            ni_lr = ni_lr/real(nImp_GF,dp)
            ni_lr_p = ni_lr_p/real(nImp_GF,dp)
            ni_lr_h = ni_lr_h/real(nImp_GF,dp)
            ResponseFn = ResponseFn/real(nImp_GF,dp)
            AvResFn_p = AvResFn_p/real(nImp_GF,dp)
            AvResFn_h = AvResFn_h/real(nImp_GF,dp)
            AvdNorm_p = AvdNorm_p/real(nImp_GF,dp)
            AvdNorm_h = AvdNorm_h/real(nImp_GF,dp)

            Diff_GF = zero  !real(abs(ResponseFn - ni_lr))
            do i=1,nImp_GF
                do j=1,nImp_GF
                    Diff_GF = Diff_GF + real((ni_lr_Mat(j,i)-ResponseFn_Mat(j,i))*dconjg(ni_lr_Mat(j,i)-ResponseFn_Mat(j,i)))
                enddo
            enddo
        
            call halt_timer(LR_EC_GF_SolveLR)

            if(.not.tFirst) then
                SpectralWeight = SpectralWeight + Omega_Step*(Prev_Spec-aimag(ResponseFn))/(2.0_dp*pi)
                Prev_Spec = -aimag(ResponseFn)
            endif

            write(iunit,"(18G22.10,2I7)") Omega,real(ResponseFn),-aimag(ResponseFn), &
                real(AvResFn_p),-aimag(AvResFn_p),real(AvResFn_h),-aimag(AvResFn_h),    &
                HL_Energy,GFChemPot,abs(AvdNorm_p),abs(AvdNorm_h),real(ni_lr),-aimag(ni_lr),real(ni_lr_p),    &
                -aimag(ni_lr_p),real(ni_lr_h),-aimag(ni_lr_h),SpectralWeight,iters_p,iters_h
            call flush(iunit)
            call flush(6)

            if(tFirst) tFirst = .false.

            Omega = Omega + Omega_Step
        enddo   !End loop over omega

        write(6,"(A,G22.10)") "Total integrated spectral weight: ",SpectralWeight

        deallocate(LinearSystem_p,LinearSystem_h,LinearSystem_p_inds,LinearSystem_h_inds,Psi1_p,Psi1_h)
        deallocate(Cre_0,Ann_0,Psi_0,SchmidtPertGF_Cre_Ket,SchmidtPertGF_Ann_Ket)
        deallocate(NI_LRMat_Cre,NI_LRMat_Ann,ResponseFn_p,ResponseFn_h,ResponseFn_Mat)
        deallocate(ni_lr_Mat,SchmidtPertGF_Cre_Bra,SchmidtPertGF_Ann_Bra)
        close(iunit)
        if(tMinRes_NonDir) then
            deallocate(RHS)
            if(tPrecond_MinRes) deallocate(Precond_Diag)
        endif
        if(tGMRES_NonDir.or.tMinRes_NonDir) close(minres_unit)

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

        !Stored intermediates
        deallocate(NFCIHam_cmps,Nm1FCIHam_cmps,Np1FCIHam_cmps)
        deallocate(NFCIHam_inds,Nm1FCIHam_inds,Np1FCIHam_inds)
        deallocate(Coup_Create,Coup_Ann,Coup_Create_inds,Coup_Ann_inds)
        deallocate(Coup_Create_cum,Coup_Ann_cum,Coup_Create_T,Coup_Ann_T)
        deallocate(Coup_Create_inds_T,Coup_Ann_inds_T,Coup_Create_cum_T,Coup_Ann_cum_T)
        deallocate(FockSchmidt_SE,FockSchmidt_SE_VV,FockSchmidt_SE_CC)
        deallocate(FockSchmidt_SE_VX,FockSchmidt_SE_CX,FockSchmidt_SE_XV,FockSchmidt_SE_XC)
        deallocate(Gc_a_F_ax_Bra,Gc_a_F_ax_Ket,Gc_b_F_ab,Ga_i_F_xi_Bra,Ga_i_F_xi_Ket,Ga_i_F_ij)

    end subroutine NonIntExCont_TDA_MCLR_Charged_Cmprs

    !Calculate linear response for charged excitations - add the hole creation to particle creation
    subroutine NonIntExCont_TDA_MCLR_Charged()
        use mat_tools, only: add_localpot_comp_inplace
        use utils, only: get_free_unit,append_ext_real,append_ext
        use DetBitOps, only: DecodeBitDet,SQOperator,CountBits
        use Davidson, only: Comp_NonDir_Davidson
        use fitting, only: Fit_SE
        use zminresqlpModule, only: MinresQLP  
        use DetToolsData
        use DetTools, only: tospat,umatind,gendets
        use solvers, only: CreateIntMats
        implicit none
        complex(dp), allocatable :: NFCIHam(:,:),Np1FCIHam_alpha(:,:),Nm1FCIHam_beta(:,:)
        real(dp), allocatable :: W(:),dNorm_p(:),dNorm_h(:)
        complex(dp), allocatable , target :: LinearSystem_p(:,:),LinearSystem_h(:,:)
        complex(dp), allocatable :: Cre_0(:,:),Ann_0(:,:),ResponseFn_p(:,:),ResponseFn_h(:,:)
        complex(dp), allocatable :: Gc_a_F_ax_Bra(:,:),Gc_a_F_ax_Ket(:,:),Gc_b_F_ab(:,:),GSHam(:,:)
        complex(dp), allocatable :: ResponseFn_Mat(:,:),Ga_i_F_xi_Ket(:,:)
        complex(dp), allocatable :: Psi1_p(:),Psi1_h(:),Ga_i_F_xi_Bra(:,:),Ga_i_F_ij(:,:),ni_lr_Mat(:,:)
        complex(dp), allocatable :: temp_vecc(:),Work(:),Psi_0(:),RHS(:),SE_Change(:,:)
        complex(dp), allocatable :: NI_LRMat_Cre(:,:),NI_LRMat_Ann(:,:)
        integer, allocatable :: Coup_Ann_alpha(:,:,:),Coup_Create_alpha(:,:,:)
        integer :: i,a,j,k,ActiveEnd,ActiveStart,CoreEnd,DiffOrb,gam,ierr,info,iunit,nCore
        integer :: nLinearSystem,nOrbs,nVirt,tempK,VIndex,VirtStart,VirtEnd
        integer :: orbdum(1),nLinearSystem_h,nGSSpace,Np1GSInd,Nm1GSInd,lWork,minres_unit
        integer :: maxminres_iter,nImp_GF,pertsite,SE_Fit_Iter,nNR_Iters
        integer(ip) :: nLinearSystem_ip,minres_unit_ip,info_ip,maxminres_iter_ip,iters_p,iters_h
        real(dp) :: Omega,GFChemPot,mu,SpectralWeight,Prev_Spec,AvdNorm_p,AvdNorm_h,Var_SE,Error_GF
        real(dp) :: Diff_GF
        complex(dp) :: ResponseFn,tempel,ni_lr,ni_lr_p,ni_lr_h,AvResFn_p,AvResFn_h
        complex(dp) :: zdotc,VNorm,CNorm,Prev_SE_Saved(nImp,nImp),TR_SE
        logical :: tParity,tFirst,tSCFConverged,tSCFFailed
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='NonIntExCont_TDA_MCLR_Charged'

        maxminres_iter = iMinRes_MaxIter
        iters_p = 0
        iters_h = 0

        SpectralWeight = 0.0_dp
        Prev_Spec = 0.0_dp
        tFirst = .true.

        call set_timer(LR_EC_GF_Precom)

        write(6,"(A)") "Calculating non-interacting EC MR-TDA LR system for *hole* & *particle* alpha spin-orbital perturbations..."
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(tMinRes_NonDir) then
            if(tPreCond_MinRes) then
                write(6,"(A)") "Solving linear system with iterative non-direct preconditioned MinRes-QLP algorithm"
            else
                write(6,"(A)") "Solving linear system with iterative non-direct MinRes-QLP algorithm"
            endif
            write(6,"(A,G22.10)") "Tolerance for solution of linear system: ",rtol_LR
            write(6,"(A,G22.10)") "Maximum iterations for each solution: ",maxminres_iter
        else
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
        endif
        if(tLR_ReoptGS) then
            if(tNonDirDavidson) then
                write(6,"(A)") "Reoptimizing ground state eigenfunction in the space of V*|1> " &
     &              //"at each point with non-direct davidson"
            else
                write(6,"(A)") "Reoptimizing ground state eigenfunction in the space of V*|1> " &
     &              //"at each point with complete diagonalizer"
            endif
        endif
        if(tSC_LR.and.(mod(nOcc,nImp).ne.0)) then
            ! This criteria must be satisfied, because otherwise, the orbitals resulting from the diagonalization of
            ! h + Self_energy will mix the occupied and virtual spaces. This will result in the NI greens function
            ! spanning a different space to previous iterations
            call stop_all(t_r,'Cannot tile the self-energy through the occupied and virtual space seperately without mixing them')
        endif
        !Surely if we are doing self-consistent LR, we can only have the GS in the n-electron space, since the others
        !aren't consistent across the different greens functions.
        !Would be ok with nImp = 1 I guess, since there is still then a unique GS.
        if(tSC_LR.and.tLR_ReoptGS) then
            call stop_all(t_r,"Reoptimizing GS with self-consistency not sorted yet - probably shouldn't happen")
        endif

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
        call Fill_N_Np1_Nm1b_FCIHam(Elec,1,1,1,NHam=NFCIHam,Np1Ham=Np1FCIHam_alpha,Nm1Ham=Nm1FCIHam_beta)

        nLinearSystem = nNp1FCIDet + nFCIDet
        nLinearSystem_h = nNm1bFCIDet + nFCIDet

        if(nLinearSystem.ne.nLinearSystem_h) then
            call stop_all(t_r,'Should change restriction on the Linear system size being the same for particle/hole addition')
        endif
        
        write(6,"(A,F14.6,A)") "Memory required for the LR system: ",2*(real(nLinearSystem,dp)**2)*16/1048576.0_dp," Mb"

        !If doing full optimization of the GS problem
        nGSSpace = nFCIDet + nNp1FCIDet + nNm1bFCIDet
        Np1GSInd = nFCIDet + 1
        Nm1GSInd = np1GSInd + nNp1FCIDet
        if(tLR_ReoptGS) then
            write(6,"(A,I9)") "Size of reoptimized ground state basis: ",nGSSpace
            write(6,"(A,F14.6,A)") "Memory required for reoptimization of GS: ",(real(nGSSpace,dp)**2)*16/1048576.0_dp," Mb"
        endif
        
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
        if(tSC_LR) then
            write(iunit,"(A)") "# 1.Frequency     2.GF_LinearResponse(Re)    3.GF_LinearResponse(Im)    " &
                & //"4.ParticleGF(Re)   5.ParticleGF(Im)   6.HoleGF(Re)   7.HoleGF(Im)    8.Old_GS    9.New_GS   " &
                & //"10.Particle_Norm  11.Hole_Norm   12.NI_GF(Re)   13.NI_GF(Im)  14.NI_GF_Part(Re)   15.NI_GF_Part(Im)   " &
                & //"16.NI_GF_Hole(Re)   17.NI_GF_Hole(Im)  18.SpectralWeight  19.Iters_p   20.Iters_h    " &
                & //"21.Tr[SE](Re)  22.Tr[SE](Im) "
        else
            write(iunit,"(A)") "# 1.Frequency     2.GF_LinearResponse(Re)    3.GF_LinearResponse(Im)    " &
                & //"4.ParticleGF(Re)   5.ParticleGF(Im)   6.HoleGF(Re)   7.HoleGF(Im)    8.Old_GS    9.New_GS   " &
                & //"10.Particle_Norm  11.Hole_Norm   12.NI_GF(Re)   13.NI_GF(Im)  14.NI_GF_Part(Re)   15.NI_GF_Part(Im)   " &
                & //"16.NI_GF_Hole(Re)   17.NI_GF_Hole(Im)  18.SpectralWeight  19.Iters_p   20.Iters_h    " 
        endif
                
        if(tMinRes_NonDir) then
            minres_unit = get_free_unit()
            open(minres_unit,file='zMinResQLP.txt',status='unknown',position='append')
            allocate(RHS(nLinearSystem))
            if(tPrecond_MinRes) then
                allocate(Precond_Diag(nLinearSystem))
            endif
        endif
        
        if(tAllImp_LR) then
            !We calculate all greens functions between all impurity sites - both in the non-interacting and interacting case
            nImp_GF = nImp 
        else
            nImp_GF = 1
        endif
        allocate(SchmidtPertGF_Cre_Ket(VirtStart:nSites,nImp_GF))
        allocate(SchmidtPertGF_Cre_Bra(VirtStart:nSites,nImp_GF))
        allocate(SchmidtPertGF_Ann_Ket(1:CoreEnd,nImp_GF))
        allocate(SchmidtPertGF_Ann_Bra(1:CoreEnd,nImp_GF))
        allocate(NI_LRMat_Cre(nImp_GF,nImp_GF))
        allocate(NI_LRMat_Ann(nImp_GF,nImp_GF))
        allocate(ni_lr_Mat(nImp_GF,nImp_GF))
        allocate(ResponseFn_p(nImp_GF,nImp_GF))
        allocate(ResponseFn_h(nImp_GF,nImp_GF))
        allocate(ResponseFn_Mat(nImp_GF,nImp_GF))
        allocate(dNorm_p(nImp_GF))
        allocate(dNorm_h(nImp_GF))
        
        !Allocate memory for hamiltonian in this system:
        allocate(LinearSystem_p(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(LinearSystem_h(nLinearSystem,nLinearSystem),stat=ierr)
        if(tLR_ReoptGS) then
            allocate(GSHam(nGSSpace,nGSSpace),stat=ierr)
            allocate(W(nGSSpace),stat=ierr)
        endif
        allocate(Psi_0(nGSSpace),stat=ierr)   !If tLR_ReoptGS = .F. , then only the first nFCIDet elements will be used
        allocate(Psi1_p(nLinearSystem),stat=ierr)
        allocate(Psi1_h(nLinearSystem),stat=ierr)
        Psi1_p = zzero
        Psi1_h = zzero
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        i = nGSSpace*3 + nLinearSystem*5
        write(6,"(A,F12.5,A)") "Memory required for wavefunctions: ",real(i,dp)*ComptoMb, " Mb"
        
        !Since this does not depend at all on the new basis, since the ground state has a completely
        !disjoint basis to |1>, the RHS of the equations can be precomputed
        allocate(Cre_0(nLinearSystem,nImp_GF))
        allocate(Ann_0(nLinearSystem,nImp_GF))
        if(.not.tLR_ReoptGS) then
            Psi_0(:) = zzero 
            Psi_0(1:nFCIDet) = HL_Vec(:)    
            GFChemPot = HL_Energy
            do i = 1,nImp_GF
                call ApplySP_PertGS_EC(Psi_0,nGSSpace,Cre_0(:,i),Ann_0(:,i),nLinearSystem,i)
            enddo
        endif
        
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

        if(tSC_LR) then
            allocate(SelfEnergy_Imp(nImp,nImp)) !The self-consistently determined self-energy correction to match the interacting and non-interacting greens functions
            SelfEnergy_Imp(:,:) = zzero
            Prev_SE_Saved(:,:) = zzero
            if(iSE_Constraints.eq.1) then
                nVarSE = nImp*nImp     !All elements of self energy are independent
                write(6,"(A,I5)") "Complete flexibility in self-energy fitting. Number of independent variables: ",nVarSE
            else
                nVarSE = (nImp*(nImp+1))/2
                write(6,"(A,I5)") "Self-energy constrained to be off-diagonal hermitian. Number of independent variables: ",nVarSE
            endif
            allocate(SE_Change(nImp,nImp))  !The change in self energy each iteration
            SE_Change(:,:) = zzero
            allocate(Emb_h0v_SE(EmbSize,EmbSize))   !The 1-electron hamiltonian (with self-energy correction) over the embedded basis
            if(allocated(h0v_se)) deallocate(h0v_se)
            allocate(h0v_se(nSites,nSites))     !1 electron hamiltonian with self-energy and correlation potential striped through space
            h0v_se(:,:) = dcmplx(h0v(:,:))              !Initially set it so h0v
            call add_localpot_comp_inplace(h0v_se,SelfEnergy_Imp,.false.)
        endif
        
        !Space for useful intermediates
        !For particle addition
        allocate(Gc_a_F_ax_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Gc_a_F_ax_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Gc_b_F_ab(VirtStart:VirtEnd,nImp_GF),stat=ierr)
        !For hole addition
        allocate(Ga_i_F_xi_Bra(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Ga_i_F_xi_Ket(ActiveStart:ActiveEnd,nImp_GF),stat=ierr)
        allocate(Ga_i_F_ij(1:nCore,nImp_GF),stat=ierr)

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

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            call set_timer(LR_EC_GF_HBuild)
        
            write(6,*) "Calculating linear response for frequency: ",Omega
            if(tMinRes_NonDir) write(minres_unit,*) "Iteratively solving for frequency: ",Omega
!            call flush(6)

            tSCFConverged = .false.
            if(tSC_LR) then
                SE_Fit_Iter = 0
                tSCFFailed = .false.
                if(iReuse_SE.eq.1) then
                    !Start from previous converged self-energy
                    SelfEnergy_Imp(:,:) = Prev_SE_Saved(:,:)
                else
                    SelfEnergy_Imp(:,:) = zzero     !Set self-energy to zero again
                endif
                h0v_se(:,:) = zzero
                do j=1,nSites
                    do i=1,nSites
                        h0v_se(i,j) = dcmplx(h0v(i,j))
                    enddo
                enddo
                call add_localpot_comp_inplace(h0v_se,SelfEnergy_Imp,.false.)   !Subtract self-energy through space
            endif
            do while(.not.tSCFConverged)
                !Self consistently calculate a self-energy function over the impurity sites to match the 
                !mean-field and HL GF calculations.
                if(tSC_LR) then
                    !Find the mean-field solution with an added self-energy. This will also return the
                    !non-interacting greens functions over the requested impurity sites
                    !In addition, this will calculate the 1 electron hamiltonian in the schmidt basis
                    !with the additional self-energy contribution: Emb_h0v_SE
                    call FindNI_Charged(Omega+mu,NI_LRMat_Cre,NI_LRMat_Ann)

                    !If we have a self-consistent self-energy contribution, we also need to update
                    !the hamiltonians which have been stored.
                    !Update tmat with the new one-electron hamiltonian, with self-energy contribution
                    tmat_comp(:,:) = Emb_h0v_SE(:,:)
                    if(tChemPot) tmat_comp(1,1) = tmat_comp(1,1) - dcmplx(U/2.0_dp,0.0_dp)
                    call Fill_N_Np1_Nm1b_FCIHam(Elec,1,1,1,NHam=NFCIHam,Np1Ham=Np1FCIHam_alpha,Nm1Ham=Nm1FCIHam_beta)
                else
                    !First, find the non-interacting solution expressed in the schmidt basis
                    !This will only calculate it over the first impurity site
                    call FindSchmidtPert_Charged(Omega+mu,NI_LRMat_Cre,NI_LRMat_Ann)
                endif

                !Construct useful intermediates, using zgemm
                !sum_a Gc_a^* F_ax (Creation)
                call ZGEMM('T','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_VX(VirtStart:VirtEnd,ActiveStart:ActiveEnd),  &
                    nVirt,SchmidtPertGF_Cre_Bra(VirtStart:VirtEnd,:),nVirt,zzero,Gc_a_F_ax_Bra,EmbSize)
                call ZGEMM('N','N',EmbSize,nImp_GF,nVirt,zone,FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:VirtEnd),  &
                    EmbSize,SchmidtPertGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Gc_a_F_ax_Ket,EmbSize)
                !sum_b Gc_b F_ab  (Creation)
                call ZGEMM('N','N',nVirt,nImp_GF,nVirt,zone,FockSchmidt_SE_VV(VirtStart:VirtEnd,VirtStart:VirtEnd),    &
                    nVirt,SchmidtPertGF_Cre_Ket(VirtStart:VirtEnd,:),nVirt,zzero,Gc_b_F_ab,nVirt)

                !sum_i Ga_i^bra F_xi (Annihilation)
                call ZGEMM('N','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd),EmbSize,    &
                    SchmidtPertGF_Ann_Bra(1:nCore,:),nCore,zzero,Ga_i_F_xi_Bra,EmbSize)
!                call writematrixcomp(FockSchmidt_SE_XC,'FockSchmidt_SE_XC',.false.)
!                call writevectorcomp(SchmidtPertGF_Ann_Bra(:,1),'SchmidtPertGF_Ann_Bra')
!                call writevectorcomp(Ga_i_F_xi_Bra(:,1),'Ga_i_F_xi_Bra')
                call ZGEMM('T','N',EmbSize,nImp_GF,nCore,zone,FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd),nCore,    &
                    SchmidtPertGF_Ann_Ket(1:nCore,:),nCore,zzero,Ga_i_F_xi_Ket,EmbSize)
!                call writevectorcomp(Ga_i_F_xi_Bra(:,1),'Ga_i_F_xi_Ket')
                !sum_i Ga_i F_ij (Annihilation)
                call ZGEMM('T','N',nCore,nImp_GF,nCore,zone,FockSchmidt_SE_CC(:,:),nCore,SchmidtPertGF_Ann_Ket(:,:),nCore,zzero, &
                    Ga_i_F_ij,nCore)

                !Find the first-order wavefunction in the interacting picture for all impurity sites 
                do pertsite = 1,nImp_GF

                    if(tLR_ReoptGS) GSHam(:,:) = zzero 
                    LinearSystem_p(:,:) = zzero 
                    LinearSystem_h(:,:) = zzero 

                    !Block 1 for particle hamiltonian
                    !First, construct n + 1 (alpha) FCI space, in determinant basis
                    LinearSystem_p(1:nNp1FCIDet,1:nNp1FCIDet) = Np1FCIHam_alpha(:,:)
                    !Block 1 for hole hamiltonian
                    LinearSystem_h(1:nNm1bFCIDet,1:nNm1bFCIDet) = Nm1FCIHam_beta(:,:)
                    if(tLR_ReoptGS) then
                        !Block 1 for GS
                        GSHam(1:nFCIDet,1:nFCIDet) = nFCIHam(:,:)
                        !Block 2 (diagonal part to come)
                        GSHam(Np1GSInd:Nm1GSInd-1,Np1GSInd:Nm1GSInd-1) = Np1FCIHam_alpha(:,:)
                        !Block 4 (diagonal part to come)
                        GSHam(Nm1GSInd:nGSSpace,Nm1GSInd:nGSSpace) = Nm1FCIHam_beta(:,:)
                    endif

                    !Block 2 for particle hamiltonian
                    !Copy the N electron FCI hamiltonian to this diagonal block
                    LinearSystem_p(VIndex:nLinearSystem,VIndex:nLinearSystem) = nFCIHam(:,:)
            
                    VNorm = zzero
                    tempel = zzero 
                    do a = VirtStart,VirtEnd
                        !Calc normalization for the CV block
                        VNorm = VNorm + SchmidtPertGF_Cre_Bra(a,pertsite)*SchmidtPertGF_Cre_Ket(a,pertsite) 
                        !VNorm = VNorm + real(SchmidtPertGF_Cre(a,pertsite))**2 + aimag(SchmidtPertGF_Cre(a,pertsite))**2
                        !Calculate the diagonal correction
                        !tempel = tempel + conjg(SchmidtPertGF_Cre(a,pertsite))*Gc_b_F_ab(a,pertsite)
                        tempel = tempel + SchmidtPertGF_Cre_Bra(a,pertsite)*Gc_b_F_ab(a,pertsite)
                    enddo
                    tempel = tempel / VNorm
                    !Add diagonal virtual correction
                    do i = VIndex,nLinearSystem
                        LinearSystem_p(i,i) = LinearSystem_p(i,i) + tempel
                    enddo
                    if(tLR_ReoptGS) then
                        !Diagonal part of block 4
                        do i = Nm1GSInd,nGSSpace
                            GSHam(i,i) = GSHam(i,i) + tempel
                        enddo
                    endif

                    !Block 2 for the hole hamiltonian
                    !Copy the N electron FCI hamiltonian to this diagonal blockk
                    LinearSystem_h(VIndex:nLinearSystem,VIndex:nLinearSystem) = nFCIHam(:,:)

                    CNorm = zzero
                    tempel = zzero 
                    do i = 1,CoreEnd
                        !Calc normalization
                        !CNorm = CNorm + real(SchmidtPertGF_Ann(i,pertsite))**2 + aimag(SchmidtPertGF_Ann(i,pertsite))**2
                        CNorm = CNorm + SchmidtPertGF_Ann_Bra(i,pertsite)*SchmidtPertGF_Ann_Ket(i,pertsite)
                        !Calc diagonal correction
                        tempel = tempel + SchmidtPertGF_Ann_Bra(i,pertsite)*Ga_i_F_ij(i,pertsite)
                    enddo
                    tempel = -tempel / CNorm
                    !Add diagonal correction
                    do i = VIndex,nLinearSystem
                        LinearSystem_h(i,i) = LinearSystem_h(i,i) + tempel
                    enddo
                    if(tLR_ReoptGS) then
                        !Diagonal part of block 2
                        do i = Np1GSInd,Nm1GSInd-1
                            GSHam(i,i) = GSHam(i,i) + tempel
                        enddo
                    endif
                
                    !Block 3 for particle hamiltonian
                    do J = 1,nNp1FCIDet
                        do K = 1,nFCIDet    !VIndex,nLinearSystem
                            if(Coup_Ann_alpha(1,K,J).ne.0) then
                                !Determinants are connected via a single SQ operator
                                !First index is the spatial index, second is parity 
                                !Matrices not necessarily hermitian with self-energy term, so need to construct them seperately
                                LinearSystem_p(K+VIndex-1,J) = Gc_a_F_ax_Bra(Coup_Ann_alpha(1,K,J),pertsite)    &
                                    *Coup_Ann_alpha(2,K,J)/sqrt(VNorm)
                                !LinearSystem_p(J,K+VIndex-1) = conjg(LinearSystem_p(K+VIndex-1,J))
                                LinearSystem_p(J,K+VIndex-1) = Gc_a_F_ax_Ket(Coup_Ann_alpha(1,K,J),pertsite)    &
                                    *Coup_Ann_alpha(2,K,J)/sqrt(VNorm)
                                if(tLR_ReoptGS) then
                                    !Block 3 contribution
                                    !Assume that all the matrices are hermitian as we shouldn't be reoptimizing the GS with self-energy
                                    GSHam(K,Np1GSInd+J-1) = -Ga_i_F_xi_Ket(Coup_Ann_alpha(1,K,J),pertsite)  &
                                        *Coup_Ann_alpha(2,K,J)/sqrt(CNorm)
                                    GSHam(Np1GSInd+J-1,K) = conjg(GSHam(K,Np1GSInd+J-1))
                                endif
                            endif
                        enddo
                    enddo

                    !Block 3 for hole hamiltonian
                    do J = 1,nNm1bFCIDet
                        do K = 1,nFCIDet
                            if(Coup_Create_alpha(1,K,J).ne.0) then
                                !Determinants are connected via a single SQ operator
                                LinearSystem_h(K+VIndex-1,J) = - Ga_i_F_xi_Bra(Coup_Create_alpha(1,K,J),pertsite)   &
                                    *Coup_Create_alpha(2,K,J)/sqrt(CNorm)
                                LinearSystem_h(J,K+VIndex-1) = - Ga_i_F_xi_Ket(Coup_Create_alpha(1,K,J),pertsite)   &
                                    *Coup_Create_alpha(2,K,J)/sqrt(CNorm)
!                                LinearSystem_h(J,K+VIndex-1) = conjg(LinearSystem_h(K+VIndex-1,J))
                                if(tLR_ReoptGS) then
                                    !Block 6 contribution
                                    GSHam(K,Nm1GSInd+J-1) = Gc_a_F_ax_Ket(Coup_Create_alpha(1,K,J),pertsite) &
                                        *Coup_Create_alpha(2,K,J)/sqrt(VNorm)
                                    GSHam(Nm1GSInd+J-1,K) = conjg(GSHam(K,Nm1GSInd+J-1))
                                endif
                            endif
                        enddo
                    enddo
                    
                    !call writematrixcomp(LinearSystem_p,'LinearSystem_p',.false.)
                    !call writematrixcomp(LinearSystem_h,'LinearSystem_h',.false.)

                    if(.not.tSC_LR) then
                        !Now, check hessian is hermitian
                        do i = 1,nLinearSystem
                            do j=i,nLinearSystem
                                if(abs(LinearSystem_p(i,j)-conjg(LinearSystem_p(j,i))).gt.1.0e-8_dp) then
                                    write(6,*) "i, j: ",i,j
                                    write(6,*) "LinearSystem_p(i,j): ",LinearSystem_p(i,j)
                                    write(6,*) "LinearSystem_p(j,i): ",LinearSystem_p(j,i)
                                    call stop_all(t_r,'Particle hessian for EC-LR not hermitian')
                                endif
                                if(abs(LinearSystem_h(i,j)-conjg(LinearSystem_h(j,i))).gt.1.0e-8_dp) then
                                    write(6,*) "i, j: ",i,j
                                    write(6,*) "LinearSystem_h(i,j): ",LinearSystem_h(i,j)
                                    write(6,*) "LinearSystem_h(j,i): ",LinearSystem_h(j,i)
                                    call stop_all(t_r,'Hole hessian for EC-LR not hermitian')
                                endif
                            enddo
                        enddo

                        if(tLR_ReoptGS) then
                            !Check GS hamiltonian is hermitian
                            do i = 1,nGSSpace
                                do j = i,nGSSpace
                                    if(abs(GSHam(j,i)-conjg(GSHam(i,j))).gt.1.0e-8_dp) then
                                        call stop_all(t_r,'Reoptimized ground state hamiltonian is not hermitian')
                                    endif
                                enddo
                            enddo
                        endif
                    endif

                    !write(6,*) "Hessian constructed successfully...",Omega
                    call halt_timer(LR_EC_GF_HBuild)
                    call set_timer(LR_EC_GF_OptGS)

                    if(tLR_ReoptGS) then
                        if(tNonDirDavidson) then
                            call Comp_NonDir_Davidson(nGSSpace,GSHam,GFChemPot,Psi_0,.false.)
                        else
                            allocate(Work(max(1,3*nGSSpace-2)))
                            allocate(temp_vecc(1))
                            W(:) = 0.0_dp
                            lWork = -1
                            info = 0
                            call zheev('V','U',nGSSpace,GSHam,nGSSpace,W,temp_vecc,lWork,Work,info)
                            if(info.ne.0) call stop_all(t_r,'workspace query failed')
                            lwork = int(temp_vecc(1))+1
                            deallocate(temp_vecc)
                            allocate(temp_vecc(lwork))
                            call zheev('V','U',nGSSpace,GSHam,nGSSpace,W,temp_vecc,lWork,Work,info)
                            if(info.ne.0) call stop_all(t_r,'GS H Diag failed 1')
                            deallocate(work,temp_vecc)

                            Psi_0(:) = GSHam(:,1)
                            GFChemPot = W(1)
                        endif

                        write(6,"(A,G20.10,A,G20.10,A)") "Reoptimized ground state energy is: ",GFChemPot, &  
                            " (old = ",HL_Energy,")"
                        !call writevector(HL_Vec(:),'Old Psi_0')
                        !call writevectorcomp(Psi_0,'New Psi_0')
                    endif
                    !write(6,*) "E0: ",GFChemPot+CoreEnergy

                    call halt_timer(LR_EC_GF_OptGS)
                    call set_timer(LR_EC_GF_SolveLR)

                    !Solve particle GF to start
                    !The V|0> for particle and hole perturbations are held in Cre_0 and Ann_0
                    !If we have reoptimized the ground state, we will need to recompute these
                    !call writevectorcomp(Psi_0,'Psi_0')
                    if(tLR_ReoptGS) then
                        !Since we are only keeping these one at a time, we don't actually need to store all pertsite vectors
                        call ApplySP_PertGS_EC(Psi_0,nGSSpace,Cre_0(:,pertsite),Ann_0(:,pertsite),nLinearSystem,pertsite)
                    endif
                    !call writevectorcomp(Psi_0,'Psi_0')
                    !call writevectorcomp(Cre_0(:,1),'Ann_0')
                    !call writevectorcomp(Ann_0,'Ann_0')
                    
                    LinearSystem_p(:,:) = -LinearSystem_p(:,:)
                    !Offset matrix
                    do i = 1,nLinearSystem
                        LinearSystem_p(i,i) = LinearSystem_p(i,i) + dcmplx(Omega+mu+GFChemPot,dDelta)
                    enddo

                    !Now solve these linear equations
                    !call writevectorcomp(Psi1_p,'Cre_0')
                    if(tMinRes_NonDir) then
!                        zShift = dcmplx(-Omega-mu-GFChemPot,-dDelta)
                        zDirMV_Mat => LinearSystem_p
!                        call writematrixcomp(LinearSystem_p,'LinearSystem_p',.true.)
                        call setup_RHS(nLinearSystem,Cre_0(:,pertsite),RHS)
                        maxminres_iter_ip = int(maxminres_iter,ip)
                        minres_unit_ip = int(minres_unit,ip)
                        nLinearSystem_ip = int(nLinearSystem,ip)
                        if(tPrecond_MinRes) then
                            call FormPrecond(nLinearSystem)
                            call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_p, &
                                itnlim=maxminres_iter_ip,Msolve=zPreCond,istop=info_ip,rtol=rtol_LR,itn=iters_p, &
                                startguess=tReuse_LS)
                        else
                            call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip, &
                                itnlim=maxminres_iter_ip,rtol=rtol_LR,x=Psi1_p,istop=info_ip,itn=iters_p,startguess=tReuse_LS)
                        endif
                        info = info_ip
                        zDirMV_Mat => null()
                        if(info.gt.7) write(6,*) "info: ",info
                        if(info.eq.8) write(6,"(A,I9)") "Linear equation solver hit iteration limit: ",iters_p
                        if((info.eq.9).or.(info.eq.10).or.(info.eq.11)) then
                            call stop_all(t_r,'Input matrices to linear solver incorrect')
                        endif
                        if(info.gt.11) call stop_all(t_r,'Linear equation solver failed')
                    else
                        !Copy V|0> to another array, since if we are not reoptimizing the GS, we want to keep them.
                        Psi1_p(:) = Cre_0(:,pertsite)
                        !Psi1_p will be overwritten with the solution
                        call SolveCompLinearSystem(LinearSystem_p,Psi1_p,nLinearSystem,info)
                        if(info.ne.0) then 
                            write(6,*) "INFO: ",info
                            call warning(t_r,'Solving linear system failed for particle hamiltonian - skipping this frequency')
                            Omega = Omega + Omega_Step
                            call halt_timer(LR_EC_GF_SolveLR)
                            cycle
                        endif
                    endif
                    !call writevectorcomp(Psi1_p,'Psi1_p')

                    !Now solve the LR for the hole addition
                    do i = 1,nLinearSystem
                        LinearSystem_h(i,i) = dcmplx(Omega+mu,dDelta) + (LinearSystem_h(i,i) - dcmplx(GFChemPot,0.0_dp))
                    enddo
                    if(tMinRes_NonDir) then
                        !zShift = dcmplx(-Omega-mu+GFChemPot,-dDelta)
                        zDirMV_Mat => LinearSystem_h
!                        call writematrixcomp(LinearSystem_h,'LinearSystem_h',.true.)
                        call setup_RHS(nLinearSystem,Ann_0(:,pertsite),RHS)
                        maxminres_iter_ip = int(maxminres_iter,ip)
                        minres_unit_ip = int(minres_unit,ip)
                        nLinearSystem_ip = int(nLinearSystem,ip)
                        if(tPrecond_MinRes) then
                            call FormPrecond(nLinearSystem)
                            call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_h, &
                                itnlim=maxminres_iter_ip,Msolve=zPreCond,istop=info_ip,rtol=rtol_LR,itn=iters_h, &
                                startguess=tReuse_LS)
                        else
                            call MinResQLP(n=nLinearSystem_ip,Aprod=zDirMV,b=RHS,nout=minres_unit_ip,x=Psi1_h, &
                                itnlim=maxminres_iter_ip,istop=info_ip,rtol=rtol_LR,itn=iters_h,startguess=tReuse_LS)
                        endif
                        info = info_ip
                        zDirMV_Mat => null()
                        if(info.gt.7) write(6,*) "info: ",info
                        if(info.eq.8) write(6,"(A,I9)") "Linear equation solver hit iteration limit: ",iters_h
                        if((info.eq.9).or.(info.eq.10).or.(info.eq.11)) then
                            call stop_all(t_r,'Input matrices to linear solver incorrect')
                        endif
                        if(info.gt.11) call stop_all(t_r,'Linear equation solver failed')
                    else
                        Psi1_h(:) = Ann_0(:,pertsite)
                        call SolveCompLinearSystem(LinearSystem_h,Psi1_h,nLinearSystem,info)
                        if(info.ne.0) then 
                            write(6,*) "INFO: ",info
                            call warning(t_r,'Solving linear system failed for hole hamiltonian - skipping this frequency')
                            Omega = Omega + Omega_Step
                            call halt_timer(LR_EC_GF_SolveLR)
                            cycle
                        endif
                    endif
            
                    !Find normalization of first-order wavefunctions
                    dNorm_p(pertsite) = znrm2(nLinearSystem,Psi1_p,1)
                    dNorm_h(pertsite) = znrm2(nLinearSystem,Psi1_h,1)
                    !write(6,*) "Normalization of first-order wavefunction: ",dNorm_p,dNorm_h
                    
                    !Now, calculate all interacting greens funtions that we want, by dotting in the first-order space
                    if((nImp_GF.gt.1).and.(tLR_ReoptGS)) then
                        call stop_all(t_r,'Cannot do dotting of vectors here, since the LHS vectors are still being computed')
                    endif
                    do j = 1,nImp_GF
                        ResponseFn_p(pertsite,j) = zdotc(nLinearSystem,Cre_0(:,j),1,Psi1_p,1)
                        ResponseFn_h(pertsite,j) = zdotc(nLinearSystem,Ann_0(:,j),1,Psi1_h,1)
                        ResponseFn_Mat(pertsite,j) = ResponseFn_p(pertsite,j) + ResponseFn_h(pertsite,j)
                        ni_lr_Mat(pertsite,j) = NI_LRMat_Cre(pertsite,j) + NI_LRMat_Ann(pertsite,j) 
                    enddo

                enddo   !End do over pertsite. We now have all the greens funtions 

                !Now, calculate NI and interacting response function as trace over diagonal parts of the local greens functions
                ni_lr = zzero
                ResponseFn = zzero
                AvResFn_p = zzero
                AvResFn_h = zzero
                ni_lr_p = zzero
                ni_lr_h = zzero
                AvdNorm_p = zero
                AvdNorm_h = zero
                do i = 1,nImp_GF
                    ni_lr = ni_lr + ni_lr_Mat(i,i)
                    ResponseFn = ResponseFn + ResponseFn_Mat(i,i)
                    AvResFn_p = AvResFn_p + ResponseFn_p(i,i)
                    AvResFn_h = AvResFn_h + ResponseFn_h(i,i)
                    ni_lr_p = ni_lr_p + NI_LRMat_Cre(i,i)
                    ni_lr_h = ni_lr_h + NI_LRMat_Ann(i,i)
                    AvdNorm_p = AvdNorm_p + dNorm_p(i)
                    AvdNorm_h = AvdNorm_h + dNorm_h(i)
                enddo
                ni_lr = ni_lr/real(nImp_GF,dp)
                ni_lr_p = ni_lr_p/real(nImp_GF,dp)
                ni_lr_h = ni_lr_h/real(nImp_GF,dp)
                ResponseFn = ResponseFn/real(nImp_GF,dp)
                AvResFn_p = AvResFn_p/real(nImp_GF,dp)
                AvResFn_h = AvResFn_h/real(nImp_GF,dp)
                AvdNorm_p = AvdNorm_p/real(nImp_GF,dp)
                AvdNorm_h = AvdNorm_h/real(nImp_GF,dp)

                Diff_GF = zero  !real(abs(ResponseFn - ni_lr))
                do i=1,nImp_GF
                    do j=1,nImp_GF
                        Diff_GF = Diff_GF + real((ni_lr_Mat(j,i)-ResponseFn_Mat(j,i))*dconjg(ni_lr_Mat(j,i)-ResponseFn_Mat(j,i)))
                    enddo
                enddo
            
                call halt_timer(LR_EC_GF_SolveLR)
                call set_timer(LR_EC_GF_FitGF)
                
                if(tSC_LR) then
                    !Do the fitting of the self energy and iterate.
                    !Update SelfEnergy_Imp and h0v_se
                    SE_Fit_Iter = SE_Fit_Iter + 1
                    
                    if(SE_Fit_Iter.eq.1) then
                        write(6,"(A)") "    Iter  NR_Iters   Delta_SE         Diff_GFs            Start_Diff_GFs      Orig_NI_GF"&
     &                      //"                               HL_GF"
                    endif
                    if((tPartialSE_Fit.and.(SE_Fit_Iter.le.iPartialSE_Fit)).or.(.not.tPartialSE_Fit)) then
                        if((iReuse_SE.eq.2).and.(SE_Fit_Iter.eq.1)) then
                            !Before first fitting, send in the previous converged self-energy to start from
                            SE_Change(:,:) = Prev_SE_Saved(:,:) 
                        else
                            SE_Change(:,:) = zzero
                        endif
                        call Fit_SE(SE_Change,Var_SE,Error_GF,nNR_Iters,ResponseFn_Mat,Omega+mu,SE_Fit_Iter,ni_lr_Mat)
                        !Write out
                        write(6,"(2I7,7G20.10)") SE_Fit_Iter,nNR_Iters,Var_SE,Error_GF,Diff_GF,ni_lr,ResponseFn
                        !call writematrixcomp(SE_Change,'SE_Change',.true.)
!                        if(Var_SE.gt.1.0e4_dp) then
!                            !The self-energy has diverged
!                            write(6,*) "DIVERGENT SELF-ENERGY TERM. Skipping frequency."
!                            tSCFConverged = .true.
!                            SE_Change(:,:) = zzero
!                            tSCFFailed = .true.
!                        elseif(SE_Fit_Iter.gt.750) then
!                            write(6,*) "SELF-ENERGY TERM NOT CONVERGED. Skipping frequency."
!                            tSCFConverged = .true.
!                            tSCFFailed = .true.
!                        endif

!                        if(iGF_Fit.eq.4) then
                            !If directly calculating self energy from the inverses, only add it to the impurity block
                            !Do not update self energy as we don't want it striped through the space, or added to the high-level calculation
                            !h0v_SE(1:nImp,1:nImp) = h0v_SE(1:nImp,1:nImp) - SE_Change(:,:) 
!                            write(6,*) "SE Change: ",SE_Change
!                            write(6,*) "Self energy: ",SelfEnergy_Imp
!                        else
!                            !Update self-energy
!                            !Emb_h0v_SE and all fock matrices are updated in FindNI_Charged routine 
!                            call add_localpot_comp_inplace(h0v_se,SE_Change,.false.)
!                            SelfEnergy_Imp(:,:) = SelfEnergy_Imp(:,:) + SE_Change(:,:)
!                        endif
                    else
                        !Write out just to show the difference in the GFs
                        write(6,"(2I7,7G20.10)") SE_Fit_Iter,0,zero,Diff_GF,Diff_GF,ni_lr,ResponseFn
                        tSCFConverged = .true.
                    endif

                    if(Var_SE.lt.1.0e-8_dp) then
                        !if(Error_GF.gt.1.0e-4) then
                        !    call stop_all(t_r,"Greens functions not consistent, but self-energy no longer changing due to damping")
                        !endif
                        !Yay - converged
                        tSCFConverged = .true.
                    endif
                else
                    tSCFConverged = .true.  !We only do one cycle, and do not try to match the greens functions.
                endif

                call halt_timer(LR_EC_GF_FitGF)

            enddo   !Finish the self consistency

            if(.not.tFirst) then
                SpectralWeight = SpectralWeight + Omega_Step*(Prev_Spec-aimag(ResponseFn))/(2.0_dp*pi)
                Prev_Spec = -aimag(ResponseFn)
            endif

            if(tSC_LR) then
                !Find mean self-energy per site
                TR_SE = zzero
                do i=1,nImp_GF
                    TR_SE = TR_SE + SelfEnergy_Imp(i,i)
                enddo
                TR_SE = TR_SE/dcmplx(nImp_GF,0.0_dp)

                write(iunit,"(18G22.10,2I7,2G22.10)") Omega,real(ResponseFn),-aimag(ResponseFn), &
                    real(AvResFn_p),-aimag(AvResFn_p),real(AvResFn_h),-aimag(AvResFn_h),    &
                    HL_Energy,GFChemPot,abs(AvdNorm_p),abs(AvdNorm_h),real(ni_lr),-aimag(ni_lr),real(ni_lr_p),    &
                    -aimag(ni_lr_p),real(ni_lr_h),-aimag(ni_lr_h),SpectralWeight,iters_p,iters_h,   &
                    real(TR_SE),aimag(TR_SE)

                call writematrixcomp(SelfEnergy_Imp,'Self Energy',.true.)
                if(.not.tSCFFailed) then
                    Prev_SE_Saved(:,:) = SelfEnergy_Imp(:,:)    !Save the self-energy contribution from the previous iteration
                else
                    Prev_SE_Saved(:,:) = zzero
                endif
            else
                write(iunit,"(18G22.10,2I7)") Omega,real(ResponseFn),-aimag(ResponseFn), &
                    real(AvResFn_p),-aimag(AvResFn_p),real(AvResFn_h),-aimag(AvResFn_h),    &
                    HL_Energy,GFChemPot,abs(AvdNorm_p),abs(AvdNorm_h),real(ni_lr),-aimag(ni_lr),real(ni_lr_p),    &
                    -aimag(ni_lr_p),real(ni_lr_h),-aimag(ni_lr_h),SpectralWeight,iters_p,iters_h
            endif
            call flush(iunit)
            call flush(6)

            if(tFirst) tFirst = .false.

            Omega = Omega + Omega_Step
        enddo   !End loop over omega

        write(6,"(A,G22.10)") "Total integrated spectral weight: ",SpectralWeight

        if(tLR_ReoptGS) then
            deallocate(W,GSHam)
        endif
        deallocate(LinearSystem_p,LinearSystem_h,Psi1_p,Psi1_h)
        deallocate(Cre_0,Ann_0,Psi_0,SchmidtPertGF_Cre_Ket,SchmidtPertGF_Ann_Ket)
        deallocate(NI_LRMat_Cre,NI_LRMat_Ann,ResponseFn_p,ResponseFn_h,ResponseFn_Mat)
        deallocate(ni_lr_Mat,SchmidtPertGF_Cre_Bra,SchmidtPertGF_Ann_Bra)
        if(tSC_LR) then
            deallocate(SelfEnergy_Imp,Emb_h0v_SE,SE_Change)
        endif
        close(iunit)
        if(tMinRes_NonDir) then
            close(minres_unit)
            deallocate(RHS)
            if(tPrecond_MinRes) deallocate(Precond_Diag)
        endif

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

        !Stored intermediates
        deallocate(NFCIHam,Nm1FCIHam_beta,Np1FCIHam_alpha)
        deallocate(Coup_Create_alpha,Coup_Ann_alpha)
        deallocate(FockSchmidt_SE,FockSchmidt_SE_VV,FockSchmidt_SE_CC)
        deallocate(FockSchmidt_SE_VX,FockSchmidt_SE_CX,FockSchmidt_SE_XV,FockSchmidt_SE_XC)
        deallocate(Gc_a_F_ax_Bra,Gc_a_F_ax_Ket,Gc_b_F_ab,Ga_i_F_xi_Bra,Ga_i_F_xi_Ket,Ga_i_F_ij)

    end subroutine NonIntExCont_TDA_MCLR_Charged

    subroutine NonIntExContracted_TDA_MCLR()
        use utils, only: get_free_unit,append_ext_real,append_ext
        use DetBitOps, only: DecodeBitDet,SQOperator,CountBits
        use zminresqlpModule, only: MinresQLP  
        use DetToolsData
        use DetTools, only: tospat,umatind,GenDets,GetHElement
        use solvers, only: CreateIntMats
        implicit none
        integer :: a,i,j,k,AVInd_tmp,b,beta
        integer :: CoreEnd,VirtStart,VirtEnd,CVInd_tmp,iunit,iunit2
        integer :: DiffOrb,nOrbs,gam,gam1,gam1_ind,gam1_spat,gam2,gam2_ind,gam2_spat,ierr
        integer :: gam_spat,nLinearSystem,tempK,nSpan,ActiveEnd,ActiveStart
        integer :: orbdum(1),CAInd_tmp,lwork,info,nSize,nCore,nVirt,nFullNm1,nFullNp1,iters
        integer :: maxminres_iter,minres_unit,InterMem,HamMem,CoupMem
        logical :: tParity,tCalcResponse,tTransformSpace
        real(dp) :: Omega,CVNorm,GSEnergy
        complex(dp) :: tempel,ResponseFn,dOrthog,dNorm,testc,ni_lr
        complex(dp) , allocatable :: LinearSystem(:,:),Overlap(:,:),Projector(:,:),VGS(:),CanTrans(:,:)
        complex(dp) , allocatable :: temp_vecc(:),tempc(:,:),RHS(:),Transform(:,:),OrthogHam(:,:)
        complex(dp) , allocatable :: Psi_0(:),S_EigVec(:,:),G_ai_G_aj(:,:),G_ai_G_bi(:,:),G_xa_G_ya(:,:)
        complex(dp) , allocatable :: G_xa_F_ab(:,:),G_xa_G_yb_F_ab(:,:),FockSchmidtComp(:,:),G_ia_G_xa(:,:)
        complex(dp) , allocatable :: F_xi_G_ia_G_ya(:,:),G_xa_F_ya(:,:),G_xi_G_yi(:,:),G_ix_F_ij(:,:)
        complex(dp) , allocatable :: G_ix_G_jy_F_ji(:,:),G_ia_G_ix(:,:),F_ax_G_ia_G_iy(:,:),G_ix_F_iy(:,:)
        complex(dp) , allocatable :: SBlock(:,:),temp_vecc_2(:),S_Diag(:),temp_vecc_3(:)
        complex(dp) , allocatable , target :: LHS(:,:)
        real(dp), allocatable :: NFCIHam(:,:),Nm1FCIHam_alpha(:,:),Nm1FCIHam_beta(:,:)
        real(dp), allocatable :: Np1FCIHam_alpha(:,:),Np1FCIHam_beta(:,:),SBlock_val(:)
        real(dp), allocatable :: AVNorm(:),CANorm(:),Work(:),H_Vals(:),S_EigVal(:)
        integer, allocatable :: Coup_Ann_alpha(:,:,:),Coup_Ann_beta(:,:,:)
        integer, allocatable :: Coup_Create_alpha(:,:,:),Coup_Create_beta(:,:,:)
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='NonIntExContracted_TDA_MCLR'
        integer(ip) :: maxminres_iter_ip,minres_unit_ip,nSpan_ip,info_ip
        complex(dp) , allocatable :: RHS2(:),Psi1(:)
        logical, parameter :: tDiagHam = .false.

        call set_timer(LR_EC_TDA_Precom)
        
        maxminres_iter = iMinRes_MaxIter
        iters = 0

        write(6,*) "Calculating non-interacting EC MR-TDA LR system..."
        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
        if(tMinRes_NonDir) then
            if(tPreCond_MinRes) then
                write(6,"(A)") "Solving linear system with iterative non-direct preconditioned MinRes-QLP algorithm"
            else
                write(6,"(A)") "Solving linear system with iterative non-direct MinRes-QLP algorithm"
            endif
            write(6,"(A,G22.10)") "Tolerance for solution of linear system: ",rtol_LR
            write(6,"(A,G22.10)") "Maximum iterations for each solution: ",maxminres_iter
        else
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
        endif
        if(tRemoveGSFromH) then
            write(6,"(A)") "Explicitly removing the ground state from the hamiltonian. Hamiltonian will now lose hermiticity"
        else
            write(6,"(A)") "Hamiltonian will not have ground state explicitly removed. Should remain hermitian"
        endif

        !umat and tmat for the active space
        call CreateIntMats(tComp=.false.)

        if(tMinRes_NonDir) then
            minres_unit = get_free_unit()
            open(minres_unit,file='zMinResQLP.txt',status='unknown',position='append')
        endif
        
        !Enumerate excitations for fully coupled space
        !Seperate the lists into different Ms sectors in the N+- lists
        call GenDets(Elec,EmbSize,.true.,.true.,.true.) 
        write(6,*) "Number of determinants in {N,N+1,N-1} FCI space: ",ECoupledSpace

        !if(allocated(Nm1BitList)) call writevectorint(Nm1BitList,'Nm1BitList')
        !if(allocated(Nm1bBitList)) call writevectorint(Nm1bBitList,'Nm1bBitList')
        !if(allocated(Np1BitList)) call writevectorint(Np1BitList,'Np1BitList')
        !if(allocated(Np1bBitList)) call writevectorint(Np1bBitList,'Np1bBitList')

        !Construct FCI hamiltonians for the N, N-1_alpha, N-1_beta, N+1_alpha and N+1_beta spaces
        !N electron
        allocate(NFCIHam(nFCIDet,nFCIDet))
        NFCIHam(:,:) = 0.0_dp
        do i = 1,nFCIDet
            do j = 1,nFCIDet
                call GetHElement(FCIDetList(:,i),FCIDetList(:,j),Elec,NFCIHam(i,j))
            enddo
        enddo
        !N-1 hamiltonian
        if(nNm1FCIDet.ne.nNm1bFCIDet) call stop_all(t_r,'Cannot deal with open shell systems')
        allocate(Nm1FCIHam_alpha(nNm1FCIDet,nNm1FCIDet))
        Nm1FCIHam_alpha(:,:) = 0.0_dp
        allocate(Nm1FCIHam_beta(nNm1FCIDet,nNm1FCIDet))
        Nm1FCIHam_beta(:,:) = 0.0_dp
        do i=1,nNm1FCIDet
            do j=1,nNm1FCIDet
                call GetHElement(Nm1FCIDetList(:,i),Nm1FCIDetList(:,j),Elec-1,Nm1FCIHam_alpha(i,j))
                call GetHElement(Nm1bFCIDetList(:,i),Nm1bFCIDetList(:,j),Elec-1,Nm1FCIHam_beta(i,j))
            enddo
        enddo

        !N+1 hamiltonian
        allocate(Np1FCIHam_alpha(nNp1FCIDet,nNp1FCIDet))
        Np1FCIHam_alpha(:,:) = 0.0_dp
        allocate(Np1FCIHam_beta(nNp1FCIDet,nNp1FCIDet))
        Np1FCIHam_beta(:,:) = 0.0_dp
        do i=1,nNp1FCIDet
            do j=1,nNp1FCIDet
                call GetHElement(Np1FCIDetList(:,i),Np1FCIDetList(:,j),Elec+1,Np1FCIHam_alpha(i,j))
                call GetHElement(Np1bFCIDetList(:,i),Np1bFCIDetList(:,j),Elec+1,Np1FCIHam_beta(i,j))
            enddo
        enddo
        HamMem = nNp1FCIDet**2 + nNm1FCIDet**2 + nFCIDet**2
        write(6,'(A,F9.3,A)') "Memory required to store hamiltonians: ",real(HamMem,dp)*RealtoMb,' Mb'

        nLinearSystem = (2*nFCIDet) + nImp*2*(nNm1FCIDet+nNm1bFCIDet) + nImp*2*(nNp1FCIDet+nNp1bFCIDet) 
        
        iunit = get_free_unit()
        call append_ext_real('EC-TDA_DDResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# 1.Frequency     2.DD_LinearResponse(Re)    3.DD_LinearResponse(Im)    " &
            & //"4.Orthog    5.Norm    6.NewGS   7.OldGS  8.OrthogRHS    9.NI_LR(Re)   10.NI_LR(Im)   11.Iters"
        
        if(tDiagHam) then
            iunit2 = get_free_unit()
            call append_ext_real('EC-TDA_EValues',U,filename)
            if(.not.tHalfFill) then
                !Also append occupation of lattice to the filename
                call append_ext(filename,nOcc,filename2)
            else
                filename2 = filename
            endif
            open(unit=iunit2,file=filename2,status='unknown')
            write(iunit2,"(A)") "# Frequency     EValues..."
        endif

        !Allocate memory for hmailtonian in this system:
        write(6,'(A,F13.3,A)') "Memory required to store overlap & linear system: ",    &
            2.0_dp*(real(nLinearSystem,dp)**2)*ComptoMb,' Mb'
        allocate(LinearSystem(nLinearSystem,nLinearSystem),stat=ierr)
        allocate(Overlap(nLinearSystem,nLinearSystem),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
        
        !Set up orbital indices
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
        !And full N pm 1 space sizes
        nFullNm1 = nNm1FCIDet + nNm1bFCIDet
        nFullNp1 = nNp1FCIDet + nNp1bFCIDet

        !Set up indices for the block of the linear system
        CVIndex = nFCIDet + 1   !Beginning of EC Core virtual excitations
        AVIndex = nFCIDet + nFCIDet + 1 !Beginning of EC Active Virtual excitations
        CAIndex = AVIndex + (nImp*2)*nFullNm1 !Beginning of EC Core Active excitations

        write(6,*) "CV indices start from: ",CVIndex
        write(6,*) "AV indices start from: ",AVIndex
        write(6,*) "CA indices start from: ",CAIndex
        write(6,*) "Total size of linear sys: ",nLinearSystem
            
        !Allocate memory for normalization constants
        allocate(AVNorm(1:nImp*2))
        allocate(CANorm(1:nImp*2))

        !Allocate and precompute 1-operator coupling coefficients between the different sized spaces.
        !First number is the index of operator, and the second is the parity change when applying the operator
        CoupMem = 2*nFCIDet*nNm1bFCIDet + 2*nFCIDet*nNm1FCIDet + 2*nFCIDet*nNp1FCIDet + 2*nFCIDet*nNp1bFCIDet
        write(6,'(A,F9.3,A)') "Memory required to store coupling matrices: ",real(CoupMem,dp)*RealtoMb,' Mb'

        allocate(Coup_Create_alpha(2,nFCIDet,nNm1bFCIDet))
        allocate(Coup_Create_beta(2,nFCIDet,nNm1FCIDet))
        allocate(Coup_Ann_alpha(2,nFCIDet,nNp1FCIDet))
        allocate(Coup_Ann_beta(2,nFCIDet,nNp1bFCIDet))
        Coup_Create_alpha(:,:,:) = 0
        Coup_Create_beta(:,:,:) = 0
        Coup_Ann_alpha(:,:,:) = 0
        Coup_Ann_beta(:,:,:) = 0
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
            !Now, creation of a beta orbital
            do K = 1,nNm1FCIDet
                DiffOrb = ieor(Nm1BitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(mod(nOrbs,2).ne.1) then 
                    !There must be an odd number of orbital differences between the determinants
                    !since they are of different electron number by one
                    call stop_all(t_r,'Not odd number of electrons 2')
                endif
                if(nOrbs.ne.1) then
                    !We only want one orbital different
                    cycle
                endif
                !Now, find out what SPINorbital this one is from the bit number set
                call DecodeBitDet(orbdum,1,DiffOrb)
                gam = orbdum(1)
                if(mod(gam,2).ne.0) then
                    call stop_all(t_r,'differing orbital should be an beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = Nm1BitList(K)
                call SQOperator(tempK,gam,tParity,.false.)
                Coup_Create_beta(1,J,K) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_beta(2,J,K) = -1
                else
                    Coup_Create_beta(2,J,K) = 1
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
            !Now, annihilation of a beta orbital
            do K = 1,nNp1bFCIDet
                DiffOrb = ieor(Np1bBitList(K),FCIBitList(J))
                nOrbs = CountBits(DiffOrb)
                if(mod(nOrbs,2).ne.1) then 
                    !There must be an odd number of orbital differences between the determinants
                    !since they are of different electron number by one
                    call stop_all(t_r,'Not odd number of electrons 4')
                endif
                if(nOrbs.ne.1) then
                    !We only want one orbital different
                    cycle
                endif
                !Now, find out what SPINorbital this one is from the bit number set
                call DecodeBitDet(orbdum,1,DiffOrb)
                gam = orbdum(1)
                if(mod(gam,2).ne.0) then
                    call stop_all(t_r,'differing orbital should be a beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = Np1bBitList(K)
                call SQOperator(tempK,gam,tParity,.true.)
                Coup_Ann_beta(1,J,K) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_beta(2,J,K) = -1
                else
                    Coup_Ann_beta(2,J,K) = 1
                endif
            enddo
        enddo
            
        !Store the fock matrix in complex form, so that we can ZGEMM easily
        write(6,'(A,F9.3,A)') "Memory required to store fock matrix: ",real(nSites**2,dp)*ComptoMb,' Mb'
        write(6,'(A,F9.3,A)') "Memory required to store bath states: ",4.0_dp*real(nSites**2,dp)*ComptoMb,' Mb'
        allocate(FockSchmidtComp(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                FockSchmidtComp(j,i) = dcmplx(FockSchmidt(j,i),0.0_dp)
            enddo
        enddo
        !call R_C_Copy_2D(FockSchmidtComp(:,:),FockSchmidt(:,:),nSites,nSites)

        InterMem = nCore**2 + nVirt**2 + 8*(EmbSize**2) + 2*(EmbSize*nVirt) + 2*(nCore*EmbSize)
        write(6,'(A,F9.3,A)') "Memory required to store contracted intermediates: ",real(InterMem,dp)*ComptoMb,' Mb'
        !Space for useful intermediates
        allocate(G_ai_G_aj(nCore,nCore))    !Core,Core, contracted over virtual space
        allocate(G_ai_G_bi(VirtStart:VirtEnd,VirtStart:VirtEnd))    !Virt,Virt, contracted over core space
        allocate(G_xa_G_ya(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))    !Active,Active, contracted over virtuals
        allocate(G_xa_F_ab(ActiveStart:ActiveEnd,VirtStart:VirtEnd))
        allocate(G_xa_G_yb_F_ab(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ia_G_xa(1:CoreEnd,ActiveStart:ActiveEnd))
        allocate(F_xi_G_ia_G_ya(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_xa_F_ya(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_xi_G_yi(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ix_F_ij(ActiveStart:ActiveEnd,1:CoreEnd))
        allocate(G_ix_G_jy_F_ji(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ia_G_ix(VirtStart:VirtEnd,ActiveStart:ActiveEnd))
        allocate(F_ax_G_ia_G_iy(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))
        allocate(G_ix_F_iy(ActiveStart:ActiveEnd,ActiveStart:ActiveEnd))

        call halt_timer(LR_EC_TDA_Precom)

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            call set_timer(LR_EC_TDA_HBuild)
        
            LinearSystem(:,:) = dcmplx(0.0_dp,0.0_dp)
            Overlap(:,:) = dcmplx(0.0_dp,0.0_dp)
            write(6,*) "Calculating linear response for frequency: ",Omega
            call flush(6)
            if(tMinRes_NonDir) write(minres_unit,*) "Iteratively solving for frequency: ",Omega

            !First, find the non-interacting solution expressed in the schmidt basis
            call FindSchmidtPert(.false.,Omega,ni_lr)

            !call writematrixcomp(SchmidtPert,'SchmidtPert',.true.)
            !write(6,"(A)",advance='no') "Constructing hessian matrix..."

            !First, construct useful intermediates
            !do i=1,nSites
            !    do j=1,nSites
            !        if(abs(SchmidtPert(i,j)-SchmidtPert(j,i)).gt.1.0e-9) then
            !            call stop_all(t_r,'Perturbation not symmetric')
            !        endif
            !    enddo
            !enddo
            !sum_a G_ai^* G_aj
            call ZGEMM('C','N',nCore,nCore,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(VirtStart:VirtEnd,1:CoreEnd),nVirt, &
                SchmidtPert(VirtStart:VirtEnd,1:CoreEnd),nVirt,dcmplx(0.0_dp,0.0_dp),G_ai_G_aj,nCore)
            !sum_i G_ai^* G_bi = sum_i G_ia^* G_bi
            call ZGEMM('C','N',nVirt,nVirt,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,VirtStart:VirtEnd),nCore, &
                SchmidtPert(1:CoreEnd,VirtStart:VirtEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ai_G_bi,nVirt)
            !sum_a G_xa^* G_ya  (where (x,y) in active space)
            call ZGEMM('C','N',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),    &
                SchmidtPert(VirtStart:VirtEnd,ActiveStart:ActiveEnd),nVirt, &
                SchmidtPert(VirtStart:VirtEnd,ActiveStart:ActiveEnd),nVirt,dcmplx(0.0_dp,0.0_dp),G_xa_G_ya,EmbSize)
            !sum_a G_xa F_ab
            call ZGEMM('N','N',EmbSize,nVirt,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(ActiveStart:ActiveEnd,VirtStart:VirtEnd), &
                EmbSize,FockSchmidtComp(VirtStart:VirtEnd,VirtStart:VirtEnd),nVirt,dcmplx(0.0_dp,0.0_dp),G_xa_F_ab,EmbSize)
            !sum_ab G_xa^* G_yb F_ab
            call ZGEMM('C','T',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(VirtStart:VirtEnd,ActiveStart:ActiveEnd), &
                nVirt,G_xa_F_ab(ActiveStart:ActiveEnd,VirtStart:VirtEnd),EmbSize,dcmplx(0.0_dp,0.0_dp),G_xa_G_yb_F_ab,EmbSize)
            !sum_a G_ia^* G_xa
            call ZGEMM('C','T',nCore,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(VirtStart:VirtEnd,1:CoreEnd),nVirt,   &
                SchmidtPert(ActiveStart:ActiveEnd,VirtStart:VirtEnd),EmbSize,dcmplx(0.0_dp,0.0_dp),G_ia_G_xa,nCore)
            !sum_ia F_xi G_ia^* G_ya
            call ZGEMM('N','N',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),FockSchmidtComp(ActiveStart:ActiveEnd,1:CoreEnd),   &
                EmbSize,G_ia_G_xa,nCore,dcmplx(0.0_dp,0.0_dp),F_xi_G_ia_G_ya,EmbSize)
            !sum_a G_xa F_ya
            call ZGEMM('N','T',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),SchmidtPert(ActiveStart:ActiveEnd,VirtStart:VirtEnd), &
                EmbSize,FockSchmidtComp(ActiveStart:ActiveEnd,VirtStart:VirtEnd),EmbSize,dcmplx(0.0_dp,0.0_dp),G_xa_F_ya,EmbSize)
            !sum_i G_ix^* G_iy
            call ZGEMM('C','N',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore, &
                SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_xi_G_yi,EmbSize)
            !sum_i G_ix F_ij
            call ZGEMM('T','N',EmbSize,nCore,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore,   &
                FockSchmidtComp(1:CoreEnd,1:CoreEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ix_F_ij,EmbSize)
            !sum_ij G_ix^* G_jy F_ji
            call ZGEMM('C','T',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore, &
                G_ix_F_ij,EmbSize,dcmplx(0.0_dp,0.0_dp),G_ix_G_jy_F_ji,EmbSize)
            !sum_i G_ia^* G_ix
            call ZGEMM('C','N',nVirt,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,VirtStart:VirtEnd),nCore,   &
                SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ia_G_ix,nVirt)
            !sum_ia F_ax G_ia^* G_iy
            call ZGEMM('T','N',EmbSize,EmbSize,nVirt,dcmplx(1.0_dp,0.0_dp),    &
                FockSchmidtComp(VirtStart:VirtEnd,ActiveStart:ActiveEnd),nVirt,G_ia_G_ix,nVirt, &
                dcmplx(0.0_dp,0.0_dp),F_ax_G_ia_G_iy,EmbSize)
            !sum_i G^ix F_iy
            call ZGEMM('T','N',EmbSize,EmbSize,nCore,dcmplx(1.0_dp,0.0_dp),SchmidtPert(1:CoreEnd,ActiveStart:ActiveEnd),nCore, &
                FockSchmidtComp(1:CoreEnd,ActiveStart:ActiveEnd),nCore,dcmplx(0.0_dp,0.0_dp),G_ix_F_iy,EmbSize)

            
            !First, construct FCI space, in determinant basis
            !****************************   Block 1   **************************
            do i = 1,nFCIDet
                do j = 1,nFCIDet
                    LinearSystem(j,i) = dcmplx(NFCIHam(j,i),0.0_dp)
                enddo
            enddo
            !call R_C_Copy_2D(LinearSystem(1:nFCIDet,1:nFCIDet),NFCIHam(:,:),nFCIDet,nFCIDet)
            !LinearSystem(1:nFCIDet,1:nFCIDet) = NFCIHam(:,:)
            !call writematrixcomp(LinearSystem(1:nFCIDet,1:nFCIDet),'N-electron hamil',.true.)

            !****************************************************************************
            !*********    CORE-VIRTUAL EXCITATION BLOCK *********************************
            !****************************************************************************

            !Calc normalization for the CV block
            CVNorm = 0.0_dp
            do i=1,CoreEnd
                CVNorm = CVNorm + real(G_ai_G_aj(i,i))
            enddo
            CVNorm = CVNorm * 2.0_dp

            !*****************************   Block 2   *****************************
            !Copy the uncontracted FCI hamiltonian space
            LinearSystem(CVIndex:AVIndex-1,CVIndex:AVIndex-1) = LinearSystem(1:nFCIDet,1:nFCIDet)

            !Now alter the diagonals of this block
            tempel = dcmplx(0.0_dp,0.0_dp)
            do i=1,CoreEnd
                do j=1,CoreEnd
                    tempel = tempel - FockSchmidt(j,i)*G_ai_G_aj(j,i)
                enddo
            enddo

            do b = VirtStart,VirtEnd
                do a = VirtStart,VirtEnd
                    tempel = tempel + FockSchmidt(a,b)*G_ai_G_bi(a,b)
                enddo
            enddo
            tempel = tempel * (2.0_dp/CVNorm)

            !write(6,*) "In the CV diagonal space, the diagonals are offset by (should be +ve): ",tempel
            !Offset all diagonals of the CV space by this value
            do i = CVIndex,AVIndex-1
                LinearSystem(i,i) = LinearSystem(i,i) + tempel
            enddo

            !Now for the coupling of the CV excitations to the uncontracted space
            !***********************   Block 3   ***************************
            tempel = dcmplx(0.0_dp,0.0_dp)
            do i = 1,CoreEnd
                do a = VirtStart,VirtEnd
                    tempel = tempel + SchmidtPert(a,i)*FockSchmidt(a,i)
                enddo
            enddo
            tempel = tempel * (2.0_dp/(sqrt(CVNorm)))

            !Add these in to the diagonals of the coupling blocks
            do i = 1,nFCIDet
                LinearSystem(i,nFCIDet+i) = tempel
                LinearSystem(nFCIDet+i,i) = conjg(tempel)
            enddo

            !****************************************************************************
            !*********    ACTIVE-VIRTUAL EXCITATION BLOCK *******************************
            !****************************************************************************

            !First, get normalization constants
            AVNorm(:) = 0.0_dp
            do a = VirtStart,VirtEnd
                do gam = 1,EmbSize
                    gam_spat = gam+CoreEnd
                    !AVNorm(gam) = AVNorm(gam) + real(SchmidtPert(gam_spat,a))**2 + aimag(SchmidtPert(gam_spat,a))**2
                    AVNorm(gam) = AVNorm(gam) + real(conjg(SchmidtPert(gam_spat,a))*SchmidtPert(gam_spat,a),dp)
                enddo
            enddo
            !call writevector(AVNorm,'AV Norm')

            !This starts at 'AVIndex'
            !'Diagonal' block
            !*****************************   Block 4   *****************************
            do gam1 = 1,EmbSize
                !Run through all active space orbitals
                gam1_spat = gam1+(nOcc-nImp)

                !gam1_ind now gives the starting index of the AV diagonal block
                !for orbital gam1
                gam1_ind = AVIndex + (gam1-1)*nFullNm1

                do gam2 = 1,EmbSize
                    gam2_spat = gam2+(nOcc-nImp)

                    gam2_ind = AVIndex + (gam2-1)*nFullNm1

                    !Construct appropriate weighting factor for the hamiltonian matrix element contribution
                    !Fill with appropriate FCI hamiltonian block
                    do i = gam2_ind,gam2_ind+nNm1bFCIDet-1
                        do j = gam1_ind,gam1_ind+nNm1bFCIDet-1
                            LinearSystem(j,i) = dcmplx(Nm1FCIHam_beta(j-gam1_ind+1,i-gam2_ind+1),0.0_dp)
                        enddo
                    enddo
!                    call R_C_Copy_2D(LinearSystem(gam1_ind:gam1_ind+nNm1bFCIDet-1,gam2_ind:gam2_ind+nNm1bFCIDet-1), &
!                        Nm1FCIHam_beta(:,:),nNm1bFCIDet,nNm1bFCIDet)

                    !Now construct for the gam1_beta : gam2_beta block
                    do i = gam2_ind+nNm1bFCIDet,gam2_ind+nFullNm1-1
                        do j = gam1_ind+nNm1bFCIDet,gam1_ind+nFullNm1-1
                            LinearSystem(j,i) =     &
                                dcmplx(Nm1FCIHam_alpha(j-(gam1_ind+nNm1bFCIDet)+1,i-(gam2_ind+nNm1bFCIDet)+1),0.0_dp)
                        enddo
                    enddo
                    !call R_C_Copy_2D(LinearSystem(gam1_ind+nNm1bFCIDet:gam1_ind+nFullNm1-1,   &
                    !    gam2_ind+nNm1bFCIDet:gam2_ind+nFullNm1-1),Nm1FCIHam_alpha(:,:),nNm1FCIDet,nNm1FCIDet)
                    
                    !Multiply every element by the appropriate weighting factor
                    !This weighting factor is the same for both spin blocks, so no need to do seperately
                    LinearSystem(gam1_ind:gam1_ind+nFullNm1-1,gam2_ind:gam2_ind+nFullNm1-1) =   & 
                        LinearSystem(gam1_ind:gam1_ind+nFullNm1-1,gam2_ind:gam2_ind+nFullNm1-1) &
                        *G_xa_G_ya(gam1_spat,gam2_spat)

                    !Now for the virtual excitation term, which is diagonal in each determinant space
                    !Add this term in to the diagonals of each block
                    do i = 0,nFullNm1-1
                        LinearSystem(gam1_ind+i,gam2_ind+i) = LinearSystem(gam1_ind+i,gam2_ind+i) +     &
                            G_xa_G_yb_F_ab(gam1_spat,gam2_spat) 
                    enddo

                    !Finally, all the elements need to be normalized correctly.
                    !Again, the normalization condition is the same for both spins.
                    tempel = sqrt(AVNorm(gam1)*AVNorm(gam2))
                    LinearSystem(gam1_ind:gam1_ind+nFullNm1-1,gam2_ind:gam2_ind+nFullNm1-1) =   &
                        LinearSystem(gam1_ind:gam1_ind+nFullNm1-1,gam2_ind:gam2_ind+nFullNm1-1) / tempel

                enddo   !End gam2
            enddo   !End gam1

            !Now for the coupling to the active-virtual excitation and uncontracted blocks
            !**********************   Block 5 & 6  ***************************
            !First, when the AV excitations are both alpha operators
            do J = 1,nFCIDet
                CVInd_tmp = CVIndex + J - 1 
                do K = 1,nNm1bFCIDet
                    if(Coup_Create_alpha(1,J,K).ne.0) then
                        do beta = 1,EmbSize
                            AVInd_tmp = AVIndex + (beta-1)*nFullNm1 + K -1

                            !Block 5
                            LinearSystem(CVInd_tmp,AVInd_tmp) = LinearSystem(CVInd_tmp,AVInd_tmp) - &
                                F_xi_G_ia_G_ya(Coup_Create_alpha(1,J,K),beta+CoreEnd)*Coup_Create_alpha(2,J,K) / &
                                sqrt(CVNorm*AVNorm(beta))

                            !Hermiticity
                            LinearSystem(AVInd_tmp,CVInd_tmp) = conjg(LinearSystem(CVInd_tmp,AVInd_tmp))

                            !Block 6
                            LinearSystem(J,AVInd_tmp) = LinearSystem(J,AVInd_tmp) + &
                                G_xa_F_ya(beta+CoreEnd,Coup_Create_alpha(1,J,K))*Coup_Create_alpha(2,J,K) / &
                                sqrt(AVNorm(beta))

                            !Hermiticity
                            LinearSystem(AVInd_tmp,J) = conjg(LinearSystem(J,AVInd_tmp))
                        enddo
                    endif
                enddo
                !Now for other spin-type
                do K = 1,nNm1FCIDet
                    if(Coup_Create_beta(1,J,K).ne.0) then
                        do beta = 1,EmbSize
                            AVInd_tmp = AVIndex + (beta-1)*nFullNm1 + nNm1bFCIDet + K -1

                            !Block 5
                            LinearSystem(CVInd_tmp,AVInd_tmp) = LinearSystem(CVInd_tmp,AVInd_tmp) - &
                                F_xi_G_ia_G_ya(Coup_Create_beta(1,J,K),beta+CoreEnd)*Coup_Create_beta(2,J,K) /  &
                                sqrt(CVNorm*AVNorm(beta))

                            !Hermiticity
                            LinearSystem(AVInd_tmp,CVInd_tmp) = conjg(LinearSystem(CVInd_tmp,AVInd_tmp))

                            !Block 6
                            LinearSystem(J,AVInd_tmp) = LinearSystem(J,AVInd_tmp) + &
                                G_xa_F_ya(beta+CoreEnd,Coup_Create_beta(1,J,K))*Coup_Create_beta(2,J,K) / &
                                sqrt(AVNorm(beta))

                            !Hermiticity
                            LinearSystem(AVInd_tmp,J) = conjg(LinearSystem(J,AVInd_tmp))
                        enddo
                    endif
                enddo
            enddo

            !****************************************************************************
            !************    CORE-ACTIVE EXCITATION BLOCK *******************************
            !****************************************************************************

            !First, get normalization constants
            CANorm(:) = 0.0_dp
            do gam = 1,nImp*2
                gam_spat = gam+(nOcc-nImp)
                do i = 1,CoreEnd 
                    CANorm(gam) = CANorm(gam) + real(SchmidtPert(i,gam_spat))**2 + aimag(SchmidtPert(i,gam_spat))**2
                enddo
            enddo
            !call writevector(CANorm,'CA Norm')


            !This starts at 'CAIndex'
            !'Diagonal' block
            !*****************************   Block 7   *****************************
            !As opposed to block 4, the beta hamiltonian now corresponds to the correct block for the excitation
            do gam1 = 1,nImp*2
                !Run through all active space orbitals
                gam1_spat = gam1+(nOcc-nImp)

                !gam1_ind now gives the starting index of the AV diagonal block
                !for orbital gam1
                gam1_ind = CAIndex + (gam1-1)*nFullNp1

                do gam2 = 1,nImp*2
                    gam2_spat = gam2+(nOcc-nImp)

                    gam2_ind = CAIndex + (gam2-1)*nFullNp1

                    !Fill with appropriate FCI hamiltonian block
                    do i = gam2_ind,gam2_ind+nNp1FCIDet-1
                        do j = gam1_ind,gam1_ind+nNp1FCIDet-1
                            LinearSystem(j,i) = dcmplx(Np1FCIHam_alpha(j-gam1_ind+1,i-gam2_ind+1),0.0_dp)
                        enddo
                    enddo
                    !call R_C_Copy_2D(LinearSystem(gam1_ind:gam1_ind+nNp1FCIDet-1,gam2_ind:gam2_ind+nNp1FCIDet-1),   &
                    !    Np1FCIHam_alpha(:,:),nNp1FCIDet,nNp1FCIDet)

                    !Now construct for the gam1_beta : gam2_beta block
                    do i = gam2_ind+nNp1FCIDet,gam2_ind+nFullNp1-1
                        do j = gam1_ind+nNp1FCIDet,gam1_ind+nFullNp1-1
                            LinearSystem(j,i) = dcmplx(Np1FCIHam_beta(j-(gam1_ind+nNp1FCIDet)+1,i-(gam2_ind+nNp1FCIDet)+1),0.0_dp)
                        enddo
                    enddo
                    !call R_C_Copy_2D(LinearSystem(gam1_ind+nNp1FCIDet:gam1_ind+nFullNp1-1,    &
                    !    gam2_ind+nNp1FCIDet:gam2_ind+nFullNp1-1),Np1FCIHam_beta(:,:),nNp1bFCIDet,nNp1bFCIDet)
                    
                    !Multiply every element by the appropriate weighting factor
                    !This weighting factor is the same for both spin blocks, so no need to do seperately
                    LinearSystem(gam1_ind:gam1_ind+nFullNp1-1,gam2_ind:gam2_ind+nFullNp1-1) =   &
                        LinearSystem(gam1_ind:gam1_ind+nFullNp1-1,gam2_ind:gam2_ind+nFullNp1-1) *   &
                        G_xi_G_yi(gam1_spat,gam2_spat)

                    !Now for the occupied excitation term, which is diagonal in each determinant space
                    !Add this term in to the diagonals of each block
                    do i = 0,nFullNp1-1
                        LinearSystem(gam1_ind+i,gam2_ind+i) = LinearSystem(gam1_ind+i,gam2_ind+i) - &
                            G_ix_G_jy_F_ji(gam1_spat,gam2_spat)
                    enddo

                    !Finally, all the elements need to be normalized correctly.
                    !Again, the normalization condition is the same for both spins.
                    tempel = sqrt(CANorm(gam1)*CANorm(gam2))
                    LinearSystem(gam1_ind:gam1_ind+nFullNp1-1,gam2_ind:gam2_ind+nFullNp1-1) =   &
                        LinearSystem(gam1_ind:gam1_ind+nFullNp1-1,gam2_ind:gam2_ind+nFullNp1-1) / tempel
                enddo   !End gam2
            enddo   !End gam1

            !*************************   Block 8 is ZERO   ******************************

            !*************************   Block 9 & 10  **************************************
            do J = 1,nFCIDet
                CVInd_tmp = CVIndex + J - 1 
                do K = 1,nNp1FCIDet
                    if(Coup_Ann_alpha(1,J,K).ne.0) then
                        do beta = 1,EmbSize
                            CAInd_tmp = CAIndex + (beta-1)*nFullNp1 + K - 1

                            !Block 9
                            LinearSystem(CVInd_tmp,CAInd_tmp) = LinearSystem(CVInd_tmp,CAInd_tmp) - &
                                F_ax_G_ia_G_iy(Coup_Ann_alpha(1,J,K),beta+CoreEnd)*Coup_Ann_alpha(2,J,K) / &
                                sqrt(CVNorm*CANorm(beta))

                            !Hermiticity
                            LinearSystem(CAInd_tmp,CVInd_tmp) = conjg(LinearSystem(CVInd_tmp,CAInd_tmp))

                            !Block 10
                            LinearSystem(J,CAInd_tmp) = LinearSystem(J,CAInd_tmp) - &
                                G_ix_F_iy(beta+CoreEnd,Coup_Ann_alpha(1,J,K))*Coup_Ann_alpha(2,J,K) / &
                                sqrt(CANorm(beta))

                            !Hermiticity
                            LinearSystem(CAInd_tmp,J) = conjg(LinearSystem(J,CAInd_tmp))
                        enddo
                    endif
                enddo
                !Now for beta spin
                do K = 1,nNp1bFCIDet
                    if(Coup_Ann_beta(1,J,K).ne.0) then
                        do beta = 1,EmbSize
                            CAInd_tmp = CAIndex + (beta-1)*nFullNp1 + nNp1FCIDet + K - 1

                            !Block 9
                            LinearSystem(CVInd_tmp,CAInd_tmp) = LinearSystem(CVInd_tmp,CAInd_tmp) - &
                                F_ax_G_ia_G_iy(Coup_Ann_beta(1,J,K),beta+CoreEnd)*Coup_Ann_beta(2,J,K) / &
                                sqrt(CVNorm*CANorm(beta))

                            !Hermiticity
                            LinearSystem(CAInd_tmp,CVInd_tmp) = conjg(LinearSystem(CVInd_tmp,CAInd_tmp))

                            !Block 10
                            LinearSystem(J,CAInd_tmp) = LinearSystem(J,CAInd_tmp) - &
                                G_ix_F_iy(beta+CoreEnd,Coup_Ann_beta(1,J,K))*Coup_Ann_beta(2,J,K) / &
                                sqrt(CANorm(beta))

                            !Hermiticity
                            LinearSystem(CAInd_tmp,J) = conjg(LinearSystem(J,CAInd_tmp))
                        enddo
                    endif
                enddo
            enddo

!            !Now check that Hessian is hermitian
!            do i=1,nLinearSystem
!                do j=i,nLinearSystem
!                    if(abs(LinearSystem(i,j)-conjg(LinearSystem(j,i))).gt.1.0e-8_dp) then
!                        write(6,*) "i, j: ",i,j
!                        write(6,*) "LinearSystem(i,j): ",LinearSystem(i,j)
!                        write(6,*) "LinearSystem(j,i): ",LinearSystem(j,i)
!                        call stop_all(t_r,'Hessian for EC-LR not hermitian')
!                    endif
!                enddo
!            enddo

            !write(6,*) "Hessian constructed successfully...",Omega
            call halt_timer(LR_EC_TDA_HBuild)

            !*********************   Hessian construction finished   **********************
            !call writematrix(LinearSystem(1:nFCIDet,1:nFCIDet),'Hessian_N',.true.)

            !****************************************************************************
            !************    OVERLAP MATRIX   *******************************************
            !****************************************************************************
            call set_timer(LR_EC_TDA_SBuild)

            ! Block 1 and 2 are equal to the identity
            do i = 1,AVIndex-1
                Overlap(i,i) = dcmplx(1.0_dp,0.0_dp)
            enddo

            !Now deal with block 4
            !AV-AV overlap
            do gam1 = 1,nImp*2
                !Run through all active space orbitals
                gam1_spat = gam1+(nOcc-nImp)

                !gam1_ind now gives the starting index of the AV diagonal block
                !for orbital gam1
                gam1_ind = AVIndex + (gam1-1)*nFullNm1

                do gam2 = 1,nImp*2
                    gam2_spat = gam2+(nOcc-nImp)

                    gam2_ind = AVIndex + (gam2-1)*nFullNm1

                    !Now for the overlap, which is diagonal in each determinant space
                    tempel = G_xa_G_ya(gam1_spat,gam2_spat) / sqrt(AVNorm(gam1)*AVNorm(gam2))
                    if((gam1.eq.gam2).and.(abs(tempel-1.0_dp).gt.1.0e-7_dp)) then
                        write(6,*) "gam1,gam2: ",gam1,gam2
                        write(6,*) "tempel: ",tempel
                        call stop_all(t_r,'Error calculating overlap 1')
                    endif
                    !Add this term in to the diagonals of each block
                    do i = 0,nFullNm1-1
                        Overlap(gam1_ind+i,gam2_ind+i) = tempel
                    enddo
                enddo   !End gam2
            enddo   !End gam1

            !Now deal with block 7
            !CA-CA overlap
            do gam1 = 1,nImp*2
                !Run through all active space orbitals
                gam1_spat = gam1+(nOcc-nImp)

                !gam1_ind now gives the starting index of the AV diagonal block
                !for orbital gam1
                gam1_ind = CAIndex + (gam1-1)*nFullNp1

                do gam2 = 1,nImp*2
                    gam2_spat = gam2+(nOcc-nImp)

                    gam2_ind = CAIndex + (gam2-1)*nFullNp1

                    !Now for the occupied excitation term, which is diagonal in each determinant space
                    tempel = G_xi_G_yi(gam1_spat,gam2_spat) / sqrt(CANorm(gam1)*CANorm(gam2))
                    if((gam1.eq.gam2).and.(abs(tempel-1.0_dp).gt.1.0e-7_dp)) then
                        call stop_all(t_r,'Error calculating overlap 2')
                    endif
                    !Add this term in to the diagonals of each block
                    do i = 0,nFullNp1-1
                        Overlap(gam1_ind+i,gam2_ind+i) = tempel
                    enddo

                enddo   !End gam2
            enddo   !End gam1
!            !Check hermiticity and normalization of overlap matrix
            do i=1,nLinearSystem
                do j=i,nLinearSystem
                    if(abs(Overlap(i,j)-conjg(Overlap(j,i))).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Overlap matrix not hermitian')
                    endif
                enddo
                if(abs(Overlap(i,i)-1.0_dp).gt.1.0e-7_dp) then
                    write(6,*) "i: ",i
                    write(6,*) "Overlap(i,i): ",Overlap(i,i)
                    call stop_all(t_r,'Functions not normalized')
                endif
            enddo

!            write(6,*) "Overlap matrix constructed successfully..."
            call halt_timer(LR_EC_TDA_SBuild)
            call set_timer(LR_EC_TDA_Project)

            if(tProjectOutNull.or.tLR_ReoptGS.or.tOrthogBasis) then
                !We need eigenvalues and vectors of S
                !Beware, that by block diagonalizing, the eigenvalues are no longer in increasing order

                !Fill in orthogonal blocks (no need to diagonalize)
                allocate(S_EigVec(nLinearSystem,nLinearSystem))
                allocate(S_Eigval(nLinearSystem))
                S_EigVec(:,:) = dcmplx(0.0_dp,0.0_dp)
                do i=1,AVIndex-1
                    S_EigVec(i,i) = dcmplx(1.0_dp,0.0_dp)
                    S_Eigval(i) = 1.0_dp
                enddo

                !Diagonalize block 4
                nSize = nFullNm1 * EmbSize
                if(nSize.ne.(CAIndex-AVIndex)) then
                    write(6,*) "nSize: ",nSize
                    write(6,*) "CVIndex-AVIndex: ",CAIndex-AVIndex
                    call stop_all(t_r,'Block sizes wrong')
                endif
                allocate(SBlock(nSize,nSize))
                SBlock(:,:) = Overlap(AVIndex:CAIndex-1,AVIndex:CAIndex-1)
                allocate(Work(max(1,3*nSize-2)))
                allocate(temp_vecc(1))
                allocate(SBlock_val(nSize))
                SBlock_val(:) = 0.0_dp
                lWork = -1
                info = 0
                call zheev('V','U',nSize,SBlock,nSize,SBlock_val,temp_vecc,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'workspace query failed')
                lwork = int(temp_vecc(1))+1
                deallocate(temp_vecc)
                allocate(temp_vecc(lwork))
                call zheev('V','U',nSize,SBlock,nSize,SBlock_val,temp_vecc,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'S Diag failed 1')
                deallocate(work,temp_vecc)
                !Fill into main arrays
                S_Eigval(AVIndex:CAIndex-1) = SBlock_val(:)
                S_EigVec(AVIndex:CAIndex-1,AVIndex:CAIndex-1) = SBlock(:,:)
                deallocate(SBlock,SBlock_val)

                !Block 7
                nSize = nFullNp1 * EmbSize
                if(nSize.ne.(nLinearSystem-CAIndex+1)) call stop_all(t_r,'Block sizes wrong')
                allocate(SBlock(nSize,nSize))
                SBlock(:,:) = Overlap(CAIndex:nLinearSystem,CAIndex:nLinearSystem)
                allocate(Work(max(1,3*nSize-2)))
                allocate(temp_vecc(1))
                allocate(SBlock_val(nSize))
                SBlock_val(:) = 0.0_dp
                lWork = -1
                info = 0
                call zheev('V','U',nSize,SBlock,nSize,SBlock_val,temp_vecc,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'workspace query failed')
                lwork = int(temp_vecc(1))+1
                deallocate(temp_vecc)
                allocate(temp_vecc(lwork))
                call zheev('V','U',nSize,SBlock,nSize,SBlock_val,temp_vecc,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'S Diag failed 2')
                deallocate(work,temp_vecc)
                !Fill into main arrays
                S_Eigval(CAIndex:nLinearSystem) = SBlock_val(:)
                S_EigVec(CAIndex:nLinearSystem,CAIndex:nLinearSystem) = SBlock(:,:)
                deallocate(SBlock,SBlock_val)

                !call writevector(S_Eigval,'Overlap EVals')
            endif

            if(tProjectOutNull) then
                !Project out the null space from the equations to ensure unique representation.
                !Work in the linear span of S
                !S_EigVec contains the eigenvectors of the overlap matrix. Only keep ones with large enough eigenvalues
                nSpan = 0
                do i=1,nLinearSystem
                    if(S_Eigval(i).lt.(-1.0e-8_dp)) then
                        call stop_all(t_r,'Error - shouldnt have negative eigenvalues in overlap spectrum')
                    elseif(S_Eigval(i).gt.MinS_Eigval) then
                        !Include this eigenvector
                        nSpan = nSpan + 1
                    endif
                enddo
                if(nSpan.eq.nLinearSystem) then
                    write(6,*) "No linear dependencies found in overlap"
                    tTransformSpace = .false.
                    allocate(Transform(nLinearSystem,nLinearSystem))
                    Transform(:,:) = S_EigVec(:,:)
                    allocate(S_Diag(nSpan))
                    S_Diag(:) = S_Eigval(:)
                else
                    tTransformSpace = .true.
                    write(6,*) "Removing linear dependencies in overlap. Vectors removed: ",nLinearSystem-nSpan

                    !Find transformation matrix to the new space with linear dependencies removed before solving
                    allocate(Transform(nLinearSystem,nSpan))
                    allocate(S_Diag(nSpan))
                    Transform(:,:) = dcmplx(0.0_dp,0.0_dp)
                    S_Diag(:) = dcmplx(0.0_dp,0.0_dp)
                    nSpan = 0
                    do i=1,nLinearSystem
                        if(S_Eigval(i).gt.MinS_Eigval) then
                            nSpan = nSpan + 1
                            !Include vector
                            Transform(:,nSpan) = S_EigVec(:,i)
                            S_Diag(nSpan) = S_Eigval(i)
                        endif
                    enddo

                endif
            else
                !Do not project out null space
                tTransformSpace = .false.
                nSpan = nLinearSystem
                if(tLR_ReoptGS) then
                    allocate(Transform(nLinearSystem,nLinearSystem))
                    Transform(:,:) = S_EigVec(:,:)
                endif
            endif
            call halt_timer(LR_EC_TDA_Project)

            call set_timer(LR_EC_TDA_OptGS)
            if(tOrthogBasis) then
                allocate(Psi_0(nSpan))  !To store the ground state in the full basis
                Psi_0(:) = dcmplx(0.0_dp,0.0_dp)
            else
                allocate(Psi_0(nLinearSystem))  !To store the ground state in the full basis
                Psi_0(:) = dcmplx(0.0_dp,0.0_dp)
            endif

            if(tLR_ReoptGS) then
                !Reoptimize the GS in this new space

                !Transform to the canonically orthogonalized representation in the non-null space of S
                allocate(CanTrans(nLinearSystem,nSpan))
                allocate(tempc(nSpan,nSpan))
                tempc(:,:) = dcmplx(0.0_dp,0.0_dp)
                j = 1
                do i=1,nLinearSystem
                    if(S_EigVal(i).gt.MinS_Eigval) then
                        tempc(j,j) = dcmplx(1.0_dp/sqrt(S_Eigval(i)),0.0_dp)
                        j=j+1
                    elseif(.not.tProjectOutNull) then
                        write(6,*) "Small overlap eigenvalue: ",S_EigVal(i)
                        call stop_all(t_r,  &
                            'Cannot reoptimize ground state without projecting out null space due to small overlap eigenvals')
                    endif
                enddo
                if(j.ne.(nSpan+1)) call stop_all(t_r,'Error in indexing when reoptimizing ground state')

                call ZGEMM('N','N',nLinearSystem,nSpan,nSpan,dcmplx(1.0_dp,0.0_dp),    &
                    Transform,nLinearSystem,tempc,nSpan,dcmplx(0.0_dp,0.0_dp),CanTrans,nLinearSystem)

                !Now, transform that hamiltonian into this new basis
                allocate(OrthogHam(nSpan,nSpan))
                deallocate(tempc)
                allocate(tempc(nLinearSystem,nSpan))
                !call writematrixcomp(CanTrans,'CanTrans',.true.)
                call ZGEMM('N','N',nLinearSystem,nSpan,nLinearSystem,dcmplx(1.0_dp,0.0_dp),    &
                    LinearSystem,nLinearSystem,CanTrans,nLinearSystem,dcmplx(0.0_dp,0.0_dp),tempc,nLinearSystem)
                call ZGEMM('C','N',nSpan,nSpan,nLinearSystem,dcmplx(1.0_dp,0.0_dp),    &
                    CanTrans,nLinearSystem,tempc,nLinearSystem,dcmplx(0.0_dp,0.0_dp),OrthogHam,nSpan)
                !call writematrixcomp(OrthogHam,'OrthogHam',.true.)

                !Rediagonalize this new hamiltonian
                allocate(Work(max(1,3*nSpan-2)))
                allocate(temp_vecc(1))
                allocate(H_Vals(nSpan))
                H_Vals(:)=0.0_dp
                lWork=-1
                info=0
                call zheev('V','U',nSpan,OrthogHam,nSpan,H_Vals,temp_vecc,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'workspace query failed')
                lwork=int(abs(temp_vecc(1)))+1
                deallocate(temp_vecc)
                allocate(temp_vecc(lwork))
                call zheev('V','U',nSpan,OrthogHam,nSpan,H_Vals,temp_vecc,lWork,Work,info)
                if (info.ne.0) then
                    write(6,*) "info: ",info
                    call stop_all(t_r,"Diag failed 2")
                endif
                deallocate(work,temp_vecc)

                GSEnergy = H_Vals(1)
                write(6,"(A,G20.10,A,G20.10,A)") "Reoptimized ground state energy is: ",GSEnergy, &  
                    " (old = ",HL_Energy,")"
                if(tDiagHam) write(iunit2,*) Omega,H_Vals(:)
                if(GSEnergy-HL_Energy.gt.1.0e-8_dp) then
                    call stop_all(t_r,'Reoptimized GS energy lower than original GS energy - this should not be possible')
                endif

                if(.not.tOrthogBasis) then
                    !Store the eigenvector in the original basis
                    !However, just rotate the ground state, not all the eigenvectors
                    call ZGEMV('N',nLinearSystem,nSpan,dcmplx(1.0_dp,0.0_dp),CanTrans,nLinearSystem,   &
                        OrthogHam(:,1),1,dcmplx(0.0_dp,0.0_dp),Psi_0,1)
                else
                    !Do not rotate - keep in the orthogonal basis
                    Psi_0(:) = OrthogHam(:,1)
                endif

                deallocate(OrthogHam,H_Vals,tempc,CanTrans)
            else
                !We are not reoptimizing the GS
                GSEnergy = HL_Energy
                do i = 1,nFCIDet
                    Psi_0(i) = dcmplx(HL_Vec(i),0.0_dp)
                enddo
            endif
            call halt_timer(LR_EC_TDA_OptGS)

            call set_timer(LR_EC_TDA_BuildLR)
            !Remove the ground state from the hamiltonian
            if(tOrthogBasis) then
                allocate(Projector(nSpan,nSpan))
                allocate(temp_vecc_2(nSpan))
                do j = 1,nSpan
                    do i = 1,nSpan
                        Projector(i,j) = Psi_0(i)*conjg(Psi_0(j))*S_Diag(j)
                    enddo
                    temp_vecc_2(j) = conjg(Psi_0(j))*S_Diag(j)
                enddo
            else
                allocate(Projector(nLinearSystem,nLinearSystem))
                allocate(temp_vecc_2(nLinearSystem))
                allocate(temp_vecc_3(nLinearSystem))
                temp_vecc_3(:) = conjg(Psi_0(:))
                call ZGEMV('T',nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),Overlap,nLinearSystem,    &
                    temp_vecc_3,1,dcmplx(0.0_dp,0.0_dp),temp_vecc_2,1)
                deallocate(temp_vecc_3)
                Projector(:,:) = dcmplx(0.0_dp,0.0_dp)
                do i = 1,nLinearSystem
                    do j = 1,nLinearSystem
                        Projector(i,j) = temp_vecc_2(j)*Psi_0(i)
                    enddo
                enddo
            endif

            !Construct LHS
            allocate(LHS(nSpan,nSpan))
            if(tOrthogBasis) then
                !Rotate linear system into the orthogonal basis.
                allocate(tempc(nSpan,nLinearSystem))
                call ZGEMM('C','N',nSpan,nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),Transform,nLinearSystem,    &
                    LinearSystem,nLinearSystem,dcmplx(0.0_dp,0.0_dp),tempc,nSpan)
                call ZGEMM('N','N',nSpan,nSpan,nLinearSystem,dcmplx(1.0_dp,0.0_dp),tempc,nSpan,Transform,nLinearSystem,    &
                    dcmplx(0.0_dp,0.0_dp),LHS,nSpan)

                !Rotated overlap should be diagonal matrix
                !Check this
                allocate(CanTrans(nSpan,nSpan))
                call ZGEMM('C','N',nSpan,nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),Transform,nLinearSystem,    &
                    Overlap,nLinearSystem,dcmplx(0.0_dp,0.0_dp),tempc,nSpan)
                call ZGEMM('N','N',nSpan,nSpan,nLinearSystem,dcmplx(1.0_dp,0.0_dp),tempc,nSpan,Transform,nLinearSystem,    &
                    dcmplx(0.0_dp,0.0_dp),CanTrans,nSpan)
                do i=1,nSpan
                    do j=1,nSpan
                        if((i.eq.j).and.(abs(CanTrans(i,i)-S_Diag(i)).gt.1.0e-8_dp)) then
                            write(6,*) "S: ",i,j,CanTrans(i,j)
                            call stop_all(t_r,'Projected overlap matrix does not have unit diagonal')
                        elseif((i.ne.j).and.(abs(CanTrans(j,i)).gt.1.0e-8_dp)) then
                            write(6,*) "S: ",i,j,CanTrans(j,i)
                            call stop_all(t_r,'Projected overlap matrix has non-zero off-diagonals')
                        endif
                    enddo
                enddo
                deallocate(CanTrans,tempc)
                
                if(tRemoveGSFromH) then
                    !Remove ground state from linear system
                    LHS(:,:) = LHS(:,:) - (Projector(:,:)*GSEnergy)
                endif

                !Form LHS - easier as S is diagonal   
                do i = 1,nSpan
                    LHS(i,i) = LHS(i,i) - S_Diag(i)*(GSEnergy + dcmplx(Omega,dDelta))
                enddo
            else
                if(tRemoveGSFromH) then
                    !Remove ground state from linear system
                    LinearSystem(:,:) = LinearSystem(:,:) - (Projector(:,:)*GSEnergy)
                endif

                !Now construct the lhs of the equations
                !Now, we want to calculate H - (E_0 + Omega)S
                do i=1,nLinearSystem
                    do j=1,nLinearSystem
                        LinearSystem(j,i) = LinearSystem(j,i) - (Overlap(j,i)*(GSEnergy + dcmplx(Omega,dDelta)))
                    enddo
                enddo

                if(tTransformSpace) then
                    !Transform the LHS to the linear span of S
                    allocate(tempc(nLinearSystem,nSpan))
                    call ZGEMM('N','N',nLinearSystem,nSpan,nLinearSystem,dcmplx(1.0_dp,0.0_dp),LinearSystem,nLinearSystem, &
                        Transform,nLinearSystem,dcmplx(0.0_dp,0.0_dp),tempc,nLinearSystem)
                    call ZGEMM('C','N',nSpan,nSpan,nLinearSystem,dcmplx(1.0_dp,0.0_dp),Transform,nLinearSystem,    &
                        tempc,nLinearSystem,dcmplx(0.0_dp,0.0_dp),LHS,nSpan)
                    deallocate(tempc)
                else
                    LHS(:,:) = LinearSystem(:,:)
                endif
            endif

            !Set up RHS of linear system: -Q V |0>

            !Change projector so that it now projects *out* the ground state
            Projector(:,:) = Projector(:,:)*dcmplx(-1.0_dp,0.0_dp)
            if(tOrthogBasis) then
                do i = 1,nSpan
                    Projector(i,i) = dcmplx(1.0_dp,0.0_dp) + Projector(i,i)
                enddo
            else
                do i = 1,nLinearSystem
                    Projector(i,i) = dcmplx(1.0_dp,0.0_dp) + Projector(i,i)
                enddo
            endif
            
            !find V |0> and store it in temp_vecc
            if(tOrthogBasis) then
                allocate(temp_vecc(nSpan))
                call ApplyDensityPert_EC_Orthog(Psi_0,temp_vecc,nSpan,nLinearSystem,Transform)
            else
                allocate(temp_vecc(nLinearSystem))
                call ApplyDensityPert_EC(Psi_0,temp_vecc,nLinearSystem)
            endif
            !call writevectorcomp(Psi_0,'Psi_0')
            !call writevectorcomp(temp_vecc,'V|0>')

            !Now, calculate -QV|0> and put into VGS
            if(tOrthogBasis) then
                allocate(VGS(nSpan))
                call ZGEMV('N',nSpan,nSpan,dcmplx(-1.0_dp,0.0_dp),Projector,nSpan, &
                    temp_vecc,1,dcmplx(0.0_dp,0.0_dp),VGS,1)

                !Check that RHS is orthogonal to the ground state
                testc = dcmplx(0.0_dp,0.0_dp)
                do j = 1,nSpan
                    testc = testc + temp_vecc_2(j)*VGS(j) !Overlap(i,j)*conjg(Psi_0(i))*VGS(j)
                enddo
            else
                allocate(VGS(nLinearSystem))
                call ZGEMV('N',nLinearSystem,nLinearSystem,dcmplx(-1.0_dp,0.0_dp),Projector,nLinearSystem, &
                    temp_vecc,1,dcmplx(0.0_dp,0.0_dp),VGS,1)

                !Check that RHS is orthogonal to the ground state
                testc = dcmplx(0.0_dp,0.0_dp)
                do j = 1,nLinearSystem
                    testc = testc + temp_vecc_2(j)*VGS(j) !Overlap(i,j)*conjg(Psi_0(i))*VGS(j)
                enddo
            endif

            !We now have the RHS in 'VGS'. Project this into the linear span of S if needed
            allocate(RHS(nSpan))
            if(tTransformSpace.and.(.not.tOrthogBasis)) then
                call ZGEMV('C',nLinearSystem,nSpan,dcmplx(1.0_dp,0.0_dp),Transform,nLinearSystem,  &
                    VGS,1,dcmplx(0.0_dp,0.0_dp),RHS,1)
            else
                RHS(:) = VGS(:)
            endif
            call halt_timer(LR_EC_TDA_BuildLR)

            call set_timer(LR_EC_TDA_SolveLR)

            !Check hamiltonian is hermitian (on non-diagonals)
            if((.not.tRemoveGSFromH).and.tOrthogBasis) then
                do i=1,nSpan
                    do j=i+1,nSpan
                        if(abs(LHS(i,j)-conjg(LHS(j,i))).gt.1.0e-7_dp) then
                            write(6,*) i,j,LHS(i,j),LHS(j,i),abs(LHS(i,j)-conjg(LHS(j,i)))
                            call warning(t_r,'Linear system off-diagonal hermiticity lost')
                        endif
                    enddo
                enddo
            endif

            !Now solve these linear equations
            if(tMinRes_NonDir) then
                zDirMV_Mat => LHS
                maxminres_iter_ip = int(maxminres_iter,ip)
                minres_unit_ip = int(minres_unit,ip)
                nSpan_ip = int(nSpan,ip)
                allocate(RHS2(nSpan))
                allocate(Psi1(nSpan))
                call setup_RHS(nSpan,RHS,RHS2)
                if(tPrecond_MinRes) then
                    allocate(Precond_Diag(nSpan))
                    call FormPrecond(nSpan)
                    call MinResQLP(n=nSpan_ip,Aprod=zDirMV,b=RHS2,nout=minres_unit_ip,x=Psi1, &
                        itnlim=maxminres_iter_ip,Msolve=zPreCond,istop=info_ip,rtol=rtol_LR,itn=iters, &
                        startguess=tReuse_LS)
                else
                    call MinResQLP(n=nSpan_ip,Aprod=zDirMV,b=RHS2,nout=minres_unit_ip,x=Psi1, &
                        itnlim=maxminres_iter_ip,istop=info_ip,rtol=rtol_LR,itn=iters,startguess=tReuse_LS)
                endif
                info = info_ip
                zDirMV_Mat => null()
                if(info.gt.7) write(6,*) "info: ",info
                if(info.eq.8) write(6,"(A,I9)") "Linear equation solver hit iteration limit: ",iters
                if((info.eq.9).or.(info.eq.10).or.(info.eq.11)) then
                    call stop_all(t_r,'Input matrices to linear solver incorrect')
                endif
                if(info.gt.11) call stop_all(t_r,'Linear equation solver failed')
                info = 0
                !Copy solution to RHS
                RHS(:) = Psi1(:)
                deallocate(RHS2,Psi1)
                if(tPrecond_MinRes) deallocate(Precond_Diag)
            else
                call SolveCompLinearSystem(LHS,RHS,nSpan,info)
            endif
            if(info.ne.0) then 
                write(6,*) "INFO: ",info
                !call stop_all(t_r,'Solving Linear system failed') 
                call warning(t_r,'Solving linear system failed - skipping this frequency')
            else

                if(tTransformSpace.and.(.not.tOrthogBasis)) then
                    !Expand back out the perturbed wavefunction into the original basis
                    call ZGEMV('N',nLinearSystem,nSpan,dcmplx(1.0_dp,0.0_dp),Transform,nLinearSystem,  &
                        RHS,1,dcmplx(0.0_dp,0.0_dp),VGS,1)
                else
                    VGS(:) = RHS(:)
                endif

                if(tExplicitlyOrthog) then
                    if(tOrthogBasis) then
                        !Explicitly project out any component on psi_0
                        call ZGEMV('N',nSpan,nSpan,dcmplx(1.0_dp,0.0_dp),Projector,nSpan,VGS,1,    &
                            dcmplx(0.0_dp,0.0_dp),temp_vecc,1)
                        VGS(:) = temp_vecc(:)
                    else
                        !Explicitly project out any component on psi_0
                        call ZGEMV('N',nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),Projector,nLinearSystem,VGS,1,    &
                            dcmplx(0.0_dp,0.0_dp),temp_vecc,1)
                        VGS(:) = temp_vecc(:)
                    endif
                endif

                !Find normalization of first-order wavefunction
                !Test whether the perturbed wavefunction is orthogonal to the ground state wavefunction
                dNorm = dcmplx(0.0_dp,0.0_dp)
                dOrthog = dcmplx(0.0_dp,0.0_dp)
                if(tOrthogBasis) then
                    do j = 1,nSpan
                        dNorm = dNorm + S_Diag(j)*conjg(VGS(j))*VGS(j)
                        dOrthog = dOrthog + temp_vecc_2(j)*VGS(j)
                    enddo
                else
                    do j = 1,nLinearSystem
                        do i = 1,nLinearSystem
                            dNorm = dNorm + Overlap(i,j)*conjg(VGS(i))*VGS(j)
                        enddo
                        dOrthog = dOrthog + temp_vecc_2(j)*VGS(j)
                    enddo
                endif
!                write(6,*) "Normalization of first-order wavefunction: ",dNorm
!                write(6,*) "Orthogonality of first-order wavefunction: ",dOrthog

                tCalcResponse = .true.
                !write(20,*) U,Omega,abs(dOrthog),abs(dNorm),real(dOrthog/dNorm),aimag(dOrthog/dNorm),abs(dOrthog/dNorm)
                if(abs(dOrthog).gt.1.0e-4_dp) then
                    !call warning(t_r,'First order wavefunction not orthogonal - skipping this frequency')
                    call warning(t_r,'First order wavefunction not orthogonal')
                    !tCalcResponse=.false.
                elseif(abs(dOrthog).gt.1.0e-8_dp) then
                    call warning(t_r,'First order solution not strictly orthogonal')
                endif

                if(tCalcResponse) then
                    !Now, use temp_vecc to hold V|1>
                    ResponseFn = dcmplx(0.0_dp,0.0_dp)
                    if(tOrthogBasis) then
                        call ApplyDensityPert_EC_Orthog(VGS,temp_vecc,nSpan,nLinearSystem,Transform)
                        do j = 1,nSpan
                            ResponseFn = ResponseFn + temp_vecc(j)*temp_vecc_2(j)
                        enddo
                    else
                        call ApplyDensityPert_EC(VGS,temp_vecc,nLinearSystem)
                        do j = 1,nLinearSystem
                            ResponseFn = ResponseFn + temp_vecc(j)*temp_vecc_2(j)
                        enddo
                    endif
                    write(iunit,"(10G22.10,I8)") Omega,real(ResponseFn),-aimag(ResponseFn), &
                        abs(dOrthog),abs(dNorm),GSEnergy,HL_Energy,abs(testc),real(ni_lr),-aimag(ni_lr),iters
                    call flush(iunit)
                    deallocate(temp_vecc)
                endif

            endif

            if(allocated(Transform)) deallocate(Transform)
            if(allocated(S_Diag)) deallocate(S_Diag)
            if(allocated(S_EigVec)) then
                deallocate(S_EigVec)
                deallocate(S_EigVal)
            endif
            deallocate(LHS,RHS,Psi_0,Projector,temp_vecc_2,VGS)

            Omega = Omega + Omega_Step

            call halt_timer(LR_EC_TDA_SolveLR)
        enddo   !End loop over omega

        deallocate(LinearSystem,Overlap)
        close(iunit)
        if(tDiagHam) close(iunit2)
        if(tMinRes_NonDir) then
            close(minres_unit)
        endif

        !Deallocate determinant lists
        if(allocated(FCIDetList)) deallocate(FCIDetList)
        if(allocated(FCIBitList)) deallocate(FCIBitList)
        if(allocated(UMat)) deallocate(UMat)
        if(allocated(TMat)) deallocate(TMat)
        if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
        if(allocated(Nm1BitList)) deallocate(Nm1BitList)
        if(allocated(Np1FCIDetList)) deallocate(Np1FCIDetList)
        if(allocated(Np1BitList)) deallocate(Np1BitList)
        if(allocated(Nm1bFCIDetList)) deallocate(Nm1bFCIDetList)
        if(allocated(Nm1bBitList)) deallocate(Nm1bBitList)
        if(allocated(Np1bFCIDetList)) deallocate(Np1bFCIDetList)
        if(allocated(Np1bBitList)) deallocate(Np1bBitList)

        !Stored intermediates
        deallocate(NFCIHam,Nm1FCIHam_alpha,Nm1FCIHam_beta,Np1FCIHam_alpha,Np1FCIHam_beta)
        deallocate(AVNorm,CANorm)
        deallocate(Coup_Create_alpha,Coup_Create_beta,Coup_Ann_alpha,Coup_Ann_beta)
        deallocate(FockSchmidtComp,G_ai_G_aj,G_ai_G_bi,G_xa_G_ya,G_xa_F_ab,G_xa_G_yb_F_ab)
        deallocate(G_ia_G_xa,F_xi_G_ia_G_ya,G_xa_F_ya,G_xi_G_yi,G_ix_F_ij,G_ix_G_jy_F_ji)
        deallocate(G_ia_G_ix,F_ax_G_ia_G_iy,G_ix_F_iy)


    end subroutine NonIntExContracted_TDA_MCLR
        
    !This will apply an alpha creation operator at site pertsite to the first order interacting wavefunction
    !Returns the wavefunction Psi in the standard N-electron determinant space
    subroutine ApplyCre_FirstOrder_EC(Psi_1,nSize1,Psi,nSize0,pertsite)
        use DetToolsData
        use DetBitOps, only: SQOperator 
        implicit none
        integer, intent(in) :: nSize1,nSize0,pertsite
        complex(dp), intent(in) :: Psi_1(nSize1)
        complex(dp), intent(out) :: Psi(nSize0)
        integer :: pertsitealpha,ilut,j,i
        logical :: tParity
        character(len=*), parameter :: t_r='ApplyCre_FirstOrder_EC'

        Psi(:) = dcmplx(0.0_dp,0.0_dp)
        if(nSize0.ne.(nFCIDet+nNp1FCIDet+nNm1bFCIDet)) then
            call stop_all(t_r,'Resultant wavefunction not expressed in expected space')
        endif

        pertsitealpha = 2*pertsite-1

        !Assume that the first nNm1bFCIDet entries of Psi_1 correspond to that space
        do i = 1,nNm1bFCIDet
            if(.not.btest(Nm1bBitList(i),pertsitealpha-1)) then
                !pertsitealpha empty. We can apply a creation operator
                ilut = Nm1bBitList(i)
                call SQOperator(ilut,pertsitealpha,tParity,.false.)

                !Now we need to find which det in FCIBitList this corresponds to
                do j = 1,nFCIDet
                    if(FCIBitList(j).eq.ilut) then
                        if(tParity) then
                            Psi(j) = -Psi_1(i)
                        else
                            Psi(j) = Psi_1(i)
                        endif
                        exit
                    endif
                enddo
                if(j.gt.nFCIDet) call stop_all(t_r,'Error in finding corresponding determinant')
            endif
        enddo
        
        if(.not.tLR_ReOptGS) then
            !There is no weight on the rest of the space in the GS, so we don't need to worry about projecting onto it.
            return
        endif

        do i = 1,nFCIDet
            if(.not.btest(FCIBitList(i),pertsitealpha-1)) then
                ilut = FCIBitList(i)
                call SQOperator(ilut,pertsitealpha,tParity,.false.)

                do j = 1,nNp1FCIDet
                    if(Np1BitList(j).eq.ilut) then
                        if(tParity) then
                            Psi(nFCIDet+j) = -Psi_1(nNm1bFCIDet+i)
                        else
                            Psi(nFCIDet+j) = Psi_1(nNm1bFCIDet+i)
                        endif
                        exit
                    endif
                enddo
                if(j.gt.nNp1FCIDet) call stop_all(t_r,'Error in finding corresponding determinant 2')
            endif
        enddo

    end subroutine ApplyCre_FirstOrder_EC


    !This will apply an annihilation operator to a wavefunction (Psi_1) in the space of the first-order interacting
    !N+1 particle space
    !This space is given by
    ! nNp1FCIDet space (+) nFCIDet space with contracted creation operator applied to core
    ! The first space couples to the uncontracted N electron space of the GS, while if the GS has been reoptimized, then
    ! we also have coupling between the second space with the N-1 (third) space of the GS wavefunction
    subroutine ApplyAnn_FirstOrder_EC(Psi_1,nSize1,Psi,nSize0,pertsite)
        use DetToolsData
        use DetBitOps, only: SQOperator 
        implicit none
        integer, intent(in) :: nSize1,nSize0,pertsite
        complex(dp), intent(in) :: Psi_1(nSize1)
        complex(dp), intent(out) :: Psi(nSize0)
        integer :: pertsitealpha,ilut,j,i
        logical :: tParity
        character(len=*), parameter :: t_r='ApplyAnn_FirstOrder_EC'

        Psi(:) = dcmplx(0.0_dp,0.0_dp)

        if(nSize0.ne.(nFCIDet+nNp1FCIDet+nNm1bFCIDet)) then
            call stop_all(t_r,'Resultant wavefunction not expressed in expected space')
        endif

        pertsitealpha = 2*pertsite-1

        do i = 1,nNp1FCIDet
            if(btest(Np1BitList(i),pertsitealpha-1)) then
                !pertsitealpha occupied. We can apply an annihilation operator
                ilut = Np1BitList(i)
                call SQOperator(ilut,pertsitealpha,tParity,.true.)

                !Now need to find which det in FCIBitList this corresponds to
                do j = 1,nFCIDet
                    if(FCIBitList(j).eq.ilut) then
                        !Found corresponding determinant
                        if(tParity) then
                            Psi(j) = -Psi_1(i)
                        else
                            Psi(j) = Psi_1(i)
                        endif
                        exit
                    endif
                enddo
                if(j.gt.nFCIDet) call stop_all(t_r,'Error in finding corresponding determinant')
            endif

        enddo

        if(.not.tLR_ReOptGS) then
            !There is no weight on the rest of the space in the GS, so we don't need to worry about projecting onto it.
            return
        endif

        !Now apply the annihilation operator to the contracted part of the wavefunction
        do i = 1,nFCIDet
            if(btest(FCIBitList(i),pertsitealpha-1)) then
                ilut = FCIBitList(i)
                call SQOperator(ilut,pertsitealpha,tParity,.true.)
                do j = 1,nNm1bFCIDet
                    if(Nm1bBitList(j).eq.ilut) then
                        if(tParity) then
                            Psi(nFCIDet+nNp1FCIDet+j) = -Psi_1(nNp1FCIDet+i)
                        else
                            Psi(nFCIDet+nNp1FCIDet+j) = Psi_1(nNp1FCIDet+i)
                        endif
                        exit
                    endif
                enddo
                if(j.gt.nNm1bFCIDet) call stop_all(t_r,'Error in finding corresponding determinant 2')
            endif
        enddo

    end subroutine ApplyAnn_FirstOrder_EC
            
    !It is only trying to create/destroy the alpha orbital applied to the ground state wavefunction
    !Will return the ground state, with a particle created or destroyed in pertsite
    subroutine ApplySP_PertGS_EC(Psi_0,SizeGS,V0_Cre,V0_Ann,nSizeLR,pertsite,tSwapExcits)
        use DetToolsData
        use DetBitOps, only: SQOperator 
        implicit none
        integer, intent(in) :: SizeGS,nSizeLR,pertsite
        complex(dp), intent(in) :: Psi_0(SizeGS)
        complex(dp), intent(out) :: V0_Cre(nSizeLR),V0_Ann(nSizeLR)
        logical, intent(in), optional :: tSwapExcits
        integer :: ilut,pertsitespin,i,j,GSInd
        logical :: tParity,tSwapExcits_
        character(len=*), parameter :: t_r='ApplySP_PertGS_EC'

        !if tSwapExcits_, then the perturbation is applied to the beta orbital
        !Otherwise, it is applied to the alpha orbital
        if(present(tSwapExcits)) then
            tSwapExcits_ = tSwapExcits
        else
            tSwapExcits_ = .false.
        endif

        V0_Ann(:) = zzero
        V0_Cre(:) = zzero 

        if(SizeGS.ne.(nFCIDet+nNp1FCIDet+nNm1bFCIDet)) then
            call stop_all(t_r,'Zeroth order wavefunction not expected size')
        endif
        if(nSizeLR.ne.(nFCIDet+nNp1FCIDet)) then
            call stop_all(t_r,'Size of linear system not expected')
        endif
        if(nNp1FCIDet.ne.nNm1bFCIDet) then
            call stop_all(t_r,"V0's are different sizes for particle/hole creation - not half filled active space")
        endif

        if(tSwapExcits_) then
            !apply to beta spin 
            pertsitespin = 2*pertsite
        else
            !Apply to alpha spin
            pertsitespin = 2*pertsite-1
        endif

        !This is going to be done in a very slow N^2 loop, rather than 
        do i=1,nFCIDet
            !We want to create a particle
            if(.not.btest(FCIBitList(i),pertsitespin-1)) then
                !spin-orbital is not occupied, we can create in it

                ilut = FCIBitList(i)
                !Create it, and check parity
                call SQOperator(ilut,pertsitespin,tParity,.false.)

                !This should really be binary searched here for efficiency
                if(tSwapExcits_) then
                    do j = 1,nNp1bFCIDet
                        if(Np1bBitList(j).eq.ilut) then
                            !Found corresponding determinant
                            if(tParity) then
                                V0_Cre(j) = -Psi_0(i)
                            else
                                V0_Cre(j) = Psi_0(i)
                            endif
                            exit
                        endif
                    enddo
                    if(j.eq.(nNp1bFCIDet+1)) call stop_all(t_r,'Could not find corresponding determinant')
                else
                    do j = 1,nNp1FCIDet
                        if(Np1BitList(j).eq.ilut) then
                            !Found corresponding determinant
                            if(tParity) then
                                V0_Cre(j) = -Psi_0(i)
                            else
                                V0_Cre(j) = Psi_0(i)
                            endif
                            exit
                        endif
                    enddo
                    if(j.eq.(nNp1FCIDet+1)) call stop_all(t_r,'Could not find corresponding determinant')
                endif
            else
                !Spin-orbital occupied, we can destroy it
                ilut = FCIBitList(i)
                !Annihilate it, and check parity
                call SQOperator(ilut,pertsitespin,tParity,.true.)

                if(tSwapExcits_) then
                    !This should really be binary searched here for efficiency
                    do j = 1,nNm1FCIDet
                        if(Nm1BitList(j).eq.ilut) then
                            !Found corresponding determinant
                            if(tParity) then
                                V0_Ann(j) = -Psi_0(i)
                            else
                                V0_Ann(j) = Psi_0(i)
                            endif
                            exit
                        endif
                    enddo
                    if(j.eq.(nNm1FCIDet+1)) call stop_all(t_r,'Could not find corresponding determinant')
                else
                    !This should really be binary searched here for efficiency
                    do j = 1,nNm1bFCIDet
                        if(Nm1bBitList(j).eq.ilut) then
                            !Found corresponding determinant
                            if(tParity) then
                                V0_Ann(j) = -Psi_0(i)
                            else
                                V0_Ann(j) = Psi_0(i)
                            endif
                            exit
                        endif
                    enddo
                    if(j.eq.(nNm1bFCIDet+1)) call stop_all(t_r,'Could not find corresponding determinant')
                endif
            endif
        enddo

        if(.not.tLR_ReoptGS) then
            !We do not have GS in other parts of the space
            return
        endif

        if(tSwapExcits_) then
            !This is part of the GS which will contribute to the particle removal (V0_Ann) wavefunction
            do i = 1,nNp1bFCIDet
                if(btest(Np1bBitList(i),pertsitespin-1)) then
                    ilut = Np1bBitList(i)
                    call SQOperator(ilut,pertsitespin,tParity,.true.)
                    GSInd = nFCIDet + i
                    do j = 1,nFCIDet
                        if(FCIBitList(j).eq.ilut) then
                            if(tParity) then
                                V0_Ann(nNm1FCIDet+j) = -Psi_0(GSInd)
                            else
                                V0_Ann(nNm1FCIDet+j) = Psi_0(GSInd)
                            endif
                            exit
                        endif
                    enddo
                    if(j.gt.nFCIDet) call stop_all(t_r,'Could not find corresponding determinant 3')
                endif
            enddo

            !This is the n-1 active space which will be operated on my the particle creation operator
            do i = 1,nNm1FCIDet
                if(.not.btest(Nm1BitList(i),pertsitespin-1)) then
                    ilut = Nm1BitList(i)
                    call SQOperator(ilut,pertsitespin,tParity,.false.)
                    GSInd = nFCIDet + nNp1bFCIDet + i 
                    do j = 1,nFCIDet
                        if(FCIBitList(j).eq.ilut) then
                            if(tParity) then
                                V0_Cre(nNp1bFCIDet+j) = -Psi_0(GSInd)
                            else
                                V0_Cre(nNp1bFCIDet+j) = Psi_0(GSInd)
                            endif
                            exit
                        endif
                    enddo
                    if(j.gt.nFCIDet) call stop_all(t_r,'Could not find corresponding determinant 2')
                endif
            enddo
        else
            !This is part of the GS which will contribute to the particle removal (V0_Ann) wavefunction
            do i = 1,nNp1FCIDet
                if(btest(Np1BitList(i),pertsitespin-1)) then
                    ilut = Np1BitList(i)
                    call SQOperator(ilut,pertsitespin,tParity,.true.)
                    GSInd = nFCIDet + i
                    do j = 1,nFCIDet
                        if(FCIBitList(j).eq.ilut) then
                            if(tParity) then
                                V0_Ann(nNm1bFCIDet+j) = -Psi_0(GSInd)
                            else
                                V0_Ann(nNm1bFCIDet+j) = Psi_0(GSInd)
                            endif
                            exit
                        endif
                    enddo
                    if(j.gt.nFCIDet) call stop_all(t_r,'Could not find corresponding determinant 3')
                endif
            enddo

            !This is the n-1 active space which will be operated on my the particle creation operator
            do i = 1,nNm1bFCIDet
                if(.not.btest(Nm1bBitList(i),pertsitespin-1)) then
                    ilut = Nm1bBitList(i)
                    call SQOperator(ilut,pertsitespin,tParity,.false.)
                    GSInd = nFCIDet + nNp1FCIDet + i 
                    do j = 1,nFCIDet
                        if(FCIBitList(j).eq.ilut) then
                            if(tParity) then
                                V0_Cre(nNp1FCIDet+j) = -Psi_0(GSInd)
                            else
                                V0_Cre(nNp1FCIDet+j) = Psi_0(GSInd)
                            endif
                            exit
                        endif
                    enddo
                    if(j.gt.nFCIDet) call stop_all(t_r,'Could not find corresponding determinant 2')
                endif
            enddo
        endif

    end subroutine ApplySP_PertGS_EC


    !Apply the density perturbation in the externally contracted space which is a rotation
    !into the non-null and orthogonal basis of the overlap matrix
    subroutine ApplyDensityPert_EC_Orthog(GS,V0,nSpan,nFull,Transform)
        use DetToolsData
        implicit none
        integer , intent(in) :: nSpan,nFull
        complex(dp), intent(in) :: GS(nSpan),Transform(nFull,nSpan)
        complex(dp), intent(out) :: V0(nSpan)
        complex(dp), allocatable :: V_Full(:,:),temppsi(:),DiagV(:),tempc(:,:)
        complex(dp), allocatable :: V_Span(:,:)
        integer :: i
        !character(len=*), parameter :: t_r='ApplyDensityPert_EC_Orthog'

        allocate(V_Full(nFull,nFull))
        V_Full(:,:) = dcmplx(0.0_dp,0.0_dp)

        !V in the full basis is diagonal.
        allocate(temppsi(nFull))
        temppsi(:) = dcmplx(1.0_dp,0.0_dp)
        allocate(DiagV(nFull))

        call ApplyDensityPert_EC(temppsi,DiagV,nFull) 
        do i = 1,nFull
            V_Full(i,i) = DiagV(i)
        enddo
        deallocate(DiagV,temppsi)

        !Now, rotate into new basis
        allocate(tempc(nSpan,nFull))
        allocate(V_Span(nSpan,nSpan))
        call ZGEMM('C','N',nSpan,nFull,nFull,dcmplx(1.0_dp,0.0_dp),Transform,nFull,    &
            V_Full,nFull,dcmplx(0.0_dp,0.0_dp),tempc,nSpan)
        call ZGEMM('N','N',nSpan,nSpan,nFull,dcmplx(1.0_dp,0.0_dp),tempc,nSpan,Transform,nFull,    &
            dcmplx(0.0_dp,0.0_dp),V_Span,nSpan)

        !Now apply V in the smaller orthogonal space
        call ZGEMV('N',nSpan,nSpan,dcmplx(1.0_dp,0.0_dp),V_Span,nSpan,GS,1,dcmplx(0.0_dp,0.0_dp),V0,1)
!        call ZGEMM('N','N',nSpan,1,nSpan,dcmplx(1.0_dp,0.0_dp),V_Span,nSpan,GS,nSpan,  &
!            dcmplx(0.0_dp,0.0_dp),V0,nSpan)

        deallocate(V_Span,tempc,V_Full)

    end subroutine ApplyDensityPert_EC_Orthog

    !Apply the density perturbation in the Externally contracted space
    subroutine ApplyDensityPert_EC(GS,V0,nSize)
        use DetToolsData
        implicit none
        integer , intent(in) :: nSize
        complex(dp), intent(in) :: GS(nSize)
        complex(dp), intent(out) :: V0(nSize)
        integer :: pertsitealpha,pertsitebeta,i,j,ind,pertsite
        character(len=*), parameter :: t_r='ApplyDensityPert_EC'

        !Initially, assume that the perturbation only acts at site 1
        pertsite = 1

        V0(:) = dcmplx(0.0_dp,0.0_dp)
        pertsitealpha = 2*pertsite-1
        pertsitebeta = 2*pertsite

        !First, look through uncontracted space
        do i = 1,nFCIDet
            if(btest(FCIBitList(i),pertsitealpha-1)) then
                !The perturbation site is occupied
                V0(i) = V0(i) + GS(i)
            endif
            if(btest(FCIBitList(i),pertsitebeta-1)) then
                V0(i) = V0(i) + GS(i)
            endif
        enddo
        !Now, CV space
        do i = CVIndex,AVIndex-1
            if(btest(FCIBitList(i-nFCIDet),pertsitealpha-1)) then
                !The perturbation site is occupied
                V0(i) = V0(i) + GS(i)
            endif
            if(btest(FCIBitList(i-nFCIDet),pertsitebeta-1)) then
                V0(i) = V0(i) + GS(i)
            endif
        enddo
        !Now, AV space - first, the alpha annihilated space
        do i = 1,nImp*2
            do j = 1,nNm1bFCIDet
                ind = AVIndex + (i-1)*(nNm1FCIDet+nNm1bFCIDet)+(j-1)
                if((ind.ge.CAIndex).or.(ind.lt.AVIndex)) then
                    call stop_all(t_r,'Indexing error here')
                endif
                if(btest(Nm1bBitList(j),pertsitealpha-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
                if(btest(Nm1bBitList(j),pertsitebeta-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
            enddo
        enddo
        !Now, AV space - beta annihilated space
        do i = 1,nImp*2
            do j = 1,nNm1FCIDet
                ind = AVIndex + (i-1)*(nNm1FCIDet+nNm1bFCIDet)+nNm1bFCIDet+(j-1)
                if((ind.ge.CAIndex).or.(ind.lt.AVIndex)) then
                    call stop_all(t_r,'Indexing error here')
                endif
                if(btest(Nm1BitList(j),pertsitealpha-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
                if(btest(Nm1BitList(j),pertsitebeta-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
            enddo
        enddo
        !Now, CA space - first, the alpha created space
        do i = 1,nImp*2
            do j = 1,nNp1FCIDet
                ind = CAIndex + (i-1)*(nNp1FCIDet+nNp1bFCIDet)+(j-1)
                if((ind.gt.nSize).or.(ind.lt.CAIndex)) then
                    call stop_all(t_r,'Indexing error here')
                endif
                if(btest(Np1BitList(j),pertsitealpha-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
                if(btest(Np1BitList(j),pertsitebeta-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
            enddo
        enddo
        !Now, CA space - beta created space
        do i = 1,nImp*2
            do j = 1,nNp1bFCIDet
                ind = CAIndex + (i-1)*(nNp1FCIDet+nNp1bFCIDet)+nNp1FCIDet+(j-1)
                if((ind.gt.nSize).or.(ind.lt.CAIndex)) then
                    call stop_all(t_r,'Indexing error here')
                endif
                if(btest(Np1bBitList(j),pertsitealpha-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
                if(btest(Np1bBitList(j),pertsitebeta-1)) then
                    V0(ind) = V0(ind) + GS(ind)
                endif
            enddo
        enddo


    end subroutine ApplyDensityPert_EC
                
    subroutine CreateCmprsCouplingMatrices(NMax_Coup_Create,NMax_Coup_Ann,Coup_Create,Coup_Ann, &
            Coup_Create_inds,Coup_Ann_inds,Coup_Create_cum,Coup_Ann_cum,Coup_Create_T,Coup_Ann_T,   &
            Coup_Create_inds_T,Coup_Ann_inds_T,Coup_Create_cum_T,Coup_Ann_cum_T,tSwapSpins)
        use DetToolsData
        use DetTools, only: tospat
        use DetBitOps, only: DecodeBitDet,SQOperator,CountBits
        implicit none
        integer, intent(in) :: NMax_Coup_Create,NMax_Coup_Ann
        integer, intent(out) :: Coup_Create(2,Nmax_Coup_Create),Coup_Ann(2,Nmax_Coup_Ann)
        integer, intent(out) :: Coup_Create_inds(Nmax_Coup_Create),Coup_Ann_inds(Nmax_Coup_Ann)
        integer, intent(out) :: Coup_Create_cum(nFCIDet+1),Coup_Ann_cum(nFCIDet+1)
        integer, intent(out) :: Coup_Create_T(2,Nmax_Coup_Create),Coup_Ann_T(2,Nmax_Coup_Ann)
        integer, intent(out) :: Coup_Create_inds_T(Nmax_Coup_Create),Coup_Ann_inds_T(Nmax_Coup_Ann)
        integer, intent(out) :: Coup_Create_cum_T(nNm1FCIDet+1),Coup_Ann_cum_T(nNp1FCIDet+1)
        logical, intent(in), optional :: tSwapSpins
        logical :: tSwapSpins_,tParity
        integer :: i,j,k,i_Ann,nOrbs,DiffOrb,orbdum(1),gam,CoreEnd,tempK
        character(len=*), parameter :: t_r='CreateCmprsCouplingMatrices'
        
        CoreEnd = nOcc-nImp

        if(nNm1FCIDet.ne.nNm1bFCIDet) call stop_all(t_r,'Cannot do different spins like this if spaces different size')

        !if tSwapSpins_ = T, then want to compute the coupling matrices for the creation and annihilation of an *alpha* spin orbital.
        !Otherwise, a beta spin orbital
        if(present(tSwapSpins)) then
            tSwapSpins_ = tSwapSpins
        else
            tSwapSpins_ = .false.
        endif

        Coup_Create(:,:) = 0
        Coup_Ann(:,:) = 0
        Coup_Create_inds(:) = 0
        Coup_Ann_inds(:) = 0
        Coup_Create_cum(:) = 0
        Coup_Ann_cum(:) = 0
        Coup_Create_cum(1) = 1
        Coup_Ann_cum(1) = 1
        
        Coup_Create_T(:,:) = 0
        Coup_Ann_T(:,:) = 0
        Coup_Create_inds_T(:) = 0
        Coup_Ann_inds_T(:) = 0
        Coup_Create_cum_T(:) = 0
        Coup_Ann_cum_T(:) = 0
        Coup_Create_cum_T(1) = 1
        Coup_Ann_cum_T(1) = 1


        i = 0
        i_Ann = 0
        do J = 1,nFCIDet
            do K = 1,nNm1FCIDet
                if(tSwapSpins_) then
                    !Start with the creation of an beta orbital
                    DiffOrb = ieor(Nm1BitList(K),FCIBitList(J))
                else
                    !Start with the creation of an alpha orbital
                    DiffOrb = ieor(Nm1bBitList(K),FCIBitList(J))
                endif
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
                if((mod(gam,2).ne.1).and.(.not.tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be an alpha spin orbital')
                elseif((mod(gam,2).ne.0).and.(tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be a beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                if(tSwapSpins_) then
                    tempK = Nm1BitList(K)
                else
                    tempK = Nm1bBitList(K)
                endif
                call SQOperator(tempK,gam,tParity,.false.)

                i = i + 1
                if(i.gt.Nmax_Coup_Create) call stop_all(t_r,'Compressed array size too small')
                Coup_Create(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create(2,i) = -1
                else
                    Coup_Create(2,i) = 1
                endif
                Coup_Create_inds(i) = K
            enddo
            !Now, annihilation of an alpha orbital
            do K = 1,nNp1FCIDet
                if(tSwapSpins_) then
                    DiffOrb = ieor(Np1bBitList(K),FCIBitList(J))
                else
                    DiffOrb = ieor(Np1BitList(K),FCIBitList(J))
                endif
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
                if(mod(gam,2).ne.1.and.(.not.tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be an alpha spin orbital')
                elseif(mod(gam,2).ne.0.and.tSwapSpins_) then
                    call stop_all(t_r,'differing orbital should be a beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                if(tSwapSpins_) then
                    tempK = Np1bBitList(K)
                else
                    tempK = Np1BitList(K)
                endif
                call SQOperator(tempK,gam,tParity,.true.)

                i_Ann = i_Ann + 1
                if(i_Ann.gt.Nmax_Coup_Ann) call stop_all(t_r,'Compressed array size too small')
                Coup_Ann(1,i_Ann) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann(2,i_Ann) = -1
                else
                    Coup_Ann(2,i_Ann) = 1
                endif
                Coup_Ann_inds(i_Ann) = K
            enddo
            Coup_Create_cum(J+1) = i+1
            Coup_Ann_cum(J+1) = i_Ann + 1
        enddo

        i = 0
        do J = 1,nNm1FCIDet
            do K = 1,nFCIDet
                if(tSwapSpins_) then
                    !Start with the creation of a beta orbital
                    DiffOrb = ieor(Nm1BitList(J),FCIBitList(K))
                else
                    !Start with the creation of an alpha orbital
                    DiffOrb = ieor(Nm1bBitList(J),FCIBitList(K))
                endif
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
                if(mod(gam,2).ne.1.and.(.not.tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be an alpha spin orbital')
                elseif(mod(gam,2).ne.0.and.(tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be a beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = FCIBitList(K)
                call SQOperator(tempK,gam,tParity,.true.)

                i = i + 1
                if(i.gt.Nmax_Coup_Create) call stop_all(t_r,'Compressed array size too small')
                Coup_Create_T(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Create_T(2,i) = -1
                else
                    Coup_Create_T(2,i) = 1
                endif
                Coup_Create_inds_T(i) = K
            enddo
            Coup_Create_cum_T(J+1) = i+1
        enddo

        i = 0
        do J = 1,nNp1FCIDet
            do K = 1,nFCIDet
                if(tSwapSpins_) then
                    !Annihilation of a beta orbital
                    DiffOrb = ieor(Np1bBitList(J),FCIBitList(K))
                else
                    !Annihilation of an alpha orbital
                    DiffOrb = ieor(Np1BitList(J),FCIBitList(K))
                endif
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
                if(mod(gam,2).ne.1.and.(.not.tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be an alpha spin orbital')
                elseif(mod(gam,2).ne.0.and.(tSwapSpins_)) then
                    call stop_all(t_r,'differing orbital should be a beta spin orbital')
                endif
                !Now find out the parity change when applying this creation operator to the original determinant
                tempK = FCIBitList(K)
                call SQOperator(tempK,gam,tParity,.false.)

                i = i + 1
                if(i.gt.Nmax_Coup_Ann) call stop_all(t_r,'Compressed array size too small')
                Coup_Ann_T(1,i) = tospat(gam)+CoreEnd
                if(tParity) then
                    Coup_Ann_T(2,i) = -1
                else
                    Coup_Ann_T(2,i) = 1
                endif
                Coup_Ann_inds_T(i) = K
            enddo
            Coup_Ann_cum_T(J+1) = i+1
        enddo
        
        i = 6*Nmax_Coup_Create + 6*Nmax_Coup_Ann + 2*nFCIDet + nNm1bFCIDet + nNp1FCIDet
        write(6,"(A,F12.5,A)") "Memory required for coupling coefficient matrices: ",real(i,dp)*RealtoMb, " Mb"
        i = 2*nFCIDet*nNm1bFCIDet + 2*nFCIDet*nNp1FCIDet
        write(6,"(A,F12.5,A)") "Uncompressed memory required for coupling coefficient matrices: ",real(i,dp)*RealtoMb, " Mb"
    
    end subroutine CreateCmprsCouplingMatrices

    !Contract the basis of single excitations, by summing together all the uncontracted parts with the non-interacting LR coefficients
    !This results in requiring the solution of a system which *only* scales with the size of the active hilbert space, not the lattice
    subroutine NonIntContracted_TDA_MCLR()
        use utils, only: get_free_unit
        use DetBitOps, only: DecodeBitDet,SQOperator
        use DetToolsData
        implicit none
!        integer :: nLinearSystem,ierr,i,j,a,b,info,lwork,n,iunit
!        integer :: i_spat,a_spat,gtid,alpha,beta,gam,dj,ExcitInd
!        integer :: AS_Spin_start,AS_Spin_end,alphap,Det,ExcitInd2,p,q,r,s,umatind
!        integer :: OrbPairs,UMatSize,alpha_AS,alphap_AS,ex(2,2),TmpDet,TmpnI(elec)
!        integer :: pertsitealpha,pertsitebeta,nSize,iunit2
!        integer, allocatable :: Pivots(:)
!        logical :: tSign
!        complex(dp) :: ni_lr
!        real(dp) :: CoreCoupling,CoreVirtualNorm,Omega,Res1,Res2,EDiff,ResponseFn 
!        real(dp) :: testham,testnorm,tmp,ParityFac
!        real(dp), allocatable :: LinearSystem(:,:),temp(:,:),Work(:),W(:),Residues(:)
!        real(dp), allocatable :: RDM1(:,:),RDM2(:,:),temp_vec(:),Projector(:,:)
!        real(dp), allocatable :: Trans1RDM_bra(:,:),Trans1RDM_ket(:,:),Overlap(:,:)
!        real(dp), allocatable :: Nm1AlphaVec(:),Np1AlphaVec(:),VGS(:)
!        real(dp), allocatable :: Nm1AlphapVec(:),Np1AlphapVec(:)
!        real(dp), allocatable :: CoreActiveNorm(:),ActiveVirtualNorm(:)
!        real(dp), allocatable :: Nm1AlphaRDM(:,:),Np1AlphaRDM(:,:)
!        real(dp), allocatable :: Nm1Alpha2RDM(:,:,:,:),Np1Alpha2RDM(:,:,:,:)
        character(len=*), parameter :: t_r='NonIntContracted_TDA_MCLR'
!        logical, parameter :: tNonIntTest = .false.
!        logical, parameter :: tDiagonalize = .false.
!        logical, parameter :: tResiduesFromRDM = .false. 
!        logical, parameter :: tDebug = .false. 

        write(6,*) "Calculating non-interacting IC MR-TDA LR system..."
        call stop_all(t_r,"This routine is currently broken. Get working for complex frequencies before debugging...")

!        if(.not.tConstructFullSchmidtBasis) call stop_all(t_r,'To solve LR, must construct full schmidt basis')
!        if(.not.tCompleteDiag) call stop_all(t_r,'To solve LR, must perform complete diag')
!        !umat and tmat for the active space
!        OrbPairs = (EmbSize*(EmbSize+1))/2
!        UMatSize = (OrbPairs*(OrbPairs+1))/2
!        if(allocated(UMat)) deallocate(UMat)
!        allocate(UMat(UMatSize))
!        UMat(:) = 0.0_dp
!        do i=1,nImp
!            umat(umatind(i,i,i,i)) = U
!        enddo
!        
!        !Enumerate excitations for fully coupled space
!        call GenDets(Elec,EmbSize,.true.,.true.,.false.)
!        write(6,*) "Number of determinants in {N,N+1,N-1} FCI space: ",ECoupledSpace
!        !Calculate the size of the hamiltonian matrix
!        !This is simply the normal active space size, plus 1 fully contracted core-virtual excitation,
!        !plus 2*nImp fully contracted core-active excitation, and 2*nImp fully contracted active-virtual excitations
!        nLinearSystem = nFCIDet+1+8*nImp 
!        !Test that we reduce to the non-interacting limit
!        if(tNonIntTest) then
!            nLinearSystem = 1
!            nFCIDet = 0
!            nImp = 0
!        endif
!        
!        write(6,"(A,F14.6,A)") "Memory required for the LR hessian: ",real((nLinearSystem**2)*8,dp)/1048576.0_dp," Mb"
!        
!        iunit = get_free_unit()
!        open(unit=iunit,file='IC-TDA_DDResponse',status='unknown')
!        write(iunit,"(A)") "# Frequency     DD_LinearResponse"
!        
!        iunit2 = get_free_unit()
!        open(unit=iunit2,file='IC-TDA_EValues',status='unknown')
!        write(iunit2,"(A)") "# Frequency     EValues..."
!
!        allocate(Trans1RDM_bra(EmbSize,EmbSize))
!        allocate(Trans1RDM_ket(EmbSize,EmbSize))
!
!        allocate(Nm1AlphaVec(nNm1FCIDet))
!        allocate(Np1AlphaVec(nNp1FCIDet))
!        allocate(Nm1AlphapVec(nNm1FCIDet))
!        allocate(Np1AlphapVec(nNp1FCIDet))
!
!        allocate(Nm1AlphaRDM(EmbSize,EmbSize))
!        allocate(Np1AlphaRDM(EmbSize,EmbSize))
!        allocate(Nm1Alpha2RDM(EmbSize,EmbSize,EmbSize,EmbSize))
!        allocate(Np1Alpha2RDM(EmbSize,EmbSize,EmbSize,EmbSize))
!
!        AS_Spin_start = ((nOcc-nImp+1)*2)-1     !The starting index of the active space in spin-orbital notation
!        AS_Spin_end = (nOcc+nImp)*2     !odd = alpha, even = beta
!
!        write(6,*) "Starting spin-orbital for active space: ",AS_Spin_start
!        write(6,*) "Final spin-orbital for active space: ",AS_Spin_end
!        
!        allocate(CoreActiveNorm(AS_Spin_start:AS_Spin_end))
!        allocate(ActiveVirtualNorm(AS_Spin_start:AS_Spin_end))
!            
!        !Allocate memory for hmailtonian in this system:
!        allocate(LinearSystem(nLinearSystem,nLinearSystem),stat=ierr)
!        allocate(Overlap(nLinearSystem,nLinearSystem),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,'Error allocating')
!
!        Omega = Start_Omega
!        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
!
!            write(6,*) "Calculating linear response for frequency: ",Omega
!
!            !First, find the non-interacting solution expressed in the schmidt basis
!            call FindSchmidtPert(tNonIntTest,Omega,ni_lr)
!
!!            call writematrix(SchmidtPert,'SchmidtPert',.true.)
!!            call writematrix(FockSchmidt,'FockSchmidt',.true.)
!
!            if(.not.tDiagonalize) then
!                allocate(VGS(nLinearSystem))
!                allocate(Pivots(nLinearSystem))
!            endif
!
!
!            LinearSystem(:,:) = 0.0_dp
!            Overlap(:,:) = 0.0_dp
!
!            !write(6,"(A)",advance='no') "Constructing hessian matrix..."
!            
!            !First, construct FCI space, in determinant basis
!            !This is the first block
!            do i=1,nFCIDet
!                LinearSystem(i,i) = Spectrum(i) 
!            enddo
!            !Now transform this block back into the determinant basis
!            allocate(temp(nFCIDet,nFCIDet))
!            if(.not.tNonIntTest) then
!                call DGEMM('N','N',nFCIDet,nFCIDet,nFCIDet,1.0_dp,FullHamil,nFCIDet,LinearSystem(1:nFCIDet,1:nFCIDet),  &
!                    nFCIDet,0.0_dp,temp,nFCIDet)
!                call DGEMM('N','T',nFCIDet,nFCIDet,nFCIDet,1.0_dp,temp,nFCIDet,FullHamil,nFCIDet,0.0_dp,    &
!                    LinearSystem(1:nFCIDet,1:nFCIDet),nFCIDet)
!            endif
!            deallocate(temp)
!
!            !Now, calculate the fully IC sum of all core-virtual excitations
!            !The diagonal hamiltonian matrix element is given by: G_ai G_bi F_ab - G_ai G_aj F_ij
!            !There is no coupling to the active space via the mean-field hamiltonian, since we cant have just a single active space index
!            !The contracted excitation in orthogonal to the FCI space uncontracted determinants (since cores are always orthogonal), but 
!            !the contracted function itself is not normalized. Find the normalization, and renomalize all its contributions.
!            CoreVirtualNorm = 0.0_dp
!            do i=1,nOcc-nImp
!                do a=nOcc+nImp+1,nSites
!                    CoreVirtualNorm = CoreVirtualNorm + SchmidtPert(i,a)*SchmidtPert(a,i)
!                enddo
!            enddo
!            CoreVirtualNorm = CoreVirtualNorm * 2.0_dp  !Spin integration 
!            write(6,*) "Fully contracted core-virtual excitations have a normalization of: ",Omega,CoreVirtualNorm
!            !Diagonal term for CV excitation
!            tmp = 0.0_dp
!            do i=1,nOcc-nImp
!                do j=1,nOcc-nImp
!                    do a=nOcc+nImp+1,nSites
!                        LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1) -     &
!                            SchmidtPert(a,i)*SchmidtPert(a,j)*FockSchmidt(i,j)*2.0_dp
!                    enddo
!                enddo
!            enddo
!!            write(6,"(A,2G25.17)") "Term 1 for CV: ",Omega,LinearSystem(nFCIDet+1,nFCIDet+1)
!            tmp = 0.0_dp
!            do i=1,nOcc-nImp
!                do a=nOcc+nImp+1,nSites
!                    do b=nOcc+nImp+1,nSites
!                        !LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1) +     &
!                        tmp = tmp +     &
!                            SchmidtPert(a,i)*SchmidtPert(b,i)*FockSchmidt(a,b)*2.0_dp
!                    enddo
!                enddo
!            enddo
!!            write(6,"(A,2G25.17)") "Term 2 for CV: ",Omega,tmp
!            LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1) + tmp
!            LinearSystem(nFCIDet+1,nFCIDet+1) = LinearSystem(nFCIDet+1,nFCIDet+1)/CoreVirtualNorm  
!            write(6,"(A,2G25.17)") "Diagonal hamiltonian contribution from fully contracted core-virtual function: ",   &
!                Omega,LinearSystem(nFCIDet+1,nFCIDet+1)
!            !We are not going to add on the active space energy, since we assume that we have offset the hamiltonian by the 
!            !zeroth order energy.
!
!            !Now, we need to define the coupling to the uncontracted determinant space.
!            CoreCoupling = 0.0_dp
!            do i=1,nOcc-nImp
!                do a=nOcc+nImp+1,nSites
!                    CoreCoupling = CoreCoupling + FockSchmidt(i,a)*SchmidtPert(a,i)
!                enddo
!            enddo
!            CoreCoupling = CoreCoupling * 2.0_dp    !For the other spin type.
!            !Now we need to include the contribution from <D_j|0>
!            do i=1,nFCIDet
!                LinearSystem(i,nFCIDet+1) = CoreCoupling*FullHamil(i,1)/sqrt(CoreVirtualNorm)
!                LinearSystem(nFCIDet+1,i) = LinearSystem(i,nFCIDet+1)   !Hermiticity
!            enddo
!
!            !*****************************************************************************************
!            !Now for the active-virtual semi-internal excitations
!            !*****************************************************************************************
!            !Precompute the normalization constants for each semi-internal excitation
!            ActiveVirtualNorm(:) = 0.0_dp
!            do alpha = AS_Spin_start,AS_Spin_end
!                do a=nOcc+nImp+1,nSites
!                    ActiveVirtualNorm(alpha) = ActiveVirtualNorm(alpha) +   &
!                        SchmidtPert(a,gtid(alpha))*SchmidtPert(gtid(alpha),a)*  &
!                        HL_1RDM(gtid(alpha)-nOcc+nImp,gtid(alpha)-nOcc+nImp)/2.0_dp !We only want one spin-type now
!                enddo
!            enddo
!!            call writevector(ActiveVirtualNorm,'AV_Norm')
!            !First, creating a particle in the virtual manifold, for each annihilation in the active space
!            ExcitInd = nFCIDet + 1
!            do alpha = AS_Spin_start,AS_Spin_end
!                ExcitInd = ExcitInd + 1
!
!                alpha_AS = alpha-2*(nOcc-nImp)  !The spin-orbital label of alpha, starting from 1
!
!                Nm1AlphaVec(:) = 0.0_dp 
!                !Create 1 and 2 body symmetric RDMs for the N-1 electron system where we have annihilated spin-orbital alpha
!                !from |0>
!                do i = 1,nFCIDet
!                    if(btest(FCIBitList(i),alpha_AS-1)) then
!!                    if(IsOcc(FCIBitList(i),alpha)) then
!                        !We can annihilate it
!
!                        Det = FCIBitList(i)
!
!                        !Delete orbital from Det
!                        call SQOperator(Det,alpha_AS,tSign,.true.)
!                        !Det = ibclr(Det,alpha_AS-1)
!
!                        !Find this in the Nm1 list
!                        do j = 1,nNm1FCIDet
!                            if(Nm1BitList(j).eq.Det) then
!                                !We have found the orbital. Store the amplitude at this point
!                                if(tSign) then
!                                    Nm1AlphaVec(j) = -FullHamil(i,1)
!                                else
!                                    Nm1AlphaVec(j) = FullHamil(i,1)
!                                endif
!                                exit
!                            endif
!                        enddo
!                        if(j.gt.nNm1FCIDet) call stop_all(t_r,'Can not find appropriate determinant in N-1 space')
!
!                    endif
!                enddo
!                !We now have the wavefunction a_alpha |0> expressed in the N-1 electron FCI basis
!                !Now create the 1 and 2 body density matrices from these amplitudes for the resulting N-1 electron wavefunction
!                call CalcNm1_1RDM(Nm1AlphaVec,Nm1AlphaVec,Nm1AlphaRDM)
!                call CalcNm1_2RDM(Nm1AlphaVec,Nm1AlphaVec,Nm1Alpha2RDM)
!
!!                !Check that these correspond correctly to the higher-ranked density matrix in the N electron wavefunction
!!                do i = 1,nImp*2
!!                    do j = 1,nImp*2
!!                        if(abs(Nm1AlphaRDM(i,j)-(HL_2RDM(gtid(alpha)-nOcc+nImp,gtid(alpha)-nOcc+nImp,i,j)/2.0_dp))  &
!!                            .gt.1.0e-7_dp) then
!!                            write(6,*) "alpha: ",gtid(alpha)
!!                            write(6,*) "AS Spatial indices: ",i,j
!!                            write(6,*) "Nm1AlphaRDM(i,j) : ",Nm1AlphaRDM(i,j)
!!                            write(6,*) "HL_2RDM(alpha,alpha,i,j) : ",HL_2RDM(gtid(alpha)-nOcc+nImp,gtid(alpha)-nOcc+nImp,i,j)
!!                            write(6,*) "HL_2RDM(alpha,alpha,i,j)/2 : ",HL_2RDM(gtid(alpha)-nOcc+nImp,gtid(alpha)-nOcc+nImp,i,j)/2.0
!!                            do dj = 1,nFCIDet
!!                                write(6,*) FCIDetList(:,dj),FullHamil(dj,1)  
!!                            enddo
!!                            call writevector(Nm1AlphaVec,'Nm1AlphaVec')
!!                            call writematrix(Nm1AlphaRDM,'Nm1AlphaRDM',.true.)
!!                            call stop_all(t_r,'N-1 RDMs do not correspond to higher rank equivalent operators')
!!                        endif
!!                    enddo
!!                enddo
!                
!                !Now for the coupling to the determinant space (dj)
!                do dj = 1,nFCIDet
!                    
!                    do beta = 1,4*nImp  !Loop over AS spin-orbitals to annihilate
!                        if(mod(beta,2).ne.mod(alpha,2)) cycle
!
!                        if(btest(FCIBitList(dj),beta-1)) then
!                            !beta is occupied in dj. Annihilate it
!                            Det = FCIBitList(dj)
!                            call SQOperator(Det,beta,tSign,.true.)
!!                            Det = ibclr(Det,beta-1)     
!
!                            do i = 1,nNm1FCIDet
!                                if(Nm1BitList(i).eq.Det) then
!                                    !We have found the corresponding determinant
!                                    exit
!                                endif
!                            enddo
!                            if(i.gt.nNm1FCIDet) call stop_all(t_r,'Can not find appropriate determinant in N-1 space 2')
!                            if(abs(Nm1AlphaVec(i)).lt.1.0e-8_dp) cycle  !There is not alpha|0> component here
!
!                            if(tSign) then
!                                ParityFac = -1.0_dp
!                            else
!                                ParityFac = 1.0_dp
!                            endif
!
!!                            if(alpha_AS.ne.beta) then
!!                                !We need to find the permutation of this excitation, rather than just taking the coefficient.
!!                                !Obviously, this is a silly way of doing it...
!!                                TmpDet = Det
!!                                TmpDet = ibset(TmpDet,alpha_AS-1)  !Create alpha and decode to get the original determinant
!!                                call DecodeBitDet(TmpnI,elec,TmpDet)
!!                                ex(1,1) = 1
!!                                call getexcitation(FCIDetList(:,dj),TmpnI,elec,ex,tSign)
!!                                if(ex(1,1).ne.beta) then
!!                                    call stop_all(t_r,'error calculating parity')
!!                                elseif(ex(2,1).ne.alpha_AS) then
!!                                    call stop_all(t_r,'error calculating parity 2')
!!                                endif
!!                                if(tSign) then
!!                                    !Negative parity
!!                                    ParityFac = -1.0_dp
!!                                else
!!                                    ParityFac = 1.0_dp
!!                                endif
!!                            else
!!                                !We should have created and annihilated the same orbital
!!                                ParityFac = 1.0_dp
!!                            endif
!
!                            !Coupling matrix element is Nm1AlphaVec(i)
!                            do a = nOcc+nImp+1,nSites
!                                LinearSystem(dj,ExcitInd) = LinearSystem(dj,ExcitInd) + SchmidtPert(a,gtid(alpha))* &
!                                    FockSchmidt(gtid(beta)+nOcc-nImp,a)*Nm1AlphaVec(i)*ParityFac/sqrt(ActiveVirtualNorm(alpha))
!                            enddo
!                        endif
!                    enddo
!                    LinearSystem(ExcitInd,dj) = LinearSystem(dj,ExcitInd)
!                enddo
!
!                !Now, for the coupling to the strongly contracted excitation
!                !Minus sign because we always annihilate first with the semi-internal excitations: 
!                do i = 1,nOcc-nImp
!                    do a = nOcc+nImp+1,nSites
!                        do beta = nOcc-nImp+1,nOcc+nImp
!                            LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) - &
!                                SchmidtPert(a,i)*SchmidtPert(a,gtid(alpha))*FockSchmidt(i,beta)* &
!                                HL_1RDM(gtid(alpha_AS),beta-nOcc+nImp)
!                        enddo
!                    enddo
!                enddo
!                LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) /  &
!                    (2.0_dp*sqrt(CoreVirtualNorm*ActiveVirtualNorm(alpha)))
!                LinearSystem(nFCIDet+1,ExcitInd) = LinearSystem(ExcitInd,nFCIDet+1)
!
!                !Now for the coupling to the other semi-internal excitations
!                ExcitInd2 = nFCIDet + 1
!                do alphap = AS_Spin_start,AS_Spin_end
!                    ExcitInd2 = ExcitInd2 + 1
!                    !We shouldn't have coupling between excitations of the same spin
!                    if(mod(alphap,2).ne.mod(alpha,2)) cycle
!                    alphap_AS = alphap-2*(nOcc-nImp)
!
!                    Nm1AlphapVec(:) = 0.0_dp
!                    !Create 1 and 2 body transition RDMs for the N-1 electron system where 
!                    !we have annihilated spin-orbital alphap from |0>
!                    do i = 1,nFCIDet
!                        if(btest(FCIBitList(i),alphap_AS-1)) then
!                            !We can annihilate it
!                            Det = FCIBitList(i)
!                            call SQOperator(Det,alphap_AS,tSign,.true.)
!!                            Det = ibclr(Det,alphap_AS-1)
!
!                            do j = 1,nNm1FCIDet
!                                if(Nm1BitList(j).eq.Det) then
!                                    !We have found the orbital. Store the amplitude at this point
!                                    if(tSign) then
!                                        Nm1AlphapVec(j) = -FullHamil(i,1)
!                                    else
!                                        Nm1AlphapVec(j) = FullHamil(i,1)
!                                    endif
!                                    exit
!                                endif
!                            enddo
!                            if(j.gt.nNm1FCIDet) call stop_all(t_r,'Can not find appropritate determinant in N-1 space')
!                        endif
!                    enddo
!                    
!                    !We now have the wavefunction a_alpha |0> expressed in the N-1 electron FCI basis
!                    !Now create the 1 and 2 body density matrices from these amplitudes for the resulting N-1 electron wavefunction
!                    !TODO: I think we will need to do these for spin-orbitals, or properly spin-integrate
!                    call CalcNm1_1RDM(Nm1AlphaVec,Nm1AlphapVec,Nm1AlphaRDM)
!                    call CalcNm1_2RDM(Nm1AlphaVec,Nm1AlphapVec,Nm1Alpha2RDM)
!                    !Check that these correspond correctly to the higher-ranked density matrix in the N electron wavefunction
!!                    do i = 1,nImp*2
!!                        do j = 1,nImp*2
!!                            if(abs(Nm1AlphaRDM(i,j)-(HL_2RDM(gtid(alpha_AS),gtid(alphap_AS),i,j)/2.0)).gt.1.0e-7_dp) then
!!                                write(6,*) "ERROR: ",abs(Nm1AlphaRDM(i,j)-(HL_2RDM(gtid(alpha_AS),gtid(alphap_AS),i,j)/2.0))
!!                                write(6,*) "alpha: ",gtid(alpha_AS),alpha_AS
!!                                write(6,*) "alphap: ",gtid(alphap_AS),alphap_AS
!!                                write(6,*) "AS Spatial indices: ",i,j
!!                                write(6,*) "Nm1AlphaRDM(i,j) : ",Nm1AlphaRDM(i,j)
!!                                write(6,*) "Nm1AlphaRDM(j,i) : ",Nm1AlphaRDM(j,i)
!!                                write(6,*) "HL_2RDM(alpha,alphap,i,j) : ",HL_2RDM(gtid(alpha_AS),gtid(alphap_AS),i,j)
!!                                write(6,*) "HL_2RDM(alpha,alphap,i,j)/2 : ",HL_2RDM(gtid(alpha_AS),gtid(alphap_AS),i,j)/2.0
!!                                do dj = 1,nFCIDet
!!                                    write(6,*) FCIDetList(:,dj),FullHamil(dj,1)  
!!                                enddo
!!                                call writevector(Nm1AlphaVec,'Nm1AlphaVec')
!!                                call writematrix(Nm1AlphaRDM,'Nm1AlphaRDM',.true.)
!!                                if(alpha_AS.eq.alphap_AS) then
!!                                    !They should at least be correct when alpha = alphap?
!!                                    call stop_all(t_r,'N-1 RDMs do not correspond to higher rank equivalent operators 2')
!!                                endif
!!                            endif
!!                        enddo
!!                    enddo
!
!!                    !This is a test - we don't actually need to do the diagonals!
!!                    tmp2 = LinearSystem(ExcitInd,ExcitInd)  !Save the diagonal element to ensure that it reduces to the same thing
!!                    LinearSystem(ExcitInd,ExcitInd) = 0.0_dp
!!
!                    do a = nOcc+nImp+1,nSites
!                        do b = nOcc+nImp+1,nSites
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) +     &
!                              SchmidtPert(a,gtid(alpha))*SchmidtPert(b,gtid(alphap))*FockSchmidt(b,a)*   &
!                              HL_1RDM(gtid(alpha_AS),gtid(alphap_AS))/2.0_dp
!                        enddo
!                    enddo
!                    do a = nOcc+nImp+1,nSites
!                        do i = 1,nOcc-nImp
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) +     &
!                              FockSchmidt(i,i)*SchmidtPert(gtid(alpha),a)*SchmidtPert(gtid(alphap),a)* &
!                              HL_1RDM(gtid(alpha_AS),gtid(alphap_AS))
!                        enddo
!                    enddo
!                    do a = nOcc+nImp+1,nSites
!                        do beta = 1,2*nImp
!                            do gam = 1,2*nImp
!                                LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) + &
!                                    SchmidtPert(gtid(alpha),a)*SchmidtPert(gtid(alphap),a)* &
!                                    FockSchmidt(beta+nOcc-nImp,gam+nOcc-nImp)*  &
!                                    HL_2RDM(gtid(alpha_AS),gtid(alphap_AS),beta,gam)/2.0_dp
!                                !TODO: Should be able to comment out the below without changing the results
!!                                    SchmidtPert(gtid(alpha),a)*SchmidtPert(gtid(alphap),a)*Nm1AlphaRDM(beta,gam)*   &
!!                                    FockSchmidt(beta+nOcc-nImp,gam+nOcc-nImp)
!                            enddo
!                        enddo
!                    enddo
!                    tmp = 0.0_dp
!                    do a = nOcc+nImp+1,nSites
!                        tmp = tmp + SchmidtPert(a,gtid(alpha))*SchmidtPert(a,gtid(alphap))
!                    enddo
!                    tmp = tmp/2.0_dp
!                    !Now for the two electron component:
!                    do p = 1,EmbSize            
!                        do q = 1,EmbSize             
!                            do r = 1,EmbSize               
!                                do s = 1,EmbSize               
!                                    LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) + &
!                                        tmp*Nm1Alpha2RDM(p,q,r,s)*umat(umatind(p,r,q,s)) 
!                                enddo
!                            enddo
!                        enddo
!                    enddo
!                    LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2)/    &
!                        (sqrt(ActiveVirtualNorm(alpha)*ActiveVirtualNorm(alphap)))
!                    LinearSystem(ExcitInd2,ExcitInd) = LinearSystem(ExcitInd,ExcitInd2)
!!                    if((ExcitInd2.eq.ExcitInd).and.(abs(tmp2-LinearSystem(ExcitInd,ExcitInd)).gt.1.0e-7_dp)) then
!!                        write(6,*) "ExcitInd: ",ExcitInd
!!                        write(6,*) "Original: ",tmp2
!!                        write(6,*) "New: ",LinearSystem(ExcitInd,ExcitInd)
!!                        call stop_all(t_r,'Error in consistent diagonal elements for semi-internal excitations')
!!                    endif
!
!                enddo !alphap
!
!            enddo  !Finish looping over active-virtual semi-internal excitations
!
!            !*****************************************************************************************
!            !Now for the core-active semi-internal excitations
!            !*****************************************************************************************
!            !Precompute the normalization constants for each semi-internal excitation
!            CoreActiveNorm(:) = 0.0_dp
!            do alpha = AS_Spin_start,AS_Spin_end
!                do i=1,nOcc-nImp 
!                    CoreActiveNorm(alpha) = CoreActiveNorm(alpha) +   &
!                        SchmidtPert(i,gtid(alpha))*SchmidtPert(gtid(alpha),i)
!                enddo
!                CoreActiveNorm(alpha) = CoreActiveNorm(alpha) *     &
!                    (1.0_dp - (HL_1RDM(gtid(alpha)-nOcc+nImp,gtid(alpha)-nOcc+nImp)/2.0_dp))
!            enddo
!!            call writevector(CoreActiveNorm,'CA_Norm')
!            !First, creating a hole in the occupied manifold, for each alpha
!            ExcitInd = nFCIDet + 1 + 4*nImp
!            do alpha = AS_Spin_start,AS_Spin_end
!                ExcitInd = ExcitInd + 1
!                if(LinearSystem(ExcitInd,ExcitInd).ne.0.0_dp) call stop_all(t_r,'Indexing error')
!
!                alpha_AS = alpha-2*(nOcc-nImp)  !The spin-orbital label of alpha, starting from 1
!
!                Np1AlphaVec(:) = 0.0_dp 
!                !Create 1 and 2 body symmetric RDMs for the N+1 electron system where we have created spin-orbital alpha
!                !from |0>
!                do i = 1,nFCIDet
!                    if(.not.btest(FCIBitList(i),alpha_AS-1)) then
!                        !It is unoccupied - We can create it
!
!                        Det = FCIBitList(i)
!                        call SQOperator(Det,alpha_AS,tSign,.false.)
!                        !Include orbital from Det
!!                        Det = ibset(Det,alpha_AS-1)
!
!                        !Find this in the Nm1 list
!                        do j = 1,nNp1FCIDet
!                            if(Np1BitList(j).eq.Det) then
!                                !We have found the orbital. Store the amplitude at this point
!                                if(tSign) then
!                                    Np1AlphaVec(j) = -FullHamil(i,1)
!                                else
!                                    Np1AlphaVec(j) = FullHamil(i,1)
!                                endif
!                                exit
!                            endif
!                        enddo
!                        if(j.gt.nNp1FCIDet) call stop_all(t_r,'Can not find appropriate determinant in N+1 space')
!
!                    endif
!                enddo
!                !We now have the wavefunction a_alpha |0> expressed in the N-1 electron FCI basis
!                !Now create the 1 and 2 body density matrices from these amplitudes for the resulting N-1 electron wavefunction
!                call CalcNm1_1RDM(Np1AlphaVec,Np1AlphaVec,Np1AlphaRDM)
!                call CalcNm1_2RDM(Np1AlphaVec,Np1AlphaVec,Np1Alpha2RDM)
!
!                !Now for the coupling to the determinant space (dj)
!                do dj = 1,nFCIDet
!                    
!                    do beta = 1,4*nImp  !Loop over AS spin-orbitals
!
!                        if(.not.btest(FCIBitList(dj),beta-1)) then
!                            !beta is unoccupied in dj. create it
!                            Det = FCIBitList(dj)
!                            call SQOperator(Det,beta,tSign,.false.)
!!                            Det = ibset(Det,beta-1)     
!
!                            do i = 1,nNp1FCIDet
!                                if(Np1BitList(i).eq.Det) then
!                                    !We have found the corresponding determinant
!                                    exit
!                                endif
!                            enddo
!                            if(i.gt.nNp1FCIDet) call stop_all(t_r,'Can not find appropriate determinant in N+1 space 2')
!                            if(abs(Np1AlphaVec(i)).lt.1.0e-8_dp) cycle  !There is not alpha^+|0> component here
!
!                            if(tSign) then
!                                ParityFac = -1.0_dp
!                            else
!                                ParityFac = 1.0_dp
!                            endif
!
!!                            if(alpha_AS.ne.beta) then
!!                                !Find permutation.
!!                                TmpDet = Det
!!                                TmpDet = ibclr(TmpDet,alpha_AS-1)
!!                                call DecodeBitDet(TmpnI,elec,TmpDet)
!!                                Ex(1,1) = 1
!!                                call getexcitation(FCIDetList(:,dj),TmpnI,elec,ex,tSign)
!!                                if(ex(1,1).ne.alpha_AS) then
!!                                    write(6,*) "FCIDetList: ",FCIDetList(:,dj)
!!                                    write(6,*) "Np1Det: ",Np1FCIDetList(:,i) 
!!                                    write(6,*) "TmpDet: ",TmpnI(:) 
!!                                    write(6,*) "Np1AlphaVec(i): ",Np1AlphaVec(i)
!!                                    write(6,*) "alpha, beta: ",alpha_AS,beta
!!                                    write(6,*) "ex(1,1): ",ex(1,1)
!!                                    call stop_all(t_r,'error calculating parity 3')
!!                                elseif(ex(2,1).ne.beta) then
!!                                    call stop_all(t_r,'error calculating parity 4')
!!                                endif
!!                                !The parities are switched, since we actually want beta alpha^+, rather than the canonical ordering
!!                                if(tSign) then
!!                                    !Positive parity
!!                                    ParityFac = 1.0_dp
!!                                else
!!                                    ParityFac = -1.0_dp
!!                                endif
!!                            else
!!                                !We should have created and annihilated the same orbital
!!                                ParityFac = 1.0_dp
!!                            endif
!
!                            !Coupling matrix element is Np1AlphaVec(i)
!                            do j = 1,nOcc-nImp
!                                LinearSystem(dj,ExcitInd) = LinearSystem(dj,ExcitInd) + SchmidtPert(j,gtid(alpha))* &
!                                    FockSchmidt(j,gtid(beta)+nOcc-nImp)*Np1AlphaVec(i)*ParityFac
!                            enddo
!                        endif
!                    enddo
!
!!                    do i = 1,nOcc-nImp
!!                        LinearSystem(dj,ExcitInd) = LinearSystem(dj,ExcitInd) + SchmidtPert(i,gtid(alpha))* &
!!                            FockSchmidt(i,gtid(alpha))*FullHamil(dj,1)
!!                    enddo
!                    LinearSystem(dj,ExcitInd) = LinearSystem(dj,ExcitInd)/sqrt(CoreActiveNorm(alpha))
!                    LinearSystem(ExcitInd,dj) = LinearSystem(dj,ExcitInd)
!                enddo
!
!                !Now, for the coupling to the strongly contracted excitation
!                do i = 1,nOcc-nImp
!                    do a = nOcc+nImp+1,nSites
!                        do beta = 1,nImp*2                           
!                            LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) - &
!                                SchmidtPert(i,a)*SchmidtPert(i,gtid(alpha))*FockSchmidt(a,beta+nOcc-nImp)* &
!                                HL_1RDM(gtid(alpha_AS),beta)/2.0_dp
!                        enddo
!                    enddo
!                enddo
!                do i = 1,nOcc-nImp
!                    do a = nOcc+nImp+1,nSites
!                        LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) + &
!                            SchmidtPert(i,a)*SchmidtPert(i,gtid(alpha))*FockSchmidt(a,gtid(alpha))
!                    enddo
!                enddo
!                LinearSystem(ExcitInd,nFCIDet+1) = LinearSystem(ExcitInd,nFCIDet+1) /  &
!                    sqrt(CoreVirtualNorm*CoreActiveNorm(alpha))
!                LinearSystem(nFCIDet+1,ExcitInd) = LinearSystem(ExcitInd,nFCIDet+1)
!
!                !Now for the coupling to the other semi-internal excitations
!                ExcitInd2 = nFCIDet + 1 + 4*nImp
!                do alphap = AS_Spin_start,AS_Spin_end
!                    ExcitInd2 = ExcitInd2 + 1
!                    !We shouldn't have coupling between excitations of the same spin
!                    if(mod(alphap,2).ne.mod(alpha,2)) cycle
!                    alphap_AS = alphap-2*(nOcc-nImp)
!
!                    Np1AlphapVec(:) = 0.0_dp
!                    !Create 1 and 2 body transition RDMs for the N-1 electron system where we have annihilated spin-orbital alphap
!                    !from |0>
!                    do i = 1,nFCIDet
!                        if(.not.btest(FCIBitList(i),alphap_AS-1)) then
!                            !We can create it
!                            Det = FCIBitList(i)
!                            call SQOperator(Det,alphap_AS,tSign,.false.)
!!                            Det = ibset(Det,alphap_AS-1)
!
!                            do j = 1,nNp1FCIDet
!                                if(Np1BitList(j).eq.Det) then
!                                    !We have found the orbital. Store the amplitude at this point
!                                    if(tSign) then
!                                        Np1AlphapVec(j) = -FullHamil(i,1)
!                                    else
!                                        Np1AlphapVec(j) = FullHamil(i,1)
!                                    endif
!                                    exit
!                                endif
!                            enddo
!                            if(j.gt.nNp1FCIDet) call stop_all(t_r,'Can not find appropritate determinant in N+1 space')
!                        endif
!                    enddo
!                    
!                    !We now have the wavefunction a_alpha |0> expressed in the N-1 electron FCI basis
!                    !Now create the 1 and 2 body density matrices from these amplitudes for the resulting N-1 electron wavefunction
!                    call CalcNp1_1RDM(Np1AlphaVec,Np1AlphapVec,Np1AlphaRDM)
!                    call CalcNp1_2RDM(Np1AlphaVec,Np1AlphapVec,Np1Alpha2RDM)
!!                    !Check that these correspond correctly to the higher-ranked density matrix in the N electron wavefunction
!!                    do i = 1,nImp*2
!!                        do j = 1,nImp*2
!!                            if(abs(Np1AlphaRDM(i,j)-(HL_2RDM(gtid(alpha),gtid(alphap),i,j)/2.0))) then
!!                                call stop_all(t_r'N-1 RDMs do not correspond to higher rank equivalent operators 2')
!!                            endif
!!                        enddo
!!                    enddo
!
!                    if(alpha.eq.alphap) then
!                        do i = 1,nOcc-nImp
!                            do j = 1,nOcc-nImp
!                                !Term 1
!                                LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) -   &
!                                    SchmidtPert(j,gtid(alphap))*SchmidtPert(i,gtid(alpha))*FockSchmidt(i,j)
!                                !Term 6a
!                                LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) +   &
!                                    SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))*FockSchmidt(j,j)*2.0_dp
!                            enddo
!                            !Term 6
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) +   &
!                                SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))*One_ElecE
!!                            do beta = 1,2*nImp
!!                                do gam = 1,2*nImp
!!                                    LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) +   &
!!                                        SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))* &
!!                                        !TODO: This should just be the one-electron energy
!!                                        FockSchmidt(beta+nOcc-nImp,gam+nOcc-nImp)*HL_1RDM(beta,gam)
!!                                enddo
!!                            enddo
!                        enddo
!                    endif
!
!                    do i = 1,nOcc-nImp
!                        do j = 1,nOcc-nImp
!                            !Term 2
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) + &
!                                SchmidtPert(j,gtid(alphap))*SchmidtPert(i,gtid(alpha))*FockSchmidt(i,j)*    &
!                                HL_1RDM(gtid(alpha_AS),gtid(alphap_AS))/2.0_dp
!                            !Term 7a
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) - &
!                                SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))*FockSchmidt(j,j)*    &
!                                HL_1RDM(gtid(alpha_AS),gtid(alphap_AS))
!                        enddo
!                        !Term 3
!                        LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) + &
!                            SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))*FockSchmidt(gtid(alphap),gtid(alpha))
!                    enddo
!
!                    do i = 1,nOcc-nImp
!                        do beta = 1,2*nImp
!                            !Term 4
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) -   &
!                                SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))*  &
!                                FockSchmidt(gtid(alphap),beta+nOcc-nImp)*HL_1RDM(gtid(alpha_AS),beta)/2.0_dp
!
!                            !Term 5
!                            LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) -   &
!                                SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))* &
!                                FockSchmidt(beta+nOcc-nImp,gtid(alpha))*HL_1RDM(beta,gtid(alphap_AS))/2.0_dp
!                        enddo
!                    enddo
!
!                    do i = 1,nOcc-nImp
!                        do beta = 1,2*nImp
!                            do gam = 1,2*nImp
!                                !Term 7
!                                LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) -   &
!                                    SchmidtPert(i,gtid(alphap))*SchmidtPert(i,gtid(alpha))* &
!                                    FockSchmidt(beta+nOcc-nImp,gam+nOcc-nImp)*  &
!                                    HL_2RDM(beta,gam,gtid(alpha_AS),gtid(alphap_AS))/2.0_dp
!                            enddo
!                        enddo
!                    enddo
!                    
!                    tmp = 0.0_dp
!                    do i = 1,nOcc-nImp           
!                        tmp = tmp + SchmidtPert(i,gtid(alpha))*SchmidtPert(i,gtid(alphap))
!                    enddo
!                    !Now for the two electron component: 
!                    do p = 1,2*nImp               
!                        do q = 1,2*nImp               
!                            do r = 1,2*nImp               
!                                do s = 1,2*nImp               
!                                    LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2) + &
!                                        tmp*Np1Alpha2RDM(p,q,r,s)*umat(umatind(p,r,q,s))
!                                enddo
!                            enddo
!                        enddo
!                    enddo
!                    LinearSystem(ExcitInd,ExcitInd2) = LinearSystem(ExcitInd,ExcitInd2)/    &
!                        (sqrt(CoreActiveNorm(alpha)*CoreActiveNorm(alphap)))
!                    LinearSystem(ExcitInd2,ExcitInd) = LinearSystem(ExcitInd,ExcitInd2)
!!                    if((ExcitInd2.eq.ExcitInd).and.(abs(tmp2-LinearSystem(ExcitInd,ExcitInd)).gt.1.0e-7_dp)) then
!!                        call stop_all(t_r,'Error in consistent diagonal elements for semi-internal excitations 2')
!!                    endif
!
!                enddo !alphap
!
!                !Final block: Coupling between two types of semi-internal excitations
!                !This is zero!
!            enddo  !Finish looping over active-virtual semi-internal excitations
!
!!            call writematrix(LinearSystem(1:nFCIDet,1:nFCIDet),'FCI Hessian',.true.)
!
!            !******************************************************************************************
!            ! OVERLAP matrix elements
!            !******************************************************************************************
!
!            !Calculate overlap matrix. 
!            !All functions normalized:
!            do i=1,nLinearSystem
!                Overlap(i,i) = 1.0_dp
!            enddo
!
!            !First the active-virtual excitations
!            ExcitInd = nFCIDet + 1
!            do alpha = AS_Spin_start,AS_Spin_end
!                ExcitInd = ExcitInd + 1
!
!                alpha_AS = alpha-2*(nOcc-nImp)  !The spin-orbital label of alpha, starting from 1
!
!                ExcitInd2 = nFCIDet + 1
!                do alphap = AS_Spin_start,AS_Spin_end
!                    ExcitInd2 = ExcitInd2 + 1
!                    alphap_AS = alphap-2*(nOcc-nImp)  !The spin-orbital label of alpha, starting from 1
!                    if(mod(alpha,2).ne.mod(alphap,2)) cycle
!
!                    if(ExcitInd2.eq.ExcitInd) then
!                        cycle
!                    else
!                        do a = nOcc+nImp+1,nSites
!                            !Accumulate S
!                            Overlap(ExcitInd,ExcitInd2) = Overlap(ExcitInd,ExcitInd2) + SchmidtPert(a,gtid(alpha))*  &
!                                SchmidtPert(a,gtid(alphap))*HL_1RDM(gtid(alpha_AS),gtid(alphap_AS)) / &
!                                (2.0_dp*sqrt(ActiveVirtualNorm(alpha)*ActiveVirtualNorm(alphap))) 
!                        enddo
!                    endif
!                enddo
!            enddo
!
!            !Now for the core-active semi internals
!            ExcitInd = nFCIDet + 1 + 4*nImp
!            do alpha = AS_Spin_start,AS_Spin_end
!                ExcitInd = ExcitInd + 1
!
!                alpha_AS = alpha-2*(nOcc-nImp)  !The spin-orbital label of alpha, starting from 1
!
!                ExcitInd2 = nFCIDet + 1 + 4*nImp
!                do alphap = AS_Spin_start,AS_Spin_end
!                    ExcitInd2 = ExcitInd2 + 1
!                    alphap_AS = alphap-2*(nOcc-nImp)  !The spin-orbital label of alpha, starting from 1
!                    if(mod(alpha,2).ne.mod(alphap,2)) cycle
!
!                    if(ExcitInd2.eq.ExcitInd) then
!                        !Overlap is 1
!                        cycle
!                    else
!                        do i = 1,nOcc-nImp
!                            Overlap(ExcitInd,ExcitInd2) = Overlap(ExcitInd,ExcitInd2) -     &
!                                SchmidtPert(i,gtid(alpha))*SchmidtPert(i,gtid(alphap))* &
!                                HL_1RDM(gtid(alphap_AS),gtid(alpha_AS))/(2.0_dp*    &
!                                sqrt(CoreActiveNorm(alpha)*CoreActiveNorm(alphap)))
!                        enddo
!                    endif
!                enddo
!            enddo
!
!            !Now, multiply the whole overlap matrix by (E_0 + Omega)
!            do i=1,nFCIDet
!                !In the det space, it is just the FCI energy
!                Overlap(i,i) = Overlap(i,i)*(Spectrum(1)+Omega)
!            enddo
!            Overlap(nFCIDet+1,nFCIDet+1) = Overlap(nFCIDet+1,nFCIDet+1)*Omega   !Zeroth order energy has not been included
!            !Remove what is expected to be the `zeroth-order' energy from the internal excitations
!            tmp = 0.0_dp
!            do i=1,nOcc-nImp
!                tmp = tmp + FockSchmidt(i,i)
!            enddo
!            tmp = tmp * 2.0_dp
!            !Remove the same energy contribution from all internal excitations
!            do alpha = nFCIDet + 2,nLinearSystem
!                do alphap = nFCIDet + 2,nLinearSystem
!                    Overlap(alpha,alphap) = Overlap(alpha,alphap)*(tmp + Spectrum(1) + Omega)
!                enddo
!            enddo
!
!            !Now construct [H-(E_0 + omega)S]
!            do i=1,nLinearSystem
!                do j=1,nLinearSystem
!                    LinearSystem(i,j) = LinearSystem(i,j) - Overlap(i,j)
!                enddo
!            enddo
!
!!            !Only now do we remove the ground state from the FCI space, since this is redundant
!!            do i=1,nFCIDet
!!                do j=1,nFCIDet
!!                    LinearSystem(j,i) = LinearSystem(j,i) - Spectrum(1)*FullHamil(i,1)*FullHamil(j,1)
!!                enddo
!!                !Finally, subtract the ground state energy from the diagonals, since we want to offset it.
!!                LinearSystem(i,i) = LinearSystem(i,i) - Spectrum(1)
!!                !Also, offset by the frequency of the transition
!!                LinearSystem(i,i) = Omega - LinearSystem(i,i)
!!            enddo
!!
!!            !The CV excitation never had its zeroth-order energy included
!!            !The CV excitation also needs to be offset by the transition frequency
!!            !Its overlap is 1
!!            LinearSystem(nFCIDet+1,nFCIDet+1) = Omega - LinearSystem(nFCIDet+1,nFCIDet+1)
!
!
!            !Check hamiltonian is hermitian
!            do i=1,nLinearSystem
!                do j=1,nLinearSystem
!                    if(abs(LinearSystem(i,j)-LinearSystem(j,i)).gt.1.0e-7_dp) then
!                        call stop_all(t_r,'Hessian not hermitian')
!                    endif
!                enddo
!            enddo
!
!            if(tDebug) call writematrix(LinearSystem,'LinearSystem',.true.)
!
!            !Attempt to spin-integrate matrix?
!            ExcitInd = nFCIDet + 1
!            do alpha = AS_Spin_start,AS_Spin_end,2  
!                ExcitInd = ExcitInd + 2 !Run through beta components
!
!                if(abs(LinearSystem(ExcitInd,1)-LinearSystem(ExcitInd-1,1)).gt.1.0e-7_dp) then
!                    call stop_all(t_r,'Error 1')
!                endif
!
!            enddo
!
!            !TEST! Diagonal approximation
!!            do i=1,nLinearSystem
!!                do j=1,nLinearSystem
!!                    if(i.eq.j) then
!!                        write(6,*) "Hess, ",i,LinearSystem(i,i)
!!                        cycle
!!                    endif
!!                    LinearSystem(i,j) = 0.0_dp
!!                enddo
!!            enddo
!!
!!
!!            !TEST! Just take FCI space
!!            do i=nFCIDet+1,nLinearSystem
!!                do j=nFCIDet+1,nLinearSystem
!!                    LinearSystem(i,j) = 0.0_dp
!!                enddo
!!            enddo
!!            !Check spectrum
!            nSize = nLinearSystem
!            allocate(temp(nSize,nSize))
!            temp(:,:) = LinearSystem(1:nSize,1:nSize)
!            allocate(Work(1))
!            if(ierr.ne.0) call stop_all(t_r,"alloc err")
!            allocate(W(nSize))
!            W(:)=0.0_dp
!            lWork=-1
!            info=0
!            call dsyev('V','U',nSize,temp,nSize,W,Work,lWork,info)
!            if(info.ne.0) call stop_all(t_r,'workspace query failed')
!            lwork=int(work(1))+1
!            deallocate(work)
!            allocate(work(lwork))
!            call dsyev('V','U',nSize,temp,nSize,W,Work,lWork,info)
!            if (info.ne.0) call stop_all(t_r,"Diag failed")
!            deallocate(work)
!            call writevector(W,'Hessian spectrum')
!            write(iunit2,*) Omega,W(:)
!            deallocate(W,temp)
!
!
!            if(.not.tDiagonalize) then
!                !Do not diagonalise. Instead, solve the linear system
!                allocate(temp_vec(nFCIDet))
!                temp_vec(:) = 0.0_dp
!                pertsitealpha = 2*pertsite-1
!                pertsitebeta = 2*pertsite
!                do i = 1,nFCIDet
!                    if(btest(FCIBitList(i),pertsitealpha-1)) then
!                        !This determinant is occupied
!                        temp_vec(i) = temp_vec(i) + FullHamil(i,1)
!                    endif
!                    if(btest(FCIBitList(i),pertsitebeta-1)) then
!                        temp_vec(i) = temp_vec(i) + FullHamil(i,1)
!                    endif
!                enddo
!
!                allocate(Projector(nFCIDet,nFCIDet))
!                Projector(:,:) = 0.0_dp
!                do i=1,nFCIDet
!                    Projector(i,i) = 1.0_dp
!                enddo
!!                    LinearSystem(j,i) = LinearSystem(j,i) - Spectrum(1)*FullHamil(i,1)*FullHamil(j,1)
!                do i=1,nFCIDet
!                    do j=1,nFCIDet
!                        Projector(j,i) = Projector(j,i) - FullHamil(i,1)*FullHamil(j,1)
!                    enddo
!                enddo
!
!                VGS(:) = 0.0_dp
!                call DGEMM('N','N',nFCIDet,1,nFCIDet,-1.0_dp,Projector,nFCIDet,temp_vec,nFCIDet,0.0_dp,VGS(1:nFCIDet),nFCIDet)
!
!                deallocate(temp_vec,Projector)
!
!                !Now solve these linear equations
!                call DGESV(nLinearSystem,1,LinearSystem,nLinearSystem,Pivots,VGS,nLinearSystem,info)
!!                call DGESV(nFCIDet,1,LinearSystem(1:nFCIDet,1:nFCIDet),nFCIDet,Pivots,VGS,nFCIDet,info)
!                if(info.ne.0) call stop_all(t_r,'Solving Linear system failed') 
!
!                ResponseFn = 0.0_dp
!                do i = 1,nFCIDet
!                    if(btest(FCIBitList(i),pertsitealpha-1)) then
!                        ResponseFn = ResponseFn + FullHamil(i,1)*VGS(i)
!                    endif
!                    if(btest(FCIBitList(i),pertsitebeta-1)) then
!                        ResponseFn = ResponseFn + FullHamil(i,1)*VGS(i)
!                    endif
!                enddo
!                write(iunit,*) Omega,ResponseFn
!
!            else
!                !Now we have the full hamiltonian. Diagonalize this fully
!                allocate(Work(1))
!                if(ierr.ne.0) call stop_all(t_r,"alloc err")
!                allocate(W(nLinearSystem))
!                W(:)=0.0_dp
!                lWork=-1
!                info=0
!                call dsyev('V','U',nLinearSystem,LinearSystem,nLinearSystem,W,Work,lWork,info)
!                if(info.ne.0) call stop_all(t_r,'workspace query failed')
!                lwork=int(work(1))+1
!                deallocate(work)
!                allocate(work(lwork))
!                call dsyev('V','U',nLinearSystem,LinearSystem,nLinearSystem,W,Work,lWork,info)
!                if (info.ne.0) call stop_all(t_r,"Diag failed")
!                deallocate(work)
!
!!                write(6,*) "First 10 MR-TDA-LR transition frequencies: "
!!                highbound = min(nLinearSystem,100)
!!                call writevector(W(1:highbound),'transition frequencies')
!!                
!                write(6,*) "Calculating residues: "
!
!                !Now find the transition moments
!                allocate(Residues(nLinearSystem))
!                Residues(:) = 0.0_dp
!
!                if(tResiduesFromRDM) then
!                    !Calc the residues from the full MO transition rdms.
!
!                    allocate(RDM1(nSites,nSites))
!                    allocate(RDM2(nSites,nSites))
!                    do n=1,nLinearSystem    !loop over states
!                        if(.not.tNonIntTest) then
!                            !If we have no FCI space at all, we should not even have core contributions
!                            !First, find the transition RDM between the ground FCI state, and state n
!                            call Calc1RDM(FullHamil(:,1),LinearSystem(1:nFCIDet,n),RDM1)
!                            !...and then vice versa
!                            call Calc1RDM(LinearSystem(1:nFCIDet,n),FullHamil(:,1),RDM2)
!                        else
!                            RDM1(:,:) = 0.0_dp
!                            RDM2(:,:) = 0.0_dp
!                        endif
!
!                        !Now add to the 1RDM the conributions from the fully IC core-active excitations
!                        do i=1,nOcc-nImp
!                            do a=nOcc+nImp+1,nSites
!
!                                RDM1(i,a) = RDM1(i,a) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(i,a)/sqrt(CoreVirtualNorm)
!                                !RDM1(a,i) = RDM1(a,i) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(a,i)/sqrt(CoreVirtualNorm)
!
!                                !RDM2(i,a) = RDM2(i,a) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(a,i)/sqrt(CoreVirtualNorm)
!                                RDM2(a,i) = RDM2(a,i) + 2.0_dp*LinearSystem(nFCIDet+1,n)*SchmidtPert(a,i)/sqrt(CoreVirtualNorm)
!
!                            enddo
!                        enddo
!
!                        !Now add to the 1RDM the contributions from the semi-internal IC active-virtual excitations
!                        ExcitInd = nFCIDet + 1
!                        do alpha = AS_Spin_start,AS_Spin_end
!                            ExcitInd = ExcitInd + 1
!                            do a=nOcc+nImp+1,nSites
!                                do beta = 1,4*nImp
!                                    if(mod(beta,2).ne.mod(alpha,2)) cycle
!                                    RDM1(gtid(beta)+nOcc-nImp,a) = RDM1(gtid(beta)+nOcc-nImp,a) + SchmidtPert(a,gtid(alpha))*    &
!                                        LinearSystem(ExcitInd,n)*HL_1RDM(gtid(beta),gtid(alpha)-nOcc+nImp) /   &
!                                        (2.0_dp*sqrt(ActiveVirtualNorm(alpha)))
!                                    RDM2(a,gtid(beta)+nOcc-nImp) = RDM2(a,gtid(beta)+nOcc-nImp) + SchmidtPert(a,gtid(alpha))*    &
!                                        LinearSystem(ExcitInd,n)*HL_1RDM(gtid(alpha)-nOcc+nImp,gtid(beta)) /   &
!                                        (2.0_dp*sqrt(ActiveVirtualNorm(alpha)))
!                                enddo
!                            enddo
!                        enddo
!
!                        !Now add tot he 1RDM the contributions from the semi-internal IC active-virtual excitations
!                        ExcitInd = nFCIDet + 1 + 4*nImp
!                        do alpha = AS_Spin_start,AS_Spin_end 
!                            ExcitInd = ExcitInd + 1
!                            do i=1,nOcc-nImp
!                                RDM1(i,gtid(alpha)) = RDM1(i,gtid(alpha)) + SchmidtPert(i,gtid(alpha))* &
!                                    LinearSystem(ExcitInd,n) / sqrt(CoreActiveNorm(alpha))
!                                
!                                RDM2(gtid(alpha),i) = RDM2(gtid(alpha),i) + SchmidtPert(i,gtid(alpha))* &
!                                    LinearSystem(ExcitInd,n) / sqrt(CoreActiveNorm(alpha))
!                                do beta = 1,4*nImp
!                                    if(mod(beta,2).ne.mod(alpha,2)) cycle
!                                    RDM1(i,gtid(beta)+nOcc-nImp) = RDM1(i,gtid(beta)+nOcc-nImp) - &
!                                        SchmidtPert(i,gtid(alpha))*LinearSystem(ExcitInd,n)*  &
!                                        HL_1RDM(gtid(alpha)-nOcc+nImp,gtid(beta))/(2.0_dp*sqrt(CoreActiveNorm(alpha)))
!                                    
!                                    RDM2(gtid(beta)+nOcc-nImp,i) = RDM2(gtid(beta)+nOcc-nImp,i) - &
!                                        SchmidtPert(i,gtid(alpha))*LinearSystem(ExcitInd,n)*  &
!                                        HL_1RDM(gtid(beta),gtid(alpha)-nOcc+nImp)/(2.0_dp*sqrt(CoreActiveNorm(alpha)))
!                                enddo
!                            enddo
!                        enddo
!                                    
!                        !Now, calculate the (pertsite,pertsite) component of this in the AO basis
!                        Res1 = 0.0_dp
!                        Res2 = 0.0_dp
!                        do i = 1,nSites
!                            do j = 1,nSites
!                                Res1 = Res1 + FullSchmidtBasis(pertsite,i)*RDM1(i,j)*FullSchmidtBasis(pertsite,j)
!                                Res2 = Res2 + FullSchmidtBasis(pertsite,i)*RDM2(i,j)*FullSchmidtBasis(pertsite,j)
!                            enddo
!                        enddo
!
!                        write(6,*) "Residues for state: ",n,Res1,Res2
!                        Residues(n) = Res1*Res2*Lambda
!                    enddo
!
!                else
!                    !Calculate residues directly from the wavefunction, since we know that the perturbation acts locally to one orbital only.
!                    pertsitealpha = 2*pertsite-1
!                    pertsitebeta = 2*pertsite
!                    
!                    do n=1,nLinearSystem
!
!                        do i = 1,nFCIDet
!                            if(btest(FCIBitList(i),pertsitealpha-1)) then
!                                Residues(n) = Residues(n) + (LinearSystem(i,n)*FullHamil(i,1))**2
!                            endif
!                            if(btest(FCIBitList(i),pertsitebeta-1)) then
!                                Residues(n) = Residues(n) + (LinearSystem(i,n)*FullHamil(i,1))**2
!                            endif
!
!                        enddo
!                        Residues(n) = Residues(n)*Lambda
!
!                    enddo
!
!                endif
!
!                ResponseFn = 0.0_dp
!                do i=1,nLinearSystem
!                    ResponseFn = ResponseFn + Residues(i)/W(i)
!    !                ResponseFn = ResponseFn - Residues(i)/(Omega+EDiff)
!                enddo
!                write(iunit,*) Omega,ResponseFn
!
!                if(tNonIntTest) then
!                    !Debug comparison info
!                    if(abs(W(1)-testham).gt.1.0e-7_dp) then
!                        write(6,*) "testham: ",testham
!                        write(6,*) "Eigenvalue: ",W(1)
!                        call stop_all(t_r,'eigenvalue not as expected')
!                    endif
!                    testham = 0.0_dp
!                    do i=1,nel
!                        do a=nel+1,nSites*2
!                            if(mod(i,2).ne.mod(a,2)) cycle
!                            i_spat = gtid(i)
!                            a_spat = gtid(a)
!
!                            testham = testham + ((FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2) / &
!                                ( (Omega**2/(FullHFEnergies(a_spat)-FullHFEnergies(i_spat))) - 2*Omega + &
!                                    (FullHFEnergies(a_spat)-FullHFEnergies(i_spat)))
!
!                        enddo
!                    enddo
!                    testham = testham / testnorm
!                    testham = Omega - testham
!                    !write(6,*) "****",abs(testham-(Omega-EDiff)),abs(testham),abs(testham-(Omega-EDiff))/abs(testham)
!                    if(abs(testham-(Omega-EDiff)).gt.1.0e-8_dp) then
!                        !write(6,*) "test denominator: ",testham
!                        !write(6,*) "Calculated Denominator: ",Omega-EDiff
!                        !call stop_all(t_r,'non interacting test denominator fail')
!                    endif
!
!                    testham = 0.0_dp
!                    do i=1,nel
!                        do a=nel+1,nSites*2
!                            if(mod(i,2).ne.mod(a,2)) cycle
!                            i_spat = gtid(i)
!                            a_spat = gtid(a)
!
!                            testham = testham + ((FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2) / &
!                                (Omega - (FullHFEnergies(a_spat)-FullHFEnergies(i_spat)))
!
!                        enddo
!                    enddo
!                    testham = testham / sqrt(testnorm)
!                    testham = testham*testham
!                    if(abs(testham-Residues(1)).gt.1.0e-7_dp) then
!                        write(6,*) "test residue: ",testham
!                        write(6,*) "Calculated Residue: ",Residues(1)
!                        call stop_all(t_r,'non interacting test residue fail')
!                    endif
!                endif
!
!            endif
!
!            Omega = Omega + Omega_Step
!        
!            if(.not.tDiagonalize) then
!                deallocate(VGS,Pivots)
!            else
!                if(tResiduesFromRDM) then
!                    deallocate(RDM1,RDM2)
!                endif
!                deallocate(Residues,W)
!            endif
!
!        enddo   !Enddo loop over omega
!
!        deallocate(LinearSystem,Overlap)
!        close(iunit)
!        close(iunit2)
!        if(tNonIntTest) call stop_all(t_r,'End of NonInt test')
!
    end subroutine NonIntContracted_TDA_MCLR
            

    !Find <Dj|p^+q|0> if tGSKet=T, or the other way around if false.
    !det is the label to the single determinant function Dj
    subroutine FindDeterminantTransRDM(det,Trans1RDM,tGSKet)
        use DetToolsData, only: FCIDetList,nFCIDet
        use DetTools, only: gtid,iGetExcitLevel,GetExcitation
        implicit none
        integer, intent(in) :: det
        logical, intent(in) :: tGSKet
        real(dp), intent(out) :: Trans1RDM(EmbSize,EmbSize)
        logical :: tSign
        integer :: IC,Ex(2),k,i
        character(len=*), parameter :: t_r='FindDeterminantTransRDM'

        if(tUHF) call stop_all(t_r,"Not sure this routine works with UHF")

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
                    Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) = Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) - HL_Vec(i)
                else
                    Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) = Trans1RDM(gtid(Ex(1)),gtid(Ex(2))) + HL_Vec(i)
                endif
            elseif(IC.eq.0) then
                !Same det
                if(i.ne.det) call stop_all(t_r,'Error here')
                do k=1,Elec
                    Trans1RDM(gtid(FCIDetList(k,i)),gtid(FCIDetList(k,i))) =        &
                        Trans1RDM(gtid(FCIDetList(k,i)),gtid(FCIDetList(k,i))) + HL_Vec(i)
                enddo
            endif
        enddo

    end subroutine FindDeterminantTransRDM
    
    
    !Calculate the 2RDM for an N+1 electron wavefunction. This can be transition
    !This is only over the ACTIVE space
    !These are spin-integrated, and have the form given in Helgakker
    subroutine CalcNp1_2RDM(Bra,Ket,RDM)
        use DetToolsData, only: nNp1FCIDet,Np1BitList,Np1FCIDetList
        use DetBitOps, only: CountBits
        use DetTools, only: gtid,GetExcitation
        implicit none
        real(dp), intent(in) :: Bra(nNp1FCIDet),Ket(nNp1FCIDet)
        real(dp), intent(out) :: RDM(nImp*2,nImp*2,nImp*2,nImp*2)
        integer :: Ex(2,2),i,j,IC,orbdiffs,k,l,kel,lel,temp
        logical :: tSign
        character(len=*), parameter :: t_r='CalcNp1_2RDM'
        
        if(tUHF) call stop_all(t_r,"Not sure this routine works with UHF")

        RDM(:,:,:,:) = 0.0_dp

        do i=1,nNp1FCIDet
            do j=1,nNp1FCIDet
                !Are they singles?
                orbdiffs = ieor(Np1BitList(i),Np1BitList(j))
                IC = CountBits(orbdiffs)
                if(mod(IC,2).ne.0) call stop_all(t_r,'Error here')
                IC = IC/2
                !IC = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,j),Elec)
                if(IC.eq.2) then
                    !Connected by a double
                    Ex(1,1) = 2
                    call GetExcitation(Np1FCIDetList(:,i),Np1FCIDetList(:,j),Elec+1,Ex,tSign)
                    if(mod(Ex(1,1),2).ne.mod(Ex(1,2),2)) then
                        !We have a mixed spin excitation
                        !Ensure that the spin of i is the same as the spin of b
                        !If its not, then reverse a and b and flip the sign
                        if(mod(Ex(1,1),2).ne.mod(Ex(2,2),2)) then
                            temp = Ex(2,2)
                            Ex(2,2) = Ex(2,1)
                            Ex(2,1) = temp
                            tSign = .not.tSign
                        endif
                    endif
                    if(tSign) then
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) -  &
                                Bra(i)*Ket(j)

                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) + &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) + &
                                Bra(i)*Ket(j)

                        endif
                    else
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) +  &
                                Bra(i)*Ket(j)
                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) - &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) - &
                                Bra(i)*Ket(j)
                        endif
                    endif

                elseif(IC.eq.1) then
                    !Connected by a single
                    Ex(1,1) = 1
                    call GetExcitation(Np1FCIDetList(:,i),Np1FCIDetList(:,j),Elec+1,Ex,tSign)
                    do k=1,elec+1
                        kel = gtid(Np1FCIDetList(k,i))

                        if(Np1FCIDetList(k,i).ne.Ex(1,1)) then
                            if(tSign) then
                                if(mod(Ex(1,1),2).eq.mod(Np1FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) + &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        Bra(i)*Ket(j)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)

                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        Bra(i)*Ket(j)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)
                                endif

                            else
                                if(mod(Ex(1,1),2).eq.mod(Np1FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) - &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        Bra(i)*Ket(j)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)
                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        Bra(i)*Ket(j)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)
                                endif

                            endif

                        endif
                    enddo
                elseif(IC.eq.0) then
                    !Same det
                    Ex(1,1) = 0
                    if(i.ne.j) call stop_all(t_r,'Error here')
                    do k=1,Elec+1
                        kel = gtid(Np1FCIDetList(k,i))
                        do l=k+1,Elec+1
                            lel = gtid(Np1FCIDetList(l,i))
                            if(Np1FCIDetList(k,i).eq.Np1FCIDetList(l,i)) cycle

                            if(mod(Np1FCIDetList(l,i),2).eq.mod(Np1FCIDetList(k,i),2)) then
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,kel,kel,lel) = RDM(lel,kel,kel,lel) - &
                                    Bra(i)*Ket(j)
                                RDM(kel,lel,lel,kel) = RDM(kel,lel,lel,kel) - &
                                    Bra(i)*Ket(j)
                            else
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    Bra(i)*Ket(j)
                            endif

                        enddo
                    enddo
                endif

            enddo
        enddo

    end subroutine CalcNp1_2RDM

    !Calculate the 1RDM for an N+1 electron wavefunction. This can be transition
    !This is only over the ACTIVE space
    !These are spin-integrated
    subroutine CalcNp1_1RDM(Bra,Ket,RDM)
        use DetToolsData, only: nNp1FCIDet,Np1BitList,Np1FCIDetList
        use DetBitOps, only: CountBits
        use DetTools, only: gtid,getexcitation
        implicit none
        real(dp), intent(in) :: Bra(nNp1FCIDet),Ket(nNp1FCIDet)
        real(dp), intent(out) :: RDM(nImp*2,nImp*2)
        !real(dp) :: trace
        integer :: i,j,ex(2),k,i_spat,a_spat,k_spat,IC    !,igetexcitlevel,nCore
        integer :: orbdiffs
        logical :: tSign
        character(len=*), parameter :: t_r='CalcNp1_1RDM'

        RDM(:,:) = 0.0_dp
        do i=1,nNp1FCIDet
            do j=1,nNp1FCIDet

                !Are they singles?
                orbdiffs = ieor(Np1BitList(i),Np1BitList(j))
                IC = CountBits(orbdiffs)
                if(mod(IC,2).ne.0) call stop_all(t_r,'Error here')
                IC = IC/2

                !IC = igetexcitlevel(Np1FCIDetList(:,i),Np1FCIDetList(:,j),elec)
                if(IC.eq.1) then
                    !Calculate parity
                    ex(1) = 1
                    call getexcitation(Np1FCIDetList(:,i),Np1FCIDetList(:,j),elec+1,ex,tSign)
                    i_spat = gtid(ex(1))
                    a_spat = gtid(ex(2))
                    if(tSign) then
                        RDM(i_spat,a_spat) = RDM(i_spat,a_spat) - Bra(i)*Ket(j)
                    else
                        RDM(i_spat,a_spat) = RDM(i_spat,a_spat) + Bra(i)*Ket(j)
                    endif
                elseif(IC.eq.0) then
                    do k=1,elec+1
                        k_spat = gtid(Np1FCIDetList(k,i)) 
                        RDM(k_spat,k_spat) = RDM(k_spat,k_spat) + Bra(i)*Ket(i)
                    enddo
                endif
            enddo
        enddo
    end subroutine CalcNp1_1RDM
    
    !Calculate the 2RDM for an N-1 electron wavefunction. This can be transition
    !This is only over the ACTIVE space
    !These are spin-integrated, and have the form given in Helgakker
    subroutine CalcNm1_2RDM(Bra,Ket,RDM)
        use DetToolsData, only: nNm1FCIDet,Nm1BitList,Nm1FCIDetList
        use DetBitOps, only: CountBits
        use DetTools, only: gtid,GetExcitation
        implicit none
        real(dp), intent(in) :: Bra(nNm1FCIDet),Ket(nNm1FCIDet)
        real(dp), intent(out) :: RDM(nImp*2,nImp*2,nImp*2,nImp*2)
        integer :: Ex(2,2),i,j,IC,orbdiffs,k,l,kel,lel,temp
        logical :: tSign
        character(len=*), parameter :: t_r='CalcNm1_2RDM'

        if(tUHF) call stop_all(t_r,"Not sure this routine works with UHF")

        RDM(:,:,:,:) = 0.0_dp

        do i=1,nNm1FCIDet
            do j=1,nNm1FCIDet
                !Are they singles?
                orbdiffs = ieor(Nm1BitList(i),Nm1BitList(j))
                IC = CountBits(orbdiffs)
                if(mod(IC,2).ne.0) call stop_all(t_r,'Error here')
                IC = IC/2
                !IC = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,j),Elec)
                if(IC.eq.2) then
                    !Connected by a double
                    Ex(1,1) = 2
                    call GetExcitation(Nm1FCIDetList(:,i),Nm1FCIDetList(:,j),Elec-1,Ex,tSign)
                    if(mod(Ex(1,1),2).ne.mod(Ex(1,2),2)) then
                        !We have a mixed spin excitation
                        !Ensure that the spin of i is the same as the spin of b
                        !If its not, then reverse a and b and flip the sign
                        if(mod(Ex(1,1),2).ne.mod(Ex(2,2),2)) then
                            temp = Ex(2,2)
                            Ex(2,2) = Ex(2,1)
                            Ex(2,1) = temp
                            tSign = .not.tSign
                        endif
                    endif
                    if(tSign) then
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) -  &
                                Bra(i)*Ket(j)

                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) + &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) + &
                                Bra(i)*Ket(j)

                        endif
                    else
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) +  &
                                Bra(i)*Ket(j)
                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) - &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) - &
                                Bra(i)*Ket(j)
                        endif
                    endif

                elseif(IC.eq.1) then
                    !Connected by a single
                    Ex(1,1) = 1
                    call GetExcitation(Nm1FCIDetList(:,i),Nm1FCIDetList(:,j),Elec-1,Ex,tSign)
                    do k=1,elec-1
                        kel = gtid(Nm1FCIDetList(k,i))

                        if(Nm1FCIDetList(k,i).ne.Ex(1,1)) then
                            if(tSign) then
                                if(mod(Ex(1,1),2).eq.mod(Nm1FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) + &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        Bra(i)*Ket(j)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)

                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        Bra(i)*Ket(j)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)
                                endif

                            else
                                if(mod(Ex(1,1),2).eq.mod(Nm1FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) - &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        Bra(i)*Ket(j)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)
                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        Bra(i)*Ket(j)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)
                                endif

                            endif

                        endif
                    enddo
                elseif(IC.eq.0) then
                    !Same det
                    Ex(1,1) = 0
                    if(i.ne.j) call stop_all(t_r,'Error here')
                    do k=1,Elec-1
                        kel = gtid(Nm1FCIDetList(k,i))
                        do l=k+1,Elec-1
                            lel = gtid(Nm1FCIDetList(l,i))
                            if(Nm1FCIDetList(k,i).eq.Nm1FCIDetList(l,i)) cycle

                            if(mod(Nm1FCIDetList(l,i),2).eq.mod(Nm1FCIDetList(k,i),2)) then
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,kel,kel,lel) = RDM(lel,kel,kel,lel) - &
                                    Bra(i)*Ket(j)
                                RDM(kel,lel,lel,kel) = RDM(kel,lel,lel,kel) - &
                                    Bra(i)*Ket(j)
                            else
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    Bra(i)*Ket(j)
                            endif

                        enddo
                    enddo
                endif

            enddo
        enddo

    end subroutine CalcNm1_2RDM

    !Calculate the 1RDM for an N-1 electron wavefunction. This can be transition
    !This is only over the ACTIVE space
    !These are spin-integrated
    subroutine CalcNm1_1RDM(Bra,Ket,RDM)
        use DetToolsData, only: nNm1FCIDet,Nm1BitList,Nm1FCIDetList
        use DetBitOps, only: CountBits
        use DetTools, only: gtid,getexcitation
        implicit none
        real(dp), intent(in) :: Bra(nNm1FCIDet),Ket(nNm1FCIDet)
        real(dp), intent(out) :: RDM(nImp*2,nImp*2)
        !real(dp) :: trace
        integer :: i,j,ex(2),k,i_spat,a_spat,k_spat,IC    !,igetexcitlevel,nCore
        integer :: orbdiffs
        logical :: tSign
        character(len=*), parameter :: t_r='CalcNm1_1RDM'

        RDM(:,:) = 0.0_dp
        do i=1,nNm1FCIDet
            do j=1,nNm1FCIDet

                !Are they singles?
                orbdiffs = ieor(Nm1BitList(i),Nm1BitList(j))
                IC = CountBits(orbdiffs)
                if(mod(IC,2).ne.0) call stop_all(t_r,'Error here')
                IC = IC/2

                !IC = igetexcitlevel(Nm1FCIDetList(:,i),Nm1FCIDetList(:,j),elec)
                if(IC.eq.1) then
                    !Calculate parity
                    ex(1) = 1
                    call getexcitation(Nm1FCIDetList(:,i),Nm1FCIDetList(:,j),elec-1,ex,tSign)
                    i_spat = gtid(ex(1))
                    a_spat = gtid(ex(2))
                    if(tSign) then
                        RDM(i_spat,a_spat) = RDM(i_spat,a_spat) - Bra(i)*Ket(j)
                    else
                        RDM(i_spat,a_spat) = RDM(i_spat,a_spat) + Bra(i)*Ket(j)
                    endif
                elseif(IC.eq.0) then
                    do k=1,elec-1
                        k_spat = gtid(Nm1FCIDetList(k,i)) 
                        RDM(k_spat,k_spat) = RDM(k_spat,k_spat) + Bra(i)*Ket(i)
                    enddo
                endif
            enddo
        enddo
    end subroutine CalcNm1_1RDM

    !Calculate the (transition) 1RDM.
    !Do this very naively to start with.
    subroutine Calc1RDM(Bra,Ket,RDM)
        use DetToolsData, only: nFCIDet,FCIDetList
        use DetTools, only: gtid,igetexcitlevel,getexcitation
        implicit none
        real(dp), intent(in) :: Bra(nFCIDet),Ket(nFCIDet)
        real(dp), intent(out) :: RDM(nSites,nSites)
        !real(dp) :: trace
        integer :: nCore,i,j,ex(2),k,i_spat,a_spat,k_spat,IC
        logical :: tSign
        character(len=*), parameter :: t_r='Calc1RDM'
        
        if(tUHF) call stop_all(t_r,"Not sure this routine works with UHF")

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
        use DetTools, only: GetHElement,GetExcitation
        implicit none
        integer :: nCoreVirt,nCoreActive,nActiveVirt,nLinearSystem,ierr,info
        integer :: i,j,CoreNEl,ind2,ind1,a,x,Ex(2),b,pertsite
        logical :: tSign_a,tSign_b
        integer, allocatable :: Pivots(:),RefCore(:),Excit1_a(:),Excit1_b(:)
        integer, allocatable :: Excit2_a(:),Excit2_b(:)
        real(dp), allocatable :: LinearSystem(:,:),Overlap(:,:),Response(:)
        real(dp), allocatable :: ResponseSaved(:), LinearSystemSaved(:,:)
        real(dp) :: Omega,ResponseVal
        character(len=*), parameter :: t_r='SolveDMETResponse'

        !Assume initially that the perturbation acts only at site 1
        pertsite = 1

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
        do i = 1,nFCIDet
            do j = 1,nFCIDet
                call GetHElement(FCIDetList(:,i),FCIDetList(:,j),Elec,LinearSystem(i,j))
            enddo
        enddo
        !Finally, subtract the ground state energy from the diagonals, since we want to offset it.
        do i=1,nFCIDet
            LinearSystem(i,i) = LinearSystem(i,i) - HL_Energy 
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
                        LinearSystem(j,nFCIDet+(2*ind2)-1) = -HL_Vec(j)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2)-1,j) = -HL_Vec(j)*FockSchmidt(i,a)
                    else
                        LinearSystem(j,nFCIDet+(2*ind2)-1) = HL_Vec(j)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2)-1,j) = HL_Vec(j)*FockSchmidt(i,a)
                    endif
                    !beta-beta excitation is 2*ind2
                    if(tSign_b) then
                        LinearSystem(j,nFCIDet+(2*ind2)) = -HL_Vec(j)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2),j) = -HL_Vec(j)*FockSchmidt(i,a)
                    else
                        LinearSystem(j,nFCIDet+(2*ind2)) = HL_Vec(j)*FockSchmidt(i,a)
                        LinearSystem(nFCIDet+(2*ind2),j) = HL_Vec(j)*FockSchmidt(i,a)
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
                    Response(i) = Response(i) + Lambda*HL_Vec(i)
                elseif(FCIDetList(j,i).eq.(pertsite*2)) then
                    !beta spin of the perturbation
                    Response(i) = Response(i) + Lambda*HL_Vec(i)
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
                        ResponseVal = ResponseVal + HL_Vec(i)*Response(i)*Lambda
                    elseif(FCIDetList(j,i).eq.(pertsite*2)) then
                        ResponseVal = ResponseVal + HL_Vec(i)*Response(i)*Lambda
                    endif
                enddo
            enddo
                        
            write(6,*) "Response function: ",Omega,ResponseVal

            Omega = Omega + Omega_Step

        enddo


    end subroutine TDA_MCLR



    !Change the RHS so that we multiply by the hermitian conjugate of the shifted matrix, so that we are dealing with a purely hermitian problem
    subroutine setup_RHS(n,V0,Trans_V0)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: V0(n)
        complex(dp), intent(out) :: Trans_V0(n)
        integer :: i,j,k
        character(len=*), parameter :: t_r='setup_RHS'
        
        if(tCompressedMats.and.(.not.associated(zDirMV_Mat_cmprs))) call stop_all(t_r,'Compressed matrix not associated')
        if(tCompressedMats.and.(.not.associated(zDirMV_Mat_cmprs_inds))) call stop_all(t_r,'Compressed indices not associated')
        if((.not.tCompressedMats).and.(.not.associated(zDirMV_Mat))) call stop_all(t_r,'Matrix not associated!')

!        call writevectorint(zDirMV_Mat_cmprs_inds,'zDirMV_Mat_cmprs_inds')
!        call writevectorcomp(zDirMV_Mat_cmprs,'zDirMV_Mat_cmprs')

!        Trans_V0(:) = V0(:)
!        call ZGEMV('C',n,n,zone,zDirMV_Mat,n,V0,1,-dconjg(zShift),Trans_V0,1)
        if(tCompressedMats) then
            if(zDirMV_Mat_cmprs_inds(1).ne.(n+2)) call stop_all(t_r,'Mismatched vector and matrix')
            do i = 1,n
                Trans_V0(i) = conjg(zDirMV_Mat_cmprs(i))*V0(i)
            enddo
            do i = 1,n
                do k = zDirMV_Mat_cmprs_inds(i),zDirMV_Mat_cmprs_inds(i+1)-1
                    j = zDirMV_Mat_cmprs_inds(k)
                    Trans_V0(j) = Trans_V0(j) + conjg(zDirMV_Mat_cmprs(k))*V0(i)
                enddo
            enddo
        else
            call ZGEMV('C',n,n,zone,zDirMV_Mat,n,V0,1,zzero,Trans_V0,1)
        endif

    end subroutine setup_RHS


    !Calculate the non-interacting greens functions, with empirical self-energy correction
    !Also, find the contraction coefficients for all impurity-impurity greens functions
    !Then, calculate the new 1-electron hamiltonian (with implicit self-energy term) for 
    !the embedded basis, as well as the core hamiltonian.
    !The NI_LRMat_Cre and NI_LRMat_Ann matrices are returned as the non-interacting LR
    !The global matrices SchmidtPertGF_Cre and SchmidtPertGF_Ann are filled, as well as the
    !new one-electron hamiltonians for the interacting problem: Emb_h0v_SE and FockSchmidt_SE
    subroutine FindNI_Charged(Omega,NI_LRMat_Cre,NI_LRMat_Ann)
        use mat_tools, only: add_localpot_comp_inplace
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        implicit none
        real(dp), intent(in) :: Omega
        complex(dp), intent(out) :: NI_LRMat_Cre(nImp,nImp),NI_LRMat_Ann(nImp,nImp)
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: AO_OneE_Ham(:,:),W_Vals(:),RVec(:,:),LVec(:,:),FullSchmidtTrans_C(:,:)
        complex(dp), allocatable :: HFPertBasis_Ann_Ket(:,:),HFPertBasis_Cre_Ket(:,:),temp(:,:),temp2(:,:)
        complex(dp), allocatable :: HFPertBasis_Ann_Bra(:,:),HFPertBasis_Cre_Bra(:,:),cWork(:),EmbeddedBasis_C(:,:)
        integer :: lwork,info,i,a,pertBra,j,pertsite,nVirt,CoreEnd,VirtStart,ActiveStart,ActiveEnd,nCore
        complex(dp) :: zdotc,test
        character(len=*), parameter :: t_r='FindNI_Charged'

        if(.not.tAllImp_LR) then
            call stop_all(t_r,"Should not be in this routine if you don't want to calc all greens functions")
        endif

!        !Check hermiticity
!        do i = 1,nImp
!            do j = 1,nImp
!                if(i.eq.j) cycle
!                if(abs(dconjg(SelfEnergy_Imp(i,j))-SelfEnergy_Imp(j,i)).gt.1.0e-6_dp) then
!                    write(6,*) "Losing off-diagonal hermiticity of the self-energy"
!                    call writematrixcomp(SelfEnergy_Imp,'Self Energy',.true.)
!                endif
!            enddo
!        enddo

        allocate(AO_OneE_Ham(nSites,nSites))
        AO_OneE_Ham(:,:) = h0v_SE(:,:)
!        AO_OneE_Ham(:,:) = zzero
        !Stripe the complex (-)self-energy through the AO one-electron hamiltonian
!        call add_localpot_comp(h0v,AO_OneE_Ham,SelfEnergy_Imp,tAdd=.false.)

        !Now, diagonalize the resultant non-hermitian one-electron hamiltonian
        allocate(W_Vals(nSites))
        allocate(RVec(nSites,nSites))
        allocate(LVec(nSites,nSites))
        RVec = zzero
        LVec = zzero
        W_Vals = zzero
        allocate(Work(max(1,2*nSites)))
        allocate(cWork(1))
        lwork = -1
        info = 0
        call zgeev('V','V',nSites,AO_OneE_Ham,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
        lwork = int(abs(cWork(1)))+1
        deallocate(cWork)
        allocate(cWork(lWork))
        call zgeev('V','V',nSites,AO_OneE_Ham,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
        if(info.ne.0) call stop_all(t_r,'Diag of H - SE failed')
        deallocate(work,cWork)

        !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
        !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
        call Order_zgeev_vecs(W_Vals,LVec,RVec)
        !call writevectorcomp(W_Vals,'Eigenvalues ordered')
        !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
        call Orthonorm_zgeev_vecs(nSites,W_Vals,LVec,RVec)

        if(tCheck) then
            !*** TEST ***
            allocate(temp(nSites,nSites))
            allocate(temp2(nSites,nSites))
            !call writematrixcomp(SelfEnergy_Imp,'SelfEnergy',.true.)
            !call writevectorcomp(W_Vals,'Eigenvalues')
!            AO_OneE_Ham(:,:) = zzero
!            call add_localpot_comp(h0v,AO_OneE_Ham,SelfEnergy_Imp,tAdd=.false.)
            AO_OneE_Ham(:,:) = h0v_SE(:,:)
    
            call ZGEMM('C','N',nSites,nSites,nSites,zone,LVec,nSites,AO_OneE_Ham,nSites,zzero,temp,nSites)
            do j = 1,nSites
                do a = 1,nSites
                    if(abs(temp(j,a)-(W_Vals(j)*dconjg(LVec(a,j)))).gt.1.0e-8) then
                        write(6,*) "Eigenvector: ",j
                        write(6,*) "Component: ",a
                        write(6,*) temp(j,a),W_Vals(j)*dconjg(LVec(a,j))
                        call stop_all(t_r,'LVecs not computed correctly')
                    endif
                enddo
            enddo
    !
    !        write(6,*) "Left eigenvectors computed correctly..."
    !
            call ZGEMM('N','N',nSites,nSites,nSites,zone,AO_OneE_Ham,nSites,RVec,nSites,zzero,temp,nSites)
            do j = 1,nSites
                do a = 1,nSites
                    if(abs(temp(a,j)-(W_Vals(j)*RVec(a,j))).gt.1.0e-8) then
                        write(6,*) "Eigenvector: ",j
                        write(6,*) "Component: ",a
                        write(6,*) temp(a,j),W_Vals(j)*RVec(a,j)
                        call stop_all(t_r,'RVecs not computed correctly')
                    endif
                enddo
            enddo
    !
    !        write(6,*) "Right eigenvectors computed correctly..."
    !
            !Test normalization and biorthogonality of eigenstate pairs 
            do i = 1,nSites
                do j = 1,nSites
                    test = zdotc(nSites,LVec(:,i),1,RVec(:,j),1)
                    !write(6,*) "LVec: ",i,"RVec: ",j,test
                    if((i.eq.j).and.(abs(test-zone).gt.1.0e-8_dp)) then
                        write(6,*) "Normalization not maintained"
                        write(6,*) "LVec: ",i,"RVec: ",j,test
                        call stop_all(t_r,'Normalization error')
                    elseif((i.ne.j).and.(abs(test).gt.1.0e-8_dp)) then
                        write(6,*) "Orthogonality not maintained"
                        write(6,*) "Overlap: ",abs(test)
                        write(6,*) "LVec: ",i,"RVec: ",j,test
                        call stop_all(t_r,'Orthogonality error')
                    endif
                enddo
            enddo

    !        call writevectorcomp(W_Vals,'Eigenvalues ordered and orthonormed')
            !call writematrixcomp(AO_OneE_Ham,'h0v',.false.)
            call zGEMM('C','N',nSites,nSites,nSites,zone,LVec,nSites,AO_OneE_Ham,nSites,zzero,temp,nSites)
            call zGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RVec,nSites,zzero,temp2,nSites)
    !        call writematrixcomp(temp2,'L* H R ordered and orthonormed',.false.)
            do i = 1,nSites
                do j = 1,nSites
                    if((i.ne.j).and.abs(temp2(j,i)).gt.1.0e-9_dp) then
                        write(6,*) "i,j: ",i,j,temp2(j,i)
                        call stop_all(t_r,'L* H R does not reproduce diagonal matrix')
                    elseif((i.eq.j).and.(abs(temp2(j,i)-W_Vals(j))).gt.1.0e-9_dp) then
                        write(6,*) "i,j: ",i,j,temp2(j,i),W_Vals(j),abs(temp2(j,i)-W_Vals(j))
                        call writevectorcomp(W_Vals,'Ordered eigenvalues')
                        call stop_all(t_r,'L* H R does not reproduce eigenvalues')
                    endif
                enddo
            enddo
            deallocate(temp,temp2)
            
        endif !Endif tCheck

        NI_LRMat_Cre(:,:) = zzero
        NI_LRMat_Ann(:,:) = zzero 

!        call writevectorcomp(W_Vals,'HF energies')
!        call writematrixcomp(RVec,'Orbitals',.false.)
        
        !Schmidt basis bounds
        nVirt = nSites-nOcc-nImp   
        nCore = nOcc-nImp
        CoreEnd = nOcc-nImp
        VirtStart = nOcc+nImp+1
        ActiveStart = nOcc-nImp+1
        ActiveEnd = nOcc+nImp

        !Memory to temperarily store the first order wavefunctions of each impurity site, in the right MO basis (For the Kets)
        !and the left MO space (for the Bras)
        allocate(HFPertBasis_Ann_Bra(1:nOcc,nImp))
        allocate(HFPertBasis_Cre_Bra(nOcc+1:nSites,nImp))
        allocate(HFPertBasis_Ann_Ket(1:nOcc,nImp))
        allocate(HFPertBasis_Cre_Ket(nOcc+1:nSites,nImp))
        HFPertBasis_Ann_Bra(:,:) = zzero
        HFPertBasis_Cre_Bra(:,:) = zzero
        HFPertBasis_Ann_Ket(:,:) = zzero
        HFPertBasis_Cre_Ket(:,:) = zzero
            
        !call writematrixcomp(RVec(1:nImp,1:nSites),'RVec(1:nImp,1:nOcc) - LR',.true.)
        !Now, form the non-interacting greens functions (but with u *and* self-energy)
        do pertsite = 1,nImp
            !Form the set of non-interacting first order wavefunctions from the new one-electron h for both Bra and Ket versions
            !TODO: Check whether I want this dconj around the LVec? Similarly with the RVec, and -dDelta's
            do i = 1,nOcc
                HFPertBasis_Ann_Ket(i,pertsite) = dconjg(LVec(pertsite,i))/(dcmplx(Omega,dDelta)-W_Vals(i))
                HFPertBasis_Ann_Bra(i,pertsite) = RVec(pertsite,i)/(dcmplx(Omega,-dDelta)-W_Vals(i))
            enddo
            do a = nOcc+1,nSites
                HFPertBasis_Cre_Ket(a,pertsite) = dconjg(LVec(pertsite,a))/(dcmplx(Omega,dDelta)-W_Vals(a))
                HFPertBasis_Cre_Bra(a,pertsite) = RVec(pertsite,a)/(dcmplx(Omega,-dDelta)-W_Vals(a))
            enddo

            !Run over operators acting of the Bra in the impurity space
            do pertBra = 1,nImp
                !Now perform the set of dot products of <0|V* with |1> for all combinations of sites
                do i = 1,nOcc
                    NI_LRMat_Ann(pertsite,pertBra) = NI_LRMat_Ann(pertsite,pertBra) +   &
                        RVec(pertBra,i)*HFPertBasis_Ann_Ket(i,pertsite)
                enddo
                do a = nOcc+1,nSites
                    NI_LRMat_Cre(pertsite,pertBra) = NI_LRMat_Cre(pertsite,pertBra) +   &
                        RVec(pertBra,a)*HFPertBasis_Cre_Ket(a,pertsite)
                enddo
            enddo
        enddo

!        write(6,*) "1,1 :",NI_LRMat_Ann(1,1)+NI_LRMat_Cre(1,1)
!        write(6,*) "Ann 1,2: ",NI_LRMat_Ann(1,2)
!        write(6,*) "Ann 2,1: ",NI_LRMat_Ann(2,1)
!        write(6,*) "Cre 1,2: ",NI_LRMat_Cre(1,2)
!        write(6,*) "Cre 2,1: ",NI_LRMat_Cre(2,1)
!        write(6,*) "Sum 1,2: ",NI_LRMat_Ann(1,2) + NI_LRMat_Cre(1,2)
!        write(6,*) "Sum 2,1: ",NI_LRMat_Ann(2,1) + NI_LRMat_Cre(2,1)
!
        !Now, we need to project these NI wavefunctions into the schmidt basis
        !We want to rotate the vectors, expressed in the 'right' basis, back into that AO basis.
        !If we can do that, we can easily rotate into the Schmidt basis.
        allocate(FullSchmidtTrans_C(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                FullSchmidtTrans_C(j,i) = dcmplx(FullSchmidtBasis(j,i),0.0_dp)  !(ao,schmidt)
            enddo
        enddo

        allocate(temp(nSites,nImp))
        call ZGEMM('N','N',nSites,nImp,nOcc,zone,RVec(:,1:nOcc),nSites,HFPertBasis_Ann_Ket(1:nOcc,:),nOcc,zzero,  &
            temp,nSites)
            !temp is now the (nSites,nImp) rotated HFPertBasis_Ann into the AO basis
            !Now rotate this into the occupied schmidt basis
        call ZGEMM('T','N',nCore,nImp,nSites,zone,FullSchmidtTrans_C(:,1:CoreEnd),nSites,temp,nSites,zzero,   &
            SchmidtPertGF_Ann_Ket(1:CoreEnd,:),nCore)
        !Now do the same for the bra contraction coefficients
        allocate(temp2(nSites,1:nOcc))
        do i=1,nOcc
            do j=1,nSites
                temp2(j,i) = dconjg(LVec(j,i))
            enddo
        enddo
        call ZGEMM('N','N',nSites,nImp,nOcc,zone,temp2(:,1:nOcc),nSites,HFPertBasis_Ann_Bra(1:nOcc,:),nOcc,zzero,    &
            temp,nSites)
        deallocate(temp2)
        call ZGEMM('T','N',nCore,nImp,nSites,zone,FullSchmidtTrans_C(:,1:CoreEnd),nSites,temp,nSites,zzero,   &
            SchmidtPertGF_Ann_Bra(1:CoreEnd,:),nCore)

        !Do the same with the particle NI GF
        call ZGEMM('N','N',nSites,nImp,nSites-nOcc,zone,RVec(:,nOcc+1:nSites),nSites,  &
            HFPertBasis_Cre_Ket,nSites-nOcc,zzero,temp,nSites)
        !Now rotate into schmidt basis
        call ZGEMM('T','N',nVirt,nImp,nSites,zone,FullSchmidtTrans_C(:,VirtStart:nSites),nSites,temp,nSites,zzero,    &
            SchmidtPertGF_Cre_Ket(VirtStart:nSites,:),nVirt)
        !Now for the Bra version of the particle NI GF
        allocate(temp2(nSites,nOcc+1:nSites))
        do i = nOcc+1,nSites
            do j = 1,nSites
                temp2(j,i) = dconjg(LVec(j,i))
            enddo
        enddo
        call ZGEMM('N','N',nSites,nImp,nSites-nOcc,zone,temp2(:,nOcc+1:nSites),nSites,  &
            HFPertBasis_Cre_Bra,nSites-nOcc,zzero,temp,nSites)
        !Now rotate into schmidt basis
        call ZGEMM('T','N',nVirt,nImp,nSites,zone,FullSchmidtTrans_C(:,VirtStart:nSites),nSites,temp,nSites,zzero,    &
            SchmidtPertGF_Cre_Bra(VirtStart:nSites,:),nVirt)
        deallocate(temp,temp2)
        
        !TODO: Check that in the absence of a self-energy, the Bra and Ket are complex conjugates of each other.

        if(.not.tNoHL_SE) then
            !Now, we need to find the new 1 electron hamiltonian for the interacting problem.
            !This wants to have the self-energy projected into the schmidt basis, but only over the bath and
            !core sites, not the impurity.

            !First, deal with the embedded basis
            !Stripe the self energy through the space
            allocate(temp(nSites,nSites))
            temp(:,:) = zzero
            call add_localpot_comp_inplace(temp,SelfEnergy_Imp,tAdd=.true.)
            !temp is now just the self-energy striped through the space, in the AO basis
            !Rotate this into the embedding basis
            allocate(temp2(EmbSize,nSites))
            allocate(EmbeddedBasis_C(nSites,EmbSize))
            do i = 1,EmbSize
                do j = 1,nSites
                    EmbeddedBasis_C(j,i) = dcmplx(EmbeddedBasis(j,i),0.0_dp)
                enddo
            enddo
            call ZGEMM('T','N',EmbSize,nSites,nSites,zone,EmbeddedBasis_C,nSites,temp,nSites,zzero,temp2,EmbSize)
            call ZGEMM('N','N',EmbSize,EmbSize,nSites,zone,temp2,EmbSize,EmbeddedBasis_C,nSites,zzero,Emb_h0v_SE,EmbSize)
            deallocate(temp2,EmbeddedBasis_C)
            
            !Now, zero out the self-energy contribution over the impurity site
            Emb_h0v_SE(1:nImp,1:nImp) = zzero
            !Now subtract this self-energy correction from the normal one-electron hamiltonian in the embedded space (with corr pot)
            !The resulting one-electron embedding potential is neither real, nor hermitian
            Emb_h0v_SE(:,:) = Emb_h0v(:,:) - Emb_h0v_SE(:,:)

            !But what about the core hamiltonian?
            !Here, we do need to include
            !the effect of the self-energy on the core orbitals, since it is explicitly a coupling between the
            !(non-bath) environment and the embedded basis which we want to consider, and we are not changing
            !the bath orbital to account for it as we do in the ground state problem.
            !Calculate the new 1-electron hamiltonian

            !This now gives the one-electron hamiltonian in the AO basis, with the self-energy subtracted.
            !Rotate from the AO basis, to the Schmidt basis
            call ZGEMM('T','N',nSites,nSites,nSites,zone,FullSchmidtTrans_C,nSites,h0v_SE,nSites,zzero,temp,nSites)
            call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,FullSchmidtTrans_C,nSites,zzero,FockSchmidt_SE,nSites)

            do i=nImp+1,EmbSize
                do j=1,EmbSize
                    !The bath orbitals, and coupling to the impurity site blocks of the one-electron hamiltonain should now agree between FockSchmidt_SE and Emb_h0v_SE surely?
                    if(abs(Emb_h0v_SE(j,i)-FockSchmidt_SE(nOcc-nImp+j,nOcc-nImp+i)).gt.1.0e-8_dp) then
                        call writematrixcomp(Emb_h0v_SE,'Emb_h0v_SE',.true.)
                        call writematrixcomp(FockSchmidt_SE(nOcc-nImp+1:nOcc+nImp,nOcc-nImp+1:nOcc+nImp),   &
                            'Fock_SE Embedded system',.true.)
                        write(6,*) "The above should be the same in the bath and coupling blocks"
                        write(6,*) "j,i: ",j,i,Emb_h0v_SE(j,i),FockSchmidt_SE(nOcc-nImp+j,nOcc-nImp+i), &
                            abs(Emb_h0v_SE(j,i)-FockSchmidt_SE(nOcc-nImp+j,nOcc-nImp+i))
                        call stop_all(t_r,'One electron hamiltonians with self energy not consistent')
                    endif
                enddo
            enddo

            !Now, convert this into a complex number (for easy ZGEMM'ing), and store it in the global array
            FockSchmidt_SE_CC(1:CoreEnd,1:CoreEnd) = FockSchmidt_SE(1:CoreEnd,1:CoreEnd)
            FockSchmidt_SE_VV(VirtStart:nSites,VirtStart:nSites) = FockSchmidt_SE(VirtStart:nSites,VirtStart:nSites)
            FockSchmidt_SE_CX(1:CoreEnd,ActiveStart:ActiveEnd) = FockSchmidt_SE(1:CoreEnd,ActiveStart:ActiveEnd)
            FockSchmidt_SE_XC(ActiveStart:ActiveEnd,1:CoreEnd) = FockSchmidt_SE(ActiveStart:ActiveEnd,1:CoreEnd)
            FockSchmidt_SE_VX(VirtStart:nSites,ActiveStart:ActiveEnd) = FockSchmidt_SE(VirtStart:nSites,ActiveStart:ActiveEnd)
            FockSchmidt_SE_XV(ActiveStart:ActiveEnd,VirtStart:nSites) = FockSchmidt_SE(ActiveStart:ActiveEnd,VirtStart:nSites)
            
            !Set the impurity:impurity parts to the correct values (even though we don't access them from FockSchmidt)
            !They are different since the correlation potential is not defined over the impurity sites.
            FockSchmidt_SE(nOcc-nImp+1:nOcc+nImp,nOcc-nImp+1:nOcc+nImp) = Emb_h0v_SE(:,:)

            deallocate(temp)

        else
            !The hamiltonians use for the interacting problem remain the same.
            !That is: FockSchmidt_SE = dcmplx(FockSchmidt)  (this should have been done at initialization, and remain the same)
            ! and     Emb_h0v_SE = Emb_h0v
            do i = 1,EmbSize
                do j = 1,EmbSize
                    Emb_h0v_SE(j,i) = dcmplx(Emb_h0v(j,i),0.0_dp)
                enddo
            enddo
        endif
        
        deallocate(FullSchmidtTrans_C,AO_OneE_Ham,W_Vals,HFPertBasis_Ann_Bra,HFPertBasis_Cre_Bra)
        deallocate(HFPertBasis_Ann_Ket,HFPertBasis_Cre_Ket,LVec,RVec)

        !call writevectorcomp(SchmidtPertGF_Cre_Ket(:,1),'SchmidtPertGF_Cre')
        !call writevectorcomp(SchmidtPertGF_Ann_Ket(:,1),'SchmidtPertGF_Ann')

    end subroutine FindNI_Charged

    !Find the non-interacting solution for alpha single-particle addition and removal, and project this operator into the schmidt basis
    !This is done seperately for electron addition and removal (SchmidtPertGF_Cre and SchmidtPertGF_Add)
    ! G_0^+
    ! This is ONLY done for the local greens functions at site 1
    subroutine FindSchmidtPert_Charged(Omega,NI_LRMat_Cre,NI_LRMat_Ann,tSwapSpins)
        use mat_tools, only: add_localpot
        implicit none
        real(dp), intent(in) :: Omega
        complex(dp), intent(out) :: NI_LRMat_Cre(1,1),NI_LRMat_Ann(1,1)
        logical, intent(in), optional :: tSwapSpins
        logical :: tSwapSpins_
        complex(dp), allocatable :: HFPertBasis_Cre(:),HFPertBasis_Ann(:)
        integer :: i,j,a,b,pertsite
        character(len=*), parameter :: t_r='FindSchmidtPert_Charged'

        !tSwapSpins = T means that the non-interacting excitations will be over the *beta* space, not the alpha space
        if(present(tSwapSpins)) then
            tSwapSpins_ = tSwapSpins
        else
            tSwapSpins_ = .false.
        endif

        if(tAllImp_LR) then
            call stop_all(t_r,'Should not be in this routine if you need to calculate all nImp GFs')
        endif

        allocate(HFPertBasis_Cre(nOcc+1:nSites))
        HFPertBasis_Cre(:) = zzero 
        allocate(HFPertBasis_Ann(1:nOcc))
        HFPertBasis_Ann(:) = zzero 
        NI_LRMat_Cre(1,1) = zzero 
        NI_LRMat_Ann(1,1) = zzero 

        !Assume perturbation is local to the first impurity site (pertsite = 1).
        pertsite = 1
        if(tSwapSpins_) then
            !Assuming that MF GF operator is V_ia/w-e_a, V_ia/w-e_i+w
            do i = 1,nOcc
                HFPertBasis_Ann(i) = dcmplx(HFOrbs_b(pertsite,i),0.0_dp)/(dcmplx(Omega,dDelta)-HFEnergies_b(i))
                NI_LRMat_Ann(1,1) = NI_LRMat_Ann(1,1) + dcmplx(HFOrbs_b(pertsite,i)**2,0.0_dp) /    &
                    (dcmplx(Omega,dDelta)-HFEnergies_b(i))
            enddo
            do a = nOcc+1,nSites
                HFPertBasis_Cre(a) = dcmplx(HFOrbs_b(pertsite,a),0.0_dp)/(dcmplx(Omega,dDelta)-HFEnergies_b(a))
                NI_LRMat_Cre(1,1) = NI_LRMat_Cre(1,1) + dcmplx(HFOrbs_b(pertsite,a)**2,0.0_dp) /    &
                    (dcmplx(Omega,dDelta)-HFEnergies_b(a))
            enddo
        else
            !Assuming that MF GF operator is V_ia/w-e_a, V_ia/w-e_i+w
            do i = 1,nOcc
                HFPertBasis_Ann(i) = dcmplx(HFOrbs(pertsite,i),0.0_dp)/(dcmplx(Omega,dDelta)-HFEnergies(i))
                NI_LRMat_Ann(1,1) = NI_LRMat_Ann(1,1) + dcmplx(HFOrbs(pertsite,i)**2,0.0_dp) /  &
                    (dcmplx(Omega,dDelta)-HFEnergies(i))
            enddo
            do a = nOcc+1,nSites
                HFPertBasis_Cre(a) = dcmplx(HFOrbs(pertsite,a),0.0_dp)/(dcmplx(Omega,dDelta)-HFEnergies(a))
                NI_LRMat_Cre(1,1) = NI_LRMat_Cre(1,1) + dcmplx(HFOrbs(pertsite,a)**2,0.0_dp) /  &
                    (dcmplx(Omega,dDelta)-HFEnergies(a))
            enddo
        endif

!        call writevectorcomp(HFPertBasis_Cre,'HFPertBasis_Cre')
        SchmidtPertGF_Cre_Ket(:,1) = zzero 
        SchmidtPertGF_Ann_Ket(:,1) = zzero 

        !Transform into schmidt basis
        if(tSwapSpins_) then
            do a = nOcc+nImp+1,nSites
                do b = nOcc+1,nSites
                    SchmidtPertGF_Cre_Ket(a,1) = SchmidtPertGF_Cre_Ket(a,1) + HFtoSchmidtTransform_b(b,a)*HFPertBasis_Cre(b)
                enddo
            enddo
            do i = 1,nOcc-nImp
                do j = 1,nOcc
                    SchmidtPertGF_Ann_Ket(i,1) = SchmidtPertGF_Ann_Ket(i,1) + HFtoSchmidtTransform_b(j,i)*HFPertBasis_Ann(j)
                enddo
            enddo
        else
            do a = nOcc+nImp+1,nSites
                do b = nOcc+1,nSites
                    SchmidtPertGF_Cre_Ket(a,1) = SchmidtPertGF_Cre_Ket(a,1) + HFtoSchmidtTransform(b,a)*HFPertBasis_Cre(b)
                enddo
            enddo
            do i = 1,nOcc-nImp
                do j = 1,nOcc
                    SchmidtPertGF_Ann_Ket(i,1) = SchmidtPertGF_Ann_Ket(i,1) + HFtoSchmidtTransform(j,i)*HFPertBasis_Ann(j)
                enddo
            enddo
        endif

        SchmidtPertGF_Cre_Bra(:,:) = dconjg(SchmidtPertGF_Cre_Ket(:,:))
        SchmidtPertGF_Ann_Bra(:,:) = dconjg(SchmidtPertGF_Ann_Ket(:,:))

!        call writevectorcomp(SchmidtPertGF_Ann_Ket(:,1),'SchmidtPertGF_Ann_Ket')
!        call writevectorcomp(SchmidtPertGF_Ann_Bra(:,1),'SchmidtPertGF_Ann_Bra')
        deallocate(HFPertBasis_Ann,HFPertBasis_Cre)

    end subroutine FindSchmidtPert_Charged


    !Find the non-interacting perturbation, and project this operator into the schmidt basis of phi^0 + its virtual space
    !This is for the density perturbation 
    subroutine FindSchmidtPert(tNonIntTest,Omega,ni_lr)
        implicit none
        logical, intent(in) :: tNonIntTest  !Test to return just the perturbation in the normal HF basis
        real(dp), intent(in) :: Omega
        complex(dp) , intent(out) :: ni_lr
        complex(dp), allocatable :: HFPertBasis(:,:),temp(:,:)
        complex(dp), allocatable :: C_HFtoSTrans(:,:)
        real(dp) :: EDiff
        integer :: a,i,j,pertsite
        character(len=*), parameter :: t_r='FindSchmidtBasis'

        !Assume (for the moment) that the perturbation operates exclusively at the first site
        pertsite = 1

        if(tNonIntTest) then
            !Screw everything up by ensuring that HFOrbs = FullHFOrbs (likewise for energy eigenvalues)
            HFOrbs(:,:) = FullHFOrbs(:,:)
            HFEnergies(:) = FullHFEnergies(:)
        endif

        allocate(HFPertBasis(nSites,nSites))
        HFPertBasis(:,:) = zzero 

        ni_lr = zzero 

        !Assume perturbation is local to the first impurity site (pertsite = 1).
        if(pertsite.ne.1) call stop_all(t_r,'Perturbation is not local to the impurity site')
        !Assuming that MF DD operator is V_ia/e_a-e_i-w + V_ia/e_a-e_i+w
        do i=1,nOcc
            do a=nOcc+1,nSites
                EDiff = HFEnergies(a)-HFEnergies(i)
                HFPertBasis(i,a) = dcmplx(HFOrbs(pertsite,i),0.0_dp)*dcmplx(HFOrbs(pertsite,a),0.0_dp)*   &
                    dcmplx(Lambda,0.0_dp)*(dcmplx(1.0_dp,0.0_dp) / &
                    (dcmplx(Omega,dDelta)-dcmplx(EDiff,0.0_dp)))
                HFPertBasis(a,i) = HFPertBasis(i,a)
                ni_lr = ni_lr + dcmplx((HFOrbs(pertsite,i)*HFOrbs(pertsite,a))**2,0.0_dp) / &
                    (dcmplx(Omega,dDelta)-dcmplx(EDiff,0.0_dp))
            enddo
        enddo
        ni_lr = ni_lr * 2.0_dp

        !call writematrixcomp(HFPertBasis,'Perturbation in HF basis',.true.)
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

        allocate(C_HFtoSTrans(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                C_HFtoSTrans(j,i) = dcmplx(HFtoSchmidtTransform(j,i),0.0_dp)
            enddo
        enddo
        !call R_C_Copy_2D(C_HFtoSTrans(:,:),HFtoSchmidtTransform(:,:),nSites,nSites)

        if(allocated(SchmidtPert)) deallocate(SchmidtPert)
        allocate(SchmidtPert(nSites,nSites))
        !call DGEMM('T','N',nSites,nSites,nSites,dcmplx(1.0_dp,0.0_dp),HFtoSchmidtTransform,nSites,HFPertBasis,nSites,0.0_dp,temp,nSites)
        call ZGEMM('T','N',nSites,nSites,nSites,zone,C_HFtoSTrans,nSites, &
            HFPertBasis,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites, &
            C_HFtoSTrans,nSites,zzero,SchmidtPert,nSites)
        deallocate(temp,HFPertBasis,C_HFtoSTrans)

        !SchmidtPert is now the perturbation in the schmidt basis
        !call writematrix(SchmidtPert,'Perturbation in schmidt basis',.true.)
        !call writematrixcomp(SchmidtPert,'Perturbation in schmidt basis',.true.)
        !call stop_all('sdg','sdf')

    end subroutine FindSchmidtPert

end module LinearResponse
