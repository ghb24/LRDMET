module Davidson
    use const
    use errors, only: stop_all,warning
    use mat_tools, only: WriteVector,WriteMatrix,WriteVectorComp,WriteMatrixComp,znrm2
    implicit none

    contains
    
    !Take a complex, hermitian matrix, and find the lowest eigenvalue and eigenvector
    !Via davidson diagonalization
    !If tStartingVec, then Vec(:) contains the initial state to use
    subroutine Comp_NonDir_Davidson(nSize,Mat,Val,Vec,tStartingVec,tol,max_iter,tLowestVal,niter)
        implicit none
        integer, intent(in) :: nSize
        complex(dp), intent(in) :: Mat(nSize,nSize)
        real(dp), intent(out) :: Val
        complex(dp), intent(inout) :: Vec(nSize)
        logical, intent(in) :: tStartingVec
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        logical, intent(in), optional :: tLowestVal     !Whether to converge to lowest or highest state
        integer, intent(out), optional :: niter     !The number of iterations it took to converge

        !Local variables for optional arguments
        integer :: max_iter_
        real(dp) :: tol_
        logical :: tLowestVal_
        integer :: niter_

        complex(dp), allocatable :: SubspaceVecs(:,:),HSubspace(:,:)
        complex(dp), allocatable :: CurrVec(:)
        complex(dp), allocatable, target :: SubspaceMat_1(:,:),SubspaceMat_2(:,:)
        complex(dp), pointer :: SubspaceMat(:,:)
        real(dp), allocatable :: Eigenvals(:)
        complex(dp), allocatable :: Eigenvecs(:,:)
        integer :: ierr,iter,i
        complex(dp) :: zdotc
        real(dp) :: dConv
        real(dp), parameter :: kappa = 0.25     !Tolerance for orthogonalization procedure
        character(len=*), parameter :: t_r='Comp_NonDir_Davidson'

        if(.not.present(tol)) then
            !Set tol default
            tol_=1.0e-12_dp
        else
            tol_ = tol
        endif
        if(.not.present(max_iter)) then
            !Set max iter default
            max_iter_=100
        else
            max_iter_ = max_iter
        endif
        if(.not.present(tLowestVal)) then
            tLowestVal_ = .true.   !Whether to compute the largest or smallest eigenvalue. Currently, it is always the lowest
        else
            tLowestVal_ = tLowestVal
        endif
            
!        allocate(Eigenvecs(nSize,nSize))
!        allocate(Eigenvals(nSize))
!        Eigenvecs(:,:) = Mat(:,:)
!        allocate(CurrVec(1))
!        lWork=-1
!        ierr=0
!        call dsyev('V','U',nSize,Eigenvecs,nSize,Eigenvals,CurrVec,lWork,ierr)
!        if(ierr.ne.0) call stop_all(t_r,'Workspace queiry failed')
!        lwork=int(CurrVec(1))+1
!        deallocate(CurrVec)
!        allocate(CurrVec(lwork))
!        call dsyev('V','U',nSize,Eigenvecs,nSize,Eigenvals,CurrVec,lWork,ierr)
!        if(ierr.ne.0) call stop_all(t_r,'Diag failed')
!        write(6,*) "Exact ground state: ",Eigenvals(1)
!        deallocate(CurrVec,Eigenvecs,Eigenvals)

        !Allocate memory for subspace vectors
        !Ridiculous overuse of memory. Should really restrict this to < max_iter, and restart if gets too much
        allocate(SubspaceVecs(nSize,max_iter_),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        SubspaceVecs(:,:) = zzero 

        !Memory for H * subspace vectors, as they are created
        allocate(HSubspace(nSize,max_iter_),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        HSubspace(:,:) = zzero

        !Memory for current vector
        allocate(CurrVec(nSize),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        CurrVec(:) = zzero

        SubspaceMat => null()

        if(.not.tStartingVec) then
            !Initialize current vector - random, or just 1 in first element...
            call init_vector_comp(nSize,CurrVec)
        else
            CurrVec(:) = Vec(:)
        endif

        do iter = 1,max_iter_

            !First, orthogonalize current vector against previous vectors, using a modified Gram-Schmidt
            !Returns normalized trial vector
            call ModGramSchmidt_comp(CurrVec,nSize,SubspaceVecs(:,1:iter-1),iter-1,kappa)

            SubspaceVecs(:,iter) = CurrVec(:)

            !Apply hamiltonian to subspace vector and store all results
            call ApplyMat_comp(Mat,SubspaceVecs(:,iter),HSubspace(:,iter),nSize)

            !Form subspace matrix
            if(.not.associated(SubspaceMat)) then
                !First iteration - subspace not allocated yet
                if(iter.ne.1) call stop_all(t_r,'Error in determining subspace size')
                allocate(SubspaceMat_1(iter,iter))
                SubspaceMat => SubspaceMat_1
            else
                !Swap subspaces around
                if(allocated(SubspaceMat_1)) then
                    allocate(SubspaceMat_2(iter,iter))
                    SubspaceMat_2(1:iter-1,1:iter-1) = SubspaceMat_1(1:iter-1,1:iter-1)
                    SubspaceMat => SubspaceMat_2
                    deallocate(SubspaceMat_1)
                else
                    allocate(SubspaceMat_1(iter,iter))
                    SubspaceMat_1(1:iter-1,1:iter-1) = SubspaceMat_2(1:iter-1,1:iter-1)
                    SubspaceMat => SubspaceMat_1
                    deallocate(SubspaceMat_2)
                endif
            endif
            !Compute final column coming from new subspace vector
            do i = 1,iter
                SubspaceMat(i,iter) = zdotc(nSize,SubspaceVecs(:,i),1,HSubspace(:,iter),1)
                SubspaceMat(iter,i) = dconjg(SubspaceMat(i,iter))
!                write(6,*) i,SubspaceMat(i,iter)
            enddo

            if(allocated(Eigenvals)) then
                deallocate(Eigenvals)
                deallocate(Eigenvecs)
            endif
            allocate(Eigenvals(iter))
            allocate(Eigenvecs(iter,iter))
            Eigenvecs(:,:) = SubspaceMat(:,:)

            !Now diagonalize the subspace, returning eigenvalues and associated vectors in Eigenvecs
            call DiagSubspaceMat_comp(Eigenvecs,iter,Eigenvals)

            !Update Vec by expanding the appropriate eigenvector back into the full space
            !These are current best estimates for the final eigenvalue and vector
            if(tLowestVal_) then
                !We are after the lowest value, expand this one
                call ZGEMV('N',nSize,iter,zone,SubspaceVecs(:,1:iter),nSize,Eigenvecs(:,1), &
                    1,zzero,Vec,1)
                Val = Eigenvals(1)
            else
                !We are after the lowest value, expand this one
                call ZGEMV('N',nSize,iter,zone,SubspaceVecs(:,1:iter),nSize,Eigenvecs(:,iter), &
                    1,zzero,Vec,1)
                Val = Eigenvals(iter)
            endif

            !We also want to compute the largest vec x H
            !Either do this by applying the matrix again to CurrVec (less memory, as don't then need to store
            !HSubspace vectors), or expand out into the H x subspace space. We do the latter here, though its not
            !necessarily better
            if(tLowestVal_) then
                !We are after the lowest value, expand this one
                call ZGEMV('N',nSize,iter,zone,HSubspace(:,1:iter),nSize,Eigenvecs(:,1), &
                    1,zzero,CurrVec,1)
            else
                !We are after the lowest value, expand this one
                call ZGEMV('N',nSize,iter,zone,HSubspace(:,1:iter),nSize,Eigenvecs(:,iter), &
                    1,zzero,CurrVec,1)
            endif

!            write(6,*) "Expanded Ritz eigenvector back into full space"

            !Calculate residual vector
            !r = H * CurrVec - val * CurrVec
            !Residual vector stored in CurrVec
            CurrVec(:) = CurrVec(:) - Val*Vec(:)

            !Find norm of residual
            dConv = znrm2(nSize,CurrVec,1)
!            write(6,*) "Residual norm :",iter,dConv,Val
            if(dConv.le.tol_) then
                !Praise the lord, convergence.
                !Final vector already stored in Vec, and Val is corresponding eigenvalue
                exit
            endif

            !Now, we need to solve for the next guess vector, t.
            !Simplest new vector, though not orthogonal to previous guess. Convergence not expected to be great.
            !t = [1/(Diag(H) - val x I)] r
            !Just pairwise multiplication
            do i = 1,nSize
                if(abs(real(Mat(i,i))-Val).gt.1.0e-8_dp) then
                    CurrVec(i) = CurrVec(i) / (real(Mat(i,i)) - Val)
                endif
            enddo

        enddo

        if(iter.gt.max_iter_) call stop_all(t_r,'Davidson routine failed to converge')

        niter_ = iter
        if(present(niter)) niter = niter_
            
        !Deallocate memory as appropriate
        deallocate(CurrVec,HSubspace,SubspaceVecs,Eigenvals)
        if(allocated(SubspaceMat_1)) deallocate(SubspaceMat_1)
        if(allocated(SubspaceMat_2)) deallocate(SubspaceMat_2)
        nullify(SubspaceMat)

    end subroutine Comp_NonDir_Davidson

    !Take a real, symmetric matrix, and find the lowest eigenvalue and eigenvector
    !Via davidson diagonalization
    !If tStartingVec, then Vec(:) contains the initial state to use
    subroutine Real_NonDir_Davidson(nSize,Val,Vec,tStartingVec,Nmax,Mat,CompressMat,IndexMat,tol,max_iter,tLowestVal,niter)
        implicit none
        integer, intent(in) :: nSize
        real(dp), intent(out) :: Val
        real(dp), intent(inout) :: Vec(nSize)
        logical, intent(in) :: tStartingVec
        integer, intent(in) :: Nmax
        real(dp), intent(in), optional :: Mat(nSize,nSize)
        real(dp), intent(in), optional :: CompressMat(Nmax)
        integer, intent(in), optional :: IndexMat(Nmax)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        logical, intent(in), optional :: tLowestVal     !Whether to converge to lowest or highest state
        integer, intent(out), optional :: niter     !The number of iterations it took to converge

        !Local variables for optional arguments
        integer :: max_iter_
        real(dp) :: tol_
        logical :: tLowestVal_
        integer :: niter_
        logical :: tCompressedMat

        real(dp), allocatable :: SubspaceVecs(:,:),HSubspace(:,:)
        real(dp), allocatable :: CurrVec(:)
        real(dp), allocatable, target :: SubspaceMat_1(:,:),SubspaceMat_2(:,:)
        real(dp), pointer :: SubspaceMat(:,:)
        real(dp), allocatable :: Eigenvals(:),Eigenvecs(:,:)
        integer :: ierr,iter,i
        real(dp) :: ddot,dConv,dnrm2,DiagEl
        real(dp), parameter :: kappa = 0.25     !Tolerance for orthogonalization procedure
        character(len=*), parameter :: t_r='Real_NonDir_Davidson'

!        write(6,*) "Entered non-direct davidson routine..."
!        call flush(6)

        if(present(CompressMat)) then
            tCompressedMat = .true.
        else
            tCompressedMat = .false.
        endif

        if((.not.present(Mat)).and.(.not.present(CompressMat))) then
            call stop_all(t_r,'Neither expanded nor compressed matrix found')
        elseif(present(Mat).and.present(CompressMat)) then
            call stop_all(t_r,'Both full and compressed matrices found')
        endif
        if(tCompressedMat.and.(.not.present(CompressMat))) then
            call stop_all(t_r,'Compressed matrix not found')
        endif
        if(tCompressedMat.and.(.not.present(IndexMat))) then
            call stop_all(t_r,'Compressed index vector not found')
        endif
        if(tCompressedMat.and.(Nmax.eq.1)) then
            call stop_all(t_r,'Compressed matrix asked for, but Nmax == 1')
        endif
        if(tCompressedMat.and.(present(Mat))) then
            call stop_all(t_r,'Compressed matrices asked for, but full matrix found')
        endif

        if(.not.present(tol)) then
            !Set tol default
            tol_=1.0e-9_dp
        else
            tol_ = tol
        endif
        if(.not.present(max_iter)) then
            !Set max iter default
            max_iter_=400
        else
            max_iter_ = max_iter
        endif
        if(.not.present(tLowestVal)) then
            tLowestVal_ = .true.   !Whether to compute the largest or smallest eigenvalue. Currently, it is always the lowest
        else
            tLowestVal_ = tLowestVal
        endif
            
!        allocate(Eigenvecs(nSize,nSize))
!        allocate(Eigenvals(nSize))
!        Eigenvecs(:,:) = Mat(:,:)
!        allocate(CurrVec(1))
!        lWork=-1
!        ierr=0
!        call dsyev('V','U',nSize,Eigenvecs,nSize,Eigenvals,CurrVec,lWork,ierr)
!        if(ierr.ne.0) call stop_all(t_r,'Workspace queiry failed')
!        lwork=int(CurrVec(1))+1
!        deallocate(CurrVec)
!        allocate(CurrVec(lwork))
!        call dsyev('V','U',nSize,Eigenvecs,nSize,Eigenvals,CurrVec,lWork,ierr)
!        if(ierr.ne.0) call stop_all(t_r,'Diag failed')
!        write(6,*) "Exact ground state: ",Eigenvals(1)
!        deallocate(CurrVec,Eigenvecs,Eigenvals)

!        write(6,*) "Allocating memory"
!        call flush(6)

        !Allocate memory for subspace vectors
        !Ridiculous overuse of memory. Should really restrict this to < max_iter, and restart if gets too much
        allocate(SubspaceVecs(nSize,max_iter_),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        SubspaceVecs(:,:) = 0.0_dp

        !Memory for H * subspace vectors, as they are created
        allocate(HSubspace(nSize,max_iter_),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        HSubspace(:,:) = 0.0_dp

        !Memory for current vector
        allocate(CurrVec(nSize),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        CurrVec(:) = 0.0_dp

        SubspaceMat => null()

!        write(6,*) "Initilizing vector"
!        call flush(6)

        if(.not.tStartingVec) then
            !Initialize current vector - random, or just 1 in first element...
            call init_vector_real(nSize,CurrVec)
        else
            CurrVec(:) = Vec(:)
        endif

!        call writevector(CurrVec,'starting vector')

        do iter = 1,max_iter_

            if(mod(iter,10).eq.0) then
                write(6,*) "Starting Davidson iteration: ",iter,dConv
                call flush(6)
            endif

            !First, orthogonalize current vector against previous vectors, using a modified Gram-Schmidt
            !Returns normalized trial vector
            call ModGramSchmidt_real(CurrVec,nSize,SubspaceVecs(:,1:iter-1),iter-1,kappa)

!            write(6,*) "Orthogonalized subspace vector"

            SubspaceVecs(:,iter) = CurrVec(:)

            !Apply hamiltonian to subspace vector and store all results
            if(tCompressedMat) then
                call ApplyMat_real(SubspaceVecs(:,iter),HSubspace(:,iter),nSize,tCompressedMat,Nmax,    &
                    CompressMat=CompressMat,IndexMat=IndexMat)
            else
                call ApplyMat_real(SubspaceVecs(:,iter),HSubspace(:,iter),nSize,tCompressedMat,Nmax,Mat=Mat)
            endif

!            write(6,*) "Applied hamiltonian"
!            call writevector(HSubspace(:,iter),'H x subspace vector')

            !Form subspace matrix
            if(.not.associated(SubspaceMat)) then
                !First iteration - subspace not allocated yet
                if(iter.ne.1) call stop_all(t_r,'Error in determining subspace size')
                allocate(SubspaceMat_1(iter,iter))
                SubspaceMat => SubspaceMat_1
            else
                !Swap subspaces around
                if(allocated(SubspaceMat_1)) then
                    allocate(SubspaceMat_2(iter,iter))
                    SubspaceMat_2(1:iter-1,1:iter-1) = SubspaceMat_1(1:iter-1,1:iter-1)
                    SubspaceMat => SubspaceMat_2
                    deallocate(SubspaceMat_1)
                else
                    allocate(SubspaceMat_1(iter,iter))
                    SubspaceMat_1(1:iter-1,1:iter-1) = SubspaceMat_2(1:iter-1,1:iter-1)
                    SubspaceMat => SubspaceMat_1
                    deallocate(SubspaceMat_2)
                endif
            endif
            !Compute final column coming from new subspace vector
            do i = 1,iter
                SubspaceMat(i,iter) = ddot(nSize,SubspaceVecs(:,i),1,HSubspace(:,iter),1)
                SubspaceMat(iter,i) = SubspaceMat(i,iter)
!                write(6,*) i,SubspaceMat(i,iter)
            enddo

!            write(6,*) "Updated subspace hamiltonian"
            
            if(allocated(Eigenvals)) then
                deallocate(Eigenvals)
                deallocate(Eigenvecs)
            endif
            allocate(Eigenvals(iter))
            allocate(Eigenvecs(iter,iter))
            Eigenvecs(:,:) = SubspaceMat(:,:)

!            call writematrix(Eigenvecs,'subspace matrix',.true.)

            !Now diagonalize the subspace, returning eigenvalues and associated vectors in Eigenvecs
            call DiagSubspaceMat_real(Eigenvecs,iter,Eigenvals)

!            write(6,*) "Diagonalized subspace hamiltonian"

            !Update Vec by expanding the appropriate eigenvector back into the full space
            !These are current best estimates for the final eigenvalue and vector
            if(tLowestVal_) then
                !We are after the lowest value, expand this one
                call DGEMV('N',nSize,iter,1.0_dp,SubspaceVecs(:,1:iter),nSize,Eigenvecs(:,1), &
                    1,0.0_dp,Vec,1)
                Val = Eigenvals(1)
            else
                call DGEMV('N',nSize,iter,1.0_dp,SubspaceVecs(:,1:iter),nSize,Eigenvecs(:,iter), &
                    1,0.0_dp,Vec,1)
                Val = Eigenvals(iter)
            endif

            !We also want to compute the largest vec x H
            !Either do this by applying the matrix again to CurrVec (less memory, as don't then need to store
            !HSubspace vectors), or expand out into the H x subspace space. We do the latter here, though its not
            !necessarily better
            if(tLowestVal_) then
                !We are after the lowest value, expand this one
                call DGEMV('N',nSize,iter,1.0_dp,HSubspace(:,1:iter),nSize,Eigenvecs(:,1), &
                    1,0.0_dp,CurrVec,1)
            else
                call DGEMV('N',nSize,iter,1.0_dp,HSubspace(:,1:iter),nSize,Eigenvecs(:,iter), &
                    1,0.0_dp,CurrVec,1)
            endif

!            write(6,*) "Expanded Ritz eigenvector back into full space"

            !Calculate residual vector
            !r = H * CurrVec - val * CurrVec
            !Residual vector stored in CurrVec
            CurrVec(:) = CurrVec(:) - Val*Vec(:)

            !Find norm of residual
            dConv = dnrm2(nSize,CurrVec,1)
!            write(6,*) "Residual norm :",iter,dConv,Val
            if(dConv.le.tol_) then
                !Praise the lord, convergence.
                !Final vector already stored in Vec, and Val is corresponding eigenvalue
                exit
            endif

            !Now, we need to solve for the next guess vector, t.
            !Simplest new vector, though not orthogonal to previous guess. Convergence not expected to be great.
            !t = [1/(Diag(H) - val x I)] r
            !Just pairwise multiplication
            do i = 1,nSize
                if(tCompressedMat) then
                    !Diagonals are first elements of compressed matrix
                    DiagEl = CompressMat(i)
                else
                    DiagEl = Mat(i,i)
                endif
    
                if(abs(DiagEl-Val).gt.1.0e-8_dp) then
                    CurrVec(i) = CurrVec(i) / (DiagEl - Val)
                endif
            enddo

        enddo

        if(iter.gt.max_iter_) call stop_all(t_r,'Davidson routine failed to converge')

        niter_ = iter
        if(present(niter)) niter = niter_
            
        !Deallocate memory as appropriate
        deallocate(CurrVec,HSubspace,SubspaceVecs,Eigenvals)
        if(allocated(SubspaceMat_1)) deallocate(SubspaceMat_1)
        if(allocated(SubspaceMat_2)) deallocate(SubspaceMat_2)
        nullify(SubspaceMat)

    end subroutine Real_NonDir_Davidson

    subroutine init_vector_real(nSize,CurrVec)
        implicit none
        integer, intent(in) :: nSize
        real(dp), intent(out) :: CurrVec(nSize)

        CurrVec(:) = 0.0_dp
!        CurrVec(1) = 1.0_dp
        call Random_number(CurrVec)

    end subroutine init_vector_real
            
    subroutine init_vector_comp(nSize,CurrVec)
        implicit none
        integer, intent(in) :: nSize
        complex(dp), intent(out) :: CurrVec(nSize)
        real(dp) :: Cmpt(2)
        integer :: i

        CurrVec(:) = zzero
!        CurrVec(1) = 1.0_dp
        do i = 1,nSize
            call Random_number(Cmpt)
            CurrVec(i) = dcmplx(Cmpt(1),Cmpt(2))
        enddo

    end subroutine init_vector_comp
            
    !Diagonalize the subspace hamiltonian. Essentially just a wrapper for DSYEV, but keeps the main loop clean
    subroutine DiagSubspaceMat_real(Mat,iSize,W)
        implicit none
        integer, intent(in) :: iSize
        real(dp), intent(inout) :: Mat(iSize,iSize)
        real(dp), intent(out) :: W(iSize)

        real(dp), allocatable :: Work(:)
        integer :: lWork,info,ierr
        character(len=*), parameter :: t_r='DiagSubspaceMat_real'

        allocate(Work(1),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        W(:)=0.0_dp
        lWork=-1
        info=0
        call dsyev('V','U',iSize,Mat,iSize,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        call dsyev('V','U',iSize,Mat,iSize,W,Work,lWork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)

    end subroutine DiagSubspaceMat_real

    !Diagonalize the subspace hamiltonian. Essentially just a wrapper for ZHEEV, but keeps the main loop clean
    subroutine DiagSubspaceMat_comp(Mat,iSize,W)
        implicit none
        integer, intent(in) :: iSize
        complex(dp), intent(inout) :: Mat(iSize,iSize)
        real(dp), intent(out) :: W(iSize)

        complex(dp), allocatable :: Work(:)
        real(dp), allocatable :: RWork(:)
        integer :: lWork,info,ierr
        character(len=*), parameter :: t_r='DiagSubspaceMat_real'

        allocate(Work(1),stat=ierr)
        allocate(RWork(max(1,3*iSize-2)),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        W(:)=0.0_dp
        lWork=-1
        info=0
        call zheev('V','U',iSize,Mat,iSize,W,Work,lWork,RWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        call zheev('V','U',iSize,Mat,iSize,W,Work,lWork,RWork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work,RWork)

    end subroutine DiagSubspaceMat_comp

    !Apply the matrix to a trial vector. Essentially just a wrapper for DGEMV, but will eventually
    !want to be modified for a direct algorithm
    subroutine ApplyMat_real(Vec,ResVec,nSize,tCompressedMat,Nmax,Mat,CompressMat,IndexMat)
        implicit none
        integer, intent(in) :: nSize
        real(dp), intent(in) :: Vec(nSize)
        real(dp), intent(out) :: ResVec(nSize)
        logical, intent(in) :: tCompressedMat
        integer, intent(in) :: Nmax
        real(dp), intent(in), optional :: Mat(nSize,nSize)
        real(dp), intent(in), optional :: CompressMat(Nmax)
        integer, intent(in), optional :: IndexMat(Nmax)
        character(len=*), parameter :: t_r='ApplyMat_real'
        integer :: i,k

        if(.not.tCompressedMat) then
            if(.not.present(Mat)) call stop_all(t_r,"Uncompressed matrix not passed through")

            call DGEMV('N',nSize,nSize,1.0_dp,Mat,nSize,Vec,1,0.0_dp,ResVec,1)
        else
            if(present(Mat)) call stop_all(t_r,'Uncompressed matrix passed through')
            if(.not.present(CompressMat)) call stop_all (t_r,'Compressed matrix not passed through')
            if(.not.present(IndexMat)) call stop_all(t_r,'Compressed index array not passed through')

            !Sparse matrix multiply
            if(IndexMat(1).ne.nSize+2) call stop_all(t_r,'Mismatched vector and matrix')

            do i = 1, nSize
                ResVec(i) = CompressMat(i)*Vec(i)   !Diagonal terms
                do k = IndexMat(i),IndexMat(i+1)-1  !Off diagonal terms
                    ResVec(i) = ResVec(i) + CompressMat(k)*Vec(IndexMat(k))
                enddo
            enddo
        endif

    end subroutine ApplyMat_real
    
    !Apply the matrix to a trial vector. Essentially just a wrapper for ZGEMV, but will eventually
    !want to be modified for a direct algorithm
    subroutine ApplyMat_comp(Mat,Vec,ResVec,nSize)
        implicit none
        integer, intent(in) :: nSize
        complex(dp), intent(in) :: Mat(nSize,nSize),Vec(nSize)
        complex(dp), intent(out) :: ResVec(nSize)
        
        call ZGEMV('N',nSize,nSize,zone,Mat,nSize,Vec,1,zzero,ResVec,1)

    end subroutine ApplyMat_comp

    !Modified Gram-Schmidt orthogonalization with refinement.
    !Kappa is introduced, which ensures that the loss of orthogonality is restricted to 1/kappa
    !times machine precision. A reasonable kappa is about 0.25
    !Vec is the vector to be orthogonalized (of size nSize)
    !Subspace are the set of nSubspace vectors to be orthogonalized against
    !Returns the orthogonalized and normalized vector
    subroutine ModGramSchmidt_real(Vec,nSize,Subspace,nSubspace,kappa)
        implicit none
        integer, intent(in) :: nSize,nSubspace
        real(dp), intent(inout) :: Vec(nSize)
        real(dp), intent(in) :: Subspace(nSize,nSubspace)
        real(dp), intent(in) :: kappa

        real(dp) :: Norm_in,Norm_out,Overlap,ddot,dnrm2
        integer :: i

!        call writevector(Vec,'Vec into GS')

        !Are there even any vectors to orthogonalize against!
        if(nSubspace.eq.0) then
            Norm_out = dnrm2(nSize,Vec,1)
            Vec(:) = Vec(:)/Norm_out
            return
        endif

        !Find initial normalization of vector
        Norm_in = dnrm2(nSize,Vec,1)

        do i = 1,nSubspace
            !Project out each component of the subspace
            Overlap = ddot(nSize,Vec,1,Subspace(:,i),1)
            Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
        enddo

        !Find new normalization
        Norm_out = dnrm2(nSize,Vec,1)

        if((Norm_out/Norm_in).le.kappa) then
            !Orthogonalize again!
            do i = 1,nSubspace 
                Overlap = ddot(nSize,Vec,1,Subspace(:,i),1)
                Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
            enddo
            Norm_out = dnrm2(nSize,Vec,1)
        endif

        Vec(:) = Vec(:)/Norm_out    !Normalize

    end subroutine ModGramSchmidt_real

    !Modified Gram-Schmidt orthogonalization with refinement for complex numbers.
    !See ModGramSchmidt_real for details.
    subroutine ModGramSchmidt_comp(Vec,nSize,Subspace,nSubspace,kappa)
        implicit none
        integer, intent(in) :: nSize,nSubspace
        complex(dp), intent(inout) :: Vec(nSize)
        complex(dp), intent(in) :: Subspace(nSize,nSubspace)
        real(dp), intent(in) :: kappa

        complex(dp) :: Overlap,zdotc
        integer :: i
        real(dp) :: Norm_out,Norm_in
!        character(len=*), parameter :: t_r='ModGramSchmidt_comp'

        !Are there even any vectors to orthogonalize against!
        if(nSubspace.eq.0) then
            Norm_out = znrm2(nSize,Vec,1)
            Vec(:) = Vec(:)/Norm_out
            return
        endif

        !Find initial normalization of vector
        Norm_in = znrm2(nSize,Vec,1)

        do i = 1,nSubspace
            !Project out each component of the subspace
            Overlap = zdotc(nSize,Subspace(:,i),1,Vec,1)
            Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
        enddo

        !Find new normalization
        Norm_out = znrm2(nSize,Vec,1)

        if((Norm_out/Norm_in).le.kappa) then
            !Orthogonalize again!
            do i = 1,nSubspace 
                Overlap = zdotc(nSize,Subspace(:,i),1,Vec,1)
                Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
            enddo
            Norm_out = znrm2(nSize,Vec,1)
        endif

        Vec(:) = Vec(:)/Norm_out

    end subroutine ModGramSchmidt_comp
end module Davidson
