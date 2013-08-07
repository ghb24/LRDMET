module LRSolvers
    use const
    use globals
    use errors, only: stop_all,warning
    implicit none
    !Parameters for matrix-vector multiplications within iterative solvers
    complex(dp), pointer :: zDirMV_Mat(:,:)
    complex(dp), pointer :: zDirMV_Mat_cmprs(:)
    integer, pointer :: zDirMV_Mat_cmprs_inds(:)
    real(dp), allocatable :: Precond_Diag(:)

    contains
    
    subroutine WriteKVecHeader(iunit,KVal)
        implicit none
        integer, intent(in) :: iunit
        real(dp), intent(in) :: KVal(LatticeDim)
        integer :: i

        write(iunit,"(A,F8.4)",advance='no') '"k = ',KVal(1)
        do i = 2,LatticeDim
            write(iunit,"(A,F8.4)",advance='no') ', ',KVal(i)
        enddo
        write(iunit,"(A)") ' "'

    end subroutine WriteKVecHeader

    subroutine GetNextkVal(kPnt,tFinishedk)
        implicit none
        integer, intent(inout) :: kPnt
        logical, intent(out) :: tFinishedk

        tFinishedk = .false.
        if(kPnt.eq.0) then
            !First kpoint
            KIndex = 1
        else
            KIndex = KIndex + 1
            if((nKCalcs.ne.0).and.(KIndex.gt.nKCalcs)) then
                tFinishedk = .true.
            elseif(KIndex.gt.nKPnts) then
                tFinishedk = .true.
            endif
        endif

        !If not run through them all, find the next kpoint
        if(.not.tFinishedk) then
            if(allocated(KCalcs)) then
                kPnt = KCalcs(KIndex)
            elseif(nKCalcs.ne.0) then
                kPnt = nint(real(KIndex-1,dp)*real(nKPnts,dp)/real(nKCalcs,dp)) + 1
            else
                !Just go through all kpoints
                kPnt = kPnt + 1
            endif
        endif

    end subroutine GetNextkVal
    
    !Calculate n-electron hamiltonians
    subroutine Fill_N_Np1_Nm1b_FCIHam(nElec,Nmax_N,Nmax_Np1,Nmax_Nm1,NHam,Np1Ham,Nm1Ham,  &
            NHam_cmps,Np1Ham_cmps,Nm1Ham_cmps,NHam_inds,Np1Ham_inds,Nm1Ham_inds,tSwapExcits)
        use DetToolsData
        use DetTools, only : GetHElement_comp
        use solvers, only: StoreCompMat_comp
        use utils, only: get_free_unit
        implicit none
        integer, intent(in) :: nElec
        integer, intent(in) :: Nmax_N,Nmax_Np1,Nmax_Nm1
        complex(dp), intent(out), optional :: NHam(nFCIDet,nFCIDet)
        complex(dp), intent(out), optional :: Np1Ham(nNp1FCIDet,nNp1FCIDet)
        complex(dp), intent(out), optional :: Nm1Ham(nNm1bFCIDet,nNm1bFCIDet)
        complex(dp), intent(out), optional :: NHam_cmps(Nmax_N) 
        complex(dp), intent(out), optional :: Np1Ham_cmps(Nmax_Np1)
        complex(dp), intent(out), optional :: Nm1Ham_cmps(Nmax_Nm1)
        integer, intent(out), optional :: NHam_inds(Nmax_N)
        integer, intent(out), optional :: Np1Ham_inds(Nmax_Np1)
        integer, intent(out), optional :: Nm1Ham_inds(Nmax_Nm1)
        logical, intent(in), optional :: tSwapExcits
        real(dp) :: rtmp
        integer :: i,j,iSize,iunit_tmp
        logical :: tCmprsMats,tExist,tSwapExcits_
        character(len=*), parameter :: t_r='Fill_N_Np1_Nm1b_FCIHam'

        !If tSwapexcits = T, then we want to calculate the N+1_beta and N-1_alpha hamiltonians
        !Otherwise, calculate N+1_alpha and N-1_beta
        if(present(tSwapExcits)) then
            !Want Np1b for particle and Nm1 for hole if tBetaExcit, or vice versa otherwise
            !We have removed the distinction, and so will have to remember...
            tSwapExcits_ = tSwapExcits
        else
            tSwapExcits_ = .false.
        endif

        if(present(NHam)) then
            tCmprsMats = .false.
        else
            if(.not.present(Np1Ham_inds)) call stop_all(t_r,'Compressed matrices not present')
            tCmprsMats = .true.
        endif

        if(tCmprsMats) then
            NHam_cmps(:) = zzero
            Np1Ham_cmps(:) = zzero
            Nm1Ham_cmps(:) = zzero
            NHam_inds(:) = 0
            Np1Ham_inds(:) = 0
            Nm1Ham_inds(:) = 0

            if(tReadMats) then
                inquire(file='CompressHam_N',exist=texist)
                if(texist) then
                    write(6,*) "N-electron hamiltonian found"
                    iunit_tmp = get_free_unit()
                    open(iunit_tmp,file='CompressHam_N',status='old')
                    rewind(iunit_tmp)
                    read(iunit_tmp,*) isize
                    if(isize.gt.Nmax_N) call stop_all(t_r,'Error here')
                    do i = 1,iSize
                        read(iunit_tmp,*) rtmp,NHam_inds(i)
                        NHam_cmps(i) = dcmplx(rtmp,0.0_dp)
                    enddo
                    close(iunit_tmp)
                else
                    call StoreCompMat_comp(FCIDetList,nElec,nFCIDet,Nmax_N,NHam_cmps,NHam_inds,FCIBitList)
                endif
            else
                call StoreCompMat_comp(FCIDetList,nElec,nFCIDet,Nmax_N,NHam_cmps,NHam_inds,FCIBitList)
            endif
            if(tSwapExcits_) then
                call StoreCompMat_comp(Np1bFCIDetList,nElec+1,nNp1bFCIDet,Nmax_Np1,Np1Ham_cmps,Np1Ham_inds,Np1bBitList)
                call StoreCompMat_comp(Nm1FCIDetList,nElec-1,nNm1FCIDet,Nmax_Nm1,Nm1Ham_cmps,Nm1Ham_inds,Nm1BitList)
            else
                call StoreCompMat_comp(Np1FCIDetList,nElec+1,nNp1FCIDet,Nmax_Np1,Np1Ham_cmps,Np1Ham_inds,Np1BitList)
                call StoreCompMat_comp(Nm1bFCIDetList,nElec-1,nNm1bFCIDet,Nmax_Nm1,Nm1Ham_cmps,Nm1Ham_inds,Nm1bBitList)
            endif
        else
            NHam(:,:) = zzero
            Np1Ham(:,:) = zzero
            Nm1Ham(:,:) = zzero
            !Beware - these can be non-hermitian. Ensure that the indices are the right way around!
            do i = 1,nFCIDet
                do j = 1,nFCIDet
                    call GetHElement_comp(FCIDetList(:,i),FCIDetList(:,j),nElec,NHam(i,j),  &
                        ilutnI=FCIBitList(i),ilutnJ=FCIBitList(j))
                enddo
            enddo
            if(tSwapExcits_) then
                do i = 1,nNp1bFCIDet
                    do j = 1,nNp1bFCIDet
                        call GetHElement_comp(Np1bFCIDetList(:,i),Np1bFCIDetList(:,j),nElec+1,Np1Ham(i,j),    &
                            ilutnI=Np1bBitList(i),ilutnJ=Np1bBitList(j))
                    enddo
                enddo
                do i = 1,nNm1FCIDet
                    do j = 1,nNm1FCIDet
                        call GetHElement_comp(Nm1FCIDetList(:,i),Nm1FCIDetList(:,j),nElec-1,Nm1Ham(i,j), &
                            ilutnI=Nm1BitList(i),ilutnJ=Nm1BitList(j))
                    enddo
                enddo
            else
                do i = 1,nNp1FCIDet
                    do j = 1,nNp1FCIDet
                        call GetHElement_comp(Np1FCIDetList(:,i),Np1FCIDetList(:,j),nElec+1,Np1Ham(i,j),    &
                            ilutnI=Np1BitList(i),ilutnJ=Np1BitList(j))
                    enddo
                enddo
                do i = 1,nNm1bFCIDet
                    do j = 1,nNm1bFCIDet
                        call GetHElement_comp(Nm1bFCIDetList(:,i),Nm1bFCIDetList(:,j),nElec-1,Nm1Ham(i,j), &
                            ilutnI=Nm1bBitList(i),ilutnJ=Nm1bBitList(j))
                    enddo
                enddo
            endif
        endif

    end subroutine Fill_N_Np1_Nm1b_FCIHam

    !RHS is overwritten with the solution
    !LHS is destroyed
    subroutine SolveCompLinearSystem(LHS,RHS,nLinearSystem,info)
        use matrixops, only: z_inv 
        implicit none
        integer, intent(in) :: nLinearSystem
        integer, intent(out) :: info
        complex(dp), intent(inout) :: LHS(nLinearSystem,nLinearSystem)
        complex(dp), intent(inout) :: RHS(nLinearSystem)
        integer, allocatable :: Pivots(:)
        integer :: i,lWork
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: cWork(:),tempc(:,:),RVec(:,:),LVec(:,:),H_Valsc(:)
        complex(dp), allocatable :: CanTrans(:,:)
        character(len=*), parameter :: t_r='SolveCompLinearSystem'
    
        if(iSolveLR.eq.1) then
            !Use standard linear equation solver. Solution returned in RHS
            info = 0
            allocate(Pivots(nLinearSystem))
            call ZGESV(nLinearSystem,1,LHS,nLinearSystem,Pivots,RHS,nLinearSystem,info)
            deallocate(Pivots)
        elseif(iSolveLR.eq.2) then
            !Use advanced linear equation solver. Solution returned in RHS
            info = 0
            lWork = -1
            allocate(cWork(1))
            call ZGELS('N',nLinearSystem,nLinearSystem,1,LHS,nLinearSystem,RHS,nLinearSystem,cWork,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query for ZGELS failed')
            lwork = int(abs(cWork(1))) + 1
            deallocate(cWork)
            allocate(cWork(lWork))
            call ZGELS('N',nLinearSystem,nLinearSystem,1,LHS,nLinearSystem,RHS,nLinearSystem,cWork,lWork,info)
            deallocate(cWork)
        elseif(iSolveLR.eq.3) then
            !Solve linear equations via direct inversion

            allocate(tempc(nLinearSystem,nLinearSystem))    !The inverse of the hamiltonian system. May fail to invert at/close to poles where eigenvalues -> 0
            tempc(:,:) = dcmplx(0.0_dp,0.0_dp)
            !Directly invert the LHS
            call z_inv(LHS,tempc)

            !Multiply the inverse by the RHS of equations
            allocate(cWork(nLinearSystem))
            call ZGEMM('N','N',nLinearSystem,1,nLinearSystem,dcmplx(1.0_dp,0.0_dp),tempc,nLinearSystem,    &
                RHS,nLinearSystem,dcmplx(0.0_dp,0.0_dp),cWork,nLinearSystem)
            !Copy final solution to RHS
            RHS(:) = cWork(:)
            deallocate(cWork,tempc)
        elseif(iSolveLR.eq.4) then
            !Solve linear equations via complete diagonalization of hamiltonian
            !Beware, this matrix is not hermitian! At least the diagonals of the matrix are complex

            allocate(H_Valsc(nLinearSystem))
            allocate(RVec(nLinearSystem,nLinearSystem))
            allocate(LVec(nLinearSystem,nLinearSystem))
            RVec(:,:) = dcmplx(0.0_dp,0.0_dp)
            LVec(:,:) = dcmplx(0.0_dp,0.0_dp)
            H_Valsc(:) = dcmplx(0.0_dp,0.0_dp)
            allocate(Work(max(1,2*nLinearSystem)))
            allocate(cWork(1))
            lWork = -1
            info = 0
            call ZGEEV('V','V',nLinearSystem,LHS,nLinearSystem,H_Valsc,LVec,nLinearSystem,RVec,nLinearSystem,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Workspace query failed')
            lwork = int(cWork(1))+1
            deallocate(cWork)
            allocate(cWork(lwork))
            call ZGEEV('V','V',nLinearSystem,LHS,nLinearSystem,H_Valsc,LVec,nLinearSystem,RVec,nLinearSystem,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diag of LHS failed')
            deallocate(work,cWork)

            !Now, find inverse
            allocate(tempc(nLinearSystem,nLinearSystem)) !Temp array to build inverse in
            tempc(:,:) = dcmplx(0.0_dp,0.0_dp)
            do i=1,nLinearSystem
                if(abs(H_Valsc(i)).lt.1.0e-12_dp) then
                    write(6,*) "Eigenvalue: ",i,H_Valsc(i)
                    call warning(t_r,'VERY small/negative eigenvalue of LHS.')
                endif
                tempc(i,i) = 1.0_dp/H_Valsc(i)
            enddo
            !Rotate back into original basis
            allocate(CanTrans(nLinearSystem,nLinearSystem))
            call ZGEMM('N','N',nLinearSystem,nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),    &
                RVec,nLinearSystem,tempc,nLinearSystem,dcmplx(0.0_dp,0.0_dp),CanTrans,nLinearSystem)
            call ZGEMM('N','C',nLinearSystem,nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),    &
                CanTrans,nLinearSystem,RVec,nLinearSystem,dcmplx(0.0_dp,0.0_dp),tempc,nLinearSystem)
            deallocate(CanTrans)
            !tempc should now be the inverse
            !Multiply by RHS
            allocate(cWork(nLinearSystem))
            call ZGEMV('N',nLinearSystem,nLinearSystem,dcmplx(1.0_dp,0.0_dp),tempc,nLinearSystem,  &
                RHS,1,dcmplx(0.0_dp,0.0_dp),cWork,1)
!            call ZGEMM('N','N',nLinearSystem,1,nLinearSystem,dcmplx(1.0_dp,0.0_dp),    &
!                tempc,nLinearSystem,RHS,nLinearSystem,dcmplx(0.0_dp,0.0_dp),cWork,nLinearSystem)
            !Copy final solution to RHS
            RHS(:) = cWork(:)
            deallocate(cWork,tempc,H_Valsc,LVec,RVec)

        endif

    end subroutine SolveCompLinearSystem

    !Input: n (size of linear system)
    !       RHS (RHS of system)
    !       nout (unit)
    !       Maxiter (max iterations)
    !       rtol    (Convergence tolerance)
    !       tGuess  (Logical: Whether x is initial guess or not)
    !       tPrecond (Logical: Whether to precondition)
    !Output: iters (number of iterations required)
    !        gmres_info (0 if success)
    !InOut: x (solution / initial guess)
    subroutine GMRES_Solve(n,RHS,nout,Maxiter,rtol,tPrecond,tGuess,iters,x,error)
        implicit none
        integer, intent(in) :: n,nout,Maxiter
        real(dp) :: rtol
        logical, intent(in) :: tGuess,tPrecond
        complex(dp), intent(in) :: RHS(n)
        complex(dp), intent(inout) :: x(n)
        integer, intent(out) :: iters,error
        integer :: i,k,m
        integer :: lwork 
        integer :: icntl(8), irc(8), info(3)
        integer :: revcom, colx, coly, colz, nbscal, ierr
        complex(dp), allocatable :: work(:)
        real(dp) :: cntl(5),rinfo(2)
        integer, parameter :: finished = 0
        integer, parameter :: matvec = 1
        integer, parameter :: precondLeft = 2
        integer, parameter :: precondRight = 3
        integer, parameter :: dotprod = 4
        character(len=*), parameter :: t_r='GMRES_Solve'

        m = nKrylov  !This is the number of krylov vectors to store before restarting. Affects memory

        !Initialise control parameters to default values (note zero unit numbers means supress output)
        call init_zgmres(icntl,cntl)
        !icntl(1) = unit for error messages (6)
        !icntl(2) = unit for warnings (6)
        !icntl(3) = unit for convergence history (0)
        !icntl(4) = Location of preconditioning (0: None, 1: Left, 2: Right, 3: Double, 4: Error)
        !icntl(5) = What orthogonalization? (0: modified Gram-Schmidt (default), 1: Iterative modified GS, 2: Classical GS, 3: Iterative GS)
        !icntl(6) = Initial Guess? (1 = yes, 0 = no - use zero)
        !icntl(7) = Max interations
        !icntl(8) = How to compute residual at restart (1)

        !cntl(1) = Convergence tolerance of backward error (default = 1.E-5)
        !cntl(2:5) = normalizations (default 0)

        icntl(1) = nout
        icntl(2) = nout
        icntl(3) = nout
        icntl(7) = Maxiter
        
        if(tPrecond) then
            icntl(4) = 1    !Left preconditioning only
        else
            icntl(4) = 0
        endif

        cntl(1) = rtol

        !Determine the size of the work array
        if((icntl(5).eq.0).or.(icntl(5).eq.1)) then
            if(icntl(8).eq.1) then
                lwork = m*m + m*(n+5) + 5*n + 2
            else
                lwork = m*m + m*(n+5) + 6*n + 2
            endif
        else
            if(icntl(8).eq.1) then
                lwork = m*m + m*(n+5) + 5*n + m + 1
            else
                lwork = m*m + m*(n+5) + 6*n + m + 1
            endif
        endif

        write(6,"(A,F15.5,A)") "Memory required for GMRES solver: ",real(lwork,dp)*ComptoMb," Mb"
        call flush(6)

        allocate(work(lwork),stat=ierr)
        if(ierr.ne.0) then
            write(6,*) "Memory requirements too much. Allocation error for GMRES"
            write(6,*) "To reduce memory requirements, reduce NKRYLOV"
            call stop_all(t_r,'Alloc error')
        endif
        work(:) = zzero
        
        !Are we inputting a guess solution?
        if(tGuess) then
            icntl(6) = 1
            do i = 1,n
                work(i) = x(i)
            enddo
        else
            icntl(6) = 0
        endif

        !Set b
!        write(6,*) "RHS set to: "
        do i = 1,n
            work(i+n) = RHS(i)
!            write(6,*) RHS(i),work(i+n)
        enddo

        do while(.true.)

            call drive_zgmres(n,n,m,lwork,work,irc,icntl,cntl,info,rinfo)

            !What does the driver want us to do?
            !irc(1) determines what the gmres algorithm wants us to do
            !while irc(2:5) determines where or how to do it.
            revcom = irc(1)
            colx = irc(2)
            coly = irc(3)
            colz = irc(4)
            nbscal = irc(5)

!            write(6,*) "Return operation: ",revcom
!            call flush(6)

            if(revcom.eq.matvec) then
                !Do Matrix vector multiplicitation with vector work(colx:colx+n-1)
                !and put the result into work(colz:colz+n-1) 
                if(((colz.ge.colx).and.(colz.le.(colx+n-1))).or.    &
                    (((colz+n-1).ge.colx).and.((colz+n-1).le.(colx+n-1)))) then
                    !Check that they will not overwrite each other
                    write(6,*) "Vector indices to multiply: ",colx,'to',colx+n-1
                    write(6,*) "Output vector: ",colz,'to',colz+n-1
                    call stop_all(t_r,'vector indices overlapping')
                endif

                if(tCompressedMats) then
                    if(.not.associated(zDirMV_Mat_cmprs)) call stop_all(t_r,'Compressed matrix not associated')
!                    write(6,*) "Matrix multiplication..."
!                    call flush(6)
                    
                    !Sparse matrix multiply
                    if(zDirMV_Mat_cmprs_inds(1).ne.(n+2)) then
                        call stop_all(t_r,'Mismatched vector and matrix')
                    endif
                    do i = 1,n
                        work(colz+i-1) = zDirMV_Mat_cmprs(i)*work(colx+i-1)
                        do k = zDirMV_Mat_cmprs_inds(i),zDirMV_Mat_cmprs_inds(i+1)-1
                            work(colz+i-1) = work(colz+i-1) + zDirMV_Mat_cmprs(k)*work(colx+zDirMV_Mat_cmprs_inds(k)-1)
                        enddo
                    enddo
                else
                    if(.not.associated(zDirMV_Mat)) call stop_all(t_r,'Maxtrix not associated')
                    call ZGEMV('N',n,n,zone,zDirMV_Mat,n,work(colx),1,zzero,work(colz),1)
                endif

            elseif(revcom.eq.precondLeft) then
                !Perform left preconditioning, i.e.
                !work(colz:colz+n-1)  <--  M_1^{-1} * work(colx:colx+n-1)
                !work(colz:colz+n-1) = work(colx:colx+n-1)
                !Assume that preconditioner is just the diagonal (Jacobi preconditioning). Probably not optimal.

                if(tCompressedMats) then
                    do i = 1,n
                        if(abs(zDirMV_Mat_cmprs(i)).gt.1.0e-7_dp) then
                            work(colz+i-1) = work(colx+i-1)/zDirMV_Mat_cmprs(i)
                        else
                            work(colz+i-1) = work(colx+i-1)*1.0e7_dp
                        endif
                    enddo
                else
                    do i = 1,n
                        if(abs(zDirMV_Mat(i,i)).gt.1.0e-7_dp) then
                            work(colz+i-1) = work(colx+i-1)/zDirMV_Mat(i,i)
                        else
                            work(colz+i-1) = work(colx+i-1)*1.0e7_dp
                        endif
                    enddo
                endif

            elseif(revcom.eq.precondRight) then
                !Perform right preconditioning, i.e.
                !work(colz:colz+n-1)  <--  M_2^{-1} * work(colx:colx+n-1)
                !TODO: Implement properly!
                work(colz:colz+n-1) = work(colx:colx+n-1)
            elseif(revcom.eq.dotProd) then
                !Dot product with all previous vectors
                !work(colz:colz+n-1)  <--  work(colx:colx+n-1) * work(coly:coly+n-1)

                call zgemv('C',n,nbscal,zone,work(colx),n,work(coly),1,zzero,work(colz),1)
                !Same as
!                do i = 0,nbscal-1
!                    work(colz+i) = zdotc(n,work(colx+i*n),1,work(coly),1)
!                enddo
            elseif(revcom.eq.finished) then
                exit
            endif
        enddo

        !info(1) = 0 (Convergence reached, info(2) = number of iterations, rinfo(1) = preconditioned error, rinfo(2) = unpreconditioned error)
        !info(1) = -1 (n < 1)
        !info(1) = -2 (m < 1)
        !info(1) = -3 (lwork too small)
        !info(1) = -4 (Convergence not achieved in the number of iterations allowed)
        !info(1) = -5 (Preconditioning type not set)

        if(info(1).eq.-1) then
            call stop_all(t_r,'n < 1')
        elseif(info(1).eq.-2) then
            call stop_all(t_r,'m < 1')
        elseif(info(1).eq.-3) then
            write(6,*) "Minimal workspace required is: ",info(2)
            call stop_all(t_r,'lwork too small for current settings')
        elseif(info(1).eq.-4) then
            write(6,"(A,I10)") "Max iter: ",Maxiter
            write(6,"(A,2G20.10)") "Convergence: ",rinfo(:)
            call warning(t_r,'Convergence not achieved in number of iterations allowed')
        elseif(info(1).eq.-5) then
            call stop_all(t_r,'Preconditioning type not set')
        elseif(info(1).ne.0) then
            call stop_all(t_r,'Unspecified error')
        endif

        !Set return variables
!        write(6,*) "Final vector is: "
        do i = 1,n
            x(i) = work(i)
!            write(6,*) work(i)
        enddo
        iters = info(2)
        error = 0         

        !Final error given by rinfo(1)
        deallocate(work)

    end subroutine GMRES_Solve

    !For the actual diagonal part of the matrix which is being applied
    subroutine FormPrecond(n)
        implicit none
        integer, intent(in) :: n
!        real(dp) :: Scal
        integer :: i,k
        complex(dp) :: zdotc
        character(len=*), parameter :: t_r='FormPrecond'
        
        Precond_Diag(:) = 0.0_dp
    
        if(tCompressedMats.and.(.not.associated(zDirMV_Mat_cmprs))) call stop_all(t_r,'Compressed matrix not associated')

        if((.not.tCompressedMats).and.(.not.associated(zDirMV_Mat))) call stop_all(t_r,'Matrix not associated!')

        if(tCompressedMats) then

            do i = 1,n
                Precond_diag(i) = real(zDirMV_Mat_cmprs(i)*conjg(zDirMV_Mat_cmprs(i)),dp)
                !Now run through off-diagonals
                do k = zDirMV_Mat_cmprs_inds(i),zDirMV_Mat_cmprs_inds(i+1)-1
                    Precond_diag(i) = Precond_diag(i) + real(zDirMV_Mat_cmprs(zDirMV_Mat_cmprs_inds(k))* &
                        conjg(zDirMV_Mat_cmprs(zDirMV_Mat_cmprs_inds(k))),dp)
                enddo
            enddo

        else
!            Scal = real(zShift*dconjg(zShift))
    !The matrix diagonals are not long necessarily real
            do i = 1,n
                Precond_diag(i) = real(zdotc(n,zDirMV_Mat(:,i),1,zDirMV_Mat(:,i),1),dp)
!
!                tmp = zzero
!                do j = 1,n
!                    tmp = tmp + zDirMV_Mat(j,i)
!                enddo
!                tmp = tmp * dconjg(zShift)
!                Precond_diag(i) = Precond_diag(i) - 2.0_dp*real(tmp,dp) + Scal
            enddo
        endif
        !call writevector(Precond_diag,'Precond')
    end subroutine FormPrecond

    !Apply preconditioning, solving for y = A^(-1)x for given x 
    !Just assume that the preconditioner is the diagonal of the matrix moved just that it is positive definite
    subroutine zPreCond(n,x,y)
        use const
        integer(ip), intent(in) :: n
        complex(dp), intent(in) :: x(n)
        complex(dp), intent(out) :: y(n)
        integer :: i
!        real(dp) :: di
!        complex(dp) :: di
!        character(len=*), parameter :: t_r='zPreCond'

!        if(.not.associated(zDirMV_Mat)) call stop_all(t_r,'Matrix not associated!')

        do i = 1,n
            y(i) = x(i) / abs(Precond_diag(i))
!            di = abs(zDirMV_Mat(i,i)**2 - 2.0_dp*real(zShift,dp) + (zShift*dconjg(zShift)))
!            di = Precond_diag(i)
!            if(di.lt.0.0_dp) then
!                write(6,*) "Precond, i: ",i,Precond_diag(i)
!                call stop_all(t_r,'Preconditioned matrix should be positive definite')
!            endif
!            if(abs(di).gt.1.0e-8_dp) then
!                y(i) = x(i) / di
!            else
!                y(i) = x(i)  / abs(di)
!            endif
        enddo
!        write(6,*) "Precond: ",MinMatEl
!        write(6,*) "x: ",x

    end subroutine zPreCond
    
    !A direct matrix multiplication
    !Multiply twice, first by A-I(zShift), then by (A*-I(conjg(zShift)))
    !The shift is now included in the matrix
    subroutine zDirMV(n,x,y)
        use const
        integer(ip), intent(in) :: n
        complex(dp), intent(in) :: x(n)
        complex(dp), intent(out) :: y(n)
        integer :: i,k,j
        complex(dp) :: temp(n)
        character(len=*), parameter :: t_r='zDirMV'

        if(tCompressedMats.and.(.not.associated(zDirMV_Mat_cmprs))) call stop_all(t_r,'Compressed matrix not associated')
        if((.not.tCompressedMats).and.(.not.associated(zDirMV_Mat))) call stop_all(t_r,'Matrix not associated!')

        if(tCompressedMats) then
            !Sparse matrix multiply
            if(zDirMV_Mat_cmprs_inds(1).ne.(n+2)) then
                write(6,*) "zDirMV_Mat_cmprs_inds(1): ",zDirMV_Mat_cmprs_inds(1)
                write(6,*) "n+2: ",n+2
                call stop_all(t_r,'Mismatched vector and matrix')
            endif
            do i = 1,n
                temp(i) = zDirMV_Mat_cmprs(i)*x(i)
                do k = zDirMV_Mat_cmprs_inds(i),zDirMV_Mat_cmprs_inds(i+1)-1
                    temp(i) = temp(i) + zDirMV_Mat_cmprs(k)*x(zDirMV_Mat_cmprs_inds(k))
                enddo
            enddo

            !Now by transpose
            do i = 1,n  !Diagonals
                y(i) = conjg(zDirMV_Mat_cmprs(i))*temp(i)
            enddo
            do i = 1,n
                do k = zDirMV_Mat_cmprs_inds(i),zDirMV_Mat_cmprs_inds(i+1)-1
                    j = zDirMV_Mat_cmprs_inds(k)
                    y(j) = y(j) + conjg(zDirMV_Mat_cmprs(k))*temp(i)
                enddo
            enddo
        else
            call ZGEMV('N',n,n,zone,zDirMV_Mat,n,x,1,zzero,temp,1)
            call ZGEMV('C',n,n,zone,zDirMV_Mat,n,temp,1,zzero,y,1)
        endif

!        temp(:) = x(:)
!        call ZGEMV('N',n,n,zone,zDirMV_Mat,n,x,1,-zShift,temp,1)
!!        temp(:) = temp(:) - temp(:)*zShift
!        y(:) = temp(:)
!        call ZGEMV('C',n,n,zone,zDirMV_Mat,n,temp,1,-dconjg(zShift),y,1)
!!        y(:) = y(:) - y(:)*dconjg(zShift)

    end subroutine zDirMV
end module LRSolvers
