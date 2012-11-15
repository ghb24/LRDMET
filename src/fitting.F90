module fitting
    use const
    use errors, only: stop_all
    use globals
    use mat_tools, only: RDMErr,FromTriangularPacked,ToTriangularPacked
    implicit none

    contains

    !Fit the correlation potential so that the RDMs match. This is returned in the global vloc_change, as well as a
    !measure of the change in the potential (VarVloc) and the initial error in the RDMs (ErrRDM)
    !The two RDMs that want to match are the FCI RDM over the embedded systems (HL_1RDM), and an RDM
    !that is produced at the mean field level by diagonalizing a fock matrix in the embedding basis (Constructed).
    subroutine Fit_vloc(VarVloc,ErrRDM)
        implicit none
        real(dp), intent(out) :: VarVloc,ErrRDM
        real(dp) , allocatable :: work(:)
        real(dp) , allocatable :: EmbNatOrbs(:,:)  !The eigenvectors from the MF+corr pot over embedded system
        real(dp) , allocatable :: HLNatOrbs(:,:),HL_OccNumbers(:)  !The eigensystem from HL_1RDM over embedded system
        real(dp) :: vloc_change_sq(nImp,nImp),DiffRDM_packed(nImpCombs),DiffRDM_unpacked(nImp,nImp)
        real(dp) , allocatable :: vloc_change_packed(:)
        integer :: lWork,info,i
        character(len=*), parameter :: t_r='Fit_vloc'

        if(allocated(vloc_change)) deallocate(vloc_change)
        allocate(vloc_change(nImp,nImp))

        !Allocate memory for the v_loc fit.
        allocate(vloc_change_packed(nImpCombs))

        !Diagonalise the DM in the embedded basis
        allocate(EmbNatOrbs(EmbSize,EmbSize))
        if(allocated(MFEmbOccs)) deallocate(MFEmbOccs)
        allocate(MFEmbOccs(EmbSize))
        EmbNatOrbs(:,:) = -1.0_dp * Emb_MF_DM(:,:)
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',EmbSize,EmbNatOrbs,EmbSize,MFEmbOccs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',EmbSize,EmbNatOrbs,EmbSize,MFEmbOccs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        MFEmbOccs(:) = -1.0_dp * MFEmbOccs(:)

        do i=1,EmbSize
            if(MFEmbOccs(i).lt.-1.0e-8) call stop_all(t_r,'Negative occupation number in MF embedding space')
            if(MFEmbOccs(i).lt.0.0_dp) MFEmbOccs(i) = 0.0_dp
        enddo

        !Now find occupation numbers over embedded system with FCI density matrix
        allocate(HLNatOrbs(EmbSize,EmbSize))
        allocate(HL_OccNumbers(EmbSize))

        HLNatOrbs(:,:) = HL_1RDM(:,:)
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',EmbSize,HLNatOrbs,EmbSize,HL_OccNumbers,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',EmbSize,HLNatOrbs,EmbSize,HL_OccNumbers,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        write(6,"(A)") "Embedded system occupation numbers (CoreH + Corr Pot , High level calc): "
        do i=1,EmbSize
            write(6,"(2F12.7)") MFEmbOccs(EmbSize-(i-1)),HL_OccNumbers(EmbSize-(i-1))
        enddo

        !The aim now, is to find a vloc (over the impurity sites) which when added to the fock matrix, and the occupation 
        !numbers of the original vloc+MF rotated into the embedded basis of the new vloc, will give the same RDM as the
        !high-level calculation in the embedded system.
        !Initial guess of vloc over impurity sites
        vloc_change_packed(:) = 0.0_dp
        !Newton-raphson fit.
        !Send in initial change in vloc_change (returned as optimal change fit in packed triangular form)
        call NR_opt(vloc_change_packed)
        call FromTriangularPacked(nImp,vloc_change_packed,vloc_change) !unpack

        !We now have the optimal change in vloc to match the rdms.
        !ErrRDM is a measure of the initial error in the RDMs (just over the impurity sites)
        DiffRDM_unpacked = HL_1RDM(1:nImp,1:nImp) - Emb_MF_DM(1:nImp,1:nImp)
        call ToTriangularPacked(nImp,DiffRDM_unpacked,DiffRDM_packed)
        ErrRDM = sum(DiffRDM_packed(:)**2.0_dp)

        !Change in potential - convergence metric: Tr[delta_v^T delta_v]
        call DGEMM('T','N',nImp,nImp,nImp,1.0_dp,vloc_change,nImp,vloc_change,nImp,0.0_dp,vloc_change_sq,nImp)
        VarVloc = 0.0_dp
        do i=1,nImp
            VarVloc = VarVloc + vloc_change_sq(i,i)
        enddo

        deallocate(EmbNatOrbs,HLNatOrbs,HL_OccNumbers)
        deallocate(vloc_change_packed)

    end subroutine Fit_vloc

    !We have a function of the nImpComb v_loc_change variables, which returns a residual over EmbCombs parameters, 
    !which we want to reduce so that their summed squares is zero.
    !vloc_change is updated, which is a global variable of the correlation potential to fit.
    subroutine NR_opt(x0)
        implicit none
        real(dp), intent(inout) :: x0(nImpCombs)    !Initial guess for potential 
        real(dp) :: x(nImpCombs), r(EmbCombs)
        real(dp) :: g(EmbCombs,nImpCombs)
        real(dp) :: g2(EmbCombs+nImpCombs,nImpCombs)    !Gradient (with added stuff)
        real(dp) :: r2(EmbCombs+nImpCombs)              !Fit residuals (with added stuff)
        real(dp) :: Sing(nImpCombs),U(EmbCombs+nImpCombs,nImpCombs),VT(nImpCombs,nImpCombs)
        real(dp) :: InvDiag(nImpCombs,nImpCombs),temp(nImpCombs,nImpCombs),temp2(nImpCombs,EmbCombs+nImpCombs)
        real(dp) :: dx(nImpCombs),x_temp(nImpCombs),r_temp(EmbCombs)
        real(dp), allocatable :: work(:)
        integer :: nr,nx,lWork,info,i,it
        real(dp) :: err,err_temp,step,Opt_step,Min_val,Searchstep,dstep,f1,f2,f3,norm,s1,s2,s3
        real(dp) :: NoMoveVal,LargestDiff
        character(len=*), parameter :: t_r='NR_opt'

        nr = EmbCombs   !Number of residuals = number of triangular packed embedding sites
        nx = nImpCombs  !Number of variables = number of triangular packed impurity sites
        x(:) = x0(:)    !Starting guess for potential

        step = 3.0e-3_dp    !Step size
        do it=1,100     !NR iterations

!            call writevector(x,'inputvars')
            call RDMErr(x,r) !Update residuals (r)
            err = sum(r(:)**2)  !Error metric

!            call writevector(r,'residuals')
!            write(6,*) "Fitting iteration: ",it,err

            if(err.lt.1.0e-15_dp) exit  !Convergence satisfied
            call MakeGradMatrix(x,g,nr,nx)    !Returns the numerically calculated gradient matrix in g from the potential x
            !g is now the jacobian

            !Now get the appropriate direction to move x in
            !Normally, we would just get the solution to the equation r == g * dx for the update to x
            !However, we have some extra padding here, with a error-dependent damping for the step length.
            !I.e. the smaller the error, the more the damping to restrict step length
            g2(:,:) = 0.0_dp
            g2(1:nr,1:nx) = g(:,:)
            do i=1,nx
                g2(nr+1+(i-1),i) = 0.1*sqrt(err)
            enddo
            r2(:) = 0.0_dp
            r2(1:nr) = r(:)

            !Least fitting for overcomplete specification: SVD.
            !Ax=b where A=g2 and b = r2
            !A = U D V^T from SVD
            !x = V (1/D) U^T b
            allocate(Work(1))
            lWork=-1
            info=0
            call DGESVD('S','S',nr+nx,nx,g2,nr+nx,Sing,U,nr+nx,VT,nx,work,lwork,info)
            if(info.ne.0) call stop_all(t_r,'SVD Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call DGESVD('S','S',nr+nx,nx,g2,nr+nx,Sing,U,nr+nx,VT,nx,work,lwork,info)
            if(info.ne.0) call stop_all(t_r,'SVD failed')
            deallocate(work)

            !Now calculate x as V (1/D) U^T b
            !Calculate the matrix representation of 1/D
            InvDiag(:,:) = 0.0_dp
            do i=1,nx
                InvDiag(i,i) = 1.0_dp/Sing(i)
            enddo
            call DGEMM('T','N',nx,nx,nx,1.0_dp,VT,nx,InvDiag,nx,0.0_dp,temp,nx)
            call DGEMM('N','T',nx,nr+nx,nx,1.0_dp,temp,nx,U,nr+nx,0.0_dp,temp2,nx)
            call DGEMM('N','N',nx,1,nr+nx,1.0_dp,temp2,nx,r2,nr+nx,0.0_dp,dx,nx)
            !dx is now the direction to move in to improve the potential

            !However, we also want to find the correct distance to move in
            !Calculate the optimal step size via a crude line search
            !We could just set this stepsize to 1, but better (but more expensive) would be to search along its length for the best
            Searchstep = 0.2_dp
            step = 0.0_dp           !Initial step attempt
            Min_val = 1.0e15_dp     !minimum value of error metric
            LargestDiff = 0.0_dp
            Opt_step = step       !Optimal step size for x
            do i=1,11
                !Calculate over relatively coarse grid
                !Calculate error function
                x_temp(:) = x(:) - (step * dx(:))   
                call RDMErr(x_temp,r_temp) !Update residuals (r)
                err_temp = sum(r_temp(:)**2)  !Error metric
                if(i.eq.1) then
                    NoMoveVal = err_temp
                else
                    if(abs(NoMoveVal-err_temp).gt.abs(LargestDiff)) LargestDiff = NoMoveVal-err_temp
                endif
!                write(6,*) "***",i,err_temp
                if(err_temp.lt.Min_val) then
                    !Better value - take it
                    Opt_step = step
                    Min_val = err_temp
                endif
                step = step + Searchstep
            enddo
!            write(6,*) "OptStep, MinVal: ",Opt_Step,Min_Val
            if(abs(LargestDiff).lt.1.0e-15) then
                !We're not moving. 
                exit
            endif
            if(Opt_step.lt.0.1) then
                !i.e. the optimal stepsize is between 0 and 0.2.
                !Now do a finer bisection to nail it down completely
                !Assume that we are bracketed by 0.0 and 0.2
                write(6,*) "Entering bisection..."
                s1 = 0.0_dp
                s2 = 0.2_dp
                f1 = Min_val
                x_temp(:) = x(:) - (s2 * dx(:))   
                call RDMErr(x_temp,r_temp) !Update residuals (r)
                f2 = sum(r_temp(:)**2)  !Error metric
                dstep = 0.2 
                do while(dstep.gt.1.0e-10)
                    s3 = (s1 + s2)/2.0_dp
                    x_temp(:) = x(:) - (s3 * dx(:))   
                    call RDMErr(x_temp,r_temp) !Update residuals (r)
                    f3 = sum(r_temp(:)**2)  !Error metric
                    if(f3.gt.max(f1,f2)) then
                        write(6,"(A)") "WARNING: Bisection not bracketed..."
                        dstep = 0.2_dp
                        exit
                    endif

                    !Replace the largest f with one from s3
                    if(f1.gt.f2) then
                        !Replace f1
                        f1 = f3
                        s1 = s3
                    else
                        !Replace f2
                        f2 = f3
                        s2 = s3
                    endif
                    dstep = s1-s2
                    if(s2.lt.s1) call stop_all(t_r,'s1 and s2 have swapped sides?!')
                enddo
                Opt_step = s3
                if(Opt_step.lt.0.0_dp) call stop_all(t_r,'Optimal step size is negative?!')
            endif

            !We now have the optimal step size as Opt_step :)
            norm = 0.0_dp
            do i=1,nx
                norm = norm + dx(i)**2
            enddo
            norm = sqrt(norm)
            if(abs(step)*norm.lt.1.0e-10_dp) then
                !I think we've probably got it  - we're not moving any more
                exit
            endif

            x(:) = x(:) - Opt_step*dx(:)    !Move x

        enddo
        x0(:) = x(:)  !Return the optimal vloc.

    end subroutine NR_opt

    !Numerically construct the Jacobian matrix, which is nr by nx
    !nr is number of residuals, nx is number of dimensions in variable
    subroutine MakeGradMatrix(x,g,nr,nx)
        implicit none
        real(dp), intent(in) :: x(nx)
        real(dp), intent(out) :: g(nr,nx)  !Jacobian
        integer, intent(in) :: nr,nx
        real(dp) :: step
        real(dp) :: r_1(nr),r_2(nr),x_1(nx),x_2(nx)
        integer :: i

        g(:,:) = 0.0_dp
        step = 1.0e-5_dp
        do i=1,nx
            !x with x_i increase by differential
            x_1(:) = x(:)
            x_1(i) = x(i) + step
            call RDMErr(x_1,r_1)

            !x with x_i decreased by differential
            x_2(:) = x(:)
            x_2(i) = x(i) - step
            call RDMErr(x_2,r_2)

            g(:,i) = (r_1(:)-r_2(:))/(2.0_dp*step)
        enddo

    end subroutine MakeGradMatrix

end module fitting
