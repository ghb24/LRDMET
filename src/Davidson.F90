module Davidson
    use const
    use errors, only: stop_all,warning
    use mat_tools, only: WriteVector,WriteMatrix,WriteVectorComp,WriteMatrixComp
    implicit none


    contains

    !Take a real, symmetric matrix, and find the lowest eigenvalue and eigenvector
    !Via davidson diagonalization
    subroutine Real_NonDir_Davidson(nSize,Mat,Val,Vec,tol,max_iter)
        implicit none
        integer, intent(in) :: nSize,max_iter
        real(dp), intent(in) :: Mat(nSize,nSize)
        real(dp), intent(out) :: Val,Vec(nSize)
        real(dp), intent(in) :: tol

        real(dp), allocatable :: SubspaceVecs(:,:)
        integer :: ierr
        real(dp) :: kappa

        kappa = 0.25    !Tolerance for orthogonalization procedure

        !Allocate memory for subspace vectors
        !Ridiculous overuse of memory. Should really restrict this to < max_iter, and restart if gets too much
        allocate(SubspaceVecs(nSize,max_iter),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        SubspaceVecs(:,:) = 0.0_dp

        allocate(CurrVec(nSize),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'Memory error')
        CurrVec(:) = 0.0_dp

        !Initialize current vector - random, or just 1 in first element...
        call init_vector_real(nSize,CurrVec)

        do iter = 1,max_iter

            !First, orthogonalize current vector against previous vectors, using a modified Gram-Schmidt
            !Returns normalized
            call ModGramSchmidt_real(CurrVec,nSize,SubspaceVecs(:,1:iter-1),iter-1,kappa)




    end subroutine Real_NonDir_Davidson

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

        real(dp) :: Norm_in,Norm_out,Overlap,ddot

        !Are there even any vectors to orthogonalize against!
        if(nSubspace.eq.0) then
            Norm_out = ddot(nSize,Vec,1,Vec,1)
            Vec(:) = Vec(:)/Norm_out
            return
        endif

        !Find initial normalization of vector
        Norm_in = ddot(nSize,Vec,1,Vec,1)

        do i = 1,nSubspace
            !Project out each component of the subspace
            Overlap = ddot(nSize,Vec,1,Subspace(:,i),1)
            Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
        enddo

        !Find new normalization
        Norm_out = ddot(nSize,Vec,1,Vec,1)

        if((Norm_out/Norm_in).le.kappa) then
            !Orthogonalize again!
            do i = 1,nSizespace
                Overlap = ddot(nSize,Vec,1,Subspace(:,i),1)
                Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
            enddo
            Norm_out = zdotc(nSize,Vec,1,Vec,1)
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

        complex(dp) :: Norm_in,Norm_out,Overlap,zdotc
        real(dp) :: Norm_out_re
        character(len=*), parameter :: t_r='ModGramSchmidt_comp'

        !Are there even any vectors to orthogonalize against!
        if(nSubspace.eq.0) then
            Norm_out = zdotc(nSize,Vec,1,Vec,1)
            Norm_out_re = real(Norm_out,dp)
            Vec(:) = Vec(:)/Norm_out_re
            return
        endif

        !Find initial normalization of vector
        Norm_in = ddot(nSize,Vec,1,Vec,1)
        if(abs(aimag(Norm_in)).gt.1.0e-9_dp) call stop_all(t_r,'Norm not real')

        do i = 1,nSubspace
            !Project out each component of the subspace
            Overlap = zdotc(nSize,Subspace(:,i),1,Vec,1)
            Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
        enddo

        !Find new normalization
        Norm_out = zdotc(nSize,Vec,1,Vec,1)

        if(abs(Norm_out/Norm_in).le.kappa) then
            !Orthogonalize again!
            do i = 1,nSizespace
                Overlap = zdotc(nSize,Subspace(:,i),1,Vec,1)
                Vec(:) = Vec(:) - Overlap*Subspace(:,i)     !Remove component
            enddo
            Norm_out = zdotc(nSize,Vec,1,Vec,1)
        endif
        Norm_out_re = real(Norm_out,dp)

        Vec(:) = Vec(:)/Norm_out_re

    end subroutine ModGramSchmidt_real
end module Davidson
