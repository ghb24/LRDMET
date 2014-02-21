module SelfConsistentUtils
    use const
    use errors
    use Globals
    use utils, only: get_free_unit,append_ext
    use SC_Data
    use mat_tools, only: MakeBlockHermitian,writematrixcomp,writevector
    implicit none

    contains
    
    !From numerical recepies - Find Legendre points and weights between interval
    subroutine gauleg(x1,x2,x,w,n)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: x1,x2
        real(dp), intent(out) :: x(n),w(n)
        real(dp), parameter :: eps=3.0e-14_dp   !Relative precision
        !Given the lower and upper limits of integration x1 and x2,
        !and given n, this routine returns arrays x(1:n) and w(1:n) of length
        !n, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula
        integer :: i,j,m
        real(dp) :: p1,p2,p3,pp,xl,xm,z,z1
        
        m = (n+1)/2     !Roots are symmetric in the interval, so we only have to find half of them
        xm = 0.5_dp*(x2+x1)
        xl = 0.5_dp*(x2-x1)
        do i = 1,m      !Loop over desired roots
            z = cos(pi*(i-0.25_dp)/(n + 0.5_dp))
            !Starting with the above approximation to the ith root, we enter the
            !main loop of refinement by Newtons method
1           continue
            p1 = one
            p2 = zero
            do j = 1,n
                p3 = p2
                p2 = p1
                p1 = ((2.0_dp*j-1.0_dp)*z*p2-(j-one)*p3)/j
            enddo
            !p1 is now the desired Legendre polynomial. We next computer pp,
            !its derivative, by a standard relation involving also p2, the
            !polynomial of one lower order
            pp = n*(z*p1-p2)/(z*z-one)
            z1 = z
            z = z1-p1/pp                    !Newtons method
            if(abs(z-z1).gt.eps) goto 1
            x(i) = xm-xl*z                  !Scale the root by the desired inteval
            x(n+1-i) = xm+xl*z              !and put in its symmetric counterpart
            w(i) = 2.0_dp*xl/((one-z*z)*pp*pp)  !Compute weight
            w(n+1-i) = w(i)                     !and symmetric counterpart
        enddo

    end subroutine gauleg 

    !Scale the legendre abscissas and weights to the range \pm infty
    subroutine ScaleFreqPointsWeights(n,f,w)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(inout) :: w(n),f(n)
        integer :: i

        do i = 1,n
            w(i) = w(i)*((one + f(i)**2)/(one - f(i)**2)**2)
            f(i) = -f(i)/(one-f(i)**2)
        enddo
    end subroutine ScaleFreqPointsWeights
    
    subroutine SetChemPot(GFChemPot)
        implicit none
        real(dp), intent(out) :: GFChemPot
        character(len=*), parameter :: t_r='SetChemPot'

        if(.not.tHalfFill) call stop_all(t_r,'Chemical potential currently only set for half-filling')
        !Set chemical potential (This will only be right for half-filling)
        if(.not.tAnderson) then
            !Hubbard chemical potential
            GFChemPot = U/2.0_dp
        else
            !For anderson model, we want h0v to be h0 with -U/2 on impurity
            !Then, we want no chemical potential for the linear response.
            GFChemPot = zero
        endif

    end subroutine SetChemPot
        
    subroutine SetLatticeParams(iLatParams)
        implicit none
        integer, intent(out) :: iLatParams
        integer :: ks,params_per_k,i
            
        if(tConstrainKSym) then
            !e(k) = e(-k), and we are always using a uniform mesh
            !If gamma-centered mesh, then have nSites/2 + 1 independent parameters (we are sampling k=0 and BZ boundary which dont pair)
            !If Shifted mesh, then we have nSites/2 independent parameters
            if(tShift_Mesh) then
                ks = nKPnts/2
            else
                ks = (nKPnts/2) + 1
            endif
        else
            ks = nKPnts
        endif

        if(tConstrainphsym) then
            !First nImp is the first column
            !Second nImp-2 is the second column (taking into account herm)
            !etc...
            params_per_k = 0
            do i = nImp,0,-2
                params_per_k = params_per_k + i
            enddo
        else
            !Still constrain hermiticity. Allow for complex off-diagonal matrix elements
            params_per_k = (nImp*(nImp+1)/2) + (nImp*(nImp-1)/2)
        endif

        iLatParams = params_per_k * ks

    end subroutine SetLatticeParams

    subroutine WriteLatticeParams(iLatParams,LatParams)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: LatParams(iLatParams)
        integer :: iunit,i
            
        write(6,"(A)") "Writing lattice k-space hamiltonian..."
        iunit = get_free_unit()
        open(unit=iunit,file='LatkBlocks',status='unknown')
        write(iunit,*) iLatParams
        do i = 1,iLatParams
            write(iunit,*) LatParams(i)
        enddo
        close(iunit)

    end subroutine WriteLatticeParams

    subroutine InitLatticeParams(iLatParams,LatParams,mu,ham)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(out) :: ham(nSites,nSites)  !Return the lattice hamiltonian too
!        complex(dp) :: ham_k(nSites,nSites),temp(nSites,nSites)
        integer :: i,j,iloclatcoups,iunit
        logical :: exists
        character(len=*), parameter :: t_r='InitLatticeParams'

        if(tReadCouplings) then

            write(6,"(A)") "Reading lattice k-space hamiltonian..."
            LatParams(:) = zero

            inquire(file='LatkBlocks_Read',exist=exists)
            if(.not.exists) then
                call stop_all(t_r,'LatkBlocks_Read file does not exist')
            endif
            iunit = get_free_unit()
            open(unit=iunit,file='LatkBlocks_Read',status='old')
            read(iunit,*) iloclatcoups
            if(iloclatcoups.gt.iLatParams) then
                call stop_all(t_r,'Cannot read in from larger lattice / check inputs compatible')
            elseif(iloclatcoups.eq.iLatParams) then
                write(6,"(A)") "In lattice read, number of free parameters the same. Assuming that the "    &
                    //"k-point mesh is commensurate, and just using these values."
                do i = 1,iloclatcoups
                    read(iunit,*) LatParams(i)
                enddo
                if(abs(dShiftLatticeEvals).gt.1.0e-8_dp) then
                    call ShiftLatParams(iLatParams,LatParams,dShiftLatticeEvals)
                endif
            else
                call stop_all(t_r,'Reading in from smaller lattice? This should be possible. Bug ghb24 to code up splines')
            endif

            if(tImposephsym.or.tImposeKSym) then
                call ImposeSym_2(iLatParams,LatParams,mu)
            endif
            call LatParams_to_ham(iLatParams,LatParams,mu,ham)
        else
            !Just set to h0v
            ham = zzero
            do i = 1,nSites
                do j = 1,nSites
                    ham(j,i) = cmplx(h0v(j,i),zero,dp)
                enddo
            enddo

!            call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham,nSites,zzero,temp,nSites)
!            call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RtoK_Rot,nSites,zzero,ham_k,nSites)
!            call writematrixcomp(ham_k,'Initial k space ham',.true.)
!
!            call writematrixcomp(ham,'Initial real space ham',.true.)
            call ham_to_LatParams(iLatParams,LatParams,ham)
!            ham = zzero
!            call LatParams_to_ham(iLatParams,LatParams,mu,ham)
!            call writematrixcomp(ham,'real space ham converted from latparams',.true.)
        endif

    end subroutine InitLatticeParams

    subroutine ShiftLatParams(iLatParams,LatParams,dShift)
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(inout) :: LatParams(iLatParams)
        real(dp), intent(in) :: dShift

    end subroutine ShiftLatParams

    subroutine ham_to_KBlocks(ham,KBlocks)
        complex(dp), intent(in) :: ham(nSites,nSites) 
        complex(dp), intent(out) :: KBlocks(nImp,nImp,nKPnts)
        complex(dp), allocatable :: ctemp(:,:)
        integer :: k,ind_1,ind_2
        
        KBlocks(:,:,:) = zzero

        allocate(ctemp(nSites,nImp))
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            call ZGEMM('N','N',nSites,nImp,nSites,zone,ham,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ctemp,nSites)
            call ZGEMM('C','N',nImp,nImp,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ctemp,nSites,zzero,KBlocks(:,:,k),nImp)
        enddo
        deallocate(ctemp)

    end subroutine ham_to_KBlocks

    subroutine ham_to_LatParams(iLatParams,LatParams,ham)
        integer, intent(in) :: iLatParams
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(in) :: ham(nSites,nSites)
        complex(dp), allocatable :: KBlocks(:,:,:)

        !First, ham to k_blocks
        allocate(KBlocks(nImp,nImp,nKPnts))
        call ham_to_KBlocks(ham,KBlocks)

        call KBlocks_to_LatParams(iLatParams,LatParams,KBlocks)
        deallocate(KBlocks)

    end subroutine ham_to_LatParams

    !Convert a set of variational parameters back into a real space one-electron hamiltonian
    subroutine LatParams_to_ham(iLatParams,LatParams,mu,ham)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: LatParams(iLatParams)
        complex(dp), intent(out) :: ham(nSites,nSites)
        complex(dp), allocatable :: KBlocks(:,:,:)
        complex(dp), allocatable :: ctemp(:,:)
        integer :: k,ind_1,ind_2
        character(len=*), parameter :: t_r='LatParams_to_ham'

        allocate(KBlocks(nImp,nImp,nKPnts))
        call LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)

        ham = zzero

        !Now, rotate back into real space from kspace
        allocate(ctemp(nSites,nImp))
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            call ZGEMM('N','N',nSites,nImp,nImp,zone,RtoK_Rot(:,ind_1:ind_2),nSites,KBlocks(:,:,k),nImp,zzero,ctemp,nSites)
            call ZGEMM('N','C',nSites,nSites,nImp,zone,ctemp,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zone,ham,nSites)
        enddo
        deallocate(ctemp)
        deallocate(KBlocks)

!        do i = 1,nKPnts
!            ind_1 = (i-1)*nImp + 1
!            ind_2 = i*nImp
!            ham(ind_1:ind_2,ind_1:ind_2) = KBlocks(:,:,i)
!        enddo
!        allocate(ctemp(nSites,nSites))
!        call ZGEMM('N','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham,nSites,zzero,ctemp,nSites)
!        call ZGEMM('N','C',nSites,nSites,nSites,zone,ctemp,nSites,RtoK_Rot,nSites,zzero,ham,nSites)
!        deallocate(ctemp)

    end subroutine LatParams_to_ham

    subroutine KBlocks_to_LatParams(iLatParams,LatParams,KBlocks)
        integer, intent(in) :: iLatParams
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(in) :: KBlocks(nImp,nImp,nKPnts)
        integer :: i,j,k,ind,ks
        character(len=*), parameter :: t_r='KBlocks_to_LatParams'

        LatParams(:) = zero
        
        if(tConstrainKSym) then
            !e(k) = e(-k), and we are always using a uniform mesh
            !If gamma-centered mesh, then have nSites/2 + 1 independent parameters (we are sampling k=0 and BZ boundary which dont pair)
            !If Shifted mesh, then we have nSites/2 independent parameters
            if(tShift_Mesh) then
                ks = nKPnts/2
            else
                ks = (nKPnts/2) + 1
            endif
        else
            ks = nKPnts
        endif

        ind = 1
        do k = 1,ks
            if(tConstrainphsym) then
                do i = 1,nImp
                    do j = i,nImp-i+1
                        if(i.eq.j) then
                            if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            LatParams(ind) = real(KBlocks(j,i,k),dp)
                            ind = ind + 1
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            LatParams(ind) = real(KBlocks(j,i,k),dp)
                            LatParams(ind+1) = aimag(KBlocks(j,i,k))
                            ind = ind + 2
                        endif
                    enddo
                enddo
            else
                do i = 1,nImp
                    !Run through columns
                    do j = i,nImp
                        !Run through rows
                        if(i.eq.j) then
                            if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            LatParams(ind) = real(KBlocks(j,i,k),dp)
                            ind = ind + 1
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            LatParams(ind) = real(KBlocks(j,i,k),dp)
                            LatParams(ind+1) = aimag(KBlocks(j,i,k))
                            ind = ind + 2
                        endif
                    enddo
                enddo
            endif
        enddo

    end subroutine KBlocks_to_LatParams

!    subroutine RealSpace_EValsVecs(EVecs_k,EVals_k,EVecs_r,EVals_r)
!        complex(dp), intent(in) :: EVecs_k(nImp,nImp,nKPnts)
!        complex(dp), intent(in) :: EVals_k(nSites)
!        complex(dp), intent(out) :: EVecs_r(nSites,nSites)
!        complex(dp), intent(out) :: EVals_r(nSites)
!
!        do k = 1,nKPnts
!
!            call ZGEMM('N','N',
!
!    end subroutine RealSpace_EValsVecs

    subroutine KBlocks_to_diag(KBlocks,EVecs,EVals)
        complex(dp), intent(in) :: KBlocks(nImp,nImp,nKPnts)
        complex(dp), intent(out) :: EVecs(nImp,nImp,nKPnts)
        complex(dp), intent(out) :: EVals(nSites)
        complex(dp) :: KBlock(nImp,nImp)
        complex(dp), allocatable :: cWork(:)
        real(dp) :: Bands_k(nImp)
        real(dp), allocatable :: rWork(:)
        integer :: k,lWork,ierr,i
        character(len=*), parameter :: t_r='KBlocks_to_diag'

        EVals = zzero
        EVecs = zzero

        allocate(rWork(max(1,3*nImp-2)))
        do k = 1,nKPnts

            EVecs(:,:,k) = KBlocks(:,:,k)
            !Diagonalize each kBlock
            allocate(cWork(1))
            lWork = -1
            ierr = 0
            call zheev('V','U',nImp,EVecs(:,:,k),nImp,Bands_k,cWork,lWork,rWork,ierr)
            if(ierr.ne.0) call stop_all(t_r,'Error in diag')
            lWork = int(real(cWork(1))) + 1
            deallocate(cWork)
            allocate(cWork(lWork))
            call zheev('V','U',nImp,EVecs(:,:,k),nImp,Bands_k,cWork,lWork,rWork,ierr)
            if(ierr.ne.0) call stop_all(t_r,'Error in diag')
            deallocate(cWork)

            do i = 1,nImp
                EVals(((k-1)*nImp)+i) = cmplx(Bands_k(i),zero,dp)
            enddo

        enddo
        deallocate(rWork)

    end subroutine KBlocks_to_diag

    subroutine LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: LatParams(iLatParams)
        complex(dp), intent(out) :: KBlocks(nImp,nImp,nKPnts)
        integer :: ind,i,j,ks,k
        character(len=*), parameter :: t_r='LatParams_to_KBlocks'

        if((nImp.ne.1).and.(mod(nImp,2).ne.0)) then
            call stop_all(t_r,'Cannot deal with non-1 odd numbers of impurity sites')
        endif

        KBlocks(:,:,:) = zzero
        
        if(tConstrainKSym) then
            !e(k) = e(-k), and we are always using a uniform mesh
            !If gamma-centered mesh, then have nSites/2 + 1 independent parameters (we are sampling k=0 and BZ boundary which dont pair)
            !If Shifted mesh, then we have nSites/2 independent parameters
            if(tShift_Mesh) then
                ks = nKPnts/2
            else
                ks = (nKPnts/2) + 1
            endif
        else
            ks = nKPnts
        endif

        ind = 1
        do k = 1,ks
            if(tConstrainphsym) then
                do i = 1,nImp
                    do j = i,nImp-i+1
                        if(i.eq.j) then
                            if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            KBlocks(j,i,k) = cmplx(LatParams(ind),zero,dp)
                            ind = ind + 1
                            !For the diagonals, we want to flip the energy about the chemical potential
                            KBlocks((nImp-i)+1,(nImp-j)+1,k) = cmplx(2.0_dp*mu,zero,dp) - KBlocks(j,i,k)
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            KBlocks(j,i,k) = cmplx(LatParams(ind),LatParams(ind+1),dp)
                            ind = ind + 2
                            !Put in the same coupling
                            !Swap indices & reverse order
                            KBlocks((nImp-i)+1,(nImp-j)+1,k) = KBlocks(j,i,k)
                        endif
                    enddo
                enddo
            else
                !The variables are the lower triangle of the matrix for this kpoint
                do i = 1,nImp
                    !Run through columns
                    do j = i,nImp
                        !Run through rows
                        if(i.eq.j) then
                            if(ind.gt.iLatParams) then
                                write(6,*) "i,j,k: ",i,j,k
                                write(6,*) "ind: ",ind
                                write(6,*) "iLatParams: ",iLatParams
                                call stop_all(t_r,'Incorrect indexing 1')
                            endif
                            KBlocks(j,i,k) = cmplx(LatParams(ind),zero,dp)
                            ind = ind + 1
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing 2')
                            KBlocks(j,i,k) = cmplx(LatParams(ind),LatParams(ind+1),dp)
                            ind = ind + 2
                        endif
                    enddo
                enddo
            endif
            call MakeBlockHermitian(KBlocks(:,:,k),nImp)
        enddo

        if(tConstrainKSym) then
            !We have only filled up half the states. Fill the others.
            !Update this so that it uses a function to map to equivalent kpoints
            !Currently this is only going to work for 1D

            if(tShift_Mesh) then
                !No gamma point sampled. All k-points symmetric.
                !Mirror the k-space hamiltonian
                do i = 1,ks
                    KBlocks(:,:,i+ks) = dconjg(KBlocks(:,:,ks-i+1))
                enddo
            else
                !Mirror the kpoints, but ignore the gamma point and BZ boundary
                do i = 2,ks-1
                    KBlocks(:,:,i+ks-1) = dconjg(KBlocks(:,:,ks-i+1))
                enddo
            endif
        endif

    end subroutine LatParams_to_KBlocks

    !This will set FreqPoints and Weights. nFitPoints is the size of FreqPoints/Weights, 
    !and may refer to real or im axis based on the tFitRealFreq flag.
    subroutine SetFreqPoints(nFreq_Re,nFreq_Im,nFitPoints)
        implicit none
        integer, intent(out) :: nFreq_Re,nFreq_Im
        integer, intent(out) :: nFitPoints
        logical :: tUsePoswFreqPoints
        integer :: OmegaVal,i,j
        real(dp) :: Omega,LowFreq,HighFreq
        real(dp), allocatable :: FreqPoints_temp(:),Weights_temp(:)
        character(len=*), parameter :: t_r='SetFreqPoints'

        if(tFitRealFreq.and.tFitPoints_Legendre) then
            call stop_all(t_r,'Trying to fit the real frequency spectrum to Legendre points?')
        endif
        if(.not.tHalfFill.and.tUsePoswFreqPoints) then
            call stop_all(t_r,'Trying to just optimize the positive frequency points, but not half filled system')
        endif
        if(tHalfFill.and.(nImp.eq.1)) then
            tUsePoswFreqPoints = .true. !For ph symmetric systems, this seems to make no difference (but is faster!)
        else
            tUsePoswFreqPoints = .false.
        endif

        nFreq_Re = 0
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,.false.)
            if(OmegaVal.lt.0) exit
            nFreq_Re = nFreq_Re + 1
        enddo
        
        if(tFitPoints_Legendre) then
            !Calculate the Legendre polynomials for Legendre integration points and weights
            allocate(FreqPoints(nFreqPoints))
            allocate(Weights(nFreqPoints))
            write(6,"(A)") "Calculating Legendre frequency points..."
            LowFreq = -1.0_dp
            HighFreq = 1.0_dp
            call gauleg(LowFreq,HighFreq,FreqPoints,Weights,nFreqPoints)
            call ScaleFreqPointsWeights(nFreqPoints,FreqPoints,Weights)
            if(tUsePoswFreqPoints) then
                !Count the positive values of w
                write(6,"(A)") "Just using the positive values of omega"
                j = 0
                do i = 1,nFreqPoints
                    if(FreqPoints(i).ge.zero) j = j+1
                enddo
                allocate(FreqPoints_temp(j))
                allocate(Weights_temp(j))
                j = 0
                do i = 1,nFreqPoints
                    if(FreqPoints(i).ge.zero) then
                        j = j+1
                        FreqPoints_temp(j) = FreqPoints(i)
                        Weights_temp(j) = Weights(i)
                    endif
                enddo
                deallocate(FreqPoints,Weights)
                nFreqPoints = j
                allocate(FreqPoints(nFreqPoints))
                allocate(Weights(nFreqPoints))
                FreqPoints(:) = FreqPoints_temp(:)
                Weights(:) = Weights_temp(:)
                deallocate(FreqPoints_temp,Weights_temp)
                write(6,"(A,I7)") "New number of frequency points: ",nFreqPoints
            endif
            nFreq_Im = nFreqPoints
            write(6,"(A)") "Frequency integration points and weights: "
            do i = 1,nFreqPoints
                write(6,"(I9,2G20.13)") i,FreqPoints(i),Weights(i)
            enddo
            !Which axis do we want to do the fitting on?
            if(tFitRealFreq) then
                nFitPoints = nFreq_Re   
                tCalcRealSpectrum = .false.
                tFitMatAxis = .false.
            else
                nFitPoints = nFreq_Im
                tFitMatAxis = .true.
            endif
        
        else
            nFreq_Im = 0
            OmegaVal = 0
            do while(.true.)
                call GetNextOmega(Omega,OmegaVal,.true.)
                if(OmegaVal.lt.0) exit
                nFreq_Im = nFreq_Im + 1
            enddo
            !Which axis do we want to do the fitting on?
            if(tFitRealFreq) then
                nFitPoints = nFreq_Re   
                tCalcRealSpectrum = .false.
                tFitMatAxis = .false.
            else
                nFitPoints = nFreq_Im
                tFitMatAxis = .true.
            endif
            allocate(FreqPoints(nFitPoints))
            allocate(Weights(nFitPoints))
            OmegaVal = 0
            do while(.true.)
                call GetNextOmega(Omega,OmegaVal,tFitMatAxis)
                if(OmegaVal.lt.0) exit
                FreqPoints(OmegaVal) = Omega
                Weights(OmegaVal) = one
            enddo
        endif

        !Initially, just see if we can fit the two different Matsubara spectral functions
        write(6,*) "Number of frequency points to fit across: ",nFitPoints
        call flush(6)

    end subroutine SetFreqPoints

    !Write out the isotropic average, and individual elements of a dynamic function in the impurity space.
    subroutine writedynamicfunction(n,Func,FileRoot,tag,tCheckCausal,tCheckOffDiagHerm,tWarn,tMatbrAxis,ErrMat,FreqPoints)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: Func(nImp,nImp,n)
        character(len=*), intent(in) :: FileRoot
        integer, intent(in), optional :: tag
        logical, intent(in), optional :: tCheckCausal
        logical, intent(in), optional :: tCheckOffDiagHerm
        logical, intent(in), optional :: tWarn  !If true, don't die if non-causal
        logical, intent(in), optional :: tMatbrAxis !If true, then this correlation function is defined on the Im Axis
        complex(dp), intent(in), optional :: ErrMat(nImp,nImp,n)    !Optional errors on the function
        real(dp), intent(in), optional :: FreqPoints(n)

        character(64) :: filename
        logical :: tCheckOffDiagHerm_,tCheckCausal_,tWarn_,tMatbrAxis_
        integer :: iunit,i,j,k,imp1,imp2
        real(dp) :: Omega,Prev_Spec,SpectralWeight
        complex(dp) :: IsoAv,IsoErr
        character(len=*), parameter :: t_r='writedynamicfunction'

        if(present(tCheckCausal)) then
            tCheckCausal_ = tCheckCausal
        else
            tCheckCausal_ = .false.
        endif

        if(present(tCheckOffDiagHerm)) then
            tCheckOffDiagHerm_ = tCheckOffDiagHerm
        else
            tCheckOffDiagHerm_ = .false.
        endif

        if(present(tWarn)) then
            tWarn_ = tWarn
        else
            tWarn_ = .false.    !i.e. die by default
        endif

        if(present(tMatbrAxis)) then
            tMatbrAxis_ = tMatbrAxis
        else
            tMatbrAxis_ = .false.
        endif
        if(tMatbrAxis_) tCheckCausal_ = .false. !The imaginary part won't be negative on the matsubara axis (fn is antisym)

        if(present(tag)) then
            call append_ext(FileRoot,tag,filename)
        else
            filename = FileRoot
        endif
        iunit = get_free_unit()
        open(unit=iunit,file=filename,status='unknown')

        SpectralWeight = 0.0_dp
        Prev_Spec = 0.0_dp
        i=0
        do while(.true.)
            if(present(FreqPoints)) then
                call GetNextOmega(Omega,i,tMatbrAxis_,nFreqPoints=n,FreqPoints=FreqPoints)
            else
                call GetNextOmega(Omega,i,tMatbrAxis_)
            endif
            if(i.lt.0) exit
            if(i.gt.n) call stop_all(t_r,'Wrong number of frequency points used')
            !Find isotropic part
            IsoAv = zzero
            IsoErr = zzero
            do j = 1,nImp
                IsoAv = IsoAv + Func(j,j,i)
                if(present(ErrMat)) IsoErr = IsoErr + ErrMat(j,j,i)
                if(tCheckOffDiagHerm_) then
                    do k = 1,nImp
                        if(k.ne.j) then
                            if(abs(Func(k,j,i)-dconjg(Func(j,k,i))).gt.1.0e-5) then
                                write(6,*) "While writing file: ",filename
                                if(present(tag)) write(6,*) "Filename extension: ",tag
                                write(6,*) "Element of function: ",k,j
                                write(6,*) "Frequency point: ",Omega
                                write(6,*) "Values: ",Func(k,j,i),Func(j,k,i),abs(Func(k,j,i)-dconjg(Func(j,k,i)))
                                if(tWarn_) then
                                    call warning(t_r,'Function no longer off-diagonal hermitian')
                                else
                                    call stop_all(t_r,'Function no longer off-diagonal hermitian')
                                endif
                            endif
                        endif
                    enddo
                endif
            enddo
            IsoAv = IsoAv / real(nImp,dp)
            IsoErr = IsoErr / real(nImp,dp)
            if(tCheckCausal_.and.(aimag(IsoAv).gt.1.0e-8_dp)) then
                write(6,*) "While writing file: ",filename
                if(present(tag)) write(6,*) "Filename extension: ",tag
                write(6,*) "Frequency point: ",Omega
                write(6,*) "Isotropic value: ",IsoAv,Func(1,1,i)
                if(tWarn_) then
                    call warning(t_r,'Function not causal')
                else
                    call stop_all(t_r,'Function not causal')
                endif
            endif

            if(present(ErrMat)) then
                write(iunit,"(5G25.10)") Omega,real(IsoAv,dp),aimag(IsoAv),real(IsoErr,dp),aimag(IsoErr)
            else
                if(i.ne.1) then
                    if(tMatbrAxis_) then
                        SpectralWeight = SpectralWeight + Omega_Step_Im*(Prev_Spec-aimag(IsoAv))/(2.0_dp*pi)
                    else
                        SpectralWeight = SpectralWeight + Omega_Step*(Prev_Spec-aimag(IsoAv))/(2.0_dp*pi)
                    endif
                    Prev_Spec = -aimag(IsoAv)
                endif
                write(iunit,"(3G25.10)",advance='no') Omega,real(IsoAv,dp),aimag(IsoAv)
                do imp1=1,nImp
                    do imp2=1,nImp
                        write(iunit,"(2G25.10)",advance='no') real(Func(imp2,imp1,i),dp),aimag(Func(imp2,imp1,i))
                    enddo
                enddo
                write(iunit,*) 
            endif
        enddo
        if(tMatbrAxis_) then
            write(6,"(A,F17.10)") "Total approximate spectral weight on Matsubara axis: ",SpectralWeight
        else
            write(6,"(A,F17.10)") "Total approximate spectral weight on real axis: ",SpectralWeight
        endif
        close(iunit)

    end subroutine writedynamicfunction
            
    subroutine CheckGFsSame(n,G,G_,dTol)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: G(nImp,nImp,n)
        complex(dp), intent(in) :: G_(nImp,nImp,n)
        real(dp) , intent(in) :: dTol
        character(len=*), parameter :: t_r='CheckGFsSame'
        integer :: i,j,k

        do i = 1,n
            do j = 1,nImp
                do k = 1,nImp
                    if(abs(G(k,j,i)-G_(k,j,i)).gt.dTol) then
                        write(6,*) "Frequency point: ",i
                        write(6,*) "Element: ",j,k
                        write(6,*) "Values: ",G(k,j,i),G_(k,j,i)
                        call stop_all(t_r,'Local functions not the same')
                    endif
                enddo
            enddo
        enddo

    end subroutine CheckGFsSame
    
    !This function inverts a local greens function matrix, allowing it to be non-hermitian
    !Function is sent in as the normal greens function, n, nImp x nImp matrices, and the inverse
    !for each frequency is returned
    subroutine InvertLocalNonHermFunc(n,InvGF)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        use matrixops, only: mat_inv
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(inout) :: InvGF(nImp,nImp,n)
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: LVec(:,:),RVec(:,:),W_Vals(:),cWork(:)
        complex(dp) :: ztemp2(nImp,nImp),IGF(nImp,nImp)
        integer :: i,lwork,info,j
        character(len=*), parameter :: t_r='InvertLocalNonHermFunc'

        if(.true.) then
            do i = 1,n
                call mat_inv(InvGF(:,:,i),IGF)
                InvGF(:,:,i) = IGF(:,:)
            enddo
        else
            !Invert via diagonalization
            allocate(LVec(nImp,nImp))
            allocate(RVec(nImp,nImp))
            allocate(W_Vals(nImp))
            allocate(Work(max(1,2*nImp)))

            do i = 1,n

                !Diagonalize
                LVec(:,:) = zzero
                RVec(:,:) = zzero
                W_Vals(:) = zzero
                IGF(:,:) = InvGF(:,:,i)

                allocate(cWork(1))
                lwork = -1
                info = 0
                call zgeev('V','V',nImp,IGF,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(abs(cWork(1)))+1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zgeev('V','V',nImp,IGF,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
                deallocate(cWork)

                !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
                !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
                call Order_zgeev_vecs(W_Vals,LVec,RVec)
                !call writevectorcomp(W_Vals,'Eigenvalues ordered')
                !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
                call Orthonorm_zgeev_vecs(nImp,W_Vals,LVec,RVec)
                !Calculate greens function for this k-vector
                IGF = zzero
                do j = 1,nImp
                    IGF(j,j) = zone / W_Vals(j)
                enddo
                !Now rotate this back into the original basis
                call zGEMM('N','N',nImp,nImp,nImp,zone,RVec,nImp,IGF,nImp,zzero,ztemp2,nImp)
                call zGEMM('N','C',nImp,nImp,nImp,zone,ztemp2,nImp,LVec,nImp,zzero,IGF,nImp)

                InvGF(:,:,i) = IGF(:,:)
            enddo

            deallocate(LVec,RVec,W_Vals,Work)
        endif

    end subroutine InvertLocalNonHermFunc
    
    !From num rec
    subroutine splint(xa,ya,y2a,n,x,y)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: x, xa(n),y2a(n),ya(n)
        real(dp), intent(out) :: y
        !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa_i's in order),
        !and given the array y2a(1:n), which is the output from spline, and given a value of x, this routine
        !returns a cubic-spline interpolated value y.
        integer :: k,khi,klo
        real(dp) :: a, b, h
        character(len=*), parameter :: t_r='splint'

        klo = 1
        khi = n
1       if (khi-klo.gt.1) then
            !Bisection to find right place in table
            k = (khi+klo)/2
            if(xa(k).gt.x) then
                khi = k
            else
                klo = k
            endif
            goto 1
        endif
        !klo and khi now bracket the value of x
        h = xa(khi)-xa(klo)
        if(h.eq.zero) call stop_all(t_r,'bad xa input in splint. kpoints should be distinct')
        !Evaluate spline
        a = (xa(khi)-x)/h
        b = (x-xa(klo))/h
        y = a*ya(klo)+b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
    end subroutine splint


    !From numerical recipies
    subroutine spline(x,y,n,yp1,ypn,y2)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: yp1,ypn,x(n),y(n)
        real(dp), intent(out) :: y2(n)
        !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e. y_i = f(x_i), with
        !x_1 < x_2 < x_3 < ... < x_n, and given values yp1 and ypn for the first derivative of the
        !interpolating function at points 1 and n, respectively, this routine returns an array y2(1:n)
        !of length n which contains the second derivatives of the interpolating function at the tabulated
        !points x_i. If yp1 and/or ypn are > 1e30, the routine is signaled to set the corresponding boundary
        !condition for a natural spline, with zero second derivative on that boundary.
        integer :: i,k
        real(dp) :: p,qn,sig,un
        real(dp), allocatable :: u(:)

        allocate(u(n))
        if(yp1.gt.0.99e30_dp) then
            !The lower boundary condition is set to be 'natural'
            y2(1) = zero
            u(1) = zero
        else
            !A specified first derivative
            y2(1)=-0.5_dp
            u(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        do i=2,n-1
            !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temp storage of the decomposed factors
            sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
            p=sig*y2(i-1)+2.0_dp
            y2(i)=(sig-one)/p
            u(i) = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        enddo
        if(ypn.gt.0.99e30_dp) then
            !Upper boundary condition is 'natural'
            qn = zero
            un = zero
        else
            !Specified first derivative
            qn = 0.5_dp
            un = (3.0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)
        do k=n-1,1,-1
            !This is the backsubstitution loop of the tridiagonal algorithm
            y2(k) = y2(k)*y2(k+1)+u(k)
        enddo
        deallocate(u)
    end subroutine spline
    
    !Get next omega value for optionally real/matsubara axis
    !OmegaVal = 0 on input initializes
    !OmegaVal = -1 on exit means we have finished
    subroutine GetNextOmega(Omega,OmegaVal,tMatbrAxis,nFreqPoints,FreqPoints)
        implicit none
        real(dp), intent(inout) :: Omega
        integer, intent(inout) :: OmegaVal
        logical, intent(in) :: tMatbrAxis
        integer, intent(in), optional :: nFreqPoints
        real(dp), intent(in), optional :: FreqPoints(:)
        logical :: tUseFreqArray
        character(len=*), parameter :: t_r='GetNextOmega'

!        write(6,*) "In Get Next Omega: ",OmegaVal
        if(present(FreqPoints)) then
            tUseFreqArray = .true.
        else
            tUseFreqArray = .false.
        endif

        if(OmegaVal.eq.0) then
            if(tUseFreqArray) then
                Omega = FreqPoints(1)
            else
                if(tMatbrAxis) then
                    Omega = Start_Omega_Im
                else
                    Omega = Start_Omega
                endif
            endif
            OmegaVal = OmegaVal + 1
!            write(6,*) "Returned Omega: ",Omega
            return
        endif

        if(tUseFreqArray) then
            if(.not.present(nFreqPoints)) call stop_all(t_r,'Error here')
            if(OmegaVal.eq.nFreqPoints) then
                OmegaVal = -1
                return
            else
                Omega = FreqPoints(OmegaVal+1)
            endif
        elseif(tMatbrAxis) then
            !Omega = Omega + Omega_Step_Im
            Omega = Start_Omega_Im + OmegaVal*Omega_Step_Im
            if((Omega.gt.(max(Start_Omega_Im,End_Omega_Im)+1.0e-6_dp))  &
                    .or.(Omega.lt.(min(Start_Omega_Im,End_Omega_Im)-1.0e-6_dp))) then
                !Hit exit criterion
                OmegaVal = -1
                return
            endif
        else
            !Omega = Omega + Omega_Step
            Omega = Start_Omega + OmegaVal*Omega_Step
            if((Omega.gt.(max(Start_Omega,End_Omega)+1.0e-6_dp))  &
                    .or.(Omega.lt.(min(Start_Omega,End_Omega)-1.0e-6_dp))) then
                !Hit exit criterion
                OmegaVal = -1
                return
            endif
        endif
        OmegaVal = OmegaVal + 1
!        write(6,*) "Returned Omega: ",Omega

    end subroutine GetNextOmega
    
    subroutine Imposesym_2(iLatParams,vars,mu)
        use sort_mod, only: sort_real
        implicit none
        integer, intent(in) :: iLatParams
        real(dp), intent(inout) :: vars(iLatParams)
        real(dp), intent(in) :: mu
        integer :: i,k,lWork,iErr
        real(dp) :: DistFromMu,Bands_k(nImp)
        integer :: ks
        real(dp), allocatable :: rWork(:)
        complex(dp) :: KBlock(nImp,nImp),KBlock2(nImp,nImp),cTemp(nImp,nImp)
        complex(dp), allocatable :: KBlocks(:,:,:),cWork(:)
        character(len=*), parameter :: t_r='Imposesym_2'
            
        allocate(KBlocks(nImp,nImp,nKPnts))
        call LatParams_to_KBlocks(iLatParams,vars,mu,KBlocks)
            
        if(tShift_Mesh) then
            ks = nKPnts/2
        else
            ks = (nKPnts/2) + 1
        endif

        if(tImposeKSym) then
            !After fitting, add back in momentum inversion symmetry (if we are working in k-space)
            !This means that h(k) = h*(-k)
            !For shifted meshes, this is easy.
            !For Gamma-centered meshes, two k-points are only singly degenerate
            !We should not be able to be constraining any syms
            if(tShift_Mesh) then
                do i = 1,nKPnts/2
                    KBlock(:,:) = (KBlocks(:,:,i) + dconjg(KBlocks(:,:,nKPnts-i+1))) / 2.0_dp
                    KBlocks(:,:,i) = KBlock(:,:)
                    KBlocks(:,:,nKPnts-i+1) = dconjg(KBlock(:,:))
                enddo
            else
                do i = 2,nKPnts/2
                    KBlock(:,:) = (KBlocks(:,:,i) + dconjg(KBlocks(:,:,nKPnts-i+2))) / 2.0_dp
                    KBlocks(:,:,i) = KBlock(:,:)
                    KBlocks(:,:,nKPnts-i+2) = dconjg(KBlock(:,:))
                enddo
            endif
        endif

        if(tImposephsym) then
        
            !write(6,*) "*** IMPOSING PH SYMMETRY ***",mu
            if(mod(nImp,2).eq.0) then
                !Multiple of two bands per kpoint. Constrain them so that they are in pairs
                allocate(rWork(max(1,3*nImp-2)))
                do k = 1,nKPnts
                    KBlock(:,:) = KBlocks(:,:,k)
                    !Diagonalize each kBlock
                    allocate(cWork(1))
                    lWork = -1
                    ierr = 0
                    call zheev('V','U',nImp,KBlock,nImp,Bands_k,cWork,lWork,rWork,ierr)
                    if(ierr.ne.0) call stop_all(t_r,'Error in diag')
                    lWork = int(real(cWork(1))) + 1
                    deallocate(cWork)
                    allocate(cWork(lWork))
                    call zheev('V','U',nImp,KBlock,nImp,Bands_k,cWork,lWork,rWork,ierr)
                    if(ierr.ne.0) call stop_all(t_r,'Error in diag')
                    deallocate(cWork)

                    !write(6,*) "For k point: ",k
                    !call writematrixcomp(KBlocks(:,:,k),'k-space ham',.true.)
                    !call writevector(Bands_k,'eigenvalues at kpoint')

                    !Now pair up the bands in this kblock
                    do i = 1,nImp/2
                        DistFromMu = 0.5_dp * (- Bands_k(i) + Bands_k(nImp-i+1) )
                        !write(6,*) "Dist from mu: ",DistFromMu 
                        Bands_k(i) = mu - DistFromMu
                        Bands_k(nImp-i+1) = mu + DistFromMu
                    enddo
                    !call writevector(Bands_k,'ph symmetrized eigenvalues at kpoint')

                    !Now, rotate the block back into k-space using the same eigenvectors
                    KBlock2(:,:) = zzero
                    do i = 1,nImp
                        KBlock2(i,i) = cmplx(Bands_k(i),zero,dp)
                    enddo
                    call ZGEMM('N','N',nImp,nImp,nImp,zone,KBlock,nImp,KBlock2,nImp,zzero,cTemp,nImp)
                    call ZGEMM('N','C',nImp,nImp,nImp,zone,cTemp,nImp,KBlock,nImp,zzero,KBlocks(:,:,k),nImp)
                    !call writematrixcomp(KBlocks(:,:,k),'new k-space ham',.true.)
                enddo
                deallocate(rWork)
            elseif(nImp.eq.1) then
                !If there is only one band per kpoint, then the pairs correspond to different kpoints.
                do i = 1,ks/2
                    DistFromMu = 0.5_dp*( (mu - real(KBlocks(1,1,ks-i+1)) ) + (real(KBlocks(1,1,i)) - mu))
                    KBlocks(1,1,ks-i+1) = cmplx(mu - DistFromMu,zero,dp)
                    KBlocks(1,1,i) = cmplx(mu + DistFromMu,zero,dp)
                enddo

                !Since we are mixing kpoint, this might break k-symmetry. Reimpose it if it is being constrained
                if(tImposeKSym) then
                    !After fitting, add back in momentum inversion symmetry (if we are working in k-space)
                    !This means that e(k) = e(-k)
                    !For shifted meshes, this is easy.
                    !For Gamma-centered meshes, two k-points are only singly degenerate
                    !We should not be able to be constraining any syms
                    if(tShift_Mesh) then
                        do i = 1,nKPnts/2
                            KBlock(:,:) = (KBlocks(:,:,i) + KBlocks(:,:,nKPnts-i+1)) / 2.0_dp
                            KBlocks(:,:,i) = KBlock(:,:)
                            KBlocks(:,:,nKPnts-i+1) = KBlock(:,:)
                        enddo
                    else
                        do i = 2,nKPnts/2
                            KBlock(:,:) = (KBlocks(:,:,i) + KBlocks(:,:,nKPnts-i+2)) / 2.0_dp
                            KBlocks(:,:,i) = KBlock(:,:)
                            KBlocks(:,:,nKPnts-i+2) = KBlock(:,:)
                        enddo
                    endif
                endif
            else
                call stop_all(t_r,'Cannot impose ph symmetry if nImp != 1 or multiple of 2')
            endif
        endif
            
        call KBlocks_to_LatParams(iLatParams,vars,KBlocks)
        deallocate(KBlocks)

    end subroutine Imposesym_2

end module SelfConsistentUtils
