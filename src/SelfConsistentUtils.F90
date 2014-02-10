module SelfConsistentUtils
    use const
    use errors
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
        real(dp), intent(in) :: LatParams
            
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
        real(dp), intent(out) :: LatParams(iLatParams)
        complex(dp), intent(out) :: ham(nSites,nSites)  !Return the lattice hamiltonian too
        integer :: i,j
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
                    read(iunit,*) kval,LatParams(i)
                enddo
                if(abs(dShiftLatticeEvals).gt.1.0e-8_dp) then
                    call ShiftLatParams(iLatParams,LatParams,dShiftLatticeEvals)
                endif
            else
                call stop_all(t_r,'Reading in from smaller lattice? This should be possible. Bug ghb24 to code up splines')
            endif

            if(tImposephsym.or.tImposeKSym) then
                call ImposeSym(iLatParams,LatParams)
            endif
        else
            !Just set to h0v
            do i = 1,nSites
                do j = 1,nSites
                    ham(j,i) = cmplx(h0v(j,i),zero,dp)
                enddo
            enddo

            call ham_to_LatParams(iLatParams,LatParams,mu,ham) 
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
        complex(dp), allocatable :: ctemp(:,:),cham(:,:)
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
        integer :: i,j,k,ind
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

    subroutine LatParams_to_KBlocks(iLatParams,LatParams,mu,KBlocks)
        integer, intent(in) :: iLatParams
        real(dp), intent(in) :: mu
        real(dp), intent(in) :: LatParams(iLatParams)
        complex(dp), intent(out) :: KBlocks(nImp,nImp,nKPnts)

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
                            KBlocks((nImp-i)+1,(nImp-j)+1) = cmplx(2.0_dp*mu,zero,dp) - KBlocks(j,i,k)
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            KBlocks(j,i,k) = cmplx(LatParams(ind),LatParams(ind+1),dp)
                            ind = ind + 2
                            !Put in the same coupling
                            !Swap indices & reverse order
                            KBlocks((nImp-i)+1,(nImp-j)+1) = KBlocks(j,i,k)
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
                            if(ind.gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            KBlocks(j,i,k) = cmplx(LatParams(ind),zero,dp)
                            ind = ind + 1
                        else
                            if((ind+1).gt.iLatParams) call stop_all(t_r,'Incorrect indexing')
                            KBlocks(j,i,k) = cmplx(LatParams(ind),LatParams(ind+1),dp)
                            ind = ind + 2
                        endif
                    enddo
                enddo
            endif
            call MakeBlockHermitian(Block,nImp)
        enddo

        if(tConstrainKSym) then
            !We have only filled up half the states. Fill the others.
            !Update this so that it uses a function to map to equivalent kpoints
            !Currently this is only going to work for 1D

            if(tShift_Mesh) then
                !No gamma point sampled. All k-points symmetric.
                !Mirror the k-space hamiltonian
                do i = 1,ks
                    KBlocks(:,:,i+ks) = KBlocks(:,:,ks-i+1)
                enddo
            else
                !Mirror the kpoints, but ignore the gamma point and BZ boundary
                do i = 2,ks-1
                    KBlocks(:,:,i+ks-1) = KBlocks(:,:,ks-i+1)
                enddo
            endif
        endif

    end subroutine LatParams_to_KBlocks

    !This will set FreqPoints and Weights. nFitPoints is the size of FreqPoints/Weights, 
    !and may refer to real or im axis based on the tFitRealFreq flag.
    subroutine SetFreqPoints(nFreq_Re,nFreq_Im,nFitPoints)
        implicit none
        real(dp), intent(out) :: nFreq_Re,nFreq_Im
        integer, intent(out) :: nFitPoints
        logical :: tUsePoswFreqPoints

        if(tFitRealFreq.and.tFitPoints_Legendre) then
            call stop_all(t_r,'Trying to fit the real frequency spectrum to Legendre points?')
        endif
        if(.not.tHalfFill.and.tUsePoswFreqPoints) then
            call stop_all(t_r,'Trying to just optimize the positive frequency points, but not half filled system')
        endif
        if(tHalfFill) then
            tUsePoswFreqPoints = .true. !For ph symmetric systems, this seems to make no difference (but is faster!)
        endif

        nFreq_Re = 0
        OmegaVal = 0
        do while(.true.)
            call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.false.)
            if(OmegaVal.lt.0) exit
            nFreq_Re = nFreq_Re + 1
        enddo
        
        if(tFitPoints_Legendre) then
            !Calculate the Legendre polynomials for Legendre integration points and weights
            allocate(FreqPoints(nFreqPoints))
            allocate(IntWeights(nFreqPoints))
            write(6,"(A)") "Calculating Legendre frequency points..."
            LowFreq = -1.0_dp
            HighFreq = 1.0_dp
            call gauleg(LowFreq,HighFreq,FreqPoints,IntWeights,nFreqPoints)
            call ScaleFreqPointsWeights(nFreqPoints,FreqPoints,IntWeights)
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
                        Weights_temp(j) = IntWeights(i)
                    endif
                enddo
                deallocate(FreqPoints,IntWeights)
                nFreqPoints = j
                allocate(FreqPoints(nFreqPoints))
                allocate(IntWeights(nFreqPoints))
                FreqPoints(:) = FreqPoints_temp(:)
                IntWeights(:) = Weights_temp(:)
                deallocate(FreqPoints_temp,Weights_temp)
                write(6,"(A,I7)") "New number of frequency points: ",nFreqPoints
            endif
            nFreq_Im = nFreqPoints
            write(6,"(A)") "Frequency integration points and weights: "
            do i = 1,nFreqPoints
                write(6,"(I9,2G20.13)") i,FreqPoints(i),IntWeights(i)
            enddo
        else
            nFreq_Im = 0
            OmegaVal = 0
            do while(.true.)
                call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.true.)
                if(OmegaVal.lt.0) exit
                nFreq_Im = nFreq_Im + 1
            enddo
        endif
        call flush(6)

        !Which axis do we want to do the fitting on?
        if(tFitRealFreq) then
            nFitPoints = nFreq_Re   
            tCalcRealSpectrum = .false.
            tFitMatAxis = .false.
        else
            nFitPoints = nFreq_Im
            tFitMatAxis = .true.
        endif
        
        !Initially, just see if we can fit the two different Matsubara spectral functions
        write(6,*) "Number of frequency points to fit across: ",nFitPoints

    end subroutine SetFreqPoints


end module SelfConsistentUtils
