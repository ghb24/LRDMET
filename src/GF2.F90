module GF2
    use const
    use timing
    use errors
    use globals
    use utils, only: get_free_unit,append_ext_real,append_ext
    use GF2Data
    use matrixops, only: mat_inv
    use mat_tools, only: DiagOneeOp
    use writedata
    implicit none

    contains

    !Input: Max iterations, Convergence threshold
    !       Initial Self-energy
    !       Whether Initial self-energy is fixed
    !       Damping of self-energy
    !       Initial Chemical potential
    !       Damping of density matrix
    subroutine GF2_Hub(MaxSEIter,InitialSE_GV,InitChemPot)
        implicit none
        integer, intent(in) :: MaxSEIter
        real(dp), intent(in), optional :: InitialSE_GV(:,:,:)
        real(dp), optional :: InitChemPot
        complex(dp), allocatable :: PreviousSE(:,:,:)
        real(dp), allocatable :: PreviousP(:,:),IterStats(:,:)
        real(dp) :: LatChemPot,nExcessElec,MuThresh,GMEnergy,PreviousEnergy
        real(dp) :: MF_ChemPot  !Mean-field approximation to chemical potential
        integer :: Iter
        real(dp), parameter :: SEThresh = 1.0e-7_dp
        character(len=*), parameter :: t_r='GF2_Hub'

        write(6,"(A)") ""
        write(6,"(A)") " RUNNING GF2 CALCULATION "
        write(6,"(A)") ""

        !Set up tau grids, and matsubara grids
        call SetupGrids()
        call AllocateMem_GV()

        !Set initial self-energy if necessary
        if(present(InitialSE_GV)) then
            if(size(InitialSE_GV,3).ne.nMatsubara) then
                call stop_all(t_r,'Input initial self-energy wrong size')
            endif
            SE_Matsu_GV(:,:,:) = InitialSE_GV(:,:,:)
        endif
        
        !Construct the fock matrix and density matrix
        !Returned in globals "FockMat_GV" and "DensityMat_MF_GV"
        !Mean-field chemical potential returned in GuessChemPot
        call GetFockandP(MF_ChemPot)
        
        if(present(InitChemPot)) then
            LatChemPot = InitChemPot
        else
            LatChemPot = MF_ChemPot   
        endif
        write(6,*) "Chemical potential initially set to: ",LatChemPot
        
        !This function will also build the global GF and correlated density,
        !as well as the global number of electrons in nElec_GV
        nExcessElec = ExcessElec(LatChemPot) 
        if(abs(nElec_GV-(nExcessElec+real(NEl,dp))).gt.1.0e-7_dp) then
            call stop_all(t_r,'Error in electron number')
        endif
        write(6,*) "Number of electrons in initial greens function: ",nElec_GV

        call CalcGMEnergy(GLat_Matsu_GV,SE_Matsu_GV,DensityMat_GV,GMEnergy)

        PreviousEnergy = GMEnergy
        allocate(PreviousSE(nSites,nSites,nMatsubara))
        allocate(PreviousP(nSites,nSites))
        PreviousSE(:,:,:) = SE_Matsu_GV(:,:,:)
        PreviousP(:,:) = DensityMat_GV(:,:)

        allocate(IterStats(6,0:MaxSEIter))
        IterStats(:,:) = zero

        write(6,"(A)") "Converging Self-Energy at the level of GF2..."

        !Convergence threshold on chemical potential microiterations
        MuThresh = 1.0e-6_dp
        Iter = 0

        call FindSEChanges(IterStats,Iter,LatChemPot,GMEnergy,PreviousEnergy,   &
            PreviousSE,PreviousP,MaxSEIter)

        do while((IterStats(3,Iter).gt.SEThresh).or.(IterStats(5,Iter).gt.SEThresh)    &
                .or.(IterStats(6,Iter).gt.SEThresh).or.(Iter.le.0))
            !Converge self-energy loop
            call WriteGF2Stats(Iter,IterStats,MaxSEIter)
            Iter = Iter + 1
            if(Iter.gt.MaxSEIter) call stop_all(t_r,'Maximum global iterations hit')

            !Converge global chemical potential and Fock matrix. Returns new chemical
            !potential, and global matsubara greens function and fock matrix st.
            !electron number is correct and consistent.
            call ConvergeChemPotAndFock(LatChemPot,MuThresh)

            !FFT Greens function from iw -> tau
            call FT_GF_MatsuToImTime(GLat_Matsu_GV,GLat_Tau_GV,.true.) 

            !Build the self-energy in tau space
            call Build_SelfEnergy(GLat_Tau_GV,SE_Tau_GV)

            !FFT Self-energy from tau space to iw
            call FT_GF_ImTimeToMatsu(SE_Tau_GV,SE_Matsu_GV,.false.)

            !Build G(iw) and density
            nExcessElec = ExcessElec(LatChemPot)
        
            !Calculate energy according to Galitskii-Migdal formula
            call CalcGMEnergy(GLat_Matsu_GV,SE_Matsu_GV,DensityMat_GV,GMEnergy)

            !Test for convergence of self energy, or if hit iteration limit
            call FindSEChanges(IterStats,Iter,LatChemPot,GMEnergy,  &
                PreviousEnergy,PreviousSE,PreviousP,MaxSEIter)

        enddo

        !Write final energy and stats
        call WriteGF2Stats(Iter,IterStats,MaxSEIter)
        write(6,"(A,I6,A)") "GF2 convergence complete! Global self-energy converged in ",Iter," iterations."

        write(6,"(A,F17.10)") "Final Energy: ",GMEnergy
        write(6,"(A,F17.10)") "Final Chemical Potential: ",LatChemPot

        !TODO: Write Self-energy(iw) and restart information

        !TODO: Analytically continue G(iw) -> G(w) with Pade (and SE)

        !TODO Write G(w)

        deallocate(PreviousSE,PreviousP,IterStats)

    end subroutine GF2_Hub

    subroutine CalcGMEnergy(G_iw,SE_iw,P,Energy)
        implicit none
        complex(dp), intent(in) :: G_iw(nSites,nSites,nMatsubara)
        complex(dp), intent(in) :: SE_iw(nSites,nSites,nMatsubara)
        real(dp), intent(in) :: P(nSites,nSites)
        real(dp), intent(out) :: Energy
        integer :: i,j,k

        Energy = zero
        !One-electron contribution
        do i = 1,nSites
            do j = 1,nSites
                Energy = Energy + P(j,i)*(h0(j,i) + 0.5_dp*real(SE_iw(j,i,nMatsubara),dp))
            enddo
        enddo

        !Two-electron contribution (sum over frequencies)
        do i = 1,nMatsubara
            do j = 1,nSites
                do k = 1,nSites
                    Energy = Energy + (2.0_dp/Beta_Temp)*(real(G_iw(k,j,i))*real(SE_iw(k,j,i)) -    &
                        aimag(G_iw(k,j,i))*aimag(SE_iw(k,j,i)))
                enddo
            enddo
        enddo

    end subroutine CalcGMEnergy

    subroutine WriteGF2Stats(Iter,IterStats,MaxIter)
        implicit none
        integer, intent(in) :: Iter,MaxIter
        real(dp), intent(in) :: IterStats(6,0:MaxIter)
        integer :: i

        write(6,"(A)") ""
        write(6,"(A)") "  GF2: Global convergence update "
        write(6,"(A)") "Iteration   No.Elec   Energy   Delta_E   ChemPot   Delta_SE   Delta_P"
        do i = 0,Iter
            write(6,"(I7,6F13.7)") Iter,IterStats(1,Iter),IterStats(2,Iter),IterStats(3,Iter),  &
                IterStats(4,Iter),IterStats(5,Iter),IterStats(6,Iter)
        enddo
        write(6,"(A)") ""

    end subroutine WriteGF2Stats

    subroutine FindSEChanges(IterStats,Iter,ChemPot,E,PreviousE,PreviousSE,PreviousP,MaxIter)
        implicit none
        integer, intent(in) :: MaxIter
        real(dp), intent(inout) :: IterStats(6,0:MaxIter)
        integer, intent(in) :: Iter
        real(dp), intent(in) :: ChemPot
        real(dp), intent(in) :: E
        real(dp), intent(inout) :: PreviousE
        complex(dp), intent(inout) :: PreviousSE(nSites,nSites,nMatsubara)
        real(dp), intent(inout) :: PreviousP(nSites,nSites)
        integer :: i,j,k

        IterStats(1,Iter) = nElec_GV
        IterStats(2,Iter) = E
        IterStats(3,Iter) = E - PreviousE
        IterStats(4,Iter) = ChemPot
        IterStats(5,Iter) = zero    !Delta_SE
        IterStats(6,Iter) = zero    !Delta_P

        do i = 1,nSites
            do j = 1,nSites
                IterStats(6,Iter) = IterStats(6,Iter) + &
                    abs(DensityMat_GV(j,i) - PreviousP(j,i))
            enddo
        enddo
        do i = 1,nMatsubara
            do j = 1,nSites
                do k = 1,nSites
                    IterStats(5,Iter) = IterStats(5,Iter) + &
                        abs(SE_Matsu_GV(k,j,i) - PreviousSE(k,j,i))
                enddo
            enddo
        enddo
        IterStats(5,Iter) = IterStats(5,Iter)/real(nMatsubara,dp)

        !Set previous values to current values
        PreviousE = E
        PreviousSE(:,:,:) = SE_Matsu_GV(:,:,:)
        PreviousP(:,:) = DensityMat_GV(:,:)

    end subroutine FindSEChanges

    !This will converge s.t. the chemical potential, fock matrix, and
    !GLat_Matsu_GV are all consistent, and give the correct number of electrons
    subroutine ConvergeChemPotAndFock(LatChemPot,MuThresh)
        implicit none
        real(dp), intent(inout) :: LatChemPot       !Guess potential (should be correct for 1st iter)
        real(dp), intent(in) :: MuThresh            !Tolerance on convergence of quantities
        real(dp), allocatable :: Previous_P_MF(:,:)  !Initial MF density matrix
        real(dp), allocatable :: Previous_P(:,:)     !Initial correlated density
        real(dp), allocatable :: P_Diag(:)
        real(dp) :: Previous_Chempot
        real(dp) :: DeltaChemPot,Delta_P_MF,DeltaP,DiffP,MF_ChemPot,nExcessElec
        integer :: MicroIt,i
        character(len=*), parameter :: t_r='ConvergeChemPot'

        allocate(Previous_P_MF(nSites,nSites))
        allocate(Previous_P(nSites,nSites))
        allocate(P_Diag(nSites))

        Previous_P_MF(:,:) = DensityMat_MF_GV(:,:)
        Previous_P(:,:) = DensityMat_GV(:,:)
        Previous_ChemPot = LatChemPot
        MicroIt = 0

        call FindDensityChanges(Previous_P_MF,Previous_P,Previous_ChemPot,LatChemPot,  &
            DeltaChemPot,Delta_P_MF,DeltaP,DiffP)

        !Write out change in chemical potential, fock matrix, correlated and MF
        !P, and also the difference between the correlated and MF P.
        write(6,"(A)") "Converging chemical potential and densities..."
        write(6,"(A)") "Microiteration   No.Elec   Chempot   Delta_ChemPot   Delta_P_MF   Delta_P   Diff_P"

        do while((Delta_P_MF.gt.MuThresh).or.(DeltaP.gt.MuThresh).or.   &
                (DeltaChemPot.gt.MuThresh).or.(MicroIt.le.0))
        
            write(6,"(I5,6F13.7)") MicroIt, nElec_GV, LatChemPot, DeltaChemPot, Delta_P_MF, DeltaP, DiffP
            MicroIt = MicroIt + 1
            if(MicroIt.gt.50) call stop_all(t_r,'Maximum iteration number hit')

            !Converge the chemical potential with fixed fock matrix
            call ConvergeChemPot(LatChemPot)

            !Build fock matrix from h0 and diagonal of (correlated) density matrix
            !Updates fock matrix, mean-field density
            do i = 1,nSites
                P_Diag(i) = DensityMat_GV(i,i)
            enddo
            !Returns the chemical potential from the mean-field density
            call GetFockandP(MF_ChemPot,GuessDensity=P_Diag)

            !What is the change in the quantities (update previous ones)
            call FindDensityChanges(Previous_P_MF,Previous_P,Previous_ChemPot,LatChemPot,  &
                DeltaChemPot,Delta_P_MF,DeltaP,DiffP)

        enddo
            
        write(6,"(I5,6F13.7)") MicroIt, nElec_GV, LatChemPot, DeltaChemPot, Delta_P_MF, DeltaP, DiffP
        write(6,"(A)") "*** Fock, MF and correlated densities, and chemical potential converged ***"
    
        !Final calculation of greens function, density and number of electrons
        nExcessElec = ExcessElec(LatChemPot)

        deallocate(Previous_P_MF,Previous_P,P_Diag)

    end subroutine ConvergeChemPotAndFock

    !Find chemical potential s.t. the number of electrons in the lattice GF (for
    !fixed fock matrix) is correct
    subroutine ConvergeChemPot(LatChemPot)
        implicit none
        real(dp), intent(inout) :: LatChemPot
        real(dp) :: ChemPotLow,ChemPotHigh,OptChemPot
        real(dp), parameter :: ChemPotTol = 1.0e-7_dp

        !First, bracket the chemical potential between ChemPotLow and
        !ChemPotHigh
        call BracketChemPot(LatChemPot,ChemPotLow,ChemPotHigh)

        !Brent algorithm to converge the 1D root finding exercise of converging
        !the chemical potential
        OptChemPot = zbrent(ExcessElec,ChemPotLow,ChemPotHigh,ChemPotTol)
        !Update chemical potential
        LatChemPot = OptChemPot

    end subroutine ConvergeChemPot

    subroutine FT_GF_ImTimeToMatsu(GF_Tau,GF_iw,tGF)
        implicit none
        complex(dp), intent(in) :: GF_Tau(nSites,nSites,nImTimePoints)
        complex(dp), intent(out) :: GF_iw(nSites,nSites,nMatsubara)
        logical, intent(in) :: tGF
        integer :: i,j
        real(dp) :: Delta_Tau
        logical :: tDum

        tDum = tGF

        !Initially, do brute force - no splining
        Delta_Tau = Beta_Temp/(nImTimePoints-1)
        do i = 1,nMatsubara
            GF_iw(:,:,i) = zzero
            
            do j = 2,nImTimePoints-1
                GF_iw(:,:,i) = GF_iw(:,:,i) + exp(cmplx(zero,MatsuPoints(i)*ImTimePoints(j),dp))*GF_Tau(:,:,j)
            enddo
            GF_iw(:,:,i) = GF_iw(:,:,i) + (GF_Tau(:,:,1) - GF_Tau(:,:,nImTimePoints))/2.0_dp
            GF_iw(:,:,i) = Delta_Tau*GF_iw(:,:,i)
        enddo

    end subroutine FT_GF_ImTimeToMatsu

    !Build the self energy in imaginary-time
    subroutine Build_SelfEnergy(GF_Tau,SE_Tau)
        implicit none
        complex(dp), intent(in) :: GF_Tau(nSites,nSites,nImTimePoints)
        complex(dp), intent(out) :: SE_Tau(nSites,nSites,nImTimePoints)
        integer :: i,j,k,nImOpp
        character(len=*), parameter :: t_r='Build_SelfEnergy'

        SE_Tau(:,:,:) = zzero
        do i = 1,nImTimePoints
            nImOpp = nImTimePoints - i + 1
            if(abs(ImTimePoints(nImOpp) - (-ImTimePoints(i)+Beta_Temp)).gt.1.0e-7_dp) then
                call stop_all(t_r,'-tau not correctly sampled')
            endif

            do j = 1,nSites
                do k = 1,nSites

                    SE_Tau(k,j,i) = GF_Tau(k,j,i)*GF_Tau(k,j,i)*GF_Tau(j,k,nImOpp)*U*U

                enddo
            enddo
        enddo

    end subroutine Build_SelfEnergy

    subroutine FT_GF_MatsuToImTime(GF_iw,GF_tau,tGF)
        implicit none
        complex(dp), intent(in) :: GF_iw(nSites,nSites,nMatsubara)
        complex(dp), intent(out) :: GF_tau(nSites,nSites,nImTimePoints)
        logical, intent(in) :: tGF    !GF or SE?
        integer :: i
        character(len=*), parameter :: t_r='FT_GF_MatsuToImTime'

        write(6,"(A)") "Fourier transforming matsubara greens function to imaginary time..."

        if(.not.tGF) then
            call stop_all(t_r,'Unsure how to FT SE into imaginary time')
        endif

        do i = 1,nImTimePoints
            call MatsuToImTimeFT(GF_iw,ImTimePoints(i),GF_tau(:,:,i),tGF)
        enddo

    end subroutine FT_GF_MatsuToImTime

    !Fourier transform the antisymmetric Matsubara axis function GF_Matsu at tau = TauPoint
    subroutine MatsuToImTimeFT(GF_Matsu,TauPoint,GF_Tau,tGF)
        implicit none
        complex(dp), intent(in) :: GF_Matsu(nSites,nSites,nMatsubara)
        real(dp), intent(in) :: TauPoint
        complex(dp), intent(out) :: GF_Tau(nSites,nSites)
        logical, intent(in) :: tGF  !Is this for a greens function?
        complex(dp), allocatable :: eye(:,:)
        integer :: i
        character(len=*), parameter :: t_r='MatsuToImTimeFT'

        if(.not.tGF) then
            !This assumes that the c1 coefficient for the 1/iw term is 1
            call stop_all(t_r,'Cannot FT from iw to tau for self-energies yet')
        endif

        allocate(eye(nSites,nSites))
        eye(:,:) = zzero
        do i = 1,nSites
            eye(i,i) = zone
        enddo

        GF_Tau(:,:) = zzero
        !The simplest way first, removing the 1/iw tail analytically (from the
        !diagonals only)
        do i = 1,nMatsubara
            GF_Tau(:,:) = GF_Tau +  &
                exp(-cmplx(zero,MatsuPoints(i)*TauPoint,dp)) * &
                (GF_Matsu(:,:,i) - eye(:,:)*(one/cmplx(zero,MatsuPoints(i),dp)) )
        enddo
        !Negative Frequencies
        do i = 1,nMatsubara
            GF_Tau(:,:) = GF_Tau +  &
                exp(cmplx(zero,MatsuPoints(i)*TauPoint,dp)) * &
                (-GF_Matsu(:,:,i) + eye(:,:)*(one/cmplx(zero,MatsuPoints(i),dp)) )
        enddo

        GF_Tau(:,:) = GF_Tau(:,:) / Beta_Temp

        !Add back on the analytically FT'ed tail
        GF_Tau(:,:) = GF_Tau(:,:) - 0.5_dp*eye(:,:)

        deallocate(eye)

    end subroutine MatsuToImTimeFT

    !P = -2 x G(tau = Beta)
    !Beta will always be the last imaginary time point 
    subroutine GetGSDensityFromMatsuGF(GF_Matsu)
        implicit none
        complex(dp), optional :: GF_Matsu(nSites,nSites,nMatsubara)
        complex(dp), allocatable :: TempMat(:,:)
        integer :: i,j
        character(len=*), parameter :: t_r='GetGSDensityFromMatsuGF'

        allocate(TempMat(nSites,nSites))

        if(abs(ImTimePoints(nImTimePoints)-Beta_Temp).gt.1.0e-8_dp) then
            call stop_all(t_r,'Error in grids')
        endif
        if(present(GF_Matsu)) then
            call MatsuToImTimeFT(GF_Matsu,ImTimePoints(nImTimePoints),TempMat,.true.)
        else
            call MatsuToImTimeFT(GLat_Matsu_GV,ImTimePoints(nImTimePoints),TempMat,.true.)
        endif

        do i = 1,nSites
            do j = i,nSites
                if(abs(aimag(TempMat(i,j))).gt.1.0e-7_dp) then
                    !Though perhaps this is ok other than the diagonals?
                    call stop_all(t_r,'Density matrix complex')
                elseif(abs(real(TempMat(i,j))-real(TempMat(j,i))).gt.1.0e-7_dp) then
                    call stop_all(t_r,'Density matrix not hermitian')
                endif
            enddo
        enddo

        DensityMat_GV(:,:) = -2.0_dp*real(TempMat(:,:),dp)
        deallocate(TempMat)

    end subroutine GetGSDensityFromMatsuGF

    !Get mean field density matrix and fock matrix
    subroutine GetFockandP(GuessChemPot,GuessDensity)
        implicit none
        real(dp), intent(out) :: GuessChemPot
        real(dp), allocatable :: EigenVals(:)
        real(dp), intent(in), optional :: GuessDensity(nSites)
        real(dp) :: InputDensity(nSites)
        real(dp), allocatable :: EigenVecs(:,:),OccOrbs(:,:)
        integer :: i

        FockMat_GV(:,:) = h0(:,:)
        if(present(GuessDensity)) then
            InputDensity(:) = GuessDensity(:)
        else
            InputDensity(:) = real(NEl,dp)/real(nSites,dp)
        endif
        do i = 1,nSites
            FockMat_GV(i,i) = FockMat_GV(i,i) + U*InputDensity(i)/2.0_dp
        enddo
        allocate(EigenVecs(nSites,nSites))
        allocate(EigenVals(nSites))
        EigenVecs(:,:) = FockMat_GV(:,:)
        call DiagOneEOp(EigenVecs,EigenVals,1,nSites,.false.,.true.)

        !Calculate Density
        allocate(OccOrbs(nSites,nOcc))
        OccOrbs(:,:) = EigenVecs(:,1:nOcc)
        call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,OccOrbs,nSites,OccOrbs,nSites,zero,DensityMat_MF_GV,nSites)

        GuessChemPot = (EigenVals(nOcc) + EigenVals(nOcc+1))/2.0_dp
        deallocate(OccOrbs,EigenVals,EigenVecs)

    end subroutine GetFockandP

    !Setup the imaginary-time and matsubara grids
    !Input parameters:  No. Matsubara points
    !                   Temperature
    !                   
    !The number of time points is hardcoded to be nFacTau x nMatsu * beta / pi
    !nFacTau is hardcoded initially to be 20.
    
    !The grids are returned in:
    !   MatsuPoints(1:nMatsubara)       [All positive frequencies only]
    !   ImTimePoints(1:nImTimePoints)   [1 is 0 time and nImTimePoints is Beta]

    !TODO:  Ultimately, we don't want a uniform grid in imaginary time
    !       We want to sample low and high times more
    subroutine SetupGrids()
        implicit none
        integer :: i
        real(dp) :: Delta_Tau

        !Matsubara grids are defined by (2n+1)pi/beta for n = 0 -> infty
        !MatsuPoints(0) = beta/pi
        if(allocated(MatsuPoints)) then
            write(6,"(A)") "Matsubara grid already set up. Resetting up..."
            deallocate(MatsuPoints)
        endif
        if(allocated(ImTimePoints)) then
            write(6,"(A)") "Tau grid already set up. Resetting up..."
            deallocate(ImTimePoints)
        endif
        allocate(MatsuPoints(nMatsubara))
        do i = 0,nMatsubara-1
            MatsuPoints(i+1) = (2*i+1)*pi / Beta_Temp
        enddo
        write(6,"(A,I7)") "Number of matsubara frequency points set to: ",nMatsubara
        write(6,"(A,F20.10)") "Highest matsubara frequency point: ",MatsuPoints(nMatsubara)

        !The number of Im Time points is controlled by ScaleImTime which
        !oversamples to avoid aliasing
        nImTimePoints = ScaleImTime * real(nMatsubara,dp) * Beta_Temp / pi
        allocate(ImTimePoints(nImTimePoints))
        write(6,"(A,I7)") "Number of imag-time points set to: ",nImTimePoints

        Delta_Tau = Beta_Temp/(nImTimePoints-1)
        write(6,"(A,G16.5)") "Imaginary time step in grid: ",Delta_Tau
        ImTimePoints(:) = zero
        do i = 2,nImTimePoints
            ImTimePoints(i) = (i-1)*Delta_Tau
        enddo

    end subroutine SetupGrids

    !TODO: Everything should be built and stored in kspace
    subroutine AllocateMem_GV()
        implicit none
        integer :: nCompMem

        if(allocated(SE_Matsu_GV)) deallocate(SE_Matsu_GV)
        if(allocated(SE_Tau_GV)) deallocate(SE_Tau_GV)
        if(allocated(GLat_Matsu_GV)) deallocate(GLat_Matsu_GV)
        if(allocated(GLat_Tau_GV)) deallocate(GLat_Tau_GV)

        ! x 3 because we allocate a 'previous' SE on the matsubara grid in a bit
        nCompMem = (nSites**2)*(nMatsubara*3+nImTimePoints*2)
        write(6,"(A,F10.3,A)") "Total memory required for SE and GFs: ",nCompMem*ComptoGb," Gb"

        allocate(SE_Matsu_GV(nSites,nSites,nMatsubara))
        allocate(SE_Tau_GV(nSites,nSites,nImTimePoints))
        allocate(GLat_Matsu_GV(nSites,nSites,nMatsubara))
        allocate(GLat_Tau_GV(nSites,nSites,nImTimePoints))

        SE_Matsu_GV(:,:,:) = zzero
        GLat_Matsu_GV(:,:,:) = zzero

        if(allocated(DensityMat_GV)) deallocate(DensityMat_GV)
        if(allocated(DensityMat_MF_GV)) deallocate(DensityMat_MF_GV)
        if(allocated(FockMat_GV)) deallocate(FockMat_GV)

        allocate(DensityMat_GV(nSites,nSites))
        allocate(DensityMat_MF_GV(nSites,nSites))
        allocate(FockMat_GV(nSites,nSites))

    end subroutine AllocateMem_GV

    !TODO: Build this in k-space
    subroutine BuildMatsubaraGF(LatChemPot)
        implicit none
        real(dp), intent(in) :: LatChemPot
        complex(dp), allocatable :: SingleFreqMat(:,:)
        integer :: i,j

        GLat_Matsu_GV(:,:,:) = zzero
!$OMP PARALLEL DO PRIVATE(j,SingleFreqMat)
        do i = 1,nMatsubara
            allocate(SingleFreqMat(nSites,nSites))
            !Build matrix
            SingleFreqMat(:,:) = - FockMat_GV(:,:) - SE_Matsu_GV(:,:,i)
            do j = 1,nSites
                SingleFreqMat(j,j) = SingleFreqMat(j,j) + cmplx(LatChemPot,MatsuPoints(i),dp)
            enddo
            !Invert
            call mat_inv(SingleFreqMat,GLat_Matsu_GV(:,:,i),nSites)
            deallocate(SingleFreqMat)
        enddo
!$OMP END PARALLEL DO
    end subroutine BuildMatsubaraGF
    
    !This function returns the overpopulation of the system with electrons <N> - N_desired,
    !for a given chemical potential
    function ExcessElec(ChemPot) result(res)
        implicit none
        real(dp), intent(in) :: ChemPot
        real(dp) :: res
        integer :: i

        !Build the Matsubara greens function
        call BuildMatsubaraGF(ChemPot)

        !First, find the density matrix from the Matsubara greens function
        call GetGSDensityFromMatsuGF()

        !Then find the number of electrons from the density
        nElec_GV = zero
        do i = 1,nSites
            nElec_GV = nElec_GV + DensityMat_GV(i,i)
        enddo

        res = nElec_GV - real(NEl,dp)

    end function ExcessElec
        
    !Find convergence criteria on the densities and chemical potential, and
    !then update the matrices
    subroutine FindDensityChanges(PreviousP_MF,PreviousP,Previous_ChemPot,LatChemPot,DeltaMu,DeltaPMF,DeltaP,DiffP)
        implicit none
        real(dp), intent(inout) :: PreviousP_MF(nSites,nSites)
        real(dp), intent(inout) :: PreviousP(nSites,nSites)
        real(dp), intent(inout) :: Previous_ChemPot
        real(dp), intent(in) :: LatChemPot
        real(dp), intent(out) :: DeltaMu,DeltaPMF,DeltaP,DiffP
        integer :: i,j

        DeltaPMF = zero
        DeltaP = zero
        DiffP = zero
        do i = 1,nSites
            do j = 1,nSites
                DeltaPMF = DeltaPMF + abs(DensityMat_MF_GV(j,i) - PreviousP_MF(j,i))
                DeltaP = DeltaP + abs(DensityMat_GV(j,i) - PreviousP(j,i))
                DiffP = DiffP + abs(DensityMat_GV(j,i) - DensityMat_MF_GV(j,i))
            enddo
        enddo
        
        DeltaMu = LatChemPot - Previous_ChemPot

        Previous_ChemPot = LatChemPot
        PreviousP_MF(:,:) = DensityMat_MF_GV(:,:)
        PreviousP(:,:) = DensityMat_GV(:,:)

    end subroutine FindDensityChanges

    
    !For given fock matrix, find a chemical potential range which brackets the
    !desired number of electrons
    subroutine BracketChemPot(LatChemPot,ChemPotLow,ChemPotHigh)
        implicit none
        real(dp), intent(in) :: LatChemPot
        real(dp), intent(out) :: ChemPotLow,ChemPotHigh
        real(dp) :: x1,x2,f1,f2
        integer :: i
        real(dp), parameter :: fac = 1.5_dp
        character(len=*), parameter :: t_r='BracketChemPot'

        !Set initial guess of 0.05t higher and lower than guess
        ChemPotLow = LatChemPot - 0.05
        ChemPotHigh = LatChemPot + 0.05

        x1 = ChemPotLow
        x2 = ChemPotHigh

        f1 = ExcessElec(x1)
        f2 = ExcessElec(x2)
        !write(6,*) "ChemPot1, f1: ",x1,f1
        !write(6,*) "ChemPot2, f2: ",x2,f2
        do i = 1,50
            if(f1*f2.lt.zero) then
                ChemPotLow = x1
                ChemPotHigh = x2
                return
            endif
            if(abs(f1).lt.abs(f2)) then
                x1 = x1+fac*(x1-x2)
                f1 = ExcessElec(x1)
            else
                x2 = x2+fac*(x2-x1)
                f2 = ExcessElec(x2)
            endif
        enddo
        call stop_all(t_r,'Couldnt manage to bracket the chemical potential')

    end subroutine BracketChemPot

    !zbrent from numerical recipies
    function zbrent(func,x1,x2,tol)
        implicit none
        real(dp), intent(in) :: x1,x2,tol
        real(dp) :: zbrent
        interface
            function func(x)
                use const, only: dp 
                implicit none
                real(dp), intent(in) :: x
                real(dp) :: func
            end function func
        end interface
        integer, parameter :: ITMAX=100
        real(dp), parameter :: EPS=epsilon(x1)
        !Using Brent’s method, find the root of a function
        !func known to lie between x1 and x2.
        !The root, returned as zbrent, will be refined until its accuracy is tol
        !Parameters: Maximum allowed number of iterations, and machine floating-point precision.
        integer :: iter
        real(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        character(len=*), parameter :: t_r='zbrent'
        a=x1
        b=x2
        fa=func(a)
        fb=func(b)
        if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
            call stop_all(t_r,'Root must be bracketed for zbrent')
        c=b
        fc=fb
        do iter=1,ITMAX
            if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
                c=a
                !Rename a, b, c and adjust bounding interval d
                fc=fa
                d=b-a
                e=d
            end if
            if (abs(fc) < abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            end if
            tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol
            !Convergence check.
            xm=0.5_sp*(c-b)
            if (abs(xm) <= tol1 .or. fb == 0.0) then
                zbrent=b
                RETURN
            end if
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                s=fb/fa
                !Attempt inverse quadratic interpolation.
                if (a == c) then
                    p=2.0_sp*xm*s
                    q=1.0_sp-s
                else
                    q=fa/fc
                    r=fb/fc
                    p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
                    q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
                end if
                if (p > 0.0) q=-q
                !Check whether in bounds.
                p=abs(p)
                if (2.0_sp*p < min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
                    e=d
                    !Accept interpolation.
                    d=p/q
                else
                    d=xm
                    !Interpolation failed; use bisection.
                    e=d
                end if
            else
                !Bounds decreasing too slowly; use bisection
                d=xm
                e=d
            end if
            a=b
            !Move last best guess to a.
            fa=fb
            b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
            !Evaluate new trial root.
            fb=func(b)
        end do
        call stop_all(t_r,'Exceeded maximum iterations')
        zbrent=b

    end function zbrent

end module GF2
