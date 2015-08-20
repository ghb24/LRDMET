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
        !Finds the uncorrelated (t0) tails, and chemical potential
        if(tFitTails) then
            allocate(C2_Coeffs_GF0(nSites,nSites))
            allocate(C3_Coeffs_GF0(nSites,nSites))
        endif
        call FindG0ChemPot(GF0_ChemPot,1.0e-8_dp)
        call AllocateMem_GV()

        !Set initial self-energy if necessary
        if(present(InitialSE_GV)) then
            if(size(InitialSE_GV,3).ne.nMatsubara) then
                call stop_all(t_r,'Input initial self-energy wrong size')
            endif
            SE_GV%Matsu(:,:,:) = InitialSE_GV(:,:,:)
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
        nExcessElec = ExcessElec(GLat_GV,LatChemPot,.true.) 
        if(abs(nElec_GV-(nExcessElec+real(NEl,dp))).gt.1.0e-7_dp) then
            call stop_all(t_r,'Error in electron number')
        endif
        write(6,*) "Number of electrons in initial greens function: ",nElec_GV

        call CalcGMEnergy(GLat_GV,SE_GV,DensityMat_GV,GMEnergy)

        PreviousEnergy = GMEnergy
        allocate(PreviousSE(nSites,nSites,nMatsubara))
        allocate(PreviousP(nSites,nSites))
        PreviousSE(:,:,:) = SE_GV%Matsu(:,:,:)
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
            call FT_GF_MatsuToImTime(GLat_GV) 

            !Build the self-energy in tau space
            call Build_SelfEnergy(GLat_GV,SE_GV)

            !FFT Self-energy from tau space to iw
            call FT_GF_ImTimeToMatsu(SE_GV)

            !Build G(iw) and density
            nExcessElec = ExcessElec(GLat_GV,LatChemPot,.true.)
        
            !Calculate energy according to Galitskii-Migdal formula
            call CalcGMEnergy(GLat_GV,SE_GV,DensityMat_GV,GMEnergy)

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

        !Deallocate memory
        deallocate(PreviousSE,PreviousP,IterStats)

    end subroutine GF2_Hub

    !TODO: Use tails for this
    subroutine CalcGMEnergy(GF,SE,P,Energy)
        implicit none
        type(GreensFunc), intent(in) :: GF
        type(GreensFunc), intent(inout) :: SE
        real(dp), intent(in) :: P(nSites,nSites)
        real(dp), intent(out) :: Energy
        complex(dp), allocatable :: SE_tail(:,:),GF_tail(:,:)
        integer :: i,j,k,i_end
        real(dp) :: w

        Energy = zero
        if(tFitTails) then
            !Sigma_infty = C2 - C2_0
            SE%C0_Coeffs(:,:) = GF%C2_Coeffs(:,:) - C2_Coeffs_GF0(:,:)
            !One-electron contribution
            do i = 1,nSites
                do j = 1,nSites
                    Energy = Energy + P(j,i)*(h0(j,i) + 0.5_dp*SE%C0_Coeffs(j,i))
                    !write(6,*) SE%C0_Coeffs(j,i) - real(SE%Matsu(j,i,nMatsubara),dp),SE%C0_Coeffs(j,i)
                enddo
            enddo
        else
            !One-electron contribution
            do i = 1,nSites
                do j = 1,nSites
                    Energy = Energy + P(j,i)*(h0(j,i) + 0.5_dp*real(SE%Matsu(j,i,nMatsubara),dp))
                enddo
            enddo
        endif

        !Two-electron contribution (sum over frequencies)
        do i = 1,nMatsubara
            do j = 1,nSites
                do k = 1,nSites
                    Energy = Energy + (2.0_dp/Beta_Temp)*(real(GF%Matsu(k,j,i))*real(SE%Matsu(k,j,i)) -    &
                        aimag(GF%Matsu(k,j,i))*aimag(SE%Matsu(k,j,i)))
                enddo
            enddo
        enddo

        if(tFitTails.and.(MatsuEnergySumFac.gt.1.0_dp)) then
            !Set C1 coefficients for SE
            SE%C1_Coeffs(:,:) = C2_Coeffs_GF0(:,:)**2 - GF%C2_Coeffs(:,:)**2 + &
                GF%C3_Coeffs(:,:) - C3_Coeffs_GF0(:,:)
            !Fit C2 and C3 coefficients for SE
            call FitGFTail(SE)

            !Sum additional contributions from long-range tails to energy
            allocate(GF_tail(nSites,nSites))
            allocate(SE_tail(nSites,nSites))
            i_end = nint(MatsuEnergySumFac*nMatsubara/2.0_dp)
            !Positive Frequencies
            do i = (nMatsubara/2),i_end
                w = (2*i+1)*pi / Beta_Temp
                call EvalTails(w,GF,SE,GF_tail,SE_tail)
                do j = 1,nSites
                    do k = 1,nSites
                        Energy = Energy + (2.0_dp/Beta_Temp)*(real(GF_tail(k,j))*real(SE_tail(k,j)) -    &
                            aimag(GF_tail(k,j))*aimag(SE_tail(k,j)))
                    enddo
                enddo
            enddo
            deallocate(SE_tail,GF_tail)
        endif

    end subroutine CalcGMEnergy

    pure subroutine EvalTails(w,GF,SE,GF_tail,SE_tail)
        implicit none
        real(dp), intent(in) :: w
        type(GreensFunc), intent(in) :: GF,SE
        complex(dp), intent(out) :: GF_tail(nSites,nSites)
        complex(dp), intent(out) :: SE_tail(nSites,nSites)
        integer :: i

        GF_tail(:,:) = zzero
        do i = 1,nSites
            GF_tail(i,i) = -cmplx(zero,w,dp)
        enddo
        !C2 and C3
        GF_tail(:,:) = GF_tail(:,:) - GF%C2_Coeffs(:,:)/cmplx(w**2,zero,dp) + &
            GF%C3_Coeffs(:,:)*cmplx(zero,one/w**3,dp)

        !SE
        SE_tail(:,:) = SE%C0_Coeffs(:,:) - SE%C1_Coeffs(:,:)*cmplx(zero,one/w,dp) - &
            SE%C2_Coeffs(:,:)/cmplx(w**2,zero,dp) + &
            SE%C3_Coeffs(:,:)*cmplx(zero,one/w**3,dp)

    end subroutine EvalTails

    subroutine WriteGF2Stats(Iter,IterStats,MaxIter)
        implicit none
        integer, intent(in) :: Iter,MaxIter
        real(dp), intent(in) :: IterStats(6,0:MaxIter)
        integer :: i

        write(6,"(A)") ""
        write(6,"(A)") "  GF2: Global convergence update "
        write(6,"(A)") "Iteration   No.Elec      Energy       Delta_E      ChemPot      Delta_SE     Delta_P"
        do i = 0,Iter
            write(6,"(I7,6F13.7)") i,IterStats(1,i),IterStats(2,i),IterStats(3,i),  &
                IterStats(4,i),IterStats(5,i),IterStats(6,i)
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
                        abs(SE_GV%Matsu(k,j,i) - PreviousSE(k,j,i))
                enddo
            enddo
        enddo
        IterStats(5,Iter) = IterStats(5,Iter)/real(nMatsubara,dp)

        !Set previous values to current values
        PreviousE = E
        PreviousSE(:,:,:) = SE_GV%Matsu(:,:,:)
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
        real(dp) :: ChemPotTol
        integer :: MicroIt,i
        character(len=*), parameter :: t_r='ConvergeChemPot'

        allocate(Previous_P_MF(nSites,nSites))
        allocate(Previous_P(nSites,nSites))
        allocate(P_Diag(nSites))

        Previous_P_MF(:,:) = DensityMat_MF_GV(:,:)
        Previous_P(:,:) = DensityMat_GV(:,:)
        Previous_ChemPot = LatChemPot
        MicroIt = 0
        ChemPotTol = MuThresh/5.0_dp

        call FindDensityChanges(Previous_P_MF,Previous_P,Previous_ChemPot,LatChemPot,  &
            DeltaChemPot,Delta_P_MF,DeltaP,DiffP)

        !Write out change in chemical potential, fock matrix, correlated and MF
        !P, and also the difference between the correlated and MF P.
        write(6,"(A)") "Converging chemical potential and densities..."
        write(6,"(A)") "MicroIter No.Elec      Chempot   Delta_ChemPot  Delta_P_MF    Delta_P      Diff_P"

        do while((Delta_P_MF.gt.MuThresh).or.(DeltaP.gt.MuThresh).or.   &
                (DeltaChemPot.gt.MuThresh).or.(MicroIt.le.0))
        
            write(6,"(I5,6F13.7)") MicroIt, nElec_GV, LatChemPot, DeltaChemPot, Delta_P_MF, DeltaP, DiffP
            MicroIt = MicroIt + 1
            if(MicroIt.gt.500) call stop_all(t_r,'Maximum iteration number hit')

            !Converge the chemical potential with fixed fock matrix
            call ConvergeChemPot(LatChemPot,GLat_GV,.true.,ChemPotTol)

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
        nExcessElec = ExcessElec(GLat_GV,LatChemPot,.true.)

        deallocate(Previous_P_MF,Previous_P,P_Diag)

    end subroutine ConvergeChemPotAndFock

    !Find the uncorrelated chemical potetial, and tails to the GF.
    !These will be used to constrain the tails of the self-energy.
    subroutine FindG0ChemPot(h0ChemPot,ChemPotTol)
        implicit none
        real(dp), intent(out) :: h0ChemPot
        real(dp), intent(in) :: ChemPotTol
        type(GreensFunc) :: GF0
        real(dp), allocatable :: EigenVecs(:,:), EigenVals(:)

        write(6,"(A)") "Converging non-interacting G0 greens function..."

        !Diagonalize h0 for guess chemical potential
        allocate(EigenVecs(nSites,nSites))
        allocate(EigenVals(nSites))
        EigenVecs(:,:) = h0(:,:)
        call DiagOneEOp(EigenVecs,EigenVals,1,nSites,.false.,.true.)
        h0ChemPot = (EigenVals(nOcc) + EigenVals(nOcc+1))/2.0_dp
        write(6,"(A,F15.7)") "Guess non-interacting chemical potential: ",h0ChemPot
        deallocate(EigenVals,EigenVecs)

        !Now refine, and obtain the tail coefficients
        call allocateGF(GF0,nSites,.true.)
        call ConvergeChemPot(h0ChemPot,GF0,.false.,ChemPotTol)
        write(6,"(A,F16.10)") "Converged non-interacting greens function chemical potential: ",h0ChemPot
        call BuildMatsubaraGF(GF0,h0ChemPot,.false.)
        if(tFitTails) then
            C2_Coeffs_GF0(:,:) = GF0%C2_Coeffs(:,:)
            C3_Coeffs_GF0(:,:) = GF0%C3_Coeffs(:,:)
            write(6,"(A)") "Tails of non-interacting greens function found."
        endif
        call DeallocateGF(GF0)

    end subroutine FindG0ChemPot

    !Find chemical potential s.t. the number of electrons in the lattice GF (for
    !fixed fock matrix) is correct
    subroutine ConvergeChemPot(LatChemPot,GF,tCorr,ChemPotTol)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        logical, intent(in) :: tCorr
        real(dp), intent(inout) :: LatChemPot
        real(dp), intent(in) :: ChemPotTol 
        real(dp) :: ChemPotLow,ChemPotHigh,OptChemPot

        !First, bracket the chemical potential between ChemPotLow and
        !ChemPotHigh
        call BracketChemPot(ExcessElec,LatChemPot,ChemPotLow,ChemPotHigh,GF,tCorr)

        !Brent algorithm to converge the 1D root finding exercise of converging
        !the chemical potential
        OptChemPot = zbrent(ExcessElec,ChemPotLow,ChemPotHigh,ChemPotTol,GF,tCorr)
        !Update chemical potential
        LatChemPot = OptChemPot

    end subroutine ConvergeChemPot

    !TODO: Fit splines
    subroutine FT_GF_ImTimeToMatsu(GF)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        integer :: i,j
        real(dp) :: Delta_Tau

        !Initially, do brute force - no splining
        Delta_Tau = Beta_Temp/(nImTimePoints-1)
        do i = 1,nMatsubara
            GF%Matsu(:,:,i) = zzero
            
            do j = 2,nImTimePoints-1
                GF%Matsu(:,:,i) = GF%Matsu(:,:,i) + exp(cmplx(zero,MatsuPoints(i)*ImTimePoints(j),dp))*GF%Tau(:,:,j)
            enddo
            GF%Matsu(:,:,i) = GF%Matsu(:,:,i) + (GF%Tau(:,:,1) - GF%Tau(:,:,nImTimePoints))/2.0_dp
            GF%Matsu(:,:,i) = Delta_Tau*GF%Matsu(:,:,i)
        enddo

    end subroutine FT_GF_ImTimeToMatsu

    !Build the self energy in imaginary-time
    subroutine Build_SelfEnergy(GF,SE)
        implicit none
        type(GreensFunc), intent(in) :: GF
        type(GreensFunc), intent(inout) :: SE
        integer :: i,j,k,nImOpp
        character(len=*), parameter :: t_r='Build_SelfEnergy'

        SE%Tau(:,:,:) = zzero
        do i = 1,nImTimePoints
            nImOpp = nImTimePoints - i + 1
            if(abs(ImTimePoints(nImOpp) - (-ImTimePoints(i)+Beta_Temp)).gt.1.0e-7_dp) then
                call stop_all(t_r,'-tau not correctly sampled')
            endif

            do j = 1,nSites
                do k = 1,nSites

                    SE%Tau(k,j,i) = GF%Tau(k,j,i)*GF%Tau(k,j,i)*GF%Tau(j,k,nImOpp)*U*U

                enddo
            enddo
        enddo

    end subroutine Build_SelfEnergy

    subroutine FT_GF_MatsuToImTime(GF)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        complex(dp), allocatable :: GF_tau(:,:)
        integer :: i
        character(len=*), parameter :: t_r='FT_GF_MatsuToImTime'

        write(6,"(A)") "Fourier transforming matsubara greens function to imaginary time..."

        if(.not.GF%tGF) then
            call stop_all(t_r,'Unsure how to FT SE into imaginary time')
        endif

        allocate(GF_tau(nSites,nSites))
        do i = 1,nImTimePoints
            call MatsuToImTimeFT(GF,ImTimePoints(i),GF_tau(:,:))
            GF%Tau(:,:,i) = GF_tau(:,:)
        enddo
        deallocate(GF_tau)

    end subroutine FT_GF_MatsuToImTime
    
    !Remove the GF tail as specified in C2_Coeff and C3_Coeff for
    !a given matsubara frequency of the greens function
    pure function GFwoTail(G_iw,w,C2,C3) result(GF)
        implicit none
        real(dp), intent(in) :: w
        complex(dp), intent(in) :: G_iw(nSites,nSites)
        real(dp), intent(in) :: C2(nSites,nSites),C3(nSites,nSites)
        complex(dp) :: GF(nSites,nSites)
        integer :: i

        GF(:,:) = G_iw(:,:)
        !Remove C1
        do i = 1,nSites
            GF(i,i) = GF(i,i) + cmplx(zero,one/w,dp)
        enddo
        !Remove C2
        GF(:,:) = GF(:,:) + (zone/(w**2))*C2(:,:)
        !Remove C3
        GF(:,:) = GF(:,:) - cmplx(zero,one/(w**3),dp)*C3(:,:)

    end function GFwoTail

    !Fourier transform the antisymmetric Matsubara axis function GF_Matsu at tau = TauPoint
    subroutine MatsuToImTimeFT(GF,TauPoint,GF_Tau)
        implicit none
        type(GreensFunc), intent(in) :: GF
        real(dp), intent(in) :: TauPoint
        complex(dp), intent(out) :: GF_Tau(nSites,nSites)
        complex(dp), allocatable :: eye(:,:)
        integer :: i
        character(len=*), parameter :: t_r='MatsuToImTimeFT'

        if(.not.GF%tGF) then
            !This assumes that the c1 coefficient for the 1/iw term is 1
            call stop_all(t_r,'Cannot FT from iw to tau for self-energies yet')
        endif

        GF_Tau(:,:) = zzero
        if(tFitTails) then
            do i = 1,nMatsubara
                GF_Tau(:,:) = GF_Tau(:,:) +  &
                    exp(-cmplx(zero,MatsuPoints(i)*TauPoint,dp)) * &
                    GFwoTail(GF%Matsu(:,:,i),MatsuPoints(i),GF%C2_Coeffs,GF%C3_Coeffs)
            enddo
            GF_Tau(:,:) = GF_Tau(:,:) / Beta_Temp

            !Analytically add back on contribution from tail
            GF_Tau(:,:) = GF_Tau(:,:) +     &
                GF%C2_Coeffs(:,:)*cmplx( ((TauPoint/2.0_dp) - Beta_Temp/4.0_dp),zero,dp) + &
                GF%C3_Coeffs(:,:)*cmplx( (Beta_Temp*TauPoint - TauPoint**2)/4.0_dp,zero,dp)
            !Diagonal C1 part
            do i = 1,nSites
                GF_Tau(i,i) = GF_Tau(i,i) - 0.5_dp
            enddo
        else

            allocate(eye(nSites,nSites))
            eye(:,:) = zzero
            do i = 1,nSites
                eye(i,i) = zone
            enddo

            !The simplest way first, removing the 1/iw tail analytically (from the
            !diagonals only)
            do i = 1,nMatsubara
                GF_Tau(:,:) = GF_Tau +  &
                    exp(-cmplx(zero,MatsuPoints(i)*TauPoint,dp)) * &
                    (GF%Matsu(:,:,i) - eye(:,:)*(one/cmplx(zero,MatsuPoints(i),dp)) )
            enddo

            GF_Tau(:,:) = GF_Tau(:,:) / Beta_Temp

            !Add back on the analytically FT'ed tail
            GF_Tau(:,:) = GF_Tau(:,:) - 0.5_dp*eye(:,:)

            deallocate(eye)
        endif

    end subroutine MatsuToImTimeFT

    !P = -2 x G(tau = Beta)
    !Beta will always be the last imaginary time point 
    subroutine GetGSDensityFromMatsuGF(GF,PMat)
        implicit none
        type(GreensFunc), intent(in) :: GF
        real(dp), intent(out) :: PMat(nSites,nSites)
        complex(dp), allocatable :: TempMat(:,:)
        integer :: i,j
        character(len=*), parameter :: t_r='GetGSDensityFromMatsuGF'

        allocate(TempMat(nSites,nSites))

        if(abs(ImTimePoints(nImTimePoints)-Beta_Temp).gt.1.0e-8_dp) then
            call stop_all(t_r,'Error in grids')
        endif

        call MatsuToImTimeFT(GF,ImTimePoints(nImTimePoints),TempMat)
        
        do i = 1,nSites
            do j = i,nSites
                if(abs(aimag(TempMat(i,j))).gt.1.0e-7_dp) then
                    !Though perhaps this is ok other than the diagonals?
                    call stop_all(t_r,'Density matrix complex')
                elseif(abs(TempMat(i,j)-conjg(TempMat(j,i))).gt.1.0e-7_dp) then
                    call stop_all(t_r,'Density matrix not hermitian')
                endif
            enddo
        enddo

        PMat(:,:) = -2.0_dp*real(TempMat(:,:),dp)

        deallocate(TempMat)
        if(tWriteOut) call writematrix(PMat,'Density Matrix',.true.) 

        !write(6,*) "First element: ",DensityMat_GV(1,1)

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
    !nFacTau is hardcoded initially to be 10.
    
    !The grids are returned in:
    !   MatsuPoints(1:nMatsubara)       [First 1:nMatsu/2 negative, nMatsu/2+1:nMatsu positive]
    !   ImTimePoints(1:nImTimePoints)   [1 is 0 time and nImTimePoints is Beta]

    !TODO:  Ultimately, we don't want a uniform grid in imaginary time
    !       We want to sample low and high times more
    subroutine SetupGrids()
        implicit none
        integer :: i,j
        real(dp) :: Delta_Tau
        character(len=*), parameter :: t_r='SetupGrids'

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
        !Number of frequency points should always be even
        !in order to get the same final point at both \pm freqs.
        if(mod(nMatsubara,2).eq.1) nMatsubara = nMatsubara + 1
        allocate(MatsuPoints(nMatsubara))
        !The first 1:nMatsu/2 are negative
        !The next nMatsu/2+1:nMatsu are positive
        j = 0
        do i = -nMatsubara/2,(nMatsubara/2)-1
            j = j + 1
            MatsuPoints(j) = (2*i+1)*pi / Beta_Temp
        enddo
        if(j.ne.nMatsubara) call stop_all(t_r,'Error in grids')

        write(6,"(A,I7)") "Number of matsubara frequency points set to: ",nMatsubara
        write(6,"(A,F20.10)") "Highest matsubara frequency point: ",MatsuPoints(nMatsubara)

        if(tFitTails) then
            write(6,"(A,F11.4)") "Fitting tails of greens functions beyond matsubara frequency: ",TailStart
            !Find fit points
            do i = 1,nMatsubara
                if(MatsuPoints(i).gt.-TailStart) then
                    iTailNeg = i-1
                    exit
                endif
            enddo
            do i = nMatsubara,1,-1
                if(MatsuPoints(i).lt.TailStart) then
                    iTailPos = i+1
                    exit
                endif
            enddo
            if(iTailNeg.ne.(nMatsubara-iTailPos+1)) call stop_all(t_r,'Error with tails')
            write(6,"(A,I7)") "Number of frequency points to be fit to tail of greens function: ",2*iTailNeg
            if(iTailNeg.lt.4) then
                call stop_all(t_r,'Number of frequency points to fit tail too small. Increase number of frequencies')
            endif

            !Calculate denominators for fitting of tails
            c2_denom = zero
            c3_denom = zero
            do i = 1,iTailNeg
                c2_denom = c2_denom + one/(MatsuPoints(i)**4)
                c3_denom = c3_denom + one/(MatsuPoints(i)**6)
            enddo
            c2_denom = c2_denom * 2.0_dp
            c3_denom = c3_denom * 2.0_dp

        endif

        !The number of Im Time points is controlled by ScaleImTime which
        !oversamples to avoid aliasing
        nImTimePoints = nint(ScaleImTime * real(nMatsubara,dp) * Beta_Temp / pi)
        allocate(ImTimePoints(nImTimePoints))
        write(6,"(A,I7)") "Number of imag-time points set to: ",nImTimePoints

        Delta_Tau = Beta_Temp/(nImTimePoints-1)
        write(6,"(A,G16.5)") "Imaginary time step in grid: ",Delta_Tau
        ImTimePoints(:) = zero
        do i = 2,nImTimePoints
            ImTimePoints(i) = (i-1)*Delta_Tau
        enddo

    end subroutine SetupGrids
    
    !For a given greens function, fit the tail and fill the global
    !data on the coefficients for the high-frequency expansion
    !Will work for SE and GF (providing the SE C0 and C1 coeffs have been
    !generated)
    subroutine FitGFTail(GF)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        real(dp), allocatable :: c2_num(:,:),c3_num(:,:),eye(:,:)
        integer :: a,i,j
        character(len=*), parameter :: t_r='FitGFTail'

        allocate(c2_num(nSites,nSites))
        allocate(c3_num(nSites,nSites))

        c2_num(:,:) = zero
        c3_num(:,:) = zero

        if(GF%tGF) then
            allocate(eye(nSites,nSites))    !For the C1 coeffs
            eye(:,:) = zero
            do a = 1,nSites
                eye(a,a) = one
            enddo
            do a = 1,iTailNeg
                c2_num(:,:) = c2_num(:,:) - real(GF%Matsu(:,:,a))/(MatsuPoints(a)**2)
                c3_num(:,:) = c3_num(:,:) + (aimag(GF%Matsu(:,:,a)) + &
                    eye(:,:)/MatsuPoints(a)) / (MatsuPoints(a)**3)
            enddo
            do a = iTailPos,nMatsubara
                c2_num(:,:) = c2_num(:,:) - real(GF%Matsu(:,:,a))/(MatsuPoints(a)**2)
                c3_num(:,:) = c3_num(:,:) + (aimag(GF%Matsu(:,:,a)) + &
                    eye(:,:)/MatsuPoints(a)) / (MatsuPoints(a)**3)
            enddo
            deallocate(eye)
        else
            !SE
            do a = 1,iTailNeg
                c2_num(:,:) = c2_num(:,:) - (real(GF%Matsu(:,:,a)) - GF%C0_Coeffs(:,:)) &
                    / (MatsuPoints(a)**2)
                c3_num(:,:) = c3_num(:,:) + (aimag(GF%Matsu(:,:,a)) + &
                    GF%C1_Coeffs(:,:)/MatsuPoints(a)) / (MatsuPoints(a)**3)
            enddo
            do a = iTailPos,nMatsubara
                c2_num(:,:) = c2_num(:,:) - (real(GF%Matsu(:,:,a)) - GF%C0_Coeffs(:,:)) &
                    / (MatsuPoints(a)**2)
                c3_num(:,:) = c3_num(:,:) + (aimag(GF%Matsu(:,:,a)) + &
                    GF%C1_Coeffs(:,:)/MatsuPoints(a)) / (MatsuPoints(a)**3)
            enddo
        endif

        GF%C2_Coeffs(:,:) = c2_num(:,:) / c2_denom
        GF%C3_Coeffs(:,:) = c3_num(:,:) / c3_denom

        do i = 1,nSites
            do j = i+1,nSites
                if(abs(GF%C2_Coeffs(j,i)-GF%C2_Coeffs(i,j)).gt.1.0e-7_dp) then
                    call stop_all(t_r,'C2 coefficient matrix is not symmetric')
                endif
                if(abs(GF%C3_Coeffs(j,i)-GF%C3_Coeffs(i,j)).gt.1.0e-7_dp) then
                    call stop_all(t_r,'C3 coefficient matrix is not antisymmetric')
                endif
            enddo
        enddo

        if(tWriteOut) then
            call writematrix(GF%C2_Coeffs,'C2 Coeffs for GF tails',.true.)
            call writematrix(GF%C3_Coeffs,'C3 Coeffs for GF tails',.true.)
        endif

        deallocate(c2_num,c3_num)

    end subroutine FitGFTail

    subroutine DeallocateGF(GF)
        implicit none
        type(GreensFunc), intent(inout) :: GF

        if(allocated(GF%Matsu)) deallocate(GF%Matsu)
        if(allocated(GF%Tau)) deallocate(GF%Tau)
        if(allocated(GF%C0_Coeffs)) deallocate(GF%C0_Coeffs)
        if(allocated(GF%C1_Coeffs)) deallocate(GF%C1_Coeffs)
        if(allocated(GF%C2_Coeffs)) deallocate(GF%C2_Coeffs)
        if(allocated(GF%C3_Coeffs)) deallocate(GF%C3_Coeffs)

    end subroutine DeallocateGF

    subroutine AllocateGF(GF,n,tGF)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        integer, intent(in) :: n
        logical, intent(in) :: tGF

        call deallocateGF(GF)
        allocate(GF%Matsu(n,n,nMatsubara))
        allocate(GF%Tau(n,n,nImTimePoints))
        GF%Matsu(:,:,:) = zzero
        GF%Tau(:,:,:) = zzero
        GF%tGF = tGF
        if(tFitTails) then
            allocate(GF%C2_Coeffs(n,n))
            allocate(GF%C3_Coeffs(n,n))
            GF%C2_Coeffs(:,:) = zero
            GF%C3_Coeffs(:,:) = zero
            if(.not.tGF) then
                !A self energy
                allocate(GF%C0_Coeffs(n,n))
                allocate(GF%C1_Coeffs(n,n))
                GF%C0_Coeffs(:,:) = zero
                GF%C1_Coeffs(:,:) = zero
            endif
        endif

    end subroutine AllocateGF

    !TODO: Everything should be built and stored in kspace
    subroutine AllocateMem_GV()
        implicit none
        integer :: nCompMem

        ! x 3 because we allocate a 'previous' SE on the matsubara grid in a bit
        nCompMem = (nSites**2)*(nMatsubara*3+nImTimePoints*2)
        write(6,"(A,F10.3,A)") "Total memory required for SE and GFs: ",nCompMem*ComptoGb," Gb"

        call AllocateGF(GLat_GV,nSites,.true.)
        call AllocateGF(SE_GV,nSites,.false.)

        if(allocated(DensityMat_GV)) deallocate(DensityMat_GV)
        if(allocated(DensityMat_MF_GV)) deallocate(DensityMat_MF_GV)
        if(allocated(FockMat_GV)) deallocate(FockMat_GV)

        allocate(DensityMat_GV(nSites,nSites))
        allocate(DensityMat_MF_GV(nSites,nSites))
        allocate(FockMat_GV(nSites,nSites))

    end subroutine AllocateMem_GV

    !TODO: Build this in k-space
    subroutine BuildMatsubaraGF(GF,LatChemPot,tCorr)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        real(dp), intent(in) :: LatChemPot
        logical, intent(in) :: tCorr
        complex(dp), allocatable :: SingleFreqMat(:,:)
        integer :: i,j

        GF%Matsu(:,:,:) = zzero
!$OMP PARALLEL DO PRIVATE(j,SingleFreqMat)
        do i = 1,nMatsubara
            allocate(SingleFreqMat(nSites,nSites))
            !Build matrix
            if(tCorr) then
                SingleFreqMat(:,:) = - FockMat_GV(:,:) - SE_GV%Matsu(:,:,i)
            else
                SingleFreqMat(:,:) = - h0(:,:)
            endif
            do j = 1,nSites
                SingleFreqMat(j,j) = SingleFreqMat(j,j) + cmplx(LatChemPot,MatsuPoints(i),dp)
            enddo
            !Invert
            call mat_inv(SingleFreqMat,GF%Matsu(:,:,i),nSites)
            deallocate(SingleFreqMat)
        enddo
!$OMP END PARALLEL DO

        if(tFitTails) call FitGFTail(GF)

    end subroutine BuildMatsubaraGF
    
    !This function returns the overpopulation of the system with electrons <N> - N_desired,
    !for a given chemical potential.
    function ExcessElec(GF,ChemPot,tCorr) result(res)
        implicit none
        type(GreensFunc), intent(inout) :: GF
        real(dp), intent(in) :: ChemPot
        logical, intent(in) :: tCorr
        real(dp) :: res
        real(dp) :: TotElec
        real(dp), allocatable :: GF_Density(:,:)
        integer :: i

        !Build the Matsubara greens function
        call BuildMatsubaraGF(GF,ChemPot,tCorr)

        !First, find the density matrix from the Matsubara greens function
        TotElec = zero
        if(tCorr) then
            call GetGSDensityFromMatsuGF(GF,DensityMat_GV)
            !Then find the number of electrons from the density
            do i = 1,nSites
                TotElec = TotElec + DensityMat_GV(i,i)
            enddo
        else
            allocate(GF_Density(nSites,nSites))
            call GetGSDensityFromMatsuGF(GF,GF_Density)

            !Then find the number of electrons from the density
            do i = 1,nSites
                TotElec = TotElec + GF_Density(i,i)
            enddo
            deallocate(GF_Density)
        endif

        res = TotElec - real(NEl,dp)

        if(tCorr) nElec_GV = TotElec

        !write(6,*) "ChemPot, NElec: ",ChemPot,nElec_GV

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
        DeltaPMF = DeltaPMF / real(nSites**2,dp)
        DeltaP = DeltaP / real(nSites**2,dp)
        DiffP = DiffP / real(nSites**2,dp)
        
        DeltaMu = LatChemPot - Previous_ChemPot

        Previous_ChemPot = LatChemPot
        PreviousP_MF(:,:) = DensityMat_MF_GV(:,:)
        PreviousP(:,:) = DensityMat_GV(:,:)

    end subroutine FindDensityChanges

    
    !For given fock matrix, find a chemical potential range which brackets the
    !desired number of electrons
    subroutine BracketChemPot(func,LatChemPot,ChemPotLow,ChemPotHigh,GF,tCorr)
        implicit none
        real(dp), intent(in) :: LatChemPot
        real(dp), intent(out) :: ChemPotLow,ChemPotHigh
        type(GreensFunc), intent(inout) :: GF
        logical, intent(in) :: tCorr
        interface
            function func(GF,x,tCorr)
                use const, only: dp 
                use GF2Data, only: GreensFunc
                implicit none
                type(GreensFunc), intent(inout) :: GF
                real(dp), intent(in) :: x
                logical, intent(in) :: tCorr
                real(dp) :: func
            end function func
        end interface
        real(dp) :: x1,x2,f1,f2
        integer :: i
        real(dp), parameter :: fac = 1.5_dp
        character(len=*), parameter :: t_r='BracketChemPot'

        !Set initial guess of 0.05t higher and lower than guess
        ChemPotLow = LatChemPot - 0.05
        ChemPotHigh = LatChemPot + 0.05

        x1 = ChemPotLow
        x2 = ChemPotHigh

        f1 = func(GF,x1,tCorr)
        f2 = func(GF,x2,tCorr)
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
                f1 = func(GF,x1,tCorr)
            else
                x2 = x2+fac*(x2-x1)
                f2 = ExcessElec(GF,x2,tCorr)
            endif
        enddo
        call stop_all(t_r,'Couldnt manage to bracket the chemical potential')

    end subroutine BracketChemPot

    subroutine SplineBeta(GF,SecDerivs,n,Tau,GFVal)
        implicit none
        type(GreensFunc), intent(in) :: GF
        integer, intent(in) :: n
        complex(dp), intent(in) :: SecDerivs(n,n,nImTimePoints)
        real(dp), intent(in) :: Tau
        complex(dp), intent(out) :: GFVal(n,n)
        !local
        integer :: klo,khi,k
        real(dp) :: a,b,c,d

        GFVal(:,:) = zzero
        !Bisection to find the right place
        klo = 1
        khi = n
        do while (khi-klo.gt.1)
            k = (khi+klo)/2
            if(ImTimePoints(k).gt.Tau) then
                khi = k
            else
                klo = k
            endif
        enddo
        !klo and khi now bracket the value of tau
        a = (ImTimePoints(khi)-tau)/(ImTimePoints(khi)-ImTimePoints(klo))
        b = (tau - ImTimePoints(klo))/(ImTimePoints(khi)-ImTimePoints(klo))
        c = (one/6.0_dp)*(a**3 - a)*((ImTimePoints(khi)-ImTimePoints(klo))**2)
        d = (one/6.0_dp)*(b**3 - b)*((ImTimePoints(khi)-ImTimePoints(klo))**2)

        GFVal(:,:) = a*GF%Tau(:,:,klo) + b*GF%Tau(:,:,khi) + &
            c*SecDerivs(:,:,klo) + d*SecDerivs(:,:,khi)

    end subroutine SplineBeta

    subroutine FindSplineCoeffs(GF,SecDerivs,n)
        implicit none
        type(GreensFunc), intent(in) :: GF
        integer, intent(in) :: n
        complex(dp), intent(out) :: SecDerivs(n,n,nImTimePoints)
        !local
        real(dp) :: DeltaTau
        complex(dp) :: TauDerivZero,TauDerivBeta,TauSecDerivZero,TauSecDerivBeta
        integer :: i,j,k,info
        integer, allocatable :: ipiv(:)
        complex(dp), allocatable :: b(:),CoeffMat(:,:)
        integer, parameter :: iBoundaryOrder = 2
        character(len=*), parameter :: t_r='FindSplineCoeffs'

        allocate(b(n))
        allocate(CoeffMat(n,n))
        allocate(ipiv(n))
        SecDerivs(:,:,:) = zero

        DeltaTau = ImTimePoints(2) - ImTimePoints(1)

        do i = 1,n
            do j = 1,n

                !First, find the difference in the second derivatives at tau=0 and beta
                !as these are the boundary conditions
                if(iBoundaryOrder.eq.2) then
                    !Second order: First four points (note this is the accuracy of the
                    !derivative in the splines...)
                    TauSecDerivZero = 2.0_dp*GF%Tau(j,i,1) - 5.0_dp*GF%Tau(j,i,2) + 4.0_dp*GF%Tau(j,i,3) - GF%Tau(j,i,4)
                    TauSecDerivBeta = 2.0_dp*GF%Tau(j,i,nImTimePoints) - 5.0_dp*GF%Tau(j,i,nImTimePoints-1) + &
                        4.0_dp*GF%Tau(j,i,nImTimePoints-2) - GF%Tau(j,i,nImTimePoints-3)
                    TauDerivZero = (-3.0_dp/2.0_dp)*GF%Tau(j,i,1) + 2.0_dp*GF%Tau(j,i,2) - (1.0_dp/2.0_dp)*GF%Tau(j,i,3)
                    TauDerivBeta = (3.0_dp/2.0_dp)*GF%Tau(j,i,nImTimePoints) - 2.0_dp*GF%Tau(j,i,nImTimePoints-1) + &
                        (1.0_dp/2.0_dp)*GF%Tau(j,i,nImTimePoints-2)
                elseif(iBoundaryOrder.eq.3) then
                    !Third order: First five points
                    TauSecDerivZero = (35.0_dp/12.0_dp)*GF%Tau(j,i,1) - (26.0_dp/3.0_dp)*GF%Tau(j,i,2) + &
                        (19.0_dp/2.0_dp)*GF%Tau(j,i,3) - (14.0_dp/3.0_dp)*GF%Tau(j,i,4) + (11.0_dp/12.0_dp)*GF%Tau(j,i,5)
                    TauSecDerivBeta = (35.0_dp/12.0_dp)*GF%Tau(j,i,nImTimePoints) - (26.0_dp/3.0_dp)*GF%Tau(j,i,nImTimePoints-1) + &
                        (19.0_dp/2.0_dp)*GF%Tau(j,i,nImTimePoints-2) - (14.0_dp/3.0_dp)*GF%Tau(j,i,nImTimePoints-3) + &
                        (11.0_dp/12.0_dp)*GF%Tau(j,i,nImTimePoints-4)
                    TauDerivZero = (-11.0_dp/6.0_dp)*GF%Tau(j,i,1) + 3.0_dp*GF%Tau(j,i,2) - &
                        (3.0_dp/2.0_dp)*GF%Tau(j,i,3) + (1.0_dp/3.0_dp)*GF%Tau(j,i,4)
                    TauDerivBeta = (11.0_dp/6.0_dp)*GF%Tau(j,i,nImTimePoints) - 3.0_dp*GF%Tau(j,i,nImTimePoints-1) + &
                        (3.0_dp/2.0_dp)*GF%Tau(j,i,nImTimePoints-2) - (1.0_dp/3.0_dp)*GF%Tau(j,i,nImTimePoints-3)
                elseif(iBoundaryOrder.eq.4) then
                    !Fourth order: First six points
                    TauSecDerivZero = (15.0_dp/4.0_dp)*GF%Tau(j,i,1) - (77.0_dp/6.0_dp)*GF%Tau(j,i,2) + &
                        (107.0_dp/6.0_dp)*GF%Tau(j,i,3) - 13.0_dp*GF%Tau(j,i,4) + (61.0_dp/12.0_dp)*GF%Tau(j,i,5) - &
                        (5.0_dp/6.0_dp)*GF%Tau(j,i,6)
                    TauSecDerivBeta = (15.0_dp/4.0_dp)*GF%Tau(j,i,nImTimePoints) - (77.0_dp/6.0_dp)*GF%Tau(j,i,nImTimePoints-1) + &
                        (107.0_dp/6.0_dp)*GF%Tau(j,i,nImTimePoints-2) - 13.0_dp*GF%Tau(j,i,nImTimePoints-3) + &
                        (61.0_dp/12.0_dp)*GF%Tau(j,i,nImTimePoints-4) - (5.0_dp/6.0_dp)*GF%Tau(j,i,nImTimePoints-5)
                    TauDerivZero = (-25.0_dp/12.0_dp)*GF%Tau(j,i,1) + 4.0_dp*GF%Tau(j,i,2) - &
                        3.0_dp*GF%Tau(j,i,3) + (4.0_dp/3.0_dp)*GF%Tau(j,i,4) - (1.0_dp/4.0_dp)*GF%Tau(j,i,5)
                    TauDerivBeta = (25.0_dp/12.0_dp)*GF%Tau(j,i,nImTimePoints) - 4.0_dp*GF%Tau(j,i,nImTimePoints-1) + &
                        3.0_dp*GF%Tau(j,i,nImTimePoints-2) - (4.0_dp/3.0_dp)*GF%Tau(j,i,nImTimePoints-3) + &
                        (1.0_dp/4.0_dp)*GF%Tau(j,i,nImTimePoints-4)
                elseif(iBoundaryOrder.eq.5) then
                    !Fifth order: First seven points
                    TauSecDerivZero = (203.0_dp/45.0_dp)*GF%Tau(j,i,1) - (87.0_dp/5.0_dp)*GF%Tau(j,i,2) + &
                        (117.0_dp/4.0_dp)*GF%Tau(j,i,3) - (254.0_dp/9.0_dp)*GF%Tau(j,i,4) + (33.0_dp/2.0_dp)*GF%Tau(j,i,5) - &
                        (27.0_dp/5.0_dp)*GF%Tau(j,i,6) + (137.0_dp/180.0_dp)*GF%Tau(j,i,7)
                    TauSecDerivBeta = (203.0_dp/45.0_dp)*GF%Tau(j,i,nImTimePoints) - &
                        (87.0_dp/5.0_dp)*GF%Tau(j,i,nImTimePoints-1) + &
                        (117.0_dp/4.0_dp)*GF%Tau(j,i,nImTimePoints-2) - (254.0_dp/9.0_dp)*GF%Tau(j,i,nImTimePoints-3) + &
                        (33.0_dp/2.0_dp)*GF%Tau(j,i,nImTimePoints-4) - &
                        (27.0_dp/5.0_dp)*GF%Tau(j,i,nImTimePoints-5) + (137.0_dp/180.0_dp)*GF%Tau(j,i,nImTimePoints-6)
                    TauDerivZero = (-137.0_dp/60.0_dp)*GF%Tau(j,i,1) + 5.0_dp*GF%Tau(j,i,2) - &
                        5.0_dp*GF%Tau(j,i,3) + (10.0_dp/3.0_dp)*GF%Tau(j,i,4) - &
                        (5.0_dp/4.0_dp)*GF%Tau(j,i,5) + (1.0_dp/5.0_dp)*GF%Tau(j,i,6)
                    TauDerivBeta = (137.0_dp/60.0_dp)*GF%Tau(j,i,nImTimePoints) - 5.0_dp*GF%Tau(j,i,nImTimePoints-1) + &
                        5.0_dp*GF%Tau(j,i,nImTimePoints-2) - (10.0_dp/3.0_dp)*GF%Tau(j,i,nImTimePoints-3) + &
                        (5.0_dp/4.0_dp)*GF%Tau(j,i,nImTimePoints-4) - (1.0_dp/5.0_dp)*GF%Tau(j,i,nImTimePoints-5)
                endif
                TauSecDerivZero = TauSecDerivZero / (DeltaTau**2)
                TauSecDerivBeta = TauSecDerivBeta / (DeltaTau**2)
                TauDerivZero = TauDerivZero / DeltaTau
                TauDerivBeta = TauDerivBeta / DeltaTau

                CoeffMat(:,:) = zzero
                do k = 2,n-1
                    CoeffMat(k,k) = cmplx(4.0_dp,zero,dp)
                    CoeffMat(k,k-1) = zone
                    CoeffMat(k,k+1) = zone
                enddo
                CoeffMat(1,1) = cmplx(6.0_dp/DeltaTau,zero,dp)
                CoeffMat(1,n) = cmplx(6.0_dp/DeltaTau,zero,dp)
                CoeffMat(n,1) = cmplx(-2.0_dp,zero,dp)
                CoeffMat(n,2) = cmplx(-1.0_dp,zero,dp)
                CoeffMat(n,n) = cmplx(2.0_dp,zero,dp)
                CoeffMat(n,n-1) = cmplx(1.0_dp,zero,dp)

                do k = 2,n-1
                    b(k) = GF%Tau(j,i,k+1) - 2.0_dp*GF%Tau(j,i,k) + GF%Tau(j,i,k-1)
                enddo
                b(1) = -(-TauSecDerivBeta - TauSecDerivZero)*DeltaTau
                b(n) = -(GF%Tau(j,i,2)-GF%Tau(j,i,1) + GF%Tau(j,i,nImTimePoints) - &
                    GF%Tau(j,i,nImTimePoints-1)) + (TauDerivBeta + TauDerivZero)*DeltaTau
                b(:) = b(:)*(6.0_dp/(DeltaTau**2))

                !Solve for spline second derivatives
                call ZGESV(n,1,CoeffMat,n,ipiv,b,n,info)
                if(info.ne.0) then
                    call stop_all(t_r,'Finding splines failed')
                endif
                SecDerivs(j,i,:) = b(:)

            enddo
        enddo
        deallocate(ipiv,CoeffMat,b)

    end subroutine FindSplineCoeffs
                


    !zbrent from numerical recipies
    function zbrent(func,x1,x2,tol,GF,tCorr)
        implicit none
        real(dp), intent(in) :: x1,x2,tol
        type(GreensFunc), intent(inout) :: GF
        logical, intent(in) :: tCorr
        real(dp) :: zbrent
        interface
            function func(GF,x,tCorr)
                use const, only: dp 
                use GF2Data, only: GreensFunc
                implicit none
                type(GreensFunc), intent(inout) :: GF
                real(dp), intent(in) :: x
                logical, intent(in) :: tCorr
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
        fa=func(GF,a,tCorr)
        fb=func(GF,b,tCorr)
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
            fb=func(GF,b,tCorr)
        end do
        call stop_all(t_r,'Exceeded maximum iterations')
        zbrent=b

    end function zbrent

end module GF2
