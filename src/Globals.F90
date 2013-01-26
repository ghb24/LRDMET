module Globals
    use const
    use timing, only: timer
    implicit none
    save

    integer :: LatticeDim   !Dimensionality
    integer :: N_occs   !Number of CS orbital occupations to loop over
    integer :: nSites   !The number of sites in the full system
    integer :: nSites_x   !The number of sites in the x direction
    integer :: nSites_y   !The number of sites in the y direction
    integer :: nImp     !The number of impurity sites
    integer :: nSys     !The number of bath sites
    real(dp) :: U       !Hubbard U
    logical :: tPeriodic !Use PBEs
    logical :: tAntiPeriodic !Use Anti-PBEs
    real(dp) :: ChemPot !The chemical potential of the system
    real(dp) :: HLGap   !The Homo-lumo gap of the system
    integer :: NEl      !The number of electrons in the entire system
    integer :: Elec    !The number of electrons in the embedded system
    integer :: nOcc     !The number of CS orbitals in the entire system
    integer :: EmbSize  !The total size of the embedded system
    real(dp) :: StartU,EndU,UStep   !The range and increment of U
    logical :: tRampDownOcc     !Whether to go up or down in filling fraction
    logical :: tConstructFullSchmidtBasis   !Whether to construct the full Schmidt basis or just the embedding basis
    real(dp) :: HL_Energy   !The energy of the embedded system from the solver
    real(dp) :: One_ElecE,Two_ElecE !The one and two-body contributions to the embedded system total energy from contractions with the HL 1RDM
    real(dp) :: TotalE_Imp,One_ElecE_Imp,Two_ElecE_Imp,CoupE_Imp !Energy contributions per impurity site for 1,2 electron and coupling to bath
    real(dp) :: Fillingerror,Actualfilling_Imp,Targetfilling_Imp    !We know what the filling should be, and this is what it actually is from the HL calc.
    integer :: nImpCombs,EmbCombs   !The size of triangular packed arrays (over impurity sites and embedding sites respectively)
    logical :: tHalfFill        !Half filling only
    logical :: tDiagFullSystem  !Diagonalize full system before DMET
    logical :: tMFResponse      !Calculate any mean-field response
    logical :: tNIResponse      !Calculate NI response
    logical :: tTDAResponse      !Calculate TDA response
    logical :: tRPAResponse      !Calculate RPA response
    logical :: tLR_DMET     !Attempt linear response based on partitioning of the perturbation into the schmidt basis of phi^0
    logical :: tEC_TDA_Response !Externally contracted response of DMET
    logical :: tIC_TDA_Response !Internall contracted response of DMET
    logical :: tCompleteDiag    !Complete rather than iterative diagonalization of the embedded system
!    real(dp) :: Omega=1.0_dp           !Perturbation frequency
    real(dp) :: Lambda=1.0_dp          !Strength of perturbation
    real(dp) :: Start_Omega,End_Omega,Omega_Step    !Parameters for Omega sweep
    integer :: pertsite=1          !Site of the density perturbation
    real(dp) :: ZerothBathNorm  !Normalization of the original bath orbital (required for correct normalization in the linear response)
    logical :: tDumpFCIDUMP
    logical :: tAnderson        !Whether to do anderson model, rather than hubbard model
    logical :: tChemPot         !Whether to include a chemical potential of U/2 at the impurity site of the anderson model
                                !Note that this potential only acts on the impurity site, and only acts on the interacting system.
                                !At half-filling, the system is naturally correct, so the chemical potential only wants to be added to the
                                !interacting case to stop the electrons fleeing the impurity site.
    logical :: tProjectOutNull  !For the LR - whether to attempt to remove linear dependencies in the basis before solving the equations
    logical :: tLR_ReoptGS      !For the LR - whether to reoptimize the ground state in the full space
    real(dp) :: MinS_Eigval     !For the LR - the smallest eigenvalue of S to keep
    logical :: tExplicitlyOrthog    !For the LR - explicitly orthogonalize the first-order solution
    logical :: tOrthogBasis     !For the LR - explicit calculate V and Q matrices, and do all calculations, in the orthogonal linear span of S
    integer :: iSolveLR         !For the LR - which routine to use to solve the LR equations.
                                ! 1   ZGESV   standard linear solver
                                ! 2   ZGELS   Advanced linear solver - should be better if hamiltonian nearly singular
                                ! 3   Direct inversion
                                ! 4   Complete diagonalization

    real(dp) :: HFEnergy    !Calculated HF energy
    real(dp) :: dDelta      !Broadening for spectral functions

    integer , allocatable :: allowed_occs(:)   !The list of CS occupations for the mean-field solution
    real(dp) , allocatable :: v_loc(:,:)    !The local correlation potential over the impurity sites
    real(dp) , allocatable :: h0(:,:)       !The mean-field core hamiltonian
    real(dp) , allocatable :: h0v(:,:)      !The mean-field core hamiltonian with local correlation potential striped across it
    real(dp) , allocatable :: HFOrbs(:,:)   !The eigenvectors of the mean-field solution
    real(dp) , allocatable :: HFEnergies(:)    !The eigenvalues of the mean-field hamiltonian
    real(dp) , allocatable :: FullHFOrbs(:,:)   !The true HF orbitals, including mean-field on site repulsion
    real(dp) , allocatable :: FullHFEnergies(:)   !The true fock eigenvalues, including MF onsite repulsion
    real(dp) , allocatable :: MeanFieldDM(:,:) !The 1e density matrix in the AO basis from the mean-field calculation
    real(dp) , allocatable :: EmbeddedBasis(:,:)    !The embedded basis orbitals of bath + impurity
    real(dp) , allocatable :: FullSchmidtBasis(:,:) !The full schmidt basis including core + virtual orbtials. Ordered: core, imp, bath, virt
    real(dp) , allocatable :: HFtoSchmidtTransform(:,:) !Transformation from HF basis to full Schmidt basis
    real(dp) , allocatable :: FockSchmidt(:,:)  !The fock operator in the full schmidt basis
    real(dp) , allocatable :: Emb_h0(:,:)        !Core hamiltonian in the embedded basis 
    real(dp) , allocatable :: Emb_h0v(:,:)        !Core hamiltonian in the embedded basis with correlation potential
    real(dp) , allocatable :: Emb_MF_DM(:,:)    !Mean-field density matrix in the embedded basis
    real(dp) , allocatable :: Emb_FockPot(:,:)  !The fock potential (i.e. stuff not in core hamiltnian, eg 2 electron terms) in emb basis.
                                                !Note that here we will set this to zero. We treat the fock matrix as just the core hamiltonian, and
                                                !ignore the 2e terms (U) which is diagonal in the basis. The effect of the U will just be captured in
                                                !the self-consistent correlation potential fitting, where U is used directly in the HL calculation over imp sites.
    real(dp) , allocatable :: Emb_CorrPot(:,:)  !The local potential which is added to the fock matrix to simulate the correlation effects in emb basis
    real(dp) , allocatable :: Emb_h1(:,:)       !The response of the bath orbital to the perturbation in the embedding basis
    real(dp) , allocatable :: Emb_Pert(:,:)     !The perturbation in the embedding basis
    real(dp) , allocatable :: HL_1RDM(:,:)      !The high-level calculation of the 1RDM over the embedded system
    real(dp) , allocatable :: HL_2RDM(:,:,:,:)  !The high-level calculation of the 2RDM over the embedded system
    real(dp) , allocatable :: Emb_Fock(:,:)     !The fock matrix in the embedded basis (h0 + v_loc (Emb_CorrPot) for hubbard)
    real(dp) , allocatable :: MFEmbOccs(:)      !The occupation numbers over the embedded system solved by the Emb_Fock
    real(dp) , allocatable :: vloc_change(:,:) !The change in the correlation potential over the impurity sites
    real(dp) , allocatable :: ResponseBasis(:,:)    !The impurity site + first order change in the bath wavefunction
    complex(dp) , allocatable :: SchmidtPert(:,:)
    real(dp) , allocatable :: FullHamil(:,:)    !In case we do a complete diagonalization
    real(dp) , allocatable :: Spectrum(:)       !Eigenvalues in case of a complete diagonalization

    !timers
    type(timer) :: Full_timer   !All routines 
    type(timer) :: FullSCF 
    type(timer) :: FCIDUMP
    type(timer) :: DiagT
    type(timer) :: ConstEmb
    type(timer) :: Trans1e
    type(timer) :: HL_Time
    type(timer) :: Fit_v_time
    !LR_SR
    type(timer) :: LR_SR_NonInt
    type(timer) :: LR_SR_TDA
    type(timer) :: LR_SR_RPA
    !LR_SR_EC
    type(timer) :: LR_EC_TDA_Precom !Precomputing (outside omega loop) various hamiltonians & generating det lists 
    type(timer) :: LR_EC_TDA_HBuild     !Building the hamiltonian at each omega 
    type(timer) :: LR_EC_TDA_SBuild     !Building the overlap at each omega
    type(timer) :: LR_EC_TDA_Project    !Diag S and project out null space
    type(timer) :: LR_EC_TDA_OptGS      !Diag H and 
    type(timer) :: LR_EC_TDA_BuildLR    !Construct LR equations
    type(timer) :: LR_EC_TDA_SolveLR    !Solve LR equations





end module Globals
