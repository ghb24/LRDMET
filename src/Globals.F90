module Globals
    use const
    implicit none

    integer :: LatticeDim   !Dimensionality
    integer :: N_occs   !Number of CS orbital occupations to loop over
    integer :: nSites   !The number of sites in the full system
    integer :: nImp     !The number of impurity sites
    integer :: nSys     !The number of bath sites
    real(dp) :: U       !Hubbard U
    logical :: tGSFCI   !Use Geralds FCI impurity solver
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
    logical :: tMFResponse      !Calculate mean-field response
    logical :: tLR_DMET     !Attempt linear response based on partitioning of the perturbation into the schmidt basis of phi^0
    logical :: tCompleteDiag    !Complete rather than iterative diagonalization of the embedded system
!    real(dp) :: Omega=1.0_dp           !Perturbation frequency
    real(dp) :: Lambda=1.0_dp          !Strength of perturbation
    real(dp) :: Start_Omega,End_Omega,Omega_Step    !Parameters for Omega sweep
    integer :: pertsite=1          !Site of the density perturbation
    real(dp) :: ZerothBathNorm  !Normalization of the original bath orbital (required for correct normalization in the linear response)
    logical :: tDumpFCIDUMP
    logical :: tAnderson        !Whether to do anderson model, rather than hubbard model
    logical :: tChemPot         !Whether to include a chemical potential of U/2 at the impurity site of the anderson model
    logical :: tProjectOutNull  !For the LR - whether to attempt to remove linear dependencies in the basis before solving the equations
    logical :: tLR_ReoptGS      !For the LR - whether to reoptimize the ground state in the full space
    real(dp) :: MinS_Eigval     !For the LR - the smallest eigenvalue of S to keep

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

end module Globals
