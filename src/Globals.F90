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
    integer :: nImp_x,nImp_y    !Number of impurities in each direction for the 2D system
    real(dp) :: U       !Hubbard U
    integer :: nU_Vals  !Number of explicitly specified U values to loop over
    logical :: tUHF     !UHF
    logical :: tPeriodic !Use PBEs
    logical :: tAntiPeriodic !Use Anti-PBEs
    integer :: iMaxIterDMET !Maximum iterations for DMET self-consistency
    logical :: tSCFHF   !Perform full scf hartree--fock calculation
    logical :: tWriteOut    !Write out additional info
    logical :: tCheck   !Perform extra checks in various subroutines to ensure correct working
    logical :: tReadSystem  !Read the system from files
    logical :: tThermal     !Are we calculating the thermal bath?
    logical :: tSingFiss    !For the Singlet Fission project
    real(dp) :: Temperature !k_B T In units of t
    real(dp) :: dTolDMET    !Convergence property of DMET calculation
    real(dp) :: ChemPot !The chemical potential of the system
    real(dp) :: HLGap   !The Homo-lumo gap of the system
    integer :: NEl      !The number of electrons in the entire system
    integer :: Elec    !The number of electrons in the embedded system
    integer :: nElecFill    !The input number of electrons in the entire system
    integer :: nOcc     !The number of CS orbitals in the entire system
    integer :: EmbSize  !The total size of the embedded system
    integer :: EmbSizeSpin  ! = EmbSize for RHF, and EmbSize*2 for UHF
    real(dp) :: StartU,EndU,UStep   !The range and increment of U
    logical :: tRampDownOcc     !Whether to go up or down in filling fraction
    logical :: tSaveCorrPot     !Use the previous correlation potential in the initialization of the next GS DMET calculation
    logical :: tReadInCorrPot   !Read in the correlation potential to use. No DMET SCF
    logical :: tFCIQMC          !Run FCIQMC as solver
    integer :: nNECICores        !Number of cores to run FCIQMC on
    logical :: tCoreH_EmbBasis  !Transform into the CoreH basis before writing FCIDUMP
    character(len=64) :: CorrPot_file   !File name with the correlation potentials in it 
    logical :: tConstructFullSchmidtBasis   !Whether to construct the full Schmidt basis or just the embedding basis
    logical :: tCompressedMats  !Whether in the non-direct davidson code to store compressed matrices or not
    real(dp) :: CompressThresh  !The integral threshold
    logical :: tDiag_kspace !Wheter to perform diagonalizations in k-space or not
    logical :: tContinueConvergence !Whether to continue convergence even if we have read in correlation potential
    real(dp) :: GS_Fit_Step !The step size for the finite difference gradient approximation in fitting the correlation potential
    real(dp) :: HL_Energy   !The energy of the embedded system from the solver
    real(dp) :: One_ElecE,Two_ElecE !The one and two-body contributions to the embedded system total energy from contractions with the HL 1RDM
    real(dp) :: TotalE_Imp,One_ElecE_Imp,Two_ElecE_Imp,CoupE_Imp !Energy contributions per impurity site for 1,2 electron and coupling to bath
    real(dp) :: Fillingerror,Actualfilling_Imp,Targetfilling_Imp    !We know what the filling should be, and this is what it actually is from the HL calc.
    real(dp) :: CoreEnergy      !The non-interacting core energy of the DMET wavefunction
    integer :: nImpCombs,EmbCombs   !The size of triangular packed arrays (over impurity sites and embedding sites respectively)
    logical :: tHalfFill        !Half filling only
    logical :: tDiagFullSystem  !Diagonalize full system before DMET
    logical :: tMFResponse      !Calculate any mean-field response
    logical :: tNIResponse      !Calculate NI response
    logical :: tTDAResponse      !Calculate TDA response
    logical :: tRPAResponse      !Calculate RPA response
    logical :: tLR_DMET     !Attempt linear response based on partitioning of the perturbation into the schmidt basis of phi^0
    logical :: tCorrNI_Spectra  !Calculate single-reference linear response functions based on correlated 1e hamiltonian
    logical :: tCorrNI_LocGF    !Calculate local GF based on correlated 1e hamil
    logical :: tCorrNI_LocDD    !Calculate local DD based on correlated 1e hamil
    logical :: tCorrNI_MomGF    !Calculate the k-dependent GF based on correlated 1e hamiltonian
    logical :: tEC_TDA_Response !Externally contracted response of DMET
    logical :: tIC_TDA_Response !Internall contracted response of DMET
    logical :: tCharged_MomResponse !Momentum resolved GFs
    logical :: tNoStatickBasis  !No basis of linear system coming from the schmidt basis of the static operator. Should be exact in NI limit any more.
    logical :: tCompleteDiag    !Complete rather than iterative diagonalization of the embedded system
    logical :: tNonDirDavidson  !Compute GS with a non-direct davidson algorithm
    logical :: tMinRes_NonDir   !Solve any systems of linear equations with a non-direct linear solver
    logical :: tGMRes_NonDir    !Solve any systems of linear equations with a non-direct GMRES linear solver
    logical :: tPrecond_MinRes  !Apply preconditioning to the solution of the linear equations
    integer :: nKrylov          !The number of krylov vectors to store for the GMRES algorithm
    logical :: tReuse_LS        ! Whether or not to reuse the previous frequency calculation for the solution of the linear equation
    real(dp) :: rtol_LR         !Tolerance for exit criterion for linear solver
    real(dp) :: Lambda=1.0_dp          !Strength of perturbation
    real(dp) :: Start_Omega,End_Omega,Omega_Step    !Parameters for Omega sweep
    logical :: tDumpFCIDUMP
    logical :: tAnderson        !Whether to do anderson model, rather than hubbard model
    logical :: tChemPot         !Whether to include a chemical potential of U/2 at the impurity site of the anderson model
                                !Note that this potential only acts on the impurity site, and only acts on the interacting system.
                                !At half-filling, the system is naturally correct, so the chemical potential only wants to be added to the
                                !interacting case to stop the electrons fleeing the impurity site.

    real(dp) :: HFEnergy    !Calculated HF energy
    logical :: tReadMats, tWriteMats    !Options for reading/writing the compressed matrices (N-electron only)
    integer :: iHamSize_N,iHamSize_Nm1,iHamSize_Np1     !Size of the compressed hamiltonians for N,N-1 and N+1 respectively

    integer , allocatable :: TD_Imp_Lat(:,:),TD_Imp_Phase(:,:)  !Parameterization of the orbital space for 2D hubbard
    integer , allocatable :: ImpSites(:)    !The list of site indices for the impurity sites
    integer , allocatable :: Perm_dir(:),Perm_indir(:)  !Mappings between the two indexing orders
    real(dp), allocatable :: U_Vals(:)      !The list of U_Values to loop over
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
    complex(dp) , allocatable :: SchmidtPert(:,:)
    real(dp) , allocatable, target :: HL_Vec(:)         !The ground state eigenvector
    real(dp) , allocatable, target :: FullHamil(:,:)    !In case we do a complete diagonalization
    real(dp) , allocatable :: Spectrum(:)       !Eigenvalues in case of a complete diagonalization
    real(dp) , allocatable :: J_Ints(:,:)       !Coulomb integrals over the impurity space
    real(dp) , allocatable :: X_Ints(:,:)       !Exchange integrals over the impurity space

    !Analogous Beta space arrays
    real(dp), allocatable :: v_loc_b(:,:) 
    real(dp), allocatable :: h0_b(:,:)
    real(dp), allocatable :: h0v_b(:,:)
    real(dp), allocatable :: HFEnergies_b(:)
    real(dp), allocatable :: HFOrbs_b(:,:)
    real(dp), allocatable :: MeanfieldDM_b(:,:)
    real(dp), allocatable :: FullSchmidtBasis_b(:,:) 
    real(dp), allocatable :: HFtoSchmidtTransform_b(:,:)
    real(dp), allocatable :: FockSchmidt_b(:,:)
    real(dp), allocatable :: EmbeddedBasis_b(:,:)
    real(dp), allocatable :: Emb_h0_b(:,:)
    real(dp), allocatable :: Emb_MF_DM_b(:,:)
    real(dp), allocatable :: Emb_FockPot_b(:,:)
    real(dp), allocatable :: Emb_CorrPot_b(:,:)
    real(dp), allocatable :: Emb_Fock_b(:,:)
    real(dp), allocatable :: Emb_h0v_b(:,:)
    real(dp), allocatable :: HL_1RDM_b(:,:)
    logical :: tBetaExcit   !For beta spin-orbital greens function excitations
    real(dp) :: ChemPot_b,HLGap_b

    !KSpace parameters
    integer :: nKPnts                           !Number of kpoints sampled in the 1st BZ
    logical :: tShift_Mesh                      ! = T <- Gamma centered mesh. = F <- Monkhorst-Pack mesh
    real(dp), allocatable :: RecipLattVecs(:,:) !Define the reciprocal lattice vectors (nDim,nDim)
    real(dp), allocatable :: KPnts(:,:)         !Define the kpoint vectors in the reciprocal lattice (nDim,nKPnts)
    complex(dp), allocatable :: RtoK_Rot(:,:)   !Rotation matrix from site basis to k-space
                                                !The order of the k-index is the same as the KPnts array
    complex(dp), allocatable :: HFtoKOrbs(:,:)  !The transformation matrix from HF orbitals to k-space
    complex(dp), allocatable :: HFtoKOrbs_b(:,:)    !For beta space
    logical :: tProjectHFKPnts                  !Whether to construct the rotation matrix from final HF orbitals to original plane wave k-points
    logical :: tKSpaceOrbs                      !Calculate the kspace orbitals from the final one-electron matrix
    integer, allocatable :: KVec_EMapping(:)    !This map gives the index of the orbitals in terms of energy
    integer, allocatable :: KVec_InvEMap(:)     !This takes the orbital, and returns its energetic index
    complex(dp) , allocatable :: k_vecs(:,:)    !The orbitals, ordered by kpoint ( nImp , nSites) 
    real(dp), allocatable :: k_HFEnergies(:)    !The HF Energies, ordered by kpoint
    complex(dp), allocatable :: k_HFtoSchmidtTransform(:,:) !Transform from the complex k-space eigenvectors (2nd ind) to the schmidt basis (Note: This ordering is the opposite way around to previously)
    integer :: nKCalcs                          !The number of k-space GFs to calculate
    integer :: KIndex                           !The index of the calculational kpoint we are on
    integer, allocatable :: KCalcs(:)           !The index of the kpoints to calculate the GF over

    !Linear response options
    real(dp) :: dDelta      !Broadening for spectral functions
    logical :: tDDResponse          !Calculate neutral DD response
    logical :: tChargedResponse     !The different perturbations to calculate the response for
    logical :: tProjectOutNull  !For the LR - whether to attempt to remove linear dependencies in the basis before solving the equations
    logical :: tLR_ReoptGS      !For the LR - whether to reoptimize the ground state in the full space
    real(dp) :: MinS_Eigval     !For the LR - the smallest eigenvalue of S to keep
    logical :: tExplicitlyOrthog    !For the LR - explicitly orthogonalize the first-order solution
    logical :: tOrthogBasis     !For the LR - explicit calculate V and Q matrices, and do all calculations, in the orthogonal linear span of S
    integer :: iMinRes_MaxIter  !For the LR solver - maximum iterations
    logical :: tRemoveGSFromH   !For the LR - whether to explicitly remove the GS from the hamiltonian before forming and solving the LR equations.
                                !Warning - this can remove the hermiticity of the hamiltonian
    integer :: iSolveLR         !For the LR - which routine to use to solve the LR equations.
                                ! 1   ZGESV   standard linear solver
                                ! 2   ZGELS   Advanced linear solver - should be better if hamiltonian nearly singular
                                ! 3   Direct inversion
                                ! 4   Complete diagonalization
    logical :: tAllImp_LR       ! Whether to calculate all nImp*nImp greens functions
    logical :: tSC_LR           ! Whether to self-consistently optimize the self-energy part of the LR h.
    integer :: iReuse_SE        ! =0: no memory of SE from previous frequencies, =1: Start macroiterations with previous converged SE, =2: Start NR microiterations with previous converged SE
    logical :: tNoHL_SE         ! Do not include the self-energy contribution in the construction of the hamiltonian used for the HL calculation (i.e. SE only enters through the NI contraction coefficients)
    integer :: iGF_Fit          ! Type of fitting in SC_LR: 0 - normal, 1 - DampedNR, 2 - Linesearch, 3 - Damped & linesearch
    logical :: tPartialSE_Fit   ! Do iPartialSE_Fit fits of the self-energy. Do not fit until self consistency between HL and NI GFs
    integer :: iPartialSE_Fit   ! Max number of fits to do of the self energy
    integer :: nVarSE           ! Number of independent variables in the packed self-energy matrix
    integer :: iSE_Constraints  ! Input constraints on flexibility of the self-energy matrix
    real(dp) :: DampingExponent ! Damping of self energy update
    logical :: tConvergeMicroSE ! Whether to converge the self-energy completely for each high-level calculation
    integer :: TDLat_Ni,TDLat_Nj
    
    !DMET_LR global data
    !When a non-hermitian self-energy is added, we need to seperately calculate the non-interacting wavefunctions expressed in the right-eigenvector space
    !and the left eigenvector space. Kets are in the right eigenvector space, and Bras in the left space.
    complex(dp), allocatable :: SchmidtPertGF_Cre_Ket(:,:) !The contraction coefficients (potentially for each impurity site) for the core excitations
    complex(dp), allocatable :: SchmidtPertGF_Ann_Ket(:,:) !The contraction coefficients (potentially for each impurity site) for the core excitations
    complex(dp), allocatable :: SchmidtPertGF_Cre_Bra(:,:) !The contraction coefficients (potentially for each impurity site) for the core excitations
    complex(dp), allocatable :: SchmidtPertGF_Ann_Bra(:,:) !The contraction coefficients (potentially for each impurity site) for the core excitations
    complex(dp), allocatable :: NI_LRMat_Cre(:,:)   !NI particle greens functions for each value of omega
    complex(dp), allocatable :: NI_LRMat_Ann(:,:)   !NI hole-addition greens functions for each value of omega
    complex(dp), allocatable :: SelfEnergy_Imp(:,:)    !The updated self-energy matrix over impurity sites
    complex(dp), allocatable :: Emb_h0v_SE(:,:)        !Neither this, or the selfEnergy itself, are hermitian
    complex(dp), allocatable :: h0v_SE(:,:)         !The full AO basis of one-electron + correlation potential + self-energy
    complex(dp), allocatable :: FockSchmidt_SE(:,:) !The one-electron hamiltonian over the whole space (apart from imp-imp block)
    complex(dp), allocatable :: FockSchmidt_SE_VV(:,:) !The one-electron hamiltonian over the virtual-virtual block
    complex(dp), allocatable :: FockSchmidt_SE_CC(:,:) !The one-electron hamiltonian over the core core block
    complex(dp), allocatable :: FockSchmidt_SE_VX(:,:) !The one-electron hamiltonian over the virtual:active block
    complex(dp), allocatable :: FockSchmidt_SE_CX(:,:) !The one-electron hamiltonian over the core:active block
    complex(dp), allocatable :: FockSchmidt_SE_XV(:,:) !The one-electron hamiltonian over the active:virtual block
    complex(dp), allocatable :: FockSchmidt_SE_XC(:,:) !The one-electron hamiltonian over the active:core block

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
    !LR_MR_EC
    type(timer) :: LR_EC_TDA_Precom !Precomputing (outside omega loop) various hamiltonians & generating det lists 
    type(timer) :: LR_EC_TDA_HBuild     !Building the hamiltonian at each omega 
    type(timer) :: LR_EC_TDA_SBuild     !Building the overlap at each omega
    type(timer) :: LR_EC_TDA_Project    !Diag S and project out null space
    type(timer) :: LR_EC_TDA_OptGS      !Diag H and 
    type(timer) :: LR_EC_TDA_BuildLR    !Construct LR equations
    type(timer) :: LR_EC_TDA_SolveLR    !Solve LR equations
    !LR_MR_EC_GF analogues
    type(timer) :: LR_EC_GF_Precom
    type(timer) :: LR_EC_GF_HBuild
    type(timer) :: LR_EC_GF_OptGS
    type(timer) :: LR_EC_GF_SolveLR
    type(timer) :: LR_EC_GF_FitGF

end module Globals
