module readinput
    use input
    use const
    use globals
    use LatticesData, only: CellShape
    use errors, only: stop_all,warning
    use GF2Data, only: nMatsubara, Beta_Temp, ScaleImTime, TailStart,tFitTails
    use GF2Data, only: MatsuEnergySumFac, ScaleImTimeSpline, tSpline, GF2_MaxIter
    use SC_Data, only: tReadChemPot,tStretchNILatticeHam,dStretchNILatticeHam
    use SC_Data, only: PotentialUpdateDamping,dSelfConsConv
    implicit none

    contains

    !These want to be input parameters at some point
    subroutine set_defaults()
        implicit none

        !Main options
        tDMETCalc = .false.
        tAnderson = .false. !Hubbard by default
        tChemPot = .false.
        tReadSystem = .false.
        tHalfFill = .true. 
        tUHF = .false.
        tThermal = .false.
        tSingFiss = .false.
        nSites = 24  
        LatticeDim = 1
        nImp = 1
        nElecFill = 0
        dTolDMET = 1.0e-8_dp
        StartU = 0.0_dp
        EndU = 4.1_dp
        UStep = 1.0_dp
        CellShape = 0
        tPeriodic = .false.
        tAntiPeriodic = .false. 
        iMaxIterDMET = 150
        tDiag_kspace = .false.
        GS_Fit_Step = 1.0e-5_dp
        tRampDownOcc = .true.
        tCompleteDiag = .true. 
        tSaveCorrPot = .false.
        tDumpFCIDUMP = .false.
        tDiagFullSystem = .false.
        tSCFHF = .false.
        tWriteout = .false.
        tReadInCorrPot = .false.
        CorrPot_file = 'CORRPOTS'
        tContinueConvergence = .true.
        tNonDirDavidson = .false.
        tFCIQMC = .false.
        tMP2 = .false.
        nNECICores = 0
        tCoreH_EmbBasis = .false.
        tCheck = .false.
        tCompressedMats = .false.
        CompressThresh = 1.0e-10_dp !Integral threshold for compressed matrices
        tReadMats = .false.
        tWriteMats = .false.
        iHamSize_N = 0
        iHamSize_Np1 = 0
        iHamSize_Nm1 = 0
        tCorrelatedBath = .false.
        tExactCorrPot = .false.
        tDMFTCalc = .false.
        tT1SCF = .false.

        !General LR options
        Start_Omega = 0.0_dp
        End_Omega = 4.05_dp
        Omega_Step = 0.1_dp
        dDelta = 0.01_dp

        !GF2
        tGF2 = .false.
        GF2_MaxIter = 150
        Beta_Temp = 100.0_dp
        nMatsubara = 250
        ScaleImTime = 5.0_dp
        tFitTails = .false. !Whether to fit the w^2 and w^3 terms
        MatsuEnergySumFac = zero
        ScaleImTimeSpline = zero
        tSpline = .false.

        !SR Response
        tMFResponse = .false. 
        tNIResponse = .false.
        tTDAResponse = .false.
        tRPAResponse = .false.
        tCorrNI_Spectra = .false.
        tCorrNI_LocGF = .false.
        tCorrNI_LocDD = .false.
        tCorrNI_MomGF = .false.
        tProjectHFKPnts = .false.
        tKSpaceOrbs = .false.

        !MR Response
        tDDResponse = .false.
        tChargedResponse = .false.
        tEC_TDA_Response = .false.
        tIC_TDA_Response = .false.
        tLR_DMET = .false. 
        tCharged_MomResponse = .false.
        tNoStatickBasis = .false.
        tConstructFullSchmidtBasis = .true. 
        tProjectOutNull = .false.
        tLR_ReoptGS = .false. 
        MinS_Eigval = 1.0e-9_dp
        tExplicitlyOrthog = .false.
        iSolveLR = 1
        tOrthogBasis = .false.
        tRemoveGSFromH = .false.
        tMinRes_NonDir = .false.
        tPreCond_MinRes = .false.
        nKrylov = 100
        rtol_LR = 1.0e-8_dp
        tReuse_LS = .false.
        tSC_LR = .false.
        tNoHL_SE = .true. 
        iReuse_SE = 0
        tAllImp_LR = .false.
        iGF_Fit = 0
        tPartialSE_Fit = .false.
        iPartialSE_Fit = 0
        iSE_Constraints = 1
        Damping_SE = one 
        tConvergeMicroSE = .false.
        iMinRes_MaxIter = 20000
        tBetaExcit = .false.
        nKCalcs = 30
        max_SE_iter = 0
        NIGF_WeightFac = 0.0_dp
        tRead_SelfEnergy = .false.
        tRandom_Init_SE = .false.
        tSE_Scan = .false.
        tRealSpaceSC = .true.   !By default, attempt self-consistency on the matsubara axis
        iLatticeFitType = 1     !By default, fit the greens functions (rather than inverses)
        iFitGFWeighting = 0     !By default, when doing lattice fits, fit to a flat model, rather than bias low weighted excitations
!        iLatticeCoups = 0
        iNonLocBlocks = 0
        iMaxFitMicroIter = 0
        iFitAlgo = 1            !Which fitting algorithm to use. 1 = simplex, 2 = Powell
        tReadCouplings = .false.    !Whether to read in previous lattice couplings fits
        tSkip_Lattice_Fit = .false. !Whether to skip the fitting of the lattice
        tEnvLatHam = .false.        !Whether to use the fit hamiltonian as the one-electron external terms
        tEveryOtherCoup = .false.   !Every other lattice coupling constrained to be zero
        tStaticBathFitLat = .false. !By default no fit hamiltonian in the bath space
        tOptGF_EVals = .false.      !Optimize lattice couplings, or lattice eigenvalues
        tAnalyticDerivs = .false.   !Whether to use analytic derivatives when optimizing with Lev-Mar algorithm
        dShiftLatticeEvals = zero   !Shift read in lattice eigenvalues?
        tFitRealFreq = .false.      !Fit on matsubara axis by default
        tFitPoints_Legendre = .false.   !Fit to Legendre points
        nFreqPoints = 0             !Number of points to fit
        Start_Omega_Im = zero                
        End_Omega_Im = zero
        Omega_Step_Im = one
        iMaxIter_MacroFit = 100
        tConstrainKSym = .false.    !Restrict the fit to conserve momentum
        tConstrainphSym = .false.   !Restrict the fit to conserve ph symmetry
        tImposeKSym = .false.       !Impose momentum symmetry on the fit values
        tImposephSym = .false.      !Impose ph symmetry on the fit values
        tDiagonalSC = .false.       !Diagonal GF approx for residual
        dFitTol_SC = 1.0e-5_dp      !Tolerance on the gradient of the fits
        tRemoveImpCoupsPreSchmidt = .false. !Do we want to remove impurity couplings pre schmidt decomp.?
        tRemakeStaticBath = .false.
        tFullReoptGS = .false.      !Whether to reoptimize the ground state in the static + dynamic bath space
        tSC_StartwGSCorrPot = .true.    !Whether to start the selfconsistency from h0v or h0
        iFitStyle = 1                   !1 = Direct fitting. 2 = DMFT
        tCalcRealSpectrum = .true. 
        tAugMinRes = .false.            !For solving an augmented problem with minres without reducing the condition number
        tReadChemPot = .false.          !Read in the chemical potential from previous calculation?
        tStretchNILatticeHam = .false.  !Stretch the non-interacting lattice bandwidth?
        dStretchNILatticeHam = one
        PotentialUpdateDamping = one
        dSelfConsConv = 1.0e-5_dp

    end subroutine set_defaults

    subroutine read_input()
        use utils, only: get_free_unit
        implicit none
        integer :: command_argument_count,ios
        character(len=255) :: cInp
        character(len=100) :: w
        logical :: tEOF,tExists
        character(len=*), parameter :: t_r='read_input'

        if(command_argument_count().le.0) then
            write(6,"(A)") "No input file found. Running from defaults."
            return
        endif
        call get_command_argument(1,cInp)
        write(6,*) "Reading from file: ",trim(cInp)
        inquire(file=cInp,exist=tExists)
        if(.not.texists) call stop_all(t_r,'File '//trim(cInp)//' does not exist.')
        ir = get_free_unit()
        open(ir,file=cInp,status='old',form='formatted',err=99,iostat=ios)
        call input_options(echo_lines=.true.,skip_blank_lines=.true.)
        write(6,'(/,64("*"),/)')

        
        do while(.true.)
            call read_line(tEOF)
            if(tEOF) exit
            call readu(w)
            !First, search for block
            select case(w)
            case("MODEL")
                call ModelReadInput()
            case("GF2")
                call GF2ReadInput()
            case("LINEAR_RESPONSE")
                call LRReadInput()
            case("SELF_CONSISTENCY")
                call SCReadInput()
            case("END")
                exit
            case default
                call stop_all(t_r,'Input block: '//trim(w)//' not recognized')
            end select
        enddo
        write(6,'(/,64("*"),/)')
        close(ir)
99      if(ios.ne.0) then
            call stop_all(t_r,'Problem reading input file'//trim(cInp))
        endif

    end subroutine read_input

    subroutine GF2ReadInput()
        implicit none
        logical :: teof
        character(len=100) :: w
        character(len=*), parameter :: t_r='GF2ReadInput'

        tGF2 = .true.
        GF: do
            call read_line(teof)
            if(teof) exit
            call readu(w)
            select case(w)
            case("MAX_MACROITER")
                call readi(GF2_MaxIter)
            case("MATSUBARA_POINTS")
                call readi(nMatsubara)
            case("TAU_POINTS_FACTOR")
                call readf(ScaleImTime)
            case("MATSU_TAILS")
                tFitTails = .true.
                call readf(TailStart)
            case("ENERGYSUM_TAILS")
                call readf(MatsuEnergySumFac)
            case("SPLINE")
                tSpline = .true.
                call readf(ScaleImTimeSpline)
            case("END")
                exit
            case default
                write(6,"(A)") "ALLOWED KEYWORDS IN GF2 BLOCK: "
                write(6,"(A)") "MAX_MACROITER"
                write(6,"(A)") "MATSUBARA_POINTS"
                write(6,"(A)") "TAU_POINTS_FACTOR"
                write(6,"(A)") "MATSU_TAILS"
                write(6,"(A)") "ENERGYSUM_TAILS"
                write(6,"(A)") "SPLINE"
                call stop_all(t_r,'Keyword '//trim(w)//' not recognized')
            end select
        enddo GF   

    end subroutine GF2ReadInput

    subroutine ModelReadInput()
        implicit none
        logical :: teof
        character(len=100) :: w,w2,w3
        logical :: tSpecifiedHub,tSpecifiedAnd,tMultipleOccs
        real(dp) :: UVals(100)
        integer :: i
        character(len=*), parameter :: t_r='ModelReadInput'

        tSpecifiedHub = .false.
        tSpecifiedAnd = .false.
        tMultipleOccs = .false.

        i = 0
        UVals(:) = -1.0_dp

        Model: do
            call read_line(teof)
            if(teof) exit
            call readu(w)
            select case(w)
            case("SYSTEM")
                if(item.lt.nitems) then
                    call readu(w2)
                    select case(w2)
                    case("HUBBARD")
                        tAnderson = .false.
                        tSpecifiedHub = .true.
                    case("ANDERSON")
                        tAnderson = .true.
                        tSpecifiedAnd = .true.
                        tChemPot = .true.   !By default include chemical potential in the model
                        if(item.lt.nitems) then
                            call readu(w3)
                            select case(w3)
                            case("NO_CHEMPOT")
                                tChemPot = .false.
                            case default
                                call stop_all(t_r,'Keyword after ANDERSON no recognised. Should be "NO_CHEMPOT"')
                            end select
                        endif
                    case("TWOBAND_LATTICE")
                        tSingFiss = .true.
                    case("READ")
                        tReadSystem = .true.    !Read the system from files.
                    case default
                        call stop_all(t_r,'Keyword after "SYSTEM" not recognised')
                    end select
                else
                    call stop_all(t_r,'No model specified after "SYSTEM"')
                endif
            case("SITES")
                call readi(nSites)
                if(item.lt.nitems) then
                    call readi(LatticeDim)
                    if(LatticeDim.gt.2) then
                        call stop_all(t_r,'Can only deal with lattice dimensions =< 2')
                    endif
                endif
            case("CELLSHAPE")
                !1 = Gerald, 2 = Square, 3 = Tilted
                call readi(CellShape)
            case("UHF")
                tUHF = .true.
            case("U_VALS")
                do while(item.lt.nitems) 
                    i = i+1
                    call readf(UVals(i))
                enddo
            case("U")
                call readf(StartU)
                call readf(EndU)
                call readf(UStep)
            case("REUSE_CORRPOT")
                tSaveCorrPot = .true.
            case("READ_CORRPOT")
                tReadInCorrPot = .true.
                tContinueConvergence = .false.
                if(item.lt.nitems) then
                    CorrPot_file = ''
                    call readu(CorrPot_file)
                endif
                if(item.lt.nitems) then
                    call readu(w2)
                    select case(w2)
                    case("CONT_CONV")
                        tContinueConvergence = .true.
                    case default
                        call stop_all(t_r,'Keyword after "READ_CORRPOT" not recognised')
                    end select
                endif
            case("FITTING_STEPSIZE")
                call readf(GS_Fit_Step)
            case("T1_SCF")
                tT1SCF = .true.
            case("PBC")
                tPeriodic = .true.
            case("APBC")
                tAntiPeriodic = .true.
            case("TEMPERATURE")
                !call readf(Temperature)
                !tThermal = .true.
                call readf(Beta_Temp)
            case("MAXITER_DMET")
                call readi(iMaxIterDMET)
                if(item.lt.nitems) then
                    call readf(dTolDMET)
                endif
            case("SCF_HF")
                tSCFHF = .true.
            case("MP2")
                tMP2 = .true.
            case("KSPACE_DIAG")
                tDiag_kspace = .true.
            case("IMPSITES")
                tDMETCalc = .true.
                call readi(nImp)
            case("COMPRESSMATS")
                tCompressedMats = .true.
                if(item.lt.nitems) then
                    !Optional threshold argument
                    call readf(CompressThresh)
                endif
            case("READCMPRSMATS")
                tReadMats = .true.
            case("WRITECMPRSMATS")
                tWriteMats = .true.
            case("N_CMPRSHAMSIZE")
                call readi(iHamSize_N)
            case("NM1_CMPRSHAMSIZE")
                call readi(iHamSize_Nm1)
            case("NP1_CMPRSHAMSIZE")
                call readi(iHamSize_Np1)
            case("COMPLETE_DIAG")
                tCompleteDiag = .true.
            case("DAVIDSON")
                tCompleteDiag = .false.
            case("FCIQMC")
                tFCIQMC = .true.
                tCompleteDiag = .false.
                if(item.lt.nitems) then
                    call readi(nNECICores)
                endif
            case("COREH_EMBBASIS")
                tCoreH_EmbBasis = .true.
            case("NONDIR_DAVIDSON") 
                tNonDirDavidson = .true.
                tCompleteDiag = .false.
            case("DIAG_SYSTEM")
                tDiagFullSystem = .true.
            case("HALF_FILL")
                tHalfFill = .true.
            case("EXACT_CORRPOT")
                tExactCorrPot = .true.
            case("FILLING")
                tMultipleOccs = .true.
                call readi(nElecFill)
                tHalfFill=.false.
            case("REDUCE_OCC")
                tMultipleOccs = .true.
                tRampDownOcc = .true.
                tHalfFill = .false.
            case("INCREASE_OCC")
                tMultipleOccs = .true.
                tRampDownOcc = .false.
                tHalfFill = .false.
            case("FCIDUMP")
                tDumpFCIDUMP = .true.
            case("DEBUGOUTPUT")
                tWriteOut = .true.
            case("CHECK_SANITY")
                tCheck = .true.
            case("END")
                exit
            case default
                write(6,"(A)") "ALLOWED KEYWORDS IN MODEL BLOCK: "
                write(6,"(A)") "SYSTEM"
                write(6,"(A)") "SITES"
                write(6,"(A)") "UHF"
                write(6,"(A)") "U"
                write(6,"(A)") "U_VALS"
                write(6,"(A)") "REUSE_CORRPOT"
                write(6,"(A)") "EXACT_CORRPOT"
                write(6,"(A)") "FITTING_STEPSIZE"
                write(6,"(A)") "SCF_HF"
                write(6,"(A)") "TEMPERATURE"
                write(6,"(A)") "PBC"
                write(6,"(A)") "APBC"
                write(6,"(A)") "IMPSITES"
                write(6,"(A)") "MAXITER"
                write(6,"(A)") "HALF_FILL"
                write(6,"(A)") "FILLING"
                write(6,"(A)") "KSPACE_DIAG"
                write(6,"(A)") "COMPRESSMATS"
                write(6,"(A)") "READCMPRSMATS"
                write(6,"(A)") "WRITECMPRSMATS"
                write(6,"(A)") "N_CMPRSHAMSIZE"
                write(6,"(A)") "NM1_CMPRSHAMSIZE"
                write(6,"(A)") "NP1_CMPRSHAMSIZE"
                write(6,"(A)") "COMPLETE_DIAG"
                write(6,"(A)") "DAVIDSON"
                write(6,"(A)") "FCIQMC"
                write(6,"(A)") "COREH_EMBBASIS"
                write(6,"(A)") "DIAG_SYSTEM"
                write(6,"(A)") "REDUCE_OCC"
                write(6,"(A)") "INCREASE_OCC"
                write(6,"(A)") "FCIDUMP"
                write(6,"(A)") "DEBUGOUTPUT"
                write(6,"(A)") "CHECK_SANITY"
                call stop_all(t_r,'Keyword '//trim(w)//' not recognized')
            end select
        enddo Model
                
        if(tMultipleOccs.and.tHalfFill) then
            call stop_all(t_r,'Loops over lattice occupations, as well half-filling only specified in input')
        endif

        if(tSpecifiedAnd.and.tSpecifiedHub) then
            call stop_all(t_r,'Cannot specify both HUBBARD and ANDERSON in MODEL input block')
        endif

        nU_Vals = i
        if(i.ne.0) then
            allocate(U_Vals(nU_Vals))
            U_Vals(1:nU_Vals) = UVals(1:nU_Vals)
        endif

    end subroutine ModelReadInput

    subroutine SCReadInput()
        implicit none
        logical :: teof
        character(len=100) :: w,w2
        character(len=*), parameter :: t_r='SC_ReadInput'
        
        tSC_LR = .true.
        SC: do
            call read_line(teof)
            if(teof) exit
            call readu(w)
            select case(w)
            case("CALC_RE_SPECTRUM")
                tCalcRealSpectrum = .true.
            case("FIT_STYLE")
                call readi(iFitStyle)
            case("FIT_REALFREQ")
                tFitRealFreq = .true.
            case("FIT_ALGO")
                !1 = simplex, 2 = powell
                call readi(iFitAlgo)
                if(item.lt.nitems) then
                    call readf(dFitTol_SC)
                endif
            case("MATSUBARA_FREQ")
                call readf(Start_Omega_Im)
                call readf(End_Omega_Im)
                call readf(Omega_Step_Im)
            case("LATTICE_FIT")
                !1 = normal, 2 = inverse, 3 = normal/norm
                call readi(iLatticeFitType) 
            case("LATTICE_FIT_WEIGHTING")
                call readi(iFitGFWeighting)
            case("LATTICE_OFFDIAG_COUPBLOCKS")
!                call readi(iLatticeCoups)
                call readi(iNonLocBlocks)
            case("REMOVE_COUPS_PRESCHMIDT")
                tRemoveImpCoupsPreSchmidt = .true. !Do we want to remove impurity couplings pre schmidt decomp.?
            case("REMAKE_STATIC_BATH")
                tRemakeStaticBath = .true.
            case("DYNAMIC_GS_OPT")
                tFullReoptGS = .true.
            case("STARTWITHCORRPOT")
                tSC_StartwGSCorrPot = .true.
                if(item.lt.nitems) then
                    call readu(w2)
                    select case(w2)
                    case("FALSE")
                        tSC_StartwGSCorrPot = .false.
                    case default
                        continue
                    end select
                endif
            case("MAXITER_LATTICEFIT")
                call readi(iMaxFitMicroIter)
            case("SC_CONV")
                call readf(dSelfConsConv)
            case("MAXITER_MACROFIT")
                call readi(iMaxIter_MacroFit)
                if(iMaxIter_MacroFit.eq.0) tSkip_Lattice_Fit = .true.
            case("READ_LATTICE_COUPLINGS")
                tReadCouplings = .true.
                if(item.lt.nitems) then
                    call readf(dShiftLatticeEVals)
                endif
            case("READ_CHEMPOT")
                tReadChemPot = .true.   !Read in the initial chemical potential from a file
            case("STRETCH_LAT_BANDWIDTH")
                !Stretch the initial NI hamiltonian bandwidth
                tStretchNILatticeHam = .true.
                call readf(dStretchNILatticeHam)
            case("SKIP_LATTICE_FIT")
                tSkip_Lattice_Fit = .true.
            case("EXT_HAM_FITLATTICE")
                tEnvLatHam = .true.     !Use the fit hamiltonian as the 1-e terms in the external space
            case("BATH_FIT_LAT")
                tEnvLatHam = .true. !By default also in the external space
                tStaticBathFitLat = .true.  !Include fit lattice hamiltonian in the bath space.
            case("EVERYOTHER_LATTICECOUP_ZERO")
                tEveryOtherCoup = .true.    !Every other lattice coupling constrained to be zero
            case("OPT_LATTICE_EVALS")
                tOptGF_EVals = .true.
                tRealSpaceSC = .false.
            case("USE_ANALYTIC_DERIVS")
                tAnalyticDerivs = .true.    !Calculate analytic derivatives in fitting
            case("CONSTRAIN_K_SYM")
                tConstrainKSym = .true.        !Conserve k-point symmetry = 0 in fitting
            case("IMPOSE_K_SYM")
                tImposeKSym = .true.
            case("CONSTRAIN_PH_SYM")
                tConstrainphsym = .true.
            case("IMPOSE_PH_SYM")
                tImposephsym = .true.             !Impose ph sym
            case("LEGENDRE_GRID")           
                call readi(nFreqPoints)
                tFitPoints_Legendre = .true.    !Legendre fitting
            case("REUSE_SELFENERGY")
                call readi(iReuse_SE)
            case("MANYBODY_SELFENERGY")
                tNoHL_SE = .false.
            case("SE_DAMPING")
                call readf(Damping_SE)
            case("POT_UPDATE_DAMPING")
                !Damp the QPSC Updates
                call readf(PotentialUpdateDamping)
            case("SELFENERGY_ITERATIONS")
                call readi(max_SE_iter)
            case("LATTICEGF_MULTIPLIER")
                call readf(NIGF_WeightFac)
            case("READ_SELFENERGY")
                tRead_SelfEnergy = .true.
            case("RANDOM_SELFENERGY")
                tRandom_Init_SE = .true.
            case("SCAN_SELFENERGY_RE")
                tSE_Scan = .true.
                call readf(Start_SE_Real)
                call readf(End_SE_Real)
                call readf(SE_Step_Real)
            case("SCAN_SELFENERGY_IM")
                tSE_Scan = .true.
                call readf(Start_SE_Im)
                call readf(End_SE_Im)
                call readf(SE_Step_Im)
            case("FIT_DIAGONALS")
                tDiagonalSC = .true.
            case("SELF-CONSISTENCY")
                if(item.lt.nitems) then
                    call readi(iGF_Fit)
                endif
            case("END")
                exit
            case default
                write(6,"(A)") "ALLOWED KEYWORDS IN SELF_CONSISTENCY BLOCK: "
                write(6,"(A)") "FIT_ALGO"
                write(6,"(A)") "LATTICE_FIT"
                write(6,"(A)") "LATTICE_FIT_WEIGHTING"
                write(6,"(A)") "LATTICE_COUPLINGS"
                write(6,"(A)") "MAXITER_LATTICEFIT"
                write(6,"(A)") "MAXITER_MACROFIT"
                write(6,"(A)") "READ_LATTICE_COUPLINGS"
                write(6,"(A)") "SKIP_LATTICE_FIT"
                write(6,"(A)") "EXT_HAM_FITLATTICE"
                write(6,"(A)") "EVERYOTHER_LATTICECOUP_ZERO"
                write(6,"(A)") "USE_ANALYTIC_DERIVS"
                write(6,"(A)") "IMPOSE_K_SYM"
                write(6,"(A)") "IMPOSE_PH_SYM"
                write(6,"(A)") "CONSTRAIN_K_SYM"
                write(6,"(A)") "CONSTRAIN_PH_SYM"
                write(6,"(A)") "FIT_DIAGONALS"
                write(6,"(A)") "FIT_REALFREQ"
                write(6,"(A)") "MANYBODY_SELFENERGY"
                write(6,"(A)") "MATSUBARA_FREQ"
                write(6,"(A)") "REUSE_SELFENERGY"
                write(6,"(A)") "SELFENERGY_DAMPING"
                write(6,"(A)") "SELFENERGY_ITERATIONS"
                write(6,"(A)") "LATTICEGF_MULTIPLIER"
                write(6,"(A)") "READ_SELFENERGY"
                write(6,"(A)") "RANDOM_SELFENERGY"
                call stop_all(t_r,'Keyword '//trim(w)//' not recognized')
            end select
        enddo SC

    end subroutine SCReadInput

    subroutine LRReadInput()
        implicit none
        logical :: teof
        character(len=100) :: w
        integer :: Ks(100),i
        character(len=*), parameter :: t_r='LRReadInput'

        !If MR, need full schmidt basis
        i = 0
        LR: do
            call read_line(teof)
            if(teof) exit
            call readu(w)
            select case(w)
            case("DD_RESPONSE")
                tDDResponse = .true.
            case("GF_RESPONSE")
                tChargedResponse = .true.
            case("MOM_GF_RESPONSE")
                tCharged_MomResponse = .true.
            case("KPNT_CALCS")
                call readi(nKCalcs)
            case("K_VALS")
                do while(item.lt.nitems) 
                    i = i+1
                    call readi(Ks(i))
                enddo
            case("NOSTATICBASIS")
                tNoStatickBasis = .true.
            case("NONINT")
                tNIResponse = .true.
            case("TDA")
                tTDAResponse = .true.
            case("RPA")
                tRPAResponse = .true.
            case("CORR_NI_LOCGF")
                tCorrNI_Spectra = .true.
                tCorrNI_LocGF = .true.
            case("CORR_NI_LOCDD")
                tCorrNI_Spectra = .true.
                tCorrNI_LocDD = .true.
            case("CORR_NI_MOMGF")
                tCorrNI_MomGF = .true.
                tCorrNI_Spectra = .true.
            case("EC_TDA")
                tEC_TDA_Response = .true.
                if(item.le.nitems) then
                    call readi(iSolveLR)
                endif
            case("BETA_GF")
                tBetaExcit = .true.
            case("IC_TDA")
                tIC_TDA_Response = .true.
            case("NONDIR_MINRES")
                tMinRes_NonDir = .true.
                if(item.lt.nitems) then
                    call readf(rtol_LR)
                endif
            case("NONDIR_MINRES_AUG")
                tMinRes_NonDir = .true.
                tAugMinRes = .true.
                if(item.lt.nitems) then
                    call readf(rtol_LR)
                endif
            case("NONDIR_GMRES")
                tGMRES_NonDir = .true.
                if(item.lt.nitems) then
                    call readf(rtol_LR)
                endif
            case("MINRES_MAXITER")
                call readi(iMinRes_MaxIter)
            case("REUSE_FIRSTORDER_PSI")
                tReuse_LS = .true.
            case("STORE_HERMIT_HAMIL")
                call stop_all(t_r,'The STORE_HERMIT_HAMIL option has been depreciated since it will not improve efficiency')
            case("PRECONDITION_LR")
                tPreCond_MinRes = .true.
            case("NKRYLOV")
                call readi(nKrylov)
            case("FREQ")
                call readf(Start_Omega)
                call readf(End_Omega)
                call readf(Omega_Step)
            case("BROADENING")
                call readf(dDelta)
            case("NON_NULL")
                tProjectOutNull = .true.
                if(item.lt.nitems) then
                    call readf(MinS_Eigval)
                endif
            case("REOPT_GS")
                tLR_ReoptGS = .true.
            case("EXPLICIT_ORTHOG")
                tExplicitlyOrthog = .true.
            case("WORKLINEARSPAN")
                tOrthogBasis = .true.
            case("REMOVE_GS_FROM_H")
                tRemoveGSFromH = .true.
            case("END")
                exit
            case default
                write(6,"(A)") "ALLOWED KEYWORDS IN LINEAR_RESPONSE BLOCK: "
                write(6,"(A)") "DD_RESPONSE"
                write(6,"(A)") "GF_RESPONSE"
                write(6,"(A)") "MOM_GF_RESPONSE"
                write(6,"(A)") "KPNT_CALCS"
                write(6,"(A)") "K_VALS"
                write(6,"(A)") "NOSTATICBASIS"
                write(6,"(A)") "BETA_GF"
                write(6,"(A)") "NONDIR_MINRES"
                write(6,"(A)") "NONDIR_GMRES"
                write(6,"(A)") "MINRES_MAXITER"
                write(6,"(A)") "PRECONDITION_LR"
                write(6,"(A)") "REUSE_FIRSTORDER_PSI"
                write(6,"(A)") "NKRYLOV"
                write(6,"(A)") "STORE_HERMIT_HAMIL"
                write(6,"(A)") "NONINT"
                write(6,"(A)") "TDA"
                write(6,"(A)") "RPA"
                write(6,"(A)") "CORR_NI_LOCGF"
                write(6,"(A)") "CORR_NI_LOCDD"
                write(6,"(A)") "CORR_NI_MOMGF"
                write(6,"(A)") "EC_TDA"
                write(6,"(A)") "IC_TDA"
                write(6,"(A)") "FREQ"
                write(6,"(A)") "BROADENING"
                write(6,"(A)") "NON_NULL"
                write(6,"(A)") "REOPT_GS"
                write(6,"(A)") "EXPLICIT_ORTHOG"
                write(6,"(A)") "WORKLINEARSPAN"
                write(6,"(A)") "REMOVE_GS_FROM_H"
                call stop_all(t_r,'Keyword '//trim(w)//' not recognized')
            end select
        enddo LR
        
        if(i.ne.0) then
            nKCalcs = i
            allocate(KCalcs(nKCalcs))
            KCalcs(1:nKCalcs) = Ks(1:nKCalcs)
        endif
        
    end subroutine LRReadInput

    subroutine check_input()
        implicit none
        character(len=*), parameter :: t_r='check_input'
        
        !First specify flags which should be also set by the given input
        if(tIC_TDA_Response.or.tEC_TDA_Response) then
            tLR_DMET = .true.
            tConstructFullSchmidtBasis = .true.
        else
            tLR_DMET = .false.
            tConstructFullSchmidtBasis = .false.
        endif
        if(tNIResponse.or.tTDAResponse.or.tRPAResponse) then
            tMFResponse = .true. 
            if(.not.tSCFHF) then
                call stop_all(t_r,'To calculate the full single-reference linear response functions, '  &
                //'it is necessary to calculate full SCF HF solutions')
            endif
        else
            tMFResponse = .false.
        endif
        if(tCorrNI_Spectra) tConstructFullSchmidtBasis = .true.
        if(tCorrNI_MomGF.or.tCharged_MomResponse) then
            !Project the final orbitals onto the original k-space
            !tProjectHFKPnts = .true.
            tKSpaceOrbs = .true.
            tConstructFullSchmidtBasis = .true.
        endif

        !Now check for sanity and implementation of specified options
        if(tReadSystem.or.tExactCorrPot.or.tAnderson) then
            !Ensure we don't do fitting
            tContinueConvergence = .false.
        endif
        if(tExactCorrPot.and.(.not.tHalfFill)) then
            call stop_all(t_r,'Cannot use exact correlation potential apart from 1 impurity, half-filled hubbard models')
        endif
        if(tExactCorrPot.and.tAnderson) then
            call stop_all(t_r,'Cannot use exact correlation potential apart from 1 impurity, half-filled hubbard models')
        endif
        if(tExactCorrPot.and.(nImp.gt.1)) then
            call stop_all(t_r,'Cannot use exact correlation potential apart from 1 impurity, half-filled hubbard models')
        endif
        if((mod(nSites,2).ne.0).and.(LatticeDim.eq.1)) then
            call stop_all(t_r,'Can currently only deal with closed shell systems')
        endif
        if(tOrthogBasis.and.(.not.tProjectOutNull)) then
            call stop_all(t_r,'Need to project out null space of S in order to work in the resulting basis')
        endif
        if(tPeriodic.and.tAntiPeriodic) then
            call stop_all(t_r,'Both PBCs and APBCs specified in input')
        endif
        if(tChemPot.and.(.not.tAnderson)) then
            call stop_all(t_r,'A chemical potential can only be applied to the 1-site Anderson model')
        endif
        if(tLR_DMET.and.(.not.(tDDResponse.or.tChargedResponse.or.tCharged_MomResponse))) then
            call stop_all(t_r,'DMET linear response specified, but type of perturbation not')
        endif
        if(tMFResponse.and.(.not.(tDDResponse.or.tChargedResponse))) then
            call stop_all(t_r,'Single-reference linear response specified, but type of perturbation not')
        endif
        if(tIC_TDA_Response.and.tChargedResponse) then
            call stop_all(t_r,'Linear response greens function not coded up in internally contracted form')
        endif
        if(tIC_TDA_Response.and.tDDResponse) then
            call warning(t_r,'DMET internally-contracted density response broken. Please fix me.')
        endif
        if(tMFResponse.and.(.not.tDDResponse)) then
            call stop_all(t_r,  &
                'A single-reference particle + hole response function asked for, but SR only coded up for Density response')
        endif
        if(tMinRes_NonDir.and.tGMRes_NonDir) then
            call stop_all(t_r,'Cannot specify both MinRes and GMRes algorithms')
        endif
        if(tPrecond_MinRes.and.(.not.(tMinRes_NonDir.or.tGMRES_NonDir))) then
            call stop_all(t_r,'Cannot precondition linear response matrix if not solving iteratively!')
        endif
        if(tDDResponse.and.tSC_LR) then
            call stop_all(t_r,'Self-consistency not yet implemented for density response')
        endif
        if(tAllImp_LR.and.tLR_DMET) then
            call stop_all(t_r,'Calculation of all impurity response functions for correlated spectra still buggy - bug ghb24')
        endif
        if((.not.tCompleteDiag).and.(.not.tNonDirDavidson).and.tLR_DMET) then
            call stop_all(t_r,'To solve DMET_LR, must perform complete diag or non-direct davidson, '   &
     &          //'rather than direct davidson solver')
        endif
        if((iReuse_SE.ne.0).and.(.not.tSC_LR)) then
            call stop_all(t_r,'Cannot reuse self energy if there is no self-consistency in reponse')
        endif
        if((iGF_Fit.ne.0).and.(.not.tSC_LR)) then
            call stop_all(t_r,'How was iGF_Fit set without SC_LR!')
        endif
        if(iGF_Fit.gt.4) then
            call stop_all(t_r,'iGF_Fit set to illegal value')
        endif
        if(tPartialSE_Fit.and.(.not.tSC_LR)) then
            call stop_all(t_r,'Need to specify self-consistency to use option PartialSE_Fit. I know this makes no sense')
        endif
        if(tPartialSE_Fit.and.(iPartialSE_Fit.le.0)) then
            call stop_all(t_r,'Partial self-energy fitting specified, but not the number of fits (should be argument)')
        endif
        if(tReuse_LS.and..not.(tMinRes_NonDir.or.tGMRes_NonDir)) then
            call stop_all(t_r,'Cannot reuse first-order wavefunctions if not using iterative solver')
        endif
        if(tCompressedMats.and.(.not.tNonDirDavidson)) then
            call stop_all(t_r,'Can only use compressed matrices option with non-direct davisdon solver')
        endif
        if((.not.(tMinRes_NonDir.or.tGMRES_NonDir)).and.tLR_DMET.and.tCompressedMats) then
            call stop_all(t_r,'Can only use MinRes/GMRES solver for MR response with compressed matrices')
        endif
        if(tGMRES_NonDir.and..not.tCompressedMats.and.tLR_DMET) then
            call stop_all(t_r,'Can only use GMRES solver with compressed matrices')
        endif
        if(tDDResponse.and.tRemoveGSFromH.and.tCompressedMats) then
            call stop_all(t_r,'REMOVE_GS_FROM_H option cannot be used with compressed matrix DD response')
        endif
        if(tDDResponse.and.tLR_ReoptGS.and.tCompressedMats) then
            call stop_all(t_r,'REOPT_GS option cannot be used with compressed matrix DD response')
        endif
        if(tDDResponse.and.tExplicitlyOrthog.and.tCompressedMats) then
            call stop_all(t_r,'EXPLICIT_ORTHOG option cannot be used with compressed matrix DD response')
        endif
        if(tDDResponse.and.tOrthogBasis.and.tCompressedMats) then
            call stop_all(t_r,'WORKLINEARSPAN option cannot be used with compressed matrix DD response')
        endif
        if(tDDResponse.and.tProjectOutNull.and.tCompressedMats) then
            call stop_all(t_r,'NON_NULL option cannot be used with compressed matrix DD response')
        endif
        if(tReadMats.and..not.tCompressedMats) then
            call stop_all(t_r,'Cannot read matrices if not compressing them')
        endif
        if(tWriteMats.and..not.tCompressedMats) then
            call stop_all(t_r,'Cannot write matrices if not if compressed form')
        endif
        if(tUHF.and..not.tReadSystem) then
            call stop_all(t_r,'UHF currently only working with systems which are read in')
        endif
        if(tUHF.and..not.(tCompleteDiag.or.tNonDirDavidson)) then
            call stop_all(t_r,'Cannot currently solve UHF impurity problem without complete or non-direct davidson solvers')
        endif
        if(tBetaExcit.and.(.not.tChargedResponse)) then
            call stop_all(t_r,'Can only be beta space correlators with GFs - DD response is a spatial orbital')
        endif
        if(tBetaExcit.and.(.not.tCompressedMats)) then
            call stop_all(t_r,'Cannot do beta space correlator without compressed matrices - sorry!')
        endif
        if(tBetaExcit.and.(.not.tUHF)) then
            call stop_all(t_r,'Must use UHF for beta greens functions (otherwise same as alpha space!)')
        endif
        if(tUHF.and.tDDResponse) then
            call stop_all(t_r,'UHF not yet working for DD response')
        endif
        if(tDiag_KSpace.and.tReadSystem) then
            call stop_all(t_r,'Cannot work in k-space if reading in the system')
        endif
        if(tDiag_KSpace.and.(LatticeDim.eq.2).and.(CellShape.eq.1)) then
            call stop_all(t_r,'Cannot currently do k-space diagonalizations for CellShape = 1.')
        endif
        if(tSC_LR.and.(.not.tRealSpaceSC).and.(iNonLocBlocks.ne.0)) then
            call stop_all(t_r,'Off-diagonal lattice coupling length specified, '  &
                //'but not optimizing real space ham')
        endif
!        if(tSC_LR.and.tOptGF_EVals.and.(iLatticeCoups.ne.0)) then
!            call stop_all(t_r,'Lattice eigenvalues self-consistency cannot do a smaller number of lattice sites')
!        endif
        if(tImposephsym.and.(.not.tHalfFill)) then
            call stop_all(t_r,'Can only impose ph sym if we are at half-filling')
        endif
        if(tConstrainphsym.and.(.not.tHalfFill)) then
            call stop_all(t_r,'Can only impose ph sym if we are at half-filling')
        endif
        if((tImposeksym.or.tConstrainksym).and.(.not.tDiag_kSpace)) then
            call stop_all(t_r,'Cannot constrain/impose k-sym without setting up k-sym')
        endif
        if(tImposeksym.and.tConstrainksym) then
            call stop_all(t_r,'No need to impose ksym if we are constraining it in the fit')
        endif
        if(tImposephsym.and.tConstrainphsym) then
            call stop_all(t_r,'No need to impose ph sym if we are constraining it in the fit')
        endif
!        if((tImposephsym.or.tConstrainphsym.or.tImposeksym.or.tConstrainksym).and.(.not.tOptGF_EVals)) then
!            call stop_all(t_r,'Cannot impose symmetries unless optimizing eigenvalues')
!        endif
        if(tConstrainphsym.and.(.not.tConstrainksym)) then
            call stop_all(t_r,'Sorry - you cant constrain ph sym without also constraining ksym')
        endif
!        if(tRealSpaceSC.and.(tImposeksym.or.tConstrainksym)) then
!            call stop_all(t_r,'If doing real-space optimization, then k-symmetry is constrained automatically')
!        endif
        if(tRealSpaceSC.and.(tConstrainphsym)) then
            call stop_all(t_r,'Cannot constrain ph sym if doing real-space optimization')
        endif
        if(tSC_LR.and.tRemakeStaticBath.and.(.not.tLR_ReoptGS)) then
            call stop_all(t_r,'If remaking static bath space, should really reoptimize GS too')
        endif
        if((LatticeDim.eq.2).and.(CellShape.eq.0)) then
            call stop_all(t_r,'No Unit Cell shape specified for 2D square lattice')
        endif
        if(tMP2.and..not.tSCFHF) then
            call stop_all(t_r,'To do MP2, you must do a full HF (SCF_HF option)')
        endif

    end subroutine check_input

    subroutine name_timers()
        implicit none

        !From main subroutine
        Full_timer%timer_name='Main'
        FullSCF%timer_name='FullSCF'
        FCIDUMP%timer_name='FCIDUMP'
        DiagT%timer_name='DiagHopping'
        ConstEmb%timer_name='Const_Emb'
        Trans1e%timer_name='1eTransform'
        HL_Time%timer_name='HL_Solver'
        Fit_v_time%timer_name='Fit_corrpot'

        !SR_LR
        LR_SR_NonInt%timer_name='LR_SR_NonInt'
        LR_SR_TDA%timer_name='LR_SR_TDA'
        LR_SR_RPA%timer_name='LR_SR_RPA'

        !MR_LR_EC
        !Density response
        LR_EC_TDA_Precom%timer_name='DD_EC_Precom'
        LR_EC_TDA_HBuild%timer_name='DD_EC_HBuild'
        LR_EC_TDA_SBuild%timer_name='DD_EC_SBuild'
        LR_EC_TDA_Project%timer_name='DD_EC_NullProj'
        LR_EC_TDA_OptGS%timer_name='DD_EC_SolveGS'
        LR_EC_TDA_BuildLR%timer_name='DD_EC_BuildLR'
        LR_EC_TDA_SolveLR%timer_name='DD_EC_SolveLR'
        !GF response
        LR_EC_GF_Precom%timer_name='GF_EC_Precom'
        LR_EC_GF_HBuild%timer_name='GF_EC_HBuild'
        LR_EC_GF_OptGS%timer_name='GF_EC_OptGS'
        LR_EC_GF_SolveLR%timer_name='GF_EC_SolveLR'
        LR_EC_GF_FitGF%timer_name='GF_EC_FitGF'
        !self-consistency
        SelfCon_LR%timer_name='SelfConsistency'
        FitLatHam%timer_name='FitLatHam'
        CalcLatSpectrum%timer_name='CalcSPSpectrum'
        CalcGrads%timer_name='CalcGrads'

        !GF2
        GF2_time%timer_name = 'GF2'
        BuildMatGF_time%timer_name = 'BuildLatGF'
        FT_MatToTau_time%timer_name = 'FT_MatToTau'
        FT_TauToMat_time%timer_name = 'FT_TauToMat'
        GMEnergy_time%timer_name = 'GMEnergy'
        BuildSE_time%timer_name = 'BuildSE'
        ConvergeMu_time%timer_name = 'ConvMu'
        FitSplines_time%timer_name = 'FitSplines'

    end subroutine name_timers



end module readinput
