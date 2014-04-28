Program RealHub
    use const
    use timing
    use errors, only: stop_all
    use globals
    use Solvers
    use SchmidtDecomp
    use LRDriver 
    use fitting
    use mat_tools
    use input
    implicit none

    call init_calc()
    call set_timer(full_timer)
    call run_DMETcalc()
    call halt_timer(full_timer)
    call end_calc()

    contains
    
    !These want to be input parameters at some point
    subroutine set_defaults()
        implicit none

        !Main options
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

        !General LR options
        Start_Omega = 0.0_dp
        End_Omega = 4.05_dp
        Omega_Step = 0.1_dp
        dDelta = 0.01_dp

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
        nKCalcs = 0
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
        tCalcRealSpectrum = .false.

    end subroutine set_defaults

    subroutine init_calc()
        use report, only: environment_report
        use timing, only: init_timing
        use utils, only: get_free_unit
        implicit none
        real(dp) :: U_tmp
        integer :: i,minres_unit
        logical :: exists
        character(len=*), parameter :: t_r='init_calc'

        write(6,"(A)") "***  Starting real-space hubbard/anderson calculation  ***"

        call init_timing()

        call name_timers()

        call environment_report()

        call set_defaults()

        call read_input()

        call check_input()

        call check_openmp()
        
        inquire(file='zMinResQLP.txt',exist=exists)
        if(exists) then
            minres_unit = get_free_unit()
            open(minres_unit,file='zMinResQLP.txt',status='old')
            close(minres_unit,status='delete')
        endif

        if(LatticeDim.eq.2) then
            call Setup2DLattice()
        endif

        if(tAnderson) then
            write(6,"(A)") "Running:    o Anderson Model (single site)"
        elseif(tReadSystem) then
            write(6,"(A)") "Reading unknown system from files..."
        elseif(tSingFiss) then
            write(6,"(A)") "Reading Singlet Fission parameters..."
        else
            write(6,"(A)") "Running:    o Hubbard Model"
        endif
        if(tChemPot) then
            write(6,"(A)") "            o Chemical potential of -U/2 at impurity site for interacting system"
        else
            write(6,"(A)") "            o No chemical potential applied at impurity site"
        endif
        if(.not.tReadSystem) then
            if(tPeriodic) then
                write(6,"(A)") "            o PBCs employed"
            elseif(tAntiPeriodic) then
                write(6,"(A)") "            o APBCs employed"
            elseif(.not.(tPeriodic.or.tAntiPeriodic)) then
                write(6,"(A)") "            o Open boundary conditions employed"
            endif
            if(LatticeDim.eq.2) then
                write(6,"(A)") "            o 2-dimensional model" 
                write(6,"(A,I7,A,I7,A)") "            o Size of lattice: ",TDLat_Ni,' x ',TDLat_Nj,' at 45 degrees' 
                write(6,"(A,I7)") "            o Total lattice sites: ",nSites
            elseif(LatticeDim.eq.1) then
                write(6,"(A)") "            o 1-dimensional model" 
                write(6,"(A,I7)") "            o Size of lattice: ",nSites 
            else
                call stop_all(t_r,"Cannot determine dimensionality of system")
            endif
        else
            write(6,"(A,I9)") "            o Number of lattice sites: ",nSites
        endif
        if(tUHF) then
            write(6,"(A)") "            o *Unrestricted* bath construction: (Anti-)Ferromagnetic phase"
        else
            write(6,"(A)") "            o *Restricted* bath construction: Paramagnetic phase"
        endif
        write(6,"(A)") "            o Range of U values to consider: " 
        if(nU_Vals.eq.0) then
            !Sweeping through
            U_tmp=StartU
            do while((U_tmp.lt.max(StartU,EndU)+1.0e-5_dp).and.(U_tmp.gt.min(StartU,EndU)-1.0e-5_dp))
                write(6,"(A,F10.5)") "            o U = ",U_tmp 
                U_tmp=U_tmp+UStep
            enddo
        else
            do i = 1,nU_Vals
                write(6,"(A,F10.5)") "            o U = ",U_Vals(i) 
            enddo
        endif
        if(.not.tReadSystem) then
            if(tReadInCorrPot) then
                write(6,"(A,A)") "            o Correlation potentials for system will be read from file: ", trim(CorrPot_file)
                write(6,"(A)") "            o No DMET self-consistency of correlation potential"
            else
                if(.not.tAnderson) then
                    write(6,"(A,I8)") "            o Maximum iterations for DMET self-consistency: ", iMaxIterDMET
                    if(tSaveCorrPot) then
                        if(tHalfFill) then
                            write(6,"(A,I8)") "            o The correlation potential from the previous U value will " &
                                //"be used as a starting point for self-consistency" 
                        else
                            write(6,"(A,I8)") "            o The correlation potential from the previous electron number will " &
                                //"be used as a starting point for self-consistency"
                        endif
                    endif
                endif
            endif
        endif
        if(tHalfFill) then
            write(6,"(A)") "            o Only half filling to be considered"
        else
            if(tRampDownOcc) then
                write(6,"(A)") "            o Ramping down from half filling occupation of lattice"
            else
                write(6,"(A)") "            o Ramping up to half filling occupation of lattice"
            endif
        endif
        if(tSCFHF.and..not.tReadSystem) then
            write(6,"(A)") "            o Full self-consistent Hartree--Fock orbitals will be calculated"
        endif
        write(6,"(A,I7)") "            o Number of impurity sites: ",nImp 
        if(tDumpFCIDump) then
            write(6,"(A)") "            o Creating FCIDUMPs for system" 
        endif
        if(tFCIQMC) then
            if(nNECICores.eq.0) then
                write(6,"(A)") "            o Impurity solver: FCIQMC via call to serial NECI code"
            else
                write(6,"(A,I6)") "            o Impurity solver: FCIQMC via call to parallel NECI code. nCores: ",nNECICores
            endif
        elseif(tNonDirDavidson) then
            write(6,"(A)") "            o Impurity solver: Non-direct iterative Davidson diagonalizer" 
        elseif(tCompleteDiag) then
            write(6,"(A)") "            o Impurity solver: Complete diagonalization"
        else
            write(6,"(A)") "            o Impurity solver: Direct Davidson diagonalizer via call to FCI code"
        endif
        if(tNIResponse) then
            write(6,"(A)") "            o Calculating non-interacting linear response function" 
        endif
        if(tTDAResponse) then
            write(6,"(A)") "            o Calculating TDA linear response function" 
        endif
        if(tRPAResponse) then
            write(6,"(A)") "            o Calculating RPA linear response function" 
        endif
        if(tEC_TDA_Response) then
            write(6,"(A)") "            o Calculating externally-contracted MC-TDA DMET linear response function" 
            if(tDDResponse) then
                write(6,"(A)") "                o Local density response calculated" 
            endif
            if(tChargedResponse) then
                write(6,"(A)") "                o Local Greens function calculated" 
            endif
            if(tCharged_MomResponse) then
                write(6,"(A)") "                o Momentum-resolved Greens functions calculated" 
            endif

        endif
        if(tIC_TDA_Response) then
            write(6,"(A)") "            o Calculating internally-contracted MC-TDA DMET linear response function" 
        endif
        if(tMFResponse.or.tLR_DMET) then
            write(6,"(A,F13.8)") "            o Spectral broadening for linear response functions: ",dDelta
        endif
        write(6,"(A)") ""
        write(6,"(A)") ""


        
    end subroutine init_calc

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
            case("PBC")
                tPeriodic = .true.
            case("APBC")
                tAntiPeriodic = .true.
            case("TEMPERATURE")
                call readf(Temperature)
                tThermal = .true.
            case("MAXITER_DMET")
                call readi(iMaxIterDMET)
                if(item.lt.nitems) then
                    call readf(dTolDMET)
                endif
            case("SCF_HF")
                tSCFHF = .true.
            case("KSPACE_DIAG")
                tDiag_kspace = .true.
            case("IMPSITES")
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
        integer :: Ks(100),i
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
            case("MAXITER_MACROFIT")
                call readi(iMaxIter_MacroFit)
            case("READ_LATTICE_COUPLINGS")
                tReadCouplings = .true.
                if(item.lt.nitems) then
                    call readf(dShiftLatticeEVals)
                endif
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
            case("DAMPING")
                call readf(Damping_SE)
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
        if(tReadSystem.or.tExactCorrPot) then
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
        if(tCoreH_EmbBasis.and.(tNonDirDavidson.or.tCompleteDiag)) then
            call stop_all(t_r,'Cannot use CoreH_EmbBasis transform with complete diag or nondir-davidson')
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
        if(tDiag_KSpace.and.LatticeDim.eq.2) then
            call stop_all(t_r,'Cannot currently do k-space diagonalizations for 2D models. Fix me!')
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

    end subroutine name_timers

    subroutine end_calc()
        use timing, only: end_timing, print_timing_report
!!$      use omp_lib
        implicit none
        real(dp) :: OpenMP_wallend 
            
        call deallocate_mem()
        call end_timing()
        call print_timing_report()
!        if(tOpenMP) then
!!$          OpenMP_wallend = OMP_get_wtime()
!            write(6,"(A)") ""
!            write(6,"(A,F14.2)") "Total calculational walltime per OpenMP thread: ",  &
!                OpenMP_wallend - OpenMP_wallstart
!        endif

    end subroutine end_calc

    !Deallocate memory which is still allocated
    subroutine deallocate_mem()
        implicit none

        if(allocated(U_Vals)) deallocate(U_Vals)
        if(allocated(TD_Imp_Lat)) deallocate(TD_Imp_Lat,TD_Imp_Phase)

    end subroutine deallocate_mem

    subroutine Setup2DLattice()
        implicit none
        real(dp) :: dWidth
        integer :: TDLat_Width,x,y,dx,dy,ci,cj,site_imp
        integer :: i,j,k
        character(len=*), parameter :: t_r='Setup2DLattice'

        !Work out an appropriate width, and an actual nSites
        dWidth = sqrt(nSites*2.0_dp)
        TDLat_Width = 2 * nint(dWidth/2.0_dp)

        !There are two lattices, an x, y lattice, and an i, j lattice. These have had their axes rotated by 45 degrees.
        !2DLat_Ni is the width of each lattice (there are two interlocking in the (i,j) representation
        TDLat_Ni = TDLat_Width / 2
        !2DLat_Nj is the width of the lattice in the (x,y) representation
        TDLat_Nj = TDLat_Width

        !Actual nSites. 
        nSites = TDLat_Ni * TDLat_Nj

        write(6,*) "Updated number of sites in the 2D hubbard model will be: ",nSites
        
        if(mod(TDLat_Width/2,2).eq.0) then
            !Use anti-periodic boundary conditions
            !HF ground state only unique if Width=2*odd_number (direct PBC)
            !or if Width=2*even_number (anti-PBC)
            tAntiperiodic = .true.
            tPeriodic = .false.
        else
            tAntiperiodic = .false.
            tPeriodic = .true.
        endif
        if(tPeriodic) then
            write(6,*) "Periodic boundary conditions now in use"
        else
            write(6,*) "Anti-Periodic boundary conditions now in use"
        endif

        !Now to set up the impurity
        if(nImp.eq.1) then
            nImp_x = 1
            nImp_y = 1
        elseif(nImp.eq.2) then
            nImp_x = 1
            nImp_y = 2
        elseif(nImp.eq.4) then
            nImp_x = 2
            nImp_y = 2
        else
            call stop_all(t_r,'Cannot deal with impurities > 4')
        endif

        !Find the x,y coordinates for the middle of the array. This will be used to
        !define the corner site of the impurity
        call ij2xy(TDLat_Ni/2,TDLat_Nj/2,x,y)

        !Setup the impurity space, and how to tile the impurity through the space.
        !This creates the matrices TD_Imp_Lat and TD_Imp_Phase
        !If the correlation potential is a matrix of nImp x nImp, then TD_Imp_Lat
        !gives the index of that correlation potential which corresponds to the tiled
        !correlation potential through the space.
        !(TD_Imp_Lat(site,site)-1)/nImp + 1 gives the impurity index of that site. 

        call MakeVLocIndices()

        !Create the impurity cluster
        allocate(ImpSites(nImp))
        ImpSites = 0
        do dx = 0,nImp_x-1
            do dy = 0,nImp_y-1
                call xy2ij(x+dx,y+dy,ci,cj)
                !Remember - our sites are 1 indexed
                site_imp = ci + TDLat_Ni*cj + 1     !No need to take mods. We certainly shouldn't exceed the bounds of the ij lattice
                !write(6,*) "***",dx,dy,site_imp,(TD_Imp_Lat(site_imp,site_imp)),((TD_Imp_Lat(site_imp,site_imp)-1)/nImp) + 1
                ImpSites(((TD_Imp_Lat(site_imp,site_imp)-1)/nImp) + 1) = site_imp
            enddo
        enddo

        write(6,*) "Impurity sites defined as: ",ImpSites(:)

        !We now want to define a mapping, from the standard site indexing to an impurity ordering of the sites, such that
        ! Perm_indir(site) = Imp_site_index
        ! Perm_dir(Imp_site_index) = site       The first index maps you onto the original impurity sites
        ! Therefore, Perm_indir(site)%nImp = impurity site it maps to which repeat of the impurity you are on

        !In the impurity ordering (indirect), the impurity sites are first, followed by the repetitions of the striped
        !impurity space.
        !The direct space is the normal lattice ordering
        allocate(Perm_indir(nSites))
        allocate(Perm_dir(nSites))
        Perm_indir(:) = 0
        Perm_dir(:) = 0
        Perm_dir(1:nImp) = ImpSites(:)
        k = nImp+1
        loop: do i = 1,nSites
            do j = 1,nImp 
                if(i.eq.ImpSites(j)) then
                    cycle loop
                endif
            enddo
            Perm_dir(k) = i
            k = k + 1
        enddo loop
        if(k.ne.nSites+1) call stop_all(t_r,"Error here")

        do i=1,nSites
            Perm_indir(Perm_dir(i)) = i
        enddo
        !write(6,*) "Perm_dir: ",Perm_dir(:)
        !write(6,*) "Perm_indir: ",Perm_indir(:)

    end subroutine Setup2DLattice

    subroutine run_DMETcalc()
        implicit none
        real(dp) :: ElecPerSite,FillingFrac,VarVloc,ErrRDM,mean_vloc
        integer :: i,it,Occ,CurrU,DMETfile
        logical :: tFinishedU,t2RDM
        character(len=*), parameter :: t_r="run_DMETcalc"

        !Set up initial conditions, i.e. starting potential
        allocate(v_loc(nImp,nImp))
        allocate(h0(nSites,nSites))     !The core hamiltonian
        allocate(h0v(nSites,nSites))    !The core hamiltonian with the local potential
        allocate(HFEnergies(nSites))    !The fock energies
        v_loc(:,:) = zero 
        h0(:,:) = zero 
        h0v(:,:) = zero 
        HFEnergies(:) = zero 
        if(tUHF) then
            allocate(v_loc_b(nImp,nImp))
            allocate(h0_b(nSites,nSites))
            allocate(h0v_b(nSites,nSites))
            allocate(HFEnergies_b(nSites))
            v_loc_b(:,:) = zero 
            h0_b(:,:) = zero 
            h0v_b(:,:) = zero 
            HFEnergies_b(:) = zero 
        endif

        CurrU = 0
        do while(.true.)

            !Find the next U value
            call GetNextUVal(CurrU,tFinishedU)
            if(tFinishedU) exit
            if(.not.tSingFiss) write(6,*) "Running DMET calculation with U = ",U
        
            !Calculate the core hamiltonian based on the hopping matrix of the hubbard model in real space
            !If reading in the hopping matrix, it is done here and stored in h0
            call make_hop_mat()

            if((tDiag_KSpace.or.tProjectHFKPnts.or.tKSpaceOrbs).and.(.not.allocated(KPnts))) then
                call setup_kspace()
            endif

            !Diagonalize the mean-field hamiltonian
            !Get occupations with unique GS
            call find_occs()

            !Loop over occupation numbers 
            do Occ=1,N_Occs

                allocate(MeanFieldDM(nSites,nSites))    !DM from mean-field
                MeanFieldDM(:,:) = zero
                if(tUHF) then
                    allocate(MeanFieldDM_b(nSites,nSites))
                    MeanFieldDM_b(:,:) = zero
                endif

                call OpenDMETFile(DMETfile)
                write(DMETfile,"(A)") " #Iteration  E_DMET/Imp   E_HL   d[V]   Initial_Err[1RDM]   "    &
                    //"Filling   Filling_Err   mean_diag_correlation"

                !These occupations refer to number of closed shell orbitals, so total electrons is 2 x nOcc
                nOcc = allowed_occs(Occ)    !Number of occupied orbitals in this iteration
                NEl = 2*nOcc    !Total electrons in system
                ElecPerSite = (2.0_dp*nOcc)/real(nSites,dp) !Average number of electrons per site
                FillingFrac = ElecPerSite/2.0_dp    !A filling fraction of 1/2 is half filling (duh)

                !Write out some stats
                write(6,*) 
                write(6,"(A,F8.3,A,I5,A,I5,A)") "Electrons per site:   ",ElecPerSite," (in ", &
                    nOcc," doubly occupied orbitals on ",nSites," sites)"
                write(6,"(A,F10.5)")          "Filling Fraction:     ",FillingFrac
                write(6,"(A,F8.3)")           "Hubbard U:            ",U
                write(6,"(A,I5,A)")           "Embedded system size: ",nImp," sites"
                if(tAnderson) then
                    write(6,"(A,I5,A)")           "Anderson lattice of ",nSites," sites"
                else
                    write(6,"(A,I5,A)")           "1D Hubbard lattice of ",nSites," sites"
                endif
                write(6,*) 
                
                if(tDiagFullSystem) then
                    !Find all eigenvalues of system
                    call DiagFullSystem()
                endif

                !Calculate full hf, including mean-field on-site repulsion (which is included in correlation potential in DMET
                call set_timer(FullSCF)
                if(tSCFHF) then
                    call run_true_hf()
                else
                    call run_hf(0)
                endif
                call halt_timer(FullSCF)

                if(tDumpFCIDUMP) then
                    call set_timer(FCIDUMP)
                    call DumpFCIDUMP()
                    call halt_timer(FCIDUMP)
                endif

                !Calculate single reference linear response - non-interacting, TDA and RPA
                if(tMFResponse) then
                    call SR_LinearResponse()
                endif

                if(tReadInCorrPot.or.tReadSystem) then
                    !Read in the correlation potential from another source
                    call read_in_corrpot()
                elseif(tExactCorrPot) then
                    !We know the exact correlation potential for the hubbard model: U/2
                    v_loc(:,:) = U/2.0_dp
                endif

                !At this point, we have h0, U and a set of system sites (the first nImp indices), as well as a local potential
                do it=1,iMaxIterDMET

                    !Write out header file for the convergence
                    write(6,"(A,I6)") "Iteration: ",it

                    !Do iMaxIterDMET microiterations to converge the DMET for this occupation number
                    call add_localpot(h0,h0v,v_loc,core_b=h0_b,core_v_b=h0v_b,CorrPot_b=v_loc_b)

                    !Now run a HF calculation by constructing and diagonalizing the fock matrix
                    !This will also return the RDM in the AO basis
                    call set_timer(DiagT)
                    if(tReadSystem) then
                        call read_orbitals()
                    else
                        call run_hf(it)
                    endif
                    call halt_timer(DiagT)

                    !Construct the embedded basis
                    call set_timer(ConstEmb)
                    if(tThermal) then
                        call ConstructThermalEmbedding()
                    else
                        if(tConstructFullSchmidtBasis) then
                            call ConstructFullSchmidtBasis()
                        else
                            call CalcEmbedding()
                        endif
                    endif
                    if(tWriteOut) then
                        call writematrix(EmbeddedBasis,'EmbeddedBasis',.true.)
                        if(tUHF) call writematrix(EmbeddedBasis_b,'BetaEmbeddedBasis',.true.)
                    endif
                    call halt_timer(ConstEmb)
                    
                    !Now transform the 1 electron quantities into the embedded basis
                    !This should be exactly correct, i.e. we can now diagonalize the fock matrix in this basis
                    !to get the same result. We could also check that the number of electrons adds to the correct thing
                    call set_timer(Trans1e)
                    call Transform1e()
                    call halt_timer(Trans1e)
                    
                    call set_timer(HL_Time)
                    !Construct the two electron integrals in the system, and solve embedded system with high-level method
                    t2RDM = .false.
                    if(tFCIQMC) t2RDM = .false.
                    call SolveSystem(t2RDM)
                    call halt_timer(HL_Time)

                    if((.not.tAnderson).and.tContinueConvergence) then
                        !Fit new potential
                        !vloc_change (global) is updated in here to reflect the optimal change
                        !VarVloc is a meansure of the change in the potential
                        !ErrRDM is a measure of the initial difference in the RDMs
                        call set_timer(Fit_v_time)
                        call Fit_vloc(VarVloc,ErrRDM)

                        if(tDebug) call writematrix(vloc_change,'vloc_change',.true.)

                        !Mean vloc is actually for the old vloc for consistency with Geralds code
                        mean_vloc = 0.0_dp
                        do i=1,nImp
                            mean_vloc = mean_vloc + v_loc(i,i)
                        enddo
                        mean_vloc = mean_vloc/real(nImp)
                        call halt_timer(Fit_v_time)

                        !Write out stats:
                        !   Iter    E/Site  d[V]    Initial_ERR[RDM]    ERR[Filling]    mean[corr_pot]      Some RDM stuff...?
                        write(6,"(I7,5G22.10)") it,TotalE_Imp,VarVloc,ErrRDM,FillingError,mean_vloc
                        write(DMETfile,"(I7,7G22.10)") it,TotalE_Imp,HL_Energy,VarVloc,ErrRDM,  &
                            Actualfilling_Imp,FillingError,mean_vloc

                        !Update vloc
                        v_loc(:,:) = v_loc(:,:) + vloc_change(:,:)

                        if(VarVloc.lt.dTolDMET) then
                            write(6,"(A)") "...correlation potential converged" 
                            exit
                        endif
                    else
                        !Write out stats:
                        !   Iter    E/Site  d[V]    ERR[RDM]    ERR[Filling]    mean[corr_pot]      Some RDM stuff...?
                        VarVloc = 0.0_dp
                        ErrRDM = 0.0_dp
                        mean_vloc = 0.0_dp

                        write(6,"(I7,5G22.10)") it,TotalE_Imp,VarVloc,ErrRDM,FillingError,mean_vloc
                        write(DMETfile,"(I7,7G22.10)") it,TotalE_Imp,HL_Energy,VarVloc,ErrRDM,  &
                            Actualfilling_Imp,FillingError,mean_vloc
                        call flush(6)

                        exit    !Anderson model, so we do not want to iterate
                    endif

                enddo   !DMET convergence

                if(it.gt.iMaxIterDMET) call warning(t_r,'DMET Convergence failed - try increasing MAXITER_DMET ?')
                    
                !Potentially run FCI again now to get correlation functions from 2RDMs?
                write(6,"(A,F10.4,A,G20.10)") "FINAL energy per site for U=",U,' is: ',TotalE_Imp
                call flush(6)
                
                if(.not.tAnderson) then
                    close(DMETfile)
                    if(tHalfFill) then
                        call WriteCorrPot()
                    endif
                    call writematrix(v_loc,'Converged Correlation Potential',.true.)
                endif
        
                deallocate(MeanFieldDM)
                if(tUHF) deallocate(MeanFieldDM_b)

                if(tProjectHFKPnts) then
                    call ProjectHFontoK()
                endif

                if(tKSpaceOrbs) then
                    call GetKSpaceOrbs()
                endif

                if(tCorrNI_Spectra) then
                    !Calculates single reference spectral functions using the correlated 1-electron hamiltonian
                    call Correlated_SR_LR()
                endif

                if(tLR_DMET) then
                    !Perform linear response on the resulting DMET state
                    call MR_LinearResponse()
                endif
                deallocate(HFOrbs)
                if(tUHF) deallocate(HFOrbs_b)

                !Set potential for the next occupation number, or wipe it?
                if(.not.tSaveCorrPot) then
                    v_loc(:,:) = zero
                    if(tUHF) v_loc_b(:,:) = zero
                endif

            enddo   !Loop over occupations

            if(.not.tHalfFill) then
                !Wipe correlation potential if we have ramped through occupations
                !We can potentially keep it though if we are just doing half-filling
                v_loc(:,:) = zero
                if(tUHF) v_loc_b(:,:) = zero
            endif

        enddo   !Loop over U values

    end subroutine run_DMETcalc

    !CurrU is zero on entry for the first time, and is returned as non-zero
    !tFinished = true when we have run through all U values.
    subroutine GetNextUVal(CurrU,tFinished)
        implicit none
        logical, intent(out) :: tFinished
        integer, intent(inout) :: CurrU
        character(len=*), parameter :: t_r='GetNextUVal'

        tFinished = .false.

        if((CurrU.ne.0).and.(tSingFiss)) then
            tFinished = .true.
            return !We don't want to loop through U!
        endif

        if(nU_Vals.eq.0) then
            !We are sweeping through U values, rather than specifying them individually
            if(CurrU.eq.0) then
                !First value
                U = StartU
                CurrU = 1   !So we don't go into this block again
            else
                !Carry on sweeping through.
                !Increment from the last value
                U = U + UStep
                !Find out if we are still in range
                if((U.gt.max(StartU,EndU)+1.0e-5_dp).or.(U.lt.min(StartU,EndU)-1.0e-5_dp)) then
                    !We have finished
                    tFinished = .true.
                endif
            endif
        else
            !We are running through specified U values
            if(.not.allocated(U_Vals)) call stop_all(t_r,'U_Vals array not allocated')
            CurrU = CurrU + 1
            if(CurrU.gt.nU_Vals) then
                !We have run through all U values we want
                tFinished = .true.
            else
                U = U_Vals(CurrU)
            endif
        endif

    end subroutine GetNextUVal

    subroutine WriteCorrPot()
        use utils, only: get_free_unit,append_ext_real
        implicit none
        integer :: iunit
        character(64) :: filename
!        character(len=*), parameter :: t_r='WriteCorrPot'

        write(6,*) "Writing out converged correlation potential..."

        call append_ext_real("CORRPOTS",U,filename)
        iunit = get_free_unit()
        open(iunit,file=filename,status='unknown')

        write(iunit,*) U,nOcc,v_loc(:,:)
        close(iunit)

    end subroutine WriteCorrPot

    subroutine read_in_corrpot()
        use utils, only: get_free_unit
        use DetTools, only: tospat
        implicit none
        integer :: iunit,Occ_val,ios,i,j,k
        logical :: texist,tFoundCorrPot
        real(dp) :: U_val,CorrPot_tmp(nImp*nImp)
        real(dp), allocatable :: temp(:)
        character(len=*), parameter :: t_r='read_in_corrpot'

        write(6,*) "Reading in correlation potential..."

        if(tReadSystem) then

            inquire(file='FinalVCorr.dat',exist=texist)
            if(.not.texist) call stop_all(t_r,'Correlation potential file cannot be found') 
            iunit = get_free_unit()
            open(iunit,file='FinalVCorr.dat',status='old',action='read')
            if(tUHF) then
                allocate(temp(nImp*2))
                do i = 1,nImp*2
                    read(iunit,*) temp(:)
                    do j = 1,nImp*2
                        if(mod(i,2).eq.1) then
                            !i is alpha spin-orbital
                            if((mod(j,2).eq.0).and.(abs(temp(j)).gt.1.0e-8_dp)) then
                                call stop_all(t_r,'Coupling between different spin types in correlation potential?!')
                            elseif(mod(j,2).eq.1) then
                                !j is also alpha
                                v_loc(tospat(i),tospat(j)) = temp(j)
                            endif
                        else
                            !i is beta spin-orbital
                            if((mod(j,2).eq.1).and.(abs(temp(j)).gt.1.0e-8)) then
                                call stop_all(t_r,'Coupling between different spin types in correlation potential?!')
                            elseif(mod(j,2).eq.0) then
                                !j is also beta
                                v_loc_b(tospat(i),tospat(j)) = temp(j)
                            endif
                        endif
                    enddo
                enddo
                deallocate(temp)
            else
                do i = 1,nImp
                    read(iunit,*) v_loc(i,:)
                enddo
            endif
            close(iunit)
        else
            iunit = get_free_unit()
            inquire(file=CorrPot_file,exist=texist)
            if(.not.texist) then
                write(6,*) "correlation potential filename: ",CorrPot_file
                call stop_all(t_r,'Expecting to read in a file with a converged '    &
     &              //'correlation potential, but unable to find appropriate file')
            endif

            open(iunit,file=CorrPot_file,status='old')

            tFoundCorrPot = .false.
            do while(.true.)
                read(iunit,*,iostat=ios) U_val,Occ_val,CorrPot_tmp(1:nImp*nImp)
                if(ios.gt.0) call stop_all(t_r,'Error reading in correlation potential')
                if(ios.lt.0) exit   !EOF
                if((abs(U_val-U).lt.1.0e-7_dp).and.(Occ_val.eq.nOcc)) then
                    tFoundCorrPot = .true.
                    exit
                endif
            enddo

            if(.not.tFoundCorrPot.and.(.not.tContinueConvergence)) then
                call stop_all(t_r,'Did not read in correlation potential corresponding to this run')
            else
                k=1
                do i=1,nImp
                    do j=1,nImp
                        v_loc(j,i) = CorrPot_tmp(k)
                        k = k+1
                    enddo
                enddo
            endif

            close(iunit)
        endif

        write(6,"(A)") "Read in correlation potential: "
        call writematrix(v_loc,"v_loc",.true.)
        do i=1,nImp
            do j=1,nImp
                if(abs(v_loc(i,j)-v_loc(j,i)).gt.1.0e-6_dp) then
                    call stop_all(t_r,'correlation potential not symmetric')
                endif
            enddo
        enddo
        if(tUHF) then
            !Test beta component too
            call writematrix(v_loc_b,"v_loc_b",.true.)
            do i=1,nImp
                do j=1,nImp
                    if(abs(v_loc_b(i,j)-v_loc_b(j,i)).gt.1.0e-6_dp) then
                        call stop_all(t_r,'beta correlation potential not symmetric')
                    endif
                enddo
            enddo
        endif

    end subroutine read_in_corrpot

    !Open an output file for the DMET convergence
    subroutine OpenDMETFile(iunit)
        use utils, only: get_free_unit,append_ext,append_ext_real
        implicit none
        integer, intent(out) :: iunit
        character(64) :: filename,filename2

        iunit = get_free_unit()
        call append_ext_real("DMET",U,filename)
        if(.not.tHalfFill) then
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')

    end subroutine OpenDMETFile
    
    subroutine DumpFCIDUMP()
        use utils, only: get_free_unit,append_ext_real
        use DetTools, only: GetHFInt_spinorb
        implicit none
        integer :: iunit,i,j,k,l,A,B,ex(2,2)
        real(dp) :: hel
        real(dp), allocatable :: temp(:,:),h0HF(:,:)
        character(len=64) :: filename
        
        iunit = get_free_unit()
        call append_ext_real('FCIDUMP',U,filename)
        open(unit=iunit,file=filename,status='unknown')
        write(iunit,'(2A6,I3,A9,I3,A6,I2,A)') '&FCI ','NORB=',nSites,', NELEC=',NEl,', MS2=',0,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,nSites
            write(iunit,'(I1,A1)',advance='no') 1,','
        enddo
        write(iunit,*) ""
        WRITE(iunit,'(A7,I1)') 'ISYM=',1
        WRITE(iunit,'(A5)') '&END'

        !Calculate hopping matrix in MO basis
        if(allocated(FullHFOrbs)) then
            write(6,'(A)') "Writing out FCIDUMP file corresponding to hartree--fock orbitals"

            allocate(temp(nSites,nSites))
            allocate(h0HF(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,h0HF,nSites)
            deallocate(temp)

            do i=1,nSites
                do j=1,nSites
                    A=(i*(i-1))/2+j
                    DO k=1,nSites
                        DO l=1,nSites
                            B=(k*(k-1))/2+l

                            !IF(B.lt.A) CYCLE
                            !IF((i.lt.j).and.(k.lt.l)) CYCLE
                            !IF((i.gt.j).and.(k.lt.l)) CYCLE

                            ex(1,1) = 2*i
                            ex(1,2) = 2*k
                            ex(2,1) = 2*j
                            ex(2,2) = 2*l
                            hel = GetHFInt_spinorb(ex,FullHFOrbs)
                            if(abs(hel).gt.1.0e-8_dp) then
                                WRITE(iunit,'(1X,G20.14,4I3)') hel,i,j,k,l
                            endif
                        enddo
                    enddo
                enddo
            enddo
            do i=1,nSites
                do j=1,i
                    if(abs(h0HF(i,j)).gt.1.0e-8_dp) then
                        WRITE(iunit,'(1X,G20.14,4I3)') h0HF(i,j),i,j,0,0
                    endif
                enddo
            enddo
            do i=1,nSites
                WRITE(iunit,'(1X,G20.14,4I3)') FullHFEnergies(i),i,0,0,0
            enddo
            WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
            close(iunit)
            deallocate(h0HF)
        else
            write(6,'(A)') "Writing out FCIDUMP file corresponding to non-interacting orbitals"
            allocate(temp(nSites,nSites))
            allocate(h0HF(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,HFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,HFOrbs,nSites,0.0_dp,h0HF,nSites)
            deallocate(temp)

            do i=1,nSites
                do j=1,nSites
                    A=(i*(i-1))/2+j
                    DO k=1,nSites
                        DO l=1,nSites
                            B=(k*(k-1))/2+l

                            !IF(B.lt.A) CYCLE
                            !IF((i.lt.j).and.(k.lt.l)) CYCLE
                            !IF((i.gt.j).and.(k.lt.l)) CYCLE

                            ex(1,1) = 2*i
                            ex(1,2) = 2*k
                            ex(2,1) = 2*j
                            ex(2,2) = 2*l
                            hel = GetHFInt_spinorb(ex,HFOrbs)
                            if(abs(hel).gt.1.0e-8_dp) then
                                WRITE(iunit,'(1X,G20.14,4I3)') hel,i,j,k,l
                            endif
                        enddo
                    enddo
                enddo
            enddo
            do i=1,nSites
                do j=1,i
                    if(abs(h0HF(i,j)).gt.1.0e-8_dp) then
                        WRITE(iunit,'(1X,G20.14,4I3)') h0HF(i,j),i,j,0,0
                    endif
                enddo
            enddo
            do i=1,nSites
                WRITE(iunit,'(1X,G20.14,4I3)') HFEnergies(i),i,0,0,0
            enddo
            WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
            close(iunit)
            deallocate(h0HF)

        endif

        !Now dump in the AO basis
        call append_ext_real('AO_FCIDUMP',U,filename)
        open(unit=iunit,file=filename,status='unknown')
        write(iunit,'(2A6,I3,A9,I3,A6,I2,A)') '&FCI ','NORB=',nSites,', NELEC=',NEl,', MS2=',0,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,nSites
            write(iunit,'(I1,A1)',advance='no') 1,','
        enddo
        write(iunit,*) ""
        WRITE(iunit,'(A7,I1)') 'ISYM=',1
        WRITE(iunit,'(A5)') '&END'
        if(tAnderson) then
            !Only first site correlated
            write(iunit,'(1X,G20.14,4I3)') U,1,1,1,1
        else
            do i=1,nSites
                write(iunit,'(1X,G20.14,4I3)') U,i,i,i,i
            enddo
        endif
        do i=1,nSites
            do j=1,i
!                if(tChemPot.and.(i.eq.j).and.(i.eq.1)) then
!                    write(iunit,'(1X,G20.14,4I3)') h0(i,j)-(U/2.0_dp),i,j,0,0
!                elseif(abs(h0(i,j)).gt.1.0e-8_dp) then
                if(abs(h0(i,j)).gt.1.0e-8_dp) then
                    WRITE(iunit,'(1X,G20.12,4I3)') h0(i,j),i,j,0,0
                endif
            enddo
        enddo
        WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine DumpFCIDUMP

    subroutine check_openmp
!$      use omp_lib
        implicit none
!$      integer(kind=OMP_integer_kind) :: maxthreads

        tOpenMP = .false.
!$      tOpenMP = .true.

        if(tOpenMP) then
            write(6,"(A)") " *** OpenMP compile detected *** "
!$          maxthreads = OMP_get_max_threads()
!$          max_omp_threads = maxthreads
!$          write(6,"(A,I7)") "Maxiumum number of threads to be used: ",max_omp_threads
!!$          OpenMP_wallstart = OMP_get_wtime()
!           Do not allow nested threads to start
!$          call OMP_set_nested(.false.)
        else
            write(6,"(A)") " No OpenMP optimizations"
            max_omp_threads = 1
        endif

        !Lets test it out
!$OMP PARALLEL
!$      write(6,*) "My thread number is: ",OMP_get_thread_num()
!$OMP END PARALLEL

        write(6,*) ""

    end subroutine check_openmp

End Program RealHub
