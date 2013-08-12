Program RealHub
    use const
    use timing
    use errors, only: stop_all
    use globals
    use Solvers
    use LRDriver 
    use fitting
    use mat_tools
    use input
    implicit none
    logical, parameter :: tDebug = .false.

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
        tNoHL_SE = .false.
        iReuse_SE = 0
        tAllImp_LR = .false.
        iGF_Fit = 0
        tPartialSE_Fit = .false.
        iPartialSE_Fit = 0
        iSE_Constraints = 1
        DampingExponent = huge(0.0_dp)
        tConvergeMicroSE = .false.
        iMinRes_MaxIter = 20000
        tBetaExcit = .false.
        nKCalcs = 0

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
            case("SELF-CONSISTENCY")
                tSC_LR = .true.
                if(item.lt.nitems) then
                    call readi(iGF_Fit)
                endif
            case("REUSE_SELFENERGY")
                call readi(iReuse_SE)
            case("PARTIAL_SELFENERGY_FIT")
                tPartialSE_Fit = .true.
                if(item.le.nitems) then
                    call readi(iPartialSE_Fit)
                endif
            case("SELF_ENERGY_CONSTRAINTS")
                call readi(iSE_Constraints)
            case("NO_MANYBODY_SELFENERGY")
                tNoHL_SE = .true.
            case("RESPONSE_ALLIMP")
                tAllImp_LR = .true.
            case("SELFENERGY_DAMPING")
                call readf(DampingExponent)
            case("CONVERGE_MICROITER_SE")
                tConvergeMicroSE = .true.
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
                write(6,"(A)") "SELF-CONSISTENCY"
                write(6,"(A)") "NO_MANYBODY_SELFENERGY"
                write(6,"(A)") "REUSE_SELFENERGY"
                write(6,"(A)") "PARTIAL_SELFENERGY_FIT"
                write(6,"(A)") "SELF_ENERGY_CONSTRAINTS"
                write(6,"(A)") "RESPONSE_ALLIMP"
                write(6,"(A)") "SELFENERGY_DAMPING"
                write(6,"(A)") "CONVERGE_MICROITER_SE"
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
        if(tReadSystem) then
            !Ensure we don't do fitting
            tContinueConvergence = .false.
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
        if(tSC_LR.and.tLR_ReoptGS.and.tLR_DMET) then
            call stop_all(t_r,"Reoptimizing ground state not sorted yet for self-consistent response "  &
     &          //"calculations - probably shouldn't happen")
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
        if(tNoHL_SE.and.(.not.tSC_LR)) then
            call stop_all(t_r,'Can only specify no self-energy in many-body hamiltonian if we are doing self-consistency')
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

    end subroutine name_timers

    subroutine end_calc()
        use timing, only: end_timing, print_timing_report
        implicit none
            
        call deallocate_mem()
        call end_timing()
        call print_timing_report()

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
        
            allocate(MeanFieldDM(nSites,nSites))    !DM from mean-field
            MeanFieldDM(:,:) = zero
            if(tUHF) then
                allocate(MeanFieldDM_b(nSites,nSites))
                MeanFieldDM_b(:,:) = zero
            endif

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
    
    !Construct thermal embedding basis. Core and virtual set not included initially
    subroutine ConstructThermalEmbedding()
        implicit none
        real(dp), allocatable :: ProjOverlap(:,:),ProjOverlapEVals(:),Work(:)
        real(dp), allocatable :: RotOccOrbs(:,:),ImpurityOrbs(:,:)
        real(dp) :: norm,DDOT,Overlap,ThermoPot
        integer :: lwork,info,i,j,imp,ispin,jspin,ni,nj
        character(len=*), parameter :: t_r='ConstructThermalEmbedding'
        integer :: nbath

        write(6,*) "Constructing schmidt basis of thermal zeroth order wavefunction"

        !First, calculate the thermodynamic potential
        !Assume that the chemical potential is zero
        ThermoPot = zero
        do i = 1,nSites
            ThermoPot = ThermoPot + log(one + exp(-HFEnergies(i)/Temperature))
        enddo
        ThermoPot = -ThermoPot * Temperature

        write(6,*) "Thermodynamic potential (assuming mu=0) : ",ThermoPot

        allocate(ProjOverlap(nSites,nSites))
        do i = 1,nSites
            do ispin = 1,2  !Run over spins of i
                do j = 1,nSites
                    do jspin = 1,2  !Run over spins of j    (Shouldn't actually make any difference
                        if(mod(ispin,2).ne.mod(jspin,2)) cycle
                        do ni = 0,1
                            do nj = 0,1
                                do imp = 1,nImp
                                    ProjOverlap(i,j) = ProjOverlap(i,j) + exp(-HFEnergies(i)*ni)*HFOrbs(imp,i)  &
                                        *HFOrbs(imp,j)*exp(-HFEnergies(j)*nj)
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        ProjOverlap(:,:) = ProjOverlap(:,:) * exp(ThermoPot/Temperature)

        !Diagonalize this
        allocate(ProjOverlapEVals(nSites))
        ProjOverlapEVals(:)=zero
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,ProjOverlap,nSites,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,ProjOverlap,nSites,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        write(6,*) "Largest 2 x nImp Projected overlap eigenvalues:"
        do i=nSites,nSites-2*nImp,-1
            write(6,*) ProjOverlapEVals(i)
        enddo

        if(tWriteOut) then
            call writevector(ProjOverlapEVals,'Projected overlap eigenvalues alpha')
        endif

        !We should only have nImp non-zero eigenvalues
        nbath = 0
        do i=1,nSites
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.gt.nImp) call stop_all(t_r,'error here')

        !Now rotate original occupied orbitals into this new orthogonal basis
        allocate(RotOccOrbs(nSites,nSites))
        call DGEMM('N','N',nSites,nSites,nSites,one,HFOrbs,nSites,ProjOverlap,nSites,zero,RotOccOrbs,nSites)

        !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
        !Only the last nImp orbitals will have any coupling to the impurity on them.
        !These RotOccOrbs constitute a legitamate wavefunction, are orthonormal to all other orbitals. Just simple rotation.

        !Construct bath states by projecting out component on impurity and renormalizing
        !Assume last nImp states are the states with overlap with impurity only
        !Also normalize these orbitals
        do i=nSites,nSites-nImp+1,-1
            RotOccOrbs(1:nImp,i) = zero
            norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
            norm = sqrt(norm)
            RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
        enddo

        !These states are now the bath states.
!        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

        allocate(ImpurityOrbs(nSites,nImp))
        ImpurityOrbs(:,:) = zero 
        do i=1,nImp
            ImpurityOrbs(i,i) = one 
        enddo

        !We now have all the orbitals. Which are orthogonal to which?
        write(6,*) "All alpha impurity/bath orbitals orthogonal by construction"
        do i=nOcc,nOcc-nImp+1,-1
            do j=nOcc,nOcc-nImp+1,-1
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(i.eq.j) then
                    if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not normalized set')
                    endif
                else
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not orthogonal set')
                    endif
                endif
            enddo
        enddo
        write(6,*) "All alpha bath/bath orbitals orthonormal"

        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,2*nImp))
        EmbeddedBasis(:,:) = zero
        EmbeddedBasis(:,1:nImp) = ImpurityOrbs(:,:)
        EmbeddedBasis(:,nImp+1:2*nImp) = RotOccOrbs(:,nOcc+1:nOcc+nImp)

        EmbSize = 2*nImp      !This is the total size of the embedded system with which to do the high-level calculation on 
        EmbSizeSpin = EmbSize
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(ProjOverlap,ProjOverlapEVals,RotOccOrbs,ImpurityOrbs)
    end subroutine ConstructThermalEmbedding
            
    !Construct full embedding basis, along with orthogonal core and virtual set, and check orthonormality of orbital space
    subroutine ConstructFullSchmidtBasis()
        implicit none
        real(dp), allocatable :: ProjOverlap(:,:),ProjOverlapEVals(:),Work(:)
        real(dp), allocatable :: RotOccOrbs(:,:),ImpurityOrbs(:,:),ProjOverlapVirt(:,:)
        real(dp), allocatable :: OverlapVirt(:,:),VirtSpace(:,:),temp(:,:) 
        real(dp) :: norm,DDOT,Overlap
        integer :: lwork,info,i,j,nVirt
        character(len=*), parameter :: t_r='ConstructFullSchmidtBasis'
        integer :: nbath

        write(6,*) "Constructing full schmidt basis of HF det, including core + virt spaces"

        !HFOrbs defines HF basis.
        !Construct full projected overlap basis and diagonalize to see redundancies
        !Just do this for the occupied orbitals, and check orthogonality to the core HF and virtual orbitals
        allocate(ProjOverlap(nOcc,nOcc))
        call DGEMM('T','N',nOcc,nOcc,nImp,1.0_dp,HFOrbs(1:nImp,1:nOcc),nImp,HFOrbs(1:nImp,1:nOcc),nImp,0.0_dp,ProjOverlap,nOcc)

        !Diagonalize this
        allocate(ProjOverlapEVals(nOcc))
        ProjOverlapEVals(:)=0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        if(tWriteOut) then
            call writevector(ProjOverlapEVals,'Projected overlap eigenvalues alpha')
        endif

        !We should only have nImp non-zero eigenvalues
        nbath = 0
        do i=1,nOcc
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.gt.nImp) call stop_all(t_r,'error here')

        !Now rotate original occupied orbitals into this new orthogonal basis
        allocate(RotOccOrbs(nSites,nOcc))
        call DGEMM('N','N',nSites,nOcc,nOcc,1.0_dp,HFOrbs(1:nSites,1:nOcc),nSites,ProjOverlap,nOcc,0.0_dp,RotOccOrbs,nSites)

        !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
        !Only the last nImp orbitals will have any coupling to the impurity on them.
        !These RotOccOrbs constitute a legitamate HF wavefunction, are orthonormal to all other orbitals. Just simple rotation.
!        call writematrix(RotOccOrbs,'Occupied Orbs in schmidt basis',.true.)

        !Construct bath states by projecting out component on impurity and renormalizing
        !Assume last nImp states are the states with overlap with impurity only
        !Also normalize these orbitals
        do i=nOcc,nOcc-nImp+1,-1
            RotOccOrbs(1:nImp,i) = 0.0_dp
            norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
            norm = sqrt(norm)
            RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
        enddo

        !These states are now the bath states.
!        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

        allocate(ImpurityOrbs(nSites,nImp))
        ImpurityOrbs(:,:) = 0.0_dp
        do i=1,nImp
            ImpurityOrbs(i,i) = 1.0_dp
        enddo

        !We now have all the orbitals. Which are orthogonal to which?
        write(6,*) "All alpha impurity/bath orbitals orthogonal by construction"
        do i=nOcc,nOcc-nImp+1,-1
            do j=nOcc,nOcc-nImp+1,-1
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(i.eq.j) then
                    if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not normalized set')
                    endif
                else
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not orthogonal set')
                    endif
                endif
            enddo
        enddo
        write(6,*) "All alpha bath/bath orbitals orthonormal"

        !Now, consider the rotated core orbitals
        !Are they orthonormal to the bath orbitals (they are othogonal to impurity, since they have no component on them)
        !They must also be orthogonal to the bath orbitals, since they never had any component on the impurity, 
        !which is the only bit which has changed in constructing the bath (norm doesn't affect)
        do i=1,nOcc-nImp
            do j=nOcc,nOcc-nImp+1,-1
                !Overlap of core (i) with bath (j)
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'bath orbitals with core not orthogonal set')
                endif
            enddo
        enddo

        !However, the virtual space is *not* orthogonal to the embedded system (though it is wrt the core).
        !Now create orthogonal set of orbitals from the virtual space. There will now be a redundancy.
        !Calculate the overlap of the virtual space with a projection onto the embedded system.
        !From the diagonalization of this, we expect exactly *nImp* non-zero eigenvalues, which are the redundant orbitals.
        !Remove these, and the rest are the now non-canonical virtual orbital space.
        nVirt = nSites - nOcc
        allocate(ProjOverlapVirt(nVirt,nVirt))

        !This array is used to calculate the overlap of each virtual space function with each impurity
        !function
        allocate(OverlapVirt(2*nImp,nVirt))
        !The first nImp correspond to the impurity orbital, and the next two correspond to the bath orbitals
        OverlapVirt(:,:) = 0.0_dp
        do i=nOcc+1,nSites  !run through virtual space
            do j=1,nImp     !run through impurity orbitals
                OverlapVirt(j,i-nOcc) = HFOrbs(j,i)
            enddo
        enddo

        !Now calculate overlap with bath orbitals
        call DGEMM('T','N',nImp,nVirt,nSites,1.0_dp,RotOccOrbs(:,nOcc-nImp+1:nOcc),nSites,HFOrbs(:,nOcc+1:nSites), &
            nSites,0.0_dp,OverlapVirt(nImp+1:2*nImp,1:nVirt),2*nImp-nImp)

        !Combine overlaps to get full projected overlap matrix
        call DGEMM('T','N',nVirt,nVirt,2*nImp,1.0_dp,OverlapVirt,2*nImp,OverlapVirt,2*nImp,0.0_dp,ProjOverlapVirt,nVirt)

        !Diagonalize this
        deallocate(ProjOverlapEVals)
        allocate(ProjOverlapEVals(nVirt))
        ProjOverlapEVals(:)=0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

!        call writevector(ProjOverlapEVals,'virtual overlap eigenvalues')

        nbath = 0   !Count the number of virtual functions which span the space of the embedded system
        do i=1,nVirt
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.ne.nImp) then
            call stop_all(t_r,'Virtual space redundancy not as expected')
        endif

        !Now rotate orbitals such that they are orthogonal, while deleting the redundant ones
        !Assume the last nImp are redundant
        allocate(VirtSpace(nSites,nVirt-nImp))
        call DGEMM('N','N',nSites,nVirt-nImp,nVirt,1.0_dp,HFOrbs(:,nOcc+1:nSites),nSites,   &
            ProjOverlapVirt(1:nVirt,1:nVirt-nImp),nVirt,0.0_dp,VirtSpace,nSites)

        !Check now that the virtual space is now orthogonal to all occupied HF orbitals, as well as the impurity and bath sites

        !First against the impurity sites
        do i=1,nImp
            do j=1,nVirt-nImp
                Overlap = DDOT(nSites,ImpurityOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'virtual orbitals not orthogonal to impurity orbitals')
                endif
            enddo
        enddo

        do i=1,nOcc
            do j=1,nVirt-nImp
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    if(i.gt.(nOcc-nImp)) then
                        call stop_all(t_r,'virtual orbitals not orthogonal to bath orbitals')
                    else
                        call stop_all(t_r,'virtual orbitals not orthogonal to core orbitals')
                    endif
                endif
            enddo
        enddo

        if(allocated(FullSchmidtBasis)) deallocate(FullSchmidtBasis)
        allocate(FullSchmidtBasis(nSites,nSites))   ! (Atomicbasis,SchmidtBasis)
        FullSchmidtBasis(:,:) = 0.0_dp
        FullSchmidtBasis(:,1:nOcc-nImp) = RotOccOrbs(:,1:nOcc-nImp)    !The first orbitals consist of core  
        FullSchmidtBasis(:,nOcc-nImp+1:nOcc) = ImpurityOrbs(:,:)    !Impurity Orbs
        FullSchmidtBasis(:,nOcc+1:nOcc+nImp) = RotOccOrbs(:,nOcc-nImp+1:nOcc)   !Bath
        FullSchmidtBasis(:,nOcc+nImp+1:nSites) = VirtSpace(:,:)     !Virtual space

        !call writematrix(FullSchmidtBasis,'FullSchmidtBasis',.true.)

        !Construct unitary basis transformation matrix from HF to embedded basis
        if(allocated(HFtoSchmidtTransform)) deallocate(HFtoSchmidtTransform)
        allocate(HFtoSchmidtTransform(nSites,nSites))
        call DGEMM('T','N',nSites,nSites,nSites,1.0_dp,HFOrbs,nSites,FullSchmidtBasis,nSites,0.0_dp,HFtoSchmidtTransform,nSites)

        allocate(temp(nSites,nSites))
        !Check that this operator is unitary
        !call DGEMM('N','T',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,HFtoSchmidtTransform,nSites,0.0_dp,temp,nSites)
        !call writematrix(temp,'Test of unitarity of HF to Schmidt Transform',.true.)
        !do i=1,nSites
        !    do j=1,nSites
        !        if((i.ne.j).and.(abs(temp(i,j)).gt.1.0e-7_dp)) then
        !            call stop_all(t_r,'Transformation matrix not unitary')
        !        elseif((i.eq.j).and.(abs(temp(i,j)-1.0_dp).gt.1.0e-7)) then
        !            call stop_all(t_r,'Transformation matrix not unitary')
        !        endif
        !    enddo
        !enddo

        !Use this to rotate the fock operator into this new basis
        if(allocated(FockSchmidt)) deallocate(FockSchmidt)
        allocate(FockSchmidt(nSites,nSites))
        !Set up FockSchmidt to temperarily to be the HF basis fock operator (i.e. diagonal)
        FockSchmidt(:,:) = zero
        do i=1,nSites
            FockSchmidt(i,i) = HFEnergies(i)
        enddo
        call DGEMM('T','N',nSites,nSites,nSites,one,HFtoSchmidtTransform,nSites,FockSchmidt,nSites,zero,temp,nSites)
        call DGEMM('N','N',nSites,nSites,nSites,one,temp,nSites,HFtoSchmidtTransform,nSites,zero,FockSchmidt,nSites)
            
        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,2*nImp))
        EmbeddedBasis(:,:) = 0.0_dp
        EmbeddedBasis(:,1:nImp) = ImpurityOrbs(:,:)
        EmbeddedBasis(:,nImp+1:2*nImp) = FullSchmidtBasis(:,nOcc+1:nOcc+nImp)

        if(tUHF) then
            !Do all the same for the beta orbital space
            write(6,"(A)") "Constructing beta bath space"
            !HFOrbs defines HF basis.
            !Construct full projected overlap basis and diagonalize to see redundancies
            !Just do this for the occupied orbitals, and check orthogonality to the core HF and virtual orbitals
            call DGEMM('T','N',nOcc,nOcc,nImp,one,HFOrbs_b(1:nImp,1:nOcc),nImp,HFOrbs_b(1:nImp,1:nOcc),nImp,zero,ProjOverlap,nOcc)

            !Diagonalize this
            deallocate(ProjOverlapEVals)
            allocate(ProjOverlapEVals(nOcc))
            ProjOverlapEVals(:) = zero
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            if(tWriteOut) then
                call writevector(ProjOverlapEVals,'Projected overlap eigenvalues beta')
            endif

            !We should only have nImp non-zero eigenvalues
            nbath = 0
            do i=1,nOcc
                if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                    nbath = nbath + 1
                endif
            enddo
            if(nbath.gt.nImp) call stop_all(t_r,'error here')

            !Now rotate original beta occupied orbitals into this new orthogonal basis
            call DGEMM('N','N',nSites,nOcc,nOcc,one,HFOrbs_b(1:nSites,1:nOcc),nSites,ProjOverlap,nOcc,zero,RotOccOrbs,nSites)

            !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
            !Only the last nImp orbitals will have any coupling to the impurity on them.
            !These RotOccOrbs constitute a legitamate HF wavefunction, are orthonormal to all other orbitals. Just simple rotation.
    !        call writematrix(RotOccOrbs,'Occupied Orbs in schmidt basis',.true.)

            !Construct bath states by projecting out component on impurity and renormalizing
            !Assume last nImp states are the states with overlap with impurity only
            !Also normalize these orbitals
            do i=nOcc,nOcc-nImp+1,-1
                RotOccOrbs(1:nImp,i) = zero
                norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
                norm = sqrt(norm)
                RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
            enddo

            !These states are now the bath states.
    !        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

            ImpurityOrbs(:,:) = zero 
            do i=1,nImp
                ImpurityOrbs(i,i) = one 
            enddo

            !We now have all the orbitals. Which are orthogonal to which?
            write(6,*) "All beta impurity/bath orbitals orthogonal by construction"
            do i=nOcc,nOcc-nImp+1,-1
                do j=nOcc,nOcc-nImp+1,-1
                    Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                    if(i.eq.j) then
                        if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                            call stop_all(t_r,'beta bath orbitals not normalized set')
                        endif
                    else
                        if(abs(Overlap).gt.1.0e-7_dp) then
                            call stop_all(t_r,'beta bath orbitals not orthogonal set')
                        endif
                    endif
                enddo
            enddo
            write(6,*) "All beta bath/bath orbitals orthonormal"

            !Now, consider the rotated core orbitals
            !Are they orthonormal to the bath orbitals (they are othogonal to impurity, since they have no component on them)
            !They must also be orthogonal to the bath orbitals, since they never had any component on the impurity, 
            !which is the only bit which has changed in constructing the bath (norm doesn't affect)
            do i=1,nOcc-nImp
                do j=nOcc,nOcc-nImp+1,-1
                    !Overlap of core (i) with bath (j)
                    Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'beta bath orbitals with core not orthogonal set')
                    endif
                enddo
            enddo

            !However, the virtual space is *not* orthogonal to the embedded system (though it is wrt the core).
            !Now create orthogonal set of orbitals from the virtual space. There will now be a redundancy.
            !Calculate the overlap of the virtual space with a projection onto the embedded system.
            !From the diagonalization of this, we expect exactly *nImp* non-zero eigenvalues, which are the redundant orbitals.
            !Remove these, and the rest are the now non-canonical virtual orbital space.
            nVirt = nSites - nOcc
            !allocate(ProjOverlapVirt(nVirt,nVirt))

            !This array is used to calculate the overlap of each virtual space function with each impurity
            !function
            !allocate(OverlapVirt(2*nImp,nVirt))
            !The first nImp correspond to the impurity orbital, and the next two correspond to the bath orbitals
            OverlapVirt(:,:) = zero
            do i=nOcc+1,nSites  !run through virtual space
                do j=1,nImp     !run through impurity orbitals
                    OverlapVirt(j,i-nOcc) = HFOrbs_b(j,i)
                enddo
            enddo

            !Now calculate overlap with bath orbitals
            call DGEMM('T','N',nImp,nVirt,nSites,one,RotOccOrbs(:,nOcc-nImp+1:nOcc),nSites,HFOrbs_b(:,nOcc+1:nSites), &
                nSites,zero,OverlapVirt(nImp+1:2*nImp,1:nVirt),2*nImp-nImp)

            !Combine overlaps to get full projected overlap matrix
            call DGEMM('T','N',nVirt,nVirt,2*nImp,one,OverlapVirt,2*nImp,OverlapVirt,2*nImp,zero,ProjOverlapVirt,nVirt)

            !Diagonalize this
            deallocate(ProjOverlapEVals)
            allocate(ProjOverlapEVals(nVirt))
            ProjOverlapEVals(:) = zero
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nVirt,ProjOverlapVirt,nVirt,ProjOverlapEVals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

    !        call writevector(ProjOverlapEVals,'virtual overlap eigenvalues')

            nbath = 0   !Count the number of virtual functions which span the space of the embedded system
            do i=1,nVirt
                if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                    nbath = nbath + 1
                endif
            enddo
            if(nbath.ne.nImp) then
                call stop_all(t_r,'Virtual beta space redundancy not as expected')
            endif

            !Now rotate orbitals such that they are orthogonal, while deleting the redundant ones
            !Assume the last nImp are redundant
            !allocate(VirtSpace(nSites,nVirt-nImp))
            call DGEMM('N','N',nSites,nVirt-nImp,nVirt,one,HFOrbs_b(:,nOcc+1:nSites),nSites,   &
                ProjOverlapVirt(1:nVirt,1:nVirt-nImp),nVirt,zero,VirtSpace,nSites)

            !Check now that the virtual space is now orthogonal to all occupied HF orbitals, as well as the impurity and bath sites

            !First against the impurity sites
            do i=1,nImp
                do j=1,nVirt-nImp
                    Overlap = DDOT(nSites,ImpurityOrbs(:,i),1,VirtSpace(:,j),1)
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'virtual beta orbitals not orthogonal to impurity orbitals')
                    endif
                enddo
            enddo

            do i=1,nOcc
                do j=1,nVirt-nImp
                    Overlap = DDOT(nSites,RotOccOrbs(:,i),1,VirtSpace(:,j),1)
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        if(i.gt.(nOcc-nImp)) then
                            call stop_all(t_r,'virtual beta orbitals not orthogonal to bath orbitals')
                        else
                            call stop_all(t_r,'virtual beta orbitals not orthogonal to core orbitals')
                        endif
                    endif
                enddo
            enddo

            if(allocated(FullSchmidtBasis_b)) deallocate(FullSchmidtBasis_b)
            allocate(FullSchmidtBasis_b(nSites,nSites))   ! (Atomicbasis,SchmidtBasis)
            FullSchmidtBasis_b(:,:) = zero 
            FullSchmidtBasis_b(:,1:nOcc-nImp) = RotOccOrbs(:,1:nOcc-nImp)    !The first orbitals consist of core  
            FullSchmidtBasis_b(:,nOcc-nImp+1:nOcc) = ImpurityOrbs(:,:)    !Impurity Orbs
            FullSchmidtBasis_b(:,nOcc+1:nOcc+nImp) = RotOccOrbs(:,nOcc-nImp+1:nOcc)   !Bath
            FullSchmidtBasis_b(:,nOcc+nImp+1:nSites) = VirtSpace(:,:)     !Virtual space

            !call writematrix(FullSchmidtBasis,'FullSchmidtBasis',.true.)

            !Construct unitary basis transformation matrix from HF to embedded basis
            if(allocated(HFtoSchmidtTransform_b)) deallocate(HFtoSchmidtTransform_b)
            allocate(HFtoSchmidtTransform_b(nSites,nSites))
            call DGEMM('T','N',nSites,nSites,nSites,one,HFOrbs_b,nSites,FullSchmidtBasis_b, &
                nSites,zero,HFtoSchmidtTransform_b,nSites)

            !Check that this operator is unitary
            !call DGEMM('N','T',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,HFtoSchmidtTransform,nSites,0.0_dp,temp,nSites)
            !call writematrix(temp,'Test of unitarity of HF to Schmidt Transform',.true.)
            !do i=1,nSites
            !    do j=1,nSites
            !        if((i.ne.j).and.(abs(temp(i,j)).gt.1.0e-7_dp)) then
            !            call stop_all(t_r,'Transformation matrix not unitary')
            !        elseif((i.eq.j).and.(abs(temp(i,j)-1.0_dp).gt.1.0e-7)) then
            !            call stop_all(t_r,'Transformation matrix not unitary')
            !        endif
            !    enddo
            !enddo

            !Use this to rotate the fock operator into this new basis
            if(allocated(FockSchmidt_b)) deallocate(FockSchmidt_b)
            allocate(FockSchmidt_b(nSites,nSites))
            !Set up FockSchmidt to temperarily to be the HF basis fock operator (i.e. diagonal)
            FockSchmidt_b(:,:) = zero 
            do i=1,nSites
                FockSchmidt_b(i,i) = HFEnergies_b(i)
            enddo
            call DGEMM('T','N',nSites,nSites,nSites,one,HFtoSchmidtTransform_b,nSites,FockSchmidt_b,nSites,zero,temp,nSites)
            call DGEMM('N','N',nSites,nSites,nSites,one,temp,nSites,HFtoSchmidtTransform_b,nSites,zero,FockSchmidt_b,nSites)
        
            if(allocated(EmbeddedBasis_b)) deallocate(EmbeddedBasis_b)
            allocate(EmbeddedBasis_b(nSites,2*nImp))
            EmbeddedBasis_b(:,:) = zero
            EmbeddedBasis_b(:,1:nImp) = ImpurityOrbs(:,:)
            EmbeddedBasis_b(:,nImp+1:2*nImp) = FullSchmidtBasis_b(:,nOcc+1:nOcc+nImp)
        endif
        deallocate(temp)

        !do i=1,nSites
        !    write(6,*) "FOCKSCHMIDT: ",i,FockSchmidt(i,i)
        !enddo
        !FockSchmidt has been overwritten with the fock matrix in the schmidt basis
!        call writematrix(FockSchmidt,'Fock in schmidt basis',.true.)

!       Calculate the non-interacting core energy of the DMET wavefunction
        CoreEnergy = zero
        do i = 1,nOcc-nImp
            CoreEnergy = CoreEnergy + FockSchmidt(i,i)
            if(tUHF) CoreEnergy = CoreEnergy + FockSchmidt_b(i,i)
        enddo
        if(.not.tUHF) CoreEnergy = CoreEnergy * 2.0_dp
        write(6,*) "Non-interacting core energy for DMET wavefunction is: ",CoreEnergy

        EmbSize = 2*nImp      !This is the total size of the embedded system with which to do the high-level calculation on 
        if(tUHF) then
            EmbSizeSpin = 2*EmbSize
        else
            EmbSizeSpin = EmbSize
        endif
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(ProjOverlap,ProjOverlapEVals,RotOccOrbs,ImpurityOrbs,ProjOverlapVirt,OverlapVirt,VirtSpace)
    end subroutine ConstructFullSchmidtBasis


    !Transform the one-electron integrals into the embedded basis by unitary transformation
    subroutine Transform1e()
        implicit none
        real(dp), allocatable :: CorrPotential(:,:),temp(:,:),CorrPotential_b(:,:)
!        real(dp), allocatable :: FockPotential(:,:)
        integer :: i,j
!        character(len=*), parameter :: t_r="Transform1e"
        
        if(allocated(Emb_Fock)) deallocate(Emb_Fock)
        if(allocated(Emb_h0)) deallocate(Emb_h0)                 !Core hamiltonian
        if(allocated(Emb_MF_DM)) deallocate(Emb_MF_DM)        !Mean field RDM
        if(allocated(Emb_FockPot)) deallocate(Emb_FockPot)      !The fock potential transforms h0v into the fock matrix in the AO basis
        if(allocated(Emb_CorrPot)) deallocate(Emb_CorrPot)      !The correlation potential is the block diagonal v_pot in the AO basis
        !Emb_h0v is the core hamiltonian in the embedded basis, with correlation potential (not over the impurity sites)
        if(allocated(Emb_h0v)) deallocate(Emb_h0v)

        !Transform all of these quantities into the embedded basis
        allocate(Emb_h0(EmbSize,EmbSize))
        allocate(Emb_MF_DM(EmbSize,EmbSize))
        allocate(Emb_FockPot(EmbSize,EmbSize))
        allocate(Emb_CorrPot(EmbSize,EmbSize))
        allocate(Emb_Fock(EmbSize,EmbSize))     !For hub, this is just the core + corrPot
        allocate(Emb_h0v(EmbSize,EmbSize))

        !Rotate them
        allocate(temp(EmbSize,nSites))
        !Core hamiltonian
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_h0,EmbSize)

        if(tDebug) call writematrix(Emb_h0,'Embedded Core hamil',.true.)

        !Mean field RDM
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,MeanFieldDM,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_MF_DM,EmbSize)
        
        if(tDebug) call writematrix(Emb_MF_DM,'Embedded RDM',.true.)

        !Since we don't actually store the fock potential (the diagonal repulsion/2 electron terms), construct it here in the AO basis
        !This is the bit of the fock matrix not in the core hamiltonian
!        allocate(FockPotential(nSites,nSites))
!        FockPotential(:,:) = 0.0_dp
!        do i=1,nSites
!            !Include the on-site repulsion
!            FockPotential(i,i) = FockPotential(i,i) + U * 0.5_dp * (NEl/real(nSites))
!        enddo
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,FockPotential,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_FockPot,EmbSize)
!        deallocate(FockPotential)
        Emb_FockPot(:,:)=0.0_dp     !We don't include the 2e terms in the fock matrix here, instead, leaving them to be captured by the correlation potential
!        call writematrix(Emb_FockPot,"Embedded Fock potential",.true.)

        !We also do not store the "Correlation potential", which is the potential which is added to the fock matrix to make the DMs match
        allocate(CorrPotential(nSites,nSites))
        CorrPotential(:,:) = 0.0_dp
        !It is just h0v - h0
        do i=1,nSites
            do j=1,nSites
                CorrPotential(i,j) = h0v(i,j) - h0(i,j)
            enddo
        enddo
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,CorrPotential,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_CorrPot,EmbSize)
        deallocate(CorrPotential)
        deallocate(temp)

        if(tDebug) call writematrix(Emb_CorrPot,"Embedded correlation potential",.true.)

        !The 2e terms left out of fock matrix everywhere, although included basically in CorrPot
        Emb_Fock(:,:) = Emb_h0(:,:) + Emb_CorrPot(:,:) ! + Emb_FockPot  
        
        !Set the correlation potential to zero over the impurity, since we do not need to include the fitted potential over the 
        !orbitals that we are going to do the high-level calculation on.
        Emb_h0v(:,:) = Emb_CorrPot(:,:)
        Emb_h0v(1:nImp,1:nImp) = 0.0_dp     !Set correlation potential over the impurity sites to zero
        Emb_h0v(:,:) = Emb_h0v(:,:) + Emb_h0(:,:)    !Add the embedded core hamiltonian over all sites

        if(tUHF) then
            !Now for the beta orbitals!!
            if(allocated(Emb_Fock_b)) deallocate(Emb_Fock_b)
            if(allocated(Emb_h0_b)) deallocate(Emb_h0_b)                 !Core hamiltonian
            if(allocated(Emb_MF_DM_b)) deallocate(Emb_MF_DM_b)        !Mean field RDM
            if(allocated(Emb_FockPot_b)) deallocate(Emb_FockPot_b)      !The fock potential transforms h0v into the fock matrix in the AO basis
            if(allocated(Emb_CorrPot_b)) deallocate(Emb_CorrPot_b)      !The correlation potential is the block diagonal v_pot in the AO basis
            !Emb_h0v is the core hamiltonian in the embedded basis, with correlation potential (not over the impurity sites)
            if(allocated(Emb_h0v_b)) deallocate(Emb_h0v_b)

            !Transform all of these quantities into the embedded basis
            allocate(Emb_h0_b(EmbSize,EmbSize))
            allocate(Emb_MF_DM_b(EmbSize,EmbSize))
            allocate(Emb_FockPot_b(EmbSize,EmbSize))
            allocate(Emb_CorrPot_b(EmbSize,EmbSize))
            allocate(Emb_Fock_b(EmbSize,EmbSize))     !For hub, this is just the core + corrPot
            allocate(Emb_h0v_b(EmbSize,EmbSize))

            !Rotate them
            allocate(temp(EmbSize,nSites))
            !Core hamiltonian
            call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis_b,nSites,h0_b,nSites,0.0_dp,temp,EmbSize)
            call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis_b,nSites,0.0_dp,Emb_h0_b,EmbSize)

            if(tDebug) call writematrix(Emb_h0_b,'Embedded Core beta hamil',.true.)

            !Mean field RDM
            call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis_b,nSites,MeanFieldDM_b,nSites,0.0_dp,temp,EmbSize)
            call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis_b,nSites,0.0_dp,Emb_MF_DM_b,EmbSize)
            
            if(tDebug) call writematrix(Emb_MF_DM_b,'Embedded beta RDM',.true.)

            !Since we don't actually store the fock potential (the diagonal repulsion/2 electron terms), construct it here in the AO basis
            !This is the bit of the fock matrix not in the core hamiltonian
    !        allocate(FockPotential(nSites,nSites))
    !        FockPotential(:,:) = 0.0_dp
    !        do i=1,nSites
    !            !Include the on-site repulsion
    !            FockPotential(i,i) = FockPotential(i,i) + U * 0.5_dp * (NEl/real(nSites))
    !        enddo
    !        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,FockPotential,nSites,0.0_dp,temp,EmbSize)
    !        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_FockPot,EmbSize)
    !        deallocate(FockPotential)
            Emb_FockPot_b(:,:)=0.0_dp     !We don't include the 2e terms in the fock matrix here, instead, leaving them to be captured by the correlation potential
    !        call writematrix(Emb_FockPot,"Embedded Fock potential",.true.)

            !We also do not store the "Correlation potential", which is the potential which is added to the fock matrix to make the DMs match
            allocate(CorrPotential_b(nSites,nSites))
            CorrPotential_b(:,:) = 0.0_dp
            !It is just h0v - h0
            do i=1,nSites
                do j=1,nSites
                    CorrPotential_b(i,j) = h0v_b(i,j) - h0_b(i,j)
                enddo
            enddo
            call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis_b,nSites,CorrPotential_b,nSites,0.0_dp,temp,EmbSize)
            call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis_b,nSites,0.0_dp,Emb_CorrPot_b,EmbSize)
            deallocate(CorrPotential_b)
            deallocate(temp)

            if(tDebug) call writematrix(Emb_CorrPot_b,"Embedded beta correlation potential",.true.)

            !The 2e terms left out of fock matrix everywhere, although included basically in CorrPot
            Emb_Fock_b(:,:) = Emb_h0_b(:,:) + Emb_CorrPot_b(:,:) ! + Emb_FockPot  
            
            !Set the correlation potential to zero over the impurity, since we do not need to include the fitted potential over the 
            !orbitals that we are going to do the high-level calculation on.
            Emb_h0v_b(:,:) = Emb_CorrPot_b(:,:)
            Emb_h0v_b(1:nImp,1:nImp) = 0.0_dp     !Set correlation potential over the impurity sites to zero
            Emb_h0v_b(:,:) = Emb_h0v_b(:,:) + Emb_h0_b(:,:)    !Add the embedded core hamiltonian over all sites

        endif
        
    end subroutine Transform1e

    !Calculate the embedded basis, and transform the 1e operators into this new basis
    subroutine CalcEmbedding()
        implicit none
        real(dp) :: ImpurityOverlap(nImp,nImp),OverlapEVs(nImp)
        real(dp), allocatable :: RDMonImp(:,:),Work(:),SminHalf(:,:),temp(:,:)
        integer :: info,lWork,i,j,nDelete,nSys
        character(len=*), parameter :: t_r="CalcEmbedding"
            
        !Take density matrix over impurity: nSites-nImp x nImp
        allocate(RDMonImp(nSites-nImp,nImp))
        do i=nImp+1,nSites
            do j=1,nImp
                RDMonImp(i-nImp,j) = MeanFieldDM(i,j)
            enddo
        enddo

!        call writematrix(RDMonImp,'RDMonImp',.true.)

        !Now, we want to Lowdin orthogonalize these orbitals, i.e. |psi_alpha> = |psi_beta> S_beta,alpha^(-1/2)

        !Find overlap of this density over the impurity (contract out nSites-nImp)
        call DGEMM('T','N',nImp,nImp,nSites-nImp,1.0_dp,RDMonImp,nSites-nImp,RDMonImp,nSites-nImp,0.0_dp,ImpurityOverlap,nImp)
!        call writematrix(ImpurityOverlap,'s',.true.)
        !Now find the eigenbasis of this overlap in order to raise it to the power of -1/2
        !Diagonalize the system over the impurity sites
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        if(nImp.eq.1) then
            if(abs(ImpurityOverlap(1,1)-1.0_dp).gt.1.0e-7_dp) call stop_all(t_r,'Diag error')
        endif
        !call writevector(OverlapEVs,'S_EVs')

        nSys = 0    !nSys is the number of orbitals in the bath having removed linear dependencies
        nDelete = 0
        do i=1,nImp
            if(OverlapEVs(i).lt.1.0e-10_dp) then
                write(6,*) "Warning: Bath basis linearly dependent"
                write(6,*) "Overlap eigenvalue: ",OverlapEVs(i)
                nDelete = nDelete + 1
            else
                nSys = nSys + 1
            endif
        enddo

        if(nSys.ne.nImp) write(6,*) "Bath orbitals removed due to linear dependency: ",nDelete
        write(6,*) "Total bath orbitals: ",nSys
        if(nDelete.ne.(nImp-nSys)) call stop_all(t_r,'nDelete.ne.(nImp-nSys)')

        !Now, we need to construct the full embedded basis including orbitals on the impurity sites
        !The orbitals to remove if necessary will be first.
        !Construct S^(-1/2), and then rotate back into the original basis
        allocate(SminHalf(nSys,nSys))
        SminHalf = 0.0_dp
        do i=1,nSys
            SminHalf(i,i) = OverlapEVs(i+nDelete)**(-0.5_dp)
        enddo
!        call writematrix(ImpurityOverlap(1:nImp,nDelete+1:nImp),'ImpOverlap',.true.)
        !Now rotate back into original basis as U D U^T (Having removed the first nDelete orbitals from the Impurity overlap eigenvectors
        allocate(temp(nImp,nSys))
        call DGEMM('N','N',nImp,nSys,nSys,1.0_dp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,SminHalf,nSys,0.0_dp,temp,nImp)
        call DGEMM('N','T',nImp,nImp,nSys,1.0_dp,temp,nImp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,0.0_dp,SminHalf,nImp)
        deallocate(temp)
        !SminHalf in now in the original basis
!        call writematrix(SminHalf,'SminHalf',.true.)

        !Now, we need to multiply the original bath orbitals extracted from the RDM on the impurity by this matrix to complete the Lowdin orthogonalization
        allocate(temp(nSites-nImp,nImp))
        call DGEMM('N','N',nSites-nImp,nImp,nImp,1.0_dp,RDMonImp,nSites-nImp,SminHalf,nImp,0.0_dp,temp,nSites-nImp)
        !temp now holds the orthogonalized bath orbitals
!        call writematrix(temp,'Orthog bath',.true.)

        !Now include the impurity sites along with the bath orbitals to make the final embedded basis
        !The system orbitals are represented by the first nImp orbitals, and the bath orbitals are the last nSys
        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,nImp+nSys))
        EmbeddedBasis(:,:) = 0.0_dp
        do i=1,nImp
            EmbeddedBasis(i,i) = 1.0_dp
        enddo

        if(nImp.ne.nSys) call stop_all(t_r,'You need to check in the code here to make sure we have the right bath orbitals')
        !Here, we only want to add orbitals on the bath which haven't been removed due to being linearly dependent
        !Which orbitals in temp are the removed linear dependent ones? TODO: Check this is being done correctly
        !Currently, we are just taking the first nSys of them
        EmbeddedBasis(nImp+1:nSites,nImp+1:2*nImp) = temp(:,1:nSys)

        if(tUHF) then
            !Now for beta orbitals
            !Take density matrix over impurity: nSites-nImp x nImp
            RDMonImp = zero
            do i=nImp+1,nSites
                do j=1,nImp
                    RDMonImp(i-nImp,j) = MeanFieldDM_b(i,j)
                enddo
            enddo

    !        call writematrix(RDMonImp,'RDMonImp',.true.)

            !Now, we want to Lowdin orthogonalize these orbitals, i.e. |psi_alpha> = |psi_beta> S_beta,alpha^(-1/2)

            !Find overlap of this density over the impurity (contract out nSites-nImp)
            call DGEMM('T','N',nImp,nImp,nSites-nImp,one,RDMonImp,nSites-nImp,RDMonImp,nSites-nImp,zero,ImpurityOverlap,nImp)
    !        call writematrix(ImpurityOverlap,'s',.true.)
            !Now find the eigenbasis of this overlap in order to raise it to the power of -1/2
            !Diagonalize the system over the impurity sites
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
            if(nImp.eq.1) then
                if(abs(ImpurityOverlap(1,1)-1.0_dp).gt.1.0e-7_dp) call stop_all(t_r,'Diag error')
            endif
            !call writevector(OverlapEVs,'S_EVs')

            nSys = 0    !nSys is the number of orbitals in the bath having removed linear dependencies
            nDelete = 0
            do i=1,nImp
                if(OverlapEVs(i).lt.1.0e-10_dp) then
                    write(6,*) "Warning: Bath basis linearly dependent"
                    write(6,*) "Overlap eigenvalue: ",OverlapEVs(i)
                    nDelete = nDelete + 1
                else
                    nSys = nSys + 1
                endif
            enddo

            if(nSys.ne.nImp) write(6,*) "Bath orbitals removed due to linear dependency: ",nDelete
            write(6,*) "Total beta bath orbitals: ",nSys
            if(nDelete.ne.(nImp-nSys)) call stop_all(t_r,'nDelete.ne.(nImp-nSys)')

            !Now, we need to construct the full embedded basis including orbitals on the impurity sites
            !The orbitals to remove if necessary will be first.
            !Construct S^(-1/2), and then rotate back into the original basis
            SminHalf = zero 
            do i=1,nSys
                SminHalf(i,i) = OverlapEVs(i+nDelete)**(-0.5_dp)
            enddo
    !        call writematrix(ImpurityOverlap(1:nImp,nDelete+1:nImp),'ImpOverlap',.true.)
            !Now rotate back into original basis as U D U^T (Having removed the first nDelete orbitals from the Impurity overlap eigenvectors
            deallocate(temp)
            allocate(temp(nImp,nSys))
            call DGEMM('N','N',nImp,nSys,nSys,1.0_dp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,SminHalf,nSys,0.0_dp,temp,nImp)
            call DGEMM('N','T',nImp,nImp,nSys,1.0_dp,temp,nImp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,0.0_dp,SminHalf,nImp)
            deallocate(temp)
            !SminHalf in now in the original basis
    !        call writematrix(SminHalf,'SminHalf',.true.)

            !Now, we need to multiply the original bath orbitals extracted from the RDM on the impurity by this matrix to complete the Lowdin orthogonalization
            allocate(temp(nSites-nImp,nImp))
            call DGEMM('N','N',nSites-nImp,nImp,nImp,1.0_dp,RDMonImp,nSites-nImp,SminHalf,nImp,0.0_dp,temp,nSites-nImp)
            !temp now holds the orthogonalized bath orbitals
    !        call writematrix(temp,'Orthog bath',.true.)

            !Now include the impurity sites along with the bath orbitals to make the final embedded basis
            !The system orbitals are represented by the first nImp orbitals, and the bath orbitals are the last nSys
            if(allocated(EmbeddedBasis_b)) deallocate(EmbeddedBasis_b)
            allocate(EmbeddedBasis_b(nSites,nImp+nSys))
            EmbeddedBasis_b(:,:) = zero 
            do i=1,nImp
                EmbeddedBasis_b(i,i) = one 
            enddo

            if(nImp.ne.nSys) call stop_all(t_r,'You need to check in the code here to make sure we have the right bath orbitals')
            !Here, we only want to add orbitals on the bath which haven't been removed due to being linearly dependent
            !Which orbitals in temp are the removed linear dependent ones? TODO: Check this is being done correctly
            !Currently, we are just taking the first nSys of them
            EmbeddedBasis_b(nImp+1:nSites,nImp+1:2*nImp) = temp(:,1:nSys)

        endif
        
        EmbSize = 2*nImp      !This is the total size of the embedded system with which to do the high-level calculation on 
        if(tUHF) then
            EmbSizeSpin = 2*EmbSize
        else
            EmbSizeSpin = EmbSize
        endif
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(temp,SminHalf)
        
    end subroutine CalcEmbedding

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
                if(tChemPot.and.(i.eq.j).and.(i.eq.1)) then
                    write(iunit,'(1X,G20.14,4I3)') h0(i,j)-(U/2.0_dp),i,j,0,0
                elseif(abs(h0(i,j)).gt.1.0e-8_dp) then
                    WRITE(iunit,'(1X,G20.12,4I3)') h0(i,j),i,j,0,0
                endif
            enddo
        enddo
        WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine DumpFCIDUMP

End Program RealHub
