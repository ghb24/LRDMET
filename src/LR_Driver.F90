module LRDriver
    use const
    use errors, only: stop_all
    use SingRefLR
    use LinearResponse
    use MomSpectra 
    use globals
    implicit none

    contains

    subroutine SC_Mom_LR_NoDyson()
        use utils, only: get_free_unit,append_ext_real,append_ext
        implicit none
        real(dp) :: Omega,GFChemPot,OmegaVal
        integer :: nESteps,iter,i,k,l,iunit,j
        complex(dp), allocatable :: SE(:,:,:),G_Mat(:,:,:),SE_Old(:,:,:),LocalMomGF(:,:,:)
        real(dp) :: dChangeSE_Tol,dFuncTol,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF
        character(64) :: filename,filename2
        logical, allocatable :: tOmegaConv(:)
        logical :: tSuccess
        character(len=*), parameter :: t_r='SC_Mom_LR_NoDyson'

        dChangeSE_Tol = 1.0e-9_dp
        dFuncTol = 1.0e-8_dp
        
        !How many frequency points are there exactly?
        Omega = Start_Omega
        nESteps = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            nESteps = nESteps + 1
            Omega = Omega + Omega_Step
        enddo
        
        !Find the k-independent self-consistent self-energy
        allocate(SE(nImp,nImp,nESteps))
        allocate(SE_Old(nImp,nImp,nESteps))
        allocate(G_Mat(nImp,nImp,nESteps))
        allocate(LocalMomGF(nImp,nImp,nESteps))

        !Array of all frequency points, and whether they have correctly converged or not
        allocate(tOmegaConv(nESteps))

        call Init_SE(nESteps,SE)

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp
        iter = 0

        do while(.true.)

            iter = iter + 1

            write(6,*) "Starting iteration: ",iter

            !First, do high-level calculation, with SE for bath construction
            call SchmidtGF_wSE(G_Mat,GFChemPot,SE,nESteps)
            write(6,*) "1"
            call flush(6)
            call writedynamicfunction(nESteps,G_Mat,'G_Imp',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)
            write(6,"(A)") "Impurity Greens function causal"

            if(tSE_Scan) then
                call PlotSigmaSurface(nESteps,G_Mat,GFChemPot)
                call stop_all(t_r,'Finished Scanning Self-energy')
            endif

            !Now, numerically calculate the self energy, st. the lattice and impurity greens function match
            SE_Old(:,:,:) = SE(:,:,:)
            call ConvergeGlobalCoupling(nESteps,G_Mat,SE,GFChemPot,dFuncTol,dChangeSE_Tol,tOmegaConv,tSuccess,150)
            write(6,*) "2"
            call flush(6)
            !Now, simply set the values of the frequency outside the causal range to simply be equal to the largest frequency causal self-energy.
            call fit_noncausal_SE(nESteps,SE,tOmegaConv,tSuccess)
            !if(.not.tSuccess) then
                !Turn to a simplex algorithm to converge the points that failed with newton-raphson.
                !call ConvergeGlobalCoupling_Direct(nESteps,G_Mat,SE,GFChemPot,dFuncTol,dChangeSE_Tol,tOmegaConv,tSuccess)
            !endif

            !Now check whether we are actually converged, i.e. does the lattice greens function match the impurity greens function
            !Construct FT of k-space GFs and take zeroth part.
            call FindLocalMomGF(nESteps,SE,LocalMomGF)
            !call FindRealSpaceLocalMomGF(nESteps,SE,LocalMomGF)
            write(6,*) "3"
            call flush(6)
            call writedynamicfunction(nESteps,LocalMomGF,'LocalMomGF',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)
            write(6,"(A)") "Lattice greens function causal"
            !call CheckGFsSame(nESteps,G_Mat,LocalMomGF,1.0e-8_dp)
            !write(6,"(A)") "Lattice GF = Impurity GF"

            call writedynamicfunction(nESteps,SE,'SelfEnergy',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)
            write(6,"(A)") "Self energy causal"

            !Find the change in self-energy (and greens functions)
            call FindSEDiffs(SE,SE_Old,G_Mat,LocalMomGF,nESteps,MaxDiffSE,MinDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF)

            write(6,"(I7,4G15.7)") iter,MaxDiffSE,MeanDiffSE,MaxDiffGF,MeanDiffGF

            if(MaxDiffSE.lt.dChangeSE_Tol) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",dChangeSE_Tol
                exit
            endif

            !if(iter.eq.1) call stop_all(t_r,'SELF ENERGY CAUSAL (in first iteration)! YAY!')
        enddo

        write(6,"(A)") "Writing out converged self-energy"
        iunit = get_free_unit()
        call append_ext_real('Converged_SE',U,filename)
        if(.not.tHalfFill) then
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            !Write out *lower* triangle
            write(iunit,"(G25.10)",advance='no') Omega
            do k = 1,nImp
                do l = k,nImp
                    if((k.eq.nImp).and.(l.eq.nImp)) then
                        exit
                    endif
                    write(iunit,"(2G25.10)",advance='no') real(SE(l,k,i),dp),aimag(SE(l,k,i))
                enddo
            enddo
            write(iunit,"(2G25.10)") real(SE(nImp,nImp,i),dp),aimag(SE(nImp,nImp,i))
            Omega = Omega + Omega_Step
        enddo
        close(iunit)

        deallocate(SE,SE_Old,G_Mat,LocalMomGF,tOmegaConv)

    end subroutine SC_Mom_LR_NoDyson
    
    !Attempt to get k-space spectral functions by self-consistenly calculating k-independent hybridization and self-energy contributions
    !As opposed to the routine below, this one attempts to use a global hybridization in the construction of the bath, rather than the 
    !self-energy
    subroutine SC_Mom_LR_Z()
        use utils, only: get_free_unit,append_ext_real,append_ext
        implicit none
        real(dp) :: Omega,GFChemPot,OmegaVal,MaxDiffSE,DiffSE,reSE,imSE,r(3)
        complex(dp) :: DiffMatSE(nImp,nImp)
        integer :: nESteps,iter,i,k,l,iunit,j
        complex(dp), allocatable :: SE(:,:,:),G_Mat(:,:,:),SE_Old(:,:,:)
        complex(dp), allocatable :: Hybrid(:,:,:),InvLocalMomGF(:,:,:),LocalMomGF(:,:,:)
        complex(dp), allocatable :: LocalCoupFn(:,:,:),GlobalCoup(:,:,:)
        logical, allocatable :: tOmegaConv(:)
        logical :: tSuccess
        character(64) :: filename,filename2
        logical :: exists
        character(len=*), parameter :: t_r='SC_Mom_LR_Z'
        
        !How many frequency points are there exactly?
        Omega = Start_Omega
        nESteps = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            nESteps = nESteps + 1
            Omega = Omega + Omega_Step
        enddo
        
        !Find the k-independent self-consistent self-energy
        allocate(SE(nImp,nImp,nESteps))
        allocate(SE_Old(nImp,nImp,nESteps))
        allocate(G_Mat(nImp,nImp,nESteps))
        allocate(Hybrid(nImp,nImp,nESteps))
        allocate(InvLocalMomGF(nImp,nImp,nESteps))
        allocate(LocalMomGF(nImp,nImp,nESteps))
        allocate(LocalCoupFn(nImp,nImp,nESteps))    !X'
        allocate(GlobalCoup(nImp,nImp,nESteps))     !Z

        Hybrid(:,:,:) = zzero
        InvLocalMomGF(:,:,:) = zzero
        LocalMomGF(:,:,:) = zzero
        LocalCoupFn(:,:,:) = zzero
        GlobalCoup(:,:,:) = zzero

        allocate(tOmegaConv(nESteps))
        
        call Init_SE(nESteps,SE)

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp

        iter = 0

        do while(.true.)

            iter = iter + 1

            write(6,*) "Starting iteration: ",iter
            !Construct FT of k-space GFs and take zeroth part.
            !write(6,*) "FT'ing mean-field greens functions for local part"
            call FindLocalMomGF(nESteps,SE,LocalMomGF)
            !call FindRealSpaceLocalMomGF(nESteps,SE,LocalMomGF)
            write(6,*) "1"
            call flush(6)
            call writedynamicfunction(nESteps,LocalMomGF,'LocalMomGF',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)

            !Invert the matrix of non-interacting local greens functions.
            !write(6,*) "Inverting Local greens function"
            InvLocalMomGF(:,:,:) = LocalMomGF(:,:,:)
            write(6,*) "2"
            call flush(6)
            call InvertLocalNonHermGF(nESteps,InvLocalMomGF)

            !Now find hybridization
            !This is given by (omega + mu + idelta - e_0 - SE - InvGF)
            !write(6,*) "Calculating Hybridization function to local greens function"
            write(6,*) "3"
            call flush(6)
            call FindHybrid(nESteps,InvLocalMomGF,SE,Hybrid)
            call writedynamicfunction(nESteps,Hybrid,'Hybrid',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)

            !Check that G^0 = G^0' now.
            write(6,*) "4"
            call flush(6)
            call CheckNIGFsSame(nESteps,LocalMomGF,SE,Hybrid,GFChemPot)

            !Now calculate X', the local coupling function (LocalCoupFn), as
            ![omega + mu + i delta - h00 - Delta]^-1
            write(6,*) "5"
            call flush(6)
            call CalcLocalCoupling(nESteps,Hybrid,LocalCoupFn,GFChemPot)
            call writedynamicfunction(nESteps,LocalCoupFn,'LocalCouplingFn',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.)

            !Iteratively converge the global, k-independent coupling quantity 'GlobalCoup' (Z), 
            !which mimics the effect of the hybridization on the whole lattice.
            write(6,*) "6"
            call flush(6)
            tOmegaConv(:) = .false.
            call ConvergeGlobalCoupling(nESteps,LocalCoupFn,GlobalCoup,GFChemPot,1.0e-9_dp,1.0e-10_dp,tOmegaConv,tSuccess,150)

            !Check here that the two coupling functions are identical
            write(6,*) "7"
            call flush(6)
            !call CheckNICoupFnsSame(nESteps,LocalCoupFn,GlobalCoup,GFChemPot)
            call writedynamicfunction(nESteps,GlobalCoup,'GlobalCoupling_Z',tag=iter,   &
                tCheckCausal=.true.,tCheckOffDiagHerm=.true.,tWarn=.false.)

            !Construct interacting greens function from the global coupling
            !calculate G_00
            write(6,"(A)") "Calculating high-level correlation function..."
            write(6,*) "8"
            call flush(6)
            !TEST!
            !GlobalCoup(:,:,:) = zzero
            call SchmidtGF_wSE(G_Mat,GFChemPot,GlobalCoup,nESteps)
            write(6,"(A)") "High-level correlation function obtained."

            !Now calculate the change in the self energy from a damped application
            !of the dyson equation. This is not iterative, and the non-interating
            !and iteractive greens function do not necessarily agree from the outset.
            !This takes the old SE, and outputs the new one.
            SE_Old(:,:,:) = SE(:,:,:)
            call Calc_SE(nESteps,Hybrid,G_Mat,GFChemPot,SE,MaxDiffSE)
            call writedynamicfunction(nESteps,SE,'SelfEnergy',tag=iter,tCheckCausal=.true.,tCheckOffDiagHerm=.true.,tWarn=.false.)

            !Finally, should we do this all in a larger self-consistency, 
            !such that the self energy is used for the frequency dependent bath?
            if(MaxDiffSE.lt.1.0e-4_dp) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",1.0e-4_dp
                exit
            endif

            !if(iter.eq.1) call stop_all(t_r,'SELF ENERGY CAUSAL (in first iteration)! YAY!')

        enddo

        write(6,"(A)") "Writing out converged self-energy"
        iunit = get_free_unit()
        call append_ext_real('Converged_SE',U,filename)
        if(.not.tHalfFill) then
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            !Write out *lower* triangle
            write(iunit,"(G25.10)",advance='no') Omega
            do k = 1,nImp
                do l = k,nImp
                    if((k.eq.nImp).and.(l.eq.nImp)) then
                        exit
                    endif
                    write(iunit,"(2G25.10)",advance='no') real(SE(l,k,i),dp),aimag(SE(l,k,i))
                enddo
            enddo
            write(iunit,"(2G25.10)") real(SE(nImp,nImp,i),dp),aimag(SE(nImp,nImp,i))
            Omega = Omega + Omega_Step
        enddo
        close(iunit)

        deallocate(LocalCoupFn,SE,SE_Old,G_Mat,Hybrid,InvLocalMomGF,LocalMomGF,GlobalCoup,tOmegaConv)

    end subroutine SC_Mom_LR_Z

    !This is the high level routine to work out how we want to do the linear response
    !These functions are for single-reference (generally non-interacting) spectral functions, but using the correlated one-electron potential in their calculation.
    !Therefore, they should be pretty good around the ground state.
    subroutine Correlated_SR_LR()
        implicit none
        !character(len=*), parameter :: t_r='Correlated_SR_LR'

        if(tCorrNI_LocGF) then
            call CorrNI_LocalGF()
        endif
        if(tCorrNI_LocDD) then
            call CorrNI_LocalDD()
        endif
        if(tCorrNI_MomGF) then
            call CorrNI_MomGF()
        endif

    end subroutine Correlated_SR_LR


    !Attempt to get k-space spectral functions by self-consistenly calculating k-independent hybridization and self-energy contributions
    subroutine SC_Mom_LR()
        use utils, only: get_free_unit,append_ext_real,append_ext
        implicit none
        real(dp) :: Omega,GFChemPot,OmegaVal,MaxDiffSE,DiffSE,reSE,imSE,r(3)
        complex(dp) :: DiffMatSE(nImp,nImp)
        integer :: nESteps,iter,i,k,l,iunit,j
        complex(dp), allocatable :: SE(:,:,:),G_Mat(:,:,:),SE_Old(:,:,:)
        character(64) :: filename,filename2
        logical :: exists
        character(len=*), parameter :: t_r='SC_Mom_LR'
        
        !How many frequency points are there exactly?
        Omega = Start_Omega
        nESteps = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            nESteps = nESteps + 1
            Omega = Omega + Omega_Step
        enddo
        
        !Find the k-independent self-consistent self-energy
        allocate(SE(nImp,nImp,nESteps))
        allocate(SE_Old(nImp,nImp,nESteps))
        allocate(G_Mat(nImp,nImp,nESteps))

        if(tRead_SelfEnergy) then
            write(6,"(A)") "Reading in self-energy..."
            call append_ext_real('Converged_SE',U,filename)
            if(.not.tHalfFill) then
                call append_ext(filename,nOcc,filename2)
            else
                filename2 = filename
            endif
            inquire(file=filename2,exist=exists)
            if(.not.exists) then
                write(6,*) "Converged self-energy does not exist, filename: ",filename2
                call stop_all(t_r,'Cannot find file with converged selfenergy')
            endif
            
            iunit = get_free_unit()
            open(unit=iunit,file=filename2,status='unknown')
            Omega = Start_Omega
            i = 0
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
                i = i + 1
                !Read in *lower* triangle
                read(iunit,"(G25.10)",advance='no') OmegaVal
                if(abs(OmegaVal-Omega).gt.1.0e-8_dp) then
                    write(6,*) OmegaVal,Omega
                    call stop_all(t_r,'Omega values do not match up for read-in self-energy')
                endif
                do k = 1,nImp
                    do l = k,nImp
                        if((k.eq.nImp).and.(l.eq.nImp)) then
                            exit
                        endif
                        read(iunit,"(2G25.10)",advance='no') reSE,imSE
                        SE(l,k,i) = dcmplx(reSE,imSE)
                    enddo
                enddo
                read(iunit,"(2G25.10)") reSE,imSE
                SE(nImp,nImp,i) = dcmplx(reSE,imSE)
                Omega = Omega + Omega_Step
            enddo
            close(iunit)
        else
            SE(:,:,:) = zzero
        endif

        !Set chemical potential (This will only be right for half-filling)
        GFChemPot = U/2.0_dp

        iter = 0

        do while(.true.)

            iter = iter + 1

            !calculate G_00
!            call NonIntExCont_TDA_MCLR_Charged_Cmprs()
!            call NonIntExCont_TDA_MCLR_Charged()
            write(6,"(A)") "Calculating high-level correlation function..."
            call SchmidtGF_wSE(G_Mat,GFChemPot,SE,nESteps)
            write(6,"(A)") "High-level correlation function obtained."
            write(6,"(A,I8,A)") "Now attempting self-consistent determination of self-energy function from ", &
                nESteps," frequency points."
            call flush(6)

            if((iter.eq.1).and.tRandom_Init_SE) then
                !Randomize the initial SE for the microiterations on the first iteration
                CALL random_seed()
                do i = 1,nESteps
                    do j = 1,nImp
                        do k = 1,nImp
                            CALL Random_Number(r(:))
                            if(r(1).gt.0.5_dp) then
                                SE(k,j,i) = dcmplx(r(2),-r(3))
                            else
                                SE(k,j,i) = dcmplx(-r(2),-r(3))
                            endif
                        enddo
                    enddo
                enddo
            endif
            !tester
            !call testNIGFs(SE,nESteps)

            !Now calculate the (hybridization and) self-energy self-consistently
            !This will read back in the greens function
            !The returned self-energy is k-independent, but will reproduce the correlated local greens function
            SE_Old(:,:,:) = SE(:,:,:)
            call Converge_SE(SE,nESteps)
            !call Converge_SE_NoHybrid(G_Mat,SE,nESteps,iter)

            !What is the overall change in self-energy
            MaxDiffSE = zero
            do i = 1,nESteps
                DiffMatSE(:,:) = SE(:,:,i) - SE_Old(:,:,i)
                DiffSE = zero
                do j = 1,nImp
                    do k = 1,nImp
                        DiffSE = DiffSE + real(dconjg(DiffMatSE(k,j))*DiffMatSE(k,j),dp)
                    enddo
                enddo
                if(DiffSE.gt.MaxDiffSE) MaxDiffSE = DiffSE
            enddo

            !Finally, should we do this all in a larger self-consistency, 
            !such that the self energy is used for the frequency dependent bath?
            if(MaxDiffSE.lt.1.0e-4_dp) then
                write(6,"(A,G15.8)") "Self-energy macroiteration converged to: ",1.0e-4_dp
                exit
            endif

        enddo

        write(6,"(A)") "Writing out converged self-energy"
        iunit = get_free_unit()
        call append_ext_real('Converged_SE',U,filename)
        if(.not.tHalfFill) then
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        Omega = Start_Omega
        i = 0
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            i = i + 1
            !Write out *lower* triangle
            write(iunit,"(G25.10)",advance='no') Omega
            do k = 1,nImp
                do l = k,nImp
                    if((k.eq.nImp).and.(l.eq.nImp)) then
                        exit
                    endif
                    write(iunit,"(2G25.10)",advance='no') real(SE(l,k,i),dp),aimag(SE(l,k,i))
                enddo
            enddo
            write(iunit,"(2G25.10)") real(SE(nImp,nImp,i),dp),aimag(SE(nImp,nImp,i))
            Omega = Omega + Omega_Step
        enddo
        close(iunit)

    end subroutine SC_Mom_LR

    !TODO:
    !   Code up all of these options!
    !   Really work out difference between non-interacting LR, TDA and RPA, and look at how quality changes in response functions as U increased
    !   Look at difference in quality between TDA-type and RPA-type MCLR methods 
    !   Look at difference in quality between full MCLR and the fully contracted type
    !   Work out self-consistency condition to optimise both the full MCLR and the fully contracted type - response of correlation potential?
    !   Consider using different methods to obtain contraction coefficients in the fully contracted methods - TDA/RPA. Does this improve things?
    !   In the fully contracted case, should we split between impurity and environment, rather than embedded system and core? More semi-internal yuckiness...
    !   Perhaps look at CC2 
    subroutine MR_LinearResponse()
        implicit none
        
        if(tIC_TDA_Response) then
            !Create contracted single excitation space using the non-interacting reference for the contractions
            !The matrix is then created in a CI fashion
            !The difference between this and the one below is that it does not include any operators in the active space, and therefore relies on coupling between
            !the N and N+1 and N-1 active spaces.
            call NonIntContracted_TDA_MCLR()
        endif

        if(tEC_TDA_Response) then
            !Externally contracted
            if(tDDResponse) then
                if(tCompressedMats) then
                    call NonIntExCont_TDA_MCLR_DD_Cmprs()
                else
                    call NonIntExContracted_TDA_MCLR()
                endif
            endif
            if(tChargedResponse) then
                if(tCompressedMats) then
                    call NonIntExCont_TDA_MCLR_Charged_Cmprs()
                else
                    call NonIntExCont_TDA_MCLR_Charged()
                endif
            endif
            if(tCharged_MomResponse) then
                call MomGF_Ex()
            endif
        endif

        !Create contracted single excitation space using the non-interacting reference for the contractions
        !The matrix is then created in an RPA fashion
        !call NonIntContracted_RPA_MCLR()

        !Full MCLR, creating excitations in a CI fashion, rather than with commutators. Should reduce to TDA in single reference limit
!        call TDA_MCLR()

        !Full MCLR, with excitations and deexcitations. Should reduce to RPA in single reference limit
!        call RPA_MCLR()

    end subroutine MR_LinearResponse

    !Run single reference linear response calculations, based on true HF calculation.
    subroutine SR_LinearResponse()
        implicit none
        character(len=*), parameter :: t_r='SR_LinearResponse'
        
        if(tNIResponse) then
            !Non-interacting linear response
            if(tChargedResponse) then
                call warning(t_r,'NonInteracting LR not yet set up for charged perturbations. Skipping...')
            else
                call set_timer(LR_SR_NonInt) 
                call NonInteractingLR()
                call halt_timer(LR_SR_NonInt)
            endif
        endif
        if(tTDAResponse) then
            !Single reference TDA
            if(tChargedResponse) then
                call warning(t_r,'TDA LR not yet set up for charged perturbations. Skipping...')
            else
                call set_timer(LR_SR_TDA) 
                call TDA_LR()
                call halt_timer(LR_SR_TDA)
            endif
        endif
        if(tRPAResponse) then
            !Single reference RPA
            if(tChargedResponse) then
                call warning(t_r,'RPA LR not yet set up for charged perturbations. Skipping...')
            else
                call set_timer(LR_SR_RPA) 
                call RPA_LR()
                call halt_timer(LR_SR_RPA)
            endif
        endif

    end subroutine SR_LinearResponse
    
end module LRDriver
