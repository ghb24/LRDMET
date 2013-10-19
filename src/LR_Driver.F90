module LRDriver
    use const
    use errors, only: stop_all
    use SingRefLR
    use LinearResponse
    use MomSpectra 
    use globals
    implicit none

    contains

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
            call Converge_SE(SE,nESteps)
            SE_Old(:,:,:) = SE(:,:,:)
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
