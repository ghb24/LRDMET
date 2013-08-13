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
    subroutine SC_SRMom_LR()
        implicit none

        !First, calculate G_00
        call NonIntExCont_TDA_MCLR_Charged_Cmprs()

        !Now calculate the hybridization and self-energy self-consistently
        !This will read back in the greens function
        call SC_Mom_LR()

    end subroutine SC_SRMom_LR

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
