module GF2
    use const
    use timing
    use errors
    use globals
    use utils, only: get_free_unit,append_ext_real,append_ext
    use GF2Data
    use matrixops, only: mat_inv
    use writedata
    implicit none

    contains

    !Input: Max iterations, Convergence threshold
    !       Initial Self-energy
    !       Whether Initial self-energy is fixed
    !       Damping of self-energy
    !       Initial Chemical potential
    !       Damping of density matrix
    subroutine GF2_Hub()
        implicit none
        integer, intent(in) :: Max


        !Construct the fock matrix and density matrix
        !Returned in globals "FockMat" and "DensityMat"
        !Mean-field chemical potential returned in GuessChemPot
        call GetFockandP(GuessChemPot)
        
        if(present(InitChemPot)) then
            LatChemPot = InitChemPot
        else
            LatChemPot = GuessChemPot
        endif
        write(6,*) "Chemical potential initially set to: ",
        
        !Set up tau grids, and matsubara grids

        !Set initial self-energy if necessary

        !Initial memory allocation

        

        !Build G(iw) (Should be done in kspace)
        call BuildMatsubaraGF()

        do while(.true.)
            !Converge self-energy loop

            !Converge global chemical potential. Returns new chemical
            !potential and new greens function
            call ConvergeChemicalPotential()

            !FFT Greens function from iw -> tau

            !Build the self-energy in tau space

            !FFT Self-energy from tau space to iw

            !Build G(iw)
        
            !Calculate energy

            !Test for convergence of self energy, or if hit iteration limit
        enddo

        !Write final energy and stats

        !Write Self-energy(iw)

        !Analytically continue G(iw) -> G(w) with Pade

    end subroutine GF2_Hub

    subroutine ConvergeChemicalPotential()
        implicit none

        do while(.true.)

            !Find Density from matsubara greens function
            call GetGSDensityFromMatsuGF()

            !Find Number of electrons from density

            !Test if global number of electrons correct

            !Build fock matrix from h0 and diagonal of density matrix

            !Update G(iw) [TODO: Do we really need to remake the whole matrix?]

        enddo

    end subroutine ConvergeChemicalPotential

    !P = -2 x G(tau = Beta)
    subroutine GetGSDensityFromMatsuGF()
        implicit none

        if(Temperature.lt.1e-6_dp) then
            !Zero temperature
            Beta = ImTimePoints(nImTimePoints)
        else
            Beta = one/Temperature
        endif
        DensityMat(:,:) = zero
        do i = 1,nMatsuFreqs
            DensityMat(:,:) = DensityMat(:,:) + MatsuGF(:,:,i)*exp(cmplx(zero,-MatsuFreqs(i)*Beta))
        enddo
        DensityMat(:,:) = -(2.0_dp/Beta)*DensityMat(:,:)

    end subroutine GetGSDensityFromMatsuGF

    !Get mean field density matrix and fock matrix
    subroutine GetFockandP(GuessChemPot)
        implicit none
        real(dp), intent(out) :: GuessChemPot
        real(dp), allocatable :: EigenVals(:)
        integer :: i

        if(.not.allocated(FockMat)) then
            allocate(FockMat(nSites,nSites))
        endif
        FockMat(:,:) = h0(:,:)
        do i = 1,nSites
            FockMat(i,i) = FockMat(i,i) + U*real(NEl,dp)/(2.0_dp*nSites)
        enddo
        allocate(EigenVecs(nSites,nSites))
        allocate(EigenVals(nSites))
        EigenVecs(:,:) = FockMat(:,:)
        call DiagOneEOp(EigenVecs,EigenVals,1,nSites,.false.,.true.)

        !Calculate Density
        if(.not.allocated(DensityMat))
            allocate(DensityMat(nSites,nSites))
        endif
        allocate(OccOrbs(nSites,nOcc))
        OccOrbs(:,:) = EigenVecs(:,1:nOcc)
        call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,OccOrbs,nSites,OccOrbs,nSites,zero,DensityMat,nSites)

        GuessChemPot = (EigenVals(nOcc) + EigenVals(nOcc+1))/2.0_dp
        deallocate(OccOrbs,EigenVals,EigenVecs)

    end subroutine GetFockandP

    !TODO: Build this in k-space
    subroutine BuildMatsubaraGF()
        implicit none
        complex(dp), allocatable :: SingleFreqMat(:,:)
        integer :: i,j

        MatsuGF(:,:,:) = zzero
!$OMP PARALLEL DO PRIVATE(j,SingleFreqMat)
        do i = 1,nMatsuFreqs
            allocate(SingleFreqMat(nSites,nSites))
            !Build matrix
            SingleFreqMat(:,:) = - FockMat(:,:) - LatSelfEnergy(:,:,i)
            do j = 1,nSites
                SingleFreqMat(j,j) = SingleFreqMat(j,j) + (LatChemPot + MatsuFreqs(i))
            enddo
            !Invert
            call mat_inv(SingleFreqMat,MatsuGF(:,:,i))
            deallocate(SingleFreqMat)
        enddo
!$OMP END PARALLEL DO
    end subroutine BuildMatsubaraGF

end module GF2
