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

        !Construct the fock matrix and density
            
        !Build G(iw)

        do while(.true.)
            !Converge self-energy

            !Converge global chemical potential. Returns new chemical
            !potential and new greens function

            !FFT Greens function from iw -> tau

            !Build the self-energy in tau space

            !FFT Self-energy from tau space to iw

            !Build G(iw)
        
            !Calculate energy

            !Test for convergence of self energy
        enddo

        !Write final energy and stats

        !Write Self-energy(iw)

        !Analytically continue G(iw) -> G(w) with Pade

    end subroutine GF2_Hub

end module GF2
