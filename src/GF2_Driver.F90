module GF2Driver
    use GF2, only: GF2_Hub
    implicit none

    contains

    subroutine GF2_Driver()
        implicit none

        !Run GF2 on the full Hubbard model
        call GF2_Hub(100)

    end subroutine GF2_Driver

end module GF2Driver
