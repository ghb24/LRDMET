module GF2Driver
    implicit none

    contains

    subroutine GF2_Driver()
        implicit none

        !Run GF2 on the full Hubbard model
        call GF2_Hub()

    end subroutine GF2_Driver

end module GF2Driver
