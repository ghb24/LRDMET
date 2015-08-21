module GF2Driver
    use GF2, only: GF2_Hub
    use GF2Data, only: GF2_MaxIter 
    implicit none

    contains

    subroutine GF2_Driver()
        implicit none

        !Run GF2 on the full Hubbard model
        call GF2_Hub(GF2_MaxIter)

    end subroutine GF2_Driver

end module GF2Driver
