module T1SCF
    use const
    use errors, only: stop_all
    use globals
    implicit none

    contains

    !Rotate the previous iterations Slater determinant into the new Slater determinant
    !over the entire lattice..
    !This is achieved by rotating the imp blocks of the T1 matrix into the site basis, 
    !before applying the rotation to all sites.
    subroutine RotDet_T1(iter)
        implicit none
        integer, intent(in) :: iter
        character(len=*), parameter :: t_r='RotDet_T1'

        call stop_all(t_r,'Not yet implemented')


    end subroutine RotDet_T1

end module T1SCF
