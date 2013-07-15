module const

    !Module containing parameters relating to constants in the code.

    implicit none


    ! selected_int_kind(0): equivalent to a byte.  -128 <= int_i0 <= 127.
    ! selected_int_kind(3): equivalent to a 16 bit integer.  -32768 <= int_i0 <= 32767.
    ! selected_int_kind(6): equivalent to a 32 bit integer.  -2147483648 <= int_i0 <= 2147483647.
    ! selected_int_kind(10): equivalent to a 64 bit integer. -9223372036854775808 <= int_i0 <= 9223372036854775807.

    integer, parameter :: i1 = selected_int_kind(0)     !A single byte integer
    integer, parameter :: i2 = selected_int_kind(3)     !int*2
    integer, parameter :: i4 = selected_int_kind(6)     !int*4
    integer, parameter :: i8 = selected_int_kind(10)    !int*8
    integer, parameter :: ip = i8

    integer, parameter :: sp = selected_real_kind(6,37)     !For single precision real numbers
    integer, parameter :: dp = selected_real_kind(15,307)   !For double precision real numbers

    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
    complex(dp), parameter :: zzero = cmplx(0.0_dp,0.0_dp,dp)
    complex(dp), parameter :: zone = cmplx(1.0_dp,0.0_dp,dp)

    real(dp), parameter :: eps = epsilon(zero)
    real(dp), parameter :: realmin = tiny(one), realmax = huge(one)

    real(dp), parameter :: pi = 3.1415926535897931_dp

    real(dp), parameter :: RealtoMb = 8.0_dp/1048576.0_dp
    real(dp), parameter :: ComptoMb = RealtoMb*2.0_dp

    ! depsilon is the precision used to compare floating point numbers.
    real(dp), parameter :: depsilon = 1.e-8

    integer :: qwerty
    integer, parameter :: n_int = kind(qwerty)
    integer, parameter :: bits_n_int = bit_size(qwerty)
    integer, parameter :: end_n_int = bits_n_int - 1
    integer, parameter :: size_n_int = bits_n_int/8

    logical, public :: debug = .false.  !This is just for the MINRESSLQ code

end module const
