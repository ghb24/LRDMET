module DetBitOps
    use const, only: n_int,bits_n_int,end_n_int,size_n_int
    implicit none

    contains

    pure subroutine EncodeBitDet(nI,nel,ilut)
      integer , intent(out) :: ilut
      integer , intent(in) :: nel
      integer , intent(in) :: nI(nel)
      integer :: i
      ilut = 0

      do i=1,nel
          ilut=ibset(ilut,nI(i)-1)
      enddo

    end subroutine EncodeBitDet

    function CountBits(ilut) result(nbits)
        integer, intent(in) :: ilut
        integer :: nbits,tmp
        integer :: m1,m2,m3,m4
        character(len=*), parameter :: t_r='CountBits'

        if(bits_n_int.eq.64) then

            m1 = 6148914691236517205        !Z'5555555555555555'
            m2 = 3689348814741910323        !Z'3333333333333333'
            m3 = 1085102592571150095        !Z'0f0f0f0f0f0f0f0f'
            m4 = 72340172838076673          !Z'0101010101010101'

            tmp = ilut - iand(ishft(ilut,-1), m1)
            tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
            tmp = iand(tmp, m3) + iand(ishft(tmp,-4), m3)
            nbits = int(ishft(tmp*m4, -56),size_n_int)

        elseif(bits_n_int.eq.32) then

            m1 = 1431655765          !Z'55555555'
            m2 = 858993459           !Z'33333333'
            m3 = 252645135           !Z'0F0F0F0F'
            m4 = 16843009            !Z'01010101'

            tmp = ilut - iand(ishft(ilut,-1), m1)
            tmp = iand(tmp, m2) + iand(ishft(tmp, -2), m2)
            tmp = iand((tmp+ishft(tmp, -4)), m3) * m4
            nbits = ishft(tmp, -24)
           ! write(6,"(B32.32,i7)") ilut,nbits
        else
            call stop_all(t_r,'Error in integer type')
        endif

    end function CountBits

end module DetBitOps
