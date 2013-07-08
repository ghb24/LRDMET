module DetBitOps
    use const, only: n_int,bits_n_int,end_n_int,size_n_int
    use errors, only: stop_all 
    implicit none

    contains
    
    pure subroutine DecodeBitDet(nI,nel,ilut)
      integer , intent(in) :: ilut
      integer , intent(in) :: nel
      integer , intent(out) :: nI(nel)
      integer :: elec,j

      nI(:) = 0
      elec = 0
      do j=0,end_n_int
          if(btest(iLut,j)) then
              !An electron is at this orbital
              elec=elec+1
              nI(elec)=j+1
              if (elec == nel) exit
          endif
      enddo

    end subroutine DecodeBitDet

    !Routine to create/annihilate orbital 'orb' from bit representation ilut
    !ilut is returned as the new occupation vector
    !tAnn = .T. - annihilate
    !tAnn = .F. - create
    !tSign is the change in parity of the vector
    subroutine SQOperator(ilut,orb,tSign,tAnn)
        integer, intent(inout) :: ilut
        integer, intent(in) :: orb
        logical, intent(in) :: tAnn
        logical, intent(out) :: tSign
        integer :: i,setorbs
        character(len=*), parameter :: t_r='SQOperator'

        if(tAnn.and.(.not.btest(ilut,orb-1))) then
            call stop_all(t_r,'Orbital not occupied for annihilation')
        endif
        if((.not.tAnn).and.(btest(ilut,orb-1))) then
            call stop_all(t_r,'Orbital not unoccupied for creation')
        endif

        !orb is bit 'orb-1'
        !Calculate parity
        setorbs = 0
        do i=0,orb-2
            if(btest(ilut,i)) then
                setorbs = setorbs + 1
            endif
        enddo
        if(mod(setorbs,2).eq.0) then
            !Even number of set spin-orbitals before desired orb. No change in parity
            tSign = .false.
        else
            tSign = .true.
        endif

        if(tAnn) then
            ilut = ibclr(ilut,orb-1)
        else
            ilut = ibset(ilut,orb-1)
        endif

    end subroutine SQOperator

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

    function FindBitExcitLevel(ilut1,ilut2) result(IC)
        implicit none
        integer, intent(in) :: ilut1,ilut2
        integer :: IC,tmp

        tmp = ieor(ilut1,ilut2)
        tmp = iand(ilut1,tmp)
        IC = CountBits(tmp)

    end function FindBitExcitLevel

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
