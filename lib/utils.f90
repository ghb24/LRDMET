module utils

! Various utilities and tools...

implicit none

contains

    elemental function binom_i(m, n) result(binom)

        ! ACM Algorithm 160 translated to Fortran.
        ! Returns the binomial coefficient ^mC_n: the number of
        ! combinations of m things taken n at a time.

        integer :: binom
        integer, intent(in) :: m, n
        integer :: p,i, n1

        n1 = n
        p = m - n1
        if (n1 < p) then
            p = n1
            n1 = m - p
        end if
        binom = n1 + 1
        if (p == 0) binom = 1
        do i = 2, p
            binom = (binom*(n1+i))/i
        end do

    end function binom_i

    elemental function binom_r(m, n) result(binom)

        ! ACM Algorithm 160 translated to Fortran.
        ! Returns the binomial coefficient ^mC_n: the number of
        ! combinations of m things taken n at a time.

        use const, only: dp

        real(dp) :: binom
        integer, intent(in) :: m, n
        integer :: p,i, n1

        n1 = n
        p = m - n1
        if (n1 < p) then
            p = n1
            n1 = m - p
        end if
        binom = n1 + 1
        if (p == 0) binom = 1
        do i = 2, p
            binom = (binom*(n1+i))/i
        end do

    end function binom_r

    function get_free_unit() result(free_unit)

        ! Returns:
        !    The first free file unit above 10 and less than or equal to
        !    the paramater max_unit (currently set to 200).

        use errors, only: stop_all

        integer, parameter :: max_unit = 100
        integer :: free_unit
        integer :: i
        logical :: t_open, t_exist

        do i = 10, max_unit
            inquire(unit=i, opened=t_open, exist=t_exist)
            if (.not.t_open .and. t_exist) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below max_unit.')

    end function get_free_unit

    elemental function char_in_int(i) result(chars)
        integer, intent(in) :: i
        integer :: chars

        if(i < 10) then
            chars = 1
        elseif(i < 100) then
            chars = 2
        elseif(i < 1000) then
            chars = 3
        elseif(i < 10000) then
            chars = 4
        elseif(i < 100000) then
            chars = 5
        elseif(i < 1000000) then
            chars = 6
        elseif(i < 10000000) then
            chars = 7
        elseif(i < 100000000) then
            chars = 8
        endif

    end function char_in_int

    elemental function int_fmt(i, padding) result(fmt1)
    
        ! In:
        !    i: an integer
        !    padding (optional): amount of padding to add to format statement.
        !        The default amount is 2.  The padding is used to include the
        !        sign if i is negative.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.
        
        ! This does take i/o formatting to a slightly OCD level addmittedly...

        character(4) :: fmt1
        integer, intent(in) :: i
        integer, intent(in), optional :: padding
        integer :: p
        real :: r

        if (present(padding)) then
            p = padding
        else
            p  = 2
        end if

        if (i == 0 .or. i==1) then
            r = 1.0
        else
            r = log10(real(abs(i)+1))
        end if
        p = ceiling(r+p)

        if (p < 10) then
            write (fmt1,'("i",i1)') p
        else if (p < 100) then
            write (fmt1,'("i",i2)') p
        else
            ! By this point we'll have hit integer overflow anyway...
            write (fmt1,'("i",i3)') p
        end if

    end function int_fmt

    subroutine append_ext(stem, n, s)

        ! Returns stem.n in s.

        character(*), intent(in) :: stem
        integer, intent(in) :: n
        character(*), intent(out) :: s
        character(10) :: ext
        !character(len=*), parameter :: t_r='append_ext'

!        i=1
!        do while(.true.)
!            if(stem(i:i).eq.' ') then
!                exit
!            else
!                s(i:i) = stem(i:i)
!            endif
!            i=i+1
!            if(i.gt.256) stop 'error in append_ext'
!        enddo
!        !i is first empty slot
!
!        s(i:i) = '_'
!        i=i+1
!
!        write (ext,'('//int_fmt(n,0)//')') n
!        do k=1,10
!            write(6,*) 'sss',ext,k,ext(k:k),char_in_int(n)
!        enddo
!        
!        do k=i,i+char_in_int(n)
!            s(k:k) = ext(k-i+1:k-i+1)
!        enddo
!
!
!        do i=1,64
!            write(6,*) stem,i,stem(i:i)
!        enddo
!
!        s = trim(s)
        ext = ''
!        !s = ''
!        !write(ext,'(I6)') n
!        !write(6,*) "S:",s
        write (ext,'('//int_fmt(n,0)//')') n
        s = adjustr(trim(stem))//'_'//adjustl(ext)
!        write(6,*) "s: ",s,int_fmt(n,0),ext

    end subroutine append_ext

    subroutine append_ext_real(stem, f, s)
        use const

        ! Returns stem.f in s

        implicit none
        character(*), intent(in) :: stem
        real(dp), intent(in) :: f
        character(*), intent(out) :: s

        s = ''
        write(s,'(F6.2)') f
        s = stem//'_'//adjustl(s)

    end subroutine append_ext_real

   subroutine get_unique_filename(stem, tnext, istart, filename)

        ! Find a filename which is either the "newest" or the next to be used.
        ! The filename is assumed to be stem.x, where x is an integer.

        ! In:
        !    stem: stem of the filename.
        !    tnext: the next unused filename is found if true, else the
        !        filename is set to be stem.x where stem.x exists and stem.x+1
        !        doesn't and x is greater than istart.
        !    istart: the integer of the first x value to check.
        !        If istart is negative, then the filename is set to be stem.x,
        !        where x = |istart+1|.  This overrides everything else.
        ! Out:
        !    filename.

        character(*), intent(in) :: stem
        logical, intent(in) :: tnext
        integer, intent(in) :: istart
        character(*), intent(out) :: filename

        integer :: i
        logical :: exists

        i = istart
        exists = .true.
        do while (exists)
            call append_ext(stem, i, filename)
            inquire(file=filename,exist=exists)
            i = i + 1
        end do

        if (.not.tnext) then
            ! actually want the last file which existed.
            ! this will return stem.istart if stem.istart doesn't exist.
            i = max(istart,i - 2)
            call append_ext(stem, i, filename)
        end if

        ! Have been asked for a specific file.
        if (istart < 0) then
            call append_ext(stem, abs(istart+1), filename)
        end if

    end subroutine get_unique_filename

end module utils
