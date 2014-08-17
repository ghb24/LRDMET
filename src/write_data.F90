module writedata
    use const
    implicit none

    interface WriteMatrix
        module procedure WriteMatrix_r
        module procedure WriteMatrix_z
        module procedure WriteMatrix_i
    end interface

    interface WriteVector
        module procedure WriteVector_r
        module procedure WriteVector_z
        module procedure WriteVector_i
    end interface

    contains

    subroutine WriteMatrix_z(mat,matname,tOneLine)
        implicit none
        complex(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        logical, intent(in) :: tOneLine
        integer :: i,j

        write(6,*) "Writing out complex matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(2G18.7)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,2G18.7)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrix_z

    subroutine WriteMatrix_r(mat,matname,tOneLine)
        implicit none
        real(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        logical, intent(in) :: tOneLine
        integer :: i,j

        write(6,*) "Writing out real matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(G25.10)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,G25.10)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrix_r

    subroutine WriteMatrix_i(mat,matname,tOneLine)
        implicit none
        integer, intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        logical, intent(in) :: tOneLine
        integer :: i,j

        write(6,*) "Writing out integer matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(I11)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,I11)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrix_i

    subroutine WriteVector_r(vec,vecname)
        implicit none
        real(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out real vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVector_r
    
    subroutine WriteVector_z(vec,vecname)
        implicit none
        complex(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out complex vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(2G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVector_z
    
    subroutine WriteVector_i(vec,vecname)
        implicit none
        integer, intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out integer vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(I12)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVector_i

end module writedata
