program openmpbughunt

    !
    ! ifort openmp/threadprivate pointer testcase.
    !
    ! To build:
    !     ifort -o test -openmp test.f90
    !
    !     (tested on ifort 14.0.2 20140120, 
    !
    ! To run:
    !     OMP_NUM_THREADS=4 ./test
    !
    ! Expected:
    !     The output should loop over the various loop elements 1-10 (in
    !     differing orders). For each iteration, the expected output is (for
    !     example):
    ! 
    !  >> Loop element:  3, on thread:   0
    !  >> Pointer T  1  5
    !  >> Pointer val 1234
    !
    !
    ! Observed:
    !     The "associated intrinsic" does not appear to work on any thread
    !     other than thread 0. HOWEVER, the pointer appears to be correctly
    !     associated, the correct arry bounds can be read, and the correct
    !     value is read:
    !
    !  >> Loop element: 10, on thread:   3
    !  >> Pointer F  1  5
    !  >> Pointer val 1234
    !
    !
    ! Additional test:
    !     If array bounds checking is enabled (-CB), then sometimes the array
    !     bounds do not get correctly assigned, causing a runtime failure. This
    !     is inconsistent, and does not always occur.

    !  >> Loop element:  7, on thread:   2
    !  >> Pointer F  1  0
    !  >> forrtl: severe (408): fort: (2): Subscript #1 of the array PTR has value 1 which is greater than the upper bound of -1

    !  >> Image              PC                Routine            Line        Source             
    !  >> test               0000000000403A15  Unknown               Unknown  Unknown
    !  >> libiomp5.so        00007F6F13A2B233  Unknown               Unknown  Unknown
    !
    ! ifort versions:
    !     12.1.3 20120212 - Only the additional test fails
    !     13.0.0 20120731 - Only the additional test fails
    !     13.1.1 20130313 - Only the additional test fails
    !     13.1.3 20130607 - Only the additional test fails
    !     14.0.0 20130728 - Both tests fail
    !     14.0.1 20131008 - Both tests fail
    !     14.0.2 20140120 - Both tests fail
    !
    !     gfortran - Behaves as expected.

!$  use omp_lib
    implicit none
    integer :: i
    integer, allocatable, target :: arr(:)

    ! This must be a threadprivate pointer. The same failure is obsevrved if it
    ! is a module level threadprivate variable
    !
    ! If this is declared as a normal pointer, and passed into the parallel do
    ! loop via a "private(ptr)" declaration, then the failure is not observed.
    integer, save, pointer :: ptr(:)
!$OMP threadprivate(ptr)

    ! This failure only appears to manifest itself with an allocatable array.
    allocate(arr(5))
    arr = 1234


!$OMP parallel do default(shared)
    do i = 1, 10

        write(6,'(a,i3,a,i3)') 'Loop element:', i, ', on thread: ', &
                               OMP_get_thread_num()

        ptr => arr
        write(6,'(a,l,2i3)') "Pointer", associated(ptr), &
                                   lbound(ptr, 1), ubound(ptr, 1)
        write(6,'(a,i5)') "Pointer val", ptr(1)
        ptr => null()

    enddo

!$OMP end parallel do

end

