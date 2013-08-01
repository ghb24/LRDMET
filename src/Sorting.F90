#define srt_ind(i) ((i)+1)
subroutine sort_int(arr, length)
    ! sort needs auxiliary storage of length 2*log_2(n).
    use errors, only: stop_all
    implicit none
    integer, parameter :: nStackMax = 50
    integer, parameter :: nInsertionSort = 7              

    integer, intent(in) :: length
    integer, intent(inout) :: arr(length)                         

    ! Oh, how lovely it would be to be able to use push/pop and not need
    ! to guess a size of the stack to start with       
    integer :: stack(nStackMax), nstack 
    integer :: pivot, lo, hi, n, i, j                                
    ! n.b. This size statement is removed if type1 is scalar ...
    integer :: a               

    ! Temporary variables for swapping  
    integer :: tmp1                                  
    character(*), parameter :: this_routine = 'sort'        
                                                                      
    ! Initialise temporary variables if required      
    
    ! *** HACK ***                                                
    ! Workaround for gfortran compiler bug
    ! n.b. This will produce spurious warnings if run in valgrind, as
    !      a is not initialised. Not a problem in optimised build.
    ! if (custom_lt(a, a)) i = i
    ! if (custom_gt(a, a)) i = i


    ! The size of the array to sort. 
    ! N.B. Use zero based indices. To obtain ! the entry into the actual 
    ! array, multiply indices by nskip and add ! 1 (hence zero based)
    ! **** See local macro srt_ind() ******
    lo = lbound(arr, 1) - 1
    n = (ubound(arr, 1) - lo -1) + 1
    hi = lo + n - 1

    nstack = 0
    do while (.true.)
        ! If the section/partition we are looking at is smaller than
        ! nInsertSort then perform an insertion sort 
        if (hi - lo < nInsertionSort) then
            do j = lo + 1, hi
                a = arr(srt_ind(j))
                do i = j - 1, 0, -1
                    if (arr(srt_ind(i)) < a) exit
                    arr(srt_ind(i+1)) = arr(srt_ind(i))
                enddo
                arr(srt_ind(i+1)) = a
            enddo

            if (nstack == 0) exit
            hi = stack(nstack)
            lo = stack(nstack-1)
            nstack = nstack - 2

        ! Otherwise start partitioning with quicksort. 
        else
            ! Pick a partitioning element, and arrange such that
            ! arr(lo) <= arr(lo+1) <= arr(hi) 
            pivot = (lo + hi) / 2
            tmp1 = arr(srt_ind(pivot))
            arr(srt_ind(pivot)) = arr(srt_ind(lo + 1))
            arr(srt_ind(lo + 1)) = tmp1

            if (arr(srt_ind(lo)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
            endif
            if (arr(srt_ind(lo+1)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
            endif
            if (arr(srt_ind(lo)) > arr(srt_ind(lo+1))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = tmp1
            endif

            i = lo + 1
            j = hi
            a = arr(srt_ind(lo + 1)) !! a is the pivot value
            do while (.true.)
                ! Scand down list to find element > a 
                i = i + 1
                do while (arr(srt_ind(i)) < a)
                    i = i + 1
                enddo

                ! Scan down list to find element < a 
                j = j - 1
                do while (arr(srt_ind(j)) > a)
                    j = j - 1
                enddo

                ! When the pointers crossed, partitioning is complete. 
                if (j < i) exit
                ! Swap the elements, so that all elements < a end up
                ! in lower indexed variables. 
                tmp1 = arr(srt_ind(i))
                arr(srt_ind(i)) = arr(srt_ind(j))
                arr(srt_ind(j)) = tmp1
            enddo;

            ! Insert partitioning element 
            arr(srt_ind(lo + 1)) = arr(srt_ind(j))
            arr(srt_ind(j)) = a
            ! Push the larger of the partitioned sections onto the stack
            ! of sections to look at later.
            ! --> need fewest stack elements. 
            nstack = nstack + 2
            if (nstack > nStackMax) then
                    call stop_all (this_routine, &
                                    "parameter nStackMax too small")
            endif
            if (hi - i + 1 >= j - lo) then
                stack(nstack) = hi
                stack(nstack-1) = i
                hi = j - 1
            else
                stack(nstack) = j - 1
                stack(nstack-1) = lo
                lo = i
            endif
        endif
    enddo
end subroutine sort_int

subroutine sort_real(arr, length)
    ! sort needs auxiliary storage of length 2*log_2(n).
    use errors, only: stop_all
    use const
    implicit none
    integer, parameter :: nStackMax = 50
    integer, parameter :: nInsertionSort = 7              

    integer, intent(in) :: length
    real(dp), intent(inout) :: arr(length)                         

    ! Oh, how lovely it would be to be able to use push/pop and not need
    ! to guess a size of the stack to start with       
    integer :: stack(nStackMax), nstack 
    integer :: pivot, lo, hi, n, i, j                                
    ! n.b. This size statement is removed if type1 is scalar ...
    real(dp) :: a               

    ! Temporary variables for swapping  
    real(dp) :: tmp1                                  
    character(*), parameter :: this_routine = 'sort_real'        
                                                                      
    ! Initialise temporary variables if required      
    
    ! *** HACK ***                                                
    ! Workaround for gfortran compiler bug
    ! n.b. This will produce spurious warnings if run in valgrind, as
    !      a is not initialised. Not a problem in optimised build.
    ! if (custom_lt(a, a)) i = i
    ! if (custom_gt(a, a)) i = i


    ! The size of the array to sort. 
    ! N.B. Use zero based indices. To obtain ! the entry into the actual 
    ! array, multiply indices by nskip and add ! 1 (hence zero based)
    ! **** See local macro srt_ind() ******
    lo = lbound(arr, 1) - 1
    n = (ubound(arr, 1) - lo -1) + 1
    hi = lo + n - 1

    nstack = 0
    do while (.true.)
        ! If the section/partition we are looking at is smaller than
        ! nInsertSort then perform an insertion sort 
        if (hi - lo < nInsertionSort) then
            do j = lo + 1, hi
                a = arr(srt_ind(j))
                do i = j - 1, 0, -1
                    if (arr(srt_ind(i)) < a) exit
                    arr(srt_ind(i+1)) = arr(srt_ind(i))
                enddo
                arr(srt_ind(i+1)) = a
            enddo

            if (nstack == 0) exit
            hi = stack(nstack)
            lo = stack(nstack-1)
            nstack = nstack - 2

        ! Otherwise start partitioning with quicksort. 
        else
            ! Pick a partitioning element, and arrange such that
            ! arr(lo) <= arr(lo+1) <= arr(hi) 
            pivot = (lo + hi) / 2
            tmp1 = arr(srt_ind(pivot))
            arr(srt_ind(pivot)) = arr(srt_ind(lo + 1))
            arr(srt_ind(lo + 1)) = tmp1

            if (arr(srt_ind(lo)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
            endif
            if (arr(srt_ind(lo+1)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
            endif
            if (arr(srt_ind(lo)) > arr(srt_ind(lo+1))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = tmp1
            endif

            i = lo + 1
            j = hi
            a = arr(srt_ind(lo + 1)) !! a is the pivot value
            do while (.true.)
                ! Scand down list to find element > a 
                i = i + 1
                do while (arr(srt_ind(i)) < a)
                    i = i + 1
                enddo

                ! Scan down list to find element < a 
                j = j - 1
                do while (arr(srt_ind(j)) > a)
                    j = j - 1
                enddo

                ! When the pointers crossed, partitioning is complete. 
                if (j < i) exit
                ! Swap the elements, so that all elements < a end up
                ! in lower indexed variables. 
                tmp1 = arr(srt_ind(i))
                arr(srt_ind(i)) = arr(srt_ind(j))
                arr(srt_ind(j)) = tmp1
            enddo;

            ! Insert partitioning element 
            arr(srt_ind(lo + 1)) = arr(srt_ind(j))
            arr(srt_ind(j)) = a
            ! Push the larger of the partitioned sections onto the stack
            ! of sections to look at later.
            ! --> need fewest stack elements. 
            nstack = nstack + 2
            if (nstack > nStackMax) then
                    call stop_all (this_routine, &
                                    "parameter nStackMax too small")
            endif
            if (hi - i + 1 >= j - lo) then
                stack(nstack) = hi
                stack(nstack-1) = i
                hi = j - 1
            else
                stack(nstack) = j - 1
                stack(nstack-1) = lo
                lo = i
            endif
        endif
    enddo
end subroutine sort_real

subroutine sort_real2(arr, length, length2)
    ! sort needs auxiliary storage of length 2*log_2(n).
    use errors, only: stop_all
    use const
    implicit none
    integer, parameter :: nStackMax = 50
    integer, parameter :: nInsertionSort = 7              

    integer, intent(in) :: length
    integer, intent(in) :: length2
    real(dp), intent(inout) :: arr(length,length2) 

    ! Oh, how lovely it would be to be able to use push/pop and not need
    ! to guess a size of the stack to start with       
    integer :: stack(nStackMax), nstack 
    integer :: pivot, lo, hi, n, i, j                                
    ! n.b. This size statement is removed if type1 is scalar ...
    real(dp) :: a(length2)

    ! Temporary variables for swapping  
    real(dp) :: tmp1(length2)                       
    character(*), parameter :: this_routine = 'sort_real'        
                                                                      
    ! Initialise temporary variables if required      
    
    ! *** HACK ***                                                
    ! Workaround for gfortran compiler bug
    ! n.b. This will produce spurious warnings if run in valgrind, as
    !      a is not initialised. Not a problem in optimised build.
    ! if (custom_lt(a, a)) i = i
    ! if (custom_gt(a, a)) i = i


    ! The size of the array to sort. 
    ! N.B. Use zero based indices. To obtain ! the entry into the actual 
    ! array, multiply indices by nskip and add ! 1 (hence zero based)
    ! **** See local macro srt_ind() ******
    lo = lbound(arr, 1) - 1
    n = (ubound(arr, 1) - lo -1) + 1
    hi = lo + n - 1

    nstack = 0
    do while (.true.)
        ! If the section/partition we are looking at is smaller than
        ! nInsertSort then perform an insertion sort 
        if (hi - lo < nInsertionSort) then
            do j = lo + 1, hi
                a = arr(srt_ind(j),:)
                do i = j - 1, 0, -1
                    if (arr(srt_ind(i),1) < a(1)) exit
                    arr(srt_ind(i+1),:) = arr(srt_ind(i),:)
                enddo
                arr(srt_ind(i+1),:) = a(:)
            enddo

            if (nstack == 0) exit
            hi = stack(nstack)
            lo = stack(nstack-1)
            nstack = nstack - 2

        ! Otherwise start partitioning with quicksort. 
        else
            ! Pick a partitioning element, and arrange such that
            ! arr(lo) <= arr(lo+1) <= arr(hi) 
            pivot = (lo + hi) / 2
            tmp1(:) = arr(srt_ind(pivot),:)
            arr(srt_ind(pivot),:) = arr(srt_ind(lo + 1),:)
            arr(srt_ind(lo + 1),:) = tmp1(:)

            if (arr(srt_ind(lo),1) > arr(srt_ind(hi),1)) then
                tmp1(:) = arr(srt_ind(lo),:)
                arr(srt_ind(lo),:) = arr(srt_ind(hi),:)
                arr(srt_ind(hi),:) = tmp1(:)
            endif
            if (arr(srt_ind(lo+1),1) > arr(srt_ind(hi),1)) then
                tmp1(:) = arr(srt_ind(lo+1),:)
                arr(srt_ind(lo+1),:) = arr(srt_ind(hi),:)
                arr(srt_ind(hi),:) = tmp1(:)
            endif
            if (arr(srt_ind(lo),1) > arr(srt_ind(lo+1),1)) then
                tmp1(:) = arr(srt_ind(lo),:)
                arr(srt_ind(lo),:) = arr(srt_ind(lo+1),:)
                arr(srt_ind(lo+1),:) = tmp1(:)
            endif

            i = lo + 1
            j = hi
            a(:) = arr(srt_ind(lo + 1),:) !! a is the pivot value
            do while (.true.)
                ! Scand down list to find element > a 
                i = i + 1
                do while (arr(srt_ind(i),1) < a(1))
                    i = i + 1
                enddo

                ! Scan down list to find element < a 
                j = j - 1
                do while (arr(srt_ind(j),1) > a(1))
                    j = j - 1
                enddo

                ! When the pointers crossed, partitioning is complete. 
                if (j < i) exit
                ! Swap the elements, so that all elements < a end up
                ! in lower indexed variables. 
                tmp1(:) = arr(srt_ind(i),:)
                arr(srt_ind(i),:) = arr(srt_ind(j),:)
                arr(srt_ind(j),:) = tmp1(:)
            enddo;

            ! Insert partitioning element 
            arr(srt_ind(lo + 1),:) = arr(srt_ind(j),:)
            arr(srt_ind(j),:) = a(:)
            ! Push the larger of the partitioned sections onto the stack
            ! of sections to look at later.
            ! --> need fewest stack elements. 
            nstack = nstack + 2
            if (nstack > nStackMax) then
                    call stop_all (this_routine, &
                                    "parameter nStackMax too small")
            endif
            if (hi - i + 1 >= j - lo) then
                stack(nstack) = hi
                stack(nstack-1) = i
                hi = j - 1
            else
                stack(nstack) = j - 1
                stack(nstack-1) = lo
                lo = i
            endif
        endif
    enddo
end subroutine sort_real2
    
subroutine sort_d_i (arr, arr2, length)
    use errors, only: stop_all
    use const
    implicit none

    ! Perform a quicksort on an array of integers, arr(n). Uses the 
    ! sample code in NumericalRecipies as a base.
    ! Optionally sort arr2 in parallel (in the routines it is enabled)

    ! Custom comparisons. Use the operator .locallt. & .localgt.
    ! interface operator(.locallt.)
    !     pure function custom_lt (b, c) result (ret)
    !         use constants
    !         real(dp), intent(in) :: b(:)
    !         real(dp), intent(in) :: c(:)
    !         logical :: ret
    !     end function 
    ! end interface     
    ! interface operator(.localgt.)
    !     pure function custom_gt (b, c) result (ret)
    !         use constants
    !         real(dp), intent(in) :: b(:)
    !         real(dp), intent(in) :: c(:)
    !         logical :: ret
    !     end function 
    ! end interface     


    ! sort needs auxiliary storage of length 2*log_2(n).
    integer, parameter :: nStackMax = 50
    integer, parameter :: nInsertionSort = 7
    integer, parameter :: nskip_ = 1
    integer, intent(in) :: length

    real(dp), intent(inout) :: arr(length)
    integer, intent(inout) :: arr2(length)

    ! Oh, how lovely it would be to be able to use push/pop and not need
    ! to guess a size of the stack to start with
    integer :: stack(nStackMax), nstack
    integer :: pivot, lo, hi, n, i, j
    ! n.b. This size statement is removed if type1 is scalar ...
    real(dp) :: a
    integer :: a2
    ! Temporary variables for swapping
    real(dp) :: tmp1
    integer :: tmp2
    character(*), parameter :: this_routine = 'sort_d_i'

    ! The size of the array to sort. 
    ! N.B. Use zero based indices. To obtain ! the entry into the actual 
    ! array, multiply indices by nskip and add ! 1 (hence zero based)
    ! **** See local macro srt_ind() ******
    lo = lbound(arr, 1) - 1
    n = ((ubound(arr, 1) - lo -1)/nskip_) + 1
    hi = lo + n - 1

    nstack = 0
    do while (.true.)
        ! If the section/partition we are looking at is smaller than
        ! nInsertSort then perform an insertion sort 
        if (hi - lo < nInsertionSort) then
            do j = lo + 1, hi
                a = arr(srt_ind(j))
                 a2 = arr2(srt_ind(j))
                do i = j - 1, 0, -1
                    if (arr(srt_ind(i)) < a) exit
                    arr(srt_ind(i+1)) = arr(srt_ind(i))
                     arr2(srt_ind(i+1)) = arr2(srt_ind(i))
                enddo
                arr(srt_ind(i+1)) = a
                 arr2(srt_ind(i+1)) = a2
            enddo

            if (nstack == 0) exit
            hi = stack(nstack)
            lo = stack(nstack-1)
            nstack = nstack - 2

        ! Otherwise start partitioning with quicksort. 
        else
            ! Pick a partitioning element, and arrange such that
            ! arr(lo) <= arr(lo+1) <= arr(hi) 
            pivot = (lo + hi) / 2
            tmp1 = arr(srt_ind(pivot))
            arr(srt_ind(pivot)) = arr(srt_ind(lo + 1))
            arr(srt_ind(lo + 1)) = tmp1
             tmp2 = arr2(srt_ind(pivot))
             arr2(srt_ind(pivot)) = arr2(srt_ind(lo+1))
             arr2(srt_ind(lo+1)) = tmp2

            if (arr(srt_ind(lo)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
                 tmp2 = arr2(srt_ind(lo))
                 arr2(srt_ind(lo)) = arr2(srt_ind(hi))
                 arr2(srt_ind(hi)) = tmp2
            endif
            if (arr(srt_ind(lo+1)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
                 tmp2 = arr2(srt_ind(lo+1))
                 arr2(srt_ind(lo+1)) = arr2(srt_ind(hi))
                 arr2(srt_ind(hi)) = tmp2
            endif
            if (arr(srt_ind(lo)) > arr(srt_ind(lo+1))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = tmp1
                 tmp2 = arr2(srt_ind(lo))
                 arr2(srt_ind(lo)) = arr2(srt_ind(lo+1))
                 arr2(srt_ind(lo+1)) = tmp2
            endif

            i = lo + 1
            j = hi
            a = arr(srt_ind(lo + 1)) !! a is the pivot value
             a2 = arr2(srt_ind(lo + 1))
            do while (.true.)
                ! Scand down list to find element > a 
                i = i + 1
                do while (arr(srt_ind(i)) < a)
                    i = i + 1
                enddo

                ! Scan down list to find element < a 
                j = j - 1
                do while (arr(srt_ind(j)) > a)
                    j = j - 1
                enddo

                ! When the pointers crossed, partitioning is complete. 
                if (j < i) exit

                ! Swap the elements, so that all elements < a end up
                ! in lower indexed variables. 
                tmp1 = arr(srt_ind(i))
                arr(srt_ind(i)) = arr(srt_ind(j))
                arr(srt_ind(j)) = tmp1
                 tmp2 = arr2(srt_ind(i))
                 arr2(srt_ind(i)) = arr2(srt_ind(j))
                 arr2(srt_ind(j)) = tmp2
            enddo;

            ! Insert partitioning element 
            arr(srt_ind(lo + 1)) = arr(srt_ind(j))
            arr(srt_ind(j)) = a
             arr2(srt_ind(lo + 1)) = arr2(srt_ind(j))
             arr2(srt_ind(j)) = a2

            ! Push the larger of the partitioned sections onto the stack
            ! of sections to look at later.
            ! --> need fewest stack elements. 
            nstack = nstack + 2
            if (nstack > nStackMax) then
                    call stop_all (this_routine, &
                                    "parameter nStackMax too small")
            endif
            if (hi - i + 1 >= j - lo) then
                stack(nstack) = hi
                stack(nstack-1) = i
                hi = j - 1
            else
                stack(nstack) = j - 1
                stack(nstack-1) = lo
                lo = i
            endif
        endif
    enddo

end subroutine sort_d_i

subroutine sort_d_a_c (arr, arr2, length, length2)
    use errors, only: stop_all
    use const
    implicit none

    ! sort needs auxiliary storage of length 2*log_2(n).
    integer, parameter :: nStackMax = 50
    integer, parameter :: nInsertionSort = 7
    integer, parameter :: nskip_ = 1
    integer, intent(in) :: length, length2

    real(dp), intent(inout) :: arr(length)
    complex(dp), intent(inout) :: arr2(length2,length)

    ! Oh, how lovely it would be to be able to use push/pop and not need
    ! to guess a size of the stack to start with
    integer :: stack(nStackMax), nstack
    integer :: pivot, lo, hi, n, i, j
    ! n.b. This size statement is removed if type1 is scalar ...
    real(dp) :: a
    complex(dp) :: a2(size(arr2(:,1)))

    ! Temporary variables for swapping
    real(dp) :: tmp1
    complex(dp) :: tmp2(size(arr2(:,1)))
    character(*), parameter :: this_routine = 'sort_d_a_c'


    ! The size of the array to sort. 
    ! N.B. Use zero based indices. To obtain ! the entry into the actual 
    ! array, multiply indices by nskip and add ! 1 (hence zero based)
    ! **** See local macro srt_ind() ******
    lo = lbound(arr, 1) - 1
    n = ((ubound(arr, 1) - lo -1)/nskip_) + 1
    hi = lo + n - 1

    nstack = 0
    do while (.true.)
        ! If the section/partition we are looking at is smaller than
        ! nInsertSort then perform an insertion sort 
        if (hi - lo < nInsertionSort) then
            do j = lo + 1, hi
                a = arr(srt_ind(j))
                 a2 = arr2(:,srt_ind(j))
                do i = j - 1, 0, -1
                    if (arr(srt_ind(i)) < a) exit
                    arr(srt_ind(i+1)) = arr(srt_ind(i))
                     arr2(:,srt_ind(i+1)) = arr2(:,srt_ind(i))
                enddo
                arr(srt_ind(i+1)) = a
                 arr2(:,srt_ind(i+1)) = a2
            enddo

            if (nstack == 0) exit
            hi = stack(nstack)
            lo = stack(nstack-1)
            nstack = nstack - 2

        ! Otherwise start partitioning with quicksort. 
        else
            ! Pick a partitioning element, and arrange such that
            ! arr(lo) <= arr(lo+1) <= arr(hi) 
            pivot = (lo + hi) / 2
            tmp1 = arr(srt_ind(pivot))
            arr(srt_ind(pivot)) = arr(srt_ind(lo + 1))
            arr(srt_ind(lo + 1)) = tmp1
             tmp2 = arr2(:,srt_ind(pivot))
             arr2(:,srt_ind(pivot)) = arr2(:,srt_ind(lo+1))
             arr2(:,srt_ind(lo+1)) = tmp2

            if (arr(srt_ind(lo)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
                 tmp2 = arr2(:,srt_ind(lo))
                 arr2(:,srt_ind(lo)) = arr2(:,srt_ind(hi))
                 arr2(:,srt_ind(hi)) = tmp2
            endif
            if (arr(srt_ind(lo+1)) > arr(srt_ind(hi))) then
                tmp1 = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = arr(srt_ind(hi))
                arr(srt_ind(hi)) = tmp1
                 tmp2 = arr2(:,srt_ind(lo+1))
                 arr2(:,srt_ind(lo+1)) = arr2(:,srt_ind(hi))
                 arr2(:,srt_ind(hi)) = tmp2
            endif
            if (arr(srt_ind(lo)) > arr(srt_ind(lo+1))) then
                tmp1 = arr(srt_ind(lo))
                arr(srt_ind(lo)) = arr(srt_ind(lo+1))
                arr(srt_ind(lo+1)) = tmp1
                 tmp2 = arr2(:,srt_ind(lo))
                 arr2(:,srt_ind(lo)) = arr2(:,srt_ind(lo+1))
                 arr2(:,srt_ind(lo+1)) = tmp2
            endif

            i = lo + 1
            j = hi
            a = arr(srt_ind(lo + 1)) !! a is the pivot value
             a2 = arr2(:,srt_ind(lo + 1))
            do while (.true.)
                ! Scand down list to find element > a 
                i = i + 1
                do while (arr(srt_ind(i)) < a)
                    i = i + 1
                enddo

                ! Scan down list to find element < a 
                j = j - 1
                do while (arr(srt_ind(j)) > a)
                    j = j - 1
                enddo

                ! When the pointers crossed, partitioning is complete. 
                if (j < i) exit

                ! Swap the elements, so that all elements < a end up
                ! in lower indexed variables. 
                tmp1 = arr(srt_ind(i))
                arr(srt_ind(i)) = arr(srt_ind(j))
                arr(srt_ind(j)) = tmp1
                 tmp2 = arr2(:,srt_ind(i))
                 arr2(:,srt_ind(i)) = arr2(:,srt_ind(j))
                 arr2(:,srt_ind(j)) = tmp2
            enddo;

            ! Insert partitioning element 
            arr(srt_ind(lo + 1)) = arr(srt_ind(j))
            arr(srt_ind(j)) = a
             arr2(:,srt_ind(lo + 1)) = arr2(:,srt_ind(j))
             arr2(:,srt_ind(j)) = a2

            ! Push the larger of the partitioned sections onto the stack
            ! of sections to look at later.
            ! --> need fewest stack elements. 
            nstack = nstack + 2
            if (nstack > nStackMax) then
                    call stop_all (this_routine, &
                                    "parameter nStackMax too small")
            endif
            if (hi - i + 1 >= j - lo) then
                stack(nstack) = hi
                stack(nstack-1) = i
                hi = j - 1
            else
                stack(nstack) = j - 1
                stack(nstack-1) = lo
                lo = i
            endif
        endif
    enddo

end subroutine sort_d_a_c




!Ensure that all LVecs are orthonormal with respect to the RVec vectors
!LVec vectors need to be complex conjugated.
!This will involve gram-schmidt rotation of the vectors in degenerate sets
subroutine Orthonorm_zgeev_vecs(N,W,LVec,RVec)
    use errors, only: stop_all 
    use const
    implicit none
    integer, intent(in) :: N
    complex(dp), intent(in) :: W(N)
    complex(dp), intent(inout) :: LVec(N,N),RVec(N,N)
    integer :: i,StartingInd,R,j,k
    complex(dp) :: zdotc,norm,overlap

    i = 1
    do while(i.le.N)

        StartingInd = i

        if(i.lt.N) then
            do while(abs(W(i+1)-W(i)).lt.1.0e-9_dp)
                !We are in a degenerate set
                i = i+1
                if(i.eq.N) exit !We have found the last degenerate block
            enddo
        endif
!        write(6,*) "Degenerate block from ",StartingInd,' to ',i
!        call flush(6)
        
        do R = StartingInd,i
            !Now, normalize the vectors
            !Remember that the square root of a complex number will have two roots, given by +- w
            norm = zzero
            do j = 1,N
                norm = norm + dconjg(LVec(j,R)) * RVec(j,R)
            enddo
            !norm = zdotc(N,LVec(:,R),1,RVec(:,R),1)
            norm = sqrt(norm)
            !write(6,*) "Norm: ",R,norm
            RVec(:,R) = RVec(:,R) / norm
            LVec(:,R) = LVec(:,R) / dconjg(norm)
        enddo

        !call writematrixcomp(LVec(:,StartingInd:i),'LVecs',.true.)
        !call writematrixcomp(RVec(:,StartingInd:i),'RVecs',.true.)

!        !The degenerate set goes from StartingInd to i
!        !Orthogonalize the R vectors against the L vectors
!        do R = StartingInd,i    !Run through the R vectors
!            do L = StartingInd,i    !Run through the L vectors
!                if(R.eq.L) cycle    !We want to maintain overlap with the corresponding vector
!
!                !Gram-Schmidt
!                overlap = zdotc(N,LVec(:,L),1,RVec(:,R),1)
!                RVec(:,R) = RVec(:,R) - overlap * dconjg(LVec(:,L))
!                write(6,*) "Overlap: ",L,R,overlap
!            enddo
!        enddo

!Perform two-sided modified gram-schmidt to bi-orthogonalize the degenerate vectors
        do j = StartingInd,i
            do k = StartingInd,j-1    !Bi-orthogonalize against the previous vectors
                overlap = zdotc(N,LVec(:,k),1,RVec(:,j),1)
                RVec(:,j) = RVec(:,j) - RVec(:,k)*overlap
                overlap = zdotc(N,RVec(:,k),1,LVec(:,j),1)
                LVec(:,j) = LVec(:,j) - LVec(:,k)*overlap
            enddo
        enddo

!        do j = StartingInd,i
!            do k = StartingInd,i
!                if(j.eq.k) cycle
!                overlap = zdotc(N,LVec(:,j),1,RVec(:,k),1) 
!                if(abs(overlap).gt.1.0e-8_dp) then
!                    write(6,*) "i,j: ",j,k,Overlap
!                    call stop_all('Orthonorm_zgeev_vecs','non-zero overlap')
!                endif
!            enddo
!        enddo
!            
        do R = StartingInd,i
            !Now, normalize the vectors
            !Remember that the square root of a complex number will have two roots, given by +- w
            norm = zdotc(N,LVec(:,R),1,RVec(:,R),1)
            norm = sqrt(norm)
            !write(6,*) "Norm: ",R,norm
            RVec(:,R) = RVec(:,R) / norm
            LVec(:,R) = LVec(:,R) / dconjg(norm)
        enddo

        i = i+1

    enddo

end subroutine Orthonorm_zgeev_vecs


module sort_mod_c_a_c_a_c
    implicit none

    ! Private operator for sorting complex numbers by their real values
    interface operator(.gt.)
        module procedure cmplx_gt_c_a_c_a_c
    end interface

    interface operator(.lt.)
        module procedure cmplx_lt_c_a_c_a_c
    end interface

    interface cmplx_gt
        module procedure cmplx_gt_c_a_c_a_c
    end interface

    interface cmplx_lt
        module procedure cmplx_lt_c_a_c_a_c
    end interface

contains

    subroutine Order_zgeev_vecs(arr, arr2, arr3,nskip)
        use errors, only: stop_all
        use const
        implicit none
        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: nStackMax = 50
        integer, parameter :: nInsertionSort = 7
        integer, intent(in), optional :: nskip

        complex(dp), intent(inout) :: arr(:)
        complex(dp), intent(inout) :: arr2(:,:)
        complex(dp), intent(inout) :: arr3(:,:)

        ! Oh, how lovely it would be to be able to use push/pop and not need
        ! to guess a size of the stack to start with
        integer :: stack(nStackMax), nstack, nskip_
        integer :: pivot, lo, hi, n, i, j
        ! n.b. This size statement is removed if type1 is scalar ...
        complex(dp) :: a
        complex(dp) :: a2(size(arr2(:,1)))
        complex(dp) :: a3(size(arr3(:,1)))

        ! Temporary variables for swapping
        complex(dp) :: tmp1
        complex(dp) :: tmp2(size(arr2(:,1)))
        complex(dp) :: tmp3(size(arr3(:,1)))
        character(*), parameter :: this_routine = 'Order_zgeev_vecs'

        if (present(nskip)) then
            nskip_ = nskip
        else
            nskip_ = 1
        endif

        ! The size of the array to sort. 
        ! N.B. Use zero based indices. To obtain ! the entry into the actual 
        ! array, multiply indices by nskip and add ! 1 (hence zero based)
        ! **** See local macro srt_ind() ******
        lo = lbound(arr, 1) - 1
        n = ((ubound(arr, 1) - lo -1)/nskip_) + 1
        hi = lo + n - 1

        nstack = 0
        do while (.true.)
            ! If the section/partition we are looking at is smaller than
            ! nInsertSort then perform an insertion sort 
            if (hi - lo < nInsertionSort) then
                do j = lo + 1, hi
                    a = arr(srt_ind(j))
                     a2 = arr2(:,srt_ind(j))
                     a3 = arr3(:,srt_ind(j))
                    do i = j - 1, 0, -1
                        if (arr(srt_ind(i)) < a) exit
                        arr(srt_ind(i+1)) = arr(srt_ind(i))
                         arr2(:,srt_ind(i+1)) = arr2(:,srt_ind(i))
                         arr3(:,srt_ind(i+1)) = arr3(:,srt_ind(i))
                    enddo
                    arr(srt_ind(i+1)) = a
                     arr2(:,srt_ind(i+1)) = a2
                     arr3(:,srt_ind(i+1)) = a3
                enddo

                if (nstack == 0) exit
                hi = stack(nstack)
                lo = stack(nstack-1)
                nstack = nstack - 2

            ! Otherwise start partitioning with quicksort. 
            else
                ! Pick a partitioning element, and arrange such that
                ! arr(lo) <= arr(lo+1) <= arr(hi) 
                pivot = (lo + hi) / 2
                tmp1 = arr(srt_ind(pivot))
                arr(srt_ind(pivot)) = arr(srt_ind(lo + 1))
                arr(srt_ind(lo + 1)) = tmp1
                 tmp2 = arr2(:,srt_ind(pivot))
                 arr2(:,srt_ind(pivot)) = arr2(:,srt_ind(lo+1))
                 arr2(:,srt_ind(lo+1)) = tmp2
                 tmp3 = arr3(:,srt_ind(pivot))
                 arr3(:,srt_ind(pivot)) = arr3(:,srt_ind(lo+1))
                 arr3(:,srt_ind(lo+1)) = tmp3

                if (arr(srt_ind(lo)) > arr(srt_ind(hi))) then
                    tmp1 = arr(srt_ind(lo))
                    arr(srt_ind(lo)) = arr(srt_ind(hi))
                    arr(srt_ind(hi)) = tmp1
                     tmp2 = arr2(:,srt_ind(lo))
                     arr2(:,srt_ind(lo)) = arr2(:,srt_ind(hi))
                     arr2(:,srt_ind(hi)) = tmp2
                     tmp3 = arr3(:,srt_ind(lo))
                     arr3(:,srt_ind(lo)) = arr3(:,srt_ind(hi))
                     arr3(:,srt_ind(hi)) = tmp3
                endif
                if (arr(srt_ind(lo+1)) > arr(srt_ind(hi))) then
                    tmp1 = arr(srt_ind(lo+1))
                    arr(srt_ind(lo+1)) = arr(srt_ind(hi))
                    arr(srt_ind(hi)) = tmp1
                     tmp2 = arr2(:,srt_ind(lo+1))
                     arr2(:,srt_ind(lo+1)) = arr2(:,srt_ind(hi))
                     arr2(:,srt_ind(hi)) = tmp2
                     tmp3 = arr3(:,srt_ind(lo+1))
                     arr3(:,srt_ind(lo+1)) = arr3(:,srt_ind(hi))
                     arr3(:,srt_ind(hi)) = tmp3
                endif
                if (arr(srt_ind(lo)) > arr(srt_ind(lo+1))) then
                    tmp1 = arr(srt_ind(lo))
                    arr(srt_ind(lo)) = arr(srt_ind(lo+1))
                    arr(srt_ind(lo+1)) = tmp1
                     tmp2 = arr2(:,srt_ind(lo))
                     arr2(:,srt_ind(lo)) = arr2(:,srt_ind(lo+1))
                     arr2(:,srt_ind(lo+1)) = tmp2
                     tmp3 = arr3(:,srt_ind(lo))
                     arr3(:,srt_ind(lo)) = arr3(:,srt_ind(lo+1))
                     arr3(:,srt_ind(lo+1)) = tmp3
                endif

                i = lo + 1
                j = hi
                a = arr(srt_ind(lo + 1)) !! a is the pivot value
                 a2 = arr2(:,srt_ind(lo + 1))
                 a3 = arr3(:,srt_ind(lo + 1))
                do while (.true.)
                    ! Scand down list to find element > a 
                    i = i + 1
                    do while (arr(srt_ind(i)) < a)
                        i = i + 1
                    enddo

                    ! Scan down list to find element < a 
                    j = j - 1
                    do while (arr(srt_ind(j)) > a)
                        j = j - 1
                    enddo

                    ! When the pointers crossed, partitioning is complete. 
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables. 
                    tmp1 = arr(srt_ind(i))
                    arr(srt_ind(i)) = arr(srt_ind(j))
                    arr(srt_ind(j)) = tmp1
                     tmp2 = arr2(:,srt_ind(i))
                     arr2(:,srt_ind(i)) = arr2(:,srt_ind(j))
                     arr2(:,srt_ind(j)) = tmp2
                     tmp3 = arr3(:,srt_ind(i))
                     arr3(:,srt_ind(i)) = arr3(:,srt_ind(j))
                     arr3(:,srt_ind(j)) = tmp3
                enddo;

                ! Insert partitioning element 
                arr(srt_ind(lo + 1)) = arr(srt_ind(j))
                arr(srt_ind(j)) = a
                 arr2(:,srt_ind(lo + 1)) = arr2(:,srt_ind(j))
                 arr3(:,srt_ind(lo + 1)) = arr3(:,srt_ind(j))
                 arr2(:,srt_ind(j)) = a2
                 arr3(:,srt_ind(j)) = a3

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements. 
                nstack = nstack + 2
                if (nstack > nStackMax) then
                        call stop_all (this_routine, &
                                        "parameter nStackMax too small")
                endif
                if (hi - i + 1 >= j - lo) then
                    stack(nstack) = hi
                    stack(nstack-1) = i
                    hi = j - 1
                else
                    stack(nstack) = j - 1
                    stack(nstack-1) = lo
                    lo = i
                endif
            endif
        enddo

    end subroutine Order_zgeev_vecs

    ! A private comparison. This should not be available outside of the
    ! module!
    elemental function cmplx_gt_c_a_c_a_c (a, b) result (bGt)
        use const
        complex(dp), intent(in) :: a, b
        logical :: bGt

        bGt = real(a, dp) > real(b, dp)
    end function

    elemental function cmplx_lt_c_a_c_a_c (a, b) result (bLt)
        use const
        complex(dp), intent(in) :: a, b
        logical :: bLt

        bLt = real(a, dp) < real(b, dp)
    end function
end module
