#define srt_ind(i) ((i)+1)

!Wrapper routine to calculate FCI determinant list
!This only calculates it for one spin-type. 
subroutine GenDets(NEl,SpatOrbs)
    use DetToolsData
    implicit none
    integer, intent(in) :: NEl,SpatOrbs     !Number of electrons in active space, number of spatial orbitals in active space
    integer :: nDetAlpha
    integer, allocatable :: alphaDetList(:,:)
    integer :: nAlphaElec,AlphaDetLoop(NEl/2),i,j,k,iEl

    nFCIDet = 0
    nDetAlpha = 0
    nAlphaElec = NEl/2

    !Count number of alpha determinants
    write(6,"(A,I6,A,I6,A)") "Calculating FCI space for ",NEl," electrons in ",SpatOrbs," orbitals."

    call GenDets_R(AlphaDetLoop,nAlphaElec,SpatOrbs,1,nDetAlpha,.true.,alphaDetList)
    write(6,*) "Total number of determinants is: ",nDetAlpha*nDetAlpha

    allocate(alphaDetList(nAlphaElec,nDetAlpha))
    nDetAlpha = 0
    !Now actually calculate all the alpha determinants
    call GenDets_R(AlphaDetLoop,nAlphaElec,SpatOrbs,1,nDetAlpha,.false.,alphaDetList)
!
!    write(6,*) "Alpha determinant list is: "
!    do i=1,nDetAlpha
!        write(6,*) alphaDetList(:,i)
!    enddo

    !Now construct the full list as the product of the two lists
    nFCIDet = nDetAlpha**2
    allocate(FCIDetList(NEl,nFCIDet))
    k = 0
    do i=1,nDetAlpha
        !Alpha string (2n-1)
        do j=1,nDetAlpha
            !Beta string (even orbital indices)
            k = k + 1   !Determinant index
            do iEl=1,nAlphaElec
                FCIDetList(2*iEl-1,k) = alphaDetList(iEl,i)*2-1
                FCIDetList(2*iEl,k) = alphaDetList(iEl,j)*2
            enddo
            !sort the list
            call sort_int(FCIDetList(:,k),NEl)
        enddo
    enddo

    write(6,*) "FCI determinant list: "
    do i=1,nFCIDet
        write(6,*) FCIDetList(:,i)
    enddo

end subroutine GenDets

!iElec is the current electron
recursive subroutine GenDets_R(Det,NEl,SpatOrbs,iElec,nDetCurr,tCount,DetList)
    implicit none
    integer, intent(in) :: NEl,SpatOrbs     !Number of electrons in active space, number of spin-orbitals in active space
    integer :: iElec    !Current electron
    integer :: nDetCurr !Current number of completed determinants
    logical, intent(in) :: tCount
    integer, intent(out) :: DetList(NEl,*)
    integer :: Det(NEl)          !Current determinant being constructed
    integer :: iEl,iStartOrb

    iStartOrb = 1
    if(iElec.gt.1) iStartOrb = Det(iElec-1)+1  !If we have already got an electron in here, start from one more than the orbital occupied by the previous electron.

    do iEl=iStartOrb,SpatOrbs  !Loop through orbitals from starting electron orbtials to final orbital
        Det(iElec)=iEl   !Set current electron to this orbtial
!        write(6,*) "Setting electron number ",iElec," equal to ",iEl
        if(iElec.eq.NEl) then
            !This is the last electron
            nDetCurr = nDetCurr + 1
!            write(6,*) "Found determinant: ",Det(:)
            if(.not.tCount) then
                DetList(:,nDetCurr) = Det(:)
            endif
        else
            call GenDets_R(Det,NEl,SpatOrbs,iElec+1,nDetCurr,tCount,DetList)
        endif
    enddo
    return

end subroutine GenDets_R

subroutine GetHElement(nI,nJ,NEl,HEl)
    use const
    implicit none
    integer, intent(in) :: NEl
    integer, intent(in) :: nI(NEl),nJ(NEl)
    real(dp), intent(out) :: HEl
    real(dp) :: sltcnd_0,sltcnd_1,sltcnd_2
    integer :: Ex(2,2),iGetExcitLevel,IC
    logical :: tSign

    IC = IGETEXCITLEVEL(NI,NJ,NEL)
    Ex(1,1) = IC
    if(Ex(1,1).le.2) then
        call GetExcitation(nI,nJ,NEl,Ex,tSign)
    else
        HEl = 0.0_dp
        return
    endif

    if(IC.eq.0) then
        !Diagonal hamiltonian matrix element
        !Sum in 1 electron terms
        hel = sltcnd_0(nI,NEl)
    elseif(IC.eq.1) then
        !Single excitation
        hel = sltcnd_1(nI, Ex, tSign, NEl)
    elseif(IC.eq.2) then
        hel = sltcnd_2(Ex, tSign)
    endif

end subroutine GetHElement

function sltcnd_0 (nI,NEl) result(hel)
    use const
    use DetToolsData, only: TMat
    implicit none
    integer, intent(in) :: NEl
    integer, intent(in) :: nI(NEl)
    integer :: id(nel),gtid,i,idX,idN,j
    real(dp) :: hel,umatel

    do i=1,NEl
        id(i) = gtID(nI(i))   !turn to spatial orbital representation
    enddo

    hel = 0.0_dp
    do i=1,NEl
        hel = hel + tmat(id(i),id(i))
    enddo

    do i=1,nel-1
        do j=i+1,nel
            idX = max(id(i), id(j))
            idN = min(id(i), id(j))
            hel = hel + umatel(idN, idX, idN, idX)
        enddo
    enddo

    do i=1,nel-1
        do j=i+1,nel
            if(mod(nI(i),2).eq.mod(nI(j),2)) then
                idX = max(id(i), id(j))
                idN = min(id(i), id(j))
                hel = hel - umatel(idN, idX, idX, idN)
            endif
        enddo
    enddo

end function sltcnd_0

function sltcnd_1 (nI, ex, tSign, nel) result(hel)
    use const
    use DetToolsData, only: TMat
    implicit none
    integer, intent(in) :: nel,nI(nel),ex(2)
    logical, intent(in) :: tSign
    real(dp) :: hel,umatel
    integer :: id_ex(2), id, i,gtid

    do i=1,2
        id_ex(i) = gtid(ex(i))
    enddo

    hel = 0.0_dp
    if(mod(ex(1),2).eq.mod(ex(2),2)) then
        do i=1,nel
            if(ex(1).ne.nI(i)) then
                id = gtID(nI(i))
                hel = hel + umatel(id_ex(1),id,id_ex(2),id)
            endif
        enddo
    endif

    if(mod(ex(1),2).eq.mod(ex(2),2)) then
        do i=1,nel
            if(ex(1).ne.nI(i)) then
                if(mod(ex(1),2).eq.mod(nI(i),2)) then
                    id = gtid(nI(i))
                    hel = hel - umatel(id_ex(1),id,id,id_ex(2))
                endif
            endif
        enddo
    endif

    if(mod(ex(1),2).eq.mod(ex(2),2)) then
        hel = hel + tmat(id_ex(1),id_ex(2))
    endif
    if(tSign) hel = -hel

end function sltcnd_1

function sltcnd_2 (ex, tSign) result(hel)
    use const
    implicit none
    integer, intent(in) :: ex(2,2)
    logical, intent(in) :: tSign
    real(dp) :: hel,umatel
    integer :: id(2,2),gtid,i,j

    do i=1,2
        do j=1,2
            id(j,i) = gtid(ex(j,i))
        enddo
    enddo

    if((mod(ex(1,1),2).eq.mod(ex(2,1),2)).and.    &
       (mod(ex(1,2),2).eq.mod(ex(2,2),2))) then
        hel = umatel(id(1,1),id(1,2),id(2,1),id(2,2))
    else
        hel = 0.0_dp
    endif

    if((mod(ex(1,1),2).eq.mod(ex(2,2),2)).and.    &
       (mod(ex(1,2),2).eq.mod(ex(2,1),2))) then
        hel = hel - umatel(id(1,1),id(1,2),id(2,2),id(2,1))
    endif
    if(tSign) hel = -hel

end function sltcnd_2

!In physical notation and spatial orbitals
function umatind(i,j,k,l) result(ind)
    implicit none
    integer, intent(in) :: i,j,k,l
    integer :: ind,a,b
    
    !Combine indices I and K, ensuring I>K
     IF(I.GT.K) THEN
         A=(I*(I-1))/2+K
     ELSE
         A=(K*(K-1))/2+I
     ENDIF                                                                                         
     
     !Combine indices J and L, ensuring J>K
     IF(J.GT.L) THEN
         B=(J*(J-1))/2+L
     ELSE
         B=(L*(L-1))/2+J
     ENDIF                                                                                         
     
     !Combine (IK) and (JL) in a unique way  (k > l or if k = l then i > j)
     IF(A.GT.B) THEN
         Ind=(A*(A-1))/2+B
     ELSE
         Ind=(B*(B-1))/2+A 
     ENDIF

end function umatind

!In physical notation and spatial orbitals
function umatel(i,j,k,l) result(hel)
    use const
    use DetToolsData, only: UMat
    implicit none
    integer, intent(in) :: i,j,k,l
    real(dp) :: hel
    integer :: umatind
     hel = umat(umatind(i,j,k,l))
end function umatel

function gtid(gind) result(id)
    implicit none
    integer, intent(in) :: gind
    integer :: id
    id = (gind-1)/2 + 1
end function gtid
    
!Get antisymmetrized integral in the spinorbital basis given by a set of spatial orbitals from the AO basis
!This orbital transformation matrix is given by AOBasisTransform
! < ex(1,1) ex(1,2) || ex(2,1) ex(2,2) >
function GetHFAntisymInt_spinorb(ex,AOBasisTrans) result(HEl)
    use const
    use Globals, only: U, nSites
    implicit none
    integer, intent(in) :: ex(2,2)
    real(dp), intent(in) :: AOBasisTrans(nSites,nSites)
    real(dp) :: HEl,HEl_coul,HEl_exch
    integer :: i,j,k,l,i_spat,j_spat,k_spat,l_spat,gtid,alpha

    i = ex(1,1)
    j = ex(1,2)
    k = ex(2,1)
    l = ex(2,2)
    i_spat = gtid(i)
    j_spat = gtid(j)
    k_spat = gtid(k)
    l_spat = gtid(l)

    HEl_coul = 0.0_dp
    HEl_exch = 0.0_dp

    !First, calculate <ij|kl> if spin allowed
    if((mod(i,2).eq.mod(k,2)).and.(mod(j,2).eq.mod(l,2))) then
        !Integral is allowed
        
        do alpha = 1,nSites
            HEl_coul = HEl_coul + AOBasisTrans(alpha,i_spat)*AOBasisTrans(alpha,j_spat)*    &
                AOBasisTrans(alpha,k_spat)*AOBasisTrans(alpha,l_spat)
        enddo
        HEl_coul = HEl_coul * U
    endif

    !Now, the exchange component if spin allowed <i j | l k>
    if((mod(i,2).eq.mod(l,2)).and.(mod(j,2).eq.mod(k,2))) then
        !Integral is allowed. It will of course be the same as the coulomb term
        do alpha = 1,nSites
            HEl_exch = HEl_exch + AOBasisTrans(alpha,i_spat)*AOBasisTrans(alpha,j_spat)*    &
                AOBasisTrans(alpha,k_spat)*AOBasisTrans(alpha,l_spat)
        enddo
        HEl_exch = HEl_exch * U
    endif

    HEl = HEl_coul - HEl_exch

end function GetHFAntisymInt_spinorb

!Get integral in the spinorbital basis given by a set of spatial orbitals from the AO basis
!This orbital transformation matrix is given by AOBasisTransform
! < ex(1,1) ex(1,2) || ex(2,1) ex(2,2) >
function GetHFInt_spinorb(ex,AOBasisTrans) result(HEl)
    use const
    use Globals, only: U, nSites
    implicit none
    integer, intent(in) :: ex(2,2)
    real(dp), intent(in) :: AOBasisTrans(nSites,nSites)
    real(dp) :: HEl
    integer :: i,j,k,l,i_spat,j_spat,k_spat,l_spat,gtid,alpha

    i = ex(1,1)
    j = ex(1,2)
    k = ex(2,1)
    l = ex(2,2)
    i_spat = gtid(i)
    j_spat = gtid(j)
    k_spat = gtid(k)
    l_spat = gtid(l)

    HEl = 0.0_dp

    !First, calculate <ij|kl> if spin allowed
    if((mod(i,2).eq.mod(k,2)).and.(mod(j,2).eq.mod(l,2))) then
        !Integral is allowed
        
        do alpha = 1,nSites
            HEl = HEl + AOBasisTrans(alpha,i_spat)*AOBasisTrans(alpha,j_spat)*    &
                AOBasisTrans(alpha,k_spat)*AOBasisTrans(alpha,l_spat)
        enddo
        HEl = HEl * U
    endif

end function GetHFInt_spinorb
    

! Get the orbitals which are excited in going from I to J
! EX(1,*) are in I, and EX(2,*) are in J
! TSIGN is set to the parity of the permutations required to line up the orbitals
!   TRUE means ODD.
! If there are too many excitations to fit, then we put -excitlevel in EX(1,1) and EX(2,1)
! EX(1,1) is the max number of excitations (passed in as a parameter)
SUBROUTINE GETEXCITATION(NI,NJ,NEL,EX,TSIGN)
    IMPLICIT NONE
    INTEGER, intent(in) :: NEl,nI(NEl),nJ(NEl)
    INTEGER, intent(out) :: EX(2,*)
    LOGICAL, intent(out) :: tSign
    INTEGER iMaxExcit
    INTEGER I,J,IPAR
    INTEGER IC1,IC2
    iMaxExcit=EX(1,1)
    EX(1:2,1:iMaxExcit)=0
    IC1=0
    IC2=0
    I=1
    J=1
    IPAR=0
!    CALL WRITEDET(6,NI,NEL,.TRUE.)
!    CALL WRITEDET(6,NJ,NEL,.TRUE.)
    DO WHILE(I.LE.NEL.AND.J.LE.NEL)
!.. Differences from I to J
!   WRITE(6,*) "GE",I,J
        DO WHILE(I.LE.NEL)
           if (NI(I) >= NJ(J)) exit
           IC1=IC1+1
           IF(IC1.LE.iMaxExcit) THEN
              EX(1,IC1)=NI(I)
              IPAR=IPAR+I
           ENDIF
           I=I+1
        ENDDO
!.. Differences from J to I
        DO WHILE(I.LE.NEL.AND.J.LE.NEL)
           if (NI(I) <= NJ(J)) exit
           IC2=IC2+1
           IF(IC2.LE.iMaxExcit) THEN
              EX(2,IC2)=NJ(J)
              IPAR=IPAR+J
           ENDIF
           J=J+1
        ENDDO
        IF(I.LE.NEL.AND.J.LE.NEL) then
           if (NI(I) == NJ(J)) then
               I=I+1
               J=J+1
           endif
        ENDIF
    ENDDO
!.. Deal with remaining I
    DO WHILE(I.LE.NEL)
        IC1=IC1+1
        IF(IC1.LE.iMaxExcit) THEN
           IPAR=IPAR+I
           EX(1,IC1)=NI(I)
        ENDIF
        I=I+1
    ENDDO
!.. Deal with remaining J
    DO WHILE(J.LE.NEL)
        IC2=IC2+1
        IF(IC2.LE.iMaxExcit) THEN
           EX(2,IC2)=NJ(J)
           IPAR=IPAR+J
        ENDIF
        J=J+1
    ENDDO
    IF(iC1.GT.iMaxExcit) THEN
!.. we actually needed more space.  Just list the excitation counts (-ve)
       DO i=1,iMaxExcit
          IF (i.EQ.1) THEN
             EX(1,1)=-iC1
             EX(2,1)=-iC2
          ELSE
             EX(1,i)=0
             EX(2,i)=0
          ENDIF
       ENDDO
    ELSEIF(iC1.EQ.0) THEN
       EX(1,1)=0
       EX(2,1)=0
    ENDIF
    TSIGN=BTEST(IPAR,0)
    RETURN
END

INTEGER FUNCTION IGETEXCITLEVEL(NI,NJ,NEL)
    IMPLICIT NONE
    INTEGER I,J,NEL,NI(NEL),NJ(NEL),IC
    IC=0
    DO I=1,NEL
       DO J=1,NEL
          IF(NI(I).EQ.NJ(J)) THEN
              IC=IC+1
              EXIT
          ENDIF
       ENDDO
    ENDDO
    IGETEXCITLEVEL=NEL-IC
    RETURN
END FUNCTION iGetExcitLevel



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

