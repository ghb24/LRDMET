#define srt_ind(i) ((i)+1)

!Wrapper routine to calculate FCI determinant list
!This only calculates it for one spin-type. 
!If tCoupledSpaces is set, then it will also create the N+1 and N-1 spaces (both Ms = +- 1/2)
!tSplitMs will split the N+1 and N-1 spaces into seperate Ms components
subroutine GenDets(NEl,SpatOrbs,tCoupledSpaces,tCreateBitRep,tSplitMs)
    use DetToolsData
    use DetBitOps, only: EncodeBitDet
    use errors, only: stop_all 
    implicit none
    integer, intent(in) :: NEl,SpatOrbs     !Number of electrons in active space, number of spatial orbitals in active space
    logical, intent(in) :: tCoupledSpaces,tCreateBitRep,tSplitMs
    integer :: nDetAlpha,Nm1DetAlpha,Np1DetAlpha
    !The list of spatial orbital occupations for one spin type of the N-electron space
    integer, allocatable :: NspinDetList(:,:)   
    integer, allocatable :: NspinDetLoop(:) !Current Det for the N-electron space
    integer, allocatable :: Nm1spinDetList(:,:)   
    integer, allocatable :: Np1spinDetList(:,:)   
    integer :: nSpinElec,Nm1SpinElec,Np1SpinElec,i,j,k,iEl
    character(len=*), parameter :: t_r='GenDets'

    nFCIDet = 0
    nDetAlpha = 0
    nSpinElec = NEl/2
    allocate(NspinDetLoop(nSpinElec))
    NspinDetLoop(:) = 0

    !Count number of alpha determinants
    write(6,"(A,I6,A,I6,A)") "Calculating FCI space for ",NEl," electrons in ",SpatOrbs," orbitals."

    call GenDets_R(NspinDetLoop,nSpinElec,SpatOrbs,1,nDetAlpha,.true.,NspinDetList)
    write(6,*) "Total number of determinants is: ",nDetAlpha*nDetAlpha

    allocate(NspinDetList(nSpinElec,nDetAlpha))
    nDetAlpha = 0
    !Now actually calculate all the alpha determinants
    call GenDets_R(NspinDetLoop,nSpinElec,SpatOrbs,1,nDetAlpha,.false.,NspinDetList)
    deallocate(NspinDetLoop)
!
!    write(6,*) "Alpha determinant list is: "
!    do i=1,nDetAlpha
!        write(6,*) NspinDetList(:,i)
!    enddo

    !Now construct the full list as the product of the two lists
    nFCIDet = nDetAlpha**2
    if(allocated(FCIDetList)) deallocate(FCIDetList)
    allocate(FCIDetList(NEl,nFCIDet))
    k = 0
    do i=1,nDetAlpha
        !Alpha string (2n-1)
        do j=1,nDetAlpha
            !Beta string (even orbital indices)
            k = k + 1   !Determinant index
            do iEl=1,nSpinElec
                FCIDetList(2*iEl-1,k) = NspinDetList(iEl,i)*2-1
                FCIDetList(2*iEl,k) = NspinDetList(iEl,j)*2
            enddo
            !sort the list
            call sort_int(FCIDetList(:,k),NEl)
        enddo
    enddo

!    write(6,*) "FCI determinant list: "
!    do i=1,nFCIDet
!        write(6,*) FCIDetList(:,i)
!    enddo

    if(tCreateBitRep) then
        if(allocated(FCIBitList)) deallocate(FCIBitList)
        allocate(FCIBitList(nFCIDet))
        do i=1,nFCIDet
            call EncodeBitDet(FCIDetList(:,i),NEl,FCIBitList(i))
        enddo
    endif
    DetListStorage = DetListStorage + ((NEl+1)*nFCIDet)

    if(.not.tCoupledSpaces) then
        deallocate(NspinDetList)
        return  !We do not want to create the N-1 and N+1 spaces
    endif

    !Now, create the N-1 space
    write(6,*) "Calculating N-1 space"
    write(6,"(A,I6,A,I6,A)") "Calculating FCI space for ",NEl-1," electrons in ",SpatOrbs," orbitals."

    !The NspinDetList is now for the N eletron space. Make something equivalent for the N-1 electron space
    nNm1FCIDet = 0
    Nm1DetAlpha = 0
    Nm1SpinElec = (NEl/2)-1

    if(Nm1SpinElec.gt.0) then
        allocate(NspinDetLoop(Nm1SpinElec))
        NspinDetLoop(:) = 0
        call GenDets_R(NspinDetLoop,Nm1SpinElec,SpatOrbs,1,Nm1DetAlpha,.true.,Nm1spinDetList)
        write(6,*) "Total number of determinants in N-1 space is: ",Nm1DetAlpha*nDetAlpha*2
        allocate(Nm1spinDetList(Nm1SpinElec,Nm1DetAlpha))
        Nm1DetAlpha = 0
        !Now actually calculate all the alpha determinants
        call GenDets_R(NspinDetLoop,Nm1SpinElec,SpatOrbs,1,Nm1DetAlpha,.false.,Nm1spinDetList)
        deallocate(NspinDetLoop)
        if(tSplitMs) then
            nNm1FCIDet = Nm1DetAlpha*nDetAlpha 
            nNm1bFCIDet = Nm1DetAlpha*nDetAlpha 
            if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
            allocate(Nm1FCIDetList(NEl-1,nNm1FCIDet))
            if(allocated(Nm1bFCIDetList)) deallocate(Nm1bFCIDetList)
            allocate(Nm1bFCIDetList(NEl-1,nNm1bFCIDet))

            k=1
            do i=1,nDetAlpha
                do j=1,Nm1DetAlpha
                    do iEl=1,Nm1SpinElec
                        Nm1bFCIDetList(iEl,k) = Nm1spinDetList(iEl,j)*2-1    !N-1 alpha orbitals
                        Nm1FCIDetList(iEl,k) = Nm1spinDetList(iEl,j)*2    !N-1 beta orbitals
                    enddo
                    do iEl=1,nSpinElec
                        Nm1bFCIDetList(Nm1SpinElec+iEl,k) = NspinDetList(iEl,i)*2
                        Nm1FCIDetList(Nm1SpinElec+iEl,k) = NspinDetList(iEl,i)*2-1
                    enddo
                    call sort_int(Nm1bFCIDetList(:,k),NEl-1)
                    call sort_int(Nm1FCIDetList(:,k),NEl-1)
                    k=k+1
                enddo
            enddo
            if(k.ne.nNm1FCIDet+1) call stop_all(t_r,'Error filling split det lists')
        else
            nNm1FCIDet = Nm1DetAlpha*nDetAlpha*2    !Because we need both Ms values
            if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
            allocate(Nm1FCIDetList(NEl-1,nNm1FCIDet))

            k=1
            do i=1,nDetAlpha
                do j=1,Nm1DetAlpha
                    do iEl=1,Nm1SpinElec
                        Nm1FCIDetList(iEl,k) = Nm1spinDetList(iEl,j)*2-1    !N-1 alpha orbitals
                        Nm1FCIDetList(iEl,k+1) = Nm1spinDetList(iEl,j)*2    !N-1 beta orbitals
                    enddo
                    do iEl=1,nSpinElec
                        Nm1FCIDetList(Nm1SpinElec+iEl,k) = NspinDetList(iEl,i)*2
                        Nm1FCIDetList(Nm1SpinElec+iEl,k+1) = NspinDetList(iEl,i)*2-1
                    enddo
                    call sort_int(Nm1FCIDetList(:,k),NEl-1)
                    call sort_int(Nm1FCIDetList(:,k+1),NEl-1)
                    k=k+2
                enddo
            enddo
        endif
        deallocate(Nm1spinDetList)
    else
        !Only 1 electron
        Nm1DetAlpha = 0
        write(6,*) "Total number of determinants in N-1 space is: ",nDetAlpha*2
        if(tSplitMs) then
            nNm1FCIDet = nDetAlpha
            nNm1bFCIDet = nDetAlpha
            if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
            allocate(Nm1FCIDetList(NEl-1,nNm1FCIDet))
            if(allocated(Nm1bFCIDetList)) deallocate(Nm1bFCIDetList)
            allocate(Nm1bFCIDetList(NEl-1,nNm1bFCIDet))
            k=1
            do i=1,nDetAlpha
                !We only have alpha/beta orbtials at a time (1 electron)
                Nm1bFCIDetList(1,k) = NspinDetList(1,i)*2
                Nm1FCIDetList(1,k) = NspinDetList(1,i)*2-1
                k=k+1
            enddo
        else
            Nm1DetAlpha = 0
            write(6,*) "Total number of determinants in N-1 space is: ",nDetAlpha*2
            nNm1FCIDet = nDetAlpha*2
            if(allocated(Nm1FCIDetList)) deallocate(Nm1FCIDetList)
            allocate(Nm1FCIDetList(NEl-1,nNm1FCIDet))
            k=1
            do i=1,nDetAlpha
                !We only have alpha/beta orbtials at a time (1 electron)
                Nm1FCIDetList(1,k) = NspinDetList(1,i)*2
                Nm1FCIDetList(1,k+1) = NspinDetList(1,i)*2-1
                k=k+2
            enddo
        endif

    endif

    if(tCreateBitRep) then
        if(allocated(Nm1BitList)) deallocate(Nm1BitList)
        allocate(Nm1BitList(nNm1FCIDet))
        if(tSplitMs) then
            if(allocated(Nm1bBitList)) deallocate(Nm1bBitList)
            allocate(Nm1bBitList(nNm1bFCIDet))
        endif
        do i=1,nNm1FCIDet
            call EncodeBitDet(Nm1FCIDetList(:,i),NEl-1,Nm1BitList(i))
        enddo
        if(tSplitMs) then
            do i=1,nNm1bFCIDet
                call EncodeBitDet(Nm1bFCIDetList(:,i),NEl-1,Nm1bBitList(i))
            enddo
        endif
    endif
    DetListStorage = DetListStorage + (NEl*nNm1FCIDet)
    if(tSplitMs) DetListStorage = DetListStorage + (NEl*nNm1bFCIDet)

    !Now for the N+1 space
    write(6,*) "Calculating N+1 space"
    write(6,"(A,I6,A,I6,A)") "Calculating FCI space for ",NEl+1," electrons in ",SpatOrbs," orbitals."

    !The NspinDetList is now for the N eletron space. Make something equivalent for the N+1 electron space
    nNp1FCIDet = 0
    Np1DetAlpha = 0
    Np1SpinElec = (NEl/2)+1
    allocate(NspinDetLoop(Np1SpinElec))
    NspinDetLoop(:) = 0

    call GenDets_R(NspinDetLoop,Np1SpinElec,SpatOrbs,1,Np1DetAlpha,.true.,Np1spinDetList)
    write(6,*) "Total number of determinants in N+1 space is: ",Np1DetAlpha*nDetAlpha*2
    allocate(Np1spinDetList(Np1SpinElec,Np1DetAlpha))
    Np1DetAlpha = 0
    !Now actually calculate all the alpha determinants
    call GenDets_R(NspinDetLoop,Np1SpinElec,SpatOrbs,1,Np1DetAlpha,.false.,Np1spinDetList)
    deallocate(NspinDetLoop)

    if(tSplitMs) then
        nNp1FCIDet = Np1DetAlpha*nDetAlpha
        nNp1bFCIDet = Np1DetAlpha*nDetAlpha
        if(allocated(Np1FCIDetList)) deallocate(Nm1FCIDetList)
        allocate(Np1FCIDetList(NEl+1,nNp1FCIDet))
        if(allocated(Np1bFCIDetList)) deallocate(Nm1bFCIDetList)
        allocate(Np1bFCIDetList(NEl+1,nNp1bFCIDet))

        k=1
        do i=1,nDetAlpha
            do j=1,Np1DetAlpha
                do iEl=1,Np1SpinElec
                    Np1FCIDetList(iEl,k) = Np1spinDetList(iEl,j)*2-1    !N-1 alpha orbitals
                    Np1bFCIDetList(iEl,k) = Np1spinDetList(iEl,j)*2    !N-1 beta orbitals
                enddo
                do iEl=1,nSpinElec
                    Np1FCIDetList(Np1SpinElec+iEl,k) = NspinDetList(iEl,i)*2
                    Np1bFCIDetList(Np1SpinElec+iEl,k) = NspinDetList(iEl,i)*2-1
                enddo
                call sort_int(Np1FCIDetList(:,k),NEl+1)
                call sort_int(Np1bFCIDetList(:,k),NEl+1)
                k=k+1
            enddo
        enddo
        if(k.ne.nNp1FCIDet+1) call stop_all(t_r,'Error counting Ms Dets N+1')
        deallocate(Np1spinDetList,NspinDetList)
        if(tCreateBitRep) then
            if(allocated(Np1BitList)) deallocate(Np1BitList)
            allocate(Np1BitList(nNp1FCIDet))
            if(allocated(Np1bBitList)) deallocate(Np1bBitList)
            allocate(Np1bBitList(nNp1bFCIDet))
            do i=1,nNp1FCIDet
                call EncodeBitDet(Np1FCIDetList(:,i),NEl+1,Np1BitList(i))
                call EncodeBitDet(Np1bFCIDetList(:,i),NEl+1,Np1bBitList(i))
            enddo
        endif
    else
        nNp1FCIDet = Np1DetAlpha*nDetAlpha*2    !Because we need both Ms values
        if(allocated(Np1FCIDetList)) deallocate(Nm1FCIDetList)
        allocate(Np1FCIDetList(NEl+1,nNp1FCIDet))

        k=1
        do i=1,nDetAlpha
            do j=1,Np1DetAlpha
                do iEl=1,Np1SpinElec
                    Np1FCIDetList(iEl,k) = Np1spinDetList(iEl,j)*2-1    !N-1 alpha orbitals
                    Np1FCIDetList(iEl,k+1) = Np1spinDetList(iEl,j)*2    !N-1 beta orbitals
                enddo
                do iEl=1,nSpinElec
                    Np1FCIDetList(Np1SpinElec+iEl,k) = NspinDetList(iEl,i)*2
                    Np1FCIDetList(Np1SpinElec+iEl,k+1) = NspinDetList(iEl,i)*2-1
                enddo
                call sort_int(Np1FCIDetList(:,k),NEl+1)
                call sort_int(Np1FCIDetList(:,k+1),NEl+1)
                k=k+2
            enddo
        enddo
        deallocate(Np1spinDetList,NspinDetList)
        if(tCreateBitRep) then
            if(allocated(Np1BitList)) deallocate(Np1BitList)
            allocate(Np1BitList(nNp1FCIDet))
            do i=1,nNp1FCIDet
                call EncodeBitDet(Np1FCIDetList(:,i),NEl+1,Np1BitList(i))
            enddo
        endif
    endif
    DetListStorage = DetListStorage + ((NEl+2)*nNp1FCIDet)
    if(tSplitMs) DetListStorage = DetListStorage + ((NEl+2)*nNp1bFCIDet)

    if(tSplitMs) then
        ECoupledSpace = nFCIDet + nNm1FCIDet + nNm1bFCIDet + nNp1FCIDet + nNp1bFCIDet
    else
        ECoupledSpace = nFCIDet + nNm1FCIDet + nNp1FCIDet
    endif

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
    integer :: Ex(2,2),iGetExcitLevel,IC,DeltaSpin
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
        DeltaSpin = mod(Ex(1,1),2) + mod(Ex(2,1),2)
        if(mod(DeltaSpin,2).eq.1) then
            !We have changed Ms - forbidden excitation
            HEl = 0.0_dp
            return
        endif
        hel = sltcnd_1(nI, Ex, tSign, NEl)
    elseif(IC.eq.2) then
        DeltaSpin = mod(Ex(1,1),2) + mod(Ex(2,1),2) + mod(Ex(1,2),2) + mod(Ex(2,2),2)
        if(mod(DeltaSpin,2).eq.1) then
            !We have changed Ms - forbidden excitation
            HEl = 0.0_dp
            return
        endif
        hel = sltcnd_2(Ex, tSign)
    endif

end subroutine GetHElement

subroutine GetHElement_comp(nI,nJ,NEl,HEl)
    use const
    implicit none
    integer, intent(in) :: NEl
    integer, intent(in) :: nI(NEl),nJ(NEl)
    complex(dp), intent(out) :: HEl
    complex(dp) :: sltcnd_0_comp,sltcnd_1_comp,sltcnd_2_comp
    integer :: Ex(2,2),iGetExcitLevel,IC,DeltaSpin
    logical :: tSign

    IC = IGETEXCITLEVEL(NI,NJ,NEL)
    Ex(1,1) = IC
    if(Ex(1,1).le.2) then
        call GetExcitation(nI,nJ,NEl,Ex,tSign)
    else
        HEl = zzero
        return
    endif

    if(IC.eq.0) then
        !Diagonal hamiltonian matrix element
        !Sum in 1 electron terms
        hel = sltcnd_0_comp(nI,NEl)
    elseif(IC.eq.1) then
        !Single excitation
        DeltaSpin = mod(Ex(1,1),2) + mod(Ex(2,1),2)
        if(mod(DeltaSpin,2).eq.1) then
            !We have changed Ms - forbidden excitation
            HEl = zzero
            return
        endif
        hel = sltcnd_1_comp(nI, Ex, tSign, NEl)
    elseif(IC.eq.2) then
        DeltaSpin = mod(Ex(1,1),2) + mod(Ex(2,1),2) + mod(Ex(1,2),2) + mod(Ex(2,2),2)
        if(mod(DeltaSpin,2).eq.1) then
            !We have changed Ms - forbidden excitation
            HEl = zzero
            return
        endif
        hel = sltcnd_2_comp(Ex, tSign)
    endif

end subroutine GetHElement_comp

function sltcnd_0_comp (nI,NEl) result(hel)
    use const
    use DetToolsData, only: TMat_comp
    implicit none
    integer, intent(in) :: NEl
    integer, intent(in) :: nI(NEl)
    integer :: id(nel),gtid,i,idX,idN,j
    complex(dp) :: hel
    real(dp) :: umatel

    do i=1,NEl
        id(i) = gtID(nI(i))   !turn to spatial orbital representation
    enddo

    hel = zzero
    do i=1,NEl
        hel = hel + tmat_comp(id(i),id(i))
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

end function sltcnd_0_comp

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

function sltcnd_1_comp (nI, ex, tSign, nel) result(hel)
    use const
    use DetToolsData, only: TMat_comp
    implicit none
    integer, intent(in) :: nel,nI(nel),ex(2)
    logical, intent(in) :: tSign
    real(dp) :: umatel
    complex(dp) :: hel
    integer :: id_ex(2), id, i,gtid

    do i=1,2
        id_ex(i) = gtid(ex(i))
    enddo

    hel = zzero
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
        hel = hel + tmat_comp(id_ex(1),id_ex(2))
    endif
    if(tSign) hel = -hel

end function sltcnd_1_comp

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

function sltcnd_2_comp (ex, tSign) result(hel)
    use const
    implicit none
    integer, intent(in) :: ex(2,2)
    logical, intent(in) :: tSign
    real(dp) :: umatel
    complex(dp) :: hel
    integer :: id(2,2),gtid,i,j

    do i=1,2
        do j=1,2
            id(j,i) = gtid(ex(j,i))
        enddo
    enddo

    if((mod(ex(1,1),2).eq.mod(ex(2,1),2)).and.    &
       (mod(ex(1,2),2).eq.mod(ex(2,2),2))) then
        hel = dcmplx(umatel(id(1,1),id(1,2),id(2,1),id(2,2)),0.0_dp)
    else
        hel = zzero
    endif

    if((mod(ex(1,1),2).eq.mod(ex(2,2),2)).and.    &
       (mod(ex(1,2),2).eq.mod(ex(2,1),2))) then
        hel = hel - dcmplx(umatel(id(1,1),id(1,2),id(2,2),id(2,1)),0.0_dp)
    endif
    if(tSign) hel = -hel

end function sltcnd_2_comp

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
    use Globals, only: U, nSites, tAnderson
    implicit none
    integer, intent(in) :: ex(2,2)
    real(dp), intent(in) :: AOBasisTrans(nSites,nSites)
    real(dp) :: HEl,HEl_coul,HEl_exch
    integer :: i,j,k,l,i_spat,j_spat,k_spat,l_spat,gtid,alpha
    integer :: Corr_sites

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
        
        if(tAnderson) then
            Corr_sites = 1
        else
            !Hubbard model
            Corr_sites = nSites
        endif
        do alpha = 1,Corr_sites
            HEl_coul = HEl_coul + AOBasisTrans(alpha,i_spat)*AOBasisTrans(alpha,j_spat)*    &
                AOBasisTrans(alpha,k_spat)*AOBasisTrans(alpha,l_spat)
        enddo
        HEl_coul = HEl_coul * U
    endif

    !Now, the exchange component if spin allowed <i j | l k>
    if((mod(i,2).eq.mod(l,2)).and.(mod(j,2).eq.mod(k,2))) then
        !Integral is allowed. It will of course be the same as the coulomb term
        if(tAnderson) then
            Corr_sites = 1
        else
            !Hubbard model
            Corr_sites = nSites
        endif
        do alpha = 1,Corr_sites
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
    use Globals, only: U, nSites, tAnderson
    implicit none
    integer, intent(in) :: ex(2,2)
    real(dp), intent(in) :: AOBasisTrans(nSites,nSites)
    real(dp) :: HEl
    integer :: i,j,k,l,i_spat,j_spat,k_spat,l_spat,gtid,alpha
    integer :: Corr_sites

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

        if(tAnderson) then
            Corr_sites = 1
        else
            !Hubbard model
            Corr_sites = nSites
        endif
        do alpha = 1,Corr_sites
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
subroutine getexcitation(ni,nj,nel,ex,tsign)
    implicit none
    integer, intent(in) :: nel,ni(nel),nj(nel)
    integer, intent(inout) :: ex(2,*)
    logical, intent(out) :: tsign
    integer imaxexcit
    integer i,j,ipar
    integer ic1,ic2
    imaxexcit=ex(1,1)
    ex(1:2,1:imaxexcit)=0
    ic1=0
    ic2=0
    i=1
    j=1
    ipar=0
!    call writedet(6,ni,nel,.true.)
!    call writedet(6,nj,nel,.true.)
    do while(i.le.nel.and.j.le.nel)
!.. differences from i to j
!   write(6,*) "ge",i,j
        do while(i.le.nel)
           if (ni(i) >= nj(j)) exit
           ic1=ic1+1
           if(ic1.le.imaxexcit) then
              ex(1,ic1)=ni(i)
              ipar=ipar+i
           endif
           i=i+1
        enddo
!.. differences from j to i
        do while(i.le.nel.and.j.le.nel)
           if (ni(i) <= nj(j)) exit
           ic2=ic2+1
           if(ic2.le.imaxexcit) then
              ex(2,ic2)=nj(j)
              ipar=ipar+j
           endif
           j=j+1
        enddo
        if(i.le.nel.and.j.le.nel) then
           if (ni(i) == nj(j)) then
               i=i+1
               j=j+1
           endif
        endif
    enddo
!.. deal with remaining i
    do while(i.le.nel)
        ic1=ic1+1
        if(ic1.le.imaxexcit) then
           ipar=ipar+i
           ex(1,ic1)=ni(i)
        endif
        i=i+1
    enddo
!.. deal with remaining j
    do while(j.le.nel)
        ic2=ic2+1
        if(ic2.le.imaxexcit) then
           ex(2,ic2)=nj(j)
           ipar=ipar+j
        endif
        j=j+1
    enddo
    if(ic1.gt.imaxexcit) then
!.. we actually needed more space.  just list the excitation counts (-ve)
       do i=1,imaxexcit
          if (i.eq.1) then
             ex(1,1)=-ic1
             ex(2,1)=-ic2
          else
             ex(1,i)=0
             ex(2,i)=0
          endif
       enddo
    elseif(ic1.eq.0) then
       ex(1,1)=0
       ex(2,1)=0
    endif
    tsign=btest(ipar,0)
    return
end subroutine

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

!Ensure that all LVecs are orthonormal with respect to the RVec vectors
!LVec vectors need to be complex conjugated.
!This will involve gram-schmidt rotation of the vectors in degenerate sets
subroutine Orthonorm_zgeev_vecs(N,W,LVec,RVec)
    use const
    implicit none
    integer, intent(in) :: N
    complex(dp), intent(in) :: W(N)
    complex(dp), intent(inout) :: LVec(N,N),RVec(N,N)
    integer :: i,StartingInd,R,L
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

        !The degenerate set goes from StartingInd to i
        !Orthogonalize the R vectors against the L vectors
        do R = StartingInd,i    !Run through the R vectors
            do L = StartingInd,i    !Run through the L vectors
                if(R.eq.L) cycle    !We want to maintain overlap with the corresponding vector

                !Gram-Schmidt
                overlap = zdotc(N,LVec(:,L),1,RVec(:,R),1)
                RVec(:,R) = RVec(:,R) - overlap * dconjg(LVec(:,L))
            enddo
        enddo
            
        do R = StartingInd,i
            !Now, normalize the vectors
            !Remember that the square root of a complex number will have two roots, given by +- w
            norm = zdotc(N,LVec(:,R),1,RVec(:,R),1)
            norm = sqrt(norm)
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
