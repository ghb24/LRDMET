program matsusum
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)   !For double precision real numbers
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
    complex(dp), parameter :: zzero = cmplx(0.0_dp,0.0_dp,dp)
    complex(dp), parameter :: zone = cmplx(1.0_dp,0.0_dp,dp)
    real(dp), parameter :: pi = 3.1415926535897931_dp

    integer :: n=2    
    real(dp) :: ScaleImTime = 10.0_dp
    integer :: nMatsu = 500
    real(8) :: beta = 100.0_8

    integer :: nImTimePoints
    real(dp), allocatable :: MatsuPoints(:),ImTimePoints(:)
    real(dp) :: Delta_Tau
    real(dp) :: h0(2,2),NEl
    real(dp), allocatable :: W(:),Orbs(:,:),Work(:)
    complex(dp) :: eye(2,2)
    integer :: lWork,info
    real(dp) :: ChemPot
    complex(dp), allocatable :: G_iw(:,:,:),SingleFreq(:,:),GF_Tau(:,:,:)
    integer :: i,j,a

    if(mod(nMatsu,2).eq.1) nMatsu = nMatsu + 1  !Should always be even to get same final point on both \pm w
    !The first 1:nMatsu/2 are the negative frequencies.
    !The next nMatsu/2+1:nMatsu are the positive frequencies
    allocate(MatsuPoints(nMatsu))
    j = 0
    do i = -nMatsu/2,(nMatsu/2)-1
        j = j + 1
        MatsuPoints(j) = (2*i+1)*pi / beta
        write(6,*) j, MatsuPoints(j) 
    enddo
    if(j.ne.nMatsu) stop 'error here'
    
    !do i = 0,nMatsu-1
    !    MatsuPoints(i+1) = (2*i+1)*pi / beta
    !enddo
    write(6,"(A,I7)") "Number of matsubara frequency points set to: ",nMatsu
    write(6,"(A,F20.10)") "Highest matsubara frequency point: ",MatsuPoints(nMatsu)
        
    nImTimePoints = ScaleImTime * real(nMatsu,dp) * Beta / pi
    allocate(ImTimePoints(nImTimePoints))
    write(6,"(A,I7)") "Number of imag-time points set to: ",nImTimePoints

    Delta_Tau = Beta/(nImTimePoints-1)
    write(6,"(A,G16.5)") "Imaginary time step in grid: ",Delta_Tau
    ImTimePoints(:) = zero
    do i = 2,nImTimePoints
        ImTimePoints(i) = (i-1)*Delta_Tau
    enddo

    !Form hubbard model
    h0(1,1) = zero
    h0(2,2) = zero
    h0(1,2) = -one
    h0(2,1) = -one

    !diagonalize
    allocate(Orbs(n,n))
    Orbs(:,:) = h0(:,:)
    allocate(W(n))
    W(:) = 0.0_dp
    allocate(Work(1))
    lWork=-1
    info=0
    call dsyev('V','L',n,Orbs,n,W,Work,lWork,info)
    if(info.ne.0) stop 'workspace' 
    lwork=int(work(1))+1
    deallocate(work)
    allocate(work(lwork))
    call dsyev('V','L',n,Orbs,n,W,Work,lWork,info)
    if(info.ne.0) stop 'Diag failed'
    deallocate(work)

    write(6,*) "Eigenvalues: ",W(1),W(2)

    !Half filled
    ChemPot = (W(1)+W(2))/2.0_dp
    !Third of way between gap
    !ChemPot = (2.0*W(1)/3.0) + (1.0*W(2)/3.0)
    write(6,*) "Chemical potential: ",ChemPot

    allocate(G_iw(2,2,nMatsu))
    G_iw(:,:,:) = zzero
    !Form G0(iw)
    do i = 1,nMatsu
        allocate(SingleFreq(2,2))
        SingleFreq(:,:) = -h0(:,:)
        do j = 1,n
            SingleFreq(j,j) = SingleFreq(j,j) + cmplx(ChemPot,MatsuPoints(i),dp)
        enddo
        !invert
        call z_inv2(SingleFreq,G_iw(:,:,i),n)
        deallocate(SingleFreq)
    enddo

    open(88,file='G_iw',status='unknown')
    do i = 1,nMatsu
        write(88,*) MatsuPoints(i),real(G_iw(1,1,i)),aimag(G_iw(1,1,i)),real(G_iw(2,1,i)),  &
            aimag(G_iw(2,1,i)),real(G_iw(2,2,i)),aimag(G_iw(2,2,i)),-one/MatsuPoints(i),    &
            aimag(G_iw(1,1,i)-(zone/cmplx(zero,MatsuPoints(i),dp)))
    enddo
    close(88)

    !Now, try to FT
    eye = zzero
    do i = 1,n
        eye(i,i) = zone
    enddo
    allocate(GF_Tau(n,n,nImTimePoints))
    GF_Tau(:,:,:) = zzero
    do a = 1,nImTimePoints

        do i = 1,nMatsu
            GF_Tau(:,:,a) = GF_Tau(:,:,a) + exp(-cmplx(zero,MatsuPoints(i)*ImTimePoints(a),dp)) * &
                (G_iw(:,:,i) - eye(:,:)*(one/cmplx(zero,MatsuPoints(i),dp)))
        enddo
        !do i = 1,nMatsu
        !    GF_Tau(:,:,a) = GF_Tau(:,:,a) + &
        !        exp(cmplx(zero,MatsuPoints(i)*ImTimePoints(a),dp)) * &
        !        (-G_iw(:,:,i) + eye(:,:)*(one/cmplx(zero,MatsuPoints(i),dp)))
        !enddo

        GF_Tau(:,:,a) = GF_Tau(:,:,a) / Beta

        GF_Tau(:,:,a) = GF_Tau(:,:,a) - 0.5_dp*eye(:,:)
    enddo

    open(88,file='G_tau',status='unknown')
    do i = 1,nImTimePoints
        write(88,*) ImTimePoints(i),real(GF_Tau(1,1,i)),aimag(GF_Tau(1,1,i)),real(GF_Tau(2,1,i)),   &
            aimag(GF_Tau(2,1,i)),real(GF_Tau(2,2,i)),aimag(GF_Tau(2,2,i))
    enddo
    close(88)

    NEl = zero
    do i = 1,n
        NEl = NEl + GF_Tau(i,i,nImTimePoints)
    enddo
    write(6,*) "Number of electrons in system: ",-2.0_dp*NEl

    !Now, transform back and see if it matches
    G_iw(:,:,:) = zzero
    do i = 1,nMatsu
        do j = 2,nImTimePoints-1
            G_iw(:,:,i) = G_iw(:,:,i) + exp(cmplx(zero,MatsuPoints(i)*ImTimePoints(j),dp))*GF_Tau(:,:,j)
        enddo
        G_iw(:,:,i) = G_iw(:,:,i) + (GF_Tau(:,:,1) - GF_Tau(:,:,nImTimePoints))/2.0_dp
        G_iw(:,:,i) = Delta_Tau*G_iw(:,:,i)
    enddo

    open(88,file='G_iw_Trans',status='unknown')
    do i = 1,nMatsu
        write(88,*) MatsuPoints(i),real(G_iw(1,1,i)),aimag(G_iw(1,1,i)),real(G_iw(2,1,i)),  &
            aimag(G_iw(2,1,i)),real(G_iw(2,2,i)),aimag(G_iw(2,2,i)),-one/MatsuPoints(i),    &
            aimag(G_iw(1,1,i)-(zone/cmplx(zero,MatsuPoints(i),dp)))
    enddo
    close(88)

end program

  SUBROUTINE z_inv2(mat,matinv,dimen)
      implicit none
    integer, parameter :: dp = selected_real_kind(15,307)   !For double precision real numbers
    integer, intent(in) :: dimen
    complex(dp), intent(in) :: mat(dimen,dimen)
    complex(dp), dimension(dimen,dimen), intent(out) :: matinv
    integer, dimension(dimen) :: ipiv
    integer :: msize,nsize,lwork,info
    complex(dp), allocatable :: cWork(:)

    msize=dimen
    nsize=dimen
    if(msize.ne.nsize) stop 'error in z_inv2'
    matinv(:,:) = mat(:,:)

    info=0
    call ZGETRF(msize,nsize,matinv,nsize,ipiv,info)
    IF (INFO /= 0) then
        write(6,*) "info: ",info
        STOP 'Error with z_inv2 1'
    endif
    allocate(cWork(1))
    lwork = -1
    call ZGETRI(msize,matinv,msize,ipiv,cwork,lwork,info)
    if(info.ne.0) stop 'error with workspace query in z_inv2'
    lWork = int(real(cWork(1))) + 1
    deallocate(cWork)
    allocate(cWork(lWork))
    call ZGETRI(msize,matinv,msize,ipiv,cwork,lwork,info)
    if(info.ne.0) stop 'error with inversion in z_inv2'
    deallocate(cWork)

  end subroutine z_inv2
